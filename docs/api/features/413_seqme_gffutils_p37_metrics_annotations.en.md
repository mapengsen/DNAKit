# seqme sequence metrics and gffutils annotations

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: seqme 0.5.1 was tested from a wheel extracted into scratch. RNA-FM/GENA-LM weights and hosted resources were not downloaded. The upstream FlyBase conversion script is excluded because its exon loop writes stale coordinates/transcript identifiers.

## `dnakit.predictions.seqme_p37_functions.seqme_authenticity`

**Purpose:** Authenticity (AuthPct): the fraction of generated sequences whose nearest training neighbour is closer to another training sequence than to the generated one (Alaa et al. 2022).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_authenticity(generated: Sequence[str], train_set: Sequence[str], *, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Authenticity (AuthPct): the fraction of generated sequences whose nearest training neighbour is closer to another training sequence than to the generated one (Alaa et al. 2022).

Higher is better; low values mean the generator copies its training data. k-mer embeddings of ``kmer_size`` are used.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_clipped_density`

**Purpose:** Clipped density, an outlier-robust fidelity score of generated sequences against a reference manifold built from k-nearest-neighbour radii (Salvy et al. 2024). Range 0 to 1, higher is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_clipped_density(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Clipped density, an outlier-robust fidelity score of generated sequences against a reference manifold built from k-nearest-neighbour radii (Salvy et al. 2024). Range 0 to 1, higher is better.

``strict=True`` requires equally many generated and reference sequences.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_clipped_coverage`

**Purpose:** Clipped coverage, an outlier-robust diversity score: how much of the reference manifold the generated sequences cover. Range 0 to 1, higher is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_clipped_coverage(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Clipped coverage, an outlier-robust diversity score: how much of the reference manifold the generated sequences cover. Range 0 to 1, higher is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_precision`

**Purpose:** Improved precision (Kynkaanniemi et al. 2019): the share of generated sequences inside the k-NN manifold of the reference set (fidelity). Higher is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_precision(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Improved precision (Kynkaanniemi et al. 2019): the share of generated sequences inside the k-NN manifold of the reference set (fidelity). Higher is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_recall`

**Purpose:** Improved recall (Kynkaanniemi et al. 2019): the share of reference sequences inside the k-NN manifold of the generated set (coverage). Higher is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_recall(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Improved recall (Kynkaanniemi et al. 2019): the share of reference sequences inside the k-NN manifold of the generated set (coverage). Higher is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_conformity_score`

**Purpose:** Distributional conformity score: how well the property distribution of the generated set fits a kernel density estimate of the reference set, with cross-validation (Gaussian KDE).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_conformity_score(generated: Sequence[str], reference: Sequence[str], *, property_names: Sequence[str]=('gc_content', 'length'), n_splits: int=5, kde_bandwidth: str | float='silverman', seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Distributional conformity score: how well the property distribution of the generated set fits a kernel density estimate of the reference set, with cross-validation (Gaussian KDE).

``property_names`` are built-in sequence properties (gc_content, length, max_homopolymer_length, cpg_fraction).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_fkea_diversity`

**Purpose:** FKEA kernel-entropy diversity (Vendi / RKE score): the effective number of distinct modes of the generated k-mer embeddings under a Gaussian kernel (Friedman and Dieng 2023; Ospanov et al. 2024).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_fkea_diversity(generated: Sequence[str], *, bandwidth: float=1.0, alpha: float=2, n_random_fourier_features: int | None=2048, kmer_size: int=3, seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
FKEA kernel-entropy diversity (Vendi / RKE score): the effective number of distinct modes of the generated k-mer embeddings under a Gaussian kernel (Friedman and Dieng 2023; Ospanov et al. 2024).

``alpha=2`` gives the RKE score, other values the Vendi-alpha score. ``n_random_fourier_features=None`` uses the exact kernel.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_hit_rate`

**Purpose:** Hit rate: the proportion of generated sequences whose property lies in ``[min_value, max_value]`` (a bound left as None is open), with its standard error.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_hit_rate(generated: Sequence[str], property_name: str, *, min_value: float | None=None, max_value: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Hit rate: the proportion of generated sequences whose property lies in ``[min_value, max_value]`` (a bound left as None is open), with its standard error.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_hypervolume`

**Purpose:** Multi-objective hypervolume of the property values of the generated set (all properties maximised): the hypervolume indicator (``method='hvi'``) or the convex-hull volume (``'convex-hull'``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_hypervolume(generated: Sequence[str], property_names: Sequence[str], *, method: str='hvi', nadir: Sequence[float] | None=None, ideal: Sequence[float] | None=None, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Multi-objective hypervolume of the property values of the generated set (all properties maximised): the hypervolume indicator (``method='hvi'``) or the convex-hull volume (``'convex-hull'``).

``nadir`` (worst value per property, default zeros) and ``ideal`` (best value, used for scaling to [0, 1]) have one number per property.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_property_mean`

**Purpose:** Mean of a sequence property over the generated set with its deviation (the seqme ``ID`` metric); ``deviation`` is ``std``, ``se`` or ``var``, ``estimate`` ``biased`` or ``unbiased``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_property_mean(generated: Sequence[str], property_name: str, *, objective: str='maximize', deviation: str='se', estimate: str='unbiased', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Mean of a sequence property over the generated set with its deviation (the seqme ``ID`` metric); ``deviation`` is ``std``, ``se`` or ``var``, ``estimate`` ``biased`` or ``unbiased``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_property_kl_divergence`

**Purpose:** KL divergence between the property distribution of the generated set and the reference set, estimated with Gaussian KDEs and ``n_draws`` Monte Carlo samples. Lower is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_property_kl_divergence(generated: Sequence[str], reference: Sequence[str], property_name: str, *, n_draws: int=10000, kde_bandwidth: str | float='silverman', seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
KL divergence between the property distribution of the generated set and the reference set, estimated with Gaussian KDEs and ``n_draws`` Monte Carlo samples. Lower is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_mmd`

**Purpose:** Maximum mean discrepancy with a Gaussian RBF kernel between the k-mer embeddings of the generated and reference sets, multiplied by ``scale``. Lower is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_mmd(generated: Sequence[str], reference: Sequence[str], *, estimate: str='biased', sigma: float=10, scale: float=1000, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Maximum mean discrepancy with a Gaussian RBF kernel between the k-mer embeddings of the generated and reference sets, multiplied by ``scale``. Lower is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_kid`

**Purpose:** Kernel inception distance: MMD with a polynomial kernel (``degree``, ``coef0``) between k-mer embeddings (Binkowski et al. 2018). Lower is better.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_kid(generated: Sequence[str], reference: Sequence[str], *, estimate: str='biased', degree: int=3, coef0: float=1.0, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Kernel inception distance: MMD with a polynomial kernel (``degree``, ``coef0``) between k-mer embeddings (Binkowski et al. 2018). Lower is better.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_threshold_fraction`

**Purpose:** Fraction of generated sequences whose property passes a threshold: at least ``min_value`` and/or at most ``max_value`` (give one or both).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_threshold_fraction(generated: Sequence[str], property_name: str, *, min_value: float | None=None, max_value: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Fraction of generated sequences whose property passes a threshold: at least ``min_value`` and/or at most ``max_value`` (give one or both).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_embed_rna_fm`

**Purpose:** Embed RNA/mRNA sequences with RNA-FM (``mRNA`` codon-aligned CDS model, 1280 dimensions, or ``ncRNA``, 640 dimensions) by averaging token embeddings of ``layer`` (last layer is 12).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_embed_rna_fm(sequences: Sequence[str], *, model_name: str='mRNA', layer: int=12, batch_size: int=32, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Embed RNA/mRNA sequences with RNA-FM (``mRNA`` codon-aligned CDS model, 1280 dimensions, or ``ncRNA``, 640 dimensions) by averaging token embeddings of ``layer`` (last layer is 12).

Needs the ``rna-fm`` package and downloads the model weights on first use. Sequences are RNA (U) or DNA (T) letters as RNA-FM expects.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_embed_genalm`

**Purpose:** Embed DNA sequences with a GENA-LM language-model checkpoint (mean of the token embeddings of ``layer``, default the last).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_embed_genalm(sequences: Sequence[str], *, checkpoint: str='bert_base_t2t', layer: int=-1, batch_size: int=16, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Embed DNA sequences with a GENA-LM language-model checkpoint (mean of the token embeddings of ``layer``, default the last).

``checkpoint`` is one of bert_base_t2t, bert_base_t2t_lastln_t2t, bert_base_t2t_multi, bert_large_t2t or bigbird_base_t2t.
Needs ``transformers`` and ``einops`` and downloads the model from Hugging Face on first use (network).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_classify_genalm`

**Purpose:** Classify DNA sequences with a fine-tuned GENA-LM model and return class probabilities.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_classify_genalm(sequences: Sequence[str], *, checkpoint: str='bert_base_t2t_promoters', batch_size: int=16, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Classify DNA sequences with a fine-tuned GENA-LM model and return class probabilities.

``checkpoint`` is a promoter model (bert_base_t2t_promoters, bert_large_t2t_promoters at 300 bp, bert_large_t2t_promoters2 at 2000 bp;
classes: no promoter, promoter) or a splice-site model (bert_base_t2t_splice_site, bert_large_t2t_splice_site at 15000 bp; classes:
neither, splice acceptor, splice donor). Needs ``transformers`` and downloads the model from Hugging Face on first use (network).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_project_embeddings`

**Purpose:** Project embeddings to low dimension by PCA, t-SNE or UMAP (``seqme.utils.pca``, ``tsne``, ``umap``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_project_embeddings(embeddings: Sequence[Sequence[float]], *, method: str='pca', seed: int | None=0, reference_embeddings: Sequence[Sequence[float]] | None=None, n_components: int=2, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Project embeddings to low dimension by PCA, t-SNE or UMAP (``seqme.utils.pca``, ``tsne``, ``umap``).

``method`` is ``pca``, ``tsne`` or ``umap``; each returns 2-D points (t-SNE with PCA initialisation, perplexity and
learning rate as in seqme). With ``method='pca'`` and ``reference_embeddings``, PCA is fitted on the reference
rows and applied to ``embeddings`` with ``n_components`` components (``seqme.models.PCA``); the explained variance ratio is returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_spearman_alignment_score`

**Purpose:** Spearman correlation between the pairwise-distance matrices of two embeddings of the same sequences (alignment of two embedding spaces, -1 to 1).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_spearman_alignment_score(embeddings_a: Sequence[Sequence[float]], embeddings_b: Sequence[Sequence[float]], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Spearman correlation between the pairwise-distance matrices of two embeddings of the same sequences (alignment of two embedding spaces, -1 to 1).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_knn_alignment_score`

**Purpose:** k-NN feature alignment score of an embedding model: the share of each sequence's ``n_neighbors`` nearest neighbours that carry the same label (0 to 1).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_knn_alignment_score(embeddings: Sequence[Sequence[float]], labels: Sequence[Any], *, n_neighbors: int=5, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
k-NN feature alignment score of an embedding model: the share of each sequence's ``n_neighbors`` nearest neighbours that carry the same label (0 to 1).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_plot_embeddings`

**Purpose:** Scatter-plot 2-D embedding projections, one coloured group per key (``seqme.utils.plot_embeddings``), and save the image (png, svg or pdf).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_plot_embeddings(projections: Mapping[str, Sequence[Sequence[float]]], output_path: str | Path, *, title: str | None=None, point_size: float=20, alpha: float=0.6, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Scatter-plot 2-D embedding projections, one coloured group per key (``seqme.utils.plot_embeddings``), and save the image (png, svg or pdf).

``projections`` maps a group name to its list of ``[x, y]`` points (for example from ``seqme_project_embeddings``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.predictions.seqme_p37_functions.seqme_plot_knn_alignment_score`

**Purpose:** Plot the k-NN feature alignment score against the number of neighbours (``seqme.utils.plot_knn_alignment_score``) and save the image; the scores are returned too.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def seqme_plot_knn_alignment_score(embeddings: Sequence[Sequence[float]], labels: Sequence[Any], output_path: str | Path, *, n_neighbors: Sequence[int]=(1, 3, 5, 10), label: str='embedding', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Plot the k-NN feature alignment score against the number of neighbours (``seqme.utils.plot_knn_alignment_score``) and save the image; the scores are returned too.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.fetch_gffutils_feature`

**Purpose:** Fetch one annotation feature by its ID from a GFF/GTF file or gffutils database (``gffutils-cli fetch``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def fetch_gffutils_feature(source: str | Path, feature_id: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Fetch one annotation feature by its ID from a GFF/GTF file or gffutils database (``gffutils-cli fetch``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.list_gffutils_features`

**Purpose:** List annotation features, all of them or one feature type, optionally on one sequence.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def list_gffutils_features(source: str | Path, *, featuretype: str | None=None, seqid: str | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List annotation features, all of them or one feature type, optionally on one sequence.

``featuretype=None`` walks ``FeatureDB.all_features()``; a type name walks ``features_of_type()``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.group_gffutils_parent_children`

**Purpose:** Group every feature of ``parent_featuretype`` with its child features (``FeatureDB.iter_by_parent_childs``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def group_gffutils_parent_children(source: str | Path, *, parent_featuretype: str='gene', level: int | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Group every feature of ``parent_featuretype`` with its child features (``FeatureDB.iter_by_parent_childs``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.fetch_gffutils_parents`

**Purpose:** Fetch the parent features of a feature (``FeatureDB.parents``, ``gffutils-cli parents``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def fetch_gffutils_parents(source: str | Path, feature_id: str, *, level: int | None=None, featuretype: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Fetch the parent features of a feature (``FeatureDB.parents``, ``gffutils-cli parents``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.get_gffutils_interfeatures`

**Purpose:** Return the gaps between consecutive features of one type (``FeatureDB.interfeatures``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def get_gffutils_interfeatures(source: str | Path, *, featuretype: str='gene', new_featuretype: str | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Return the gaps between consecutive features of one type (``FeatureDB.interfeatures``).

Features are sorted by sequence and start; gaps are made only between neighbours on the same sequence.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.delete_gffutils_features`

**Purpose:** Delete features by ID or all features of one type; the edited database is written to a new ``.db`` file.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def delete_gffutils_features(source: str | Path, output_path: str | Path, *, feature_ids: Sequence[str] | None=None, featuretype: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Delete features by ID or all features of one type; the edited database is written to a new ``.db`` file.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.update_gffutils_features`

**Purpose:** Add GFF lines to an annotation database (``FeatureDB.update``); the result is a new ``.db`` file.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def update_gffutils_features(source: str | Path, additions: Sequence[str], output_path: str | Path, *, merge_strategy: str='error', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Add GFF lines to an annotation database (``FeatureDB.update``); the result is a new ``.db`` file.

``merge_strategy`` is a gffutils strategy: ``error``, ``warning``, ``merge``, ``create_unique`` or ``replace``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.add_gffutils_relation`

**Purpose:** Add a parent-child relation of a given level between two features (``FeatureDB.add_relation``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def add_gffutils_relation(source: str | Path, parent_id: str, child_id: str, level: int, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Add a parent-child relation of a given level between two features (``FeatureDB.add_relation``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.create_gffutils_introns`

**Purpose:** Derive intron features from the gaps between the exons of each transcript (``FeatureDB.create_introns``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def create_gffutils_introns(source: str | Path, *, exon_featuretype: str='exon', grandparent_featuretype: str | None='gene', parent_featuretype: str | None=None, new_featuretype: str='intron', merge_attributes: bool=True, numeric_sort: bool=False, output_path: str | Path | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Derive intron features from the gaps between the exons of each transcript (``FeatureDB.create_introns``).

The new features are returned and, when ``output_path`` (a new file) is given, written as GFF3.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.create_gffutils_splice_sites`

**Purpose:** Derive 5' and 3' splice-site features from the exon boundaries of each transcript (``FeatureDB.create_splice_sites``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def create_gffutils_splice_sites(source: str | Path, *, exon_featuretype: str='exon', grandparent_featuretype: str | None='gene', parent_featuretype: str | None=None, merge_attributes: bool=True, numeric_sort: bool=False, output_path: str | Path | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Derive 5' and 3' splice-site features from the exon boundaries of each transcript (``FeatureDB.create_splice_sites``).

Every exon must carry an ``ID`` attribute, as gffutils builds the new IDs from them.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.gffutils_children_bp`

**Purpose:** Total base pairs of the children of one type of a feature, for example exonic bp of a gene (``FeatureDB.children_bp``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def gffutils_children_bp(source: str | Path, feature_id: str, *, child_featuretype: str='exon', merge: bool=False, ignore_strand: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Total base pairs of the children of one type of a feature, for example exonic bp of a gene (``FeatureDB.children_bp``).

With ``merge=True`` overlapping children are merged first (same sequence, strand and type).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.extract_gffutils_feature_sequence`

**Purpose:** Extract the genome sequence of a feature from a FASTA file (``Feature.sequence``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def extract_gffutils_feature_sequence(source: str | Path, fasta: str | Path, feature_id: str, *, use_strand: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Extract the genome sequence of a feature from a FASTA file (``Feature.sequence``).

With ``use_strand=True`` a minus-strand feature is reverse-complemented, as gffutils does.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.write_gffutils_gff`

**Purpose:** Write features as GFF records with ``GFFWriter`` (``write_rec`` / ``write_recs``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def write_gffutils_gff(source: str | Path, output_path: str | Path, *, featuretype: str | None=None, canonical_gene_order: bool=False, gene_featuretype: str='gene', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Write features as GFF records with ``GFFWriter`` (``write_rec`` / ``write_recs``).

``canonical_gene_order=True`` instead writes each gene followed by its mRNAs (longest first), their
exons and the remaining children (``write_gene_recs``); ``featuretype`` is then ignored.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.remove_gffutils_duplicates`

**Purpose:** Remove duplicate GFF records by merging records that share an ID (``gffutils-cli rmdups``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def remove_gffutils_duplicates(source: str | Path, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Remove duplicate GFF records by merging records that share an ID (``gffutils-cli rmdups``).

The database is built with ``merge_strategy="merge"`` and every record is written once to ``output_path``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.gffutils_p37_functions.plot_gffutils_gene_model`

**Purpose:** Draw a gene model (transcript line, UTRs, CDS) with ``gffutils.contrib.plotting.Gene`` and save the image.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def plot_gffutils_gene_model(source: str | Path, gene_id: str, output_path: str | Path, *, transcripts: Sequence[str]=('mRNA',), utrs: Sequence[str]=("3'UTR", "5'UTR"), cds: Sequence[str]=('CDS',), verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Draw a gene model (transcript line, UTRs, CDS) with ``gffutils.contrib.plotting.Gene`` and save the image.

``utrs`` depends on the annotation: GTF uses ``3'UTR``/``5'UTR``, FlyBase GFF ``three_prime_UTR``/``five_prime_UTR``.
Needs pybedtools and matplotlib. The result holds the image path and the drawn transcripts.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.predictions.seqme_p37_functions import seqme_property_mean
r = seqme_property_mean(['ACGT', 'GGGG', 'AAAA'], 'gc_content', verbose=False)
print(r['value'], r['generated_count'])
```

**Example output:**

```text
0.5 3
```
