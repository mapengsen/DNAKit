# seqme 序列指标与 gffutils 注释

本页记录 P37 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：seqme 0.5.1 was tested from a wheel extracted into scratch. RNA-FM/GENA-LM weights and hosted resources were not downloaded. The upstream FlyBase conversion script is excluded because its exon loop writes stale coordinates/transcript identifiers.

## `dnakit.predictions.seqme_p37_functions.seqme_authenticity`

**作用：** Authenticity (AuthPct): the fraction of generated sequences whose nearest training neighbour is closer to another training sequence than to the generated one (Alaa et al. 2022).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_authenticity(generated: Sequence[str], train_set: Sequence[str], *, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Authenticity (AuthPct): the fraction of generated sequences whose nearest training neighbour is closer to another training sequence than to the generated one (Alaa et al. 2022).

Higher is better; low values mean the generator copies its training data. k-mer embeddings of ``kmer_size`` are used.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_clipped_density`

**作用：** Clipped density, an outlier-robust fidelity score of generated sequences against a reference manifold built from k-nearest-neighbour radii (Salvy et al. 2024). Range 0 to 1, higher is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_clipped_density(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Clipped density, an outlier-robust fidelity score of generated sequences against a reference manifold built from k-nearest-neighbour radii (Salvy et al. 2024). Range 0 to 1, higher is better.

``strict=True`` requires equally many generated and reference sequences.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_clipped_coverage`

**作用：** Clipped coverage, an outlier-robust diversity score: how much of the reference manifold the generated sequences cover. Range 0 to 1, higher is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_clipped_coverage(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Clipped coverage, an outlier-robust diversity score: how much of the reference manifold the generated sequences cover. Range 0 to 1, higher is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_precision`

**作用：** Improved precision (Kynkaanniemi et al. 2019): the share of generated sequences inside the k-NN manifold of the reference set (fidelity). Higher is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_precision(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Improved precision (Kynkaanniemi et al. 2019): the share of generated sequences inside the k-NN manifold of the reference set (fidelity). Higher is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_recall`

**作用：** Improved recall (Kynkaanniemi et al. 2019): the share of reference sequences inside the k-NN manifold of the generated set (coverage). Higher is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_recall(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int=5, kmer_size: int=3, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Improved recall (Kynkaanniemi et al. 2019): the share of reference sequences inside the k-NN manifold of the generated set (coverage). Higher is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_conformity_score`

**作用：** Distributional conformity score: how well the property distribution of the generated set fits a kernel density estimate of the reference set, with cross-validation (Gaussian KDE).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_conformity_score(generated: Sequence[str], reference: Sequence[str], *, property_names: Sequence[str]=('gc_content', 'length'), n_splits: int=5, kde_bandwidth: str | float='silverman', seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Distributional conformity score: how well the property distribution of the generated set fits a kernel density estimate of the reference set, with cross-validation (Gaussian KDE).

``property_names`` are built-in sequence properties (gc_content, length, max_homopolymer_length, cpg_fraction).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_fkea_diversity`

**作用：** FKEA kernel-entropy diversity (Vendi / RKE score): the effective number of distinct modes of the generated k-mer embeddings under a Gaussian kernel (Friedman and Dieng 2023; Ospanov et al. 2024).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_fkea_diversity(generated: Sequence[str], *, bandwidth: float=1.0, alpha: float=2, n_random_fourier_features: int | None=2048, kmer_size: int=3, seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
FKEA kernel-entropy diversity (Vendi / RKE score): the effective number of distinct modes of the generated k-mer embeddings under a Gaussian kernel (Friedman and Dieng 2023; Ospanov et al. 2024).

``alpha=2`` gives the RKE score, other values the Vendi-alpha score. ``n_random_fourier_features=None`` uses the exact kernel.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_hit_rate`

**作用：** Hit rate: the proportion of generated sequences whose property lies in ``[min_value, max_value]`` (a bound left as None is open), with its standard error.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_hit_rate(generated: Sequence[str], property_name: str, *, min_value: float | None=None, max_value: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Hit rate: the proportion of generated sequences whose property lies in ``[min_value, max_value]`` (a bound left as None is open), with its standard error.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_hypervolume`

**作用：** Multi-objective hypervolume of the property values of the generated set (all properties maximised): the hypervolume indicator (``method='hvi'``) or the convex-hull volume (``'convex-hull'``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_hypervolume(generated: Sequence[str], property_names: Sequence[str], *, method: str='hvi', nadir: Sequence[float] | None=None, ideal: Sequence[float] | None=None, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Multi-objective hypervolume of the property values of the generated set (all properties maximised): the hypervolume indicator (``method='hvi'``) or the convex-hull volume (``'convex-hull'``).

``nadir`` (worst value per property, default zeros) and ``ideal`` (best value, used for scaling to [0, 1]) have one number per property.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_property_mean`

**作用：** Mean of a sequence property over the generated set with its deviation (the seqme ``ID`` metric); ``deviation`` is ``std``, ``se`` or ``var``, ``estimate`` ``biased`` or ``unbiased``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_property_mean(generated: Sequence[str], property_name: str, *, objective: str='maximize', deviation: str='se', estimate: str='unbiased', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Mean of a sequence property over the generated set with its deviation (the seqme ``ID`` metric); ``deviation`` is ``std``, ``se`` or ``var``, ``estimate`` ``biased`` or ``unbiased``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_property_kl_divergence`

**作用：** KL divergence between the property distribution of the generated set and the reference set, estimated with Gaussian KDEs and ``n_draws`` Monte Carlo samples. Lower is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_property_kl_divergence(generated: Sequence[str], reference: Sequence[str], property_name: str, *, n_draws: int=10000, kde_bandwidth: str | float='silverman', seed: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
KL divergence between the property distribution of the generated set and the reference set, estimated with Gaussian KDEs and ``n_draws`` Monte Carlo samples. Lower is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_mmd`

**作用：** Maximum mean discrepancy with a Gaussian RBF kernel between the k-mer embeddings of the generated and reference sets, multiplied by ``scale``. Lower is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_mmd(generated: Sequence[str], reference: Sequence[str], *, estimate: str='biased', sigma: float=10, scale: float=1000, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Maximum mean discrepancy with a Gaussian RBF kernel between the k-mer embeddings of the generated and reference sets, multiplied by ``scale``. Lower is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_kid`

**作用：** Kernel inception distance: MMD with a polynomial kernel (``degree``, ``coef0``) between k-mer embeddings (Binkowski et al. 2018). Lower is better.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_kid(generated: Sequence[str], reference: Sequence[str], *, estimate: str='biased', degree: int=3, coef0: float=1.0, kmer_size: int=3, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Kernel inception distance: MMD with a polynomial kernel (``degree``, ``coef0``) between k-mer embeddings (Binkowski et al. 2018). Lower is better.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_threshold_fraction`

**作用：** Fraction of generated sequences whose property passes a threshold: at least ``min_value`` and/or at most ``max_value`` (give one or both).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_threshold_fraction(generated: Sequence[str], property_name: str, *, min_value: float | None=None, max_value: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fraction of generated sequences whose property passes a threshold: at least ``min_value`` and/or at most ``max_value`` (give one or both).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_embed_rna_fm`

**作用：** Embed RNA/mRNA sequences with RNA-FM (``mRNA`` codon-aligned CDS model, 1280 dimensions, or ``ncRNA``, 640 dimensions) by averaging token embeddings of ``layer`` (last layer is 12).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_embed_rna_fm(sequences: Sequence[str], *, model_name: str='mRNA', layer: int=12, batch_size: int=32, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Embed RNA/mRNA sequences with RNA-FM (``mRNA`` codon-aligned CDS model, 1280 dimensions, or ``ncRNA``, 640 dimensions) by averaging token embeddings of ``layer`` (last layer is 12).

Needs the ``rna-fm`` package and downloads the model weights on first use. Sequences are RNA (U) or DNA (T) letters as RNA-FM expects.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_embed_genalm`

**作用：** Embed DNA sequences with a GENA-LM language-model checkpoint (mean of the token embeddings of ``layer``, default the last).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_embed_genalm(sequences: Sequence[str], *, checkpoint: str='bert_base_t2t', layer: int=-1, batch_size: int=16, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Embed DNA sequences with a GENA-LM language-model checkpoint (mean of the token embeddings of ``layer``, default the last).

``checkpoint`` is one of bert_base_t2t, bert_base_t2t_lastln_t2t, bert_base_t2t_multi, bert_large_t2t or bigbird_base_t2t.
Needs ``transformers`` and ``einops`` and downloads the model from Hugging Face on first use (network).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_classify_genalm`

**作用：** Classify DNA sequences with a fine-tuned GENA-LM model and return class probabilities.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_classify_genalm(sequences: Sequence[str], *, checkpoint: str='bert_base_t2t_promoters', batch_size: int=16, cache_dir: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Classify DNA sequences with a fine-tuned GENA-LM model and return class probabilities.

``checkpoint`` is a promoter model (bert_base_t2t_promoters, bert_large_t2t_promoters at 300 bp, bert_large_t2t_promoters2 at 2000 bp;
classes: no promoter, promoter) or a splice-site model (bert_base_t2t_splice_site, bert_large_t2t_splice_site at 15000 bp; classes:
neither, splice acceptor, splice donor). Needs ``transformers`` and downloads the model from Hugging Face on first use (network).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_project_embeddings`

**作用：** Project embeddings to low dimension by PCA, t-SNE or UMAP (``seqme.utils.pca``, ``tsne``, ``umap``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_project_embeddings(embeddings: Sequence[Sequence[float]], *, method: str='pca', seed: int | None=0, reference_embeddings: Sequence[Sequence[float]] | None=None, n_components: int=2, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Project embeddings to low dimension by PCA, t-SNE or UMAP (``seqme.utils.pca``, ``tsne``, ``umap``).

``method`` is ``pca``, ``tsne`` or ``umap``; each returns 2-D points (t-SNE with PCA initialisation, perplexity and
learning rate as in seqme). With ``method='pca'`` and ``reference_embeddings``, PCA is fitted on the reference
rows and applied to ``embeddings`` with ``n_components`` components (``seqme.models.PCA``); the explained variance ratio is returned.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_spearman_alignment_score`

**作用：** Spearman correlation between the pairwise-distance matrices of two embeddings of the same sequences (alignment of two embedding spaces, -1 to 1).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_spearman_alignment_score(embeddings_a: Sequence[Sequence[float]], embeddings_b: Sequence[Sequence[float]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Spearman correlation between the pairwise-distance matrices of two embeddings of the same sequences (alignment of two embedding spaces, -1 to 1).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_knn_alignment_score`

**作用：** k-NN feature alignment score of an embedding model: the share of each sequence's ``n_neighbors`` nearest neighbours that carry the same label (0 to 1).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_knn_alignment_score(embeddings: Sequence[Sequence[float]], labels: Sequence[Any], *, n_neighbors: int=5, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
k-NN feature alignment score of an embedding model: the share of each sequence's ``n_neighbors`` nearest neighbours that carry the same label (0 to 1).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_plot_embeddings`

**作用：** Scatter-plot 2-D embedding projections, one coloured group per key (``seqme.utils.plot_embeddings``), and save the image (png, svg or pdf).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_plot_embeddings(projections: Mapping[str, Sequence[Sequence[float]]], output_path: str | Path, *, title: str | None=None, point_size: float=20, alpha: float=0.6, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Scatter-plot 2-D embedding projections, one coloured group per key (``seqme.utils.plot_embeddings``), and save the image (png, svg or pdf).

``projections`` maps a group name to its list of ``[x, y]`` points (for example from ``seqme_project_embeddings``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.predictions.seqme_p37_functions.seqme_plot_knn_alignment_score`

**作用：** Plot the k-NN feature alignment score against the number of neighbours (``seqme.utils.plot_knn_alignment_score``) and save the image; the scores are returned too.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def seqme_plot_knn_alignment_score(embeddings: Sequence[Sequence[float]], labels: Sequence[Any], output_path: str | Path, *, n_neighbors: Sequence[int]=(1, 3, 5, 10), label: str='embedding', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Plot the k-NN feature alignment score against the number of neighbours (``seqme.utils.plot_knn_alignment_score``) and save the image; the scores are returned too.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.fetch_gffutils_feature`

**作用：** Fetch one annotation feature by its ID from a GFF/GTF file or gffutils database (``gffutils-cli fetch``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_gffutils_feature(source: str | Path, feature_id: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fetch one annotation feature by its ID from a GFF/GTF file or gffutils database (``gffutils-cli fetch``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.list_gffutils_features`

**作用：** List annotation features, all of them or one feature type, optionally on one sequence.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_gffutils_features(source: str | Path, *, featuretype: str | None=None, seqid: str | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
List annotation features, all of them or one feature type, optionally on one sequence.

``featuretype=None`` walks ``FeatureDB.all_features()``; a type name walks ``features_of_type()``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.group_gffutils_parent_children`

**作用：** Group every feature of ``parent_featuretype`` with its child features (``FeatureDB.iter_by_parent_childs``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def group_gffutils_parent_children(source: str | Path, *, parent_featuretype: str='gene', level: int | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Group every feature of ``parent_featuretype`` with its child features (``FeatureDB.iter_by_parent_childs``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.fetch_gffutils_parents`

**作用：** Fetch the parent features of a feature (``FeatureDB.parents``, ``gffutils-cli parents``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fetch_gffutils_parents(source: str | Path, feature_id: str, *, level: int | None=None, featuretype: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fetch the parent features of a feature (``FeatureDB.parents``, ``gffutils-cli parents``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.get_gffutils_interfeatures`

**作用：** Return the gaps between consecutive features of one type (``FeatureDB.interfeatures``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def get_gffutils_interfeatures(source: str | Path, *, featuretype: str='gene', new_featuretype: str | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Return the gaps between consecutive features of one type (``FeatureDB.interfeatures``).

Features are sorted by sequence and start; gaps are made only between neighbours on the same sequence.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.delete_gffutils_features`

**作用：** Delete features by ID or all features of one type; the edited database is written to a new ``.db`` file.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def delete_gffutils_features(source: str | Path, output_path: str | Path, *, feature_ids: Sequence[str] | None=None, featuretype: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Delete features by ID or all features of one type; the edited database is written to a new ``.db`` file.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.update_gffutils_features`

**作用：** Add GFF lines to an annotation database (``FeatureDB.update``); the result is a new ``.db`` file.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def update_gffutils_features(source: str | Path, additions: Sequence[str], output_path: str | Path, *, merge_strategy: str='error', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Add GFF lines to an annotation database (``FeatureDB.update``); the result is a new ``.db`` file.

``merge_strategy`` is a gffutils strategy: ``error``, ``warning``, ``merge``, ``create_unique`` or ``replace``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.add_gffutils_relation`

**作用：** Add a parent-child relation of a given level between two features (``FeatureDB.add_relation``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def add_gffutils_relation(source: str | Path, parent_id: str, child_id: str, level: int, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Add a parent-child relation of a given level between two features (``FeatureDB.add_relation``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.create_gffutils_introns`

**作用：** Derive intron features from the gaps between the exons of each transcript (``FeatureDB.create_introns``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def create_gffutils_introns(source: str | Path, *, exon_featuretype: str='exon', grandparent_featuretype: str | None='gene', parent_featuretype: str | None=None, new_featuretype: str='intron', merge_attributes: bool=True, numeric_sort: bool=False, output_path: str | Path | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Derive intron features from the gaps between the exons of each transcript (``FeatureDB.create_introns``).

The new features are returned and, when ``output_path`` (a new file) is given, written as GFF3.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.create_gffutils_splice_sites`

**作用：** Derive 5' and 3' splice-site features from the exon boundaries of each transcript (``FeatureDB.create_splice_sites``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def create_gffutils_splice_sites(source: str | Path, *, exon_featuretype: str='exon', grandparent_featuretype: str | None='gene', parent_featuretype: str | None=None, merge_attributes: bool=True, numeric_sort: bool=False, output_path: str | Path | None=None, max_rows: int=_DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Derive 5' and 3' splice-site features from the exon boundaries of each transcript (``FeatureDB.create_splice_sites``).

Every exon must carry an ``ID`` attribute, as gffutils builds the new IDs from them.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.gffutils_children_bp`

**作用：** Total base pairs of the children of one type of a feature, for example exonic bp of a gene (``FeatureDB.children_bp``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def gffutils_children_bp(source: str | Path, feature_id: str, *, child_featuretype: str='exon', merge: bool=False, ignore_strand: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Total base pairs of the children of one type of a feature, for example exonic bp of a gene (``FeatureDB.children_bp``).

With ``merge=True`` overlapping children are merged first (same sequence, strand and type).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.extract_gffutils_feature_sequence`

**作用：** Extract the genome sequence of a feature from a FASTA file (``Feature.sequence``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def extract_gffutils_feature_sequence(source: str | Path, fasta: str | Path, feature_id: str, *, use_strand: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Extract the genome sequence of a feature from a FASTA file (``Feature.sequence``).

With ``use_strand=True`` a minus-strand feature is reverse-complemented, as gffutils does.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.write_gffutils_gff`

**作用：** Write features as GFF records with ``GFFWriter`` (``write_rec`` / ``write_recs``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_gffutils_gff(source: str | Path, output_path: str | Path, *, featuretype: str | None=None, canonical_gene_order: bool=False, gene_featuretype: str='gene', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write features as GFF records with ``GFFWriter`` (``write_rec`` / ``write_recs``).

``canonical_gene_order=True`` instead writes each gene followed by its mRNAs (longest first), their
exons and the remaining children (``write_gene_recs``); ``featuretype`` is then ignored.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.remove_gffutils_duplicates`

**作用：** Remove duplicate GFF records by merging records that share an ID (``gffutils-cli rmdups``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def remove_gffutils_duplicates(source: str | Path, output_path: str | Path, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Remove duplicate GFF records by merging records that share an ID (``gffutils-cli rmdups``).

The database is built with ``merge_strategy="merge"`` and every record is written once to ``output_path``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.gffutils_p37_functions.plot_gffutils_gene_model`

**作用：** Draw a gene model (transcript line, UTRs, CDS) with ``gffutils.contrib.plotting.Gene`` and save the image.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def plot_gffutils_gene_model(source: str | Path, gene_id: str, output_path: str | Path, *, transcripts: Sequence[str]=('mRNA',), utrs: Sequence[str]=("3'UTR", "5'UTR"), cds: Sequence[str]=('CDS',), verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Draw a gene model (transcript line, UTRs, CDS) with ``gffutils.contrib.plotting.Gene`` and save the image.

``utrs`` depends on the annotation: GTF uses ``3'UTR``/``5'UTR``, FlyBase GFF ``three_prime_UTR``/``five_prime_UTR``.
Needs pybedtools and matplotlib. The result holds the image path and the drawn transcripts.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.predictions.seqme_p37_functions import seqme_property_mean
r = seqme_property_mean(['ACGT', 'GGGG', 'AAAA'], 'gc_content', verbose=False)
print(r['value'], r['generated_count'])
```

**示例输出：**

```text
0.5 3
```
