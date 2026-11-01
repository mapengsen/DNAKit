"""Quality, diversity and novelty metrics for generated DNA sets, embedding models and projections through seqme (pack P37).

seqme scores a set of generated sequences against a reference (training) set. Its metrics take an embedding
function or a property predictor; the functions here build both from plain arguments:

- ``kmer_size`` selects a k-mer frequency embedder (all ``4**k`` upper-case DNA k-mers, the seqme
  ``KmerFrequencyEmbedding``). Sequences use the letters A, C, G, T and N (N is ignored by the embedder); k-mers
  are case sensitive, so lower-case input is rejected instead of silently scoring zero.
- ``property_name`` selects a built-in sequence property: ``gc_content`` (fraction of G and C), ``length``,
  ``max_homopolymer_length`` or ``cpg_fraction`` (CG dinucleotides per dinucleotide position).

Every metric function returns ``value`` and ``deviation`` (None when seqme gives none) and the optimisation
direction. seqme is imported only when a function runs; a missing package raises ``BackendUnavailableError`` with
the ``dnakit[dna-generation-metrics]`` hint. Computation runs on CPU. Projection and plot functions take embedding matrices
(lists of equal-length number lists) as seqme does.
"""

from __future__ import annotations

import itertools
import re
from collections.abc import Callable, Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "dna-generation-metrics"
_ALPHABET = frozenset("ACGTN")
_MAX_SEQUENCES = 100_000
_MAX_LENGTH = 10_000_000
_PROPERTIES = ("gc_content", "length", "max_homopolymer_length", "cpg_fraction")
_PROJECTIONS = ("pca", "tsne", "umap")
_GENALM_EMBEDDING = (
    "bert_base_t2t", "bert_base_t2t_lastln_t2t", "bert_base_t2t_multi", "bert_large_t2t", "bigbird_base_t2t",
)  # fmt: skip
_GENALM_CLASSIFIERS = (
    "bert_base_t2t_promoters", "bert_large_t2t_promoters", "bert_large_t2t_promoters2",
    "bert_base_t2t_splice_site", "bert_large_t2t_splice_site",
)  # fmt: skip
_IMAGE_SUFFIXES = (".png", ".svg", ".pdf")


def _provenance() -> dict[str, Any]:
    try:
        installed: str | None = metadata.version("seqme")
    except metadata.PackageNotFoundError:
        installed = None
    return {"package": "seqme", "package_version": installed}


def _seqme(module: str = "seqme") -> Any:
    return require(module, _EXTRA)


def _integer(value: object, name: str, minimum: int, maximum: int = 10**9) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    if positive and value <= 0:
        raise ConfigurationError(f"{name} must be positive.")
    return float(value)


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _sequences(value: object, name: str, *, minimum: int = 1) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) < minimum:
        raise ConfigurationError(f"{name} must be a list of at least {minimum} DNA sequence(s).")
    if len(value) > _MAX_SEQUENCES:
        raise ConfigurationError(f"{name} may hold at most {_MAX_SEQUENCES} sequences.")
    for item in value:
        if not isinstance(item, str) or not item or len(item) > _MAX_LENGTH or set(item) - _ALPHABET:
            raise ConfigurationError(f"every sequence in {name} must be a non-empty upper-case string of A, C, G, T and N.")
    return list(value)


def _matrix(value: object, name: str, *, minimum_rows: int = 2) -> Any:
    np = require("numpy", _EXTRA)
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a matrix of numbers.") from exc
    if array.ndim != 2 or array.shape[0] < minimum_rows or array.shape[1] < 1:
        raise ConfigurationError(f"{name} must be a matrix with at least {minimum_rows} rows and one column.")
    if not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must hold finite numbers.")
    return array


def _embedder(kmer_size: int) -> Callable[[list[str]], Any]:
    size = _integer(kmer_size, "kmer_size", 1, 6)
    models = _seqme("seqme.models")
    kmers = ["".join(item) for item in itertools.product("ACGT", repeat=size)]
    return models.KmerFrequencyEmbedding(kmers)


def _property(name: str) -> Callable[[list[str]], Any]:
    _choice(name, _PROPERTIES, "property_name")
    np = require("numpy", _EXTRA)

    def predictor(sequences: list[str]) -> Any:
        values = []
        for sequence in sequences:
            if name == "gc_content":
                values.append((sequence.count("G") + sequence.count("C")) / len(sequence))
            elif name == "length":
                values.append(len(sequence))
            elif name == "max_homopolymer_length":
                values.append(max(len(m.group(0)) for m in re.finditer(r"(.)\1*", sequence)))
            else:
                values.append(sequence.count("CG") / max(len(sequence) - 1, 1))
        return np.asarray(values, dtype=float)

    return predictor


def _properties(names: Sequence[str], minimum: int) -> list[Callable[[list[str]], Any]]:
    if isinstance(names, (str, bytes)) or len(names) < minimum:
        raise ConfigurationError(f"property_names must list at least {minimum} property name(s) from {list(_PROPERTIES)}.")
    return [_property(item) for item in names]


def _result(operation: str, metric: Any, outcome: Any, generated: Sequence[str], **extra: object) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "metric": metric.name,
        "objective": metric.objective,
        "value": float(outcome.value),
        "deviation": None if outcome.deviation is None else float(outcome.deviation),
        "generated_count": len(generated),
        **extra,
        "provenance": _provenance(),
    }


def _run(operation: str, build: Callable[[], Any], generated: Sequence[str], **extra: object) -> dict[str, Any]:
    try:
        metric = build()
        outcome = metric(list(generated))
    except ConfigurationError:
        raise
    except (ValueError, TypeError, IndexError, KeyError) as exc:
        raise ConfigurationError(f"seqme could not compute the metric: {exc}") from exc
    return _result(operation, metric, outcome, generated, **extra)


def seqme_authenticity(generated: Sequence[str], train_set: Sequence[str], *, kmer_size: int = 3, verbose: bool = True) -> dict[str, Any]:
    """Authenticity (AuthPct): the fraction of generated sequences whose nearest training neighbour is closer to another training sequence than to the generated one (Alaa et al. 2022).

    Higher is better; low values mean the generator copies its training data. k-mer embeddings of ``kmer_size`` are used.
    """

    seqs, train = _sequences(generated, "generated"), _sequences(train_set, "train_set", minimum=2)
    metrics = _seqme("seqme.metrics")
    embedder = _embedder(kmer_size)
    return _run("seqme_authenticity", lambda: metrics.AuthPct(train, embedder), seqs, kmer_size=kmer_size)


def _manifold_metric(operation: str, name: str, generated: Sequence[str], reference: Sequence[str], n_neighbors: int, kmer_size: int, strict: bool) -> dict[str, Any]:
    seqs = _sequences(generated, "generated")
    ref = _sequences(reference, "reference", minimum=2)
    neighbors = _integer(n_neighbors, "n_neighbors", 1, 1000)
    strict = _flag(strict, "strict")
    metrics = _seqme("seqme.metrics")
    embedder = _embedder(kmer_size)
    cls = getattr(metrics, name)
    return _run(operation, lambda: cls(neighbors, ref, embedder, strict=strict), seqs, n_neighbors=neighbors, kmer_size=kmer_size)


def seqme_clipped_density(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int = 5, kmer_size: int = 3, strict: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Clipped density, an outlier-robust fidelity score of generated sequences against a reference manifold built from k-nearest-neighbour radii (Salvy et al. 2024). Range 0 to 1, higher is better.

    ``strict=True`` requires equally many generated and reference sequences.
    """

    return _manifold_metric("seqme_clipped_density", "ClippedDensity", generated, reference, n_neighbors, kmer_size, strict)


def seqme_clipped_coverage(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int = 5, kmer_size: int = 3, strict: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Clipped coverage, an outlier-robust diversity score: how much of the reference manifold the generated sequences cover. Range 0 to 1, higher is better."""

    return _manifold_metric("seqme_clipped_coverage", "ClippedCoverage", generated, reference, n_neighbors, kmer_size, strict)


def seqme_precision(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int = 5, kmer_size: int = 3, strict: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Improved precision (Kynkaanniemi et al. 2019): the share of generated sequences inside the k-NN manifold of the reference set (fidelity). Higher is better."""

    return _manifold_metric("seqme_precision", "Precision", generated, reference, n_neighbors, kmer_size, strict)


def seqme_recall(generated: Sequence[str], reference: Sequence[str], *, n_neighbors: int = 5, kmer_size: int = 3, strict: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Improved recall (Kynkaanniemi et al. 2019): the share of reference sequences inside the k-NN manifold of the generated set (coverage). Higher is better."""

    return _manifold_metric("seqme_recall", "Recall", generated, reference, n_neighbors, kmer_size, strict)


def seqme_conformity_score(
    generated: Sequence[str],
    reference: Sequence[str],
    *,
    property_names: Sequence[str] = ("gc_content", "length"),
    n_splits: int = 5,
    kde_bandwidth: str | float = "silverman",
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distributional conformity score: how well the property distribution of the generated set fits a kernel density estimate of the reference set, with cross-validation (Gaussian KDE).

    ``property_names`` are built-in sequence properties (gc_content, length, max_homopolymer_length, cpg_fraction).
    """

    seqs, ref = _sequences(generated, "generated"), _sequences(reference, "reference", minimum=2)
    predictors = _properties(property_names, 1)
    metrics = _seqme("seqme.metrics")
    bandwidth = kde_bandwidth if isinstance(kde_bandwidth, str) else _number(kde_bandwidth, "kde_bandwidth", positive=True)
    return _run(
        "seqme_conformity_score",
        lambda: metrics.ConformityScore(ref, predictors, n_splits=_integer(n_splits, "n_splits", 2, 100), kde_bandwidth=bandwidth, seed=seed),
        seqs, property_names=list(property_names),
    )  # fmt: skip


def seqme_fkea_diversity(generated: Sequence[str], *, bandwidth: float = 1.0, alpha: float = 2, n_random_fourier_features: int | None = 2048, kmer_size: int = 3, seed: int = 0, verbose: bool = True) -> dict[str, Any]:
    """FKEA kernel-entropy diversity (Vendi / RKE score): the effective number of distinct modes of the generated k-mer embeddings under a Gaussian kernel (Friedman and Dieng 2023; Ospanov et al. 2024).

    ``alpha=2`` gives the RKE score, other values the Vendi-alpha score. ``n_random_fourier_features=None`` uses the exact kernel.
    """

    seqs = _sequences(generated, "generated")
    metrics = _seqme("seqme.metrics")
    embedder = _embedder(kmer_size)
    fourier = None if n_random_fourier_features is None else _integer(n_random_fourier_features, "n_random_fourier_features", 1, 10**6)
    np = require("numpy", _EXTRA)

    def single_precision(items: list[str]) -> Any:
        # FKEA multiplies the embeddings with float32 random features, so float64 k-mer frequencies must be cast.
        return np.asarray(embedder(items), dtype=np.float32)

    return _run(
        "seqme_fkea_diversity",
        lambda: metrics.FKEA(single_precision, _number(bandwidth, "bandwidth", positive=True), alpha=_number(alpha, "alpha", positive=True), n_random_fourier_features=fourier, seed=seed),
        seqs, kmer_size=kmer_size,
    )  # fmt: skip


def seqme_hit_rate(generated: Sequence[str], property_name: str, *, min_value: float | None = None, max_value: float | None = None, verbose: bool = True) -> dict[str, Any]:
    """Hit rate: the proportion of generated sequences whose property lies in ``[min_value, max_value]`` (a bound left as None is open), with its standard error."""

    seqs = _sequences(generated, "generated")
    predictor = _property(property_name)
    if min_value is None and max_value is None:
        raise ConfigurationError("give min_value, max_value or both.")
    low = None if min_value is None else _number(min_value, "min_value")
    high = None if max_value is None else _number(max_value, "max_value")
    np = require("numpy", _EXTRA)
    metrics = _seqme("seqme.metrics")

    def condition(items: list[str]) -> Any:
        values = predictor(items)
        keep = np.ones(len(values), dtype=bool)
        if low is not None:
            keep &= values >= low
        if high is not None:
            keep &= values <= high
        return keep

    return _run("seqme_hit_rate", lambda: metrics.HitRate(condition, name=f"Hit-rate {property_name}"), seqs, property_name=property_name, min_value=low, max_value=high)


def seqme_hypervolume(
    generated: Sequence[str],
    property_names: Sequence[str],
    *,
    method: str = "hvi",
    nadir: Sequence[float] | None = None,
    ideal: Sequence[float] | None = None,
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multi-objective hypervolume of the property values of the generated set (all properties maximised): the hypervolume indicator (``method='hvi'``) or the convex-hull volume (``'convex-hull'``).

    ``nadir`` (worst value per property, default zeros) and ``ideal`` (best value, used for scaling to [0, 1]) have one number per property.
    """

    seqs = _sequences(generated, "generated")
    predictors = _properties(property_names, 2)
    _choice(method, ("hvi", "convex-hull"), "method")
    np = require("numpy", _EXTRA)
    metrics = _seqme("seqme.metrics")
    worst = None if nadir is None else np.asarray([_number(x, "nadir") for x in nadir])
    best = None if ideal is None else np.asarray([_number(x, "ideal") for x in ideal])
    return _run(
        "seqme_hypervolume",
        lambda: metrics.Hypervolume(predictors, method=method, nadir=worst, ideal=best, strict=_flag(strict, "strict")),
        seqs, property_names=list(property_names), method=method,
    )  # fmt: skip


def seqme_property_mean(generated: Sequence[str], property_name: str, *, objective: str = "maximize", deviation: str = "se", estimate: str = "unbiased", verbose: bool = True) -> dict[str, Any]:
    """Mean of a sequence property over the generated set with its deviation (the seqme ``ID`` metric); ``deviation`` is ``std``, ``se`` or ``var``, ``estimate`` ``biased`` or ``unbiased``."""

    seqs = _sequences(generated, "generated")
    predictor = _property(property_name)
    _choice(objective, ("maximize", "minimize"), "objective")
    _choice(deviation, ("std", "se", "var"), "deviation")
    _choice(estimate, ("biased", "unbiased"), "estimate")
    metrics = _seqme("seqme.metrics")
    return _run("seqme_property_mean", lambda: metrics.ID(predictor, property_name, objective, deviation=deviation, estimate=estimate), seqs, property_name=property_name)


def seqme_property_kl_divergence(generated: Sequence[str], reference: Sequence[str], property_name: str, *, n_draws: int = 10_000, kde_bandwidth: str | float = "silverman", seed: int = 0, verbose: bool = True) -> dict[str, Any]:
    """KL divergence between the property distribution of the generated set and the reference set, estimated with Gaussian KDEs and ``n_draws`` Monte Carlo samples. Lower is better."""

    seqs, ref = _sequences(generated, "generated", minimum=2), _sequences(reference, "reference", minimum=2)
    predictor = _property(property_name)
    metrics = _seqme("seqme.metrics")
    bandwidth = kde_bandwidth if isinstance(kde_bandwidth, str) else _number(kde_bandwidth, "kde_bandwidth", positive=True)
    return _run(
        "seqme_property_kl_divergence",
        lambda: metrics.KLDivergence(ref, predictor, n_draws=_integer(n_draws, "n_draws", 1, 10**7), kde_bandwidth=bandwidth, seed=seed),
        seqs, property_name=property_name,
    )  # fmt: skip


def seqme_mmd(generated: Sequence[str], reference: Sequence[str], *, estimate: str = "biased", sigma: float = 10, scale: float = 1000, kmer_size: int = 3, verbose: bool = True) -> dict[str, Any]:
    """Maximum mean discrepancy with a Gaussian RBF kernel between the k-mer embeddings of the generated and reference sets, multiplied by ``scale``. Lower is better."""

    seqs, ref = _sequences(generated, "generated"), _sequences(reference, "reference", minimum=2)
    _choice(estimate, ("biased", "unbiased"), "estimate")
    metrics = _seqme("seqme.metrics")
    embedder = _embedder(kmer_size)
    return _run("seqme_mmd", lambda: metrics.MMD(ref, embedder, estimate=estimate, sigma=_number(sigma, "sigma", positive=True), scale=_number(scale, "scale", positive=True)), seqs, kmer_size=kmer_size)


def seqme_kid(generated: Sequence[str], reference: Sequence[str], *, estimate: str = "biased", degree: int = 3, coef0: float = 1.0, kmer_size: int = 3, verbose: bool = True) -> dict[str, Any]:
    """Kernel inception distance: MMD with a polynomial kernel (``degree``, ``coef0``) between k-mer embeddings (Binkowski et al. 2018). Lower is better."""

    seqs, ref = _sequences(generated, "generated"), _sequences(reference, "reference", minimum=2)
    _choice(estimate, ("biased", "unbiased"), "estimate")
    metrics = _seqme("seqme.metrics")
    embedder = _embedder(kmer_size)
    return _run("seqme_kid", lambda: metrics.KID(ref, embedder, estimate=estimate, degree=_integer(degree, "degree", 1, 20), coef0=_number(coef0, "coef0")), seqs, kmer_size=kmer_size)


def seqme_threshold_fraction(generated: Sequence[str], property_name: str, *, min_value: float | None = None, max_value: float | None = None, verbose: bool = True) -> dict[str, Any]:
    """Fraction of generated sequences whose property passes a threshold: at least ``min_value`` and/or at most ``max_value`` (give one or both)."""

    seqs = _sequences(generated, "generated")
    predictor = _property(property_name)
    low = None if min_value is None else _number(min_value, "min_value")
    high = None if max_value is None else _number(max_value, "max_value")
    if low is None and high is None:
        raise ConfigurationError("give min_value, max_value or both.")
    metrics = _seqme("seqme.metrics")
    return _run("seqme_threshold_fraction", lambda: metrics.Threshold(predictor, property_name, min_value=low, max_value=high), seqs, property_name=property_name, min_value=low, max_value=high)


# Embedding models and projections -----------------------------------------------------------------------


def seqme_embed_rna_fm(sequences: Sequence[str], *, model_name: str = "mRNA", layer: int = 12, batch_size: int = 32, verbose: bool = True) -> dict[str, Any]:
    """Embed RNA/mRNA sequences with RNA-FM (``mRNA`` codon-aligned CDS model, 1280 dimensions, or ``ncRNA``, 640 dimensions) by averaging token embeddings of ``layer`` (last layer is 12).

    Needs the ``rna-fm`` package and downloads the model weights on first use. Sequences are RNA (U) or DNA (T) letters as RNA-FM expects.
    """

    _choice(model_name, ("mRNA", "ncRNA"), "model_name")
    seqs = _sequences_loose(sequences)
    models = _seqme("seqme.models")
    try:
        model = models.RNAFM(model_name=model_name, batch_size=_integer(batch_size, "batch_size", 1, 4096), device="cpu")
        matrix = model.embed(seqs, layer=_integer(layer, "layer", 0, 12))
    except (ValueError, TypeError, RuntimeError) as exc:
        raise ConfigurationError(f"RNA-FM could not embed the sequences: {exc}") from exc
    return _embedding_result("seqme_embed_rna_fm", matrix, model_name=model_name, layer=layer)


def _sequences_loose(value: object) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value or len(value) > _MAX_SEQUENCES:
        raise ConfigurationError("sequences must be a non-empty list of sequence strings.")
    if not all(isinstance(item, str) and item and re.fullmatch(r"[ACGTUNacgtun]+", item) for item in value):
        raise ConfigurationError("every sequence must be a non-empty string of nucleotide letters.")
    return list(value)


def _embedding_result(operation: str, matrix: Any, **extra: object) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "embeddings": [[float(x) for x in row] for row in matrix],
        **extra,
        "provenance": _provenance(),
    }


def seqme_embed_genalm(sequences: Sequence[str], *, checkpoint: str = "bert_base_t2t", layer: int = -1, batch_size: int = 16, cache_dir: str | None = None, verbose: bool = True) -> dict[str, Any]:
    """Embed DNA sequences with a GENA-LM language-model checkpoint (mean of the token embeddings of ``layer``, default the last).

    ``checkpoint`` is one of bert_base_t2t, bert_base_t2t_lastln_t2t, bert_base_t2t_multi, bert_large_t2t or bigbird_base_t2t.
    Needs ``transformers`` and ``einops`` and downloads the model from Hugging Face on first use (network).
    """

    _choice(checkpoint, _GENALM_EMBEDDING, "checkpoint")
    seqs = _sequences(sequences, "sequences")
    models = _seqme("seqme.models")
    try:
        model = models.GENALM(getattr(models.GENALMCheckpoint, checkpoint), device="cpu", batch_size=_integer(batch_size, "batch_size", 1, 4096), cache_dir=cache_dir)
        matrix = model.embed(seqs, layer=_integer(layer, "layer", -64, 64))
    except (ValueError, TypeError, RuntimeError, OSError) as exc:
        raise ConfigurationError(f"GENA-LM could not embed the sequences: {exc}") from exc
    return _embedding_result("seqme_embed_genalm", matrix, checkpoint=checkpoint, layer=layer)


def seqme_classify_genalm(sequences: Sequence[str], *, checkpoint: str = "bert_base_t2t_promoters", batch_size: int = 16, cache_dir: str | None = None, verbose: bool = True) -> dict[str, Any]:
    """Classify DNA sequences with a fine-tuned GENA-LM model and return class probabilities.

    ``checkpoint`` is a promoter model (bert_base_t2t_promoters, bert_large_t2t_promoters at 300 bp, bert_large_t2t_promoters2 at 2000 bp;
    classes: no promoter, promoter) or a splice-site model (bert_base_t2t_splice_site, bert_large_t2t_splice_site at 15000 bp; classes:
    neither, splice acceptor, splice donor). Needs ``transformers`` and downloads the model from Hugging Face on first use (network).
    """

    _choice(checkpoint, _GENALM_CLASSIFIERS, "checkpoint")
    seqs = _sequences(sequences, "sequences")
    models = _seqme("seqme.models")
    try:
        model = models.GENALM(getattr(models.GENALMCheckpoint, checkpoint), device="cpu", batch_size=_integer(batch_size, "batch_size", 1, 4096), cache_dir=cache_dir)
        probabilities = model.classify(seqs)
    except (ValueError, TypeError, RuntimeError, OSError) as exc:
        raise ConfigurationError(f"GENA-LM could not classify the sequences: {exc}") from exc
    classes = ["no_promoter", "promoter"] if "promoters" in checkpoint else ["neither", "splice_acceptor", "splice_donor"]
    return {
        "status": "completed",
        "operation": "seqme_classify_genalm",
        "checkpoint": checkpoint,
        "classes": classes,
        "probabilities": [[float(x) for x in row] for row in probabilities],
        "provenance": _provenance(),
    }


def seqme_project_embeddings(
    embeddings: Sequence[Sequence[float]],
    *,
    method: str = "pca",
    seed: int | None = 0,
    reference_embeddings: Sequence[Sequence[float]] | None = None,
    n_components: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project embeddings to low dimension by PCA, t-SNE or UMAP (``seqme.utils.pca``, ``tsne``, ``umap``).

    ``method`` is ``pca``, ``tsne`` or ``umap``; each returns 2-D points (t-SNE with PCA initialisation, perplexity and
    learning rate as in seqme). With ``method='pca'`` and ``reference_embeddings``, PCA is fitted on the reference
    rows and applied to ``embeddings`` with ``n_components`` components (``seqme.models.PCA``); the explained variance ratio is returned.
    """

    _choice(method, _PROJECTIONS, "method")
    matrix = _matrix(embeddings, "embeddings")
    utils = _seqme("seqme.utils")
    np = require("numpy", _EXTRA)
    extra: dict[str, Any] = {}
    try:
        if reference_embeddings is not None:
            if method != "pca":
                raise ConfigurationError("reference_embeddings is only used with method='pca'.")
            reference = _matrix(reference_embeddings, "reference_embeddings")
            models = _seqme("seqme.models")
            fitted = models.PCA(lambda items: np.asarray(items, dtype=float), reference, _integer(n_components, "n_components", 1, 1000))
            projected = fitted(matrix)
            extra["explained_variance_ratio"] = [float(x) for x in fitted.variance_explained]
        else:
            projected = getattr(utils, "pca" if method == "pca" else method)(matrix, seed=seed)
    except ConfigurationError:
        raise
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"seqme could not project the embeddings: {exc}") from exc
    return {
        "status": "completed",
        "operation": "seqme_project_embeddings",
        "method": method,
        "shape": [int(projected.shape[0]), int(projected.shape[1])],
        "projection": [[float(x) for x in row] for row in projected],
        **extra,
        "provenance": _provenance(),
    }


def seqme_spearman_alignment_score(embeddings_a: Sequence[Sequence[float]], embeddings_b: Sequence[Sequence[float]], *, verbose: bool = True) -> dict[str, Any]:
    """Spearman correlation between the pairwise-distance matrices of two embeddings of the same sequences (alignment of two embedding spaces, -1 to 1)."""

    a, b = _matrix(embeddings_a, "embeddings_a", minimum_rows=3), _matrix(embeddings_b, "embeddings_b", minimum_rows=3)
    if a.shape[0] != b.shape[0]:
        raise ConfigurationError("both embeddings must describe the same number of sequences.")
    utils = _seqme("seqme.utils")
    return {
        "status": "completed",
        "operation": "seqme_spearman_alignment_score",
        "score": float(utils.spearman_alignment_score(a, b)),
        "sequence_count": int(a.shape[0]),
        "provenance": _provenance(),
    }


def _labels(labels: object, rows: int) -> Any:
    np = require("numpy", _EXTRA)
    if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence) or len(labels) != rows:
        raise ConfigurationError("labels must be a list with one label per embedding row.")
    return np.asarray(list(labels))


def seqme_knn_alignment_score(embeddings: Sequence[Sequence[float]], labels: Sequence[Any], *, n_neighbors: int = 5, verbose: bool = True) -> dict[str, Any]:
    """k-NN feature alignment score of an embedding model: the share of each sequence's ``n_neighbors`` nearest neighbours that carry the same label (0 to 1)."""

    matrix = _matrix(embeddings, "embeddings")
    k = _integer(n_neighbors, "n_neighbors", 1, 1000)
    if k >= matrix.shape[0]:
        raise ConfigurationError("n_neighbors must be smaller than the number of sequences.")
    utils = _seqme("seqme.utils")
    return {
        "status": "completed",
        "operation": "seqme_knn_alignment_score",
        "n_neighbors": k,
        "score": float(utils.knn_alignment_score(matrix, _labels(labels, matrix.shape[0]), k)),
        "provenance": _provenance(),
    }


def _new_image(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.suffix.lower() not in _IMAGE_SUFFIXES:
        raise ConfigurationError(f"output_path must end with one of {list(_IMAGE_SUFFIXES)}.")
    if path.exists():
        raise ConfigurationError(f"output_path already exists: {path.name}; choose a new file.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of output_path does not exist: {path.parent}.")
    return path


def seqme_plot_embeddings(
    projections: Mapping[str, Sequence[Sequence[float]]],
    output_path: str | Path,
    *,
    title: str | None = None,
    point_size: float = 20,
    alpha: float = 0.6,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scatter-plot 2-D embedding projections, one coloured group per key (``seqme.utils.plot_embeddings``), and save the image (png, svg or pdf).

    ``projections`` maps a group name to its list of ``[x, y]`` points (for example from ``seqme_project_embeddings``).
    """

    target = _new_image(output_path)
    if not isinstance(projections, Mapping) or not projections:
        raise ConfigurationError("projections must be a non-empty mapping from group name to 2-D points.")
    groups = {str(name): _matrix(points, f"projections[{name}]", minimum_rows=1) for name, points in projections.items()}
    for name, array in groups.items():
        if array.shape[1] != 2:
            raise ConfigurationError(f"projections[{name}] must have two columns.")
    utils = _seqme("seqme.utils")
    plt = require("matplotlib.pyplot", _EXTRA)
    figure, axes = plt.subplots(figsize=(5, 3.5))
    try:
        utils.plot_embeddings(list(groups.values()), values=list(groups), title=title, point_size=_number(point_size, "point_size", positive=True), alpha=_number(alpha, "alpha", positive=True), ax=axes)
        figure.savefig(target, dpi=150, bbox_inches="tight")
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "seqme_plot_embeddings",
        "image_path": str(target),
        "groups": {name: int(array.shape[0]) for name, array in groups.items()},
    }


def seqme_plot_knn_alignment_score(
    embeddings: Sequence[Sequence[float]],
    labels: Sequence[Any],
    output_path: str | Path,
    *,
    n_neighbors: Sequence[int] = (1, 3, 5, 10),
    label: str = "embedding",
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the k-NN feature alignment score against the number of neighbours (``seqme.utils.plot_knn_alignment_score``) and save the image; the scores are returned too."""

    target = _new_image(output_path)
    matrix = _matrix(embeddings, "embeddings")
    if isinstance(n_neighbors, (str, bytes)) or not n_neighbors:
        raise ConfigurationError("n_neighbors must be a non-empty list of integers.")
    ks = [_integer(k, "n_neighbors", 1, 1000) for k in n_neighbors]
    if max(ks) >= matrix.shape[0]:
        raise ConfigurationError("every n_neighbors value must be smaller than the number of sequences.")
    array = _labels(labels, matrix.shape[0])
    utils = _seqme("seqme.utils")
    plt = require("matplotlib.pyplot", _EXTRA)
    figure, axes = plt.subplots(figsize=(5, 3.5))
    try:
        utils.plot_knn_alignment_score(matrix, array, ks, label=label, legend_loc="best", ax=axes)
        figure.savefig(target, dpi=150, bbox_inches="tight")
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "seqme_plot_knn_alignment_score",
        "image_path": str(target),
        "n_neighbors": ks,
        "scores": [float(utils.knn_alignment_score(matrix, array, k)) for k in ks],
    }


__all__ = [
    "seqme_authenticity",
    "seqme_classify_genalm",
    "seqme_clipped_coverage",
    "seqme_clipped_density",
    "seqme_conformity_score",
    "seqme_embed_genalm",
    "seqme_embed_rna_fm",
    "seqme_fkea_diversity",
    "seqme_hit_rate",
    "seqme_hypervolume",
    "seqme_kid",
    "seqme_knn_alignment_score",
    "seqme_mmd",
    "seqme_plot_embeddings",
    "seqme_plot_knn_alignment_score",
    "seqme_precision",
    "seqme_project_embeddings",
    "seqme_property_kl_divergence",
    "seqme_property_mean",
    "seqme_recall",
    "seqme_spearman_alignment_score",
    "seqme_threshold_fraction",
]

_register_public_verbose(globals(), tuple(__all__))
