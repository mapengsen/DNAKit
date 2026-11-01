"""gReLU operations as typed DNAKit facades: genome and interval utilities, model-zoo and W&B look-ups, tomtom motif
comparison, and plots of predictions and attributions.

The user installs grelu (``pip install grelu``) with its genome tools (genomepy) for genome downloads. Each function
calls gReLU's own function for the operation with the caller's arguments; nothing is re-implemented. Missing optional
packages raise BackendUnavailableError with an install hint. Network-backed functions (genome downloads, the
Hugging Face and Weights & Biases model zoos, JASPAR) reach the network only when called; the tests replace them with
fakes. Plot functions save an image file and return its path together with the plotted data.
"""

from __future__ import annotations

import functools
import importlib
from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_HINT = "Install gReLU in the same environment as DNAKit: python -m pip install grelu."
_GENOME_HINT = "Install genomepy (python -m pip install genomepy) to download and read genomes."
_HUB_HINT = "Install huggingface_hub (python -m pip install huggingface_hub) to use the gReLU model zoo."
_WANDB_HINT = "Install wandb (python -m pip install wandb) to use the gReLU Weights & Biases model zoo."
_PLOT_HINT = "Install matplotlib (python -m pip install matplotlib) to draw plots."
_MEMELITE_HINT = "Install memelite (python -m pip install memelite) for the tomtom motif comparison."
_IMAGE_SUFFIXES = (".png", ".pdf", ".svg")
_WANDB_HOST = "https://api.wandb.ai"


def _upstream(module: str, hint: str = _HINT) -> Any:
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(f"{module} could not be imported: {exc}", hint=hint) from exc


@functools.lru_cache(maxsize=1)
def _version() -> str:
    try:
        return metadata.version("grelu")
    except metadata.PackageNotFoundError:
        return "unknown"


def _text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def _existing(path: str | Path, name: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return str(resolved)


def _new_image(path: str | Path, name: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if resolved.suffix.lower() not in _IMAGE_SUFFIXES:
        raise ConfigurationError(f"{name} must end with .png, .pdf or .svg.")
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return str(resolved)


def _call(module: str, function: str, hint: str, *args: Any, **kwargs: Any) -> Any:
    """Call gReLU's function; upstream failures are raised as BackendExecutionError with the original message."""
    target = getattr(_upstream(module, hint), function)
    try:
        return target(*args, **kwargs)
    except Exception as exc:
        raise BackendExecutionError(f"gReLU {function} failed: {exc}") from exc


def grelu_check_chrom_ends(
    data: Any,
    genome: str | None = None,
    *,
    verbose: bool = True,
) -> None:
    """Check that genomic intervals do not extend past the chromosome ends of a genome (gReLU check_chrom_ends; genome sizes via genomepy)."""
    return _call("grelu.data.preprocess", "check_chrom_ends", _GENOME_HINT, data, genome=genome)


def grelu_filter_chrom_ends(
    data: Any,
    genome: str | None = None,
    *,
    pad: int = 0,
    inplace: bool = False,
    verbose: bool = True,
) -> Any:
    """Drop intervals that extend beyond the chromosome ends, optionally after padding them (gReLU filter_chrom_ends)."""
    return _call("grelu.data.preprocess", "filter_chrom_ends", _GENOME_HINT, data, genome=genome, pad=pad,
                 inplace=inplace)


def grelu_get_gc_matched_intervals(
    intervals: Any,
    genome: str,
    *,
    binwidth: float = 0.1,
    chroms: str = "autosomes",
    blacklist: str | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> Any:
    """Sample background intervals whose GC content matches the given intervals (gReLU get_gc_matched_intervals; genome via genomepy)."""
    _text(genome, "genome")
    return _call("grelu.data.preprocess", "get_gc_matched_intervals", _GENOME_HINT, intervals, genome,
                 binwidth=binwidth, chroms=chroms, blacklist=blacklist, seed=seed)


def grelu_make_insertion_bigwig(
    frag_file: str,
    genome: str,
    *,
    out_prefix: str | None = None,
    plus_shift: int = 0,
    minus_shift: int = 0,
    chroms: Sequence[str] | str | None = None,
    tmp_dir: str = "./",
    out_dir: str = "./",
    verbose: bool = True,
) -> str:
    """Convert a fragment file of Tn5 insertions into a strand-specific insertion bigWig for a genome (gReLU make_insertion_bigwig)."""
    _text(genome, "genome")
    source = _existing(frag_file, "frag_file")
    return _call("grelu.data.preprocess", "make_insertion_bigwig", _GENOME_HINT, source, genome,
                 out_prefix=out_prefix, plus_shift=plus_shift, minus_shift=minus_shift, chroms=chroms,
                 tmp_dir=tmp_dir, out_dir=out_dir)


def grelu_read_sizes(genome: str = "hg38", *, verbose: bool = True) -> Any:
    """Read the chromosome sizes of a genome as a table (gReLU read_sizes; the genome is installed through genomepy)."""
    return _call("grelu.io.genome", "read_sizes", _GENOME_HINT, _text(genome, "genome"))


def grelu_get_genome(genome: str, *, verbose: bool = True) -> Any:
    """Install a genome through genomepy and load it as a gReLU/genomepy Genome object (network download on first use)."""
    return _call("grelu.io.genome", "get_genome", _GENOME_HINT, _text(genome, "genome"))


def grelu_read_gtf(genome: str, features: Sequence[str] | str | None = None, *, verbose: bool = True) -> Any:
    """Install a genome annotation through genomepy and read its gene features as a table (gReLU read_gtf)."""
    return _call("grelu.io.genome", "read_gtf", _GENOME_HINT, _text(genome, "genome"), features)


def grelu_get_jaspar(
    release: str = "JASPAR2024",
    collection: str = "CORE",
    *,
    tax_group: str | None = None,
    species: str | None = None,
    verbose: bool = True,
) -> Any:
    """Download motifs from the JASPAR database as position frequency matrices (gReLU get_jaspar; network)."""
    return _call("grelu.io.motifs", "get_jaspar", _HINT, release=release, collection=collection,
                 tax_group=tax_group, species=species)


def grelu_intervals_to_strings(intervals: Any, genome: str, *, verbose: bool = True) -> Any:
    """Extract the DNA sequences of genomic intervals from a genome (gReLU intervals_to_strings; genome via genomepy)."""
    return _call("grelu.sequence.format", "intervals_to_strings", _GENOME_HINT, intervals, _text(genome, "genome"))


def grelu_variant_to_seqs(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    genome: str,
    *,
    seq_len: int = 1,
    verbose: bool = True,
) -> tuple[str, str]:
    """Return the reference and alternate sequences around a single-nucleotide variant of a genome (gReLU variant_to_seqs)."""
    return _call("grelu.variant", "variant_to_seqs", _GENOME_HINT, _text(chrom, "chrom"), pos, ref, alt,
                 _text(genome, "genome"), seq_len=seq_len)


def grelu_check_reference(
    variants: Any,
    *,
    genome: str = "hg38",
    null_string: str = "-",
    verbose: bool = True,
) -> None:
    """Check that the reference alleles of a variant table match the reference genome (gReLU check_reference)."""
    return _call("grelu.variant", "check_reference", _GENOME_HINT, variants, genome=_text(genome, "genome"),
                 null_string=null_string)


def grelu_run_tomtom(
    motifs: Mapping[str, Sequence[Sequence[float]]],
    meme_file: str,
    *,
    verbose: bool = True,
) -> Any:
    """Compare query motifs to a MEME reference database with the tomtom algorithm (gReLU run_tomtom, memelite implementation).

    The motifs are probability matrices (rows = positions, columns = A, C, G, T). ``meme_file`` is a MEME-format file
    or the name of one bundled with gReLU ("hocomoco_v12", "hocomoco_v13", "consensus"). The comparison is pure Python
    through memelite; the MEME Suite tomtom executable is not used.
    """
    np = _upstream("numpy")
    _upstream("memelite", _MEMELITE_HINT)
    if not motifs:
        raise ConfigurationError("motifs must contain at least one probability matrix.")
    arrays = {}
    for name, matrix in motifs.items():
        array = np.asarray(matrix, dtype=float)
        if array.ndim != 2 or array.shape[1] != 4 or not np.isfinite(array).all():
            raise ConfigurationError(f"Motif {name!r} must be a finite matrix with four columns (A, C, G, T).")
        arrays[str(name)] = array
    if meme_file not in ("hocomoco_v12", "hocomoco_v13", "consensus"):
        meme_file = _existing(meme_file, "meme_file")
    return _call("grelu.interpret.motifs", "run_tomtom", _MEMELITE_HINT, arrays, meme_file)


def grelu_list_models(*, verbose: bool = True) -> list[str]:
    """List the model repository IDs of the gReLU model zoo on Hugging Face (network)."""
    return _call("grelu.resources", "list_models", _HUB_HINT)


def grelu_list_datasets(*, verbose: bool = True) -> list[str]:
    """List the dataset repository IDs of the gReLU model zoo on Hugging Face (network)."""
    return _call("grelu.resources", "list_datasets", _HUB_HINT)


def grelu_download_model(repo_id: str, filename: str = "model.ckpt", *, verbose: bool = True) -> str:
    """Download a model checkpoint from the gReLU model zoo on Hugging Face and return its local path (network)."""
    return _call("grelu.resources", "download_model", _HUB_HINT, _text(repo_id, "repo_id"), filename=filename)


def grelu_download_dataset(repo_id: str, filename: str = "data.h5ad", *, verbose: bool = True) -> str:
    """Download a dataset file from the gReLU model zoo on Hugging Face and return its local path (network)."""
    return _call("grelu.resources", "download_dataset", _HUB_HINT, _text(repo_id, "repo_id"), filename=filename)


def grelu_load_model(
    repo_id: str | None = None,
    filename: str = "model.ckpt",
    *,
    device: str | int = "cpu",
    project: str | None = None,
    model_name: str | None = None,
    verbose: bool = True,
) -> Any:
    """Load a gReLU model from the Hugging Face model zoo (repo_id) or from Weights & Biases (project and model_name)."""
    return _call("grelu.resources", "load_model", _HUB_HINT, repo_id=repo_id, filename=filename, device=device,
                 project=project, model_name=model_name)


def grelu_get_model_info(repo_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Read the full model card metadata of a model repository on Hugging Face (network)."""
    return _call("grelu.resources", "get_model_info", _HUB_HINT, _text(repo_id, "repo_id"))


def grelu_get_dataset_info(repo_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Read the full dataset card metadata of a dataset repository on Hugging Face (network)."""
    return _call("grelu.resources", "get_dataset_info", _HUB_HINT, _text(repo_id, "repo_id"))


def grelu_get_datasets_by_model(repo_id: str, *, verbose: bool = True) -> list[str]:
    """List the datasets linked to a model in its model card (network)."""
    return _call("grelu.resources", "get_datasets_by_model", _HUB_HINT, _text(repo_id, "repo_id"))


def grelu_get_base_models(repo_id: str, *, verbose: bool = True) -> list[str]:
    """List the base models a model was fine-tuned from, from its model card (network)."""
    return _call("grelu.resources", "get_base_models", _HUB_HINT, _text(repo_id, "repo_id"))


def grelu_get_models_by_dataset(repo_id: str, *, verbose: bool = True) -> list[str]:
    """List the models of the gReLU collection that were trained on a dataset (network)."""
    return _call("grelu.resources", "get_models_by_dataset", _HUB_HINT, _text(repo_id, "repo_id"))


def grelu_wandb_projects(*, host: str = _WANDB_HOST, verbose: bool = True) -> list[str]:
    """List the projects of the gReLU model zoo on Weights & Biases (network)."""
    return _call("grelu.resources.wandb", "projects", _WANDB_HINT, host=host)


def grelu_wandb_artifacts(
    project: str,
    *,
    host: str = _WANDB_HOST,
    type_is: str | None = None,
    type_contains: str | None = None,
    verbose: bool = True,
) -> list[str]:
    """List the artifacts of a Weights & Biases project of the gReLU model zoo (network)."""
    return _call("grelu.resources.wandb", "artifacts", _WANDB_HINT, project=_text(project, "project"), host=host,
                 type_is=type_is, type_contains=type_contains)


def grelu_wandb_models(project: str, *, host: str = _WANDB_HOST, verbose: bool = True) -> list[str]:
    """List the models of a Weights & Biases project of the gReLU model zoo (network)."""
    return _call("grelu.resources.wandb", "models", _WANDB_HINT, project=_text(project, "project"), host=host)


def grelu_wandb_datasets(project: str, *, host: str = _WANDB_HOST, verbose: bool = True) -> list[str]:
    """List the datasets of a Weights & Biases project of the gReLU model zoo (network)."""
    return _call("grelu.resources.wandb", "datasets", _WANDB_HINT, project=_text(project, "project"), host=host)


def grelu_wandb_runs(
    project: str,
    *,
    host: str = _WANDB_HOST,
    field: str = "id",
    filters: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> list[Any]:
    """List one attribute of every run of a Weights & Biases project of the gReLU model zoo (network)."""
    return _call("grelu.resources.wandb", "runs", _WANDB_HINT, project=_text(project, "project"), host=host,
                 field=field, filters=dict(filters) if filters is not None else None)


def grelu_wandb_get_artifact(
    name: str,
    project: str,
    *,
    host: str = _WANDB_HOST,
    alias: str = "latest",
    verbose: bool = True,
) -> Any:
    """Retrieve a Weights & Biases artifact of the gReLU model zoo by name and alias (network download)."""
    return _call("grelu.resources.wandb", "get_artifact", _WANDB_HINT, name=_text(name, "name"),
                 project=_text(project, "project"), host=host, alias=alias)


def grelu_wandb_get_dataset_by_model(
    model_name: str,
    project: str,
    *,
    host: str = _WANDB_HOST,
    alias: str = "latest",
    verbose: bool = True,
) -> list[str]:
    """List the datasets linked to a model of the gReLU model zoo on Weights & Biases (network)."""
    return _call("grelu.resources.wandb", "get_dataset_by_model", _WANDB_HINT, model_name=_text(model_name, "model_name"),
                 project=_text(project, "project"), host=host, alias=alias)


def grelu_wandb_get_model_by_dataset(
    dataset_name: str,
    project: str,
    *,
    host: str = _WANDB_HOST,
    alias: str = "latest",
    verbose: bool = True,
) -> list[str]:
    """List the models trained on a dataset of the gReLU model zoo on Weights & Biases (network)."""
    return _call("grelu.resources.wandb", "get_model_by_dataset", _WANDB_HINT,
                 dataset_name=_text(dataset_name, "dataset_name"), project=_text(project, "project"), host=host,
                 alias=alias)


def grelu_wandb_load_model(
    project: str,
    model_name: str,
    *,
    device: str | int = "cpu",
    host: str = _WANDB_HOST,
    alias: str = "latest",
    checkpoint_file: str = "model.ckpt",
    weights_only: bool = False,
    verbose: bool = True,
) -> Any:
    """Load a gReLU model from a Weights & Biases project of the model zoo (network download)."""
    return _call("grelu.resources.wandb", "load_model", _WANDB_HINT, project=_text(project, "project"),
                 model_name=_text(model_name, "model_name"), device=device, host=host, alias=alias,
                 checkpoint_file=checkpoint_file, weights_only=weights_only)


_MAX_VALUES = 100_000
_MAX_ROWS = 1000


def _describe(value: Any) -> Any:
    """JSON-friendly description of a plotted input: its values when small, otherwise its shape."""
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if hasattr(value, "columns") and hasattr(value, "to_dict"):
        table = value.head(_MAX_ROWS)
        return {"columns": [str(c) for c in value.columns], "rows": int(len(value)), "values": table.to_dict("list")}
    if hasattr(value, "tolist"):
        size = int(getattr(value, "size", 0))
        shape = [int(s) for s in getattr(value, "shape", ())]
        if size > _MAX_VALUES:
            return {"shape": shape, "truncated": True}
        return {"shape": shape, "values": value.tolist()}
    if isinstance(value, (list, tuple)):
        return list(value) if len(value) <= _MAX_VALUES else {"length": len(value), "truncated": True}
    return str(type(value).__name__)


def _numeric(value: Any, name: str, ndim: int | tuple[int, ...]) -> Any:
    np = _upstream("numpy")
    allowed = (ndim,) if isinstance(ndim, int) else ndim
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a numeric array.") from exc
    if array.ndim not in allowed or array.size == 0 or not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must be a finite array with {' or '.join(map(str, allowed))} dimensions.")
    return array


def _plot(function: str, output_path: str, inputs: Mapping[str, Any], options: Mapping[str, Any]) -> dict[str, Any]:
    """Draw one gReLU plot, save it as an image and return the image path with the plotted data."""
    image = _new_image(output_path, "output_path")
    plt = _upstream("matplotlib.pyplot", _PLOT_HINT)
    visualize = _upstream("grelu.visualize", _HINT)
    before = set(plt.get_fignums())
    try:
        getattr(visualize, function)(**inputs, **options)
        plt.gcf().savefig(image, bbox_inches="tight")
    except Exception as exc:
        raise BackendExecutionError(f"gReLU {function} failed: {exc}") from exc
    finally:
        for number in set(plt.get_fignums()) - before:
            plt.close(number)
    return {
        "status": "completed",
        "backend": "grelu",
        "version": _version(),
        "function": f"grelu.visualize.{function}",
        "image_path": image,
        "data": {key: _describe(value) for key, value in inputs.items()},
    }


def grelu_plot_distribution(
    values: Sequence[float],
    output_path: str,
    *,
    title: str = "metric",
    method: str = "histogram",
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a histogram, density or box plot of a list of values (gReLU plot_distribution)."""
    return _plot("plot_distribution", output_path, {"values": list(values)},
                 {"title": title, "method": method, "figsize": figsize})


def grelu_plot_pred_distribution(
    preds: Any,
    labels: Any,
    output_path: str,
    *,
    tasks: Sequence[int] | Sequence[str] | None = None,
    bins: Sequence[int] | None = None,
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the distribution of model predictions for the positive and negative labels of each task (gReLU plot_pred_distribution)."""
    return _plot("plot_pred_distribution", output_path, {"preds": _numeric(preds, "preds", 3),
                                                         "labels": _numeric(labels, "labels", 3)},
                 {"tasks": tasks, "bins": bins, "figsize": figsize})


def grelu_plot_pred_scatter(
    preds: Any,
    labels: Any,
    output_path: str,
    *,
    tasks: Sequence[int] | Sequence[str] | None = None,
    bins: Sequence[int] | None = None,
    density: bool = False,
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Scatter plot of predictions against measured labels for each task (gReLU plot_pred_scatter)."""
    return _plot("plot_pred_scatter", output_path, {"preds": _numeric(preds, "preds", 3),
                                                    "labels": _numeric(labels, "labels", 3)},
                 {"tasks": tasks, "bins": bins, "density": density, "figsize": figsize})


def grelu_plot_binary_preds(
    preds: Any,
    labels: Any,
    output_path: str,
    *,
    tasks: Sequence[int] | Sequence[str] | None = None,
    bins: Sequence[int] | None = None,
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot predicted probabilities of binary labels as box plots per class and task (gReLU plot_binary_preds)."""
    return _plot("plot_binary_preds", output_path, {"preds": _numeric(preds, "preds", 3),
                                                    "labels": _numeric(labels, "labels", 3)},
                 {"tasks": tasks, "bins": bins, "figsize": figsize})


def grelu_plot_calibration_curve(
    probs: Any,
    labels: Any,
    output_path: str,
    *,
    tasks: Sequence[int] | Sequence[str] | None = None,
    bins: Sequence[int] | None = None,
    aggregate: bool = True,
    figsize: tuple[float, float] = (4, 3),
    show_legend: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the calibration curve of predicted probabilities against observed binary labels (gReLU plot_calibration_curve)."""
    return _plot("plot_calibration_curve", output_path, {"probs": _numeric(probs, "probs", (2, 3)),
                                                         "labels": _numeric(labels, "labels", 3)},
                 {"tasks": tasks, "bins": bins, "aggregate": aggregate, "figsize": figsize,
                  "show_legend": show_legend})


def grelu_plot_evolution(
    df: Any,
    output_path: str,
    *,
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the change in scores and predictions over the rounds of a directed-evolution run (gReLU plot_evolution)."""
    return _plot("plot_evolution", output_path, {"df": df}, {"figsize": figsize})


def grelu_plot_gc_match(
    positives: Any,
    negatives: Any,
    output_path: str,
    *,
    binwidth: float = 0.1,
    genome: str = "hg38",
    figsize: tuple[float, float] = (4, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of GC content of positive and GC-matched negative intervals (gReLU plot_gc_match; genome via genomepy)."""
    return _plot("plot_gc_match", output_path, {"positives": positives, "negatives": negatives},
                 {"binwidth": binwidth, "genome": genome, "figsize": figsize})


def grelu_plot_attributions(
    attrs: Any,
    output_path: str,
    *,
    start_pos: int = 0,
    end_pos: int = -1,
    figsize: tuple[float, float] = (20, 2),
    ticks: int = 10,
    highlight_centers: Sequence[int] | None = None,
    highlight_width: int | Sequence[int] = 5,
    highlight_positions: Sequence[int] | None = None,
    ylim: tuple[float, float] | None = None,
    facecolor: str | None = "yellow",
    edgecolor: str | None = None,
    alpha: float | None = 0.15,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot per-base attribution scores of one sequence as a bar track with optional highlights (gReLU plot_attributions)."""
    return _plot("plot_attributions", output_path, {"attrs": _numeric(attrs, "attrs", 2)},
                 {"start_pos": start_pos, "end_pos": end_pos, "figsize": figsize, "ticks": ticks,
                  "highlight_centers": highlight_centers, "highlight_width": highlight_width,
                  "highlight_positions": highlight_positions, "ylim": ylim, "facecolor": facecolor,
                  "edgecolor": edgecolor, "alpha": alpha})


def grelu_plot_ism(
    ism_preds: Any,
    output_path: str,
    *,
    start_pos: int | None = None,
    end_pos: int | None = None,
    figsize: tuple[float, float] = (8, 1.5),
    method: str = "heatmap",
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot an in-silico mutagenesis table as a heatmap or a line plot (gReLU plot_ISM)."""
    return _plot("plot_ISM", output_path, {"ism_preds": ism_preds},
                 {"start_pos": start_pos, "end_pos": end_pos, "figsize": figsize, "method": method})


def grelu_plot_tracks(
    tracks: Any,
    output_path: str,
    *,
    start_pos: int = 0,
    end_pos: int | None = None,
    titles: Sequence[str] | None = None,
    figsize: tuple[float, float] = (20, 1.5),
    highlight_intervals: Any = None,
    facecolor: str | None = "yellow",
    edgecolor: str | None = None,
    alpha: float | None = 0.15,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot one or more prediction or signal tracks over a sequence window (gReLU plot_tracks)."""
    return _plot("plot_tracks", output_path, {"tracks": _numeric(tracks, "tracks", 2)},
                 {"start_pos": start_pos, "end_pos": end_pos, "titles": titles, "figsize": figsize,
                  "highlight_intervals": highlight_intervals, "facecolor": facecolor, "edgecolor": edgecolor,
                  "alpha": alpha})


def grelu_plot_attention_matrix(
    attn: Any,
    output_path: str,
    *,
    start_pos: int = 0,
    end_pos: int | None = None,
    highlight_intervals: Any = None,
    figsize: tuple[float, float] = (5, 4),
    verbose: bool = True,
) -> dict[str, Any]:
    """Heatmap of an attention matrix of a transformer model over a sequence window (gReLU plot_attention_matrix)."""
    return _plot("plot_attention_matrix", output_path, {"attn": _numeric(attn, "attn", 2)},
                 {"start_pos": start_pos, "end_pos": end_pos, "highlight_intervals": highlight_intervals,
                  "figsize": figsize})


def grelu_plot_motif(
    motif: Any,
    output_path: str,
    *,
    name: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a sequence logo of a motif given as a position matrix with four columns (gReLU plot_motif)."""
    return _plot("plot_motif", output_path, {"motif": _numeric(motif, "motif", 2)}, {"name": name})


def grelu_plot_position_effect(
    preds: Any,
    positions: Sequence[int],
    output_path: str,
    *,
    title: str | None = None,
    xlab: str | None = None,
    figsize: tuple[float, float] = (6, 3),
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot how the model output changes with the position of an inserted motif (gReLU plot_position_effect)."""
    return _plot("plot_position_effect", output_path,
                 {"preds": _numeric(preds, "preds", 2), "positions": list(positions)},
                 {"title": title, "xlab": xlab, "figsize": figsize})


__all__ = [
    "grelu_check_chrom_ends",
    "grelu_check_reference",
    "grelu_download_dataset",
    "grelu_download_model",
    "grelu_filter_chrom_ends",
    "grelu_get_base_models",
    "grelu_get_dataset_info",
    "grelu_get_datasets_by_model",
    "grelu_get_gc_matched_intervals",
    "grelu_get_genome",
    "grelu_get_jaspar",
    "grelu_get_model_info",
    "grelu_get_models_by_dataset",
    "grelu_intervals_to_strings",
    "grelu_list_datasets",
    "grelu_list_models",
    "grelu_load_model",
    "grelu_make_insertion_bigwig",
    "grelu_plot_attention_matrix",
    "grelu_plot_attributions",
    "grelu_plot_binary_preds",
    "grelu_plot_calibration_curve",
    "grelu_plot_distribution",
    "grelu_plot_evolution",
    "grelu_plot_gc_match",
    "grelu_plot_ism",
    "grelu_plot_motif",
    "grelu_plot_position_effect",
    "grelu_plot_pred_distribution",
    "grelu_plot_pred_scatter",
    "grelu_plot_tracks",
    "grelu_read_gtf",
    "grelu_read_sizes",
    "grelu_run_tomtom",
    "grelu_variant_to_seqs",
    "grelu_wandb_artifacts",
    "grelu_wandb_datasets",
    "grelu_wandb_get_artifact",
    "grelu_wandb_get_dataset_by_model",
    "grelu_wandb_get_model_by_dataset",
    "grelu_wandb_load_model",
    "grelu_wandb_models",
    "grelu_wandb_projects",
    "grelu_wandb_runs",
]
_register_public_verbose(globals())
