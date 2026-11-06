"""gReLU training, data and scoring operations (pack P36): augmentation, datasets, model building, training,
fine-tuning, evaluation, prediction, embeddings, metrics, label and prediction transforms, resources.

Each function builds the gReLU object that carries the operation (``Augmenter``, ``LabeledSeqDataset`` and its
subclasses, ``LightningModel``, ``LightningModelEnsemble``, the metric and transform classes) and calls it with
the caller's arguments and gReLU's own defaults. Models are passed as local checkpoint paths (written by
``build_sequence_model`` or by training) and datasets as a JSON-friendly ``dataset`` mapping, so every
argument and result is serialisable:

* ``{"table": "x.tsv"}`` - a delimited table, sequences (or ``chrom, start, end`` intervals) in the first
  column(s) and one label column per task (``DFSeqDataset``; ``.csv`` is comma separated, other suffixes tab);
* ``{"sequences": [...], "labels": [[[...]]], "tasks": [...]}`` - arrays, labels shaped tasks x length
  (``LabeledSeqDataset``);
* ``{"intervals": "x.bed", "bigwigs": ["a.bw"], "genome": "hg38 or a FASTA with .fai/.sizes"}`` -
  ``BigWigSeqDataset``;
* ``{"anndata": "x.h5ad", "label_key": null}`` - ``AnnDataSeqDataset``.

Extra keys (``seq_len``, ``genome``, ``end``, ``rc``, ``max_seq_shift``, ``max_pair_shift``, ``label_len``,
``label_aggfunc``, ``bin_size``, ``min_label_clip``, ``max_label_clip``, ``label_transform_func``, ``seed``,
``augment_mode``) are passed on. Training runs on the CPU unless ``device`` says otherwise. Nothing is
downloaded here except by the ``*_pretrained`` architectures, which fetch published weights. Global NumPy,
torch and Python random state is restored after every call. gReLU is imported when a function runs; a
missing package raises ``BackendUnavailableError`` with the ``bio-grelu`` hint.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import io
import random
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-grelu"
_ERRORS = (
    AssertionError,
    AttributeError,
    IndexError,
    KeyError,
    NotImplementedError,
    OSError,
    RuntimeError,
    TypeError,
    ValueError,
    ImportError,
)
_ARCHITECTURES = {
    "conv": "ConvModel",
    "dilated_conv": "DilatedConvModel",
    "conv_gru": "ConvGRUModel",
    "conv_transformer": "ConvTransformerModel",
    "conv_mlp": "ConvMLPModel",
    "borzoi": "BorzoiModel",
    "borzoi_pretrained": "BorzoiPretrainedModel",
    "enformer": "EnformerModel",
    "enformer_pretrained": "EnformerPretrainedModel",
}
_DATASET_KEYS = (
    "seq_len",
    "genome",
    "end",
    "rc",
    "max_seq_shift",
    "max_pair_shift",
    "label_len",
    "label_aggfunc",
    "bin_size",
    "min_label_clip",
    "max_label_clip",
    "label_transform_func",
    "seed",
    "augment_mode",
)
_MAX_VALUES = 500_000
_LOG_TAIL = 4000


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _provenance(operation: str) -> dict[str, str]:
    return {
        "package": "grelu",
        "version": importlib.metadata.version("grelu"),
        "operation": operation,
    }


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return str(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def _dna_list(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of DNA strings.")
    if len(values) > 1_000_000:
        raise ConfigurationError(f"{name} holds more than 1000000 sequences.")
    out = []
    for item in values:
        if not isinstance(item, str) or not item or set(item.upper()) - set("ACGTN"):
            raise ConfigurationError(f"{name} must contain DNA strings (ACGTN) only.")
        out.append(item.upper())
    return out


def _new_file(value: object, name: str, suffixes: Sequence[str]) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.suffix.lower() not in suffixes:
        raise ConfigurationError(f"{name} must end with {' or '.join(suffixes)}.")
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}")
    return path


def _new_dir(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a folder path.")
    path = Path(value).expanduser().resolve()
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ConfigurationError(f"{name} must be a new or empty folder: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The parent folder of {name} does not exist: {path.parent}")
    return path


@contextlib.contextmanager
def _preserved() -> Iterator[None]:
    """Restore NumPy, torch and Python random state after an upstream call."""

    np = require("numpy", _EXTRA)
    torch = require("torch", _EXTRA)
    states = (np.random.get_state(), torch.random.get_rng_state(), random.getstate())
    try:
        yield
    finally:
        np.random.set_state(states[0])
        torch.random.set_rng_state(states[1])
        random.setstate(states[2])


def _run(operation: str, function: Any, *args: Any, **kwargs: Any) -> tuple[Any, str]:
    """Call gReLU with its console output captured; upstream failures become BackendExecutionError."""

    buffer = io.StringIO()
    try:
        with _preserved(), contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
            result = function(*args, **kwargs)
    except _ERRORS as exc:
        raise BackendExecutionError(f"gReLU {operation} failed: {exc}") from exc
    return result, buffer.getvalue()[-_LOG_TAIL:]


def _grelu() -> tuple[Any, Any]:
    return require("grelu.lightning", _EXTRA), require("grelu.data.dataset", _EXTRA)


def _load_model(checkpoint: object) -> Any:
    lightning = require("grelu.lightning", _EXTRA)
    path = input_path(checkpoint)  # type: ignore[arg-type]
    model, _ = _run(
        "load_from_checkpoint",
        lightning.LightningModel.load_from_checkpoint,
        str(path),
        map_location="cpu",
        weights_only=False,
    )
    return model


def _array_result(array: Any, output_path: object, name: str) -> dict[str, Any]:
    np = require("numpy", _EXTRA)
    if output_path is not None:
        target = _new_file(output_path, "output_path", (".npy",))
        np.save(target, array)
        return {"shape": list(array.shape), "output": str(target)}
    if array.size > _MAX_VALUES:
        raise ConfigurationError(
            f"The {name} array has {array.size} values; give output_path (.npy) for results above {_MAX_VALUES}."
        )
    return {"shape": list(array.shape), name: array.tolist()}


def _device(value: object) -> str:
    if value == "cpu":
        return "cpu"
    if (isinstance(value, str) and value.startswith("cuda")) or (
        isinstance(value, str) and value.isdigit()
    ):
        return str(value)
    raise ConfigurationError("device must be 'cpu', 'cuda', 'cuda:N' or a GPU index.")


def _frame(path: object, name: str) -> Any:
    pandas = require("pandas", _EXTRA)
    resolved = input_path(path)  # type: ignore[arg-type]
    sep = "," if resolved.name.endswith((".csv", ".csv.gz")) else "\t"
    try:
        return pandas.read_csv(resolved, sep=sep, dtype={"chrom": str})
    except _ERRORS as exc:
        raise BackendExecutionError(f"{name} could not be read: {exc}") from exc


def _dataset(spec: object, name: str = "dataset") -> Any:
    """Build a gReLU labelled dataset from the ``dataset`` mapping."""

    if not isinstance(spec, Mapping):
        raise ConfigurationError(f"{name} must be a mapping (see the module description).")
    np = require("numpy", _EXTRA)
    _, data = _grelu()
    extra = {k: spec[k] for k in _DATASET_KEYS if k in spec}
    unknown = (
        set(spec)
        - set(_DATASET_KEYS)
        - {"table", "sequences", "labels", "tasks", "intervals", "bigwigs", "anndata", "label_key"}
    )
    if unknown:
        raise ConfigurationError(f"{name} has unknown keys: {', '.join(sorted(unknown))}.")
    kinds = [k for k in ("table", "sequences", "intervals", "anndata") if k in spec]
    if len(kinds) != 1:
        raise ConfigurationError(
            f"{name} needs exactly one of table, sequences, intervals, anndata."
        )
    kind = kinds[0]
    if kind == "table":
        for key in (
            "label_len",
            "label_aggfunc",
            "bin_size",
            "min_label_clip",
            "max_label_clip",
            "label_transform_func",
            "max_pair_shift",
        ):
            extra.pop(key, None)
        built, _log = _run(
            "DFSeqDataset", data.DFSeqDataset, _frame(spec["table"], "table"), **extra
        )
    elif kind == "sequences":
        seqs = _dna_list(spec["sequences"], f"{name}.sequences")
        labels = np.asarray(spec.get("labels"), dtype="float32")
        if labels.ndim != 3 or labels.shape[0] != len(seqs):
            raise ConfigurationError(f"{name}.labels must have shape (sequences, tasks, length).")
        tasks = spec.get("tasks") or [f"task{i}" for i in range(labels.shape[1])]
        if len(tasks) != labels.shape[1]:
            raise ConfigurationError(f"{name}.tasks must name every label task.")
        built, _log = _run(
            "LabeledSeqDataset",
            data.LabeledSeqDataset,
            seqs,
            labels,
            [str(t) for t in tasks],
            **extra,
        )
    elif kind == "intervals":
        pandas = require("pandas", _EXTRA)
        bigwigs = spec.get("bigwigs")
        if isinstance(bigwigs, (str, bytes)) or not isinstance(bigwigs, Sequence) or not bigwigs:
            raise ConfigurationError(f"{name}.bigwigs must be a non-empty list of bigWig paths.")
        bed = pandas.read_csv(
            input_path(spec["intervals"]),
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            comment="#",
            dtype={"chrom": str},
        )
        built, _log = _run(
            "BigWigSeqDataset",
            data.BigWigSeqDataset,
            bed,
            [str(input_path(p)) for p in bigwigs],
            tasks=spec.get("tasks"),
            **extra,
        )
    else:
        anndata = require("anndata", _EXTRA)
        for key in (
            "label_len",
            "label_aggfunc",
            "bin_size",
            "min_label_clip",
            "max_label_clip",
            "label_transform_func",
            "max_pair_shift",
        ):
            extra.pop(key, None)
        adata, _log = _run("read_h5ad", anndata.read_h5ad, str(input_path(spec["anndata"])))
        built, _log = _run(
            "AnnDataSeqDataset",
            data.AnnDataSeqDataset,
            adata,
            label_key=spec.get("label_key"),
            **extra,
        )
    return built


def _model_summary(model: Any) -> dict[str, Any]:
    return {
        "model_type": model.model_params.get("model_type"),
        "n_tasks": int(model.model_params["n_tasks"]),
        "parameters": int(model.count_params()),
        "tasks": list(model.data_params.get("tasks", {}).get("name", []))
        if model.data_params
        else [],
    }


# ---------------------------------------------------------------------------------------------
# augmentation and datasets
# ---------------------------------------------------------------------------------------------


def augment_sequences(
    sequences: Sequence[str],
    *,
    labels: Sequence[Sequence[Sequence[float]]] | None = None,
    seq_len: int | None = None,
    label_len: int | None = None,
    rc: bool = False,
    max_seq_shift: int = 0,
    max_pair_shift: int = 0,
    n_mutated_seqs: int = 0,
    n_mutated_bases: int | None = None,
    protect: Sequence[int] = (),
    seed: int | None = None,
    mode: str = "serial",
    n_draws: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate the training augmentations of DNA sequences (and labels): shift (jitter), reverse complement
    and random mutation (gReLU ``Augmenter``, ``shift``, ``rc_seq``, ``rc_label``).

    ``rc=True`` adds the reverse-complement copy (and reverses the label along its length); ``max_seq_shift``
    shifts the sequence window by up to that many bases either way; ``max_pair_shift`` shifts sequence and label
    together; ``n_mutated_seqs`` with ``n_mutated_bases`` makes randomly mutated copies (``protect`` lists
    positions to keep). Input sequences must be ``seq_len + 2*(max_seq_shift + max_pair_shift)`` bases long
    (``seq_len`` defaults to the unshifted length) and ``labels`` (optional) are tasks x ``label_len +
    2*max_pair_shift``. ``mode="serial"`` returns every combination for every sequence, ``mode="random"``
    returns ``n_draws`` random draws per sequence. Each output record gives the source index, the augmented
    sequence and, when labels were given, the augmented label.
    """

    np = require("numpy", _EXTRA)
    augment = require("grelu.data.augment", _EXTRA)
    fmt = require("grelu.sequence.format", _EXTRA)
    seqs = _dna_list(sequences, "sequences")
    if len({len(s) for s in seqs}) != 1:
        raise ConfigurationError("sequences must all have the same length.")
    shift_total = 2 * (
        _int(max_seq_shift, "max_seq_shift", 0, 10_000)
        + _int(max_pair_shift, "max_pair_shift", 0, 10_000)
    )
    out_len = (
        _int(seq_len, "seq_len", 1, 10_000_000)
        if seq_len is not None
        else len(seqs[0]) - shift_total
    )
    if out_len < 1 or out_len + shift_total != len(seqs[0]):
        raise ConfigurationError(
            "Sequence length must equal seq_len + 2*(max_seq_shift + max_pair_shift)."
        )
    n_mut = _int(n_mutated_seqs, "n_mutated_seqs", 0, 1000)
    n_bases = (
        None if n_mutated_bases is None else _int(n_mutated_bases, "n_mutated_bases", 0, 10_000)
    )
    if n_mut > 0 and not n_bases:
        raise ConfigurationError("n_mutated_bases (> 0) is required when n_mutated_seqs > 0.")
    kind = _choice(mode, "mode", ("serial", "random"))
    draws = _int(n_draws, "n_draws", 1, 10_000)
    label_array = None
    if labels is not None:
        label_array = np.asarray(labels, dtype="float32")
        if label_array.ndim != 3 or label_array.shape[0] != len(seqs):
            raise ConfigurationError("labels must have shape (sequences, tasks, length).")
    pair_len = None
    if label_array is not None:
        pair_len = (
            _int(label_len, "label_len", 1, 10_000_000)
            if label_len is not None
            else label_array.shape[2] - 2 * max_pair_shift
        )
    augmenter, log = _run(
        "Augmenter",
        augment.Augmenter,
        rc=_flag(rc, "rc"),
        max_seq_shift=max_seq_shift,
        max_pair_shift=max_pair_shift,
        n_mutated_seqs=n_mut,
        n_mutated_bases=n_bases,
        protect=[int(p) for p in protect],
        seq_len=out_len,
        label_len=pair_len,
        seed=seed,
        mode=kind,
    )
    n_variants = len(augmenter) if kind == "serial" else draws
    if len(seqs) * n_variants > 100_000:
        raise ConfigurationError("The augmentation would produce more than 100000 records.")
    records: list[dict[str, Any]] = []
    for index, text in enumerate(seqs):
        encoded = fmt.strings_to_indices(text)
        for variant in range(n_variants):
            try:
                with _preserved():
                    if label_array is None:
                        seq = augmenter(variant, encoded)
                        label = None
                    else:
                        seq, label = augmenter(variant, encoded, label_array[index])
            except _ERRORS as exc:
                raise BackendExecutionError(f"gReLU augmentation failed: {exc}") from exc
            row: dict[str, Any] = {
                "source_index": index,
                "variant": variant,
                "sequence": fmt.indices_to_strings(seq),
            }
            if label is not None:
                row["label"] = np.asarray(label).tolist()
            records.append(row)
    return {
        "status": "completed",
        "operation": "augment_sequences",
        "mode": kind,
        "variants_per_sequence": n_variants,
        "count": len(records),
        "records": records,
        "log_tail": log if verbose else "",
        "provenance": _provenance("Augmenter"),
    }


def build_labeled_sequence_dataset(
    dataset: Mapping[str, Any],
    *,
    indices: Sequence[int] = (0,),
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a labelled sequence dataset and report it, with selected items (gReLU ``LabeledSeqDataset``,
    ``DFSeqDataset``, ``BigWigSeqDataset``, ``AnnDataSeqDataset``).

    ``dataset`` is described in the module description: a ``table`` (DataFrame-style, one label column per
    task), ``sequences`` with ``labels`` arrays, ``intervals`` with ``bigwigs`` coverage, or an ``anndata``
    file. The result gives the number of sequences, tasks, augmented length and the decoded sequence and
    label of every requested ``indices`` entry (indices run over augmented items).
    """

    np = require("numpy", _EXTRA)
    fmt = require("grelu.sequence.format", _EXTRA)
    built = _dataset(dataset)
    picks = [_int(i, "indices item", 0, 10**9) for i in indices]
    if not picks or len(picks) > 1000:
        raise ConfigurationError("indices must hold between 1 and 1000 entries.")
    total = len(built)
    items = []
    for index in picks:
        if index >= total:
            raise ConfigurationError(
                f"index {index} is out of range for a dataset of {total} items."
            )
        one_hot, label = built[index]
        items.append(
            {
                "index": index,
                "sequence": fmt.indices_to_strings(np.asarray(one_hot.argmax(0))),
                "label": np.asarray(label).tolist(),
            }
        )
    return {
        "status": "completed",
        "operation": "build_labeled_sequence_dataset",
        "class": type(built).__name__,
        "sequences": int(built.n_seqs),
        "tasks": [str(t) for t in built.tasks.index]
        if hasattr(built.tasks, "index")
        else list(built.tasks),
        "n_tasks": int(built.n_tasks),
        "sequence_length": int(built.seq_len),
        "label_length": int(built.label_len),
        "augmented_items": int(total),
        "items": items,
        "provenance": _provenance(type(built).__name__),
    }


def scan_motif_insertions(
    sequences: Sequence[str],
    motifs: Sequence[str],
    *,
    positions: Sequence[int] | None = None,
    max_records: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert each motif at every (or the given) position of background sequences (gReLU ``MotifScanDataset``).

    Motifs overwrite the bases at the position, so the sequence length does not change. The result lists, for
    every sequence, position and motif, the new sequence; pair it with a model prediction
    (``predict_sequence_model``) to get a positional motif-effect scan.
    """

    np = require("numpy", _EXTRA)
    fmt = require("grelu.sequence.format", _EXTRA)
    _, data = _grelu()
    seqs = _dna_list(sequences, "sequences")
    pieces = _dna_list(motifs, "motifs")
    if len({len(s) for s in seqs}) != 1:
        raise ConfigurationError("sequences must all have the same length.")
    cap = _int(max_records, "max_records", 1, 1_000_000)
    where = None if positions is None else [_int(p, "positions item", 0, 10**9) for p in positions]
    scan, _log = _run("MotifScanDataset", data.MotifScanDataset, seqs, pieces, positions=where)
    total = len(scan)
    if total > cap:
        raise ConfigurationError(f"The scan has {total} sequences; raise max_records (now {cap}).")
    records = []
    for index in range(total):
        seq_idx, pos_idx, motif_idx = _split(index, (scan.n_seqs, scan.n_positions, scan.n_alleles))
        position = pos_idx if where is None else where[pos_idx]
        encoded = scan[index]
        records.append(
            {
                "source_index": seq_idx,
                "position": int(position),
                "motif_index": motif_idx,
                "sequence": fmt.indices_to_strings(np.asarray(encoded.argmax(0))),
            }
        )
    return {
        "status": "completed",
        "operation": "scan_motif_insertions",
        "count": total,
        "records": records,
        "provenance": _provenance("MotifScanDataset"),
    }


def _split(index: int, sizes: Sequence[int]) -> list[int]:
    out = []
    for position in range(len(sizes)):
        stride = 1
        for later in sizes[position + 1 :]:
            stride *= later
        out.append(int((index // stride) % sizes[position]))
    return out


def filter_cell_types(
    anndata_path: str | Path,
    output_path: str | Path,
    *,
    cutoff: int = 1000,
    count_key: str = "n_cells",
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop cell types (tasks) composed of few cells from an AnnData file (gReLU ``filter_cells``).

    ``anndata_path`` is an ``.h5ad`` file with intervals in ``.var`` and per-cell-type cell counts in
    ``.obs[count_key]``; cell types with fewer than ``cutoff`` cells are removed and the result is written to
    the new ``output_path``.
    """

    anndata = require("anndata", _EXTRA)
    preprocess = require("grelu.data.preprocess", _EXTRA)
    source = input_path(anndata_path)
    target = _new_file(output_path, "output_path", (".h5ad",))
    limit = _int(cutoff, "cutoff", 0, 10**12)
    key = _text(count_key, "count_key")
    adata, _log = _run("read_h5ad", anndata.read_h5ad, str(source))
    if key not in adata.obs:
        raise ConfigurationError(f"count_key {key!r} is not a column of .obs.")
    before = list(map(str, adata.obs_names))
    kept, _log = _run("filter_cells", preprocess.filter_cells, adata, cutoff=limit, count_key=key)
    kept.write_h5ad(target)
    after = list(map(str, kept.obs_names))
    return {
        "status": "completed",
        "operation": "filter_cell_types",
        "output": str(target),
        "cutoff": limit,
        "cell_types_before": len(before),
        "cell_types_after": len(after),
        "removed": [name for name in before if name not in set(after)],
        "provenance": _provenance("filter_cells"),
    }


# ---------------------------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------------------------


def build_sequence_model(
    architecture: str,
    n_tasks: int,
    output_path: str | Path,
    *,
    model_params: Mapping[str, Any] | None = None,
    train_params: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create an untrained gReLU sequence model and save it as a Lightning checkpoint (``LightningModel``).

    ``architecture`` selects the network (gReLU ``grelu.model.models``): ``conv`` (ConvModel),
    ``dilated_conv`` (DilatedConvModel), ``conv_gru`` (ConvGRUModel), ``conv_transformer``
    (ConvTransformerModel), ``conv_mlp`` (ConvMLPModel, needs ``seq_len`` in ``model_params``), ``borzoi``
    (BorzoiModel), ``enformer`` (EnformerModel), and ``borzoi_pretrained`` / ``enformer_pretrained``
    (BorzoiPretrainedModel / EnformerPretrainedModel, which download the published weights and are not run
    by the offline tests). ``model_params`` carries the architecture's own arguments (channels, kernel sizes,
    ``crop_len``, ``final_pool_func`` ...); ``train_params`` the training settings (``task``
    binary/multiclass/regression, ``loss``, ``lr``, ``batch_size``, ``max_epochs`` ...). The checkpoint is
    the input of the training, evaluation, prediction and embedding functions.
    """

    name = _ARCHITECTURES.get(_choice(architecture, "architecture", tuple(_ARCHITECTURES)))
    tasks = _int(n_tasks, "n_tasks", 1, 100_000)
    target = _new_file(output_path, "output_path", (".ckpt",))
    lightning, _ = _grelu()
    pl = require("pytorch_lightning", _EXTRA)
    params = {"model_type": name, "n_tasks": tasks, **dict(model_params or {})}
    model, log = _run(
        "LightningModel", lightning.LightningModel, params, train_params=dict(train_params or {})
    )

    def save() -> None:
        trainer = pl.Trainer(
            accelerator="cpu",
            devices=1,
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=False,
            enable_model_summary=False,
        )
        trainer.strategy.connect(model)
        trainer.save_checkpoint(str(target))

    _run("save_checkpoint", save)
    return {
        "status": "completed",
        "operation": "build_sequence_model",
        "architecture": architecture,
        "checkpoint": str(target),
        "model": _model_summary(model),
        "log_tail": log if verbose else "",
        "provenance": _provenance(f"models.{name}"),
    }


def replace_model_head(
    checkpoint: str | Path,
    n_tasks: int,
    output_path: str | Path,
    *,
    final_pool_func: str = "avg",
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the output head of a model by a new randomly initialised head with ``n_tasks`` outputs
    (``LightningModel.change_head``) and save the result as a new checkpoint.

    ``final_pool_func`` (``avg``, ``max`` or ``attn``) is the pooling applied to the new head; the trunk
    weights are kept, as the first step of transfer learning.
    """

    model = _load_model(checkpoint)
    tasks = _int(n_tasks, "n_tasks", 1, 100_000)
    pool = _choice(final_pool_func, "final_pool_func", ("avg", "max", "attn"))
    target = _new_file(output_path, "output_path", (".ckpt",))
    before = int(model.model_params["n_tasks"])
    _run("change_head", model.change_head, n_tasks=tasks, final_pool_func=pool)
    pl = require("pytorch_lightning", _EXTRA)

    def save() -> None:
        trainer = pl.Trainer(
            accelerator="cpu",
            devices=1,
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=False,
            enable_model_summary=False,
        )
        trainer.strategy.connect(model)
        trainer.save_checkpoint(str(target))

    _run("save_checkpoint", save)
    return {
        "status": "completed",
        "operation": "replace_model_head",
        "checkpoint": str(target),
        "tasks_before": before,
        "tasks_after": tasks,
        "final_pool_func": pool,
        "provenance": _provenance("LightningModel.change_head"),
    }


def train_sequence_model(
    checkpoint: str | Path,
    train_dataset: Mapping[str, Any],
    val_dataset: Mapping[str, Any],
    output_dir: str | Path,
    *,
    tune: bool = False,
    freeze_embedding: bool = False,
    final_pool_func: str = "avg",
    max_epochs: int | None = None,
    batch_size: int | None = None,
    learning_rate: float | None = None,
    device: str = "cpu",
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train a model on a dataset, or fine-tune a pretrained model (``LightningModel.train_on_dataset`` /
    ``tune_on_dataset``).

    ``checkpoint`` is a model from ``build_sequence_model`` (or an earlier training). With ``tune=False`` the
    model is trained as it is. With ``tune=True`` a new head with the dataset's number of tasks is built
    first (``final_pool_func``), optionally keeping the trunk fixed (``freeze_embedding``), and then trained.
    ``max_epochs``, ``batch_size`` and ``learning_rate`` override the checkpoint's training settings. Lightning
    writes its checkpoints and logs under ``output_dir`` (new or empty); the best and last checkpoint paths
    and the final metrics are returned.
    """

    model = _load_model(checkpoint)
    train = _dataset(train_dataset, "train_dataset")
    val = _dataset(val_dataset, "val_dataset")
    out = _new_dir(output_dir, "output_dir")
    out.mkdir(exist_ok=True)
    flag, freeze = _flag(tune, "tune"), _flag(freeze_embedding, "freeze_embedding")
    pool = _choice(final_pool_func, "final_pool_func", ("avg", "max", "attn"))
    model.train_params["save_dir"] = str(out)
    model.train_params["devices"] = _device(device)
    model.train_params["num_workers"] = 0
    model.train_params["logger"] = "csv"
    model.train_params["checkpoint"] = True
    if max_epochs is not None:
        model.train_params["max_epochs"] = _int(max_epochs, "max_epochs", 1, 100_000)
    if batch_size is not None:
        model.train_params["batch_size"] = _int(batch_size, "batch_size", 1, 1_000_000)
    if learning_rate is not None:
        if (
            isinstance(learning_rate, bool)
            or not isinstance(learning_rate, (int, float))
            or learning_rate <= 0
        ):
            raise ConfigurationError("learning_rate must be a positive number.")
        model.train_params["lr"] = float(learning_rate)
    if seed is not None:
        require("pytorch_lightning", _EXTRA).seed_everything(
            _int(seed, "seed", 0, 2**32 - 1), workers=True
        )
    if flag:
        trainer, log = _run(
            "tune_on_dataset",
            model.tune_on_dataset,
            train,
            val,
            final_pool_func=pool,
            freeze_embedding=freeze,
        )
    else:
        trainer, log = _run("train_on_dataset", model.train_on_dataset, train, val)
    callback = trainer.checkpoint_callback
    metrics = {k: float(v) for k, v in trainer.callback_metrics.items()}
    return {
        "status": "completed",
        "operation": "train_sequence_model",
        "mode": "fine-tune" if flag else "train",
        "output_dir": str(out),
        "best_checkpoint": callback.best_model_path or None,
        "last_checkpoint": callback.last_model_path or None,
        "epochs_run": int(trainer.current_epoch),
        "metrics": metrics,
        "model": _model_summary(model),
        "log_tail": log if verbose else "",
        "provenance": _provenance(
            "LightningModel.tune_on_dataset" if flag else "LightningModel.train_on_dataset"
        ),
    }


def evaluate_sequence_model(
    checkpoint: str | Path,
    dataset: Mapping[str, Any],
    *,
    device: str = "cpu",
    batch_size: int = 256,
    output_checkpoint: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate a model on a labelled test dataset and return one row of metrics per task
    (``LightningModel.test_on_dataset``).

    The metrics follow the model's ``task`` setting (for example ``test_mse`` and ``test_pearson`` for
    regression, AUROC/average precision/accuracy for binary tasks). ``output_checkpoint`` (optional, new
    ``.ckpt`` path) saves a copy of the model that stores the test data parameters and performance.
    """

    model = _load_model(checkpoint)
    data = _dataset(dataset)
    target = (
        None
        if output_checkpoint is None
        else _new_file(output_checkpoint, "output_checkpoint", (".ckpt",))
    )
    frame, log = _run(
        "test_on_dataset",
        model.test_on_dataset,
        data,
        devices=_device(device),
        num_workers=0,
        batch_size=_int(batch_size, "batch_size", 1, 1_000_000),
        write_path=None if target is None else str(target),
    )
    rows = [
        {"task": str(name), **{k: float(v) for k, v in row.items()}}
        for name, row in frame.iterrows()
    ]
    return {
        "status": "completed",
        "operation": "evaluate_sequence_model",
        "tasks": len(rows),
        "metrics": rows,
        "checkpoint": None if target is None else str(target),
        "log_tail": log if verbose else "",
        "provenance": _provenance("LightningModel.test_on_dataset"),
    }


def predict_sequence_model(
    checkpoints: str | Path | Sequence[str | Path],
    sequences: Sequence[str],
    *,
    model_names: Sequence[str] | None = None,
    rc: bool = False,
    max_seq_shift: int = 0,
    device: str = "cpu",
    batch_size: int = 256,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict with one model, or with several models as an ensemble, on DNA sequences
    (``LightningModel.predict_on_dataset`` with ``SeqDataset``; ``LightningModelEnsemble`` for several).

    ``checkpoints`` is one checkpoint path or a list. For a list the predictions of all models are
    concatenated along the task axis in the given order and the tasks are prefixed with ``model_names``
    (default ``model0``, ``model1`` ...). The output has shape sequences x tasks x length. ``rc`` and
    ``max_seq_shift`` average over the reverse-complement and shifted copies (gReLU's ``augment_aggfunc``,
    default mean). ``output_path`` (``.npy``) saves large results instead of returning them.
    """

    lightning, data = _grelu()
    paths = [checkpoints] if isinstance(checkpoints, (str, Path)) else list(checkpoints)
    if not paths or len(paths) > 100:
        raise ConfigurationError("checkpoints must hold between 1 and 100 paths.")
    models = [_load_model(p) for p in paths]
    seqs = _dna_list(sequences, "sequences")
    seq_set, _ = _run(
        "SeqDataset", data.SeqDataset, seqs, rc=_flag(rc, "rc"), max_seq_shift=max_seq_shift
    )
    settings = {
        "devices": _device(device),
        "num_workers": 0,
        "batch_size": _int(batch_size, "batch_size", 1, 1_000_000),
    }
    if len(models) == 1:
        predictor: Any = models[0]
        task_names = _model_summary(models[0])["tasks"]
    else:
        names = list(model_names) if model_names is not None else None
        if names is not None and len(names) != len(models):
            raise ConfigurationError("model_names needs one name per checkpoint.")
        predictor, _ = _run(
            "LightningModelEnsemble", lightning.LightningModelEnsemble, models, names
        )
        task_names = list(predictor.data_params["tasks"].get("name", []))
    array, log = _run("predict_on_dataset", predictor.predict_on_dataset, seq_set, **settings)
    result = _array_result(array, output_path, "predictions")
    return {
        "status": "completed",
        "operation": "predict_sequence_model",
        "models": len(models),
        "sequences": len(seqs),
        "tasks": task_names,
        **result,
        "log_tail": log if verbose else "",
        "provenance": _provenance("LightningModel.predict_on_dataset"),
    }


def embed_sequences(
    checkpoint: str | Path,
    sequences: Sequence[str],
    *,
    device: str = "cpu",
    batch_size: int = 256,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the model's embeddings (the output of its convolutional/attention trunk, before the head) for
    DNA sequences (``LightningModel.embed_on_dataset``).

    The embedding array has shape sequences x channels x length (length after pooling and cropping).
    ``output_path`` (``.npy``) saves large results instead of returning them.
    """

    model = _load_model(checkpoint)
    _, data = _grelu()
    seqs = _dna_list(sequences, "sequences")
    seq_set, _ = _run("SeqDataset", data.SeqDataset, seqs)
    dev: Any = _device(device)
    array, log = _run(
        "embed_on_dataset",
        model.embed_on_dataset,
        seq_set,
        device=int(dev) if dev.isdigit() else dev,
        num_workers=0,
        batch_size=_int(batch_size, "batch_size", 1, 1_000_000),
    )
    result = _array_result(array, output_path, "embeddings")
    return {
        "status": "completed",
        "operation": "embed_sequences",
        "sequences": len(seqs),
        **result,
        "log_tail": log if verbose else "",
        "provenance": _provenance("LightningModel.embed_on_dataset"),
    }


# ---------------------------------------------------------------------------------------------
# metrics, transforms, scores, resources
# ---------------------------------------------------------------------------------------------


def compute_prediction_metric(
    metric: str,
    predictions: Sequence[Sequence[Sequence[float]]],
    targets: Sequence[Sequence[Sequence[float]]],
    *,
    average: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a gReLU evaluation metric on predictions and labels, both shaped sequences x tasks x length
    (gReLU ``MSE`` and ``BestF1`` metric classes with their ``update`` and ``compute``).

    ``metric`` is ``mse`` (mean squared error per task) or ``best_f1`` (best F1 over all thresholds per task;
    targets must be 0/1 and predictions probabilities). With ``average=True`` the mean over tasks is
    returned in ``value``; ``per_task`` always lists the task values.
    """

    np = require("numpy", _EXTRA)
    torch = require("torch", _EXTRA)
    metrics = require("grelu.lightning.metrics", _EXTRA)
    kind = _choice(metric, "metric", ("mse", "best_f1"))
    pred = np.asarray(predictions, dtype="float32")
    true = np.asarray(targets, dtype="float32")
    if pred.ndim != 3 or pred.shape != true.shape or pred.size == 0 or pred.size > 5_000_000:
        raise ConfigurationError("predictions and targets must have the same non-empty 3-D shape.")
    flag = _flag(average, "average")
    tasks = pred.shape[1]
    if kind == "mse":
        calc, _ = _run("MSE", metrics.MSE, num_outputs=tasks, average=False)
    else:
        calc, _ = _run("BestF1", metrics.BestF1, num_labels=tasks, average=False)
    _run("update", calc.update, torch.tensor(pred), torch.tensor(true))
    value, _ = _run("compute", calc.compute)
    per_task = [float(v) for v in torch.as_tensor(value).flatten()]
    return {
        "status": "completed",
        "operation": "compute_prediction_metric",
        "metric": kind,
        "tasks": tasks,
        "per_task": per_task,
        "value": float(np.mean(per_task)) if flag else None,
        "provenance": _provenance("MSE" if kind == "mse" else "BestF1"),
    }


def transform_labels(
    labels: Sequence[Sequence[Sequence[float]]],
    *,
    min_clip: float | None = None,
    max_clip: float | None = None,
    transform_func: str | None = None,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip and/or transform label values (gReLU ``LabelTransform``).

    ``labels`` is sequences x tasks x length. Values are first clipped to ``[min_clip, max_clip]``, then
    ``transform_func`` (``log`` or ``log1p``) is applied.
    """

    np = require("numpy", _EXTRA)
    transforms = require("grelu.transforms.label_transforms", _EXTRA)
    array = np.asarray(labels, dtype="float64")
    if array.ndim != 3 or array.size == 0:
        raise ConfigurationError("labels must be a non-empty (sequences, tasks, length) array.")
    if transform_func is not None:
        _choice(transform_func, "transform_func", ("log", "log1p"))
    transform, _ = _run(
        "LabelTransform",
        transforms.LabelTransform,
        min_clip=min_clip,
        max_clip=max_clip,
        transform_func=transform_func,
    )
    out, _ = _run("forward", transform.forward, array)
    return {
        "status": "completed",
        "operation": "transform_labels",
        "min_clip": min_clip,
        "max_clip": max_clip,
        "transform_func": transform_func,
        **_array_result(np.asarray(out), output_path, "labels"),
        "provenance": _provenance("LabelTransform"),
    }


def score_prediction_specificity(
    predictions: Sequence[Sequence[Sequence[float]]],
    on_tasks: Sequence[int],
    *,
    off_tasks: Sequence[int] | None = None,
    on_aggfunc: str = "mean",
    off_aggfunc: str = "mean",
    off_weight: float = 1.0,
    off_thresh: float | None = None,
    positions: Sequence[int] | None = None,
    length_aggfunc: str = "sum",
    compare_func: str = "divide",
    verbose: bool = True,
) -> dict[str, Any]:
    """Cell-type specificity of model predictions: foreground against background tasks (gReLU ``Specificity``).

    ``predictions`` is sequences x tasks x length. The ``on_tasks`` (indices) and ``off_tasks`` (default: all
    others) are each aggregated over tasks (``on_aggfunc`` / ``off_aggfunc``: sum, mean, min, max) after
    ``length_aggfunc`` over ``positions``; the background is multiplied by ``off_weight`` (only values above
    ``off_thresh`` when given) and combined with ``compare_func`` (``divide``, ``subtract`` or ``log2FC``).
    One score per sequence is returned.
    """

    np = require("numpy", _EXTRA)
    transforms = require("grelu.transforms.prediction_transforms", _EXTRA)
    array = np.asarray(predictions, dtype="float64")
    if array.ndim != 3 or array.size == 0:
        raise ConfigurationError(
            "predictions must be a non-empty (sequences, tasks, length) array."
        )
    aggs = ("sum", "mean", "min", "max")
    on = [_int(t, "on_tasks item", 0, array.shape[1] - 1) for t in on_tasks]
    off = (
        None
        if off_tasks is None
        else [_int(t, "off_tasks item", 0, array.shape[1] - 1) for t in off_tasks]
    )
    if not on:
        raise ConfigurationError("on_tasks must not be empty.")
    if off is None:
        off = [i for i in range(array.shape[1]) if i not in on]
    if not off:
        raise ConfigurationError(
            "off_tasks is empty: the background needs at least one other task."
        )
    spec, _ = _run(
        "Specificity",
        transforms.Specificity,
        on,
        off_tasks=off,
        on_aggfunc=_choice(on_aggfunc, "on_aggfunc", aggs),
        off_aggfunc=_choice(off_aggfunc, "off_aggfunc", aggs),
        off_weight=off_weight,
        off_thresh=off_thresh,
        positions=None if positions is None else [int(p) for p in positions],
        length_aggfunc=_choice(length_aggfunc, "length_aggfunc", aggs),
        compare_func=_choice(compare_func, "compare_func", ("divide", "subtract", "log2FC")),
    )
    scores, _ = _run("compute", spec.compute, array)
    return {
        "status": "completed",
        "operation": "score_prediction_specificity",
        "on_tasks": on,
        "off_tasks": off,
        "scores": np.asarray(scores).reshape(array.shape[0], -1)[:, 0].tolist(),
        "provenance": _provenance("Specificity"),
    }


def score_sequence_patterns(
    sequences: Sequence[str],
    patterns: Sequence[str],
    *,
    weights: Sequence[float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weighted count of (overlapping) occurrences of subsequences in each DNA sequence (gReLU ``PatternScore``).

    Patterns are regular expressions (for example ``ACGT`` or ``TA[CT]A``) matched on the forward strand with
    overlaps allowed; ``weights`` (default 1 each) multiply the per-pattern counts, which are summed per
    sequence.
    """

    transforms = require("grelu.transforms.seq_transforms", _EXTRA)
    seqs = _dna_list(sequences, "sequences")
    if isinstance(patterns, (str, bytes)) or not isinstance(patterns, Sequence) or not patterns:
        raise ConfigurationError("patterns must be a non-empty list.")
    pats = [_text(p, "patterns item") for p in patterns]
    if weights is not None and len(weights) != len(pats):
        raise ConfigurationError("weights needs one value per pattern.")
    scorer, _ = _run(
        "PatternScore",
        transforms.PatternScore,
        pats,
        None if weights is None else [float(w) for w in weights],
    )
    scores, _ = _run("forward", scorer.forward, seqs)
    return {
        "status": "completed",
        "operation": "score_sequence_patterns",
        "scores": [float(v) for v in scores],
        "provenance": _provenance("PatternScore"),
    }


def score_sequence_motifs(
    sequences: Sequence[str],
    motifs: Mapping[str, Sequence[Sequence[float]]] | str | Path,
    *,
    names: Sequence[str] | None = None,
    weights: Sequence[float] | None = None,
    pthresh: float = 1e-3,
    rc: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weighted count of PWM motif occurrences (FIMO-style hits) in each DNA sequence (gReLU ``MotifScore``).

    ``motifs`` is a MEME file path or a mapping from motif name to a 4 x L matrix (rows A, C, G, T);
    ``names`` selects motifs from a MEME file; hits need a p-value at or below ``pthresh`` and both strands
    are scanned when ``rc`` is true. The score of a sequence is the sum over motifs of
    ``weight * number of hits`` (weights default to 1). gReLU's ``MotifScore.forward`` returns hit counts per
    (sequence, motif) row in value-count order; here the hits are aggregated per input sequence, which is
    the documented intent of the score.
    """

    np = require("numpy", _EXTRA)
    seq_transforms = require("grelu.transforms.seq_transforms", _EXTRA)
    seqs = _dna_list(sequences, "sequences")
    if isinstance(motifs, (str, Path)):
        source: Any = str(input_path(motifs))
    else:
        if not isinstance(motifs, Mapping) or not motifs:
            raise ConfigurationError("motifs must be a MEME file path or a non-empty mapping.")
        source = {str(k): np.asarray(v, dtype="float64") for k, v in motifs.items()}
        for key, value in source.items():
            if value.ndim != 2 or value.shape[0] != 4:
                raise ConfigurationError(f"motif {key} must be a 4 x L matrix.")
    scorer, _ = _run(
        "MotifScore",
        seq_transforms.MotifScore,
        source,
        names=None if names is None else [str(n) for n in names],
        weights=None if weights is None else [float(w) for w in weights],
        pthresh=pthresh,
        rc=_flag(rc, "rc"),
    )
    sites, _ = _run(
        "scan_sequences",
        seq_transforms.scan_sequences,
        seqs,
        motifs=scorer.motifs,
        pthresh=pthresh,
        rc=scorer.rc,
    )
    scores = [0.0] * len(seqs)
    hits = [0] * len(seqs)
    if len(sites):
        for row in sites.to_dict(orient="records"):
            index = int(row["seq_idx"])
            weight = 1.0 if scorer.weights is None else float(scorer.weights[row["motif"]])
            scores[index] += weight
            hits[index] += 1
    return {
        "status": "completed",
        "operation": "score_sequence_motifs",
        "hits": hits,
        "scores": scores,
        "provenance": _provenance("MotifScore"),
    }


def validate_genomic_file(
    path: str | Path,
    kind: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that a path is a valid bigWig or FASTA file by gReLU's rules (``check_bigwig`` / ``check_fasta``).

    ``kind`` is ``bigwig`` (extension ``.bw``, ``.bigwig`` or ``.bigWig`` and an existing file) or ``fasta``
    (``.fa``, ``.fasta``, ``.fa.gz``, ``.fasta.gz``, ``.fa.bgz``, ``.fasta.bgz`` and an existing file). The
    check is by name and existence, as in gReLU; it does not parse the content.
    """

    which = _choice(kind, "kind", ("bigwig", "fasta"))
    text = str(Path(_text(str(path), "path")).expanduser())
    if which == "bigwig":
        module = require("grelu.io.bigwig", _EXTRA)
        valid, _ = _run("check_bigwig", module.check_bigwig, text)
    else:
        module = require("grelu.io.fasta", _EXTRA)
        valid, _ = _run("check_fasta", module.check_fasta, text)
    return {
        "status": "completed",
        "operation": "validate_genomic_file",
        "kind": which,
        "path": text,
        "valid": bool(valid),
        "provenance": _provenance(f"check_{which}"),
    }


def get_grelu_resource_path(
    kind: str,
    name: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate a data file bundled with gReLU (``get_meme_file_path`` / ``get_blacklist_file``).

    ``kind="motifs"`` with ``name`` ``hocomoco_v12``, ``hocomoco_v13`` or ``consensus`` (JASPAR 2024
    consensus motifs) returns the MEME file; ``kind="blacklist"`` with ``name`` ``hg19``, ``hg38`` or
    ``mm10`` returns the ENCODE blacklist BED file.
    """

    which = _choice(kind, "kind", ("motifs", "blacklist"))
    resources = require("grelu.resources.utils", _EXTRA)
    if which == "motifs":
        label = _choice(name, "name", ("hocomoco_v12", "hocomoco_v13", "consensus"))
        path, _ = _run("get_meme_file_path", resources.get_meme_file_path, label)
    else:
        label = _choice(name, "name", ("hg19", "hg38", "mm10"))
        path, _ = _run("get_blacklist_file", resources.get_blacklist_file, label)
    location = Path(str(path))
    return {
        "status": "completed",
        "operation": "get_grelu_resource_path",
        "kind": which,
        "name": label,
        "path": str(location),
        "exists": location.is_file(),
        "size_bytes": location.stat().st_size if location.is_file() else None,
        "provenance": _provenance("resources.utils"),
    }


__all__ = [
    "augment_sequences",
    "build_labeled_sequence_dataset",
    "build_sequence_model",
    "compute_prediction_metric",
    "embed_sequences",
    "evaluate_sequence_model",
    "filter_cell_types",
    "get_grelu_resource_path",
    "predict_sequence_model",
    "replace_model_head",
    "scan_motif_insertions",
    "score_prediction_specificity",
    "score_sequence_motifs",
    "score_sequence_patterns",
    "train_sequence_model",
    "transform_labels",
    "validate_genomic_file",
]
_register_public_verbose(globals())
