"""genomic_benchmarks wrappers: named benchmark datasets, sequence downloads and the Zenodo cloud
cache.

The package downloads benchmark data from its Zenodo cloud cache (or from the GitHub datasets
repository) and stores it under ``~/.genomic_benchmarks`` unless a destination is given. The
functions here are network wrappers: each one calls the package's own download or getter function
through the module attribute, so tests can replace it. Each result is JSON-serialisable and holds
the
local directory, the record count, the label counts and the first record. The dataset objects that
the PyTorch getters build are not returned, because they are not JSON data.

genomic_benchmarks, torch and numpy are imported only when a function runs. Nothing else is fetched
from the network.
"""

from __future__ import annotations

import importlib
from collections import Counter
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "bio-genomic-benchmarks"
_PYTORCH = "genomic_benchmarks.dataset_getters.pytorch_datasets"
_LOC2SEQ = "genomic_benchmarks.loc2seq.loc2seq"
_CLOUD = "genomic_benchmarks.loc2seq.cloud_caching"
_PATHS = "genomic_benchmarks.utils.paths"
_SPLITS = ("train", "test")
_SERVICE_PYTORCH = "genomic_benchmarks.dataset_getters.pytorch_datasets"
_SERVICE_LOC2SEQ = "genomic_benchmarks.loc2seq"
_SERVICE_CLOUD = "genomic_benchmarks.loc2seq.cloud_caching"


def _pytorch_module() -> Any:
    """Import the PyTorch getters. genomic-benchmarks does not install torch, so it is checked
    first."""

    try:
        importlib.import_module("torch")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"torch could not be imported: {exc}",
            hint="The PyTorch dataset getters need PyTorch, which the genomic-benchmarks extra "
            "does not install: python -m pip install torch",
        ) from exc
    return require(_PYTORCH, _EXTRA)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _split(value: object) -> str:
    if not isinstance(value, str) or value not in _SPLITS:
        raise ConfigurationError("split must be 'train' or 'test'.")
    return value


def _version(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError("version must be a non-negative integer or null.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _directory_summary(path: Path) -> dict[str, Any]:
    """Summarise a full-sequence dataset directory laid out as <split>/<label>/<id>.txt."""

    if not path.is_dir():
        raise BackendExecutionError(
            f"genomic_benchmarks did not leave a dataset directory at {path}."
        )
    split_counts: dict[str, int] = {}
    labels: Counter[str] = Counter()
    first: dict[str, Any] | None = None
    for split in _SPLITS:
        split_dir = path / split
        if not split_dir.is_dir():
            continue
        total = 0
        for label_dir in sorted(item for item in split_dir.iterdir() if item.is_dir()):
            files = sorted(label_dir.glob("*.txt"))
            total += len(files)
            labels[label_dir.name] += len(files)
            if first is None and files:
                sequence = files[0].read_text()
                first = {
                    "split": split,
                    "label": label_dir.name,
                    "length": len(sequence),
                    "sequence": sequence,
                }
        split_counts[split] = total
    return {
        "record_count": sum(split_counts.values()),
        "split_counts": split_counts,
        "label_counts": dict(sorted(labels.items())),
        "first_record": first,
    }


def _dataset_summary(dataset: Any, split: str) -> dict[str, Any]:
    """Summarise a loaded PyTorch dataset (``len``, ``[index]`` returning ``(text, label)``)."""

    count = len(dataset)
    all_labels = getattr(dataset, "all_labels", None)
    labels = all_labels if all_labels is not None else [dataset[index][1] for index in range(count)]
    counts = Counter(str(label) for label in labels)
    first: dict[str, Any] | None = None
    if count:
        text, label = dataset[0]
        first = {"split": split, "label": str(label), "length": len(text), "sequence": text}
    return {
        "record_count": count,
        "split_counts": {split: count},
        "label_counts": dict(sorted(counts.items())),
        "first_record": first,
    }


def _pytorch_source(dataset_name: str, use_cloud_cache: bool, version: int | None) -> str:
    cloud = require(_CLOUD, _EXTRA)
    paths = require(_PATHS, _EXTRA)
    key = (dataset_name, 0 if version is None else version)
    if use_cloud_cache and key in cloud.CLOUD_CACHE:
        return str(cloud.CLOUD_CACHE[key])
    return str(paths.DATASET_URL_PATH)


def _pytorch_result(
    dataset_name: str,
    split: str,
    dataset: Any,
    *,
    use_cloud_cache: bool,
    version: int | None,
) -> dict[str, Any]:
    pytorch = _pytorch_module()
    return {
        "service": _SERVICE_PYTORCH,
        "source": _pytorch_source(dataset_name, use_cloud_cache, version),
        "result": {
            "dataset": dataset_name,
            "split": split,
            "local_path": str(Path(pytorch.CACHE_PATH) / dataset_name),
            "summary": _dataset_summary(dataset, split),
        },
    }


def _pytorch_options(
    force_download: object, version: object, use_cloud_cache: object
) -> dict[str, Any]:
    return {
        "force_download": _flag(force_download, "force_download"),
        "version": _version(version),
        "use_cloud_cache": _flag(use_cloud_cache, "use_cloud_cache"),
    }


def get_genomic_benchmark_dataset(
    dataset_name: str,
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a named genomic_benchmarks classification dataset for one split and describe it.

    This is a network call when the dataset is not yet on disk: the package downloads it into
    ``~/.genomic_benchmarks/<dataset_name>`` (from the Zenodo cloud cache when ``use_cloud_cache``).
    ``split`` is ``train`` or ``test``. The result holds the record count, label counts and the
    first
    record of that split; the dataset object itself is not returned.
    """

    pytorch = _pytorch_module()
    name = _text(dataset_name, "dataset_name")
    chosen = _split(split)
    options = _pytorch_options(force_download, version, use_cloud_cache)
    dataset = pytorch.get_dataset(name, chosen, **options)
    return _pytorch_result(
        name,
        chosen,
        dataset,
        use_cloud_cache=options["use_cloud_cache"],
        version=options["version"],
    )


def load_demo_coding_vs_intergenomic_seqs(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the demo coding versus intergenomic sequence dataset (network if not on disk)."""

    return _named_pytorch(
        "DemoCodingVsIntergenomicSeqs",
        "demo_coding_vs_intergenomic_seqs",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_demo_human_or_worm(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the demo human-or-worm sequence dataset (network if not on disk)."""

    return _named_pytorch(
        "DemoHumanOrWorm",
        "demo_human_or_worm",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_drosophila_enhancers_stark(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the Drosophila enhancer dataset from Stark et al. (network if not on disk).

    The upstream getter has no cloud-cache switch and always uses the cloud cache when it is
    possible.
    """

    return _named_pytorch(
        "DrosophilaEnhancersStark",
        "drosophila_enhancers_stark",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=True,
        pass_cloud_switch=False,
    )


def load_demo_mouse_enhancers(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the demo mouse enhancer dataset (stored upstream as dummy_mouse_enhancers_ensembl)."""

    return _named_pytorch(
        "DemoMouseEnhancers",
        "dummy_mouse_enhancers_ensembl",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_human_enhancers_cohn(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the human enhancer dataset of Cohn et al. (network if not on disk)."""

    return _named_pytorch(
        "HumanEnhancersCohn",
        "human_enhancers_cohn",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_human_enhancers_ensembl(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the human enhancer dataset built from Ensembl regulatory regions (network if not on
    disk)."""

    return _named_pytorch(
        "HumanEnhancersEnsembl",
        "human_enhancers_ensembl",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_human_nontata_promoters(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the human TATA-less promoter dataset (network if not on disk)."""

    return _named_pytorch(
        "HumanNontataPromoters",
        "human_nontata_promoters",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def load_human_ocr_ensembl(
    split: str,
    *,
    force_download: bool = False,
    version: int | None = None,
    use_cloud_cache: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the human open-chromatin region dataset built from Ensembl (network if not on disk)."""

    return _named_pytorch(
        "HumanOcrEnsembl",
        "human_ocr_ensembl",
        split,
        force_download=force_download,
        version=version,
        use_cloud_cache=use_cloud_cache,
    )


def _named_pytorch(
    getter_name: str,
    dataset_name: str,
    split: str,
    *,
    force_download: bool,
    version: int | None,
    use_cloud_cache: bool,
    pass_cloud_switch: bool = True,
) -> dict[str, Any]:
    pytorch = _pytorch_module()
    chosen = _split(split)
    options = _pytorch_options(force_download, version, use_cloud_cache)
    if not pass_cloud_switch:
        del options["use_cloud_cache"]
    getter = getattr(pytorch, getter_name)
    dataset = getter(chosen, **options)
    return _pytorch_result(
        dataset_name,
        chosen,
        dataset,
        use_cloud_cache=bool(options.get("use_cloud_cache", True)),
        version=options["version"],
    )


def download_genomic_benchmark_sequences(
    interval_list_dataset: str,
    *,
    version: int | None = None,
    dest_path: str | None = None,
    cache_path: str | None = None,
    force_download: bool = False,
    use_cloud_cache: bool = True,
    local_repo: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the full-sequence version of an interval-list dataset and describe it (network call).

    ``interval_list_dataset`` is a benchmark name (for example ``demo_human_or_worm``) or, with
    ``local_repo=True``, a path to a local dataset folder with a ``metadata.yaml``. The upstream
    function reads the dataset metadata and the reference genome FASTA files (downloaded into
    ``cache_path`` when missing) and writes ``<dest_path>/<name>/<split>/<label>/<id>.txt``. Without
    ``use_cloud_cache`` it always rebuilds from the references; with it, a pre-built cloud copy is
    used
    when one exists. ``dest_path`` and ``cache_path`` default to the package's locations under the
    home directory.
    """

    loc2seq = require(_LOC2SEQ, _EXTRA)
    name = _text(interval_list_dataset, "interval_list_dataset")
    destination = (
        Path(dest_path).expanduser() if dest_path is not None else Path(loc2seq.CACHE_PATH)
    )
    references = (
        Path(cache_path).expanduser() if cache_path is not None else Path(loc2seq.REF_CACHE_PATH)
    )
    target = Path(
        loc2seq.download_dataset(
            name,
            version=_version(version),
            dest_path=destination,
            cache_path=references,
            force_download=_flag(force_download, "force_download"),
            use_cloud_cache=_flag(use_cloud_cache, "use_cloud_cache"),
            local_repo=_flag(local_repo, "local_repo"),
        )
    )
    return {
        "service": _SERVICE_LOC2SEQ,
        "source": name,
        "result": {
            "dataset": target.name,
            "local_path": str(target),
            "summary": _directory_summary(target),
        },
    }


def download_genomic_benchmark_cloud_cache(
    dataset_name: str,
    version: int,
    dest_path: str,
    *,
    force_download: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a pre-built full-sequence dataset from the Zenodo cloud cache and describe it
    (network call).

    Only the (name, version) pairs listed in the package's cloud cache are allowed; other names are
    rejected before any connection is made. ``dest_path`` is the dataset directory to create (the
    zip
    archive is unpacked into its parent directory). With ``force_download`` the existing directory
    is
    removed first (upstream behaviour); the default keeps existing files.
    """

    cloud = require(_CLOUD, _EXTRA)
    name = _text(dataset_name, "dataset_name")
    release = _version(version)
    if release is None:
        raise ConfigurationError("version must be a non-negative integer.")
    key = (name, release)
    if key not in cloud.CLOUD_CACHE:
        raise ConfigurationError(
            f"{name} version {release} is not in the genomic_benchmarks cloud cache."
        )
    destination = Path(_text(dest_path, "dest_path")).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    target = Path(
        cloud.download_from_cloud_cache(
            key,
            destination,
            force_download=_flag(force_download, "force_download"),
        )
    )
    return {
        "service": _SERVICE_CLOUD,
        "source": str(cloud.CLOUD_CACHE[key]),
        "result": {
            "dataset": name,
            "local_path": str(target),
            "summary": _directory_summary(target),
        },
    }


__all__ = [
    "download_genomic_benchmark_cloud_cache",
    "download_genomic_benchmark_sequences",
    "get_genomic_benchmark_dataset",
    "load_demo_coding_vs_intergenomic_seqs",
    "load_demo_human_or_worm",
    "load_demo_mouse_enhancers",
    "load_drosophila_enhancers_stark",
    "load_human_enhancers_cohn",
    "load_human_enhancers_ensembl",
    "load_human_nontata_promoters",
    "load_human_ocr_ensembl",
]

_register_public_verbose(globals(), tuple(__all__))
