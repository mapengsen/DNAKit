"""Tests for dnakit.predictions.genomic_benchmarks_gap_functions.

Downloads are replaced: the PyTorch getters, the loc2seq download and the cloud-cache transport are
fakes on the genomic_benchmarks module attributes. One test runs the real PyTorch dataset class on a
layout written under tmp_path, with its "is downloaded" check patched. A socket guard makes any
accidental connection fail. Guard tests run without genomic_benchmarks installed.
"""

from __future__ import annotations

import importlib
import inspect
import os
import socket
import sys
import urllib.request
import zipfile
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.predictions import genomic_benchmarks_gap_functions as gap

LIVE = os.environ.get("DNAKIT_LIVE_TESTS") == "1"
PYTORCH = "genomic_benchmarks.dataset_getters.pytorch_datasets"
LOC2SEQ = "genomic_benchmarks.loc2seq.loc2seq"
CLOUD = "genomic_benchmarks.loc2seq.cloud_caching"
ZENODO = "https://zenodo.org/records/16605299"
# (split, label, sequence): one record per file, named 0.txt in each label folder.
LAYOUT = [
    ("train", "human", "ACGTAC"),
    ("train", "worm", "GG"),
    ("test", "worm", "TTTT"),
    ("test", "human", "AC"),
]


@pytest.fixture(autouse=True)
def _block_network(monkeypatch, request):
    if LIVE and request.node.name.startswith("test_live_"):
        return

    original = socket.socket

    class BlockedSocket(original):
        """Subclassable so that third-party imports still work; creating a socket fails."""

        def __init__(self, *args, **kwargs):
            raise RuntimeError("network access is blocked in unit tests")

    monkeypatch.setattr(socket, "socket", BlockedSocket)


def _module(name: str):
    pytest.importorskip("genomic_benchmarks")
    pytest.importorskip("torch")
    return importlib.import_module(name)


def _write_layout(root: Path, layout=LAYOUT) -> None:
    for split, label, sequence in layout:
        folder = root / split / label
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "0.txt").write_text(sequence)


class FakeDataset:
    """Stands in for a PyTorch dataset: ``len`` and ``[index]`` returning ``(text, label)``."""

    def __init__(self, items):
        self._items = items

    def __len__(self):
        return len(self._items)

    def __getitem__(self, index):
        return self._items[index]


# ---------- guards (no genomic_benchmarks needed)


def test_missing_pytorch_module_raises_backend_unavailable(monkeypatch):
    # Without torch itself the getters report torch first; this checks the next guard.
    pytest.importorskip("torch")
    monkeypatch.setitem(sys.modules, PYTORCH, None)
    with pytest.raises(BackendUnavailableError, match="bio-genomic-benchmarks"):
        gap.load_demo_human_or_worm("train", verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-genomic-benchmarks"):
        gap.get_genomic_benchmark_dataset("demo_human_or_worm", "train", verbose=False)


def test_missing_torch_raises_backend_unavailable_with_its_own_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    with pytest.raises(BackendUnavailableError, match="python -m pip install torch"):
        gap.load_demo_human_or_worm("train", verbose=False)


def test_missing_loc2seq_module_raises_backend_unavailable(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, LOC2SEQ, None)
    with pytest.raises(BackendUnavailableError, match="bio-genomic-benchmarks"):
        gap.download_genomic_benchmark_sequences(
            "demo_human_or_worm", dest_path=str(tmp_path), cache_path=str(tmp_path), verbose=False
        )


def test_missing_cloud_cache_module_raises_backend_unavailable(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, CLOUD, None)
    with pytest.raises(BackendUnavailableError, match="bio-genomic-benchmarks"):
        gap.download_genomic_benchmark_cloud_cache(
            "demo_human_or_worm", 0, str(tmp_path / "demo_human_or_worm"), verbose=False
        )


def test_verbose_is_last_keyword_only_parameter():
    for name in gap.__all__:
        parameters = list(inspect.signature(getattr(gap, name)).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name


# ---------- PyTorch getters (fakes on the module)


NAMED = [
    (
        "load_demo_coding_vs_intergenomic_seqs",
        "DemoCodingVsIntergenomicSeqs",
        "demo_coding_vs_intergenomic_seqs",
        True,
    ),
    ("load_demo_human_or_worm", "DemoHumanOrWorm", "demo_human_or_worm", True),
    (
        "load_drosophila_enhancers_stark",
        "DrosophilaEnhancersStark",
        "drosophila_enhancers_stark",
        False,
    ),
    ("load_demo_mouse_enhancers", "DemoMouseEnhancers", "dummy_mouse_enhancers_ensembl", True),
    ("load_human_enhancers_cohn", "HumanEnhancersCohn", "human_enhancers_cohn", True),
    ("load_human_enhancers_ensembl", "HumanEnhancersEnsembl", "human_enhancers_ensembl", True),
    ("load_human_nontata_promoters", "HumanNontataPromoters", "human_nontata_promoters", True),
    ("load_human_ocr_ensembl", "HumanOcrEnsembl", "human_ocr_ensembl", True),
]


@pytest.mark.parametrize("function_name, getter_name, dataset, takes_cloud_switch", NAMED)
def test_named_getter_calls_upstream_and_describes_the_split(
    monkeypatch, tmp_path, function_name, getter_name, dataset, takes_cloud_switch
):
    pytorch = _module(PYTORCH)
    monkeypatch.setattr(pytorch, "CACHE_PATH", tmp_path)
    calls = []

    def fake_getter(split, **kwargs):
        calls.append((split, kwargs))
        return FakeDataset([("ACGTAC", 1), ("GG", 0), ("TTTT", 1)])

    monkeypatch.setattr(pytorch, getter_name, fake_getter)
    result = getattr(gap, function_name)("train", force_download=True, version=0, verbose=False)

    expected = {"force_download": True, "version": 0}
    if takes_cloud_switch:
        expected["use_cloud_cache"] = True
    assert calls == [("train", expected)]
    assert result["service"] == "genomic_benchmarks.dataset_getters.pytorch_datasets"
    assert result["source"] == ZENODO
    description = result["result"]
    assert description["dataset"] == dataset
    assert description["split"] == "train"
    assert description["local_path"] == str(tmp_path / dataset)
    assert description["summary"] == {
        "record_count": 3,
        "split_counts": {"train": 3},
        "label_counts": {"0": 1, "1": 2},
        "first_record": {"split": "train", "label": "1", "length": 6, "sequence": "ACGTAC"},
    }


def test_named_getter_passes_cloud_switch_false_through(monkeypatch, tmp_path):
    pytorch = _module(PYTORCH)
    monkeypatch.setattr(pytorch, "CACHE_PATH", tmp_path)
    seen = {}

    def fake_getter(split, **kwargs):
        seen.update(kwargs)
        return FakeDataset([("AC", 0)])

    monkeypatch.setattr(pytorch, "DemoHumanOrWorm", fake_getter)
    result = gap.load_demo_human_or_worm("test", use_cloud_cache=False, verbose=False)
    assert seen["use_cloud_cache"] is False
    assert (
        result["source"]
        == "http://raw.githubusercontent.com/ML-Bioinfo-CEITEC/genomic_benchmarks/main/datasets/"
    )


def test_generic_getter_takes_dataset_name_and_split(monkeypatch, tmp_path):
    pytorch = _module(PYTORCH)
    monkeypatch.setattr(pytorch, "CACHE_PATH", tmp_path)
    calls = []

    def fake_get_dataset(
        dataset_name, split, force_download=False, version=None, use_cloud_cache=True
    ):
        calls.append((dataset_name, split, force_download, version, use_cloud_cache))
        return FakeDataset([("GATTACA", 2)])

    monkeypatch.setattr(pytorch, "get_dataset", fake_get_dataset)
    result = gap.get_genomic_benchmark_dataset(
        "human_ocr_ensembl", "test", version=0, verbose=False
    )
    assert calls == [("human_ocr_ensembl", "test", False, 0, True)]
    assert result["result"]["dataset"] == "human_ocr_ensembl"
    assert result["result"]["summary"]["first_record"] == {
        "split": "test",
        "label": "2",
        "length": 7,
        "sequence": "GATTACA",
    }
    assert result["source"] == ZENODO


def test_split_and_version_are_checked_before_the_getter_runs(monkeypatch, tmp_path):
    pytorch = _module(PYTORCH)

    def fail(*args, **kwargs):
        raise AssertionError("upstream getter must not run for invalid arguments")

    monkeypatch.setattr(pytorch, "DemoHumanOrWorm", fail)
    monkeypatch.setattr(pytorch, "get_dataset", fail)
    with pytest.raises(ConfigurationError):
        gap.load_demo_human_or_worm("validation", verbose=False)
    with pytest.raises(ConfigurationError):
        gap.load_demo_human_or_worm("train", version=-1, verbose=False)
    with pytest.raises(ConfigurationError):
        gap.get_genomic_benchmark_dataset("demo_human_or_worm", "dev", verbose=False)


def test_real_pytorch_getter_reads_a_local_layout_without_network(monkeypatch, tmp_path):
    """The upstream dataset class runs for real; only the download check and the cache root are
    patched."""

    pytorch = _module(PYTORCH)
    monkeypatch.setattr(pytorch, "CACHE_PATH", tmp_path)
    monkeypatch.setattr(pytorch, "is_downloaded", lambda name, cache_path=None: True)
    _write_layout(tmp_path / "demo_human_or_worm")
    result = gap.load_demo_human_or_worm("train", verbose=False)
    summary = result["result"]["summary"]
    assert summary["record_count"] == 2
    assert summary["split_counts"] == {"train": 2}
    assert sorted(summary["label_counts"].values()) == [1, 1]
    assert summary["first_record"]["sequence"] in {"ACGTAC", "GG"}
    assert result["result"]["local_path"] == str(tmp_path / "demo_human_or_worm")
    test_result = gap.get_genomic_benchmark_dataset("demo_human_or_worm", "test", verbose=False)
    assert test_result["result"]["summary"]["record_count"] == 2
    assert test_result["result"]["summary"]["first_record"]["split"] == "test"


@pytest.mark.skipif(
    not LIVE, reason="live test: set DNAKIT_LIVE_TESTS=1 to download a benchmark dataset"
)
def test_live_cloud_cache_downloads_demo_dataset(tmp_path):
    _module(CLOUD)
    result = gap.download_genomic_benchmark_cloud_cache(
        "demo_human_or_worm", 0, str(tmp_path / "demo_human_or_worm"), verbose=False
    )
    assert result["result"]["summary"]["record_count"] > 0


# ---------------------------------------------------------------- sequence downloads (fake loc2seq)


def test_sequence_download_passes_arguments_and_summarises_the_output(monkeypatch, tmp_path):
    loc2seq = _module(LOC2SEQ)
    calls = []

    def fake_download(
        interval_list_dataset,
        version=None,
        dest_path=None,
        cache_path=None,
        force_download=False,
        use_cloud_cache=True,
        local_repo=False,
    ):
        calls.append(
            {
                "interval_list_dataset": interval_list_dataset,
                "version": version,
                "dest_path": Path(dest_path),
                "cache_path": Path(cache_path),
                "force_download": force_download,
                "use_cloud_cache": use_cloud_cache,
                "local_repo": local_repo,
            }
        )
        root = Path(dest_path) / "demo_human_or_worm"
        _write_layout(root)
        return root

    monkeypatch.setattr(loc2seq, "download_dataset", fake_download)
    result = gap.download_genomic_benchmark_sequences(
        "demo_human_or_worm",
        dest_path=str(tmp_path / "data"),
        cache_path=str(tmp_path / "refs"),
        verbose=False,
    )
    assert calls == [
        {
            "interval_list_dataset": "demo_human_or_worm",
            "version": None,
            "dest_path": tmp_path / "data",
            "cache_path": tmp_path / "refs",
            "force_download": False,
            "use_cloud_cache": True,
            "local_repo": False,
        }
    ]
    assert result["service"] == "genomic_benchmarks.loc2seq"
    assert result["source"] == "demo_human_or_worm"
    description = result["result"]
    assert description["dataset"] == "demo_human_or_worm"
    assert description["local_path"] == str(tmp_path / "data" / "demo_human_or_worm")
    assert description["summary"] == {
        "record_count": 4,
        "split_counts": {"train": 2, "test": 2},
        "label_counts": {"human": 2, "worm": 2},
        "first_record": {"split": "train", "label": "human", "length": 6, "sequence": "ACGTAC"},
    }


def test_sequence_download_uses_package_locations_by_default(monkeypatch, tmp_path):
    loc2seq = _module(LOC2SEQ)
    monkeypatch.setattr(loc2seq, "CACHE_PATH", tmp_path / "home")
    monkeypatch.setattr(loc2seq, "REF_CACHE_PATH", tmp_path / "home" / "fasta")
    seen = {}

    def fake_download(interval_list_dataset, **kwargs):
        seen.update(kwargs)
        root = Path(kwargs["dest_path"]) / "human_nontata_promoters"
        _write_layout(root, [("train", "pos", "ACGT")])
        return root

    monkeypatch.setattr(loc2seq, "download_dataset", fake_download)
    result = gap.download_genomic_benchmark_sequences("human_nontata_promoters", verbose=False)
    assert seen["dest_path"] == tmp_path / "home"
    assert seen["cache_path"] == tmp_path / "home" / "fasta"
    assert result["result"]["summary"]["record_count"] == 1


def test_sequence_download_rejects_a_missing_output_directory(monkeypatch, tmp_path):
    loc2seq = _module(LOC2SEQ)
    monkeypatch.setattr(loc2seq, "download_dataset", lambda name, **kwargs: tmp_path / "absent")
    with pytest.raises(BackendExecutionError, match="did not leave a dataset directory"):
        gap.download_genomic_benchmark_sequences(
            "demo_human_or_worm", dest_path=str(tmp_path), cache_path=str(tmp_path), verbose=False
        )


def test_sequence_download_rejects_bad_flags_before_downloading(monkeypatch, tmp_path):
    loc2seq = _module(LOC2SEQ)

    def fail(*args, **kwargs):
        raise AssertionError("upstream download must not run for invalid arguments")

    monkeypatch.setattr(loc2seq, "download_dataset", fail)
    with pytest.raises(ConfigurationError):
        gap.download_genomic_benchmark_sequences("demo_human_or_worm", version=-2, verbose=False)
    with pytest.raises(ConfigurationError):
        gap.download_genomic_benchmark_sequences(
            "demo_human_or_worm", local_repo="yes", verbose=False
        )


# ---------- cloud cache (fake transport, real unpacking)


def _zip_layout(archive: Path, top: str) -> None:
    with zipfile.ZipFile(archive, "w") as handle:
        for split, label, sequence in LAYOUT:
            handle.writestr(f"{top}/{split}/{label}/0.txt", sequence)


def test_cloud_cache_download_unpacks_the_archive_it_fetches(monkeypatch, tmp_path):
    _module(CLOUD)
    requested = []

    def fake_urlretrieve(url, filename, *args, **kwargs):
        requested.append(url)
        _zip_layout(Path(filename), "demo_human_or_worm")
        return filename, None

    monkeypatch.setattr(urllib.request, "urlretrieve", fake_urlretrieve)
    dest = tmp_path / "out" / "demo_human_or_worm"
    result = gap.download_genomic_benchmark_cloud_cache(
        "demo_human_or_worm", 0, str(dest), verbose=False
    )
    assert requested == [f"{ZENODO}/files/demo_human_or_worm_v0.zip?download=1"]
    assert result["service"] == "genomic_benchmarks.loc2seq.cloud_caching"
    assert result["source"] == ZENODO
    assert result["result"]["local_path"] == str(dest.resolve())
    assert result["result"]["summary"]["record_count"] == 4
    assert result["result"]["summary"]["label_counts"] == {"human": 2, "worm": 2}
    assert not (tmp_path / "out" / "demo_human_or_worm.zip").exists()


def test_cloud_cache_keeps_existing_files_unless_forced(monkeypatch, tmp_path):
    _module(CLOUD)

    def fake_urlretrieve(url, filename, *args, **kwargs):
        _zip_layout(Path(filename), "demo_human_or_worm")
        return filename, None

    monkeypatch.setattr(urllib.request, "urlretrieve", fake_urlretrieve)
    dest = tmp_path / "demo_human_or_worm"
    dest.mkdir()
    (dest / "marker.txt").write_text("keep")
    gap.download_genomic_benchmark_cloud_cache("demo_human_or_worm", 0, str(dest), verbose=False)
    assert (dest / "marker.txt").read_text() == "keep"
    gap.download_genomic_benchmark_cloud_cache(
        "demo_human_or_worm", 0, str(dest), force_download=True, verbose=False
    )
    assert not (dest / "marker.txt").exists()
    assert (dest / "train" / "human" / "0.txt").exists()


def test_cloud_cache_rejects_names_outside_the_cache_before_connecting(monkeypatch, tmp_path):
    _module(CLOUD)

    def fail(url, filename, *args, **kwargs):
        raise AssertionError("no transfer may start for a name outside the cloud cache")

    monkeypatch.setattr(urllib.request, "urlretrieve", fail)
    with pytest.raises(ConfigurationError, match="not in the genomic_benchmarks cloud cache"):
        gap.download_genomic_benchmark_cloud_cache(
            "not_a_dataset", 0, str(tmp_path / "x"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.download_genomic_benchmark_cloud_cache(
            "demo_human_or_worm", "0", str(tmp_path / "x"), verbose=False
        )
