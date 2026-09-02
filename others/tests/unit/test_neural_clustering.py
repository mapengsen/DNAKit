"""Tests for checkpoint reuse, rep extraction, and neural k-means clustering."""

from __future__ import annotations

import json
import sys
import types
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from dnakit.core import DNARecord, DNASequence, Gap
from dnakit.datasets import (
    NeuralClusteringConfig,
    neural_cluster_sequences,
)
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    SequenceError,
    UnsupportedGapOperationError,
)
from dnakit.representations import (
    RepresentationConfig,
    available_embedding_models,
    ensure_model_checkpoint,
    extract_representations,
    get_embedding_model,
)


def _record(record_id: str, sequence: str, *, alphabet: str = "strict") -> DNARecord:
    return DNARecord(DNASequence(sequence, alphabet=alphabet), record_id)


class _FakeBackend:
    def __init__(self, matrix: list[list[float]]) -> None:
        self.matrix = matrix
        self.sequences: tuple[str, ...] = ()
        self.show_progress: bool | None = None

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        self.sequences = tuple(sequences)
        self.show_progress = show_progress
        return np.asarray(self.matrix, dtype=np.float32)


def test_model_registry_covers_compared_model_families_and_aliases() -> None:
    assert available_embedding_models() == (
        "alphagenome",
        "caduceus",
        "dnabert2",
        "enformer",
        "evo2",
        "generator",
        "grover",
        "hyenadna",
        "janusdna",
        "lucaone",
        "ntv2",
    )
    assert get_embedding_model("DNABERT-2").name == "dnabert2"
    assert get_embedding_model("nt-v2").name == "ntv2"
    assert get_embedding_model("alphagenome").allow_patterns is None
    assert all(
        model.source_repository.startswith("https://")
        and model.checkpoint_url.startswith("https://")
        for model in map(get_embedding_model, available_embedding_models())
    )
    with pytest.raises(ConfigurationError) as error:
        get_embedding_model("unknown")
    assert error.value.code == "INVALID_EMBEDDING_MODEL"


def test_remote_checkpoint_code_requires_explicit_opt_in() -> None:
    assert RepresentationConfig().model == "lucaone"
    with pytest.raises(ConfigurationError) as default_error:
        extract_representations(
            [_record("a", "AAAA")],
            config=RepresentationConfig(show_progress=False),
        )
    assert default_error.value.code == "MODEL_REMOTE_CODE_NOT_ALLOWED"

    config = RepresentationConfig(model="dnabert2", show_progress=False)
    with pytest.raises(ConfigurationError) as error:
        extract_representations([_record("a", "AAAA")], config=config)
    assert error.value.code == "MODEL_REMOTE_CODE_NOT_ALLOWED"
    assert RepresentationConfig(model="dnabert2", allow_remote_code=True).model == "dnabert2"

    custom = extract_representations(
        [_record("a", "AAAA")],
        config=config,
        backend=_FakeBackend([[1.0]]),
    )
    assert custom.model_name == "dnabert2"


def test_existing_checkpoint_in_working_directory_is_reused(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    target = tmp_path / "ckpt" / "dnabert2-117m"
    target.mkdir(parents=True)
    (target / "config.json").write_text("{}", encoding="utf-8")
    (target / "pytorch_model.bin").write_bytes(b"existing")

    first = ensure_model_checkpoint("dnabert2", show_progress=False)
    second = ensure_model_checkpoint("DNABERT-2", show_progress=False)

    assert first.path == str(target)
    assert not first.downloaded
    assert second == first
    manifest = json.loads((target / ".dnakit-checkpoint.json").read_text(encoding="utf-8"))
    assert manifest["checkpoint_id"] == "zhihan1996/DNABERT-2-117M"


def test_default_lucaone_checkpoint_downloads_to_requested_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[dict[str, object]] = []
    fake_hub = types.ModuleType("huggingface_hub")

    def fake_snapshot_download(**kwargs: object) -> str:
        calls.append(dict(kwargs))
        local_dir = Path(str(kwargs["local_dir"]))
        local_dir.mkdir(parents=True, exist_ok=True)
        (local_dir / "config.json").write_text("{}", encoding="utf-8")
        (local_dir / "model.safetensors").write_bytes(b"weights")
        return str(local_dir)

    fake_hub.snapshot_download = fake_snapshot_download  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hub)

    result = ensure_model_checkpoint(
        RepresentationConfig().model,
        checkpoint_dir=tmp_path / "models",
        show_progress=True,
    )

    assert result.downloaded
    assert Path(result.path) == tmp_path / "models" / "lucaone-gene-step36-8m"
    assert calls[0]["repo_id"] == "LucaGroup/LucaOne-gene-step36.8M"
    assert "*.safetensors" in calls[0]["allow_patterns"]  # type: ignore[operator]


def test_representation_extraction_normalizes_iupac_and_returns_read_only_matrix() -> None:
    backend = _FakeBackend([[1.0, 2.0, 3.0]])
    result = extract_representations(
        [_record("iupac", "ARYN", alphabet="iupac")],
        config=RepresentationConfig(show_progress=False),
        backend=backend,
    )

    assert backend.sequences == ("ANNN",)
    assert backend.show_progress is False
    assert result.representations.shape == (1, 3)
    assert not result.representations.flags.writeable
    assert result.checkpoint_path is None
    assert result.to_dict()["representations"] == [[1.0, 2.0, 3.0]]


def test_representation_extraction_rejects_ambiguity_gaps_and_invalid_shape() -> None:
    ambiguous = _record("ambiguous", "AN", alphabet="iupac")
    with pytest.raises(SequenceError) as ambiguity_error:
        extract_representations(
            [ambiguous],
            config=RepresentationConfig(ambiguity_policy="error", show_progress=False),
            backend=_FakeBackend([[1.0]]),
        )
    assert ambiguity_error.value.code == "MODEL_AMBIGUOUS_INPUT"

    gapped = DNARecord(
        DNASequence.from_fragments(["AA", "TT"], [Gap(3)]),
        "gapped",
    )
    with pytest.raises(UnsupportedGapOperationError) as gap_error:
        extract_representations(
            [gapped],
            config=RepresentationConfig(show_progress=False),
            backend=_FakeBackend([[1.0]]),
        )
    assert gap_error.value.code == "MODEL_GAPPED_INPUT"

    with pytest.raises(BackendExecutionError) as matrix_error:
        extract_representations(
            [_record("a", "AAAA"), _record("b", "CCCC")],
            config=RepresentationConfig(show_progress=False),
            backend=_FakeBackend([[1.0, 2.0]]),
        )
    assert matrix_error.value.code == "INVALID_REPRESENTATION_MATRIX"


def test_representation_extraction_rejects_checkpoint_file(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError) as config_error:
        RepresentationConfig(checkpoint_dir=tmp_path, checkpoint_path=tmp_path)
    assert config_error.value.code == "INVALID_REPRESENTATION_CONFIG"

    checkpoint_file = tmp_path / "weights.bin"
    checkpoint_file.write_bytes(b"not-a-directory")
    with pytest.raises(ConfigurationError) as error:
        extract_representations(
            [_record("a", "AAAA")],
            config=RepresentationConfig(
                model="grover",
                checkpoint_path=checkpoint_file,
                show_progress=False,
            ),
        )
    assert error.value.code == "MODEL_CHECKPOINT_NOT_FOUND"


def test_neural_kmeans_clusters_rep_and_selects_centroid_representatives() -> None:
    pytest.importorskip("sklearn")
    records = [
        _record("a", "AAAA"),
        _record("b", "AAAT"),
        _record("c", "CCCC"),
        _record("d", "CCCG"),
    ]
    backend = _FakeBackend(
        [
            [1.0, 0.0, 0.0],
            [0.9, 0.1, 0.0],
            [0.0, 1.0, 0.0],
            [0.1, 0.9, 0.0],
        ]
    )
    config = NeuralClusteringConfig(
        representation=RepresentationConfig(show_progress=False),
        n_clusters=2,
        seed=13,
    )

    first = neural_cluster_sequences(records, config=config, backend=backend)
    second = neural_cluster_sequences(records, config=config, backend=backend)

    assert first.labels == second.labels
    assert first.representatives.ids == second.representatives.ids
    assert first.inertia == pytest.approx(second.inertia, rel=1e-12, abs=1e-15)
    assert first.labels == (0, 0, 1, 1)
    assert tuple(cluster.member_ids for cluster in first.clusters) == (("a", "b"), ("c", "d"))
    assert first.representatives.ids == ("a", "c")
    assert first.silhouette_score is not None and first.silhouette_score > 0.8
    assert first.embedding_dimension == 3
    assert first.clustering_dimension == 3
    assert json.loads(json.dumps(first.to_dict()))["model_name"] == "lucaone"


def test_neural_kmeans_supports_pca_and_validates_cluster_count() -> None:
    pytest.importorskip("sklearn")
    records = [_record("a", "AAAA"), _record("b", "CCCC"), _record("c", "GGGG")]
    backend = _FakeBackend([[1.0, 0.0, 0.0], [0.9, 0.1, 0.0], [0.0, 1.0, 1.0]])
    result = neural_cluster_sequences(
        records,
        config=NeuralClusteringConfig(
            representation=RepresentationConfig(show_progress=False),
            n_clusters=2,
            pca_components=2,
            seed=7,
        ),
        backend=backend,
    )
    assert result.clustering_dimension == 2
    assert len(result.pca_explained_variance_ratio) == 2

    with pytest.raises(ConfigurationError) as cluster_error:
        neural_cluster_sequences(
            records,
            config=NeuralClusteringConfig(
                representation=RepresentationConfig(show_progress=False),
                n_clusters=4,
            ),
            backend=backend,
        )
    assert cluster_error.value.code == "INVALID_NEURAL_CLUSTER_COUNT"

    with pytest.raises(ConfigurationError) as pca_error:
        neural_cluster_sequences(
            records[:1],
            config=NeuralClusteringConfig(
                representation=RepresentationConfig(show_progress=False),
                n_clusters=1,
                pca_components=1,
            ),
            backend=_FakeBackend([[1.0, 0.0, 0.0]]),
        )
    assert pca_error.value.code == "INVALID_NEURAL_CLUSTER_PCA"
