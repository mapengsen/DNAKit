"""DNABERT-2 / DNABERT-S embeddings with the LucaOneTasks protocol (offline; needs a local checkpoint)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import dnakit.representations.lucaonetasks_p16_functions as luca
from dnakit.exceptions import ConfigurationError


def _checkpoint() -> str:
    explicit = os.environ.get("DNAKIT_DNABERT_CHECKPOINT")
    if explicit and Path(explicit).is_dir():
        return explicit
    pytest.importorskip("torch")
    pytest.importorskip("transformers")
    from huggingface_hub import snapshot_download

    try:
        return snapshot_download("zhihan1996/DNABERT-2-117M", local_files_only=True)
    except Exception:
        pytest.skip("no local DNABERT checkpoint (set DNAKIT_DNABERT_CHECKPOINT)")


def test_requires_explicit_remote_code_opt_in() -> None:
    with pytest.raises(ConfigurationError, match="Python code"):
        luca.embed_dna_with_dnabert(["ACGT"], verbose=False)
    with pytest.raises(ConfigurationError, match="model"):
        luca.embed_dna_with_dnabert(["ACGT"], model="bert", allow_remote_code=True, verbose=False)


def test_matrix_vector_and_mean(tmp_path: Path) -> None:
    pytest.importorskip("numpy")
    checkpoint = _checkpoint()
    sequences = {
        "a": "ACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGA",
        "b": "ttttggggccccaaaattttggggccccaaaa",
    }
    common = {
        "model": "dnabert2",
        "checkpoint_path": checkpoint,
        "allow_remote_code": True,
        "device": "cpu",
        "verbose": False,
    }
    matrix = luca.embed_dna_with_dnabert(sequences, output_path=str(tmp_path / "m.npz"), **common)
    with_special = luca.embed_dna_with_dnabert(sequences, include_special_tokens=True, **common)
    vector = luca.embed_dna_with_dnabert(sequences, embedding_type="vector", **common)
    mean = luca.embed_dna_with_dnabert(sequences, embedding_type="mean", **common)
    rows = {s["id"]: s["shape"][0] for s in matrix["sequences"]}
    assert all(s["shape"][0] == rows[s["id"]] + 2 for s in with_special["sequences"])
    assert len(vector["embeddings"]["a"]) == 768 and len(mean["embeddings"]["b"]) == 768
    import numpy as np

    saved = np.load(tmp_path / "m.npz")
    assert saved["a"].shape == (rows["a"], 768) and not np.allclose(
        vector["embeddings"]["a"], mean["embeddings"]["a"]
    )
