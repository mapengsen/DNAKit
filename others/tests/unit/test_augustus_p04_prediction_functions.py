"""Tests for dnakit.annotation.augustus_p04_prediction_functions (AUGUSTUS 3.5.0 facades).

Argument-validation tests run everywhere; tests that run AUGUSTUS skip when ``augustus`` is not on PATH.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from dnakit.annotation import augustus_p04_prediction_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENOME = DATA / "genome.fa"
HAVE = shutil.which("augustus") is not None and shutil.which("etraining") is not None
requires_augustus = pytest.mark.skipif(
    not HAVE, reason="AUGUSTUS (conda install -c bioconda augustus) is not on PATH"
)


def _config_source() -> Path | None:
    value = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if value and Path(value).is_dir():
        return Path(value)
    exe = shutil.which("augustus")
    if exe:
        sibling = Path(exe).resolve().parent.parent / "config"
        if sibling.is_dir():
            return sibling
    return None


@pytest.fixture()
def mini_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A small writable config folder (templates plus the human parameter set)."""
    source = _config_source()
    if source is None:
        pytest.skip("no AUGUSTUS config folder")
    target = tmp_path / "config"
    (target / "species").mkdir(parents=True)
    for name in ("cgp", "extrinsic", "model", "parameters", "profile"):
        shutil.copytree(source / name, target / name)
    for name in ("generic", "human"):
        shutil.copytree(source / "species" / name, target / "species" / name)
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(target))
    return target


def test_missing_program_raises_install_hint(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        m.augustus_predict_genes(
            str(GENOME),
            "human",
            str(tmp_path / "o.gtf"),
            augustus_bin_dir=str(tmp_path),
            verbose=False,
        )


def test_validation_errors(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        m.augustus_predict_genes(str(GENOME), "--bad", str(tmp_path / "o"), verbose=False)
    with pytest.raises(ConfigurationError):
        m.augustus_predict_genes(str(GENOME), "human", str(tmp_path / "o"), sample=5, verbose=False)
    with pytest.raises(ConfigurationError):
        m.augustus_predict_genes(
            str(GENOME), "human", str(tmp_path / "o"), strand="left", verbose=False
        )
    existing = tmp_path / "exists"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        m.augustus_predict_genes(str(GENOME), "human", str(existing), verbose=False)
    with pytest.raises(ConfigurationError):
        m.augustus_choose_parameter_set(str(GENOME), ["human"], mode="metagenomic", verbose=False)
    with pytest.raises(ConfigurationError):
        m.augustus_join_segment_predictions([], str(tmp_path / "j"), verbose=False)


@requires_augustus
def test_predict_genes_modes(tmp_path: Path, mini_config: Path) -> None:
    plain = m.augustus_predict_genes(
        str(GENOME), "human", str(tmp_path / "plain.gtf"), verbose=False
    )
    assert plain["status"] == "completed" and plain["version"] and plain["genes"] >= 1
    assert plain["transcripts"] == plain["genes"]
    mea = m.augustus_predict_genes(
        str(GENOME),
        "human",
        str(tmp_path / "mea.gtf"),
        mea=True,
        protein_sequences=True,
        verbose=False,
    )
    assert (
        "--mea=1" in mea["command"]
        and "# protein sequence" in Path(mea["outputs"]["predictions"]).read_text()
    )
    gff3 = m.augustus_predict_genes(
        str(GENOME), "human", str(tmp_path / "p.gff3"), gff3=True, verbose=False
    )
    assert "ID=g1" in Path(gff3["outputs"]["predictions"]).read_text() and gff3["transcripts"] >= 1
    sampled = m.augustus_predict_genes(
        str(GENOME),
        "human",
        str(tmp_path / "s.gtf"),
        alternatives_from_sampling=True,
        sample=3,
        verbose=False,
    )
    assert "--sample=3" in sampled["command"]


@requires_augustus
def test_join_functions(tmp_path: Path, mini_config: Path) -> None:
    first = m.augustus_predict_genes(str(GENOME), "human", str(tmp_path / "a.gtf"), verbose=False)
    joined = m.augustus_join_segment_predictions(
        [first["outputs"]["predictions"]] * 2, str(tmp_path / "j.gtf"), verbose=False
    )
    assert joined["input_files"] == 2 and 1 <= joined["genes"] <= 2 * first["genes"]
    merged = m.augustus_join_gene_sets(
        [first["outputs"]["predictions"]], str(tmp_path / "m.gtf"), verbose=False
    )
    assert merged["transcripts"] == first["transcripts"]


@requires_augustus
def test_train_optimize_and_parameter_set(tmp_path: Path, mini_config: Path) -> None:
    from dnakit.annotation import augustus_p04_conversion_functions as conv

    gtf = tmp_path / "pred.gtf"
    m.augustus_predict_genes(str(GENOME), "human", str(gtf), verbose=False)
    genbank = tmp_path / "train.gb"
    conv.augustus_gff_to_genbank(str(gtf), str(GENOME), 500, str(genbank), verbose=False)
    trained = m.augustus_train_species(
        "toy", str(genbank), set_stop_codon_frequencies=True, verbose=False
    )
    assert "toy_exon_probs.pbl" in trained["parameter_files"]
    assert abs(sum(trained["stop_codon_frequencies"].values()) - 1.0) < 0.01
    with pytest.raises(ConfigurationError):
        m.augustus_train_species("toy", str(genbank), method="other", verbose=False)
    if os.environ.get("DNAKIT_AUGUSTUS_SLOW") == "1":  # about five minutes
        optimized = m.augustus_optimize_parameters(
            "toy",
            str(genbank),
            str(tmp_path / "opt"),
            rounds=1,
            kfold=2,
            pstep=2,
            cleanup=True,
            verbose=False,
        )
        assert optimized["status"] == "completed"
    chosen = m.augustus_choose_parameter_set(str(GENOME), ["human", "toy"], verbose=False)
    assert chosen["species"] in {"human", "toy"}
    table = m.augustus_choose_parameter_set(
        str(GENOME),
        ["human", "toy"],
        mode="metagenomic",
        classification_path=str(tmp_path / "cls.tsv"),
        verbose=False,
    )
    assert table["sequences_classified"] >= 1
    predicted = m.augustus_predict_genes(
        str(GENOME), "toy", str(tmp_path / "toy.gtf"), verbose=False
    )
    assert predicted["status"] == "completed"
