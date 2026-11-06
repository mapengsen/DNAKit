"""Tests for the P19 Bio.motifs functions (the WebLogo web service is only used with DNAKIT_LIVE_TESTS=1)."""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio")

from dnakit.exceptions import ConfigurationError
from dnakit.regulation import biopython_p19_motif_functions as m

SITES = ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]


def _motif(pseudocounts=None):
    import warnings

    from Bio import motifs
    from Bio.Seq import Seq

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        motif = motifs.create([Seq(site) for site in SITES])
        if pseudocounts is not None:
            motif.pseudocounts = pseudocounts
    return motif


def test_motif_matrices_match_biopython() -> None:
    result = m.motif_matrices(SITES, pseudocounts=0.5, verbose=False)
    motif = _motif(0.5)
    for base in "ACGT":
        assert result["counts"][base] == list(motif.counts[base])
        assert result["pwm"][base] == pytest.approx(list(motif.pwm[base]))
        assert result["pssm"][base] == pytest.approx(list(motif.pssm[base]))
    assert result["consensus"] == "TACGC"
    assert result["degenerate_consensus"] == "WACVC"
    assert result["anticonsensus"] == "CCATG"
    assert result["instance_count"] == 7


def test_motif_matrices_infinity_and_counts_input() -> None:
    plain = m.motif_matrices(SITES, verbose=False)  # Biopython default: no pseudocounts
    assert plain["pssm"]["C"][0] is None  # -inf
    counts = {base: list(values) for base, values in plain["counts"].items()}
    again = m.motif_matrices(counts=counts, verbose=False)
    assert again["pwm"] == plain["pwm"]
    background = m.motif_matrices(
        SITES, pseudocounts=1, background={"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}, verbose=False
    )
    assert background["background"]["A"] == pytest.approx(0.3)
    with pytest.raises(ConfigurationError):
        m.motif_matrices(verbose=False)
    with pytest.raises(ConfigurationError):
        m.motif_matrices(SITES, counts=counts, verbose=False)
    with pytest.raises(ConfigurationError):
        m.motif_matrices(["ACGU", "ACGT"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.motif_matrices(["ACG", "ACGT"], verbose=False)


def test_anticonsensus() -> None:
    from_sites = m.motif_anticonsensus(SITES, verbose=False)
    counts = m.motif_matrices(SITES, verbose=False)["counts"]
    assert (
        m.motif_anticonsensus(counts=counts, verbose=False)["anticonsensus"]
        == from_sites["anticonsensus"]
    )
    assert str(_motif().anticonsensus) == from_sites["anticonsensus"]


def test_pssm_search_matches_biopython_both_strands() -> None:
    from Bio.Seq import Seq

    sequence = "TACACGTACGCAACGCATTTGTAAC"
    result = m.pssm_search(sequence, SITES, threshold=3.0, verbose=False)
    pssm = _motif(0.5).pssm
    expected = [(int(p), float(s)) for p, s in pssm.search(Seq(sequence), threshold=3.0)]
    assert [h["position"] for h in result["hits"]] == [position for position, _ in expected]
    assert [h["score"] for h in result["hits"]] == pytest.approx(
        [score for _, score in expected], rel=1e-6
    )
    assert any(h["strand"] == "-" for h in result["hits"])
    for hit in result["hits"]:
        window = sequence[hit["start"] : hit["end"]]
        assert len(window) == 5
        if hit["strand"] == "+":
            assert hit["start"] == hit["position"]
    only_forward = m.pssm_search(sequence, SITES, threshold=3.0, both_strands=False, verbose=False)
    assert all(h["strand"] == "+" for h in only_forward["hits"])
    with pytest.raises(ConfigurationError):
        m.pssm_search("ACG", SITES, verbose=False)
    with pytest.raises(ConfigurationError):
        m.pssm_search("ACGTXACGT", SITES, verbose=False)
    assert m.pssm_search("T" * 20, SITES, threshold=3.0, verbose=False)["hit_count"] == 0


def test_pssm_score_statistics_match_biopython() -> None:
    result = m.pssm_score_statistics(SITES, verbose=False)
    pssm = _motif(0.5).pssm
    assert result["mean"] == pytest.approx(pssm.mean())
    assert result["std"] == pytest.approx(pssm.std())
    assert result["max"] == pytest.approx(pssm.max)
    assert result["min"] == pytest.approx(pssm.min)
    skewed = m.pssm_score_statistics(
        SITES, background={"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4}, verbose=False
    )
    assert skewed["mean"] != pytest.approx(result["mean"])
    unseen = m.pssm_score_statistics(SITES, pseudocounts=None, verbose=False)
    assert unseen["min"] is None  # -inf


def test_distribution_and_thresholds_match_biopython() -> None:
    dist = _motif(0.5).pssm.distribution(precision=1000)
    result = m.pssm_score_distribution(SITES, precision=100, verbose=False)
    assert result["bin_count"] == 500
    assert sum(result["background_density"]) == pytest.approx(1.0)
    assert sum(result["motif_density"]) == pytest.approx(1.0)
    fpr = m.pssm_threshold(SITES, method="fpr", rate=0.01, verbose=False)
    assert fpr["threshold"] == pytest.approx(dist.threshold_fpr(0.01))
    fnr = m.pssm_threshold(SITES, method="fnr", rate=0.1, verbose=False)
    assert fnr["threshold"] == pytest.approx(dist.threshold_fnr(0.1))
    balanced = m.pssm_threshold(SITES, method="balanced", verbose=False)
    threshold, rate = dist.threshold_balanced(1.0, return_rate=True)
    assert (balanced["threshold"], balanced["rate"]) == pytest.approx((threshold, rate))
    patser = m.pssm_threshold(SITES, method="patser", verbose=False)
    assert patser["threshold"] == pytest.approx(dist.threshold_patser())
    assert fnr["threshold"] < patser["threshold"] < fpr["threshold"]
    with pytest.raises(ConfigurationError):
        m.pssm_threshold(SITES, method="fpr", verbose=False)  # rate missing
    with pytest.raises(ConfigurationError):
        m.pssm_threshold(SITES, method="other", verbose=False)
    with pytest.raises(ConfigurationError):
        m.pssm_score_distribution(SITES, precision=100_000, verbose=False)


def test_jaspar_pseudocounts() -> None:
    from Bio.motifs import jaspar

    result = m.jaspar_pseudocounts(SITES, verbose=False)
    expected = jaspar.calculate_pseudocounts(_motif())
    assert result["pseudocounts"] == pytest.approx(dict(expected))
    assert result["pseudocounts"]["A"] == pytest.approx(7**0.5 * 0.25)


def test_weblogo_validation(tmp_path) -> None:
    with pytest.raises(ConfigurationError):
        m.weblogo_service_logo(SITES, tmp_path / "x.png", image_format="gif", verbose=False)
    with pytest.raises(ConfigurationError):
        m.weblogo_service_logo(
            SITES, tmp_path / "x.png", options={"not_an_option": "1"}, verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.weblogo_service_logo(
            SITES, tmp_path / "x.png", options={"show_errorbars": "yes"}, verbose=False
        )


@pytest.mark.skipif(os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live WebLogo service test")
def test_weblogo_service_live(tmp_path) -> None:
    result = m.weblogo_service_logo(
        SITES, tmp_path / "logo.png", options={"logo_title": "test"}, verbose=False
    )
    assert result["bytes"] > 1000
    assert (tmp_path / "logo.png").read_bytes()[:4] == b"\x89PNG"
