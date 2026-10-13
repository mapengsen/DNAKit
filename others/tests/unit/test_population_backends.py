"""Contract tests for optional tskit and dadi population-genetics adapters."""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population.backends import (
    build_tree_sequence,
    calculate_diversity,
    calculate_fst,
    calculate_sfs_log_likelihood,
    calculate_site_frequency_spectrum,
    calculate_tajimas_d,
    estimate_sfs_theta,
    expected_neutral_sfs,
    fold_sfs,
    read_tree_sequence,
    write_tree_sequence,
)


def _nodes() -> list[dict[str, Any]]:
    """Four sample nodes (0-3), coalescent nodes 4 and 5 at times 1 and 2, root 6 at time 3."""

    return [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}]


def _edges(pairs: list[tuple[int, int]]) -> list[dict[str, Any]]:
    return [{"left": 0, "right": 100, "parent": p, "child": c} for c, p in pairs]


def _tables() -> dict[str, Any]:
    """One tree ((0,1),(2,3)) on 0..100; singleton 0 at 10, clade {0,1} at 20, clade {2,3} at 30."""

    return {
        "sequence_length": 100,
        "nodes": _nodes(),
        "edges": _edges([(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]),
        "sites": [
            {"position": 10, "ancestral_state": "A"},
            {"position": 20, "ancestral_state": "A"},
            {"position": 30, "ancestral_state": "A"},
        ],
        "mutations": [
            {"site": 0, "node": 0, "derived_state": "G"},
            {"site": 1, "node": 4, "derived_state": "G"},
            {"site": 2, "node": 5, "derived_state": "G"},
        ],
    }


def _caterpillar_tables() -> dict[str, Any]:
    """Tree ((0,1,2),3) with one site carried by samples 0, 1 and 2 (derived count 3 of 4)."""

    return {
        "sequence_length": 100,
        "nodes": [*_nodes()[:4], {"time": 1}, {"time": 2}],
        "edges": _edges([(0, 4), (1, 4), (2, 4), (3, 5), (4, 5)]),
        "sites": [{"position": 10, "ancestral_state": "A"}],
        "mutations": [{"site": 0, "node": 4, "derived_state": "G"}],
    }


def _direct_tree_sequence(tskit: Any) -> Any:
    """The same fixture as _tables(), built with tskit calls only, for direct comparisons."""

    tables = tskit.TableCollection(sequence_length=100)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    for time in (1, 2, 3):
        tables.nodes.add_row(flags=0, time=time)
    for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]:
        tables.edges.add_row(left=0, right=100, parent=parent, child=child)
    for position in (10, 20, 30):
        tables.sites.add_row(position=position, ancestral_state="A")
    for site, node in [(0, 0), (1, 4), (2, 5)]:
        tables.mutations.add_row(site=site, node=node, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


def _tajima_d_closed_form(n: int, segregating: int, pi: Fraction) -> float:
    """Tajima (1989) with exact fractions: (pi - S/a1) / sqrt(e1 S + e2 S(S-1))."""

    a1 = sum(Fraction(1, i) for i in range(1, n))
    a2 = sum(Fraction(1, i * i) for i in range(1, n))
    b1 = Fraction(n + 1, 3 * (n - 1))
    b2 = Fraction(2 * (n * n + n + 3), 9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - Fraction(n + 2, a1 * n) + a2 / a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    variance = e1 * segregating + e2 * segregating * (segregating - 1)
    return float(pi - Fraction(segregating) / a1) / math.sqrt(float(variance))


def test_build_tree_sequence_summary_matches_direct_tskit() -> None:
    tskit = pytest.importorskip("tskit")
    result = build_tree_sequence(_tables(), verbose=False)
    direct = _direct_tree_sequence(tskit)
    assert result["sequence_length"] == 100.0
    assert result["num_trees"] == direct.num_trees == 1
    assert result["num_samples"] == direct.num_samples == 4
    assert result["num_nodes"] == direct.num_nodes == 7
    assert result["num_edges"] == direct.num_edges == 6
    assert result["num_sites"] == direct.num_sites == 3
    assert result["num_mutations"] == direct.num_mutations == 3
    assert result["provenance"]["package"].startswith("tskit ")


def test_build_tree_sequence_sorts_rows_before_building() -> None:
    pytest.importorskip("tskit")
    tables = _tables()
    # Reverse the site rows and renumber the mutations so they still point at the same positions.
    tables["sites"] = list(reversed(tables["sites"]))
    tables["mutations"] = [dict(row, site=2 - row["site"]) for row in tables["mutations"]]
    result = calculate_site_frequency_spectrum(tables, verbose=False)
    assert result["spectra"]["all"] == [0, 1, 2, 0, 0]


@pytest.mark.parametrize(
    ("change", "pattern"),
    [
        (lambda t: {**t, "extra": 1}, "unknown keys"),
        (lambda t: {k: v for k, v in t.items() if k != "nodes"}, "must include sequence_length"),
        (lambda t: {**t, "sequence_length": True}, "finite number"),
        (lambda t: {**t, "sequence_length": 0}, "greater than zero"),
        (lambda t: {**t, "nodes": []}, "max_rows|at least one row"),
        (
            lambda t: {**t, "nodes": [{"time": 0, "sample": "yes"}, *t["nodes"][1:]]},
            "true or false",
        ),
        (lambda t: {**t, "nodes": [{"sample": True}, *t["nodes"][1:]]}, "missing keys"),
        (
            lambda t: {**t, "edges": [{**t["edges"][0], "weight": 1}, *t["edges"][1:]]},
            "unknown keys",
        ),
        (
            lambda t: {**t, "mutations": [{**t["mutations"][0], "node": 99}, *t["mutations"][1:]]},
            "tskit",
        ),
        (
            lambda t: {**t, "edges": [{**t["edges"][0], "right": 200}, *t["edges"][1:]]},
            "tskit",
        ),
    ],
)
def test_build_tree_sequence_rejects_invalid_tables(change: Any, pattern: str) -> None:
    pytest.importorskip("tskit")
    with pytest.raises(ConfigurationError, match=pattern):
        build_tree_sequence(change(_tables()), verbose=False)


def test_build_tree_sequence_enforces_max_rows() -> None:
    pytest.importorskip("tskit")
    with pytest.raises(ConfigurationError, match="max_rows=2"):
        build_tree_sequence(_tables(), max_rows=2, verbose=False)


def test_write_then_read_tree_sequence_round_trip(tmp_path: Path) -> None:
    pytest.importorskip("tskit")
    path = tmp_path / "fixture.trees"
    written = write_tree_sequence(_tables(), path, verbose=False)
    read = read_tree_sequence(path, verbose=False)
    for key in ("num_trees", "num_samples", "num_nodes", "num_edges", "num_sites", "num_mutations"):
        assert written[key] == read[key]
    assert read["bytes"] == written["bytes"] == path.stat().st_size
    with pytest.raises(ConfigurationError, match="already exists"):
        write_tree_sequence(_tables(), path, verbose=False)


def test_read_and_write_reject_bad_files(tmp_path: Path) -> None:
    pytest.importorskip("tskit")
    path = tmp_path / "fixture.trees"
    write_tree_sequence(_tables(), path, verbose=False)
    with pytest.raises(ConfigurationError, match="max_bytes"):
        read_tree_sequence(path, max_bytes=10, verbose=False)
    plain = tmp_path / "plain.trees"
    plain.write_text("not a tree sequence", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="could not read"):
        read_tree_sequence(plain, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        write_tree_sequence(_tables(), tmp_path / "missing" / "out.trees", verbose=False)


def test_diversity_matches_hand_calculation_and_direct_tskit() -> None:
    tskit = pytest.importorskip("tskit")
    # Pairwise differences per site: 3 (singleton), 4 and 4, over C(4, 2) = 6 pairs and 100 bases.
    result = calculate_diversity(_tables(), verbose=False)
    direct = _direct_tree_sequence(tskit).diversity(sample_sets=[[0, 1, 2, 3]])[0]
    assert result["diversity_per_base"]["all"] == pytest.approx((3 + 4 + 4) / 6 / 100)
    assert result["diversity_per_base"]["all"] == pytest.approx(float(direct))  # per-base default
    # Within {0, 1} one pair differs at site 10; within {2, 3} no pair differs.
    within = calculate_diversity(_tables(), sample_sets={"A": [0, 1], "B": [2, 3]}, verbose=False)
    assert within["diversity_per_base"]["A"] == pytest.approx(1 / 100)
    assert within["diversity_per_base"]["B"] == pytest.approx(0.0)


def test_tajimas_d_matches_closed_form_and_direct_tskit() -> None:
    tskit = pytest.importorskip("tskit")
    # n = 4, S = 3 and pi = 11/6 (sum of pairwise differences over the three sites).
    result = calculate_tajimas_d(_tables(), verbose=False)
    direct = _direct_tree_sequence(tskit).Tajimas_D(sample_sets=[[0, 1, 2, 3]])[0]
    assert result["segregating_sites"]["all"] == 3
    assert result["tajimas_d"]["all"] == pytest.approx(_tajima_d_closed_form(4, 3, Fraction(11, 6)))
    assert result["tajimas_d"]["all"] == pytest.approx(float(direct))


def test_tajimas_d_is_none_without_segregating_sites() -> None:
    pytest.importorskip("tskit")
    tables = _tables()
    tables["sites"] = []
    tables["mutations"] = []
    result = calculate_tajimas_d(tables, verbose=False)
    assert result["segregating_sites"]["all"] == 0
    assert result["tajimas_d"]["all"] is None


def test_fst_matches_hand_calculation_and_direct_tskit() -> None:
    tskit = pytest.importorskip("tskit")
    # d_A = 1 (one pair differs), d_B = 0 and d_AB = 10 differences over 4 pairs = 2.5.
    # Fst = 1 - 2(1 + 0) / (1 + 2 * 2.5 + 0) = 2/3.
    sets = {"A": [0, 1], "B": [2, 3]}
    result = calculate_fst(_tables(), sample_sets=sets, verbose=False)
    direct = _direct_tree_sequence(tskit).Fst(sample_sets=[[0, 1], [2, 3]])
    assert result["fst"] == pytest.approx(2 / 3)
    assert result["fst"] == pytest.approx(float(direct))
    assert result["population_a"] == {"name": "A", "samples": 2}


@pytest.mark.parametrize(
    ("sample_sets", "pattern"),
    [
        ({"A": [0, 1, 2, 3]}, "exactly two sample sets"),
        ({"A": [0, 1], "B": [2, 4]}, "non-sample nodes"),
        ({"A": [0, 1], "B": [2, 2]}, "more than once"),
        ({"A": [0], "B": [2, 3]}, "at least two sample nodes"),
    ],
)
def test_fst_rejects_invalid_sample_sets(sample_sets: dict[str, list[int]], pattern: str) -> None:
    pytest.importorskip("tskit")
    with pytest.raises(ConfigurationError, match=pattern):
        calculate_fst(_tables(), sample_sets=sample_sets, verbose=False)


def test_site_frequency_spectrum_counts_and_folds_derived_alleles() -> None:
    tskit = pytest.importorskip("tskit")
    unfolded = calculate_site_frequency_spectrum(_tables(), verbose=False)
    # Singleton at index 1, two doubletons at index 2.
    assert unfolded["spectra"]["all"] == [0, 1, 2, 0, 0]
    # A site carried by three of four samples folds onto index 1.
    caterpillar = _caterpillar_tables()
    assert calculate_site_frequency_spectrum(caterpillar, verbose=False)["spectra"]["all"] == [
        0,
        0,
        0,
        1,
        0,
    ]
    folded = calculate_site_frequency_spectrum(caterpillar, folded=True, verbose=False)
    assert folded["spectra"]["all"] == [0, 1, 0, 0, 0]
    direct = _direct_tree_sequence(tskit).allele_frequency_spectrum(
        sample_sets=[[0, 1, 2, 3]], polarised=True, span_normalise=False
    )
    assert unfolded["spectra"]["all"] == [int(value) for value in direct]


def test_missing_tskit_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "tskit", None)
    with pytest.raises(BackendUnavailableError, match="bio-tskit"):
        build_tree_sequence(_tables(), verbose=False)


def test_expected_neutral_sfs_matches_theta_over_i_and_direct_dadi() -> None:
    dadi = pytest.importorskip("dadi")
    spectrum = expected_neutral_sfs(4, verbose=False)["spectrum"]
    assert spectrum[0] is None and spectrum[4] is None  # monomorphic classes are masked
    # Standard neutral model: E[xi_i] = theta / i; grid error is about 2e-5 at 1000 points.
    for index in (1, 2, 3):
        assert spectrum[index] == pytest.approx(1 / index, rel=1e-4)
    direct = dadi.Demographics1D.snm(None, (4,), 1000)
    assert spectrum[1] == pytest.approx(float(direct[1]), rel=1e-12)
    scaled = expected_neutral_sfs(4, theta=2.5, verbose=False)["spectrum"]
    assert scaled[1] == pytest.approx(2.5 * spectrum[1], rel=1e-12)


def test_folded_expected_sfs_matches_hand_fold() -> None:
    pytest.importorskip("dadi")
    # theta = 2, n = 4: minor classes are 2(1 + 1/3) = 8/3 and 2 * (1/2) = 1.
    spectrum = expected_neutral_sfs(4, theta=2.0, folded=True, verbose=False)["spectrum"]
    assert spectrum[1] == pytest.approx(8 / 3, rel=1e-4)
    assert spectrum[2] == pytest.approx(1.0, rel=1e-4)
    assert spectrum[3] is None and spectrum[4] is None


def test_fold_sfs_merges_classes_for_even_and_odd_sample_sizes() -> None:
    pytest.importorskip("dadi")
    # n = 4: classes 1 and 3 merge (2 + 1 = 3); class 2 stays (1).
    assert fold_sfs([0, 2, 1, 1, 0], verbose=False)["spectrum"] == [None, 3.0, 1.0, None, None]
    # n = 5: 1 with 4 gives 5, and 2 with 3 gives 5.
    assert fold_sfs([0, 1, 2, 3, 4, 0], verbose=False)["spectrum"] == [
        None,
        5.0,
        5.0,
        None,
        None,
        None,
    ]


def test_estimate_sfs_theta_matches_harmonic_ratio() -> None:
    pytest.importorskip("dadi")
    # theta_MLE = sum(observed classes 1..n-1) / sum(1/i) = 4 / (1 + 1/2 + 1/3) = 24/11.
    result = estimate_sfs_theta([0, 2, 1, 1, 0], verbose=False)
    assert result["theta"] == pytest.approx(24 / 11, rel=1e-4)


def test_log_likelihood_matches_poisson_formula_and_dadi_multinomial() -> None:
    dadi = pytest.importorskip("dadi")
    observed = [0, 2, 1, 1, 0]
    # Poisson terms x log m - m - log(x!) for classes 1..3 with m = 1 / i (theta = 1).
    expected = sum(
        count * math.log(1 / index) - 1 / index - math.lgamma(count + 1)
        for index, count in ((1, 2), (2, 1), (3, 1))
    )
    given = calculate_sfs_log_likelihood(observed, theta=1.0, verbose=False)
    assert given["theta_source"] == "given"
    assert given["log_likelihood"] == pytest.approx(expected, abs=1e-4)
    fitted = calculate_sfs_log_likelihood(observed, verbose=False)
    assert fitted["theta_source"] == "maximum likelihood"
    assert fitted["theta"] == pytest.approx(24 / 11, rel=1e-4)
    model = dadi.Demographics1D.snm(None, (4,), 1000)
    direct = float(dadi.Inference.ll_multinom(model, dadi.Spectrum(observed)))
    assert fitted["log_likelihood"] == pytest.approx(direct, rel=1e-9)


@pytest.mark.parametrize(
    ("call", "pattern"),
    [
        (lambda: expected_neutral_sfs(True, verbose=False), "sample_size must be an integer"),
        (lambda: expected_neutral_sfs(1, verbose=False), "sample_size must be an integer"),
        (
            lambda: expected_neutral_sfs(4, theta=0, verbose=False),
            "theta must be greater than zero",
        ),
        (lambda: expected_neutral_sfs(4, grid_points=10, verbose=False), "grid_points"),
        (lambda: expected_neutral_sfs(4, folded="yes", verbose=False), "folded must be"),
        (lambda: fold_sfs([1, 2], verbose=False), "observed must be a list"),
        (lambda: fold_sfs([0, -1, 1], verbose=False), "non-negative"),
        (lambda: fold_sfs([0, 1, 2, 1, 0], max_sample_size=3, verbose=False), "max_sample_size=3"),
        (lambda: estimate_sfs_theta([0, 0, 0, 0], verbose=False), "at least one count"),
        (lambda: calculate_sfs_log_likelihood([0, 0, 0, 0], verbose=False), "at least one count"),
        (lambda: calculate_sfs_log_likelihood([0, 1, 1, 0], theta=True, verbose=False), "finite"),
    ],
)
def test_dadi_adapters_reject_invalid_input(call: Any, pattern: str) -> None:
    pytest.importorskip("dadi")
    with pytest.raises(ConfigurationError, match=pattern):
        call()


def test_missing_dadi_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "dadi", None)
    with pytest.raises(BackendUnavailableError, match="bio-dadi"):
        expected_neutral_sfs(4, verbose=False)
