"""Tests for dnakit.population.dadi_p23_phi_functions (skipped when dadi is missing)."""

from __future__ import annotations

import pytest

dadi = pytest.importorskip("dadi")
np = pytest.importorskip("numpy")

from dnakit.exceptions import ConfigurationError
from dnakit.population import dadi_p23_phi_functions as m

PTS = 16
XX = dadi.Numerics.default_grid(PTS)

SNP_TEXT = """Context Out Allele1 A B Allele2 A B ID
ACA ACA C 4 3 T 2 3 s1
CCG CCG C 5 2 T 1 4 s2
GCT GCT G 3 3 A 3 3 s3
TCA TCA C 6 1 T 0 5 s4
ACG ACG C 2 2 T 4 4 s5
"""


def test_one_population_models_match_dadi():
    for model, ref in (
        ("dominance", dadi.PhiManip.phi_1D(XX, gamma=1.0, h=0.3)),
        ("genic", dadi.PhiManip.phi_1D_genic(XX, gamma=1.0)),
        ("neutral", dadi.PhiManip.phi_1D_snm(XX)),
        ("x_linked", dadi.PhiManip.phi_1D_X(XX, gamma=1.0, alpha=2.0)),
    ):
        kwargs = {"gamma": 1.0}
        if model == "dominance":
            kwargs["h"] = 0.3
        if model == "x_linked":
            kwargs["alpha"] = 2.0
        if model == "neutral":
            kwargs = {}
        got = m.phi_one_population(grid_points=PTS, model=model, **kwargs)
        assert np.allclose(got["phi"], ref)
        assert got["status"] == "completed"
    with pytest.raises(ConfigurationError):
        m.phi_one_population(grid_points=PTS, model="bad")


def test_events_match_dadi_and_shapes():
    phi1 = dadi.PhiManip.phi_1D(XX)
    two = m.phi_split_population(phi1.tolist(), grid_points=PTS, event="1_to_2")
    assert np.allclose(two["phi"], dadi.PhiManip.phi_1D_to_2D(XX, phi1))
    three = m.phi_split_population(two["phi"], grid_points=PTS, event="2_to_3_split_2")
    assert three["shape"] == [PTS] * 3
    ref = dadi.PhiManip.phi_2D_to_3D_split_1(XX, np.array(two["phi"]))
    assert np.allclose(
        m.phi_split_population(two["phi"], grid_points=PTS, event="2_to_3_split_1")["phi"], ref
    )
    adm = m.phi_admix_into_population(
        three["phi"], grid_points=PTS, target_population=2, proportions=[0.3, 0.2]
    )
    ref = dadi.PhiManip.phi_3D_admix_1_and_3_into_2(np.array(three["phi"]), 0.3, 0.2, XX, XX, XX)
    assert np.allclose(adm["phi"], ref)
    two_adm = m.phi_admix_into_population(
        two["phi"], grid_points=PTS, target_population=1, proportions=[0.4]
    )
    assert np.allclose(
        two_adm["phi"], dadi.PhiManip.phi_2D_admix_2_into_1(np.array(two["phi"]), 0.4, XX, XX)
    )
    new = m.phi_create_admixed_population(two["phi"], grid_points=PTS, proportions=[0.4])
    assert np.allclose(
        new["phi"], dadi.PhiManip.phi_2D_to_3D_admix(np.array(two["phi"]), 0.4, XX, XX, XX)
    )
    four = m.phi_create_admixed_population(three["phi"], grid_points=PTS, proportions=[0.2, 0.3])
    assert four["num_populations"] == 4
    into4 = m.phi_admix_into_population(
        four["phi"], grid_points=PTS, target_population=1, proportions=[0.1, 0.2, 0.3]
    )
    assert into4["shape"] == [PTS] * 4


def test_select_populations(tmp_path):
    phi1 = dadi.PhiManip.phi_1D(XX)
    two = np.array(m.phi_split_population(phi1.tolist(), grid_points=PTS, event="1_to_2")["phi"])
    three = np.array(
        m.phi_split_population(two.tolist(), grid_points=PTS, event="2_to_3_split_2")["phi"]
    )
    out = tmp_path / "phi.npy"
    removed = m.phi_select_populations(
        three.tolist(), grid_points=PTS, operation="remove", populations=[2], output_path=out
    )
    assert np.allclose(np.load(out), dadi.PhiManip.remove_pop(three, XX, 2))
    assert removed["phi_path"] == str(out)
    kept = m.phi_select_populations(str(out), grid_points=PTS, operation="keep", populations=[1])
    assert kept["num_populations"] == 1
    re = m.phi_select_populations(
        three.tolist(), grid_points=PTS, operation="reorder", populations=[3, 1, 2]
    )
    assert np.allclose(re["phi"], dadi.PhiManip.reorder_pops(three, [3, 1, 2]))
    with pytest.raises(ConfigurationError):
        m.phi_select_populations(
            three.tolist(), grid_points=PTS, operation="reorder", populations=[1, 2]
        )
    with pytest.raises(ConfigurationError):  # existing output is never overwritten
        m.phi_select_populations(
            three.tolist(), grid_points=PTS, operation="remove", populations=[1], output_path=out
        )


def test_spectrum_from_phi_and_inbreeding():
    phi = dadi.PhiManip.phi_1D(XX)
    got = m.spectrum_from_phi(phi.tolist(), grid_points=PTS, sample_sizes=[8])
    ref = dadi.Spectrum.from_phi(phi, (8,), (XX,))
    assert np.allclose(np.array(got["spectrum"][1:-1]), ref[1:-1])
    inb = m.spectrum_from_phi(
        phi.tolist(), grid_points=PTS, sample_sizes=[8], inbreeding_coefficients=[0.2], ploidies=[2]
    )
    ref = dadi.Spectrum.from_phi_inbreeding(phi, (8,), (XX,), [0.2], [2])
    assert inb["mode"] == "from_phi_inbreeding"
    assert np.allclose(np.array(inb["spectrum"][1:-1]), ref[1:-1])
    with pytest.raises(ConfigurationError):
        m.spectrum_from_phi(phi.tolist(), grid_points=PTS, sample_sizes=[8, 8])


def test_scramble_is_seeded():
    fs = dadi.Spectrum(np.arange(1.0, 26.0).reshape(5, 5))
    one = m.scramble_spectrum_pop_ids(fs.tolist(), seed=3)
    two = m.scramble_spectrum_pop_ids(fs.tolist(), seed=3)
    assert one["spectrum"] == two["spectrum"]
    assert one["shape"] == [5, 5]
    with pytest.raises(ConfigurationError):
        m.scramble_spectrum_pop_ids([1.0, 2.0, 3.0, 4.0])


def test_spectrum_from_snp_data(tmp_path):
    snp = tmp_path / "toy.txt"
    snp.write_text(SNP_TEXT)
    data = dadi.Misc.make_data_dict(str(snp))
    for polarized in (True, False):
        got = m.spectrum_from_snp_data(
            snp, pop_ids=["A", "B"], projections=[4, 4], polarized=polarized
        )
        ref = dadi.Spectrum.from_data_dict(data, ["A", "B"], [4, 4], polarized=polarized)
        assert got["folded"] == (not polarized)
        assert got["segregating_sites"] == pytest.approx(float(ref.S()))
    fux = tmp_path / "fux.txt"
    lines = [
        f"{a}{b}{c} {u} 0.01"
        for a in "ACGT"
        for b in "ACGT"
        for c in "ACGT"
        for u in "ACGT"
        if u != b
    ]
    fux.write_text("\n".join(lines) + "\n")
    cor = m.spectrum_from_snp_data(snp, pop_ids=["A"], projections=[6], misidentification_path=fux)
    ref = dadi.Spectrum.from_data_dict_corrected(data, ["A"], [6], str(fux))
    assert cor["mode"] == "from_data_dict_corrected"
    assert cor["segregating_sites"] == pytest.approx(float(ref.S()))
    with pytest.raises(ConfigurationError):
        m.spectrum_from_snp_data(snp, pop_ids=["Z"], projections=[4])
