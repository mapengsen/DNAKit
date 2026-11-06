"""Tests for dnakit.variants.allel_p32_functions (scikit-allel genotype, allele-count, position and state-run operations).

Expected values are written by hand from the definitions (a genotype call is homozygous when all alleles
agree, a dosage counts alleles, and so on) and, where scikit-allel has a method, also compared with the
direct scikit-allel result. Everything is offline.
"""

from __future__ import annotations

import inspect
import sys

import pytest

np = pytest.importorskip("numpy")
allel = pytest.importorskip("allel")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import allel_p32_functions as ap

# variants x samples x ploidy: hom-ref, het, hom-alt / het-12, missing, hom-2
G = [[[0, 0], [0, 1], [1, 1]], [[1, 2], [-1, -1], [2, 2]]]
H = [[0, 1, 1, -1], [2, 0, 0, 2]]
AC = [[2, 0, 0], [1, 1, 0], [0, 0, 0], [0, 3, 0], [1, 2, 3]]
CALL = ap.__dict__


def _run(name: str, *args, **kwargs):
    kwargs.setdefault("verbose", False)
    return getattr(ap, name)(*args, **kwargs)


def test_all_public_functions_have_keyword_verbose() -> None:
    assert len(ap.__all__) == 25
    for name in ap.__all__:
        parameter = inspect.signature(getattr(ap, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True


def test_missing_allel_raises_backend_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "allel", None)
    with pytest.raises(BackendUnavailableError, match="bio-scikit-allel"):
        _run("locate_genotype_calls", G)


@pytest.mark.parametrize(
    "kind, expected",
    [
        ("called", [[True, True, True], [True, False, True]]),
        ("missing", [[False, False, False], [False, True, False]]),
        ("het", [[False, True, False], [True, False, False]]),
        ("hom", [[True, False, True], [False, False, True]]),
        ("hom_ref", [[True, False, False], [False, False, False]]),
        ("hom_alt", [[False, False, True], [False, False, True]]),
    ],
)
def test_locate_genotype_calls_matches_definitions(kind: str, expected: list) -> None:
    result = _run("locate_genotype_calls", G, kind=kind)
    assert result["mask"] == expected
    assert result["count"] == sum(map(sum, expected))


def test_locate_genotype_calls_specific_allele_and_call() -> None:
    assert _run("locate_genotype_calls", G, kind="hom", allele=2)["mask"] == [
        [False] * 3,
        [False, False, True],
    ]
    assert _run("locate_genotype_calls", G, kind="call", call=[0, 1])["mask"] == [
        [False, True, False],
        [False] * 3,
    ]
    with pytest.raises(ConfigurationError):
        _run("locate_genotype_calls", G, kind="call")
    with pytest.raises(ConfigurationError):
        _run("locate_genotype_calls", G, kind="nonsense")


def test_locate_genotype_calls_rejects_bad_arrays() -> None:
    for bad in ([], [[0, 1]], [[[0.5, 1.0]]], [[[0, 1], [0]]]):
        with pytest.raises(ConfigurationError):
            _run("locate_genotype_calls", bad)


def test_locate_haplotype_calls() -> None:
    assert _run("locate_haplotype_calls", H, kind="ref")["mask"] == [
        [True, False, False, False],
        [False, True, True, False],
    ]
    assert _run("locate_haplotype_calls", H, kind="alt")["mask"] == [
        [False, True, True, False],
        [True, False, False, True],
    ]
    assert _run("locate_haplotype_calls", H, kind="missing")["count"] == 1
    assert _run("locate_haplotype_calls", H, kind="called")["count"] == 7
    assert _run("locate_haplotype_calls", H, kind="call", allele=2)["mask"] == [
        [False] * 4,
        [True, False, False, True],
    ]
    with pytest.raises(ConfigurationError):
        _run("locate_haplotype_calls", H, kind="call")


def test_count_genotype_calls_axes_match_allel() -> None:
    g = allel.GenotypeArray(G)
    assert _run("count_genotype_calls", G, kind="hom")["counts"] == int(g.count_hom())
    assert (
        _run("count_genotype_calls", G, kind="hom", axis=1)["counts"]
        == g.count_hom(axis=1).tolist()
        == [2, 1]
    )
    assert _run("count_genotype_calls", G, kind="hom", axis=0)["counts"] == [1, 0, 2]
    assert _run("count_genotype_calls", G, kind="het")["counts"] == 2
    assert _run("count_genotype_calls", G, kind="call", call=[2, 2])["counts"] == 1
    assert _run("count_genotype_calls", G, kind="hom_alt", axis=1)["counts"] == [1, 1]
    with pytest.raises(ConfigurationError):
        _run("count_genotype_calls", G, axis=2)


def test_count_haplotype_calls() -> None:
    assert _run("count_haplotype_calls", H, kind="ref")["counts"] == 3
    assert _run("count_haplotype_calls", H, kind="alt", axis=1)["counts"] == [2, 2]
    assert _run("count_haplotype_calls", H, kind="call", allele=2, axis=0)["counts"] == [1, 0, 0, 1]


def test_genotype_dosage_ref_and_alt() -> None:
    assert _run("genotype_dosage", G, allele="alt")["dosage"] == [[0, 1, 2], [2, 0, 2]]
    assert _run("genotype_dosage", G, allele="ref")["dosage"] == [[2, 1, 0], [0, 0, 0]]
    assert _run("genotype_dosage", G, allele="ref", fill=-1)["dosage"][1][1] == -1
    with pytest.raises(ConfigurationError):
        _run("genotype_dosage", G, allele="other")


def test_genotype_allele_counts_per_call_and_gt_strings() -> None:
    counts = _run("genotype_allele_counts_per_call", G)
    assert counts["max_allele"] == 2
    assert counts["counts"][0][1] == [1, 1, 0] and counts["counts"][1][1] == [0, 0, 0]
    assert _run("genotype_allele_counts_per_call", G, max_allele=3)["counts"][0][0] == [2, 0, 0, 0]
    assert _run("genotypes_to_gt_strings", G)["gt"] == [
        ["0/0", "0/1", "1/1"],
        ["1/2", "./.", "2/2"],
    ]


def test_remap_alleles_for_each_array_kind() -> None:
    swap = [[1, 0, 2], [0, 2, 1]]
    expected = allel.GenotypeArray(G).map_alleles(np.array(swap)).values.tolist()
    assert _run("remap_alleles", G, swap, kind="genotypes")["values"] == expected
    assert _run("remap_alleles", H, swap, kind="haplotypes")["values"][0] == [1, 0, 0, -1]
    counts = _run("remap_alleles", [[2, 1, 0], [0, 4, 5]], swap, kind="allele_counts")["values"]
    assert counts == [[1, 2, 0], [0, 5, 4]]
    with pytest.raises(ConfigurationError):
        _run("remap_alleles", G, [[0, 1, 2]], kind="genotypes")


def test_fill_masked_genotypes() -> None:
    result = _run(
        "fill_masked_genotypes", G, [[True, False, False], [False, False, True]], fill_value=-1
    )
    assert result["genotypes"][0][0] == [-1, -1] and result["genotypes"][1][2] == [-1, -1]
    assert result["genotypes"][0][1] == [0, 1]
    with pytest.raises(ConfigurationError):
        _run("fill_masked_genotypes", G, [[True, False]])


def test_haploidify_is_reproducible_and_picks_existing_alleles() -> None:
    first = _run("haploidify_genotypes", G, seed=7)["haplotypes"]
    assert first == _run("haploidify_genotypes", G, seed=7)["haplotypes"]
    for row, calls in zip(first, G, strict=False):
        for allele, call in zip(row, calls, strict=False):
            assert allele in call


def test_allele_count_site_masks_and_counts() -> None:
    expected = {
        "variant": [False, True, False, True, True],
        "non_variant": [True, False, True, False, False],
        "segregating": [False, True, False, False, True],
        "non_segregating": [True, False, True, True, False],
        "biallelic": [False, True, False, False, False],
        "biallelic_01": [False, True, False, False, False],
        "singleton": [False, True, False, False, True],
        "doubleton": [False, False, False, False, False],
    }
    # a site is a variant when a non-reference allele is observed; sites 0 and 2 (all reference / no calls) are not
    arr = allel.AlleleCountsArray(AC)
    for kind, _mask in expected.items():
        got = _run("locate_allele_count_sites", AC, kind=kind)
        method = getattr(arr, f"is_{kind}")
        direct = method(1) if kind in ("singleton", "doubleton") else method()
        assert got["mask"] == direct.tolist(), kind
        assert got["count"] == int(direct.sum())
    assert _run("locate_allele_count_sites", AC, kind="variant")["mask"] == expected["variant"]
    assert _run("locate_allele_count_sites", AC, kind="doubleton", allele=2)["mask"] == [
        False,
        False,
        False,
        False,
        False,
    ]
    assert _run("locate_allele_count_sites", AC, kind="singleton", allele=2)["mask"] == [
        False,
        False,
        False,
        False,
        False,
    ]
    for kind in (
        "variant",
        "non_variant",
        "segregating",
        "non_segregating",
        "singleton",
        "doubleton",
    ):
        direct = (
            getattr(arr, f"count_{kind}")(1)
            if kind in ("singleton", "doubleton")
            else getattr(arr, f"count_{kind}")()
        )
        assert _run("count_allele_count_sites", AC, kind=kind)["count"] == int(direct)
    with pytest.raises(ConfigurationError):
        _run("count_allele_count_sites", AC, kind="biallelic")


def test_summarize_alleles_per_variant() -> None:
    assert _run("summarize_alleles_per_variant", AC, kind="allelism")["values"] == [1, 2, 0, 1, 3]
    assert _run("summarize_alleles_per_variant", AC, kind="max_allele")["values"] == [
        0,
        1,
        -1,
        1,
        2,
    ]


def test_feature_coverage_mask() -> None:
    result = _run("feature_coverage_mask", [[1, 5], [8, 9]], 10)
    assert result["mask"] == [True] * 5 + [False, False, True, True, False]
    assert result["covered"] == 7
    with pytest.raises(ConfigurationError):
        _run("feature_coverage_mask", [[5, 1]], 10)


def test_locate_sorted_positions_single_and_many() -> None:
    assert _run("locate_sorted_positions", [10, 20, 30, 40], 30)["located"] == 2
    many = _run("locate_sorted_positions", [10, 20, 30, 40], [20, 25, 40], strict=False)
    assert many["mask"] == [False, True, False, True] and many["count"] == 2
    with pytest.raises(ConfigurationError):
        _run("locate_sorted_positions", [10, 20, 30], 25)
    with pytest.raises(ConfigurationError):
        _run("locate_sorted_positions", [10, 20, 30], [20, 25])
    with pytest.raises(ConfigurationError):
        _run("locate_sorted_positions", [30, 20], 20)
    with pytest.raises(ConfigurationError):
        _run("locate_sorted_positions", [10, 10, 20], 10)


def test_intersect_sorted_positions() -> None:
    result = _run("intersect_sorted_positions", [1, 3, 5, 7], [3, 7, 9])
    assert result["positions"] == [3, 7]
    assert result["mask_positions"] == [False, True, False, True]
    assert result["mask_other"] == [True, True, False]
    assert _run("intersect_sorted_positions", [1, 2], [5, 6])["positions"] == []


def test_select_positions_and_overlap_with_ranges() -> None:
    selected = _run("select_positions_in_ranges", [1, 3, 5, 7, 9], [[2, 5], [8, 20]])
    assert selected["positions"] == [3, 5, 9] and selected["mask"] == [
        False,
        True,
        True,
        False,
        True,
    ]
    single = _run("select_positions_in_ranges", [1, 3, 5, 7, 9], [[3, 3]])
    assert single["positions"] == [3]
    with pytest.raises(ConfigurationError):
        _run("select_positions_in_ranges", [1, 3], [[10, 20]], strict=True)
    overlap = _run("overlap_positions_with_ranges", [1, 3, 5, 7, 9], [[2, 5], [100, 200]])
    assert overlap["position_mask"] == [False, True, True, False, False]
    assert overlap["range_mask"] == [True, False]
    with pytest.raises(ConfigurationError):
        _run("overlap_positions_with_ranges", [1, 3], [[5, 2]])


def test_create_allele_index_mapping_matches_allel() -> None:
    ref, alt, alleles = ["A", "C"], [["G", "T"], ["A", ""]], [["A", "G"], ["A", "C"]]
    got = _run("create_allele_index_mapping", ref, alt, alleles)["mapping"]
    assert got == [[0, 1, -1], [1, 0, -1]]
    assert (
        got == allel.create_allele_mapping(np.array(ref), np.array(alt), np.array(alleles)).tolist()
    )
    with pytest.raises(ConfigurationError):
        _run("create_allele_index_mapping", ref, alt, [["A", "G"]])


def test_principal_coordinates_analysis_matches_allel_and_checks_input() -> None:
    d = np.array([[0, 1, 2, 3], [1, 0, 1.5, 2], [2, 1.5, 0, 1], [3, 2, 1, 0]], dtype=float)
    result = _run("principal_coordinates_analysis", d.tolist())
    coords, explained = allel.pcoa(d)
    assert np.allclose(result["coordinates"], coords) and np.allclose(
        result["explained_ratio"], explained
    )
    assert abs(sum(result["explained_ratio"]) - 1.0) < 1e-9
    with pytest.raises(ConfigurationError):
        _run("principal_coordinates_analysis", [[0, 1], [1, 0]])
    with pytest.raises(ConfigurationError):
        _run("principal_coordinates_analysis", [[0, 1, 2], [3, 0, 1], [2, 1, 0]])


def test_tabulate_state_transitions_and_blocks() -> None:
    series, states, pos = [0, 0, 1, 1, 2], [0, 1, 2], [1, 3, 5, 7, 9]
    transitions = _run("tabulate_state_transitions", series, states, pos=pos)["transitions"]
    assert [(r["lstate"], r["rstate"], r["lpos"], r["rpos"]) for r in transitions] == [
        (-1, 0, -1, 1),
        (0, 1, 3, 5),
        (1, 2, 7, 9),
        (2, -1, 9, -1),
    ]
    blocks = _run("tabulate_state_blocks", series, states, pos=pos)["blocks"]
    assert [(b["state"], b["support"]) for b in blocks] == [(0, 2), (1, 2), (2, 1)]
    no_pos = _run("tabulate_state_transitions", series, states)["transitions"]
    assert "lpos" not in no_pos[0]
    with pytest.raises(ConfigurationError):
        _run("tabulate_state_transitions", series, states, pos=[1, 2])


def test_scale_genotype_dosages_center_and_patterson() -> None:
    dosage = [[0, 1, 2, 1], [2, 2, 1, 1]]
    centred = _run("scale_genotype_dosages", dosage, method="center")["scaled"]
    assert centred[0] == [-1.0, 0.0, 1.0, 0.0]
    patterson = np.array(_run("scale_genotype_dosages", dosage, method="patterson")["scaled"])
    mean = np.mean(dosage, axis=1, keepdims=True)
    p = mean / 2
    assert np.allclose(patterson, (np.array(dosage) - mean) / np.sqrt(p * (1 - p)))
    with pytest.raises(ConfigurationError):
        _run("scale_genotype_dosages", dosage, method="zscore")


def test_fold_site_frequency_spectrum() -> None:
    assert _run("fold_site_frequency_spectrum", [0, 5, 3, 2, 1, 0], 5)["folded"] == [0, 6, 5]
    assert _run("fold_site_frequency_spectrum", [0, 4, 3, 2, 1, 5, 0, 0], 7)["folded"] == [
        0,
        4,
        8,
        3,
    ]
    assert _run("fold_site_frequency_spectrum", [0, 5, 3], 5)["folded"] == [
        0,
        5,
        3,
    ]  # zero-padded to n + 1
    with pytest.raises(ConfigurationError):
        _run("fold_site_frequency_spectrum", [0, 5, 3, 1], 2)


def test_write_sites_vcf_is_readable_and_refuses_overwrite(tmp_path) -> None:
    target = tmp_path / "sites.vcf"
    records = [
        {"chrom": "1", "pos": 5, "ref": "A", "alt": ["G", "T"], "qual": 30},
        {"chrom": "1", "pos": 9, "ref": "C", "alt": "T", "id": "rs1", "filters": ["LowQual"]},
    ]
    result = _run("write_sites_vcf", str(target), records)
    assert result["variant_count"] == 2
    rows = [
        line.split("\t") for line in target.read_text().splitlines() if not line.startswith("#")
    ]
    assert rows[0][:8] == ["1", "5", ".", "A", "G,T", "30", "PASS", "."]
    assert rows[1][:8] == ["1", "9", "rs1", "C", "T", ".", "LowQual", "."]
    parsed = allel.read_vcf(
        str(target), fields=["variants/POS", "variants/ALT", "variants/FILTER_LowQual"]
    )
    assert parsed["variants/POS"].tolist() == [5, 9]
    assert (
        parsed["variants/ALT"][0].tolist() == ["G", "T", ""]
        and parsed["variants/ALT"][1].tolist()[0] == "T"
    )
    assert parsed["variants/FILTER_LowQual"].tolist() == [False, True]
    with pytest.raises(ConfigurationError):
        _run("write_sites_vcf", str(target), records)
    with pytest.raises(ConfigurationError):
        _run(
            "write_sites_vcf",
            str(tmp_path / "x.vcf"),
            [{"chrom": "1", "pos": 0, "ref": "A", "alt": "G"}],
        )


def test_large_genomic_coordinates_and_counts_are_accepted() -> None:
    positions = [1_000_000, 2_500_000, 99_000_000]
    assert _run("locate_sorted_positions", positions, 2_500_000)["located"] == 1
    assert _run("select_positions_in_ranges", positions, [[1_000_001, 100_000_000]])[
        "positions"
    ] == [2_500_000, 99_000_000]
    assert (
        _run("tabulate_state_blocks", [0, 0, 1], [0, 1], pos=positions)["blocks"][0]["support"] == 2
    )
    assert _run("feature_coverage_mask", [[5_000, 9_000]], 12_000)["covered"] == 4_001
    assert _run("summarize_alleles_per_variant", [[250_000, 3]], kind="allelism")["values"] == [2]
    assert _run("fold_site_frequency_spectrum", [0, 5_000, 3_000], 2)["folded"] == [5_000]


def test_undefined_scaled_values_are_null() -> None:
    scaled = _run("scale_genotype_dosages", [[0, 0, 0], [0, 1, 2]], method="patterson")["scaled"]
    assert scaled[0] == [None, None, None] and scaled[1] == [-2.0, 0.0, 2.0]
