"""Contract tests for DnaChisel, DnaCauldron and pydna functions in dnakit.molecular_lab (G16a)."""

from __future__ import annotations

import inspect
import os
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import dnacauldron_functions as dnacauldron_api
from dnakit.molecular_lab import dnachisel_functions as dnachisel_api
from dnakit.molecular_lab import pydna_functions as pydna_api
from dnakit.molecular_lab.dnacauldron_functions import (
    autoselect_golden_gate_enzyme,
    crop_dna_record_with_features,
    list_type2s_digestion_overhangs,
)
from dnakit.molecular_lab.dnachisel_functions import (
    find_sequence_difference_segments,
    generate_constrained_dna_sequence,
    group_nearby_positions,
    list_common_restriction_enzymes,
)
from dnakit.molecular_lab.pydna_functions import (
    assemble_pydna_fragments,
    fetch_genbank_record,
    find_common_substrings,
    find_terminal_overlaps,
    parse_primer_sequences,
    primer_melting_temperature,
    primer_pair_hybridization,
    simulate_pcr_amplicon,
)

# Three Golden Gate parts: BsaI sites give overhangs ATGC, GCTT and CGTA.
GG_PARTS = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]
FWD = "ACGTTGCAAGGCTTAGCAT"
REV = "GGTCCAATGCTAAGCCTTGC"
DUPLEX = "ACGTTGCAAGGCTTAGCATG"


def _rc(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# DnaChisel --------------------------------------------------------------------


def test_generate_constrained_matches_upstream_and_constraints() -> None:
    dnachisel = pytest.importorskip("dnachisel")

    result = generate_constrained_dna_sequence(
        30,
        constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
        seed=1,
        verbose=False,
    )
    direct = dnachisel.random_compatible_dna_sequence(
        sequence_length=30,
        constraints=[dnachisel.AvoidPattern("GGTCTC")],
        seed=1,
        logger=None,
    )

    assert result["sequence"] == str(direct)
    assert len(result["sequence"]) == 30
    assert "GGTCTC" not in result["sequence"]
    assert result["all_constraints_pass"] is True
    assert result["constraint_names"] == ["AvoidPattern"]


def test_generate_constrained_gc_window_is_hand_checkable() -> None:
    pytest.importorskip("dnachisel")

    result = generate_constrained_dna_sequence(
        40,
        constraints=[{"name": "EnforceGCContent", "mini": 0.4, "maxi": 0.6}],
        seed=7,
        verbose=False,
    )
    gc = sum(base in "GC" for base in result["sequence"])

    assert 16 <= gc <= 24  # 40% to 60% of 40 bases
    assert result["all_constraints_pass"] is True


def test_generate_constrained_rejects_bad_input() -> None:
    pytest.importorskip("dnachisel")

    with pytest.raises(ConfigurationError, match="Unsupported DnaChisel constraint"):
        generate_constrained_dna_sequence(
            20, constraints=[{"name": "NotAConstraint"}], verbose=False
        )
    with pytest.raises(ConfigurationError, match="sequence_length"):
        generate_constrained_dna_sequence(True, verbose=False)


def test_find_sequence_difference_segments_hand_calculated() -> None:
    pytest.importorskip("dnachisel")

    result = find_sequence_difference_segments("ACGTTACGTA", "ACGAAACGTC", verbose=False)

    # Positions 3 and 4 differ (T/A, T/A) and merge; position 9 differs (A/C).
    assert result["segments"] == [[3, 5], [9, 10]]
    assert result["segment_count"] == 2
    assert result["differing_positions"] == 3


def test_find_sequence_difference_segments_requires_equal_length() -> None:
    pytest.importorskip("dnachisel")

    with pytest.raises(ConfigurationError, match="same length"):
        find_sequence_difference_segments("ACGT", "ACG", verbose=False)


def test_list_common_restriction_enzymes_matches_upstream() -> None:
    dnachisel = pytest.importorskip("dnachisel")

    result = list_common_restriction_enzymes(site_length=6, min_suppliers=5, verbose=False)
    direct = sorted(
        str(name) for name in dnachisel.list_common_enzymes(site_length=(6,), min_suppliers=5)
    )

    assert result["enzymes"] == direct
    assert result["enzyme_count"] == len(direct) == 37
    assert result["enzymes"] == sorted(result["enzymes"])


def test_group_nearby_positions_hand_calculated() -> None:
    pytest.importorskip("dnachisel")

    result = group_nearby_positions([30, 1, 11, 2, 10, 3], max_gap=2, verbose=False)

    # Gaps of 1 join; the gap 3 to 10 and 11 to 30 split the groups.
    assert result["groups"] == [[1, 2, 3], [10, 11], [30]]
    assert result["group_count"] == 3
    assert group_nearby_positions([3, 4], max_gap=2, verbose=False)["groups"] == [[3, 4]]
    assert group_nearby_positions([3, 5], max_gap=2, verbose=False)["groups"] == [[3], [5]]
    spread = group_nearby_positions([1, 2, 3, 10, 11, 30], max_group_spread=3, verbose=False)
    assert spread["groups"] == [[1, 2, 3], [10, 11], [30]]


def test_group_nearby_positions_rejects_duplicates() -> None:
    pytest.importorskip("dnachisel")

    with pytest.raises(ConfigurationError, match="unique"):
        group_nearby_positions([1, 1, 2], max_gap=2, verbose=False)


# DnaCauldron ------------------------------------------------------------------


def test_autoselect_golden_gate_enzyme_matches_upstream() -> None:
    dnacauldron = pytest.importorskip("dnacauldron")
    from dnacauldron.biotools.autoselect_enzyme import autoselect_enzyme

    result = autoselect_golden_gate_enzyme(GG_PARTS, enzymes=["BsaI", "BbsI"], verbose=False)
    records = []
    for index, text in enumerate(GG_PARTS, start=1):
        record = dnacauldron.sequence_to_biopython_record(
            text, id=f"part_{index}", name=f"part_{index}"
        )
        record.annotations["topology"] = "linear"
        records.append(record)

    assert result["enzyme"] == str(autoselect_enzyme(records, enzymes=("BsaI", "BbsI")))
    assert result["enzyme"] == "BsaI"
    assert result["part_count"] == 3


def test_autoselect_golden_gate_enzyme_rejects_unknown_enzyme() -> None:
    pytest.importorskip("dnacauldron")

    with pytest.raises(ConfigurationError, match="must be one of"):
        autoselect_golden_gate_enzyme(GG_PARTS, enzymes=["EcoRI"], verbose=False)


def test_list_type2s_digestion_overhangs_hand_calculated() -> None:
    pytest.importorskip("dnacauldron")

    result = list_type2s_digestion_overhangs(GG_PARTS, enzyme="BsaI", verbose=False)

    # Junctions: part 1 ATGC/GCTT, part 2 GCTT/CGTA, part 3 CGTA/ATGC.
    assert result["overhangs"] == ["ATGC", "CGTA", "GCTT"]
    assert result["overhang_count"] == 3
    assert result["none_placeholder_count"] == 0


def test_crop_dna_record_with_features_hand_calculated() -> None:
    pytest.importorskip("dnacauldron")
    sequence = "ACGT" * 50
    features = [
        {"type": "CDS", "start": 0, "end": 40, "label": "a"},
        {"type": "misc_feature", "start": 60, "end": 100},
        {"type": "misc_feature", "start": 150, "end": 180},
    ]

    result = crop_dna_record_with_features(
        sequence, start=10, end=120, features=features, verbose=False
    )

    assert result["sequence"] == sequence[10:120]
    assert result["sequence_length"] == 110
    # The 0-40 feature is clipped to 10-40 and shifted by 10; 150-180 lies outside.
    assert result["features"] == [
        {"type": "CDS", "start": 0, "end": 30, "label": "a"},
        {"type": "misc_feature", "start": 50, "end": 90, "label": ""},
    ]


# pydna ------------------------------------------------------------------------


def test_simulate_pcr_amplicon_hand_calculated() -> None:
    pytest.importorskip("pydna")
    template = "AAAAA" + FWD + "C" * 16 + _rc(REV) + "TTTTT"

    result = simulate_pcr_amplicon(template, FWD, REV, verbose=False)

    # 19 forward + 16 spacer + 20 reverse-primer bases.
    assert result["product_length"] == 19 + 16 + 20
    assert result["product_sequence"] == FWD + "C" * 16 + _rc(REV)


def test_simulate_pcr_amplicon_without_product_raises() -> None:
    pytest.importorskip("pydna")

    with pytest.raises(ConfigurationError, match="could not amplify"):
        simulate_pcr_amplicon("GGGGGGGGGGGGGGGGGGGGGGGG", FWD, REV, verbose=False)


def test_assemble_gibson_recovers_original_sequence() -> None:
    pytest.importorskip("pydna")
    original = (
        "ACGTTGCAAGGCTTAGCATCGATCGGATCCAAGCTTGGCCATGGAATTCGCGGCCGCTCGAGTCTAGA"
        "GGTACCGTCGACTGCAGCGGCCGCTCTAGAGGCATGCCCGGGTATGCATGCAAGCTTAAAGGAGCTC"
    )
    left, right = original[:80], original[50:]

    result = assemble_pydna_fragments([left, right], method="gibson", min_overlap=25, verbose=False)

    # The 30-base overlap (positions 50-79) is the only join; one product returns.
    assert result["product_count"] == 1
    assert result["products"][0]["sequence"] == original
    assert result["products"][0]["length"] == len(original)


def test_assemble_golden_gate_hand_calculated() -> None:
    pytest.importorskip("pydna")

    result = assemble_pydna_fragments(
        GG_PARTS, method="golden_gate", enzymes=["BsaI"], circular_only=True, verbose=False
    )

    expected = (
        "ATGC" + "TTTAAACCCGGGTTTT" + "GCTT" + "AAAACCCCGGGGAAAA" + "CGTA" + "CCCCAAAATTTTGGGG"
    )
    assert result["product_count"] == 1
    assert result["products"] == [{"sequence": expected, "length": 60, "circular": True}]


@pytest.mark.parametrize(
    ("method", "function_name"),
    [
        ("in_fusion", "in_fusion_assembly"),
        ("fusion_pcr", "fusion_pcr_assembly"),
        ("in_vivo", "in_vivo_assembly"),
    ],
)
def test_assemble_overlap_methods_match_upstream(method: str, function_name: str) -> None:
    assembly = pytest.importorskip("pydna.assembly2")
    from pydna.dseqrecord import Dseqrecord

    original = "ACGTTGCAAGGCTTAGCATCGATCGGATCCAAGCTTGGCCATGGAATTCGCGGCCGCTCGAGTC"
    fragments = [original[:50], original[20:]]
    direct = getattr(assembly, function_name)([Dseqrecord(text) for text in fragments], limit=25)
    result = assemble_pydna_fragments(fragments, method=method, min_overlap=25, verbose=False)

    assert sorted(item["sequence"] for item in result["products"]) == sorted(
        str(item.seq).upper() for item in direct
    )


def test_assemble_restriction_ligation_matches_upstream() -> None:
    assembly = pytest.importorskip("pydna.assembly2")
    restriction = pytest.importorskip("Bio.Restriction")
    from pydna.dseqrecord import Dseqrecord

    direct = assembly.restriction_ligation_assembly(
        [Dseqrecord(text) for text in GG_PARTS], enzymes=[restriction.BsaI], circular_only=True
    )
    result = assemble_pydna_fragments(
        GG_PARTS, method="restriction_ligation", enzymes=["BsaI"], circular_only=True, verbose=False
    )

    assert sorted(item["sequence"] for item in result["products"]) == sorted(
        str(item.seq).upper() for item in direct
    )
    assert result["enzymes"] == ["BsaI"]


def test_assemble_rejects_mixed_method_arguments() -> None:
    pytest.importorskip("pydna")

    with pytest.raises(ConfigurationError, match="enzymes applies only"):
        assemble_pydna_fragments(GG_PARTS[:2], method="gibson", enzymes=["BsaI"], verbose=False)
    with pytest.raises(ConfigurationError, match="method must be one of"):
        assemble_pydna_fragments(GG_PARTS[:2], method="unknown", verbose=False)


def test_find_common_substrings_hand_calculated() -> None:
    pytest.importorskip("pydna")
    x, y = "ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC"

    result = find_common_substrings(x, y, min_length=5, verbose=False)

    # x[3:19] == y[0:16] == "TTGCAAGGCTTAGCAT"
    assert result["matches"] == [{"x_start": 3, "y_start": 0, "length": 16}]
    assert find_common_substrings(x, y, min_length=17, verbose=False)["match_count"] == 0


def test_find_terminal_overlaps_hand_calculated() -> None:
    pytest.importorskip("pydna")
    x, y = "ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC"

    result = find_terminal_overlaps(x, y, min_length=5, verbose=False)

    assert result["matches"] == [{"x_start": 3, "y_start": 0, "length": 16}]
    assert result["match_count"] == 1


def test_primer_pair_hybridization_hand_calculated() -> None:
    pytest.importorskip("pydna")
    forward = "GGGGG" + DUPLEX
    reverse = _rc(DUPLEX) + "CCCCC"

    result = primer_pair_hybridization(forward, reverse, minimal_annealing=10, verbose=False)

    # The primers share a 20-base duplex with no 5' overhang; the product is the forward primer.
    assert result["overhangs"] == [0]
    assert result["products"] == [{"sequence": forward, "length": 25}]


def test_primer_melting_temperature_matches_upstream() -> None:
    tm_module = pytest.importorskip("pydna.tm")

    result = primer_melting_temperature(DUPLEX, verbose=False)

    assert result["tm_celsius"] == round(float(tm_module.tm_default(DUPLEX)), 3)
    assert result["length"] == 20


def test_parse_primer_sequences_from_file_and_text(tmp_path: Path) -> None:
    parser = pytest.importorskip("pydna.parsers")
    text = ">p1\nACGTACGTACGT\n>p2\nGGGGCCCCAAAA\n"
    path = _write(tmp_path, "primers.fa", text)

    from_file = parse_primer_sequences(path, verbose=False)
    from_text = parse_primer_sequences(text, verbose=False)

    assert from_file["primers"] == [
        {"name": "p1", "sequence": "ACGTACGTACGT", "length": 12},
        {"name": "p2", "sequence": "GGGGCCCCAAAA", "length": 12},
    ]
    assert from_text["primers"] == from_file["primers"]
    assert [item.seq for item in parser.parse_primers(text)] == [
        "ACGTACGTACGT",
        "GGGGCCCCAAAA",
    ]
    assert from_file["source"] == "primers.fa"


def test_fetch_genbank_record_is_mocked(monkeypatch: pytest.MonkeyPatch) -> None:
    genbank_module = pytest.importorskip("pydna.genbank")
    dseqrecord = pytest.importorskip("pydna.dseqrecord").Dseqrecord
    calls: list[tuple[str, str]] = []

    def fake_genbank(accession: str, *args: object, email: str | None = None, **kwargs: object):
        calls.append((accession, str(email)))
        return dseqrecord("ACGTACGTAA", name="TESTREC")

    monkeypatch.setattr(genbank_module, "genbank", fake_genbank)
    result = fetch_genbank_record("NC_000913.3", email="reader@example.org", verbose=False)

    assert calls == [("NC_000913.3", "reader@example.org")]
    assert result["accession"] == "NC_000913.3"
    assert result["length"] == 10
    assert result["sequence"] == "ACGTACGTAA"
    assert "reader@example.org" not in str(result)  # the contact address is not echoed


def test_fetch_genbank_record_validates_inputs() -> None:
    pytest.importorskip("pydna")

    with pytest.raises(ConfigurationError, match="accession"):
        fetch_genbank_record("not an accession", email="reader@example.org", verbose=False)
    with pytest.raises(ConfigurationError, match="email"):
        fetch_genbank_record("CS570233.1", email="no-at-sign", verbose=False)


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="live NCBI download; set DNAKIT_LIVE_TESTS=1 and DNAKIT_NCBI_EMAIL",
)
def test_fetch_genbank_record_live() -> None:
    email = os.environ.get("DNAKIT_NCBI_EMAIL", "")
    result = fetch_genbank_record("CS570233.1", email=email, verbose=False)
    assert result["length"] > 0


# Missing packages and public contract -----------------------------------------


@pytest.mark.parametrize(
    ("module_name", "extra", "call"),
    [
        ("dnachisel", "bio-dnachisel", lambda: list_common_restriction_enzymes(verbose=False)),
        (
            "dnacauldron",
            "bio-dnacauldron",
            lambda: list_type2s_digestion_overhangs(GG_PARTS, verbose=False),
        ),
        (
            "pydna",
            "bio-assembly",
            lambda: simulate_pcr_amplicon("ACGT" * 10, FWD, REV, verbose=False),
        ),
    ],
)
def test_missing_package_raises_backend_unavailable(
    monkeypatch: pytest.MonkeyPatch, module_name: str, extra: str, call
) -> None:
    for name in list(sys.modules):
        if name == module_name or name.startswith(module_name + "."):
            monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, module_name, None)

    with pytest.raises(BackendUnavailableError, match=extra):
        call()


def test_public_functions_take_native_keyword_only_verbose() -> None:
    for module in (dnachisel_api, dnacauldron_api, pydna_api):
        for name in module.__all__:
            signature = inspect.signature(getattr(module, name))
            parameter = signature.parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
            assert parameter.default is True, name
