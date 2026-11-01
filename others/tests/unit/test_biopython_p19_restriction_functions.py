"""Tests for the P19 Bio.Restriction and Bio.CAPS functions (no network)."""

from __future__ import annotations

import warnings

import pytest

pytest.importorskip("Bio")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.sequence import biopython_p19_restriction_functions as m  # noqa: E402

SEQUENCE = "GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC"
ENZYMES = ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"]


def _restriction():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio import Restriction

    return Restriction


def test_enzyme_properties_match_biopython() -> None:
    result = m.restriction_enzyme_properties(["EcoRI", "SmaI", "KpnI", "BsaI"], verbose=False)
    by_name = {item["name"]: item for item in result["enzymes"]}
    restriction = _restriction()
    assert by_name["EcoRI"]["elucidate"] == restriction.EcoRI.elucidate() == "G^AATT_C"
    assert by_name["EcoRI"]["end_type"] == "5_overhang" and by_name["EcoRI"]["overhang_sequence"] == "AATT"
    assert by_name["SmaI"]["end_type"] == "blunt"
    assert by_name["KpnI"]["end_type"] == "3_overhang"
    assert by_name["EcoRI"]["frequency"] == restriction.EcoRI.frequency() == 4096.0
    assert by_name["EcoRI"]["commercial"] is True
    assert "New England Biolabs" in by_name["EcoRI"]["suppliers"]
    assert by_name["BsaI"]["palindromic"] is False and by_name["BsaI"]["size"] == 6
    with pytest.raises(ConfigurationError):
        m.restriction_enzyme_properties(["NotAnEnzyme"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.restriction_enzyme_properties([], verbose=False)
    with pytest.raises(ConfigurationError):
        m.restriction_enzyme_properties(["EcoRI; rm"], verbose=False)


def test_related_enzymes() -> None:
    restriction = _restriction()
    iso = m.restriction_related_enzymes("KpnI", "isoschizomers", verbose=False)
    assert iso["enzymes"] == sorted(str(e) for e in restriction.KpnI.isoschizomers()) == ["Acc65I", "Asp718I"]
    neo = m.restriction_related_enzymes("SmaI", "neoschizomers", verbose=False)
    assert set(neo["enzymes"]) >= {"XmaI"}
    equi = m.restriction_related_enzymes("KpnI", "equischizomers", verbose=False)
    assert equi["enzymes"] == []
    compat = m.restriction_related_enzymes("EcoRI", "compatible_end", verbose=False)
    assert compat["count"] == len(restriction.EcoRI.compatible_end()) > 10
    among = m.restriction_related_enzymes("SmaI", "neoschizomers", among=["XmaI", "SmaI", "KpnI"], verbose=False)
    assert among["enzymes"] == ["XmaI"]
    with pytest.raises(ConfigurationError):
        m.restriction_related_enzymes("EcoRI", "cousins", verbose=False)


def test_suppliers() -> None:
    table = m.restriction_enzymes_by_supplier(verbose=False)
    codes = {item["code"]: item for item in table["suppliers"]}
    assert codes["N"]["name"] == "New England Biolabs" and codes["N"]["enzyme_count"] > 100
    by_code = m.restriction_enzymes_by_supplier("N", verbose=False)
    by_name = m.restriction_enzymes_by_supplier("New England Biolabs", verbose=False)
    assert by_code["enzymes"] == by_name["enzymes"] and "EcoRI" in by_code["enzymes"]
    assert by_code["enzyme_count"] == codes["N"]["enzyme_count"]
    with pytest.raises(ConfigurationError):
        m.restriction_enzymes_by_supplier("Nope", verbose=False)


def test_filter_enzymes() -> None:
    restriction = _restriction()
    blunt = m.filter_restriction_enzymes(end_type="blunt", verbose=False)
    assert blunt["count"] == len([e for e in restriction.AllEnzymes if e.is_blunt()]) > 100
    assert blunt["candidate_count"] == len(restriction.AllEnzymes)
    chosen = m.filter_restriction_enzymes(
        ENZYMES, end_type="5_overhang", palindromic=True, min_site_size=6, verbose=False
    )
    assert chosen["enzymes"] == ["BamHI", "EcoRI", "HindIII"]
    sold = m.filter_restriction_enzymes(commercial=True, supplier="N", max_site_size=4, defined_site=True, verbose=False)
    assert sold["count"] > 0 and all(restriction.AllEnzymes.get(n).size <= 4 for n in sold["enzymes"])
    assert m.filter_restriction_enzymes(ENZYMES, end_type="3_overhang", verbose=False)["enzymes"] == ["KpnI"]
    with pytest.raises(ConfigurationError):
        m.filter_restriction_enzymes(end_type="sticky", verbose=False)
    with pytest.raises(ConfigurationError):
        m.filter_restriction_enzymes(supplier="?", verbose=False)


def test_analysis_report_matches_biopython() -> None:
    restriction = _restriction()
    from Bio.Seq import Seq

    reference = restriction.Analysis(restriction.RestrictionBatch([restriction.EcoRI, restriction.BamHI, restriction.SmaI, restriction.KpnI, restriction.HindIII]), Seq(SEQUENCE), linear=True).full()
    result = m.restriction_analysis_report(SEQUENCE, ENZYMES, view="list", verbose=False)
    assert result["cuts"] == {str(k): v for k, v in reference.items() if v}
    assert result["no_cut"] == ["HindIII"]
    assert "EcoRI      :  2, 20." in result["report"]
    as_map = m.restriction_analysis_report(SEQUENCE, ENZYMES, view="map", verbose=False)
    assert "20 EcoRI" in as_map["report"] and SEQUENCE in as_map["report"]
    number = m.restriction_analysis_report(SEQUENCE, ENZYMES, view="number", verbose=False)
    assert "enzymes which cut 2 times" in number["report"]
    circular = m.restriction_analysis_report(SEQUENCE, ENZYMES, linear=False, verbose=False)
    assert circular["linear"] is False
    with pytest.raises(ConfigurationError):
        m.restriction_analysis_report(SEQUENCE, ENZYMES, view="pie", verbose=False)
    with pytest.raises(ConfigurationError):
        m.restriction_analysis_report("GAATTCXX", ENZYMES, verbose=False)


def test_filter_results_each_criterion() -> None:
    def run(**kwargs):
        return m.filter_restriction_results(SEQUENCE, ENZYMES, verbose=False, **kwargs)["results"]

    assert run(end_type="blunt") == {"SmaI": [27]}
    assert run(end_type="5_overhang") == {"BamHI": [11], "EcoRI": [2, 20], "HindIII": []}  # non-cutters stay, as in Biopython
    assert run(end_type="5_overhang", cuts="cutting") == {"BamHI": [11], "EcoRI": [2, 20]}
    assert run(end_type="3_overhang") == {"KpnI": [38]}
    assert run(cut_count=2) == {"EcoRI": [2, 20]}
    assert run(cut_count=[1, 2], end_type="5_overhang") == {"BamHI": [11], "EcoRI": [2, 20]}
    assert run(cuts="non_cutting") == {"HindIII": []}
    assert set(run(cuts="cutting")) == {"BamHI", "EcoRI", "KpnI", "SmaI"}
    assert run(names=["EcoRI", "KpnI"]) == {"EcoRI": [2, 20], "KpnI": [38]}
    assert run(site_sizes=[6], position_filter="between", start=1, end=20) == {"BamHI": [11], "EcoRI": [2, 20]}
    assert run(position_filter="show_only_between", start=15, end=38) == {"EcoRI": [20], "SmaI": [27]}
    assert run(position_filter="do_not_cut", start=15, end=38) == {"BamHI": [11], "HindIII": [], "KpnI": [38]}
    assert run(position_filter="only_outside", start=15, end=38) == {"BamHI": [11], "KpnI": [38]}
    assert run(position_filter="outside", start=15, end=38) == {"BamHI": [11], "EcoRI": [2, 20], "KpnI": [38]}
    assert run(end_type="blunt", cut_count=5) == {}  # an empty step must not widen the result again
    with pytest.raises(ConfigurationError):
        run()
    with pytest.raises(ConfigurationError):
        run(position_filter="between")
    with pytest.raises(ConfigurationError):
        run(position_filter="sideways", start=1, end=2)


def test_caps_differential_cut_sites_match_biopython() -> None:
    import warnings as warn

    from Bio.Align import MultipleSeqAlignment
    from Bio.CAPS import CAPSMap
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    alignment = {"a": "AAGAATTCGGG", "b": "AAGAATTAGGG", "c": "AAGAATTCGGG"}
    result = m.caps_differential_cut_sites(alignment, ["EcoRI"], verbose=False)
    restriction = _restriction()
    with warn.catch_warnings():
        warn.simplefilter("ignore")
        msa = MultipleSeqAlignment([SeqRecord(Seq(s), id=n) for n, s in alignment.items()])
        reference = CAPSMap(msa, [restriction.EcoRI]).dcuts
    assert [(s["start"], s["cuts_in"], s["blocked_in"]) for s in result["sites"]] == [
        (d.start, ["a", "c"], ["b"]) for d in reference
    ]
    assert m.caps_differential_cut_sites({"a": "AAGAATTCGGG", "b": "AAGAATTCGGG"}, ["EcoRI"], verbose=False)["site_count"] == 0
    with pytest.raises(ConfigurationError):
        m.caps_differential_cut_sites({"a": "AAGAATTC-GGG", "b": "AAGAATTCAGGG"}, ["EcoRI"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.caps_differential_cut_sites({"a": "AAGAATTC", "b": "AAGAATT"}, ["EcoRI"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.caps_differential_cut_sites({"a": "AAGAATTC"}, ["EcoRI"], verbose=False)
