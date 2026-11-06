"""Tests for the pydna functions of pack P34 (dnakit.molecular_lab.pydna_p34_functions)."""

from __future__ import annotations

import random
import struct

import pytest

pytest.importorskip("pydna")

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import pydna_p34_functions as m

Q = {"verbose": False}
STICKY = {"watson": "aggctaa", "crick": "tagcc"}
TEMPLATE = "ttcgatcgATGCCCAAATTTGGGCCCTAGGGATCCTTTAAAGGGCCC"
FWD, REV = "ATGCCCAAATTTGGG", "GGGCCCTTTAAAGGATCC"


def _reference(spec):
    from pydna.dseq import Dseq

    return Dseq(spec["watson"], spec["crick"], ovhg=spec.get("ovhg"))


def _same(result, expected):
    assert result["watson"] == str(expected.watson)
    assert result["crick"] == str(expected.crick)
    assert result["ovhg"] == expected.ovhg


def test_fill_in_matches_pydna_and_restricts_nucleotides():
    spec = {"watson": "caaa", "crick": "cttt"}
    out = m.pydna_fill_in_ends(spec, **Q)
    _same(out["result"], _reference(spec).fill_in())
    assert out["result"]["is_blunt"] is True
    only_a = m.pydna_fill_in_ends(spec, nucleotides="a", **Q)
    _same(only_a["result"], _reference(spec).fill_in("a"))
    with pytest.raises(ConfigurationError):
        m.pydna_fill_in_ends(spec, nucleotides="x", **Q)


def test_mung_bean_t4_and_cast_to_ds_agree_with_pydna():
    ref = _reference(STICKY)
    _same(m.pydna_mung_bean_blunt(STICKY, **Q)["result"], ref.nibble_to_blunt())
    _same(m.pydna_t4_polymerase_polish(STICKY, **Q)["result"], ref.T4())
    both = m.pydna_cast_to_double_strand(STICKY, **Q)
    _same(both["result"], ref.cast_to_ds())
    _same(m.pydna_cast_to_double_strand(STICKY, end="left", **Q)["result"], ref.cast_to_ds_left())
    _same(m.pydna_cast_to_double_strand(STICKY, end="right", **Q)["result"], ref.cast_to_ds_right())
    assert m.pydna_mung_bean_blunt(STICKY, **Q)["result"]["is_blunt"] is True


@pytest.mark.parametrize("end", ["left", "right"])
def test_exonuclease_resection_matches_pydna(end):
    ref = _reference({"watson": "aggctaa", "crick": "ttagcct"})
    five = m.pydna_exonuclease_five_prime("aggctaa", end=end, bases=2, **Q)["result"]
    _same(five, getattr(ref, f"nibble_five_prime_{end}")(2))
    three = m.pydna_exonuclease_three_prime("aggctaa", end=end, bases=2, **Q)["result"]
    _same(three, getattr(ref, f"nibble_three_prime_{end}")(2))
    with pytest.raises(ConfigurationError):
        m.pydna_exonuclease_five_prime("aggctaa", end="middle", **Q)


def test_terminal_transferase_and_user_treatment():
    from pydna.dseq import Dseq

    out = m.pydna_terminal_transferase_tail("gctgaa", **Q)
    _same(out["result"], Dseq("gctgaa").terminal_transferase("a"))
    assert out["result"]["watson"] == "gctgaaa"
    gc = m.pydna_terminal_transferase_tail("gctgaa", nucleotides="gg", **Q)
    assert gc["result"]["watson"] == "gctgaagg"
    user = m.pydna_user_enzyme_treatment("tagaagtaggUat", **Q)
    _same(user["result"], Dseq("tagaagtaggUat").user())
    assert user["result"]["watson"] == "tagaagtagg at"


def test_melt_short_strands_finds_the_gap_and_matches_pydna():
    spec = {"watson": "tagaagtagg at", "crick": "catacctacttcta", "ovhg": 0}
    out = m.pydna_melt_short_strands(spec, length=3, **Q)
    ref = _reference(spec)
    assert out["watson_sites"] == [[int(a), int(b)] for a, b in ref.get_ss_meltsites(3)[0]]
    assert out["shed"][0]["watson"] == "at"
    remaining, _strands = ref.melt_ss_dna(3)
    assert out["remaining"]["watson"] == str(remaining.watson)
    assert len(out["melt_fragments"]) == len(ref.melt(3))
    clean = m.pydna_melt_short_strands("tagaagtaggtatg", length=3, **Q)
    assert clean["shed"] == [] and clean["watson_sites"] == []


def test_shed_strand_segments_matches_pydna():
    from pydna.dseq import Dseq

    out = m.pydna_shed_strand_segments("tagaagtaggtatg", watson_cutpairs=[[6, 8]], **Q)
    new, strands = Dseq("tagaagtaggtatg").shed_ss_dna([(6, 8)], [])
    assert out["result"]["watson"] == str(new.watson) == "tagaag  ggtatg"
    assert [item["watson"] for item in out["shed"]] == [str(s.watson) for s in strands] == ["ta"]
    with pytest.raises(ConfigurationError):
        m.pydna_shed_strand_segments("tagaagtaggtatg", watson_cutpairs=[[1]], **Q)


def test_restriction_cutters_classes():
    seq = "GGATCCaaaGAATTC"
    names = ["BamHI", "EcoRI", "HindIII"]
    once = m.pydna_restriction_cutters(seq, cutter_class="once", enzymes=names, **Q)
    assert once["enzymes"] == ["BamHI", "EcoRI"]
    assert once["cut_counts"] == {"BamHI": 1, "EcoRI": 1}
    assert m.pydna_restriction_cutters(seq, cutter_class="none", enzymes=names, **Q)["enzymes"] == [
        "HindIII"
    ]
    twice = m.pydna_restriction_cutters(seq + "GGATCC", cutter_class="twice", enzymes=names, **Q)
    assert twice["enzymes"] == ["BamHI"]
    anyc = m.pydna_restriction_cutters(seq, cutter_class="any", enzymes=names, **Q)
    assert anyc["enzymes"] == ["BamHI", "EcoRI"]
    three = m.pydna_restriction_cutters(
        seq + "GGATCC" * 2, cutter_class="n", cuts=3, enzymes=names, **Q
    )
    assert three["enzymes"] == ["BamHI"]
    from pydna.dseq import Dseq

    default = m.pydna_restriction_cutters(seq, cutter_class="once", **Q)
    assert set(default["enzymes"]) == {str(e) for e in Dseq(seq).unique_cutters()}
    with pytest.raises(ConfigurationError):
        m.pydna_restriction_cutters(seq, enzymes=["NoSuchEnzyme"], **Q)


def test_restriction_ligation_assembly_products_and_partial_digest_flags():
    vector = {"sequence": "ttGGATCCaaaaaGAATTCtt", "name": "vec"}
    insert = {"sequence": "GGATCCgggggGAATTC", "name": "ins"}
    out = m.pydna_restriction_ligation_assembly([vector, insert], ["BamHI", "EcoRI"], **Q)
    assert out["product_count"] == len(out["products"]) >= 1
    circular = [p for p in out["products"] if p["circular"]]
    assert circular and circular[0]["length"] == 22
    assert all(isinstance(p["partially_digested_inputs"], list) for p in out["products"])
    only = m.pydna_restriction_ligation_assembly(
        [vector, insert], ["BamHI", "EcoRI"], circular_only=True, **Q
    )
    assert all(p["circular"] for p in only["products"])
    with pytest.raises(ConfigurationError):
        m.pydna_restriction_ligation_assembly([vector], ["Nope"], **Q)


def test_pcr_alignment_figure_equals_pydna_amplicon_figure():
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    out = m.pydna_pcr_alignment_figure(TEMPLATE, FWD, REV, **Q)
    amp = pcr(Primer(FWD), Primer(REV), Dseqrecord(TEMPLATE))
    assert out["figure"] == amp.figure()
    assert out["product"] == str(amp.seq.watson)
    assert out["product_length"] == 39
    with pytest.raises(ConfigurationError):
        m.pydna_pcr_alignment_figure(TEMPLATE, "TTTTTTTTTTTTTTTTT", REV, **Q)


FRAGMENTS = ["acgtacgtacgtaaGATCGATCGATCGgg", "GATCGATCGATCGggttccaaggttcc"]


def test_assembly_figures_compact_detailed_and_matplotlib(tmp_path):
    compact = m.pydna_assembly_figure(FRAGMENTS, limit=12, **Q)
    assert "15" in compact["figure"] and "fragment1" in compact["figure"]
    assert compact["product_length"] == 41
    detailed = m.pydna_assembly_figure(FRAGMENTS, limit=12, figure="detailed", **Q)
    assert "GATCGATCGATCG" in detailed["figure"].upper()
    pytest.importorskip("matplotlib")
    target = tmp_path / "asm.png"
    drawn = m.pydna_assembly_figure(
        FRAGMENTS, limit=12, figure="matplotlib", output_path=str(target), **Q
    )
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and drawn["image_path"] == str(target)
    with pytest.raises(ConfigurationError):
        m.pydna_assembly_figure(FRAGMENTS, limit=12, figure="matplotlib", **Q)
    with pytest.raises(ConfigurationError):
        m.pydna_assembly_figure(FRAGMENTS, limit=200, **Q)


def test_cloning_strategy_export_import_and_history_round_trip():
    built = m.pydna_build_cloning_strategy([TEMPLATE], operation="pcr", primers=[FWD, REV], **Q)
    assert built["product_names"] == ["product1"] and "PCRSource" in built["history"][0]
    records = m.pydna_cloning_strategy_to_records(built["strategy_json"], **Q)
    assert records["record_count"] == 1
    assert records["records"][0]["sequence"].upper() == "ATGCCCAAATTTGGGCCCTAGGGATCCTTTAAAGGGCCC"
    shown = m.pydna_cloning_history(built["strategy_json"], action="display", **Q)
    assert "PCRSource" in shown["records"][0]["history"] and shown["records"][0]["has_history"]
    valid = m.pydna_cloning_history(built["strategy_json"], action="validate", **Q)
    assert valid["records"][0]["valid"] is True
    norm = m.pydna_cloning_history(built["strategy_json"], action="normalize", **Q)
    assert norm["records"][0]["unchanged"] is True

    gibson = m.pydna_build_cloning_strategy(FRAGMENTS, operation="gibson", limit=12, **Q)
    assert "GibsonAssemblySource" in gibson["history"][0]
    assert m.pydna_cloning_history(gibson["strategy_json"], action="validate", **Q)["records"][0][
        "valid"
    ]

    plain = m.pydna_build_cloning_strategy([TEMPLATE], **Q)
    assert (
        m.pydna_cloning_history(plain["strategy_json"], **Q)["records"][0]["has_history"] is False
    )
    with pytest.raises(ConfigurationError):
        m.pydna_cloning_strategy_to_records("{not json", **Q)
    with pytest.raises(ConfigurationError):
        m.pydna_build_cloning_strategy([TEMPLATE], operation="restriction_ligation", **Q)


def test_dna_mass_peptide_search_and_origin_sync():
    from pydna.dseqrecord import Dseqrecord

    out = m.pydna_dna_mass_from_amount("atgcatgc", **Q)
    ref = Dseqrecord("atgcatgc")
    assert out["mass_g"] == pytest.approx(ref.m())
    assert out["molecular_weight"] == pytest.approx(ref.seq.mw())
    double = m.pydna_dna_mass_from_amount("atgcatgc", amount_mol=1e-13, **Q)
    assert double["mass_g"] == pytest.approx(2 * out["mass_g"])

    found = m.pydna_find_peptide("atgtacgatcgtatgctggttatattttag", "RML", **Q)
    assert (found["start"], found["end"], found["dna"]) == (9, 18, "cgtatgctg")
    assert m.pydna_find_peptide("atgtacgatcgtatgctggttatattttag", "WWW", **Q)["found"] is False

    synced = m.pydna_sync_circular_origin("atgaCCC", "gaat", **Q)
    assert synced["result"] == "gaCCCat" and synced["circular"] is True


def _write_abif(path, bases):
    data1, data2 = bases.encode(), bytes([30]) * len(bases)
    start = 128
    directory = start + len(data1) + len(data2)
    header = b"ABIF" + struct.pack(">H4sI2H3I", 101, b"tdir", 1, 1023, 28, 2, 56, directory)
    entries = [
        struct.pack(">4sI2H4I", b"PBAS", 2, 2, 1, len(data1), len(data1), start, 0),
        struct.pack(">4sI2H4I", b"PCON", 2, 2, 1, len(data2), len(data2), start + len(data1), 0),
    ]
    blob = header.ljust(start, b"\0") + data1 + data2 + b"".join(entries)
    path.write_bytes(blob)


def test_map_sanger_traces_on_a_synthetic_abif_file(tmp_path):
    pytest.importorskip("Bio")
    rng = random.Random(7)
    reference = "".join(rng.choice("acgt") for _ in range(120))
    _write_abif(tmp_path / "read1.ab1", reference[20:90].upper())
    _write_abif(tmp_path / "other.ab1", "".join(rng.choice("ACGT") for _ in range(60)))
    out = m.pydna_map_sanger_traces(reference, str(tmp_path / "*.ab1"), limit=20, **Q)
    assert [p.split("/")[-1] for p in out["matching_files"]] == ["read1.ab1"]
    assert (out["mapped"][0]["start"], out["mapped"][0]["end"]) == (20, 90)
    with pytest.raises(ConfigurationError):
        m.pydna_map_sanger_traces(reference, str(tmp_path / "none*.ab1"), **Q)


def test_clean_genbank_text_repairs_and_parses():
    text = (
        "LOCUS       test                      12 bp    DNA     linear   UNK 01-JAN-1980\n"
        "DEFINITION  .\nFEATURES             Location/Qualifiers\n     CDS             1..9\n"
        "                     /label=x\nORIGIN\n        1 atgaaatgat ag\n//\n"
    )
    out = m.pydna_clean_genbank_text(text, **Q)
    assert out["name"] == "test" and out["length"] == 12 and out["feature_count"] == 1
    assert "ds-DNA" in out["genbank_text"]
    with pytest.raises(ConfigurationError):
        m.pydna_clean_genbank_text("   ", **Q)


def _primer_screen_ready():
    pytest.importorskip("ahocorasick")


def _screen_data():
    from Bio.Seq import Seq

    rng = random.Random(1)
    template = "".join(rng.choice("ACGT") for _ in range(1200))
    insert = "".join(rng.choice("ACGT") for _ in range(300))
    other = template[:650] + insert + template[700:]
    rc = lambda s: str(Seq(s).reverse_complement())  # noqa: E731
    primers = [
        {"name": "F", "sequence": template[100:120]},
        {"name": "R", "sequence": rc(template[900:920])},
        {"name": "F2", "sequence": template[300:320]},
    ]
    return template, other, insert, rc, primers


def test_primer_pairs_by_size_and_flanking_pairs():
    _primer_screen_ready()
    template, _, _, _, primers = _screen_data()
    pairs = m.pydna_primer_pairs_by_size(template, primers, short=500, long=1500, **Q)
    assert {(p["forward_primer"], p["reverse_primer"], p["size"]) for p in pairs["pairs"]} == {
        ("F", "R", 820),
        ("F2", "R", 620),
    }
    narrow = m.pydna_primer_pairs_by_size(template, primers, short=700, long=1500, **Q)
    assert [p["size"] for p in narrow["pairs"]] == [820]
    flank = m.pydna_flanking_primer_pairs(template, primers, target_start=400, target_end=800, **Q)
    assert {p["size"] for p in flank["pairs"]} == {820, 620}
    tight = m.pydna_flanking_primer_pairs(template, primers, target_start=200, target_end=800, **Q)
    assert {p["size"] for p in tight["pairs"]} == {820}


def test_diagnostic_primer_pairs_and_triplets():
    _primer_screen_ready()
    template, other, insert, rc, primers = _screen_data()
    shorter = template[:500] + template[800:]
    pairs = m.pydna_diagnostic_primer_sets([template, shorter], primers, **Q)
    assert pairs["set_count"] == 1
    assert [x["product_size"] for x in pairs["sets"][0]] == [820, 520]
    strict = m.pydna_diagnostic_primer_sets(
        [template, shorter], primers, min_relative_size_difference=0.9, **Q
    )
    assert strict["set_count"] == 0
    triplet_primers = [
        {"name": "F", "sequence": template[100:120]},
        {"name": "rA", "sequence": rc(template[660:680])},
        {"name": "rB", "sequence": rc(insert[200:220])},
    ]
    trio = m.pydna_diagnostic_primer_sets(
        [template, other], triplet_primers, primers_per_set=3, short=300, long=1500, **Q
    )
    assert trio["set_count"] == 1
    assert [(x["reverse_primer"], x["product_size"]) for x in trio["sets"][0]] == [
        ("rA", 580),
        ("rB", 770),
    ]
    with pytest.raises(ConfigurationError):
        m.pydna_diagnostic_primer_sets([template], primers, **Q)


def test_recombinase_sites_found_on_both_strands():
    from pydna.dseqrecord import Dseqrecord
    from pydna.recombinase import Recombinase

    site1, site2 = "ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT"
    seq = "aaaATGCCCTAAaaTTtt"
    out = m.pydna_find_recombinase_sites(seq, site1, site2, **Q)
    expected = Recombinase(site1, site2).find(Dseqrecord(seq))
    assert out["sites"]["site1"][0]["start"] == int(expected["site1"][0].start) == 3
    assert out["annotated_features"][0]["type"] == "protein_bind"
    from Bio.Seq import Seq

    flipped = "ccc" + str(Seq("ATGCCCTAAAATT").reverse_complement()) + "ggg"
    reverse = m.pydna_find_recombinase_sites(flipped, site1, site2, **Q)
    assert reverse["sites"]["site1"][0]["strand"] == -1
    assert m.pydna_find_recombinase_sites("aaaaaaaaaaaa", site1, site2, **Q)["sites"] == {}
    with pytest.raises(ConfigurationError):
        m.pydna_find_recombinase_sites(seq, "ATGCCCTAAaaTT", "CCCCCCCCC", **Q)


def test_terminal_codon_usage_and_expression_report():
    from pydna import codon

    out = m.pydna_terminal_codon_usage("ATGAAACCCTAA", **Q)
    assert out["start_value"] == codon.start["sce"]["ATG"] == 1.0
    assert out["stop_value"] == codon.stop["sce"]["TAA"] == 0.47
    only = m.pydna_terminal_codon_usage("TTTAAACCCTGA", codon="start", **Q)
    assert only["start_value"] is None and "stop_value" not in only
    with pytest.raises(ConfigurationError):
        m.pydna_terminal_codon_usage("ATG", organism="human", **Q)

    cds = "ATGAAACCCGGGTTTAAATAA"
    report = m.pydna_cds_expression_report(cds, **Q)
    assert report["codon_count"] == 7 and report["label"] == "ATG...TAA"
    assert report["gc"] == pytest.approx(7 / 21, abs=1e-3)
    assert report["start_value"] == 1.0 and report["stop_value"] == 0.47
    assert report["n_end"] == codon.n_end["sce"]["Lys"]
    assert report["rare_codon_counts"]["CGA"] == 0
    rare = m.pydna_cds_expression_report("ATGAAACGACGACGATAA", **Q)
    assert rare["rare_codon_counts"]["CGA"] == 3 and rare["rare_fraction"] == pytest.approx(0.5)
    with pytest.raises(ConfigurationError):
        m.pydna_cds_expression_report("ATGAA", **Q)


def test_open_reading_frame_check_matches_pydna():
    from pydna.seqrecord import SeqRecord

    for text in ["atgtaa", "atgaaa", "atttaa", "atgaaataa", "atgtaatag"]:
        assert m.pydna_is_open_reading_frame(text, **Q)["is_orf"] == SeqRecord(text).isorf()
    assert m.pydna_is_open_reading_frame("atgtaa", **Q)["is_orf"] is True
    assert m.pydna_is_open_reading_frame("atgaaa", **Q)["is_orf"] is False
    for table in (1, 4, 11):
        for text in ["atgtgataa", "atgtaa", "ctgaaataa"]:
            assert m.pydna_is_open_reading_frame(text, table=table, **Q)["is_orf"] == SeqRecord(
                text
            ).isorf(table)


def test_seguid_stamp_is_deterministic_and_topology_aware():
    from pydna.dseqrecord import Dseqrecord

    out = m.pydna_seguid_stamp("atgc", timestamp="2020-01-01", **Q)
    record = Dseqrecord("atgc")
    assert out["stamp"] == record.stamp(now=lambda: "2020-01-01") == out["stamp"]
    assert out["stamp"].startswith("ldseguid=")
    assert out["comment"] == f"pydna {out['stamp']} 2020-01-01"
    circular = m.pydna_seguid_stamp({"watson": "atgc", "circular": True}, timestamp="x", **Q)
    assert circular["stamp"].startswith("cdseguid=") and circular["circular"] is True
    noted = m.pydna_seguid_stamp("atgc", timestamp="t", tool="dnakit", comment="hello", **Q)
    assert noted["comment"].startswith("dnakit ") and noted["comment"].endswith("hello")


def test_iupac_site_to_regex_matches_pydna_and_search():
    import re

    from pydna.sequence_regex import compute_regex_site

    out = m.pydna_iupac_site_to_regex("GGATCCNRY", **Q)
    assert out["regex"] == compute_regex_site("GGATCCNRY") == "(?i)GGATCC[GATC][AG][CT]"
    assert re.search(out["regex"], "ttggatccaact") is not None
    assert re.search(out["regex"], "ttggatccatt") is None
    with pytest.raises(ConfigurationError):
        m.pydna_iupac_site_to_regex("GG TT", **Q)


def test_all_public_functions_have_verbose_and_registry_friendly_signatures():
    import inspect

    for name in m.__all__:
        function = getattr(m, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
