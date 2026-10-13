"""Contract tests for the G18n review wrappers (pydna, DnaCauldron and DnaChisel).

Each wrapper is checked against the upstream call it wraps, or against a hand calculation.
All tests run offline; no live services are used.
"""

from __future__ import annotations

import random

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import dnacauldron_review_functions as dc
from dnakit.molecular_lab import dnachisel_review_functions as dchisel
from dnakit.molecular_lab import pydna_review_functions as pyd

A_SEQ = (
    "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGT"
    "GGAGAGGGTGAAGG"
)
B_SEQ = (
    "GGCGAAGGGATGCTTACCGTTGGTACACACCCCAATCTTTTACTGGAGCGTAAAGGTTACGCATGGCAATTAAATCGCTGGATTGTCGAT"
    "CCACGCCATCTGG"
)
INS1 = "ATGGCTAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGG"
INS2 = "ATGGCTAGCGTGATCAAGCCGGAAATGCGCATCCACCTGCGTAAAAGTGTTCCAGATGGCTACCACCCTGGG"


def _rc(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _part(left: str, right: str, insert: str) -> str:
    return "CGTCTCA" + left + insert + right + "A" + "GAGACG"


def _random_dna(seed: int, length: int) -> str:
    rng = random.Random(seed)
    text = "".join(rng.choice("ACGT") for _ in range(length))
    return text.replace("CGTCTC", "CGTCAC").replace("GAGACG", "GAGATG")


# pydna


def test_dna_molecules_equal_matches_pydna_and_hand_cases() -> None:
    pytest.importorskip("pydna")
    from pydna.utils import eq

    result = pyd.dna_molecules_equal(["ATGC", "GCAT"])
    assert result["equal"] is True
    assert result["equal"] == bool(eq("ATGC", "GCAT", circular=False))
    assert pyd.dna_molecules_equal(["ATGC", "TGCA"])["equal"] is False
    rotated = pyd.dna_molecules_equal(["ATGCAAGTCC", "AAGTCCATGC"], circular=True)
    assert rotated["equal"] is True
    assert pyd.dna_molecules_equal(["ATGCAAGTCC", "AAGTCCATGC"])["equal"] is False


def test_rare_codon_positions_hand_calculation() -> None:
    pytest.importorskip("pydna")
    # Codons: ATG CGA TAG AGA AAA AAA AGC TCG; only CGA is rare in the yeast table.
    result = pyd.rare_codon_positions("ATGCGATAGAGAAAAAAAAGCTCGA")
    assert result["codon_count"] == 8
    assert result["positions"] == [{"start": 3, "end": 6, "codon": "CGA"}]
    with pytest.raises(ConfigurationError):
        pyd.rare_codon_positions("ATG", organism="human")


def test_smallest_circular_rotation_hand_calculation() -> None:
    pytest.importorskip("pydna")
    result = pyd.smallest_circular_rotation("CAGT")
    assert result["sequence"] == "AGTC"
    assert result["length"] == 4


def test_generate_random_orf_is_reproducible_and_restores_random_state() -> None:
    pytest.importorskip("pydna")
    random.seed(99)
    before = random.getstate()
    first = pyd.generate_random_orf(10, seed=1)
    assert random.getstate() == before
    second = pyd.generate_random_orf(10, seed=1)
    assert first["sequence"] == second["sequence"]
    sequence = first["sequence"]
    assert len(sequence) == 30
    assert sequence.startswith("ATG")
    assert sequence[-3:] in {"TAA", "TAG", "TGA"}
    codons = [sequence[index : index + 3] for index in range(0, len(sequence) - 3, 3)]
    assert not set(codons) & {"TAA", "TAG", "TGA"}


def test_annotate_primer_binding_sites_hand_calculation() -> None:
    pytest.importorskip("pydna")
    product = "ATGCGTACGTTAGCATCGGATCCAAAGGGCCC"
    result = pyd.annotate_primer_binding_sites(product, "ATGCGTACGTTAGCAT", "GGGCCCTTT")
    assert result["product_length"] == 32
    assert result["features"] == [
        {"type": "primer_bind", "start": 0, "end": 16, "strand": 1, "label": "forward_primer"},
        {"type": "primer_bind", "start": 23, "end": 32, "strand": -1, "label": "reverse_primer"},
    ]
    with pytest.raises(ConfigurationError):
        pyd.annotate_primer_binding_sites(product, "TTTTTTTTTTTTTTTT", "GGGCCCTTT")


def test_suggest_pcr_program_matches_pydna_programme() -> None:
    pytest.importorskip("pydna")
    from pydna.amplify import pcr
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from pydna.tm import dbd_program, program

    forward = A_SEQ[:22]
    reverse = _rc(A_SEQ[-22:])
    amplicon = pcr(Primer(forward), Primer(reverse), Dseqrecord(A_SEQ), limit=14)
    taq = pyd.suggest_pcr_program(A_SEQ, forward, reverse)
    assert taq["program"] == program(amplicon)
    assert taq["product_length"] == len(amplicon)
    assert taq["program_lines"] == taq["program"].splitlines()
    dbd = pyd.suggest_pcr_program(A_SEQ, forward, reverse, polymerase="dbd")
    assert dbd["program"] == dbd_program(amplicon)
    assert taq["program"] != dbd["program"]
    with pytest.raises(ConfigurationError):
        pyd.suggest_pcr_program(A_SEQ, forward, reverse, polymerase="phi29")


def test_ligate_blunt_fragments_gives_circular_products_of_total_length() -> None:
    pytest.importorskip("pydna")
    first = "ATGCGTACGTTAGCATCGGATCCAAAGGGCCCTTTAAAGGAGATCGTTTAGCGCCC"
    second = "TTTAGCGCCCAAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGC"
    result = pyd.ligate_blunt_fragments([first, second])
    assert result["product_count"] == 2
    assert all(row["circular"] for row in result["products"])
    assert {row["length"] for row in result["products"]} == {len(first) + len(second)}


def test_gateway_recombination_products_matches_pydna() -> None:
    pytest.importorskip("pydna")
    from Bio.Data.IUPACData import ambiguous_dna_values
    from pydna.assembly2 import gateway_assembly
    from pydna.dseqrecord import Dseqrecord
    from pydna.gateway import recombinase_dict

    site = recombinase_dict["BP"]["conservative"].recombinases[0]

    def expand(text: str) -> str:
        return "".join(
            ambiguous_dna_values[ch][0] if ch in ambiguous_dna_values else ch for ch in text.upper()
        )

    att_b = expand(site.site1)
    att_p = expand(site.site2)
    linear = att_b + "TTTTACGTACGTAAAGGGCCCAAATTTGGGGCCAAGGCCTTAAGGCCAA" + att_b
    plasmid = att_p + "AAAACCCGGGTTTAAACCCGGG" + att_p
    result = pyd.gateway_recombination_products([linear, plasmid])
    upstream = gateway_assembly([Dseqrecord(linear), Dseqrecord(plasmid)], "BP", greedy=False)
    assert result["product_count"] == len(upstream) > 0
    with pytest.raises(ConfigurationError):
        pyd.gateway_recombination_products([linear, plasmid], reaction_type="XX")


def test_simulate_genome_edit_cre_lox_excision_matches_pydna_docstring() -> None:
    pytest.importorskip("pydna")
    from pydna.cre_lox import LOXP_SEQUENCE

    integrated = (
        "CCCCCCATAACTTCGTATAGCATACATTATACGAAGTTATGGGGGATAACTTCGTATAGCATACATTATACGAAGTTATAAAAA"
    )
    result = pyd.simulate_genome_edit(
        integrated, mechanism="cre_lox", operation="excision", limit=40
    )
    assert result["limit"] is None
    assert sorted(row["length"] for row in result["products"]) == [39, 45]
    assert LOXP_SEQUENCE in integrated
    integration = pyd.simulate_genome_edit(
        "CCCCCC" + LOXP_SEQUENCE + "AAAAA",
        mechanism="cre_lox",
        operation="integration",
        inserts=[LOXP_SEQUENCE + "GGGGG"],
        inserts_circular=True,
    )
    assert [row["length"] for row in integration["products"]] == [84]


def test_simulate_genome_edit_homologous_recombination_matches_pydna() -> None:
    pytest.importorskip("pydna")
    from pydna.assembly2 import homologous_recombination_excision
    from pydna.dseqrecord import Dseqrecord

    repeat = "ACGTACGTAGGCTTAACCGGTTAGCCTAGGATCGATCGGA"
    genome = (
        "CCCAAATTTGGGCCCAAATTT"
        + repeat
        + "GGGTTTAAACCCGGGTTTAAA"
        + repeat
        + "TTTGGGAAACCCTTTGGGAAA"
    )
    result = pyd.simulate_genome_edit(
        genome, mechanism="homologous_recombination", operation="excision", limit=40
    )
    upstream = homologous_recombination_excision(Dseqrecord(genome), limit=40)
    assert result["limit"] == 40
    assert sorted(row["length"] for row in result["products"]) == sorted(len(p) for p in upstream)
    assert result["product_count"] == 2


def test_simulate_genome_edit_crispr_docstring_example() -> None:
    pytest.importorskip("pydna")
    result = pyd.simulate_genome_edit(
        "aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac",
        mechanism="crispr",
        operation="integration",
        inserts=["aaccggttAAAAAAAAAttcaaagcac"],
        guides=["ttcaatgcaaacagtaatga"],
        limit=8,
    )
    assert [row["length"] for row in result["products"]] == [27]


def test_simulate_genome_edit_recombinase_excision_and_validation() -> None:
    pytest.importorskip("pydna")
    sites = ("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
    genome = "TTTTATGCCCTAAaaTTGGGGGGGGGGGGGAAaaTTTTTTTCCCTCCCC"
    result = pyd.simulate_genome_edit(
        genome, mechanism="recombinase", operation="excision", recombinase_sites=sites
    )
    assert sorted(row["length"] for row in result["products"]) == [19, 30]
    with pytest.raises(ConfigurationError):
        pyd.simulate_genome_edit(genome, mechanism="recombinase", operation="excision")
    with pytest.raises(ConfigurationError):
        pyd.simulate_genome_edit(genome, mechanism="cre_lox", operation="integration")
    with pytest.raises(ConfigurationError):
        pyd.simulate_genome_edit(genome, mechanism="crispr", operation="excision")


def test_assemble_recombinase_fragments_matches_pydna() -> None:
    pytest.importorskip("pydna")
    from pydna.assembly2 import recombinase_assembly
    from pydna.dseqrecord import Dseqrecord
    from pydna.recombinase import Recombinase

    sites = ("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")
    fragments = ["CCCCATGCCCTAAaaTTGGGGGGGG", "AAaaTTTTTTTCCCTTTTTTTTT"]
    result = pyd.assemble_recombinase_fragments(fragments, recombinase_sites=sites)
    upstream = recombinase_assembly(
        [Dseqrecord(text) for text in fragments], Recombinase(*sites), circular_only=False
    )
    assert result["product_count"] == len(upstream) == 2
    assert sorted(row["length"] for row in result["products"]) == sorted(len(p) for p in upstream)


def test_extract_protospacers_finds_spacer_before_scaffold() -> None:
    pytest.importorskip("pydna")
    from pydna.crispr import cas9

    spacer = "GACGTACGATCGTAGCTAGG"
    construct = "ATTA" + spacer + cas9.scaffold + "TTTT"
    result = pyd.extract_protospacers(construct)
    assert result["protospacers"] == [spacer]
    assert result["protospacer_length"] == 20
    assert pyd.extract_protospacers("ACGT" * 5)["protospacer_count"] == 0


def test_parse_sequence_text_matches_pydna_parser() -> None:
    pytest.importorskip("pydna")
    from pydna.parsers import embl_gb_fasta

    text = ">one\nATGCCC\n>two\nGGGAAA\n"
    result = pyd.parse_sequence_text(text)
    upstream = embl_gb_fasta(text)
    assert [row["name"] for row in result["records"]] == [record.name for record in upstream]
    assert [row["sequence"] for row in result["records"]] == [
        str(record.seq).upper() for record in upstream
    ]
    with pytest.raises(ConfigurationError):
        pyd.parse_sequence_text("   ")


def test_design_user_cloning_primers_matches_pydna_design() -> None:
    pytest.importorskip("pydna")
    from pydna.amplify import pcr
    from pydna.design import user_assembly_design
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    fragments = [
        {"template": A_SEQ, "forward_primer": A_SEQ[:22], "reverse_primer": _rc(A_SEQ[-22:])},
        {"template": B_SEQ, "forward_primer": B_SEQ[:22], "reverse_primer": _rc(B_SEQ[-22:])},
    ]
    result = pyd.design_user_cloning_primers(fragments)
    amplicons = [
        pcr(
            Primer(item["forward_primer"]),
            Primer(item["reverse_primer"]),
            Dseqrecord(item["template"]),
            limit=14,
        )
        for item in fragments
    ]
    upstream = user_assembly_design(amplicons, max_overlap=15, min_overlap=4, max_tail=50)
    assert result["fragment_count"] == 2
    assert [row["forward_primer"] for row in result["fragments"]] == [
        str(item.forward_primer.seq) for item in upstream
    ]
    with pytest.raises(ConfigurationError):
        pyd.design_user_cloning_primers([fragments[0]])


# DnaCauldron


def test_simulate_type2s_assembly_matches_dnacauldron_mix() -> None:
    pytest.importorskip("dnacauldron")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from dnacauldron.AssemblyMix import generate_type2s_restriction_mix
    from dnacauldron.biotools import set_record_topology

    parts = {"p1": _part("AATG", "GCTT", INS1), "p2": _part("GCTT", "AATG", INS2)}
    result = dc.simulate_type2s_assembly(parts, enzyme="BsmBI")
    records = []
    for name, text in parts.items():
        record = SeqRecord(Seq(text), id=name, name=name, annotations={"molecule_type": "DNA"})
        set_record_topology(record, "linear")
        records.append(record)
    upstream = list(
        generate_type2s_restriction_mix(parts=records, enzyme="BsmBI").compute_circular_assemblies()
    )
    assert result["circular_product_count"] == len(upstream) == 1
    assert result["products"][0]["sequence"] == str(upstream[0].seq).upper()
    assert result["products"][0]["length"] == len(upstream[0])
    auto = dc.simulate_type2s_assembly(parts)
    assert auto["enzyme"] == "BsmBI"
    with pytest.raises(ConfigurationError):
        dc.simulate_type2s_assembly(parts, enzyme="EcoRI")


def test_record_has_backbone_hand_lengths() -> None:
    pytest.importorskip("dnacauldron")
    backbone = _random_dna(7, 700)
    donor = _part("AATG", "GCTT", INS1) + backbone
    result = dc.record_has_backbone(donor)
    assert result["construct_length"] == len(donor) == 794
    assert result["insert_length"] == 72
    assert result["backbone_length"] == 722
    assert result["has_backbone"] is True
    assert dc.record_has_backbone(donor, min_backbone_length=800)["has_backbone"] is False


def test_substitute_type2s_overhangs_matches_dnacauldron() -> None:
    pytest.importorskip("dnacauldron")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from dnacauldron.biotools import set_record_topology
    from dnacauldron.utils.utils import substitute_overhangs

    donor = _part("AATG", "GCTT", INS1) + _random_dna(7, 700)
    result = dc.substitute_type2s_overhangs(donor, {"AATG": "CCCC"})
    record = SeqRecord(Seq(donor), id="v", name="v", annotations={"molecule_type": "DNA"})
    set_record_topology(record, "circular")
    upstream = substitute_overhangs(record, {"AATG": "CCCC"}, enzyme="BsmBI")
    assert result["product_sequence"] == str(upstream.seq).upper()
    assert result["product_length"] == 794
    assert result["applied_overhangs"] == ["AATG"]
    assert dc.substitute_type2s_overhangs(donor, {"TTTT": "CCCC"})["applied_overhangs"] == []


def test_clone_into_donor_vector_gives_backbone_plus_insert() -> None:
    pytest.importorskip("dnacauldron")
    donor = _part("AATG", "GCTT", INS1) + _random_dna(7, 700)
    insert = _part("AATG", "GCTT", INS2)
    result = dc.clone_into_donor_vector(donor, insert, enzyme="BsmBI")
    # The backbone (722 bp) plus the 72 bp insert replaces the 72 bp stuffer.
    assert result["product_length"] == 722 + 72
    assert "CCCCATGGCTAGCGTG" not in result["product_sequence"]
    with pytest.raises(ConfigurationError):
        dc.clone_into_donor_vector(donor, _part("CCCC", "GGGG", INS2), enzyme="BsmBI")


def test_assign_parts_to_backbones_reports_each_part() -> None:
    pytest.importorskip("dnacauldron")
    backbone = _part("AATG", "GCTT", _random_dna(11, 40)) + _random_dna(7, 700)
    result = dc.assign_parts_to_backbones(
        {"p1": _part("AATG", "GCTT", INS1), "p2": _part("TTAA", "CCGG", INS2)},
        {"bb": backbone},
    )
    states = {row["part"]: row["status_per_part"] for row in result["assignments"]}
    assert states == {"p1": "assigned", "p2": "no_backbone"}
    assert result["assigned_count"] == 1
    assert result["assignments"][0]["backbone"] == "bb"
    with pytest.raises(ConfigurationError):
        dc.assign_parts_to_backbones({"p1": "ACGT"}, {"bb": "ACGT"}, min_backbone_length=0)


# DnaChisel


def test_dna_window_overlap_hand_cases() -> None:
    pytest.importorskip("dnachisel")
    overlap = dchisel.dna_window_overlap([0, 10], [5, 20])
    assert overlap["overlap"] == [5, 10]
    assert overlap["overlap_length"] == 5
    assert overlap["overlaps"] is True
    touching = dchisel.dna_window_overlap([0, 5], [5, 9])
    assert touching["overlap_length"] == 0
    assert touching["overlaps"] is False
    assert dchisel.dna_window_overlap([0, 3], [5, 9])["overlap"] is None
    with pytest.raises(ConfigurationError):
        dchisel.dna_window_overlap([4, 2], [0, 9])


def test_subdivide_dna_window_hand_cases() -> None:
    pytest.importorskip("dnachisel")
    result = dchisel.subdivide_dna_window([0, 10], max_span=4)
    assert result["windows"] == [[0, 4], [4, 8], [8, 10]]
    assert dchisel.subdivide_dna_window([3, 3], max_span=2)["windows"] == []
    with pytest.raises(ConfigurationError):
        dchisel.subdivide_dna_window([0, 10], max_span=0)


def test_group_nearby_dna_segments_hand_case() -> None:
    pytest.importorskip("dnachisel")
    segments = [[0, 5], [2, 4], [10, 12], [11, 20]]
    result = dchisel.group_nearby_dna_segments(segments, max_start_gap=3)
    assert result["groups"] == [[[0, 5], [2, 4]], [[10, 12], [11, 20]]]
    assert result["group_count"] == 2


def test_codon_backtranslation_table_matches_standard_code() -> None:
    pytest.importorskip("dnachisel")
    result = dchisel.codon_backtranslation_table()
    table = result["table"]
    assert table["M"] == ["ATG"]
    assert table["*"] == ["TAA", "TAG", "TGA"]
    assert len(table["X"]) == 61
    assert result["entry_count"] == 26
    with pytest.raises(ConfigurationError):
        dchisel.codon_backtranslation_table("No such table")
