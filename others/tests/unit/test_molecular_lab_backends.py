"""Contract tests for optional molecular-lab backends: design, records, SBOL3, SnapGene, plots."""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab.backends import (
    build_feature_map,
    build_genbank_feature_map,
    build_sbol3_document,
    calculate_primer3_dimer,
    calculate_primer3_hairpin,
    calculate_primer3_tm,
    convert_snapgene_to_genbank,
    design_assembly_primers,
    design_pcr_primers,
    evaluate_dnachisel_specs,
    find_recombination_sites,
    optimize_dnachisel_record,
    read_dnacauldron_parts,
    read_genbank_genome_tracks,
    read_sbol3_sequences,
    read_snapgene_features,
    read_snapgene_sequence,
    render_feature_map_svg,
    render_genome_tracks_svg,
    reverse_translate_protein,
    validate_sbol3_document,
    write_dnacauldron_constructs,
)

LOXP = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"  # 34 bp loxP site


def _random_dna(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_genbank(path: Path, record_id: str, sequence: str, features=()) -> Path:
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqIO import write
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq(sequence), id=record_id, name=record_id, description="test part")
    record.annotations["molecule_type"] = "DNA"
    record.features.extend(
        SeqFeature(FeatureLocation(start, end, strand=strand), type=kind, qualifiers=qualifiers)
        for start, end, strand, kind, qualifiers in features
    )
    write(record, str(path), "genbank")
    return path


# DnaChisel


def test_evaluate_dnachisel_specs_reports_pattern_location_by_hand(tmp_path: Path) -> None:
    dnachisel = pytest.importorskip("dnachisel")
    sequence = "ATG" + "GGTCTC" + "ATGCCCAAA"  # GGTCTC occupies 0-based [3, 9)
    result = evaluate_dnachisel_specs(
        sequence, constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}], verbose=False
    )

    upstream = dnachisel.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[dnachisel.AvoidPattern("GGTCTC")],
        logger=None,
    )
    direct = next(iter(upstream.constraints_evaluations()))
    assert result["all_constraints_pass"] is False
    assert result["constraints"][0]["passes"] is False
    assert bool(direct.passes) is False
    assert result["constraints"][0]["locations"] == [[3, 9]]
    assert result["failed_constraints"] == [result["constraints"][0]["label"]]
    assert result["provenance"]["package"] == "dnachisel"


def test_evaluate_dnachisel_specs_rejects_unknown_specs_and_empty_input() -> None:
    pytest.importorskip("dnachisel")
    with pytest.raises(ConfigurationError, match="Unsupported DnaChisel specification"):
        evaluate_dnachisel_specs("ACGT", constraints=[{"name": "os.system"}], verbose=False)
    with pytest.raises(ConfigurationError, match="at least one constraint"):
        evaluate_dnachisel_specs("ACGT", verbose=False)
    with pytest.raises(ConfigurationError, match="DNA string"):
        evaluate_dnachisel_specs(12, constraints=[{"name": "AvoidHairpins"}], verbose=False)


def test_optimize_dnachisel_record_removes_pattern_and_keeps_length(tmp_path: Path) -> None:
    dnachisel = pytest.importorskip("dnachisel")
    source = _write_genbank(tmp_path / "part.gb", "partA", "ATG" + "GGTCTC" + "ATGCCCAAATAA")
    target = tmp_path / "optimized.gb"

    result = optimize_dnachisel_record(
        source,
        output_path=target,
        constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
        verbose=False,
    )

    assert result["sequence_length"] == 21  # 3 + 6 + 12 bases
    assert result["edits"] >= 1
    assert result["all_constraints_pass"] is True
    from Bio import SeqIO

    written = str(next(SeqIO.parse(str(target), "genbank")).seq)
    assert len(written) == 21
    assert "GGTCTC" not in written
    upstream = dnachisel.DnaOptimizationProblem(
        sequence=written, constraints=[dnachisel.AvoidPattern("GGTCTC")], logger=None
    )
    assert upstream.all_constraints_pass()


def test_optimize_dnachisel_record_refuses_to_overwrite(tmp_path: Path) -> None:
    pytest.importorskip("dnachisel")
    source = _write_genbank(tmp_path / "part.gb", "partA", "ATGGCTAAATAA")
    existing = tmp_path / "exists.gb"
    existing.write_text("occupied", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        optimize_dnachisel_record(
            source,
            output_path=existing,
            constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
            verbose=False,
        )


def test_reverse_translate_protein_uses_first_codons_by_hand() -> None:
    pytest.importorskip("dnachisel")
    # M = ATG; A first codon in table order is GCT; K = AAA; stop = TAA.
    result = reverse_translate_protein("MAK*", verbose=False)
    assert result["sequence"] == "ATGGCTAAATAA"
    assert result["sequence_length"] == 12
    assert result["randomized"] is False


def test_reverse_translate_protein_matches_upstream_and_validates_input() -> None:
    dnachisel = pytest.importorskip("dnachisel")
    protein = "MKTAYIAKQRQISFVKSHFSRQ"
    result = reverse_translate_protein(protein, verbose=False)
    assert result["sequence"] == dnachisel.reverse_translate(protein)
    assert dnachisel.translate(result["sequence"]) == protein
    with pytest.raises(ConfigurationError, match="standard residues"):
        reverse_translate_protein("MKB", verbose=False)
    with pytest.raises(ConfigurationError, match="Unknown codon table"):
        reverse_translate_protein("MK", table="NoSuchTable", verbose=False)


# pydna


def test_design_pcr_primers_spans_template_and_matches_pydna() -> None:
    pydna_design = pytest.importorskip("pydna.design")
    pytest.importorskip("pydna.dseqrecord")
    template = _random_dna(7, 300)
    result = design_pcr_primers(template, verbose=False)

    forward, reverse = result["forward_primer"], result["reverse_primer"]
    assert template.startswith(forward)
    assert template.endswith(_reverse_complement(reverse))
    assert result["product_length"] == len(template) == 300

    from pydna.dseqrecord import Dseqrecord

    direct = pydna_design.primer_design(Dseqrecord(template), target_tm=55.0, limit=13)
    assert forward == str(direct.forward_primer.seq)
    assert reverse == str(direct.reverse_primer.seq)


def test_design_pcr_primers_rejects_bad_targets() -> None:
    pytest.importorskip("pydna.design")
    with pytest.raises(ConfigurationError, match="target_tm"):
        design_pcr_primers(_random_dna(1, 100), target_tm=True, verbose=False)
    with pytest.raises(ConfigurationError, match="too short"):
        design_pcr_primers("ACGTACGT", verbose=False)


def test_design_assembly_primers_share_overlap_at_each_junction() -> None:
    pytest.importorskip("pydna.design")
    templates = [_random_dna(11, 300), _random_dna(12, 250)]
    result = design_assembly_primers(templates, overlap=20, verbose=False)

    first, second = result["fragments"]
    # The reverse primer of fragment 1 carries the 20 bases that begin fragment 2.
    assert _reverse_complement(first["reverse_primer"])[-20:] == second["forward_primer"][:20]
    assert first["forward_primer"] == templates[0][: len(first["forward_primer"])]
    assert result["overlap"] == 20


def test_design_assembly_primers_rejects_single_fragment_and_short_overlap() -> None:
    pytest.importorskip("pydna.design")
    with pytest.raises(ConfigurationError, match="between 2"):
        design_assembly_primers([_random_dna(1, 200)], verbose=False)
    with pytest.raises(ConfigurationError, match="overlap"):
        design_assembly_primers(
            [_random_dna(1, 200), _random_dna(2, 200)], overlap=5, verbose=False
        )


def test_find_recombination_sites_locates_loxp_by_hand() -> None:
    pytest.importorskip("pydna.cre_lox")
    sequence = "GGG" + LOXP + "CCC"  # loxP occupies [3, 37)
    result = find_recombination_sites(sequence, system="loxp", verbose=False)
    assert result["sites"] == [{"name": "loxP", "start": 3, "end": 37, "strand": 1}]
    assert result["sequence_length"] == 3 + 34 + 3


def test_find_recombination_sites_matches_pydna_gateway_finder() -> None:
    pytest.importorskip("pydna.gateway")
    from pydna.dseqrecord import Dseqrecord
    from pydna.gateway import find_gateway_sites

    sequence = "AAACAAGTTTGTACAAAAAAGCAGGCTCCCC"
    result = find_recombination_sites(sequence, system="gateway", verbose=False)
    direct = find_gateway_sites(Dseqrecord(sequence), greedy=False)
    expected = [int(item.start) for item in direct["attB1"]]
    assert [row["start"] for row in result["sites"] if row["name"] == "attB1"] == expected
    with pytest.raises(ConfigurationError, match="loxp"):
        find_recombination_sites(sequence, system="flp", verbose=False)


# SBOL3


SBOL_NAMESPACE = "https://dnakit.example.org/parts/"


def test_build_sbol3_document_round_trips_sequences(tmp_path: Path) -> None:
    pytest.importorskip("sbol3")
    target = tmp_path / "parts.ttl"
    built = build_sbol3_document(
        {"part_one": "ATGCGT", "part_two": "GGATCC"},
        namespace=SBOL_NAMESPACE,
        output_path=target,
        verbose=False,
    )
    assert built["component_count"] == 2
    assert built["sequence_bases"] == 12
    assert built["validation"]["error_count"] == 0
    assert built["document"] is None

    read = read_sbol3_sequences(target, verbose=False)
    by_id = {row["display_id"]: row["elements"] for row in read["sequences"]}
    assert by_id == {"part_one_sequence": "ATGCGT", "part_two_sequence": "GGATCC"}
    assert read["sequences"][0]["encoding"].endswith("format_1207")  # EDAM IUPAC DNA

    checked = validate_sbol3_document(target, verbose=False)
    assert checked["valid"] is True


def test_build_sbol3_document_returns_text_and_validates_inputs() -> None:
    pytest.importorskip("sbol3")
    built = build_sbol3_document({"part_one": "ATGCGT"}, namespace=SBOL_NAMESPACE, verbose=False)
    assert "sbol:Component" in built["document"]
    with pytest.raises(ConfigurationError, match="namespace"):
        build_sbol3_document({"part_one": "ATG"}, namespace="ftp://example.org/", verbose=False)
    with pytest.raises(ConfigurationError, match="display ID"):
        build_sbol3_document({"1bad": "ATG"}, namespace=SBOL_NAMESPACE, verbose=False)
    with pytest.raises(ConfigurationError, match="file_format"):
        build_sbol3_document(
            {"part_one": "ATG"}, namespace=SBOL_NAMESPACE, file_format="pdf", verbose=False
        )


def test_validate_sbol3_document_reports_upstream_errors(tmp_path: Path) -> None:
    sbol3 = pytest.importorskip("sbol3")
    document = sbol3.Document()
    # A Component without a type breaks the SBOL3 rule that types are required.
    document.add(sbol3.Component(SBOL_NAMESPACE + "untyped", [], namespace=SBOL_NAMESPACE))
    path = tmp_path / "untyped.ttl"
    document.write(str(path), sbol3.TURTLE)

    checked = validate_sbol3_document(path, verbose=False)
    assert checked["valid"] is False
    assert checked["error_count"] >= 1
    assert checked["errors"][0]["object"] == SBOL_NAMESPACE + "untyped"


def test_read_sbol3_sequences_rejects_non_rdf_text(tmp_path: Path) -> None:
    pytest.importorskip("sbol3")
    path = tmp_path / "junk.ttl"
    path.write_text("this is not rdf at all @@@", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="could not read"):
        read_sbol3_sequences(path, verbose=False)


# SnapGene


def test_snapgene_reader_rejects_non_snapgene_bytes(tmp_path: Path) -> None:
    pytest.importorskip("snapgene_reader")
    path = tmp_path / "not_snapgene.dna"
    path.write_text(">fasta\nACGT\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="could not parse"):
        read_snapgene_sequence(path, verbose=False)
    with pytest.raises(ConfigurationError, match="could not parse"):
        read_snapgene_features(path, verbose=False)
    with pytest.raises(ConfigurationError, match="could not convert"):
        convert_snapgene_to_genbank(path, tmp_path / "out.gb", verbose=False)
    assert not (tmp_path / "out.gb").exists()


@pytest.mark.skipif(
    not os.environ.get("DNAKIT_SNAPGENE_TEST_FILE"),
    reason="set DNAKIT_SNAPGENE_TEST_FILE to a SnapGene-written .dna file to run this test",
)
def test_snapgene_functions_match_upstream_on_real_file(tmp_path: Path) -> None:
    snapgene_reader = pytest.importorskip("snapgene_reader")
    source = Path(os.environ["DNAKIT_SNAPGENE_TEST_FILE"])
    direct = snapgene_reader.snapgene_file_to_dict(filepath=str(source))

    sequence = read_snapgene_sequence(source, verbose=False)
    assert sequence["sequence"] == direct["seq"]
    assert sequence["length"] == len(direct["seq"])
    features = read_snapgene_features(source, verbose=False)
    assert features["feature_count"] == len(direct["features"])
    converted = convert_snapgene_to_genbank(source, tmp_path / "converted.gb", verbose=False)
    assert converted["bytes"] > 0


# DnaCauldron


def test_read_dnacauldron_parts_loads_fasta_and_genbank(tmp_path: Path) -> None:
    dnacauldron = pytest.importorskip("dnacauldron")
    fasta = tmp_path / "partx.fa"
    fasta.write_text(">partX\nATGCGTACGTTAGC\n", encoding="utf-8")
    genbank = _write_genbank(tmp_path / "partb.gb", "partB", "GGATCCGGATCCGGATCCGG")

    result = read_dnacauldron_parts([fasta, genbank], verbose=False)

    by_id = {part["id"]: part for part in result["parts"]}
    assert set(by_id) == {"partX", "partB"}
    assert by_id["partX"]["length"] == 14
    assert by_id["partX"]["sequence"] == "ATGCGTACGTTAGC"
    assert by_id["partB"]["length"] == 20
    direct = dnacauldron.load_record(str(fasta), topology="default_to_linear")
    assert str(direct.seq).upper() == by_id["partX"]["sequence"]


def test_read_dnacauldron_parts_rejects_duplicate_ids(tmp_path: Path) -> None:
    pytest.importorskip("dnacauldron")
    first = _write_genbank(tmp_path / "one.gb", "same", "ATGC")
    second = _write_genbank(tmp_path / "two.gb", "same", "GGCC")
    with pytest.raises(ConfigurationError, match="Duplicate part ID"):
        read_dnacauldron_parts([first, second], verbose=False)


def test_write_dnacauldron_constructs_writes_circular_gibson_product(tmp_path: Path) -> None:
    pytest.importorskip("dnacauldron")
    rng_a = _random_dna(3, 120)
    shared_1 = _random_dna(4, 30)
    middle = _random_dna(5, 100)
    shared_2 = _random_dna(6, 30)
    tail = _random_dna(7, 90)
    # The three parts overlap in a circle, so the Gibson product is circular.
    parts = {
        "A": rng_a + shared_1,
        "B": shared_1 + middle + shared_2,
        "C": shared_2 + tail + rng_a[:30],
    }
    result = write_dnacauldron_constructs(
        parts, output_dir=tmp_path / "constructs", method="gibson", verbose=False
    )

    assert len(result["constructs"]) == 1
    construct = result["constructs"][0]
    assert construct["length"] == 120 + 30 + 100 + 30 + 90 == 370
    assert construct["topology"] == "circular"
    from Bio import SeqIO

    written = next(SeqIO.parse(str(tmp_path / "constructs" / construct["file"]), "genbank"))
    expected = rng_a + shared_1 + middle + shared_2 + tail
    assert len(written.seq) == 370
    assert expected in str(written.seq) * 2  # a circular product may start at any base


def test_write_dnacauldron_constructs_refuses_existing_directory(tmp_path: Path) -> None:
    pytest.importorskip("dnacauldron")
    (tmp_path / "taken").mkdir()
    with pytest.raises(ConfigurationError, match="already exists"):
        write_dnacauldron_constructs(
            {"A": "ATGCGT", "B": "GCGTAA"}, output_dir=tmp_path / "taken", verbose=False
        )


# DNA Features Viewer


def test_build_feature_map_assigns_lanes_by_hand() -> None:
    pytest.importorskip("dna_features_viewer")
    features = [
        {"start": 0, "end": 40, "label": "a"},
        {"start": 30, "end": 60, "strand": -1, "label": "b"},
        {"start": 70, "end": 90, "label": "c"},
    ]
    result = build_feature_map(100, features, verbose=False)
    # a and b overlap, so b moves to lane 1; c overlaps nothing and stays on lane 0.
    assert [row["lane"] for row in result["features"]] == [0.0, 1.0, 0.0]
    assert result["lane_count"] == 2
    assert result["features"][1]["strand"] == -1
    assert result["features"][0]["color"] == "#000080"


def test_build_feature_map_matches_upstream_levels() -> None:
    dnafv_levels = pytest.importorskip("dna_features_viewer.compute_features_levels")
    from dna_features_viewer import GraphicFeature

    features = [{"start": 0, "end": 40, "label": "a"}, {"start": 30, "end": 60, "label": "b"}]
    graphic = [
        GraphicFeature(start=row["start"], end=row["end"], label=row["label"]) for row in features
    ]
    direct = dnafv_levels.compute_features_levels(graphic)
    result = build_feature_map(100, features, verbose=False)
    assert [row["lane"] for row in result["features"]] == [direct[item] for item in graphic]


def test_build_feature_map_rejects_bad_coordinates() -> None:
    pytest.importorskip("dna_features_viewer")
    with pytest.raises(ConfigurationError, match="feature end"):
        build_feature_map(50, [{"start": 10, "end": 5}], verbose=False)
    with pytest.raises(ConfigurationError, match="strand"):
        build_feature_map(50, [{"start": 1, "end": 5, "strand": True}], verbose=False)
    with pytest.raises(ConfigurationError, match="#RRGGBB"):
        build_feature_map(50, [{"start": 1, "end": 5, "color": "red"}], verbose=False)


def test_render_feature_map_svg_is_deterministic(tmp_path: Path) -> None:
    pytest.importorskip("dna_features_viewer")
    features = [{"start": 0, "end": 40, "label": "a"}, {"start": 30, "end": 60, "strand": -1}]
    first = render_feature_map_svg(100, features, output_path=tmp_path / "a.svg", verbose=False)
    second = render_feature_map_svg(100, features, output_path=tmp_path / "b.svg", verbose=False)
    content = (tmp_path / "a.svg").read_bytes()
    assert content == (tmp_path / "b.svg").read_bytes()
    assert b"<svg" in content
    assert first["bytes"] == second["bytes"] == len(content)
    with pytest.raises(ConfigurationError, match=r"must end with \.svg"):
        render_feature_map_svg(100, features, output_path=tmp_path / "x.png", verbose=False)


def test_build_genbank_feature_map_returns_every_feature(tmp_path: Path) -> None:
    pytest.importorskip("dna_features_viewer")
    path = _write_genbank(
        tmp_path / "map.gb",
        "map1",
        "A" * 100,
        features=[
            (10, 40, 1, "CDS", {"gene": ["geneX"]}),
            (30, 60, -1, "gene", {"gene": ["geneY"]}),
        ],
    )
    result = build_genbank_feature_map(path, verbose=False)
    spans = sorted((row["start"], row["end"]) for row in result["features"])
    assert spans == [(10, 40), (30, 60)]
    assert result["sequence_length"] == 100
    assert result["lane_count"] == 2  # the two features overlap


# pyGenomeViz


def _genes_genbank(path: Path) -> Path:
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqIO import write
    from Bio.SeqRecord import SeqRecord

    sequence = "ATG" + "A" * 97 + "TAA" + "C" * 200  # 303 bp
    record = SeqRecord(Seq(sequence), id="contig1", name="contig1", description="genes")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(FeatureLocation(0, 100, strand=1), type="CDS", qualifiers={"gene": ["geneA"]}),
        SeqFeature(
            FeatureLocation(150, 260, strand=-1),
            type="CDS",
            qualifiers={"locus_tag": ["L2"], "product": ["hyp"]},
        ),
        SeqFeature(FeatureLocation(10, 50, strand=1), type="gene", qualifiers={}),
    ]
    write(record, str(path), "genbank")
    return path


def test_read_genbank_genome_tracks_selects_cds_with_labels(tmp_path: Path) -> None:
    pgv_genbank = pytest.importorskip("pygenomeviz.parser.genbank")
    path = _genes_genbank(tmp_path / "genes.gb")

    result = read_genbank_genome_tracks(path, feature_type="CDS", verbose=False)

    assert result["track_count"] == 1
    track = result["tracks"][0]
    assert track["length"] == 303
    assert track["features"] == [
        {"type": "CDS", "start": 0, "end": 100, "strand": 1, "label": "geneA"},
        {"type": "CDS", "start": 150, "end": 260, "strand": -1, "label": "L2"},
    ]
    direct = pgv_genbank.Genbank(str(path)).get_seqid2size()
    assert direct == {"contig1": 303}


def test_render_genome_tracks_svg_draws_adjacent_links_only(tmp_path: Path) -> None:
    pytest.importorskip("pygenomeviz")
    tracks = [
        {
            "name": "a",
            "length": 300,
            "features": [{"start": 10, "end": 100, "strand": 1, "label": "g1"}],
        },
        {"name": "b", "length": 250, "features": [{"start": 5, "end": 90, "strand": -1}]},
        {"name": "c", "length": 200, "features": []},
    ]
    link = {"track_a": "a", "start_a": 20, "end_a": 90, "track_b": "b", "start_b": 10, "end_b": 80}
    result = render_genome_tracks_svg(
        tracks[:2], output_path=tmp_path / "tracks.svg", links=[link], verbose=False
    )
    assert result["track_count"] == 2
    assert result["link_count"] == 1
    assert result["feature_count"] == 2
    assert b"<svg" in (tmp_path / "tracks.svg").read_bytes()

    with pytest.raises(ConfigurationError, match="adjacent"):
        render_genome_tracks_svg(
            tracks,
            output_path=tmp_path / "bad.svg",
            links=[{**link, "track_b": "c", "start_b": 10, "end_b": 80}],
            verbose=False,
        )
    assert not (tmp_path / "bad.svg").exists()


# primer3-py


def test_calculate_primer3_tm_matches_direct_binding() -> None:
    primer3 = pytest.importorskip("primer3")
    primer = "ACGTTGCAAGCTTGCA"
    result = calculate_primer3_tm(primer, verbose=False)
    direct = primer3.calc_tm(primer, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)
    assert result["tm_c"] == pytest.approx(direct)
    assert result["sequence_length"] == 16
    shifted = calculate_primer3_tm(primer, mv_conc=100.0, verbose=False)["tm_c"]
    assert shifted > result["tm_c"]  # more monovalent salt stabilises the duplex


def test_calculate_primer3_hairpin_reports_structure_and_none_when_absent() -> None:
    primer3 = pytest.importorskip("primer3")
    hairpin = calculate_primer3_hairpin("GGGGAAAACCCC", verbose=False)
    direct = primer3.calc_hairpin("GGGGAAAACCCC")
    assert hairpin["structure_found"] is True
    assert hairpin["tm_c"] == pytest.approx(direct.tm)
    assert hairpin["dg_cal_per_mol"] == pytest.approx(direct.dg)

    absent = calculate_primer3_hairpin("AAAAAAAAAAAA", verbose=False)
    assert absent["structure_found"] is False
    assert absent["tm_c"] is None
    assert absent["dg_cal_per_mol"] is None


def test_calculate_primer3_dimer_covers_homo_and_hetero_dimers() -> None:
    primer3 = pytest.importorskip("primer3")
    homo = calculate_primer3_dimer("ATGCGCAT", verbose=False)
    assert homo["dimer_type"] == "homodimer"
    assert homo["partner_length"] is None
    assert homo["dg_cal_per_mol"] == pytest.approx(primer3.calc_homodimer("ATGCGCAT").dg)

    hetero = calculate_primer3_dimer("ATGCGTACGT", partner="ACGTACGCAT", verbose=False)
    assert hetero["dimer_type"] == "heterodimer"
    assert hetero["partner_length"] == 10
    direct = primer3.calc_heterodimer("ATGCGTACGT", "ACGTACGCAT")
    assert hetero["tm_c"] == pytest.approx(direct.tm)


def test_primer3_functions_validate_sequences_and_conditions() -> None:
    pytest.importorskip("primer3")
    with pytest.raises(ConfigurationError, match="only A, C, G and T"):
        calculate_primer3_tm("ACGN", verbose=False)
    with pytest.raises(ConfigurationError, match="exceeds 200 bases"):
        calculate_primer3_tm("A" * 201, verbose=False)
    with pytest.raises(ConfigurationError, match="mv_conc"):
        calculate_primer3_tm("ACGT", mv_conc=-1.0, verbose=False)
    with pytest.raises(ConfigurationError, match="max_loop"):
        calculate_primer3_hairpin("ACGTACGT", max_loop=True, verbose=False)


# Missing optional packages


@pytest.mark.parametrize(
    ("module", "extra", "call"),
    [
        (
            "dnachisel",
            "bio-dnachisel",
            lambda p: evaluate_dnachisel_specs(
                "ACGT", constraints=[{"name": "AvoidHairpins"}], verbose=False
            ),
        ),
        ("dnacauldron", "bio-dnacauldron", lambda p: read_dnacauldron_parts([p], verbose=False)),
        (
            "pydna",
            "bio-assembly",
            lambda p: design_pcr_primers("ACGTACGTACGTACGTACGTACGT", verbose=False),
        ),
        ("sbol3", "bio-pysbol3", lambda p: validate_sbol3_document(p, verbose=False)),
        ("snapgene_reader", "bio-snapgene", lambda p: read_snapgene_sequence(p, verbose=False)),
        (
            "dna_features_viewer",
            "bio-dna-features-viewer",
            lambda p: build_feature_map(10, [{"start": 0, "end": 5}], verbose=False),
        ),
        ("primer3", "bio-primer3", lambda p: calculate_primer3_tm("ACGTACGTAC", verbose=False)),
        (
            "pygenomeviz",
            "bio-pygenomeviz",
            lambda p: render_genome_tracks_svg(
                [{"name": "a", "length": 10}], output_path=p.with_suffix(".svg"), verbose=False
            ),
        ),
    ],
)
def test_missing_optional_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, module: str, extra: str, call
) -> None:
    placeholder = tmp_path / "input.dat"
    placeholder.write_text("placeholder", encoding="utf-8")
    # Submodules imported by earlier tests are cached, so block every loaded name in the package.
    for name in [key for key in sys.modules if key == module or key.startswith(module + ".")]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.raises(BackendUnavailableError, match=extra):
        call(placeholder)
