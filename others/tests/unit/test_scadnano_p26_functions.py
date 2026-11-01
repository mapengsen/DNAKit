"""Tests for dnakit.genomics.scadnano_p26_functions (needs scadnano; skips when missing)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.genomics import scadnano_p26_functions as p26


def _design() -> dict:
    """Two helices: a 48-base scaffold crossing over and a named staple."""

    sc = pytest.importorskip("scadnano")
    design = sc.Design(helices=[sc.Helix(max_offset=48) for _ in range(2)], grid=sc.Grid.square)
    design.draw_strand(0, 0).to(24).cross(1).move(-24).as_scaffold()
    design.draw_strand(1, 0).to(24).with_name("stapA").with_vendor_fields("25nm", "STD", "P1", "A1")
    return json.loads(design.to_json())


def _model(data: dict):
    sc = pytest.importorskip("scadnano")
    return sc.Design.from_scadnano_json_map(json.loads(json.dumps(data)))


@pytest.fixture(autouse=True)
def _need_scadnano(request: pytest.FixtureRequest) -> None:
    if request.node.name != "test_rejects_non_mapping_design_without_backend":
        pytest.importorskip("scadnano")


@pytest.fixture()
def design() -> dict:
    return _design()


@pytest.fixture()
def sequenced(design: dict) -> dict:
    return p26.assign_dna_sequence(design, "ACGT" * 12, strand_index=0, verbose=False)["design"]


def test_rejects_non_mapping_design_without_backend() -> None:
    with pytest.raises(ConfigurationError):
        p26.compute_base_pairs("not a design", verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        p26.load_scadnano_design(verbose=False)


def test_load_three_sources(design: dict, tmp_path: Path) -> None:
    by_map = p26.load_scadnano_design(design, verbose=False)
    by_text = p26.load_scadnano_design(json_text=json.dumps(design), verbose=False)
    target = tmp_path / "d.sc"
    target.write_text(json.dumps(design))
    by_file = p26.load_scadnano_design(path=target, verbose=False)
    assert by_map["strands"] == by_text["strands"] == by_file["strands"] == 2
    assert [by_map["source_kind"], by_text["source_kind"], by_file["source_kind"]] == [
        "json_map",
        "json_string",
        "file",
    ]
    with pytest.raises(ConfigurationError):
        p26.load_scadnano_design(design, json_text="{}", verbose=False)


def test_export_scadnano_matches_to_json(design: dict, tmp_path: Path) -> None:
    out = tmp_path / "x.sc"
    result = p26.export_scadnano_design(design, output_path=out, verbose=False)
    assert result["json_text"] == _model(design).to_json()
    assert out.read_text() == result["json_text"]
    with pytest.raises(ConfigurationError):
        p26.export_scadnano_design(design, output_path=out, verbose=False)


def test_cadnano_round_trip(tmp_path: Path) -> None:
    rectangle = p26.create_origami_rectangle(num_helices=4, num_cols=10, verbose=False)["design"]
    exported = p26.export_cadnano_v2_design(rectangle, verbose=False)
    assert _model(rectangle).to_cadnano_v2_json() == exported["json_text"]
    imported = p26.import_cadnano_v2_design(exported["cadnano_v2"], verbose=False)
    assert imported["helices"] == 4
    file = tmp_path / "c.json"
    file.write_text(exported["json_text"])
    assert p26.import_cadnano_v2_design(file, verbose=False)["strands"] == imported["strands"]


def test_oxview_export(sequenced: dict, tmp_path: Path) -> None:
    out = tmp_path / "d.oxview"
    result = p26.export_oxview_design(sequenced, output_path=out, verbose=False)
    assert result["strands"] == 2
    assert result["nucleotides"] == 72
    assert json.loads(out.read_text())["systems"]


def test_idt_bulk_and_plate(sequenced: dict, tmp_path: Path) -> None:
    result = p26.export_idt_order(sequenced, verbose=False)
    assert result["text"] == _model(sequenced).to_idt_bulk_input_format()
    assert result["text"].startswith("stapA,")
    pytest.importorskip("openpyxl")
    out = tmp_path / "plate.xlsx"
    plate = p26.export_idt_order(
        sequenced, order_format="plate_excel", output_path=out, verbose=False
    )
    assert out.is_file() and plate["bytes"] > 0
    with pytest.raises(ConfigurationError):
        p26.export_idt_order(sequenced, order_format="plate_excel", verbose=False)


def test_plate_maps(design: dict) -> None:
    pytest.importorskip("tabulate")
    result = p26.build_plate_maps(design, verbose=False)
    assert result["plates"][0]["plate_name"] == "P1"
    assert result["plates"][0]["wells"] == {"A1": "stapA"}
    assert "stapA" in result["plates"][0]["table"]


def test_sequences_and_complement(design: dict, sequenced: dict) -> None:
    model = _model(sequenced)
    assert model.strands[1].dna_sequence == "ACGT" * 6
    assert p26.get_strand_sequence(sequenced, 1, verbose=False)["sequence"] == "ACGT" * 6
    delimited = p26.get_strand_sequence(sequenced, 0, delimiter="|", verbose=False)["sequence"]
    assert delimited == model.strands[0].dna_sequence_delimited("|")
    vendor = p26.get_strand_sequence(sequenced, 0, vendor_format=True, verbose=False)["sequence"]
    assert vendor == model.strands[0].vendor_dna_sequence()
    sub = p26.get_domain_subsequence(sequenced, 0, 0, 2, 5, verbose=False)["sequence"]
    assert sub == model.strands[0].domains[0].dna_sequence_in(2, 5)
    with pytest.raises(ConfigurationError):
        p26.get_strand_sequence(sequenced, 9, verbose=False)


def test_assign_complement_from_strand(design: dict) -> None:
    first = p26.assign_dna_sequence(
        design, "ACGT" * 12, strand_index=0, assign_complement=False, verbose=False
    )["design"]
    result = p26.assign_complement_sequence(first, 1, 0, verbose=False)
    assert result["assigned"]["dna_sequence"] == "ACGT" * 6
    with pytest.raises(ConfigurationError):
        p26.assign_dna_sequence(design, "A" * 5, strand_index=0, check_length=True, verbose=False)


def test_assign_helix_sequences(design: dict) -> None:
    seq = p26.assign_dna_sequence(
        design, "ACGT" * 12, strand_index=0, assign_complement=False, verbose=False
    )["design"]
    result = p26.assign_helix_domain_sequences(seq, 1, helix=0, forward=False, verbose=False)
    assert result["assigned"]["dna_sequence"] == "ACGT" * 6


def test_complementarity() -> None:
    pytest.importorskip("scadnano")
    assert p26.check_complementarity("A", "T", verbose=False)["complementary"]
    assert not p26.check_complementarity("A", "G", verbose=False)["complementary"]
    assert p26.check_complementarity("A", "?", allow_wildcard=True, verbose=False)["complementary"]
    assert p26.check_complementarity("ACG", "TGC", mode="sequences", verbose=False)["complementary"]
    with pytest.raises(ConfigurationError):
        p26.check_complementarity("AC", "T", verbose=False)
    with pytest.raises(ConfigurationError):
        p26.check_complementarity("A", "T", mode="other", verbose=False)


def test_m13(design: dict) -> None:
    sc = pytest.importorskip("scadnano")
    result = p26.get_m13_sequence(verbose=False)
    assert result["sequence"] == sc.m13() and result["length"] == 7249
    assert p26.get_m13_sequence(rotation=0, variant="p8064", verbose=False)["length"] == 8064
    assigned = p26.assign_m13_scaffold(design, verbose=False)
    scaffold = _model(assigned["design"]).strands[0]
    assert scaffold.dna_sequence == sc.m13()[:48]
    with pytest.raises(ConfigurationError):
        p26.get_m13_sequence(variant="p1", verbose=False)


def test_modifications(sequenced: dict) -> None:
    added = p26.set_strand_modification(
        sequenced, 1, location="5p", display_text="B", vendor_code="/5Biosg/", verbose=False
    )
    assert _model(added["design"]).strands[1].modification_5p.vendor_code == "/5Biosg/"
    listed = p26.list_design_modifications(added["design"], verbose=False)
    assert listed["count"] == 1
    assert (
        p26.list_design_modifications(added["design"], location="3p", verbose=False)["count"] == 0
    )
    internal = p26.set_strand_modification(
        added["design"],
        1,
        location="internal",
        display_text="T*",
        vendor_code="/iBiodT/",
        index=3,
        allowed_bases=["T"],
        verbose=False,
    )
    assert _model(internal["design"]).strands[1].modifications_int
    plain = p26.get_strand_without_modifications(internal["design"], 1, verbose=False)
    assert "5prime_modification" not in plain["strand"]
    removed = p26.remove_strand_modification(added["design"], 1, location="5p", verbose=False)
    assert _model(removed["design"]).strands[1].modification_5p is None
    with pytest.raises(ConfigurationError):
        p26.set_strand_modification(
            sequenced, 1, location="internal", display_text="x", vendor_code="y", verbose=False
        )


def test_base_pairs(sequenced: dict) -> None:
    pairs = p26.compute_base_pairs(sequenced, verbose=False)
    assert pairs["total"] == 24
    expected = _model(sequenced).base_pairs()
    assert pairs["base_pairs"] == {str(h): sorted(int(o) for o in v) for h, v in expected.items()}


def test_strand_editing(design: dict) -> None:
    helix = p26.add_design_helix(design, idx=2, max_offset=48, grid_position=[0, 2], verbose=False)
    assert helix["helices"] == 3
    added = p26.add_design_strand(
        helix["design"],
        {"domains": [{"helix": 2, "forward": True, "start": 0, "end": 10}], "name": "n1"},
        verbose=False,
    )
    assert added["added_strand_index"] == 2
    assert p26.find_strand_by_name(added["design"], "n1", verbose=False)["strand_index"] == 2
    assert not p26.find_strand_by_name(added["design"], "zz", verbose=False)["found"]
    removed = p26.remove_design_strand(added["design"], 2, verbose=False)
    assert removed["strands"] == 2
    ends = p26.list_strands_at_helix_end(design, 0, end="5p", verbose=False)
    assert [s["strand_index"] for s in ends["strands"]] == [0]
    at = p26.find_domains_at_offset(design, 0, 10, verbose=False)
    assert len(at["domains"]) == 1
    none = p26.find_domains_at_offset(design, 0, 10, forward=False, verbose=False)
    assert none["domains"] == []


def test_build_strand(design: dict) -> None:
    result = p26.build_strand(
        design,
        helix=0,
        offset=24,
        steps=[
            {"op": "to", "offset": 40},
            {"op": "loopout", "helix": 1, "length": 3},
            {"op": "move", "delta": -16},
            {"op": "with_deletions", "deletions": [30]},
            {"op": "as_circular"},
        ],
        verbose=False,
    )
    strand = _model(result["design"]).strands[result["added_strand_index"]]
    assert strand.circular and len(strand.domains) == 3
    with pytest.raises(ConfigurationError):
        p26.build_strand(design, helix=0, offset=0, steps=[{"op": "shell"}], verbose=False)
    with pytest.raises(ConfigurationError):
        p26.build_strand(
            design, helix=0, offset=0, steps=[{"op": "to", "offset": 4, "x": 1}], verbose=False
        )


def test_domain_list_edits(sequenced: dict) -> None:
    design = p26.build_strand(
        _design(),
        helix=0,
        offset=24,
        steps=[
            {"op": "to", "offset": 32},
            {"op": "cross", "helix": 1},
            {"op": "move", "delta": -8},
        ],
        verbose=False,
    )
    idx = design["added_strand_index"]
    assert design["strand"]["domains"] == 2
    removed = p26.remove_strand_domain(design["design"], idx, 1, verbose=False)
    assert removed["assigned"]["domains"] == 1
    appended = p26.insert_strand_domain(
        removed["design"],
        idx,
        {"helix": 1, "forward": False, "start": 24, "end": 32},
        verbose=False,
    )
    assert appended["assigned"]["domains"] == 2
    rotated = p26.rotate_strand_domains(appended["design"], idx, 1, verbose=False)
    first = _model(rotated["design"]).strands[idx].domains[0]
    assert first.helix == 1
    replaced = p26.set_strand_domains(
        appended["design"],
        idx,
        [{"helix": 0, "forward": True, "start": 24, "end": 30}],
        verbose=False,
    )
    assert replaced["assigned"]["domains"] == 1
    reversed_one = p26.reverse_strands(appended["design"], strand_index=idx, verbose=False)
    assert _model(reversed_one["design"]).strands[idx].domains[0].helix == 1
    reversed_all = p26.reverse_strands(appended["design"], verbose=False)
    assert reversed_all["strands"] == appended["strands"]


def test_strand_properties(design: dict) -> None:
    result = p26.set_strand_properties(
        design, 1, scaffold=True, circular=False, vendor_fields={"scale": "100nm"}, verbose=False
    )
    strand = _model(result["design"]).strands[1]
    assert strand.is_scaffold and strand.vendor_fields.scale == "100nm"
    with pytest.raises(ConfigurationError):
        p26.set_strand_properties(design, 1, vendor_fields={"bad": "x"}, verbose=False)


def test_overlap_and_indices(sequenced: dict) -> None:
    same = p26.check_domain_overlap(
        sequenced, 0, 1, domain_index=0, other_domain_index=0, verbose=False
    )
    assert same["overlaps"] is False
    whole = p26.check_domain_overlap(sequenced, 0, 1, verbose=False)
    assert whole["overlaps"] is True
    first = p26.check_domain_overlap(
        sequenced, 0, 0, domain_index=0, other_domain_index=0, verbose=False
    )
    assert first["overlap_interval"] == [0, 24] and first["overlaps_illegally"]
    assert not first["overlaps"]
    found = p26.find_overlapping_domains(sequenced, 0, 0, verbose=False)
    model = _model(sequenced)
    others = [d for st in model.strands for d in st.domains if d is not model.strands[0].domains[0]]
    expected = list(model.strands[0].domains[0].find_overlapping_ranges(others))
    assert len(found["domains"]) == len(expected) > 0
    pairs = p26.find_helix_domain_overlap_pairs(sequenced, verbose=False)
    assert pairs["count"] == 1 and pairs["pairs"][0]["helix"] == 1
    crossing = p26.domain_has_crossover(sequenced, 0, 0, 23, verbose=False)
    assert crossing["has_crossover"] is True
    mapped = p26.map_domain_offset_to_dna_index(sequenced, 0, 1, offset=10, verbose=False)
    model = _model(sequenced)
    assert mapped["start_index"] == model.strands[0].dna_index_start_domain(
        model.strands[0].domains[1]
    )
    assert mapped["dna_index"] == model.strands[0].domains[1].domain_offset_to_strand_dna_idx(
        10, True
    )


def test_crossover_listings() -> None:
    rectangle = p26.create_origami_rectangle(num_helices=4, num_cols=10, verbose=False)["design"]
    model = _model(rectangle)
    crossovers = p26.list_helix_crossovers(rectangle, 0, verbose=False)
    assert crossovers["count"] == len(model.helices[0].crossover_addresses(model.helices))
    split = p26.split_helix_strands(rectangle, 0, verbose=False)
    staples, scaffolds = model.helices[0].get_strand_sets()
    assert len(split["staple_domains"]) == len(staples)
    assert len(split["scaffold_domains"]) == len(scaffolds)
    assert p26.list_crossover_points(rectangle, verbose=False)["scaffold_low"][0] == [4, 26, 15]
    assert p26.get_helix_neighbors(rectangle, 0, verbose=False)["neighbors"] == [
        None,
        1,
        None,
        None,
    ]
    potential = p26.list_potential_crossovers(rectangle, 0, verbose=False)
    assert potential["count"] == len(model.potential_crossover_list(model.helices[0]))


def test_crossover_nick_ligate_deletion_insertion() -> None:
    sc = pytest.importorskip("scadnano")
    base = sc.Design(helices=[sc.Helix(max_offset=48) for _ in range(2)], grid=sc.Grid.square)
    base.draw_strand(0, 0).to(24)
    base.draw_strand(1, 24).move(-24)
    data = json.loads(base.to_json())
    crossed = p26.add_crossover(data, helix=0, helix2=1, offset=23, forward=True, verbose=False)
    assert crossed["strands"] == 1
    nicked = p26.add_design_nick(crossed["design"], 0, 5, True, verbose=False)
    assert nicked["strands"] == 2
    ligated = p26.ligate_design_nick(nicked["design"], 0, 5, True, verbose=False)
    assert ligated["strands"] == 1
    deleted = p26.add_design_deletion(crossed["design"], 0, 10, verbose=False)
    assert deleted["strand_lengths"][0] == crossed["strand_lengths"][0] - 1
    inserted = p26.add_design_insertion(crossed["design"], 0, 10, 2, verbose=False)
    assert inserted["strand_lengths"][0] == crossed["strand_lengths"][0] + 2
    inlined = p26.inline_deletions_insertions(inserted["design"], verbose=False)
    assert inlined["strand_lengths"] == inserted["strand_lengths"]
    assert "insertions" not in json.dumps(inlined["design"])
    full = p26.add_crossover(
        data, helix=0, helix2=1, offset=10, forward=True, full=True, verbose=False
    )
    assert full["strands"] >= 1


def test_move_strands(design: dict) -> None:
    moved = p26.move_strands(design, 3, verbose=False)
    assert _model(moved["design"]).strands[0].domains[0].start == 3
    helix_move = p26.move_strands(design, 0, axis="helix", verbose=False)
    assert helix_move["strands"] == 2
    with pytest.raises(ConfigurationError):
        p26.move_strands(design, 1, axis="z", verbose=False)


def test_origami_rectangle_and_autostaple() -> None:
    pytest.importorskip("scadnano")
    rectangle = p26.create_origami_rectangle(num_helices=4, num_cols=10, verbose=False)
    assert rectangle["helices"] == 4 and rectangle["scaffolds"] == 1
    original = json.loads(
        __import__("scadnano.origami_rectangle", fromlist=["create"])
        .create(num_helices=4, num_cols=10)
        .to_json()
    )
    assert rectangle["design"] == original
    assert rectangle["strand_lengths"][0] == 640
    twist = p26.create_origami_rectangle(
        num_helices=4, num_cols=10, twist_correction_deletion_spacing=3, verbose=False
    )
    assert sum(twist["strand_lengths"]) < sum(rectangle["strand_lengths"])
    stapled = p26.autostaple_design(rectangle["design"], verbose=False)
    assert stapled["scaffolds"] == 1 and stapled["strands"] < rectangle["strands"]
    bare = p26.create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False, verbose=False)
    added = p26.add_twist_correction_deletions(
        bare["design"],
        offset_start=0,
        deletion_spacing=2,
        deletion_start_col=1,
        deletion_offset=-1,
        num_helices=4,
        num_cols=10,
        verbose=False,
    )
    assert sum(added["strand_lengths"]) < sum(bare["strand_lengths"])
    one = p26.add_deletion_in_range(bare["design"], helix=1, start=17, end=30, verbose=False)
    assert sum(one["strand_lengths"]) == sum(bare["strand_lengths"]) - 2
    with pytest.raises(ConfigurationError):
        p26.create_origami_rectangle(num_helices=4, num_cols=10, nick_pattern="x", verbose=False)


def test_helix_roll_and_backbone() -> None:
    rectangle = p26.create_origami_rectangle(num_helices=4, num_cols=10, verbose=False)["design"]
    model = _model(rectangle)
    plan = p26.relax_helix_rolls(rectangle, compute_only=True, verbose=False)
    assert plan["applied"] is False and plan["design"] == rectangle
    expected = model.helices[0].compute_relaxed_roll_delta(
        model.helices, model.grid, model.geometry
    )
    assert plan["roll_deltas"]["0"] == pytest.approx(expected)
    applied = p26.relax_helix_rolls(rectangle, verbose=False)
    model.relax_helix_rolls()
    assert applied["rolls"]["0"] == pytest.approx(model.helices[0].roll)
    angle = p26.get_backbone_angle(rectangle, 0, 10, True, verbose=False)["angle_degrees"]
    assert angle == pytest.approx(
        _model(rectangle).helices[0].backbone_angle_at_offset(10, True, model.geometry)
    )
