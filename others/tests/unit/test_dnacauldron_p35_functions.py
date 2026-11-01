"""Tests for dnakit.molecular_lab.dnacauldron_p35_functions (records, sticky ends, mixes, graphs, plans, BASIC).

Offline: every input is a small synthetic DNA fixture generated from a fixed seed, and the outputs go to
tmp_path. Results are cross-checked against DNA Cauldron's own API called directly.
"""

from __future__ import annotations

import inspect
import random
import socket
from pathlib import Path

import pytest

pytest.importorskip("dnacauldron")
pytest.importorskip("Bio")

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import dnacauldron_p35_functions as m

PNG = b"\x89PNG\r\n\x1a\n"


class _NoConnectSocket(socket.socket):
    def connect(self, *args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    connect_ex = connect  # type: ignore[assignment]


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", _NoConnectSocket)


@pytest.fixture(autouse=True)
def _close_figures() -> None:
    plt = pytest.importorskip("matplotlib.pyplot")
    plt.close("all")


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT", "TGCA"))


class Dna:
    def __init__(self, seed: int) -> None:
        self.rng = random.Random(seed)

    def random(self, length: int) -> str:
        while True:
            text = "".join(self.rng.choice("ACGT") for _ in range(length))
            if not any(site in text for site in ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC")):
                return text


def golden_gate_parts() -> dict[str, str]:
    dna = Dna(11)

    def part(left: str, right: str) -> str:
        return dna.random(60) + "CGTCTCA" + left + dna.random(30) + right + "AGAGACG" + dna.random(60)

    return {"p1": part("AATG", "CCTA"), "p2": part("CCTA", "GGAT"), "p3": part("GGAT", "AATG")}


PARTS = golden_gate_parts()


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(m.__all__) == 28
    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


# Records.


def test_build_and_annotate_record() -> None:
    built = m.build_dnacauldron_record(
        PARTS["p1"], record_id="x", topology="circular", features=[{"start": 0, "end": 10, "label": "a"}], verbose=False
    )["record"]
    assert built["topology"] == "circular" and built["length"] == len(PARTS["p1"])
    assert built["features"] == [{"type": "misc_feature", "start": 0, "end": 10, "strand": 1, "label": "a"}]
    edited = m.annotate_dnacauldron_record(
        built, features=[{"type": "gene", "start": 5, "end": 20, "strand": -1, "label": "g"}], topology="linear", verbose=False
    )["record"]
    assert edited["topology"] == "linear" and len(edited["features"]) == 2
    assert edited["features"][1]["strand"] == -1 and built["topology"] == "circular"  # input unchanged


def test_record_validation_errors() -> None:
    with pytest.raises(ConfigurationError, match="only A, C, G and T"):
        m.build_dnacauldron_record("ACGN", verbose=False)
    with pytest.raises(ConfigurationError, match="start < end"):
        m.build_dnacauldron_record("ACGT", features=[{"start": 3, "end": 2}], verbose=False)
    with pytest.raises(ConfigurationError, match="topology"):
        m.build_dnacauldron_record("ACGT", topology="ring", verbose=False)


def test_write_records_and_read_back_from_repository(tmp_path: Path) -> None:
    record = m.build_dnacauldron_record(PARTS["p1"], record_id="part_one", topology="circular", verbose=False)["record"]
    out = tmp_path / "gb"
    written = m.write_dnacauldron_records([record], str(out), verbose=False)
    path = Path(written["files"][0]["path"])
    assert path.name == "part_one.gb" and path.read_text().startswith("LOCUS")
    fasta = m.write_dnacauldron_records([record], str(tmp_path / "fa"), output_format="fasta", verbose=False)
    assert Path(fasta["files"][0]["path"]).read_text().splitlines()[1] == PARTS["p1"][:60]
    with pytest.raises(ConfigurationError, match="already exists"):
        m.write_dnacauldron_records([record], str(out), verbose=False)
    found = m.get_dnacauldron_repository_records(["part_one"], files=[str(path)], verbose=False)
    assert found["records"][0]["sequence"] == PARTS["p1"]
    assert found["records"][0]["topology"] == "circular"
    several = m.get_dnacauldron_repository_records(["part_one"], folder=str(out), verbose=False)
    assert several["repository_ids"] == ["part_one"]
    with pytest.raises(ConfigurationError, match="not found"):
        m.get_dnacauldron_repository_records(["missing"], files=[str(path)], verbose=False)


# Sticky ends.


def test_digestion_matches_dnacauldron_directly() -> None:
    import dnacauldron as dc
    from Bio import Restriction
    from Bio.Seq import Seq

    result = m.digest_dnacauldron_sticky_ends(PARTS["p1"], ["BsmBI"], circular=True, verbose=False)
    direct = dc.StickyEndSeq.list_from_sequence_digestion(Seq(PARTS["p1"]), Restriction.BsmBI, linear=False)
    assert [f["sequence"] for f in result["fragments"]] == [str(f) for f in direct]
    assert result["fragment_count"] == 2
    insert = next(f for f in result["fragments"] if f["length"] == 30)
    assert (insert["left_end"], insert["right_end"]) == (
        {"sequence": "AATG", "strand": 1},
        {"sequence": "CCTA", "strand": -1},
    )
    by_record = m.digest_dnacauldron_sticky_ends(PARTS["p1"], ["BsmBI"], circular=True, method="record", verbose=False)
    assert by_record["fragment_count"] == 2
    none = m.digest_dnacauldron_sticky_ends("ACGT" * 10, ["BsmBI"], verbose=False)
    assert none["fragments"][0]["left_end"] is None
    with pytest.raises(ConfigurationError, match="unknown restriction enzyme"):
        m.digest_dnacauldron_sticky_ends(PARTS["p1"], ["NoSuchEnzyme"], verbose=False)


def sticky() -> dict[str, object]:
    return {
        "sequence": "GGAGCT",
        "left_end": {"sequence": "AATG", "strand": 1},
        "right_end": {"sequence": "CCTA", "strand": -1},
    }


def test_sticky_sequence_transforms() -> None:
    reverse = m.transform_dnacauldron_sticky_sequence(sticky(), verbose=False)["fragment"]
    assert reverse["sequence"] == revcomp("GGAGCT")
    assert reverse["left_end"] == {"sequence": revcomp("CCTA"), "strand": 1}
    assert reverse["right_end"] == {"sequence": revcomp("AATG"), "strand": -1}
    again = m.transform_dnacauldron_sticky_sequence(reverse, verbose=False)["fragment"]
    assert again["sequence"] == "GGAGCT" and again["left_end"] == sticky()["left_end"]
    flat = m.transform_dnacauldron_sticky_sequence(sticky(), operation="to_standard", verbose=False)
    assert flat["sequence"] == "AATGGGAGCTCCTA"
    middle = m.transform_dnacauldron_sticky_sequence(sticky(), operation="to_standard", discard_sticky_ends=True, verbose=False)
    assert middle["sequence"] == "GGAGCT"
    sliced = m.transform_dnacauldron_sticky_sequence(sticky(), operation="slice", start=1, end=4, verbose=False)["fragment"]
    assert sliced == {"sequence": "GAG", "length": 3, "left_end": None, "right_end": None}
    with pytest.raises(ConfigurationError, match="strand"):
        m.transform_dnacauldron_sticky_sequence({"sequence": "AC", "left_end": {"sequence": "A", "strand": 2}}, verbose=False)


def test_sticky_end_compatibility() -> None:
    nxt = {"sequence": "TTTT", "left_end": {"sequence": "CCTA", "strand": 1}, "right_end": {"sequence": "AATG", "strand": -1}}
    checked = m.check_dnacauldron_sticky_end_compatibility([sticky(), nxt], circular=True, verbose=False)
    assert [row["compatible"] for row in checked["pairs"]] == [True, True]
    assert checked["all_compatible"] is True
    bad = m.check_dnacauldron_sticky_end_compatibility([sticky(), sticky()], verbose=False)
    assert bad["all_compatible"] is False and bad["pairs"][0]["right_end"] == "CCTA"
    with pytest.raises(ConfigurationError, match="2 to 200"):
        m.check_dnacauldron_sticky_end_compatibility([sticky()], verbose=False)


def test_sticky_fragment_to_record() -> None:
    record = m.sticky_fragment_to_dnacauldron_record(sticky(), verbose=False)["record"]
    assert record["sequence"] == "AATGGGAGCTCCTA"
    labels = [f["label"] for f in record["features"]]
    assert labels == ["(+) strand", "(-) strand"]
    end = m.sticky_fragment_to_dnacauldron_record(sticky(), part="right_end", verbose=False)["record"]
    assert end["sequence"] == "CCTA" and end["features"][0]["label"] == "(-) strand"
    with pytest.raises(ConfigurationError, match="no left_end"):
        m.sticky_fragment_to_dnacauldron_record({"sequence": "ACGT"}, part="left_end", verbose=False)


def test_homology_check() -> None:
    same = "ACGTACGTACGTACGTAC"
    result = m.check_dnacauldron_homology(same, same, verbose=False)
    assert result["acceptable"] is True and result["tm"] == 54 and result["edit_distance"] == 0
    assert m.check_dnacauldron_homology(same, same[:-1] + "A", verbose=False)["acceptable"] is False
    assert m.check_dnacauldron_homology(same, same[:-1] + "A", max_distance=1, verbose=False)["acceptable"] is True
    assert m.check_dnacauldron_homology("ACGT", "ACGT", verbose=False)["acceptable"] is False  # shorter than 15
    assert m.check_dnacauldron_homology(same, same, min_tm=60, verbose=False)["acceptable"] is False
    assert m.check_dnacauldron_homology(same, verbose=False)["edit_distance"] is None


# Mixes.


def test_slots_and_slots_graph() -> None:
    slots = m.compute_dnacauldron_slots(PARTS, verbose=False)
    assert slots["slots"] == [
        {"left_overhang": "AATG", "right_overhang": "CCTA", "parts": ["p1"]},
        {"left_overhang": "ATCC", "right_overhang": "TAGG", "parts": ["p2"]},
        {"left_overhang": "CATT", "right_overhang": "ATCC", "parts": ["p3"]},
    ]
    graph = m.compute_dnacauldron_slots_graph(PARTS, verbose=False)
    assert graph["directed"] is True and graph["node_count"] == 6 and graph["edge_count"] == 6
    assert {n["id"] for n in graph["nodes"] if n["kind"] == "slot"} == {"AATG-CCTA", "ATCC-TAGG", "CATT-ATCC"}
    flat = m.compute_dnacauldron_slots_graph(PARTS, with_overhangs=False, verbose=False)
    assert flat["directed"] is False and flat["node_count"] == 3
    with pytest.raises(ConfigurationError, match="unknown restriction enzyme"):
        m.compute_dnacauldron_slots(PARTS, enzymes=["Bogus"], verbose=False)


def test_connections_graph_matches_dnacauldron() -> None:
    import dnacauldron as dc
    from dnacauldron.biotools import sequence_to_biopython_record, set_record_topology

    records = []
    for key, text in sorted(PARTS.items()):
        record = sequence_to_biopython_record(text, id=key)
        set_record_topology(record, "circular")
        records.append(record)
    mix = dc.generate_type2s_restriction_mix(records, enzyme="BsmBI")
    expected = sorted([int(a), int(b)] for a, b in mix.filtered_connections_graph.edges())
    result = m.compute_dnacauldron_connections_graph(PARTS, verbose=False)
    assert result["edges"] == expected and result["edge_count"] == 6
    assert result["fragment_count"] == len(mix.fragments_dict)
    first = result["fragments"][0]
    assert set(first) == {"index", "label", "part", "reversed", "length", "left_end", "right_end"}


def test_circular_assemblies_all_output_modes() -> None:
    records = m.enumerate_dnacauldron_circular_assemblies(PARTS, verbose=False)
    assert records["returned"] == 1
    construct = records["assemblies"][0]
    assert construct["parts"] == ["p1", "p2", "p3"] and construct["topology"] == "circular"
    assert construct["length"] == 102 and len(construct["sequence"]) == 102
    sets = m.enumerate_dnacauldron_circular_assemblies(PARTS, output="fragments_sets", verbose=False)
    assert sets["assemblies"] == [["p1", "p2", "p3"]]
    random_sets = m.enumerate_dnacauldron_circular_assemblies(
        PARTS, output="fragments_sets", sampling="random", random_seed=3, verbose=False
    )
    assert random_sets["assemblies"] == [["p1", "p2", "p3"]]
    random_records = m.enumerate_dnacauldron_circular_assemblies(PARTS, sampling="random", random_seed=3, verbose=False)
    assert random_records["assemblies"][0]["sequence"] == construct["sequence"]
    annotated = m.enumerate_dnacauldron_circular_assemblies(PARTS, annotate_homologies=True, verbose=False)
    assert annotated["returned"] == 1
    with pytest.raises(ConfigurationError, match="output"):
        m.enumerate_dnacauldron_circular_assemblies(PARTS, output="x", verbose=False)


def test_random_sampling_state_is_restored_and_seed_repeats() -> None:
    import numpy

    numpy.random.seed(5)
    before = numpy.random.get_state()[1][:3].tolist()
    first = m.enumerate_dnacauldron_circular_assemblies(PARTS, sampling="random", random_seed=1, verbose=False)
    second = m.enumerate_dnacauldron_circular_assemblies(PARTS, sampling="random", random_seed=1, verbose=False)
    assert first["assemblies"] == second["assemblies"]
    assert numpy.random.get_state()[1][:3].tolist() == before


def test_linear_assemblies_and_fragment_chains() -> None:
    linear = m.enumerate_dnacauldron_linear_assemblies(PARTS, verbose=False)
    assert linear["returned"] == 6
    assert all(row["left_end"] and row["right_end"] for row in linear["assemblies"])
    three = [row for row in linear["assemblies"] if len(row["parts"]) == 3]
    assert three and three[0]["length"] == len(three[0]["sequence"])
    assert {row["parts"][0] for row in three} <= {"p1", "p2", "p3", "p1_r", "p2_r", "p3_r"}
    longer = m.enumerate_dnacauldron_linear_assemblies(PARTS, min_parts=3, verbose=False)
    assert all(len(row["parts"]) >= 3 for row in longer["assemblies"])
    capped = m.enumerate_dnacauldron_linear_assemblies(PARTS, max_results=2, verbose=False)
    assert capped["returned"] == 2 and capped["limit_reached"] is True
    chains = m.transform_dnacauldron_fragment_chain(PARTS, verbose=False)
    assert chains["chain_count"] == 2
    for chain in chains["chains"]:
        assert chain["after"] == ["p1", "p2", "p3"] and chain["is_standardized"] is True
    reverse = m.transform_dnacauldron_fragment_chain(PARTS, operation="reverse_complement", verbose=False)
    for chain in reverse["chains"]:
        flipped = [name[:-2] if name.endswith("_r") else name + "_r" for name in chain["before"]][::-1]
        assert chain["after"] == flipped


def test_homology_and_ligase_cycling_mixes() -> None:
    dna = Dna(2)
    x = dna.random(80)
    y = x[-30:] + dna.random(50)
    z = y[-30:] + dna.random(50) + x[:30]
    gibson = m.enumerate_dnacauldron_circular_assemblies({"X": x, "Y": y, "Z": z}, mix_type="homology", verbose=False)
    assert gibson["returned"] == 1 and gibson["assemblies"][0]["length"] == 80 + 50 + 50
    strict = m.enumerate_dnacauldron_circular_assemblies(
        {"X": x, "Y": y, "Z": z}, mix_type="homology", min_size=40, verbose=False
    )
    assert strict["returned"] == 0
    a, b, c = dna.random(60), dna.random(60), dna.random(60)
    oligos = {"oAB": a[-20:] + b[:20], "oBC": b[-20:] + c[:20], "oCA": c[-20:] + a[:20]}
    circular = m.simulate_dnacauldron_lcr_assembly({"A": a, "B": b, "C": c}, oligos, verbose=False)
    assert circular["returned"] == 1
    assert sorted(circular["constructs"][0]["bridging_oligos"]) == ["oAB", "oBC", "oCA"]
    assert circular["constructs"][0]["length"] == 180 and sorted(circular["constructs"][0]["parts"]) == ["A", "B", "C"]
    linear = m.simulate_dnacauldron_lcr_assembly(
        {"A": a, "B": b, "C": c}, {"oAB": oligos["oAB"], "oBC": oligos["oBC"]}, circular=False, verbose=False
    )
    assert sorted(row["parts"] for row in linear["constructs"]) == [["A", "B"], ["A", "B", "C"], ["B", "C"]]
    found = m.find_dnacauldron_bridging_oligo(a, b, oligos, verbose=False)
    assert (found["found"], found["oligo"], found["start"], found["end"]) == (True, "oAB", 40, 80)
    assert m.find_dnacauldron_bridging_oligo(a, c, {"oAB": oligos["oAB"]}, verbose=False)["found"] is False
    assert m.find_dnacauldron_bridging_oligo(a, b, oligos, min_size=25, verbose=False)["found"] is False
    with pytest.raises(ConfigurationError, match="bridging_oligos"):
        m.compute_dnacauldron_connections_graph({"A": a}, mix_type="ligase_cycling", verbose=False)


def test_select_connectors() -> None:
    result = m.select_dnacauldron_connectors({"p1": PARTS["p1"], "p3": PARTS["p3"]}, {"c2": PARTS["p2"]}, verbose=False)
    assert result["selected_connectors"] == ["c2"]
    assert result["constructs"][0]["parts"] == ["p1", "c2", "p3"]
    with pytest.raises(ConfigurationError, match="connector"):
        m.select_dnacauldron_connectors({"p1": PARTS["p1"], "p3": PARTS["p3"]}, {"c1": PARTS["p1"]}, verbose=False)


# Plots.


def test_plots_write_png_and_return_data(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt

    graph = m.plot_dnacauldron_connections_graph(PARTS, str(tmp_path / "graph.png"), verbose=False)
    assert Path(graph["path"]).read_bytes()[:8] == PNG
    assert graph["data"]["node_count"] == 3 and graph["data"]["edge_count"] == 3
    slots = m.plot_dnacauldron_slots_graph(PARTS, str(tmp_path / "slots.png"), highlighted_parts=["p1"], verbose=False)
    assert Path(slots["path"]).read_bytes()[:8] == PNG
    assert slots["data"]["highlighted_parts"] == ["p1"] and len(slots["data"]["slots"]) == 3
    record = {"id": "x", "sequence": PARTS["p1"], "topology": "circular"}
    cuts = m.plot_dnacauldron_cuts(record, ["BsmBI"], str(tmp_path / "cuts.png"), verbose=False)
    assert Path(cuts["path"]).read_bytes()[:8] == PNG
    assert cuts["data"]["cuts"] == {"BsmBI": [68, 102]} and cuts["data"]["linear"] is False
    linear = m.plot_dnacauldron_cuts({**record, "topology": "linear"}, ["BsmBI", "EcoRI"], str(tmp_path / "lin.png"), verbose=False)
    assert linear["data"]["cuts"]["EcoRI"] == [] and linear["data"]["linear"] is True
    assert plt.get_fignums() == []
    with pytest.raises(ConfigurationError, match="already exists"):
        m.plot_dnacauldron_cuts(record, ["BsmBI"], str(tmp_path / "cuts.png"), verbose=False)
    with pytest.raises(ConfigurationError, match="must end with"):
        m.plot_dnacauldron_cuts(record, ["BsmBI"], str(tmp_path / "cuts.pdf"), verbose=False)


# Overhangs of a record.


def test_record_overhangs_and_insert_overhangs(tmp_path: Path) -> None:
    features = [
        {"type": "homology", "start": 10, "end": 14, "label": "CCTA"},
        {"type": "homology", "start": 2, "end": 6, "label": "AATG"},
        {"type": "homology", "start": 20, "end": 28, "label": "homology"},
    ]
    record = {"id": "c", "sequence": "A" * 40, "topology": "circular", "features": features}
    assert m.list_dnacauldron_record_overhangs(record, verbose=False)["overhangs"] == ["AATG", "CCTA"]
    located = m.list_dnacauldron_record_overhangs(record, with_locations=True, verbose=False)
    assert located["overhangs"] == [[2, "AATG"], [10, "CCTA"]]
    gb = m.write_dnacauldron_records([record], str(tmp_path), verbose=False)["files"][0]["path"]
    assert m.list_dnacauldron_record_overhangs(genbank_path=gb, verbose=False)["overhangs"] == ["AATG", "CCTA"]
    with pytest.raises(ConfigurationError, match="exactly one"):
        m.list_dnacauldron_record_overhangs(verbose=False)
    pair = m.get_dnacauldron_insert_overhangs(PARTS["p1"], verbose=False)
    assert (pair["left_overhang"], pair["right_overhang"]) == ("AATG", "CCTA")
    assert m.get_dnacauldron_insert_overhangs(revcomp(PARTS["p1"]), verbose=False)["left_overhang"] == "AATG"
    raw = m.get_dnacauldron_insert_overhangs(revcomp(PARTS["p1"]), standardize=False, verbose=False)
    assert (raw["left_overhang"], raw["right_overhang"]) == ("TAGG", "CATT")


# Oligo annealing and BASIC.


def test_oligo_pair_annealing() -> None:
    long = "AAAACCCCGGGGTTTTACGTACGT"
    short = revcomp("AAAACCCCGGGG")
    result = m.anneal_dnacauldron_oligo_pair(long, short, verbose=False)
    assert result["product_count"] == 1
    product = result["products"][0]
    assert product["sequence"] == "AAAACCCCGGGG"
    assert product["right_end"] == {"sequence": "TTTTACGTACGT", "strand": 1}
    assert result["products"] == m.anneal_dnacauldron_oligo_pair(short, long, verbose=False)["products"]
    with pytest.raises(ConfigurationError, match="No homology"):
        m.anneal_dnacauldron_oligo_pair("AAAAAAAAAA", "CCCCCCCCCC", verbose=False)


def basic_inputs() -> tuple[dict[str, str], dict[str, list[str]], list[str]]:
    dna = Dna(3)
    ends = [("CTCC", "GCTT"), ("TGAC", "AGGA")]
    parts = {
        f"part{i + 1}": dna.random(10) + "GGTCTCA" + left + dna.random(40) + right + "AGAGACC" + dna.random(10)
        for i, (left, right) in enumerate(ends)
    }
    linker_a, linker_b = dna.random(12), dna.random(12)
    cores = [dna.random(20) for _ in range(4)]

    def adapter(left: str, core: str, right: str) -> list[str]:
        return [left + core + right, revcomp(core)]

    adapters = {
        "ad1": adapter(revcomp("CTCC"), cores[0], linker_a),
        "ad2": adapter("GCTT", cores[1], linker_b),
        "ad3": adapter(revcomp("TGAC"), cores[2], revcomp(linker_b)),
        "ad4": adapter("AGGA", cores[3], revcomp(linker_a)),
    }
    return parts, adapters, ["ad1", "part1", "ad2", "ad3", "part2", "ad4"]


def test_basic_assembly_simulation() -> None:
    parts, adapters, order = basic_inputs()
    result = m.simulate_dnacauldron_basic_assembly(parts, adapters, order, verbose=False)
    assert result["construct_count"] == 1 and result["errors"] == []
    construct = result["constructs"][0]
    assert sorted(construct["parts"]) == ["part1", "part2"] and construct["topology"] == "circular"
    assert result["outcomes"] == {"cancelled_assemblies": 0, "errored_assemblies": 0, "valid_assemblies": 5}
    assert set(construct["sequence"]) <= set("ACGT") and construct["length"] > 100
    with pytest.raises(ConfigurationError, match="multiple of 3"):
        m.simulate_dnacauldron_basic_assembly(parts, adapters, order[:-1], verbose=False)
    with pytest.raises(ConfigurationError, match="neither parts nor adapters"):
        m.simulate_dnacauldron_basic_assembly(parts, adapters, ["ad1", "nope", "ad2"], verbose=False)


# Assembly plans.


def test_plan_read_write_and_tabulate(tmp_path: Path) -> None:
    plan = tmp_path / "plan.csv"
    plan.write_text("construct1,class: type2s_assembly,p1,p2,p3\nbig,class: type2s_assembly,construct1,p4,\n")
    read = m.read_dnacauldron_assembly_plan(str(plan), verbose=False)
    assert read["plan_name"] == "plan" and read["assembly_count"] == 2
    by_name = {row["name"]: row for row in read["assemblies"]}
    assert by_name["construct1"]["parts"] == ["p1", "p2", "p3"] and by_name["construct1"]["level"] == 1
    assert by_name["big"]["depends_on"] == ["construct1"] and by_name["big"]["level"] == 2
    assert by_name["big"]["class"] == "Type2sRestrictionAssembly"
    out = tmp_path / "out.csv"
    written = m.write_dnacauldron_assembly_plan(
        [{"name": "c2", "parts": ["c1", "d"]}, {"name": "c1", "parts": ["a", "b"]}], output_path=str(out), verbose=False
    )
    assert written["rows"] == [{"assembly": "c1", "parts": "a, b"}, {"assembly": "c2", "parts": "c1, d"}]
    assert out.read_text() == "construct,parts\nc2,c1,d\nc1,a,b"
    with pytest.raises(ConfigurationError, match="Multiple assemblies named"):
        m.write_dnacauldron_assembly_plan([{"name": "c", "parts": ["a"]}, {"name": "c", "parts": ["b"]}], verbose=False)
    sequences = tmp_path / "seq"
    sequences.mkdir()
    for key, text in PARTS.items():
        (sequences / f"{key}.fa").write_text(f">{key}\n{text}\n")
    gg = tmp_path / "gg.csv"
    gg.write_text("asm1,class: type2s_assembly,p1,p2,p3\n")
    table = m.tabulate_dnacauldron_plan_constructs(str(gg), str(sequences), verbose=False)
    assert table["construct_count"] == 1 and table["outcomes"]["valid_assemblies"] == 1
    row = table["constructs"][0]
    assert row["construct_id"] == "asm1" and row["number_of_parts"] == 3 and row["construct_size"] == 102
    assert row["parts"] == ["p1", "p2", "p3"] and row["assembly_level"] == 1
    with pytest.raises(ConfigurationError, match="existing file"):
        m.read_dnacauldron_assembly_plan(str(tmp_path / "missing.csv"), verbose=False)
