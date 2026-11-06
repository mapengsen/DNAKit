"""Real implementations of the microbiome, design_crispr, epigenome, regulation and genome3d candidate functions.

Every expected number below was derived by hand (or from a textbook definition) and, where scipy is available,
cross-checked against it.
"""

from __future__ import annotations

import math
import re

import pytest

from dnakit import design_crispr as crispr
from dnakit import epigenome as epi
from dnakit import genome3d as g3d
from dnakit import microbiome as micro
from dnakit import regulation as reg
from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError

np = pytest.importorskip("numpy")
scipy_stats = pytest.importorskip("scipy.stats")


def ok(envelope: dict) -> object:
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# ============================================================================ microbiome
TABLE = {
    "s1": {"a": 5, "b": 5, "c": 5, "d": 5},
    "s2": {"a": 1, "b": 1, "c": 2, "d": 6},
    "s3": {"a": 10, "b": 0, "c": 3, "d": 0},
}


def test_alpha_diversity_hand_values_even_community() -> None:
    # four equally abundant taxa: H = ln 4, Gini-Simpson = 1 - 4*(1/16) = 0.75, 1/D = 4, J = 1, no singletons
    m = ok(micro.quantify(TABLE, options={"entity": "alpha_diversity"}))["metrics"]["s1"]
    assert m["observed"] == 4
    assert m["shannon"] == pytest.approx(math.log(4))
    assert m["simpson"] == pytest.approx(0.75)
    assert m["inverse_simpson"] == pytest.approx(4.0)
    assert m["pielou_evenness"] == pytest.approx(1.0)
    assert m["chao1"] == pytest.approx(4.0)
    assert m["goods_coverage"] == pytest.approx(1.0)


def test_alpha_diversity_hand_values_uneven_community() -> None:
    # counts 1,1,2,6 (N=10): p = .1,.1,.2,.6; sum p^2 = .42; F1 = 2, F2 = 1
    m = ok(micro.quantify(TABLE, options={"entity": "alpha_diversity"}))["metrics"]["s2"]
    h = -(2 * 0.1 * math.log(0.1) + 0.2 * math.log(0.2) + 0.6 * math.log(0.6))
    assert m["shannon"] == pytest.approx(h)
    assert m["simpson"] == pytest.approx(0.58)
    assert m["inverse_simpson"] == pytest.approx(1 / 0.42)
    assert m["pielou_evenness"] == pytest.approx(h / math.log(4))
    assert m["chao1"] == pytest.approx(4 + 2 * 1 / (2 * (1 + 1)))  # S + F1(F1-1)/(2(F2+1)) = 4.5
    assert m["goods_coverage"] == pytest.approx(0.8)


def test_alpha_diversity_single_taxon_evenness_is_undefined_and_log_base() -> None:
    result = ok(micro.quantify([[7, 0, 0]], options={"entity": "alpha_diversity", "log_base": 2}))[
        "metrics"
    ]["sample1"]
    assert result["pielou_evenness"] is None and result["shannon"] == 0.0
    even = ok(micro.quantify([[1, 1, 1, 1]], options={"entity": "alpha_diversity", "log_base": 2}))[
        "metrics"
    ]["sample1"]
    assert even["shannon"] == pytest.approx(2.0)  # log2(4) bits


def test_relative_abundance_and_clr() -> None:
    rel = ok(micro.quantify([[1, 3]]))["relative_abundance"]
    assert rel == [[0.25, 0.75]]
    clr = ok(micro.quantify([[1, 3]], options={"entity": "clr", "pseudocount": 1.0}))["clr"][0]
    # counts+1 = (2, 4): clr = ln(2) - mean(ln 2, ln 4) = -ln(2)/2 and +ln(2)/2
    assert clr == pytest.approx([-math.log(2) / 2, math.log(2) / 2])


def test_beta_diversity_hand_values_and_scipy() -> None:
    table = {"x": {"a": 10, "b": 0, "c": 5}, "y": {"a": 0, "b": 10, "c": 5}}
    bc = ok(micro.compare(table))["distance_matrix"]
    assert bc[0][1] == pytest.approx(20 / 30)  # sum|a-b| / sum(a+b)
    assert bc[0][0] == 0 and bc[0][1] == bc[1][0]
    jac = ok(micro.compare(table, options={"metric": "jaccard"}))["distance_matrix"]
    assert jac[0][1] == pytest.approx(1 - 1 / 3)  # shared {c}, union {a,b,c}
    from scipy.spatial import distance

    counts = [[10, 0, 5], [0, 10, 5]]
    assert bc[0][1] == pytest.approx(distance.braycurtis(*counts))
    eu = ok(micro.compare(table, options={"metric": "euclidean"}))["distance_matrix"]
    assert eu[0][1] == pytest.approx(math.sqrt(200))


def test_aitchison_distance() -> None:
    # (1,0) and (0,1) with pseudocount 1 -> (2,1) and (1,2): clr = (+-ln2/2); distance = sqrt(2) * ln 2
    d = ok(micro.compare([[1, 0], [0, 1]], options={"metric": "aitchison"}))["distance_matrix"][0][
        1
    ]
    assert d == pytest.approx(math.sqrt(2) * math.log(2))
    # compositionally equivalent profiles (after +1 they are proportional) have distance 0
    same = ok(micro.compare([[0, 1, 3], [1, 3, 7]], options={"metric": "aitchison"}))[
        "distance_matrix"
    ][0][1]
    assert same == pytest.approx(0.0, abs=1e-12)


def test_rarefaction_is_seeded_and_conserves_depth() -> None:
    options = {"entity": "rarefaction", "depth": 10, "seed": 7}
    first = ok(micro.quantify(TABLE, options=options))
    second = ok(micro.quantify(TABLE, options=options))
    assert first == second
    for row, original in zip(
        first["counts"], [[5, 5, 5, 5], [1, 1, 2, 6], [10, 0, 3, 0]], strict=False
    ):
        assert sum(row) == 10
        assert all(r <= o for r, o in zip(row, original, strict=False))
    assert first["counts"][1] == [1, 1, 2, 6]  # depth equals the total: nothing can be dropped
    assert first["counts"][2][1] == 0 and first["counts"][2][3] == 0
    other = ok(micro.quantify(TABLE, options={**options, "seed": 8}))
    assert other["counts"][0] != first["counts"][0] or other["counts"][2] != first["counts"][2]


def test_rarefaction_drops_samples_below_depth() -> None:
    result = ok(micro.quantify(TABLE, options={"entity": "rarefaction", "depth": 15, "seed": 0}))
    # totals are 20, 10 and 13: only s1 reaches a depth of 15
    assert result["dropped_samples"] == ["s2", "s3"] and result["samples"] == ["s1"]
    assert sum(result["counts"][0]) == 15
    with pytest.raises(ConfigurationError):
        micro.quantify(TABLE, options={"entity": "rarefaction", "depth": 15, "drop_short": False})
    with pytest.raises(ConfigurationError):
        micro.quantify(TABLE, options={"entity": "rarefaction"})


def test_core_microbiome_and_prevalence() -> None:
    result = ok(micro.detect(TABLE, options={"prevalence": 0.9}))
    # a and c are present in all three samples, b in two, d in two
    assert result["core_taxa"] == ["a", "c"]
    prev = {r["taxon"]: r["prevalence"] for r in result["prevalence_table"]}
    assert prev == {"a": 1.0, "b": pytest.approx(2 / 3), "c": 1.0, "d": pytest.approx(2 / 3)}
    only = ok(micro.detect(TABLE, options={"entity": "prevalence"}))
    assert "core_taxa" not in only
    strict = ok(micro.detect(TABLE, options={"detection": 0.3, "prevalence": 0.5}))
    # relative abundance > 0.3: a in s3 (10/13), b in none (0.25, 0.1, 0), c in none, d in s2 (0.6)
    assert {r["taxon"]: r["samples_present"] for r in strict["prevalence_table"]} == {
        "a": 1,
        "b": 0,
        "c": 0,
        "d": 1,
    }


@pytest.mark.parametrize(
    ("lineage", "lowest", "genus"),
    [
        ("k__Bacteria; p__Firmicutes; c__Bacilli; g__Lactobacillus; s__", "genus", "Lactobacillus"),
        ("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria", "class", None),
        (
            "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus",
            "genus",
            "Lactobacillus",
        ),
        ("k__Bacteria|p__Firmicutes|c__|g__Unassigned", "phylum", None),
    ],
)
def test_taxonomy_parsing(lineage: str, lowest: str, genus: str | None) -> None:
    row = ok(micro.annotate(lineage))[0]
    assert row["lowest_rank"] == lowest
    assert row["ranks"]["genus"] == genus
    assert row["ranks"]["domain"] == "Bacteria"


def test_taxonomy_mapping_input() -> None:
    rows = ok(micro.annotate({"otu1": "k__Bacteria;p__Firmicutes", "otu2": "k__Archaea"}))
    assert [r["id"] for r in rows] == ["otu1", "otu2"] and rows[1]["lowest_taxon"] == "Archaea"


def test_cooccurrence_network_perfect_anticorrelation_of_two_taxa() -> None:
    # with only two taxa the relative abundances sum to 1, so r = -1 exactly (Pearson and Spearman)
    table = {
        "s1": {"a": 1, "b": 9},
        "s2": {"a": 3, "b": 7},
        "s3": {"a": 5, "b": 5},
        "s4": {"a": 8, "b": 2},
    }
    for method in ("pearson", "spearman"):
        result = ok(micro.build_cooccurrence_network(table, options={"method": method}))
        assert result["edges"] == [
            {"source": "a", "target": "b", "correlation": pytest.approx(-1.0)}
        ]


def test_cooccurrence_spearman_matches_scipy() -> None:
    counts = [[1, 4, 2, 9], [2, 3, 7, 1], [3, 9, 1, 5], [4, 1, 8, 3], [5, 7, 3, 2]]
    result = ok(
        micro.build_cooccurrence_network(
            counts, options={"min_abs_correlation": 0.0, "min_samples_present": 1}
        )
    )
    rel = np.array(counts, dtype=float)
    rel = rel / rel.sum(axis=1, keepdims=True)
    expected = scipy_stats.spearmanr(rel[:, 0], rel[:, 1]).statistic
    edge = next(e for e in result["edges"] if e["source"] == "taxon1" and e["target"] == "taxon2")
    assert edge["correlation"] == pytest.approx(expected)


def test_kmer_ani_matches_mash_formula_by_hand() -> None:
    # k=3: AAAC -> {AAA, AAC}, AACC -> {AAC, ACC}; J = 1/3; D = -ln(2J/(1+J))/3 = ln(2)/3
    pair = ok(micro.calculate_ani(["AAAC", "AACC"], options={"k": 3, "canonical": False}))["pairs"][
        0
    ]
    assert pair["shared_kmers"] == 1 and pair["jaccard"] == pytest.approx(1 / 3)
    assert pair["containment_a_in_b"] == pytest.approx(0.5)
    assert pair["mash_distance"] == pytest.approx(math.log(2) / 3)
    assert pair["ani_estimate"] == pytest.approx(1 - math.log(2) / 3)


def test_kmer_ani_identical_unrelated_and_reverse_complement() -> None:
    genome = "ACGTTGCATGCCGATAGCTTAGGCTAAGCTTCGAT"
    rc = genome[::-1].translate(str.maketrans("ACGT", "TGCA"))
    same = ok(micro.calculate_ani({"g": genome, "h": genome}, options={"k": 5}))["pairs"][0]
    assert same["ani_estimate"] == pytest.approx(1.0)
    flipped = ok(micro.calculate_ani({"g": genome, "r": rc}, options={"k": 5}))["pairs"][0]
    assert flipped["ani_estimate"] == pytest.approx(1.0)  # canonical k-mers ignore orientation
    unrelated = ok(micro.calculate_ani(["AAAAAAAAAA", "CCCCCCCCCC"], options={"k": 4}))["pairs"][0]
    assert unrelated["ani_estimate"] is None and unrelated["jaccard"] == 0


def test_differential_abundance_two_groups_matches_scipy() -> None:
    table = {
        "a": {"x": 1, "y": 9},
        "b": {"x": 2, "y": 8},
        "c": {"x": 3, "y": 7},
        "d": {"x": 8, "y": 2},
        "e": {"x": 9, "y": 1},
        "f": {"x": 7, "y": 3},
    }
    result = ok(micro.calculate_differential_abundance(table, options={"groups": list("AAABBB")}))
    rows = {r["taxon"]: r for r in result["taxa"]}
    # complete separation of 3 vs 3: exact two-sided Mann-Whitney p = 2/20 = 0.1, U(A) = 0 for x, 9 for y
    assert rows["x"]["p_value"] == pytest.approx(0.1) and rows["x"]["u_statistic"] == 0
    assert rows["y"]["u_statistic"] == 9
    assert rows["x"]["mean_a"] == pytest.approx(0.2) and rows["x"]["mean_b"] == pytest.approx(0.8)
    assert rows["x"]["log2_fold_change"] == pytest.approx(2.0, abs=1e-4)  # log2(0.8 / 0.2)
    assert rows["y"]["log2_fold_change"] == pytest.approx(-2.0, abs=1e-4)
    assert rows["x"]["q_value"] == pytest.approx(0.1)  # BH over two equal p-values
    ref = scipy_stats.mannwhitneyu([0.1, 0.2, 0.3], [0.8, 0.9, 0.7], alternative="two-sided")
    assert rows["x"]["p_value"] == pytest.approx(ref.pvalue)


def test_differential_abundance_requires_groups() -> None:
    with pytest.raises(ConfigurationError):
        micro.calculate_differential_abundance(TABLE)
    with pytest.raises(ConfigurationError):
        micro.calculate_differential_abundance(TABLE, options={"groups": ["a", "a", "a"]})


def test_dereplicate_exact_and_reverse_complement() -> None:
    seqs = {"g1": "ACGTACGTAC", "g2": "GTACGTACGT", "g3": "TTTTGGGG", "g4": "acgtacgtac"}
    clusters = ok(micro.dereplicate(seqs))["clusters"]
    assert clusters == [
        {"representative": "g1", "members": ["g1", "g2", "g4"]},
        {"representative": "g3", "members": ["g3"]},
    ]
    strict = ok(micro.dereplicate(seqs, options={"canonical": False}))
    assert strict["cluster_count"] == 3  # g2 (reverse complement of g1) is no longer merged
    by_ani = ok(
        micro.dereplicate(
            {"a": "ACGT" * 30, "b": "ACGT" * 30, "c": "GGCC" * 30},
            options={"method": "ani", "k": 8},
        )
    )
    assert by_ani["cluster_count"] == 2


def test_qc_flags_low_depth_and_empty_taxa() -> None:
    table = {"deep": {"a": 2000, "b": 0, "z": 0}, "shallow": {"a": 5, "b": 3, "z": 0}}
    result = ok(micro.qc(table, options={"min_depth": 1000}))
    assert result["flagged_samples"] == ["shallow"]
    assert result["empty_taxa"] == ["z"]
    rows = {r["sample"]: r for r in result["samples"]}
    assert rows["deep"]["observed"] == 1 and rows["deep"]["zero_fraction"] == pytest.approx(2 / 3)
    assert rows["shallow"]["total"] == 8


def test_table_read_export_roundtrip_and_orientation() -> None:
    text = "# Constructed from biom file\n#OTU ID\ts1\ts2\ntaxA\t1\t0\ntaxB\t3\t4\n"
    table = ok(micro.read(text))
    assert table == {
        "samples": ["s1", "s2"],
        "taxa": ["taxA", "taxB"],
        "counts": [[1.0, 3.0], [0.0, 4.0]],
    }
    exported = ok(micro.export(table))
    assert exported["text"] == "#OTU ID\ts1\ts2\ntaxA\t1\t0\ntaxB\t3\t4\n"
    assert ok(micro.read(exported["text"])) == table
    rows_text = "sample\ttaxA\ttaxB\ns1\t1\t3\ns2\t0\t4\n"
    assert ok(micro.read(rows_text, options={"orientation": "samples_rows"})) == table
    as_json = ok(micro.export(table, options={"format": "json"}))["text"]
    assert '"samples": ["s1", "s2"]' in as_json


def test_export_writes_file_and_refuses_overwrite(tmp_path) -> None:
    target = tmp_path / "t.tsv"
    result = ok(micro.export(TABLE, options={"output": str(target)}))
    assert target.read_text() == result["text"]
    with pytest.raises(ConfigurationError):
        micro.export(TABLE, options={"output": str(target)})


def test_plot_stacked_bars_has_one_segment_per_taxon_and_sample() -> None:
    svg = ok(micro.plot(TABLE, options={"top_n": 4}))["svg"]
    segments = re.findall(
        r'class="bar" data-sample="(\w+)" data-taxon="(\w+)" data-fraction="([\d.]+)"', svg
    )
    assert len(segments) == 4 + 4 + 2  # absent taxa are not drawn: s3 has only a and c
    s1 = [float(f) for s, _, f in segments if s == "s1"]
    assert sum(s1) == pytest.approx(1.0) and all(f == pytest.approx(0.25) for f in s1)


MICRO_BAD_INPUT = [
    (micro.quantify, {"s": {"a": -1}}),
    (micro.quantify, [[1, 2], [3]]),
    (micro.compare, [[1, 2]]),
    (micro.qc, {}),
    (micro.read, "#OTU ID\ts1\ntaxA\tabc\n"),
    (micro.annotate, 42),
]


@pytest.mark.parametrize(("function", "value"), MICRO_BAD_INPUT)
def test_microbiome_rejects_invalid_input(function, value) -> None:
    with pytest.raises(DNAKitError):
        function(value)


def test_microbiome_unknown_entity_and_metric_raise() -> None:
    with pytest.raises(ConfigurationError):
        micro.quantify(TABLE, options={"entity": "amr"})
    with pytest.raises(ConfigurationError):
        micro.compare(TABLE, options={"metric": "unifrac"})


# ============================================================================ design_crispr
GUIDE_A = "GACGTACGTCGTAAGCTTCA"
TARGET = GUIDE_A + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"


def _brute_force_spcas9(sequence: str, length: int = 20) -> set[tuple[str, str, int]]:
    """(guide, strand, start) for every NGG-adjacent protospacer on either strand."""

    comp = str.maketrans("ACGT", "TGCA")
    found = set()
    for strand, text in (("forward", sequence), ("reverse", sequence.translate(comp)[::-1])):
        for i in range(length, len(text) - 2):
            if text[i + 1 : i + 3] == "GG":
                start = i - length if strand == "forward" else len(sequence) - i
                found.add((text[i - length : i], strand, start))
    return found


def test_guide_enumeration_matches_brute_force_and_cut_geometry() -> None:
    guides = ok(
        crispr.design(
            TARGET,
            options={"min_gc": 0.0, "max_gc": 1.0, "max_homopolymer": 100, "avoid_poly_t": False},
        )
    )["guides"]
    assert {(g["sequence"], g["strand"], g["start"]) for g in guides} == _brute_force_spcas9(TARGET)
    forward = next(g for g in guides if g["sequence"] == GUIDE_A)
    assert (forward["start"], forward["end"], forward["pam"], forward["strand"]) == (
        0,
        20,
        "AGG",
        "forward",
    )
    assert forward["cut_site"] == 17  # SpCas9 cuts 3 bp upstream of the PAM
    reverse = [g for g in guides if g["strand"] == "reverse"]
    assert reverse
    for g in reverse:
        assert g["cut_site"] == g["start"] + 3  # PAM sits on the left in forward coordinates


def test_design_strand_and_nuclease_options() -> None:
    loose = {"min_gc": 0, "max_gc": 1, "max_homopolymer": 100, "avoid_poly_t": False}
    forward = ok(crispr.design(TARGET, options={**loose, "strand": "forward"}))["guides"]
    assert forward and {g["strand"] for g in forward} == {"forward"}
    cas12a = ok(
        crispr.design(
            "TTTG" + GUIDE_A + "CCCC",
            options={**loose, "nuclease": "AsCas12a", "strand": "forward"},
        )
    )["guides"]
    # Cas12a: the PAM (TTTV) comes first, the 20-nt guide follows it; no Cas9 cut geometry is claimed
    assert (
        cas12a[0]["sequence"] == GUIDE_A
        and cas12a[0]["pam"] == "TTTG"
        and cas12a[0]["cut_site"] is None
    )
    with pytest.raises(ConfigurationError):
        crispr.design(TARGET, options={"nuclease": "Cas99"})


def test_design_filters_gc_homopolymer_and_poly_t() -> None:
    target = "AAAAAAAAAAAAAAAAAAAATGGGACTGACTGACTGACTGACTGACTGACTGCAGG"
    kept = ok(crispr.design(target))["guides"]
    assert all(0.4 <= g["gc"] <= 0.7 and g["longest_homopolymer"] <= 4 for g in kept)
    assert "A" * 20 not in {g["sequence"] for g in kept}
    all_guides = ok(
        crispr.design(target, options={"min_gc": 0, "max_gc": 1, "max_homopolymer": 100})
    )["guides"]
    assert "A" * 20 in {g["sequence"] for g in all_guides}


def test_design_with_references_counts_off_target_sites() -> None:
    near = GUIDE_A[:5] + "T" + GUIDE_A[6:] + "TGG"  # one mismatch (position 6) + valid PAM
    refs = {"chrA": "CCCCC" + GUIDE_A + "AGG" + "TT" + near + "CCCCC"}
    guides = ok(
        crispr.design(
            TARGET,
            options={
                "references": refs,
                "min_gc": 0,
                "max_gc": 1,
                "max_homopolymer": 100,
                "avoid_poly_t": False,
            },
        )
    )["guides"]
    first = next(g for g in guides if g["sequence"] == GUIDE_A)
    assert first["exact_sites"] == 1 and first["near_sites"] == 1 and first["off_target_sites"] == 1


def test_donor_template_geometry() -> None:
    target = "A" * 30 + "C" * 20 + "G" * 30
    donor = ok(
        crispr.design(
            target,
            options={
                "entity": "donor_template",
                "cut_site": 50,
                "insert": "TTTT",
                "delete_length": 5,
                "left_arm": 10,
                "right_arm": 8,
            },
        )
    )
    assert (
        donor["donor"] == "C" * 10 + "TTTT" + "G" * 5 + "G" * 3
    )  # arms flank the 5 deleted bases (G*5 at 50..55)
    assert donor["left_arm"] == {"start": 40, "end": 50, "sequence": "C" * 10}
    assert donor["right_arm"]["start"] == 55 and donor["deleted"] == "GGGGG"
    assert donor["edited_sequence"] == "A" * 30 + "C" * 20 + "TTTT" + "G" * 25
    with pytest.raises(ConfigurationError):
        crispr.design(target, options={"entity": "donor_template", "cut_site": 5, "left_arm": 10})


def test_simulate_hdr_replaces_region_between_arms() -> None:
    target = "AAAACCCCGGGGTTTT"
    donor = "CCCC" + "AA" + "GGGG"
    result = ok(crispr.simulate(target, options={"donor": donor, "left_arm": 4, "right_arm": 4}))
    assert result["edited_sequence"] == "AAAACCCCAAGGGGTTTT" and result["length_change"] == 2
    assert (
        result["region"] == {"start": 8, "end": 8}
        and result["inserted"] == "AA"
        and result["deleted"] == ""
    )
    replaced = ok(
        crispr.simulate(
            "AAAACCCCTTGGGGTTTT", options={"donor": donor, "left_arm": 4, "right_arm": 4}
        )
    )
    assert (
        replaced["deleted"] == "TT"
        and replaced["edited_sequence"] == "AAAACCCCAAGGGGTTTT"
        and replaced["length_change"] == 0
    )
    with pytest.raises(InputFormatError):
        crispr.simulate("AAAACCCC", options={"donor": donor, "left_arm": 4, "right_arm": 4})


def test_simulate_roundtrips_design_donor() -> None:
    target = "ACGTTGCAAGGCTTAAGCTAGCTAGGATCCATGCAATGCCGGATATCGCGTA"
    design = ok(
        crispr.design(
            target,
            options={
                "entity": "donor_template",
                "cut_site": 25,
                "insert": "GGG",
                "left_arm": 10,
                "right_arm": 10,
            },
        )
    )
    edited = ok(
        crispr.simulate(target, options={"donor": design["donor"], "left_arm": 10, "right_arm": 10})
    )
    assert edited["edited_sequence"] == design["edited_sequence"]


def test_offtarget_hamming_search_counts_by_mismatch_and_requires_pam() -> None:
    guide = GUIDE_A
    one_mm = guide[:9] + ("A" if guide[9] != "A" else "C") + guide[10:]
    three_mm = "".join(
        ("A" if c != "A" else "C") if i in (0, 1, 2) else c for i, c in enumerate(guide)
    )
    five_mm = "".join(
        ("A" if c != "A" else "C") if i in (0, 4, 8, 12, 16) else c for i, c in enumerate(guide)
    )
    reference = "TT".join(
        [guide + "TGG", one_mm + "CGG", three_mm + "AGG", five_mm + "TGG", one_mm + "AAA"]
    )
    result = ok(
        crispr.scan_grna_offtargets(
            guide, options={"references": {"r": reference}, "max_mismatches": 4}
        )
    )
    hits = result["guides"][0]["hits"]
    # exact + 1 mismatch + 3 mismatches; five-mismatch site exceeds the budget; AAA is not an NGG PAM
    assert [h["mismatches"] for h in hits] == [0, 1, 3]
    assert result["guides"][0]["counts_by_mismatch"] == {"0": 1, "1": 1, "3": 1}
    assert hits[1]["mismatch_positions"] == [10] and hits[2]["mismatch_positions"] == [1, 2, 3]
    assert hits[0]["start"] == 0 and hits[0]["end"] == 20
    no_pam = ok(
        crispr.scan_grna_offtargets(
            guide, options={"references": {"r": reference}, "max_mismatches": 4, "pam": None}
        )
    )
    assert [h["mismatches"] for h in no_pam["guides"][0]["hits"]] == [0, 1, 1, 3]


def test_offtarget_search_reverse_strand_coordinates() -> None:
    guide = GUIDE_A
    comp = str.maketrans("ACGT", "TGCA")
    site = (guide + "TGG").translate(comp)[::-1]  # reverse-complement of protospacer+PAM
    reference = "GGGGG" + site + "GGGGG"
    hit = ok(crispr.scan_grna_offtargets(guide, options={"references": {"r": reference}}))[
        "guides"
    ][0]["hits"]
    assert len(hit) == 1 and hit[0]["strand"] == "-" and hit[0]["mismatches"] == 0
    assert (hit[0]["start"], hit[0]["end"]) == (
        5 + 3,
        5 + 3 + 20,
    )  # PAM (3 nt) comes first on the forward strand


@pytest.mark.parametrize(
    ("guide", "expected"),
    [
        ("GGGGGGGGGGGGGGGGGGGG", 100 - 40 - 15),  # GC 100 % (> 80): -40; G homopolymer: -15
        ("ACGTACGTTTTTACGTACGT", 100 - 40),  # GC 40 % ok; TTTT: -40
        ("GACGTACGTCGTAAGCTTCA", 100),  # GC 50 %, no runs
        ("ATATATATATATATATATAT", 100 - 40),  # GC 0 %
        ("GCGCGCGCGCGCGCGCGCAT", 100 - 15),  # GC 90 % -> > 80 => -40 ... see below
    ],
)
def test_heuristic_proxy_score_rules(guide: str, expected: float) -> None:
    gc = (guide.count("G") + guide.count("C")) / len(guide)
    if gc > 0.8:
        expected = 100 - 40 - (15 if re.search(r"(A{4}|C{4}|G{4})", guide) else 0)
    scored = crispr.score(guide)
    row = ok(scored)["scores"][0]
    assert row["proxy_score"] == expected
    assert scored["warnings"] and "heuristic" in scored["warnings"][0].lower()


def test_score_rejects_bad_guides() -> None:
    with pytest.raises(InputFormatError):
        crispr.score("ACGU" + "X")
    with pytest.raises(InputFormatError):
        crispr.score("ACGT")  # too short


def test_rank_orders_by_proxy_score_or_custom_field() -> None:
    guides = ["GGGGGGGGGGGGGGGGGGGG", GUIDE_A, "ACGTACGTTTTTACGTACGT"]
    ranked = ok(crispr.rank(guides))["ranked"]
    assert [g["sequence"] for g in ranked] == [
        GUIDE_A,
        "ACGTACGTTTTTACGTACGT",
        "GGGGGGGGGGGGGGGGGGGG",
    ]
    assert [g["rank"] for g in ranked] == [1, 2, 3]
    records = [
        {"id": "x", "sequence": GUIDE_A, "off_target_sites": 5},
        {"id": "y", "sequence": GUIDE_A, "off_target_sites": 1},
    ]
    low_first = ok(
        crispr.rank(records, options={"by": "off_target_sites", "ascending": True, "top_n": 1})
    )["ranked"]
    assert [g["id"] for g in low_first] == ["y"]
    with pytest.raises(InputFormatError):
        crispr.rank(records, options={"by": "missing_field"})


def test_filter_reports_reasons() -> None:
    guides = {
        "ok": GUIDE_A,
        "gc_low": "ATATATATATATATATATAT",
        "polyT": "GACTTTTCGACTGACTGACT",
        "run": "GACTGGGGGACTGACTGACT",
        "site": "GACGAATTCGTAAGCTTCAG",
    }
    result = ok(crispr.filter(guides, options={"avoid_enzymes": ["EcoRI"]}))
    reasons = {r["id"]: r["failed"] for r in result["rejected"]}
    assert [k["id"] for k in result["kept"]] == ["ok"]
    assert reasons["gc_low"] == ["gc_below_min"]
    assert "poly_t" in reasons["polyT"] and "homopolymer_run" in reasons["run"]
    assert reasons["site"] == ["excluded_motif:EcoRI"]


def test_validate_checks_pam_and_length() -> None:
    result = ok(
        crispr.validate(
            [GUIDE_A, GUIDE_A[:-1], "ACGTN" + "A" * 15, "CATCATCATCATCATCATCA"],
            options={"target": TARGET},
        )
    )
    rows = {r["id"]: r for r in result["guides"]}
    assert rows["guide1"]["valid"] and rows["guide1"]["target_sites"][0] == {
        "strand": "forward",
        "start": 0,
        "end": 20,
        "pam": "AGG",
    }
    assert "length_19_expected_20" in rows["guide2"]["errors"]
    assert rows["guide3"]["errors"] == ["invalid_characters"]
    assert "no_pam_adjacent_site_in_target" in rows["guide4"]["errors"]
    assert result["all_valid"] is False


def test_gibson_style_fragment_assembly() -> None:
    overlap = "ACGTACGTACGTACGTAAAA"
    first, second = "AAAAAC" + overlap, overlap + "CCCCCT"
    result = ok(crispr.assemble_dna_fragments([first, second], options={"min_overlap": 15}))
    assert result["product"] == "AAAAAC" + overlap + "CCCCCT" and result["length"] == 32
    assert result["junctions"][0]["overlap"] == overlap
    with pytest.raises(DNAKitError):
        crispr.assemble_dna_fragments(["AAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCC"])
    circular = ok(
        crispr.assemble_dna_fragments(
            ["GGGGGGGGGGAAAAATTTTT", "AAAAATTTTTCCCCCCCCCCGGGGG"],
            options={"min_overlap": 5, "circularize": True},
        )
    )
    assert circular["circular"] is True and circular["junctions"][-1]["overlap"] == "GGGGG"


def test_primer_dimer_runs() -> None:
    # AAAAGCGC is its own reverse complement at GCGC: 3' run = 4 (last 4 bases GCGC occur in rc = GCGCTTTT)
    row = ok(crispr.check_primer_dimer("AAAAGCGC"))["pairs"][0]
    assert row["longest_complementary_run"] == 4 and row["three_prime_run_a"] == 4
    assert row["three_prime_dimer"] is True
    # AAAAAAAA and CCCCCCCC share no complementary stretch at all
    none = next(
        p
        for p in ok(crispr.check_primer_dimer(["AAAAAAAA", "CCCCCCCC"]))["pairs"]
        if p["a"] != p["b"]
    )
    assert none["longest_complementary_run"] == 0 and none["three_prime_dimer"] is False
    cross = ok(crispr.check_primer_dimer(["AAAAGGGG", "CCCCTTTT"]))["pairs"]
    pair = next(p for p in cross if p["a"] != p["b"])
    assert pair["longest_complementary_run"] == 8  # CCCCTTTT is the reverse complement of AAAAGGGG
    assert pair["three_prime_run_a"] == 8


def test_primer_hairpin_detection() -> None:
    row = ok(crispr.check_primer_secondary_structure("GCGCAAAAGCGC"))["primers"][0]
    assert row["longest_stem"] == 4 and row["hairpins"][0]["loop_length"] == 4
    assert row["hairpins"][0]["left_start"] == 0 and row["hairpins"][0]["right_start"] == 8
    assert row["three_prime_hairpin"] is True
    none = ok(crispr.check_primer_secondary_structure("ACACACACACACACAC", options={"min_stem": 5}))[
        "primers"
    ][0]
    assert none["hairpins"] == [] and none["longest_stem"] == 0


def test_primer_specificity_counts_sites_on_both_strands() -> None:
    primer = "GAATTCGGGACG"
    template = (
        "TTT" + primer + "AAAA" + primer[::-1].translate(str.maketrans("ACGT", "TGCA")) + "CC"
    )
    row = ok(crispr.check_primer_specificity(primer, options={"templates": {"t": template}}))[
        "primers"
    ][0]
    assert row["total_hits"] == 2 and row["specific"] is False
    assert sorted(h["strand"] for h in row["hits"]) == ["forward", "reverse"]
    single = ok(
        crispr.check_primer_specificity(primer, options={"templates": {"t": "TTT" + primer + "CC"}})
    )["primers"][0]
    assert single["specific"] is True and single["hits"][0]["start"] == 3
    approx = ok(
        crispr.check_primer_specificity(
            primer,
            options={
                "templates": {"t": "TTT" + primer[:5] + "A" + primer[6:] + "CC"},
                "max_mismatches": 1,
            },
        )
    )["primers"][0]
    assert approx["hits"][0]["mismatch_positions"] == [5]
    with pytest.raises(ConfigurationError):
        crispr.check_primer_specificity(primer)


def test_remove_restriction_sites_edits_only_the_site() -> None:
    sequence = "TTTTGAATTCAAAACCCGGATCCTTTT"
    result = ok(crispr.remove_restriction_sites(sequence, options={"enzymes": ["EcoRI", "BamHI"]}))
    assert (
        result["sites_before"] == 2 and result["sites_after"] == 0 and result["removed_all"] is True
    )
    assert len(result["sequence"]) == len(sequence)
    assert "GAATTC" not in result["sequence"] and "GGATCC" not in result["sequence"]
    for change in result["changes"]:
        assert sequence[change["position"]] == change["before"]
        assert 4 <= change["position"] < 10 or 16 <= change["position"] < 22
    custom = ok(crispr.remove_restriction_sites("AAGCTTAA", options={"sites": ["AAGCTT"]}))
    assert "AAGCTT" not in custom["sequence"]
    with pytest.raises(ConfigurationError):
        crispr.remove_restriction_sites(sequence)


def test_screen_design_library() -> None:
    other = "TGCATGCCGATCGATCGGACTGGCATCGATCGGTACGCTAGCTAGGCCTAGCGGATCGATAGGC"
    result = ok(
        crispr.run_crispr_screen_design(
            {"geneA": TARGET, "geneB": other},
            options={"guides_per_target": 2, "min_cut_spacing": 5},
        )
    )
    by_target = {t: [g for g in result["library"] if g["target"] == t] for t in ("geneA", "geneB")}
    assert len(by_target["geneB"]) == 2
    for rows in by_target.values():
        scores = [g["proxy_score"] for g in rows]
        assert scores == sorted(scores, reverse=True)
        cuts = [g["cut_site"] for g in rows]
        assert all(abs(a - b) >= 5 for i, a in enumerate(cuts) for b in cuts[i + 1 :])
    assert result == ok(
        crispr.run_crispr_screen_design(
            {"geneA": TARGET, "geneB": other},
            options={"guides_per_target": 2, "min_cut_spacing": 5},
        )
    )
    greedy = ok(
        crispr.run_crispr_screen_design({"x": "ACGT" * 5}, options={"guides_per_target": 3})
    )
    assert greedy["targets_short"] == [{"target": "x", "requested": 3, "designed": 0}]


def test_export_formats(tmp_path) -> None:
    records = [
        {"id": "g1", "sequence": GUIDE_A, "strand": "forward"},
        {"id": "g2", "sequence": "ACGTACGTACGTACGTACGT", "strand": "reverse"},
    ]
    tsv = ok(crispr.export(records))
    assert tsv["text"].splitlines()[0] == "id\tsequence\tstrand" and tsv["count"] == 2
    assert (
        ok(crispr.export(records, options={"format": "fasta"}))["text"]
        == f">g1\n{GUIDE_A}\n>g2\nACGTACGTACGTACGTACGT\n"
    )
    assert (
        ok(crispr.export(records, options={"format": "csv"}))["text"].splitlines()[1]
        == f"g1,{GUIDE_A},forward"
    )
    out = tmp_path / "g.tsv"
    ok(crispr.export(records, options={"output": str(out)}))
    assert out.read_text().count("\n") == 3


def test_guide_map_plot() -> None:
    guides = ok(
        crispr.design(
            TARGET,
            options={"min_gc": 0, "max_gc": 1, "max_homopolymer": 100, "avoid_poly_t": False},
        )
    )["guides"]
    svg = ok(crispr.plot({"length": len(TARGET), "guides": guides}))["svg"]
    assert len(re.findall(r'class="guide"', svg)) == len(guides)
    assert 'data-start="0" data-end="20" data-strand="forward"' in svg
    with pytest.raises(InputFormatError):
        crispr.plot({"guides": [{"start": 5, "end": 5}]})


def test_crispr_still_exposes_the_original_specific_functions() -> None:
    assert (
        ok(crispr.evaluate_primer_pair(["ACGTACGT", "ACGTACGT"], options={"method": "wallace"}))[
            "forward"
        ]["tm_celsius"]
        == 24
    )  # Wallace rule 2*(A+T) + 4*(G+C)
    assert (
        ok(crispr.encode_data_into_dna("hi")) == "CGGACGGC"
    )  # 2 bits per base, h=0x68=01 10 10 00


# ============================================================================ epigenome
BED = "chr1\t100\t200\tp1\nchr1\t150\t250\tp2\nchr1\t250\t300\tp3\nchr2\t5\t10\tq"
COV = "chr1\t10\t10\t50.0\t5\t5\nchr1\t20\t20\t100\t10\t0\nchrL\t5\t5\t1\t1\t99\n"


def test_merge_matches_bedtools_semantics() -> None:
    merged = ok(epi.merge(BED))["intervals"]
    assert merged == [
        {"chrom": "chr1", "start": 100, "end": 300, "n_merged": 3},
        {"chrom": "chr2", "start": 5, "end": 10, "n_merged": 1},
    ]
    gapped = ok(epi.merge("chr1\t0\t10\nchr1\t15\t20\nchr1\t40\t50", options={"distance": 5}))[
        "intervals"
    ]
    assert [(m["start"], m["end"], m["n_merged"]) for m in gapped] == [(0, 20, 2), (40, 50, 1)]
    unmerged = ok(epi.merge("chr1\t0\t10\nchr1\t11\t20"))["intervals"]
    assert len(unmerged) == 2  # distance 0 merges overlapping and book-ended intervals only


def test_read_formats_and_errors() -> None:
    assert ok(epi.read(BED))["count"] == 4
    bg = ok(
        epi.read(
            "track type=bedGraph\nchr1\t0\t10\t1.5\n#c\nchr1\t10\t20\t2\n",
            options={"format": "bedgraph"},
        )
    )
    assert [r["value"] for r in bg["records"]] == [1.5, 2.0]
    cov = ok(epi.read(COV, options={"format": "bismark_cov"}))["records"]
    assert cov[0] == {"chrom": "chr1", "start": 10, "methylated": 5, "unmethylated": 5}
    for bad, fmt in [
        ("chr1\t10\t5\n", "bed"),
        ("chr1\tx\t5\n", "bed"),
        ("chr1\t0\t10\n", "bedgraph"),
        ("chr1\t1\t1\t5.0\t3\n", "bismark_cov"),
    ]:
        with pytest.raises(InputFormatError):
            epi.read(bad, options={"format": fmt})


def test_methylation_levels_per_site_and_region() -> None:
    sites = ok(epi.calculate(COV, options={"min_coverage": 11}))
    assert [(s["chrom"], s["coverage"], s["level"]) for s in sites["sites"]] == [
        ("chrL", 100, 0.01)
    ]
    assert sites["dropped_low_coverage"] == 2
    every = ok(epi.calculate(COV))["sites"]
    assert [s["level"] for s in every] == [0.5, 1.0, 0.01]
    regions = ok(
        epi.calculate(
            COV, options={"one_based": True, "regions": "chr1\t0\t15\nchr1\t0\t30\nchr1\t100\t200"}
        )
    )["regions"]
    # one-based starts 10 and 20 -> 0-based 9 and 19; region 0-15 contains site 9 only (5 of 10 reads methylated)
    assert regions[0]["n_sites"] == 1 and regions[0]["weighted_level"] == 0.5
    assert regions[1]["coverage"] == 20 and regions[1]["weighted_level"] == pytest.approx(15 / 20)
    assert regions[1]["mean_site_level"] == pytest.approx((0.5 + 1.0) / 2)
    assert regions[2]["n_sites"] == 0 and regions[2]["weighted_level"] is None


def test_weighted_level_differs_from_mean_of_site_levels() -> None:
    rows = [
        {"chrom": "c", "start": 1, "methylated": 1, "unmethylated": 1},
        {"chrom": "c", "start": 2, "methylated": 90, "unmethylated": 10},
    ]
    region = ok(epi.calculate(rows, options={"regions": "c\t0\t10"}))["regions"][0]
    assert region["weighted_level"] == pytest.approx(91 / 102)
    assert region["mean_site_level"] == pytest.approx((0.5 + 0.9) / 2)


def test_cpg_observed_expected_by_hand() -> None:
    # CGCGAT: L=6, C=2, G=2, CpG=2 -> 2*6/(2*2) = 3.0; GC fraction 4/6
    row = ok(epi.calculate("CGCGAT", options={"entity": "cpg_observed_expected"}))["sequences"][0]
    assert (row["c"], row["g"], row["cpg"]) == (2, 2, 2)
    assert row["observed_expected"] == pytest.approx(3.0) and row["gc_fraction"] == pytest.approx(
        4 / 6
    )
    none = ok(epi.calculate("AAAATT", options={"entity": "cpg_observed_expected"}))["sequences"][0]
    assert none["observed_expected"] is None
    windows = ok(
        epi.calculate(
            "CGCGATAT", options={"entity": "cpg_observed_expected", "window": 4, "step": 4}
        )
    )["sequences"][0]["windows"]
    assert [w["cpg"] for w in windows] == [2, 0]


def test_cpg_islands_found_by_delegation() -> None:
    islands = ok(epi.find("AT" * 150 + "CG" * 100 + "AT" * 150, options={"entity": "cpg_islands"}))[
        "islands"
    ]
    # CG repeats occupy [300, 500); the merged qualifying 200-bp windows (>= 50 % GC, O/E >= 0.6) span [200, 600)
    assert (islands[0]["start"], islands[0]["end"], islands[0]["length"]) == (200, 600, 400)
    assert islands[0]["observed_expected"] >= 0.6 and len(islands) == 1
    pure = ok(epi.find("CG" * 100, options={"entity": "cpg_islands"}))["islands"]
    assert (pure[0]["start"], pure[0]["end"]) == (0, 200) and pure[0][
        "observed_expected"
    ] == pytest.approx(2.0)


def test_differential_methylation_fisher_matches_scipy() -> None:
    rows = [
        {
            "chrom": "chr1",
            "start": 0,
            "end": 100,
            "methylated_a": 10,
            "unmethylated_a": 90,
            "methylated_b": 30,
            "unmethylated_b": 70,
        },
        {
            "chrom": "chr1",
            "start": 100,
            "end": 200,
            "methylated_a": 50,
            "unmethylated_a": 50,
            "methylated_b": 52,
            "unmethylated_b": 48,
        },
        {
            "chrom": "chr1",
            "start": 200,
            "end": 300,
            "methylated_a": 1,
            "unmethylated_a": 2,
            "methylated_b": 3,
            "unmethylated_b": 4,
        },
    ]
    result = ok(epi.find(rows))
    assert result["tested"] == 2 and result["excluded_low_coverage"] == 1
    first = result["regions"][0]
    ref = scipy_stats.fisher_exact([[10, 90], [30, 70]])
    assert first["p_value"] == pytest.approx(ref.pvalue, rel=1e-9) and first[
        "odds_ratio"
    ] == pytest.approx(ref.statistic)
    assert first["delta"] == pytest.approx(0.2) and first["significant"] is True
    second = result["regions"][1]
    assert second["p_value"] == pytest.approx(scipy_stats.fisher_exact([[50, 50], [52, 48]]).pvalue)
    assert second["q_value"] == pytest.approx(
        second["p_value"]
    )  # largest p keeps its value under BH
    assert first["q_value"] == pytest.approx(min(1.0, first["p_value"] * 2))


@pytest.mark.parametrize(
    "table", [[3, 1, 1, 3], [0, 5, 5, 0], [7, 2, 4, 9], [2, 8, 2, 8], [0, 0, 3, 4]]
)
def test_fisher_exact_reference_tables(table: list[int]) -> None:
    from dnakit.epigenome._impl_local import fisher_exact

    a, b, c, d = table
    assert fisher_exact(a, b, c, d)["p_value"] == pytest.approx(
        scipy_stats.fisher_exact([[a, b], [c, d]]).pvalue, rel=1e-9
    )


def test_fisher_hand_value() -> None:
    # classic 2x2 [[3,1],[1,3]]: tables with probability <= 16/70 are 0,1(=16/70),3,4 -> (1+16+16+1)/70 = 34/70
    from dnakit.epigenome._impl_local import fisher_exact

    assert fisher_exact(3, 1, 1, 3)["p_value"] == pytest.approx(34 / 70)


def test_differential_regions_from_sites_with_regions() -> None:
    sites = [
        {
            "chrom": "c",
            "start": s,
            "methylated_a": 1,
            "unmethylated_a": 9,
            "methylated_b": 8,
            "unmethylated_b": 2,
        }
        for s in (1, 2, 3)
    ]
    sites.append(
        {
            "chrom": "c",
            "start": 50,
            "methylated_a": 5,
            "unmethylated_a": 5,
            "methylated_b": 5,
            "unmethylated_b": 5,
        }
    )
    result = ok(epi.find(sites, options={"regions": "c\t0\t10\nc\t40\t60"}))
    top = result["regions"][0]
    assert (top["start"], top["coverage_a"], top["level_a"], top["level_b"]) == (
        0,
        30,
        pytest.approx(0.1),
        pytest.approx(0.8),
    )


def test_consensus_peaks_by_sample_depth() -> None:
    sets = {"s1": "chr1\t0\t100", "s2": "chr1\t50\t150", "s3": "chr1\t120\t200"}
    peaks = ok(epi.call(sets, options={"min_samples": 2}))["peaks"]
    assert peaks == [
        {"chrom": "chr1", "start": 50, "end": 100},
        {"chrom": "chr1", "start": 120, "end": 150},
    ]
    assert ok(epi.call(sets, options={"min_samples": 3}))["peaks"] == []
    union = ok(epi.call(sets, options={"min_samples": 1}))["peaks"]
    assert union == [{"chrom": "chr1", "start": 0, "end": 200}]
    one_sample_twice = ok(
        epi.call(
            {"s1": "chr1\t0\t100\nchr1\t50\t150", "s2": "chr1\t500\t600"},
            options={"min_samples": 2},
        )
    )["peaks"]
    assert one_sample_twice == []  # overlapping peaks of one sample vote once


def test_interval_overlap_jaccard() -> None:
    result = ok(
        epi.calculate(
            "chr1\t0\t100\nchr1\t200\t300",
            options={"entity": "interval_overlap", "b": "chr1\t50\t250"},
        )
    )
    assert result["overlap_bp"] == 100 and result["union_bp"] == 200 + 200 - 100
    assert (
        result["jaccard"] == pytest.approx(100 / 300) and result["a_intervals_overlapping_b"] == 2
    )


def test_signal_normalisation_hand_values() -> None:
    cpm = [r["value"] for r in ok(epi.normalize([10, 30, 60]))["records"]]
    assert cpm == [100000.0, 300000.0, 600000.0]
    bedgraph = "chr1\t0\t1000\t10\nchr1\t1000\t3000\t30"
    rpkm = [r["value"] for r in ok(epi.normalize(bedgraph, options={"method": "rpkm"}))["records"]]
    assert rpkm == [
        pytest.approx(250000.0),
        pytest.approx(375000.0),
    ]  # cpm 250000 / 1 kb, 750000 / 2 kb
    z = [r["value"] for r in ok(epi.normalize([1, 2, 3], options={"method": "zscore"}))["records"]]
    assert z == [pytest.approx(-math.sqrt(1.5)), 0.0, pytest.approx(math.sqrt(1.5))]
    assert [
        r["value"] for r in ok(epi.normalize([2, 4, 6], options={"method": "minmax"}))["records"]
    ] == [0.0, 0.5, 1.0]
    log = ok(epi.normalize([1, 1], options={"method": "log2cpm"}))["records"][0]["value"]
    assert log == pytest.approx(math.log2(500000 + 1))
    with pytest.raises(InputFormatError):
        epi.normalize([5, 5], options={"method": "zscore"})
    with pytest.raises(InputFormatError):
        epi.normalize([0, 0])


def test_chromatin_state_segments() -> None:
    result = ok(
        epi.build_chromatin_state_segments(
            {"chr1": ["E1", "E1", "E2", "E1"]}, options={"bin_size": 200}
        )
    )
    assert result["segments"] == [
        {"chrom": "chr1", "start": 0, "end": 400, "state": "E1", "n_bins": 2},
        {"chrom": "chr1", "start": 400, "end": 600, "state": "E2", "n_bins": 1},
        {"chrom": "chr1", "start": 600, "end": 800, "state": "E1", "n_bins": 1},
    ]
    bed = "chr1\t0\t100\tA\nchr1\t100\t200\tA\nchr1\t300\t400\tA\n"
    assert ok(epi.build_chromatin_state_segments(bed))["segment_count"] == 2
    assert (
        ok(epi.build_chromatin_state_segments(bed, options={"max_gap": 100}))["segment_count"] == 1
    )
    with pytest.raises(InputFormatError):
        epi.build_chromatin_state_segments("chr1\t0\t100\n")


GENES = [
    {"chrom": "chr1", "start": 5000, "end": 9000, "name": "G1", "strand": "+"},
    {"chrom": "chr1", "start": 0, "end": 1500, "name": "G2", "strand": "-"},
]


def test_enhancer_gene_links_nearest_and_window() -> None:
    nearest = ok(epi.link_enhancers_to_genes("chr1\t1000\t1200\te1", options={"genes": GENES}))[
        "links"
    ][0]
    # midpoint 1100; G2 is on the - strand so its TSS is end-1 = 1499 (distance 399); G1 TSS 5000 (3900)
    assert nearest["midpoint"] == 1100
    assert [(g["gene"], g["distance"], g["tss"]) for g in nearest["genes"]] == [("G2", 399, 1499)]
    window = ok(
        epi.link_enhancers_to_genes(
            "chr1\t1000\t1200\te1",
            options={"genes": GENES, "method": "window", "max_distance": 4000},
        )
    )["links"][0]
    assert [(g["gene"], g["distance"]) for g in window["genes"]] == [("G2", 399), ("G1", 3900)]
    far = ok(
        epi.link_enhancers_to_genes(
            "chr1\t1000\t1200\te1", options={"genes": GENES, "max_distance": 100}
        )
    )["links"][0]
    assert far["genes"] == []
    other_chrom = ok(epi.link_enhancers_to_genes("chr9\t1000\t1200", options={"genes": GENES}))[
        "links"
    ][0]
    assert other_chrom["genes"] == []
    with pytest.raises(ConfigurationError):
        epi.link_enhancers_to_genes("chr1\t1000\t1200")


def test_peak_annotation_overlap_and_nearest() -> None:
    features = "chr1\t150\t300\tgA\nchr1\t500\t600\tgB\nchr2\t0\t10\tgC"
    rows = ok(
        epi.annotate(
            "chr1\t100\t200\nchr1\t350\t400\nchr1\t100000\t100100", options={"features": features}
        )
    )["regions"]
    assert rows[0]["overlapping"] == ["gA"] and rows[0]["nearest"] == {
        "name": "gA",
        "distance": 0,
        "side": "overlap",
    }
    assert rows[1]["overlapping"] == [] and rows[1]["nearest"] == {
        "name": "gA",
        "distance": 50,
        "side": "upstream",
    }
    assert rows[2]["nearest"]["name"] == "gB" and rows[2]["nearest"]["distance"] == 100000 - 600
    capped = ok(
        epi.annotate("chr1\t100000\t100100", options={"features": features, "max_distance": 10})
    )["regions"][0]
    assert capped["nearest"] is None


def test_methylation_qc_with_control_chromosome() -> None:
    result = ok(epi.qc(COV, options={"control_chroms": ["chrL"], "min_coverage": 10}))
    assert (
        result["sites"] == 2
        and result["mean_coverage"] == 10.0
        and result["median_coverage"] == 10.0
    )
    assert result["global_methylation"] == pytest.approx(15 / 20)  # control reads excluded
    assert result["conversion_rate"] == pytest.approx(0.99)  # 1 of 100 control reads is methylated
    assert result["fraction_sites_min_coverage"] == 1.0
    low = ok(epi.qc("chr1\t1\t1\t0\t0\t3\nchr1\t2\t2\t0\t0\t7\n", options={"min_coverage": 5}))
    assert low["fraction_sites_min_coverage"] == 0.5 and low["median_coverage"] == 5.0


def test_signal_track_plot() -> None:
    svg = ok(epi.plot("chr1\t0\t10\t1\nchr1\t10\t20\t2\nchr2\t0\t5\t9", options={"chrom": "chr1"}))[
        "svg"
    ]
    bars = re.findall(r'data-start="(\d+)" data-end="(\d+)" data-value="([\d.]+)"', svg)
    assert bars == [("0", "10", "1.0"), ("10", "20", "2.0")]
    with pytest.raises(InputFormatError):
        epi.plot("chr1\t0\t10\t1", options={"chrom": "chrX"})


# ============================================================================ regulation
# a 3-column motif that prefers A, C, G: p = 0.7 for the preferred base, 0.1 otherwise, uniform background 0.25
PWM = {"A": [0.7, 0.1, 0.1], "C": [0.1, 0.7, 0.1], "G": [0.1, 0.1, 0.7], "T": [0.1, 0.1, 0.1]}
HIGH = math.log2(0.7 / 0.25)
LOW = math.log2(0.1 / 0.25)


def test_pwm_score_best_window_by_hand() -> None:
    best = ok(reg.score("TTACGTT", options={"motifs": PWM}))["sequences"][0]["motifs"][0]["best"]
    assert (best["start"], best["end"], best["strand"]) == (2, 5, "+")
    assert best["score"] == pytest.approx(3 * HIGH) and best["relative_score"] == pytest.approx(1.0)
    # TTACGTT also contains CGT = reverse complement of ACG at 3..6 with the same score; the + strand hit is reported first
    only_reverse = ok(reg.score("TTCGTTT", options={"motifs": PWM}))["sequences"][0]["motifs"][0][
        "best"
    ]
    assert only_reverse["strand"] == "-" and only_reverse["score"] == pytest.approx(3 * HIGH)
    forward_only = ok(reg.score("TTCGTTT", options={"motifs": PWM, "strand": "forward"}))[
        "sequences"
    ][0]["motifs"][0]["best"]
    assert forward_only["score"] < 3 * HIGH
    short = ok(reg.score("AC", options={"motifs": PWM}))["sequences"][0]["motifs"][0]["best"]
    assert short is None


def test_pwm_scores_agree_with_dnakit_scan_pwm() -> None:
    from dnakit.core import DNASequence
    from dnakit.patterns import PWM as PatternPWM
    from dnakit.patterns import scan_pwm

    sequence = "GATTACAACGTTCAGGACGTAC"
    reference = scan_pwm(DNASequence(sequence), PatternPWM("m", PWM), threshold=-50, verbose=False)
    mine = ok(reg.annotate(sequence, options={"motifs": PWM, "threshold": -50}))["sequences"][0][
        "hits"
    ]
    key = lambda start, strand, score: (start, strand, round(score, 9))  # noqa: E731
    assert {key(h.symbol_location.start, h.strand.value[0], h.score) for h in reference.hits} == {
        key(h["start"], "f" if h["strand"] == "+" else "r", h["score"]) for h in mine
    }


def test_annotate_motif_hits_respect_relative_threshold() -> None:
    hits = ok(reg.annotate("TTACGTTACG", options={"motifs": PWM, "min_relative_score": 0.95}))[
        "sequences"
    ][0]["hits"]
    assert [(h["start"], h["strand"], h["matched"]) for h in hits] == [
        (2, "+", "ACG"),
        (3, "-", "ACG"),
        (7, "+", "ACG"),
    ]
    loose = ok(
        reg.annotate(
            "TTACGTTACG", options={"motifs": PWM, "min_relative_score": 0.0, "strand": "forward"}
        )
    )["sequences"][0]
    assert loose["hit_count"] == 8  # every one of the 10 - 3 + 1 windows


def test_pwm_information_content_by_hand() -> None:
    result = ok(reg.calculate(PWM, options={"entity": "pwm_information_content"}))["motifs"][0]
    per_column = 0.7 * math.log2(0.7 / 0.25) + 3 * 0.1 * math.log2(0.1 / 0.25)
    assert result["per_position"] == [pytest.approx(per_column)] * 3
    assert result["total_bits"] == pytest.approx(3 * per_column)
    conserved = {"A": [1, 0], "C": [0, 0.5], "G": [0, 0.5], "T": [0, 0]}
    # column 1: one base with p = 1 -> 2 bits; column 2: two bases at 0.5 -> 1 bit (floor of 1e-6 is negligible)
    ic = ok(reg.calculate(conserved, options={"entity": "pwm_information_content"}))["motifs"][0][
        "per_position"
    ]
    assert ic == [pytest.approx(2.0, abs=1e-4), pytest.approx(1.0, abs=1e-4)]


def test_counts_matrix_with_pseudocount_and_jaspar_text() -> None:
    counts = {"A": [8, 0], "C": [1, 0], "G": [1, 10], "T": [0, 0]}
    result = ok(
        reg.score("AG", options={"motifs": counts, "pseudocount": 1.0, "strand": "forward"})
    )["sequences"][0]["motifs"][0]["best"]
    # (count + 1) / (10 + 4): A col0 = 9/14, G col1 = 11/14
    assert result["score"] == pytest.approx(math.log2(9 / 14 / 0.25) + math.log2(11 / 14 / 0.25))
    jaspar = ">MA0000.1 test\nA [ 8 0 ]\nC [ 1 0 ]\nG [ 1 10 ]\nT [ 0 0 ]\n"
    again = ok(
        reg.score("AG", options={"motifs": jaspar, "pseudocount": 1.0, "strand": "forward"})
    )["sequences"][0]["motifs"][0]["best"]
    assert again["score"] == pytest.approx(result["score"])


def test_regulatory_effect_delta_by_hand() -> None:
    delta = ok(
        reg.calculate(
            {"id": "v1", "sequence": "TTACGTT", "position": 3, "ref": "C", "alt": "A"},
            options={"motifs": PWM},
        )
    )["variants"][0]
    motif = delta["motifs"][0]
    # ACG -> AAG: the middle column now sees A (0.1): loses log2(0.7/0.1) bits
    assert motif["ref_score"] == pytest.approx(3 * HIGH) and motif["alt_score"] == pytest.approx(
        2 * HIGH + LOW
    )
    assert motif["delta"] == pytest.approx(LOW - HIGH) and delta["max_abs_delta"] == pytest.approx(
        math.log2(7)
    )
    with pytest.raises(InputFormatError):
        reg.calculate(
            {"id": "v", "sequence": "TTACGTT", "position": 3, "ref": "T", "alt": "A"},
            options={"motifs": PWM},
        )


def test_prioritize_ranks_by_absolute_delta() -> None:
    variants = [
        {"id": "neutral", "ref_sequence": "ACG", "alt_sequence": "ACG"},
        {"id": "strong", "ref_sequence": "ACG", "alt_sequence": "ACT"},
        {"id": "unscorable", "ref_sequence": "AC", "alt_sequence": "AG"},
    ]
    ranked = ok(reg.prioritize(variants, options={"motifs": PWM}))["ranked"]
    assert [(r["id"], r["rank"]) for r in ranked] == [
        ("strong", 1),
        ("neutral", 2),
        ("unscorable", 3),
    ]
    assert (
        ranked[0]["max_abs_delta"] == pytest.approx(math.log2(7))
        and ranked[1]["max_abs_delta"] == 0.0
    )
    assert ranked[2]["max_abs_delta"] is None
    assert len(ok(reg.prioritize(variants, options={"motifs": PWM, "top_n": 1}))["ranked"]) == 1


def test_in_silico_mutagenesis_by_hand() -> None:
    result = ok(reg.perform_in_silico_mutagenesis("AAACGAA", options={"motifs": PWM}))["motifs"][0]
    assert result["reference_score"] == pytest.approx(3 * HIGH)
    positions = result["positions"]
    assert len(positions) == 7 and positions[2]["ref"] == "A"
    # every substitution of the motif's first base (pos 2) turns a 0.7 column into 0.1: a loss of log2(7) bits
    for base in "CGT":
        assert positions[2]["effects"][base] == pytest.approx(LOW - HIGH)
    assert positions[2]["effects"].keys() == {"C", "G", "T"}
    # the last base lies only in windows that never reach the ACG optimum, so mutating it cannot change the best score
    assert positions[6]["effects"]["T"] == pytest.approx(0.0)
    # recomputation: brute-force best score of every mutated sequence
    from dnakit.regulation import _impl_local as impl

    pwm = impl.parse_motifs({"motifs": PWM})[0]
    for pos in range(7):
        for base, effect in positions[pos]["effects"].items():
            mutated = "AAACGAA"[:pos] + base + "AAACGAA"[pos + 1 :]
            assert effect == pytest.approx(
                impl.best_score(mutated, pwm) - impl.best_score("AAACGAA", pwm)
            )


def test_motif_comparison_orientation_and_offset() -> None:
    rc_pwm = {
        "A": list(reversed(PWM["T"])),
        "C": list(reversed(PWM["G"])),
        "G": list(reversed(PWM["C"])),
        "T": list(reversed(PWM["A"])),
    }
    pairs = ok(reg.compare({"m": PWM, "same": PWM, "rc": rc_pwm}))["pairs"]
    by_pair = {(p["a"], p["b"]): p for p in pairs}
    assert (
        by_pair[("m", "same")]["similarity"] == pytest.approx(1.0)
        and by_pair[("m", "same")]["offset"] == 0
    )
    assert (
        by_pair[("m", "rc")]["similarity"] == pytest.approx(1.0)
        and by_pair[("m", "rc")]["orientation"] == "reverse_complement"
    )
    shifted = {
        "A": [0.25] + PWM["A"],
        "C": [0.25] + PWM["C"],
        "G": [0.25] + PWM["G"],
        "T": [0.25] + PWM["T"],
    }
    sh = ok(reg.compare([PWM, shifted]))["pairs"][0]
    assert sh["offset"] == -1 and sh["similarity"] == pytest.approx(
        1.0
    )  # shifted motif's real columns start one column later
    with pytest.raises(InputFormatError):
        reg.compare(PWM)


def test_enrichment_hypergeometric_hand_and_scipy() -> None:
    # N = 6 sequences, K = 3 with the motif, all 3 foreground sequences have it: P = 1 / C(6,3) = 1/20
    counts = {
        "foreground_hits": 3,
        "foreground_total": 3,
        "background_hits": 0,
        "background_total": 3,
    }
    row = ok(reg.enrich(counts))["motifs"][0]
    assert row["p_value"] == pytest.approx(0.05) and row["fold_enrichment"] is None
    big = {
        "foreground_hits": 8,
        "foreground_total": 10,
        "background_hits": 20,
        "background_total": 100,
    }
    hyper = ok(reg.enrich(big))["motifs"][0]
    assert hyper["p_value"] == pytest.approx(scipy_stats.hypergeom.sf(7, 110, 28, 10), rel=1e-9)
    assert hyper["fold_enrichment"] == pytest.approx(4.0)
    binom = ok(reg.enrich(big, options={"test": "binomial"}))["motifs"][0]
    assert binom["p_value"] == pytest.approx(scipy_stats.binom.sf(7, 10, 0.2), rel=1e-9)
    # P(X >= 3 | n = 3, p = 0.5) = 1/8
    three = ok(
        reg.enrich(
            {
                "foreground_hits": 3,
                "foreground_total": 3,
                "background_hits": 5,
                "background_total": 10,
            },
            options={"test": "binomial"},
        )
    )
    assert three["motifs"][0]["p_value"] == pytest.approx(0.125)


def test_enrichment_from_sequences_counts_sequences_with_hits() -> None:
    foreground = ["TTACGTT", "GACGAA", "CCACGCC", "TTTTTT"]
    background = ["TTTTTT", "CCCCCC", "GGGGGG", "TATATA", "ACGTTT"]
    result = ok(
        reg.enrich(
            {"foreground": foreground, "background": background},
            options={"motifs": {"acg": PWM}, "min_relative_score": 1.0, "strand": "forward"},
        )
    )["motifs"][0]
    assert (
        result["foreground_hits"],
        result["foreground_total"],
        result["background_hits"],
        result["background_total"],
    ) == (3, 4, 1, 5)
    assert result["p_value"] == pytest.approx(scipy_stats.hypergeom.sf(2, 9, 4, 4))
    assert result["fold_enrichment"] == pytest.approx((3 / 4) / (1 / 5))
    with pytest.raises(InputFormatError):
        reg.enrich(
            {
                "foreground_hits": 5,
                "foreground_total": 3,
                "background_hits": 0,
                "background_total": 3,
            }
        )


def test_multiple_motifs_get_bh_adjusted() -> None:
    other = {"A": [0.1, 0.1, 0.1], "C": [0.1, 0.1, 0.1], "G": [0.1, 0.1, 0.1], "T": [0.7, 0.7, 0.7]}
    rows = ok(
        reg.enrich(
            {"foreground": ["TTT", "TTT"], "background": ["ACG", "CCC", "GGG"]},
            options={
                "motifs": {"acg": PWM, "ttt": other},
                "min_relative_score": 1.0,
                "strand": "forward",
            },
        )
    )["motifs"]
    assert {r["motif"] for r in rows} == {"acg", "ttt"}
    ttt = next(r for r in rows if r["motif"] == "ttt")
    assert ttt["foreground_hits"] == 2 and ttt["p_value"] == pytest.approx(
        scipy_stats.hypergeom.sf(1, 5, 2, 2)
    )
    assert all(r["q_value"] >= r["p_value"] - 1e-15 for r in rows)


def test_region_to_gene_link_and_annotate_target_gene() -> None:
    linked = ok(reg.link("chr1\t1000\t1200\te1", options={"genes": GENES}))["links"][0]
    assert linked["genes"][0]["gene"] == "G2" and linked["genes"][0]["distance"] == 399
    annotated = ok(
        reg.annotate("chr1\t1000\t1200\te1", options={"entity": "target_gene", "genes": GENES})
    )["links"][0]
    assert annotated == linked
    both = ok(
        reg.link(
            "chr1\t1000\t1200", options={"genes": GENES, "method": "window", "max_distance": 4000}
        )
    )["links"][0]["genes"]
    assert [g["gene"] for g in both] == ["G2", "G1"]
    assert both[0]["direction"] == "downstream" and both[1]["direction"] == "upstream"


def test_classify_promoter_proximal_distal() -> None:
    regions = "chr1\t4500\t4600\nchr1\t6000\t6100\nchr1\t100000\t100100\nchr1\t1450\t1600\nchr1\t1600\t1700"
    rows = ok(reg.classify(regions, options={"genes": GENES}))["regions"]
    # G1 (+, TSS 5000) promoter window [4000, 5100); G2 (-, TSS 1499) window [1400, 2500)
    assert [r["class"] for r in rows] == ["promoter", "proximal", "distal", "promoter", "promoter"]
    assert rows[0]["nearest_gene"] == "G1" and rows[0]["distance_to_tss"] == 450
    narrow = ok(
        reg.classify(
            "chr1\t1700\t1800",
            options={"genes": GENES, "promoter_upstream": 100, "promoter_downstream": 50},
        )
    )["regions"][0]
    assert narrow["class"] == "proximal"  # minus-strand window is [1450, 1600)
    assert (
        ok(reg.classify("chrZ\t1\t10", options={"genes": GENES}))["regions"][0]["class"] == "distal"
    )


def test_motif_validation_report() -> None:
    result = ok(
        reg.validate(
            {
                "good": PWM,
                "cnt": {"A": [5, 0], "C": [0, 5], "G": [0, 0], "T": [0, 0]},
                "negative": {"A": [-1, 0], "C": [1, 1], "G": [0, 0], "T": [0, 0]},
                "zero": {"A": [0], "C": [0], "G": [0], "T": [0]},
            }
        )
    )
    rows = {r["motif"]: r for r in result["motifs"]}
    assert (
        rows["good"]["valid"]
        and rows["good"]["matrix_type"] == "probabilities"
        and rows["good"]["length"] == 3
    )
    assert rows["cnt"]["valid"] and rows["cnt"]["matrix_type"] == "counts"
    assert not rows["negative"]["valid"] and not rows["zero"]["valid"]
    assert result["all_valid"] is False
    with pytest.raises(InputFormatError):
        reg.score("ACGT", options={"motifs": {"A": [1], "C": [1]}})


# ============================================================================ genome3d
TAD = [[10 if (i < 4) == (j < 4) else 1 for j in range(8)] for i in range(8)]


def test_insulation_score_and_boundary_by_hand() -> None:
    # two 4-bin domains (contacts 10 inside, 1 across), w=2: S(i)=mean M[i-2:i, i:i+2]
    # S = [-, -, 10, 5.5, 1, 5.5, 10, -]; mean(S) = 6.4; score(4) = log2(1 / 6.4)
    result = ok(g3d.calculate(TAD, options={"entity": "insulation", "window": 2}))
    assert result["mean_contact"] == [None, None, 10.0, 5.5, 1.0, 5.5, 10.0, None]
    assert result["score"][4] == pytest.approx(math.log2(1 / 6.4))
    assert result["score"][2] == pytest.approx(math.log2(10 / 6.4))
    annotated = ok(g3d.annotate(TAD, options={"window": 2, "boundary_radius": 2, "bin_size": 1000}))
    assert annotated["boundary_bins"] == [4] and annotated["boundary_bp"] == [4000]
    two = [[10 if (i // 4) == (j // 4) else 1 for j in range(12)] for i in range(12)]
    assert ok(g3d.annotate(two, options={"window": 2, "boundary_radius": 2}))["boundary_bins"] == [
        4,
        8,
    ]
    assert ok(g3d.annotate(two, options={"window": 2, "boundary_radius": 2}))["tads"] == [
        {"start_bin": 0, "end_bin": 4, "edge": True},
        {"start_bin": 4, "end_bin": 8, "edge": False},
        {"start_bin": 8, "end_bin": 12, "edge": True},
    ]


def test_distance_decay_by_hand() -> None:
    curve = ok(g3d.calculate(TAD, options={"bin_size": 5000}))["curve"]
    assert curve[0] == {"distance_bins": 0, "n_pairs": 8, "mean_contact": 10.0, "distance_bp": 0}
    assert (
        curve[1]["mean_contact"] == pytest.approx((6 * 10 + 1) / 7) and curve[1]["n_pairs"] == 7
    )  # only pair (3,4) crosses the domains
    assert curve[7]["mean_contact"] == 1.0 and curve[7]["distance_bp"] == 35000


def test_distance_decay_scaling_slope_recovers_power_law() -> None:
    n = 40
    matrix = [[100.0 / max(abs(i - j), 1) ** 1.5 for j in range(n)] for i in range(n)]
    scaling = ok(g3d.calculate(matrix, options={"s_min": 2, "s_max": 20}))["scaling"]
    assert scaling["slope"] == pytest.approx(-1.5, abs=1e-9)


def test_ice_normalisation_recovers_planted_bias_and_has_unit_row_sums() -> None:
    # circulant "true" contact map: every row sums to 9, so ICE must recover it from the biased matrix
    truth = np.array(
        [[4.0, 2.0, 1.0, 2.0], [2.0, 4.0, 2.0, 1.0], [1.0, 2.0, 4.0, 2.0], [2.0, 1.0, 2.0, 4.0]]
    )
    bias = np.array([0.5, 1.0, 2.0, 3.0])
    raw = truth * np.outer(bias, bias)
    out = ok(g3d.normalize(raw.tolist(), options={"tol": 1e-12, "max_iter": 1000}))
    balanced = np.array(out["matrix"])
    assert out["converged"] is True and np.allclose(balanced, balanced.T)
    assert np.allclose(balanced.sum(axis=1), 1.0, atol=1e-9)
    assert np.allclose(
        balanced, truth / truth.sum(axis=1)[0], atol=1e-9
    )  # raw / (b_i b_j), scaled to unit rows
    assert np.allclose(raw / np.outer(out["bias"], out["bias"]), balanced, atol=1e-9)
    scaled = ok(
        g3d.normalize(
            raw.tolist(), options={"tol": 1e-12, "max_iter": 1000, "target_row_sum": 10.0}
        )
    )
    assert np.allclose(np.array(scaled["matrix"]).sum(axis=1), 10.0, atol=1e-8)


def test_ice_masks_empty_bins_and_general_matrix_row_sums() -> None:
    rng = np.random.default_rng(0)
    a = rng.integers(1, 30, size=(6, 6)).astype(float)
    sym = a + a.T
    sym[2, :] = 0
    sym[:, 2] = 0
    out = ok(g3d.normalize(sym.tolist(), options={"tol": 1e-12, "max_iter": 5000}))
    m = np.array(out["matrix"])
    valid = [0, 1, 3, 4, 5]
    assert out["masked_bins"] == [2] and out["bias"][2] is None
    assert np.allclose(m, m.T) and np.allclose(m[valid].sum(axis=1), 1.0, atol=1e-9)
    assert np.all(m[2] == 0)
    diag_free = ok(
        g3d.normalize(
            (sym + np.diag(np.full(6, 100.0))).tolist(),
            options={"ignore_diagonals": 1, "tol": 1e-12, "max_iter": 5000},
        )
    )
    assert np.allclose(np.array(diag_free["matrix"])[valid].sum(axis=1), 1.0, atol=1e-9)
    assert np.allclose(np.diag(np.array(diag_free["matrix"]))[valid], 0.0)


def test_coverage_normalisation_by_hand() -> None:
    # row sums 3,3,6 (mean 4): weights .75,.75,1.5 -> M_ij / (w_i w_j)
    matrix = [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 6.0]]
    out = ok(g3d.normalize(matrix, options={"method": "coverage"}))
    assert out["matrix"][0][1] == pytest.approx(2 / (0.75 * 0.75)) and out["matrix"][2][
        2
    ] == pytest.approx(6 / (1.5 * 1.5))
    assert out["bias"] == [pytest.approx(0.75), pytest.approx(0.75), pytest.approx(1.5)]


def test_aggregate_conserves_contacts_and_handles_remainder() -> None:
    result = ok(g3d.aggregate(TAD, options={"factor": 4, "bin_size": 1000}))
    assert result["matrix"] == [[160.0, 16.0], [16.0, 160.0]] and result["bin_size"] == 4000
    odd = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]
    coarse = ok(g3d.aggregate(odd, options={"factor": 2}))["matrix"]
    assert coarse == [[1 + 2 + 2 + 4, 3 + 5], [3 + 5, 6]]
    assert sum(map(sum, coarse)) == sum(map(sum, odd))
    with pytest.raises(ConfigurationError):
        g3d.aggregate(TAD)


def test_read_dense_sparse_labels_and_binary_rejection(tmp_path) -> None:
    dense = ok(g3d.read("# c\nbin\tA\tB\nA\t1\t2\nB\t2\t3\n"))
    assert (
        dense["matrix"] == [[1.0, 2.0], [2.0, 3.0]]
        and dense["total_contacts"] == 6.0
        and dense["nonzero_entries"] == 4
    )
    sparse = ok(g3d.read("0\t0\t5\n0\t1\t2\n1\t1\t3\n", options={"format": "sparse"}))
    assert sparse["matrix"] == [[5.0, 2.0], [2.0, 3.0]]
    padded = ok(g3d.read("0 2 4", options={"format": "sparse", "n_bins": 4}))["matrix"]
    assert padded[2][0] == 4.0 and padded[0][2] == 4.0 and len(padded) == 4
    file = tmp_path / "m.txt"
    file.write_text("1,2\n2,3\n")
    assert ok(g3d.read(str(file)))["matrix"] == [[1.0, 2.0], [2.0, 3.0]]
    for bad in ("sample.hic", tmp_path / "x.mcool"):
        with pytest.raises(ConfigurationError):
            g3d.read(bad)


def test_export_write_roundtrip(tmp_path) -> None:
    matrix = [[5, 0, 1], [0, 2, 0], [1, 0, 7]]
    sparse = ok(g3d.export(matrix, options={"out_format": "sparse"}))
    assert sparse["text"] == "0\t0\t5\n0\t2\t1\n1\t1\t2\n2\t2\t7\n"
    assert ok(g3d.read(sparse["text"], options={"format": "sparse"}))["matrix"] == [
        [5.0, 0.0, 1.0],
        [0.0, 2.0, 0.0],
        [1.0, 0.0, 7.0],
    ]
    path = tmp_path / "out.tsv"
    written = ok(g3d.write(matrix, options={"output": str(path)}))
    assert (
        written["bytes_written"] == path.stat().st_size
        and path.read_text() == "5\t0\t1\n0\t2\t0\n1\t0\t7\n"
    )
    assert ok(g3d.read(path))["matrix"] == [[5.0, 0.0, 1.0], [0.0, 2.0, 0.0], [1.0, 0.0, 7.0]]
    with pytest.raises(ConfigurationError):
        g3d.write(matrix, options={"output": str(path)})
    with pytest.raises(ConfigurationError):
        g3d.write(matrix)


def test_compare_matrices_correlation() -> None:
    assert ok(g3d.compare({"a": TAD, "b": TAD}))["pearson"] == pytest.approx(1.0)
    scaled = [[2 * x for x in row] for row in TAD]
    both = ok(g3d.compare([TAD, scaled]))
    assert both["pearson"] == pytest.approx(1.0) and both["spearman"] == pytest.approx(1.0)
    assert both["mean_absolute_difference"] == pytest.approx(
        np.mean([x for i, row in enumerate(TAD) for j, x in enumerate(row) if j > i])
    )
    flipped = [[1 if (i < 4) == (j < 4) else 10 for j in range(8)] for i in range(8)]
    assert ok(g3d.compare([TAD, flipped]))["pearson"] == pytest.approx(-1.0)
    rng = np.random.default_rng(1)
    a = rng.random((9, 9))
    a = a + a.T
    b = rng.random((9, 9))
    b = b + b.T
    got = ok(g3d.compare([a.tolist(), b.tolist()]))
    iu = np.triu_indices(9, k=1)
    assert got["pearson"] == pytest.approx(scipy_stats.pearsonr(a[iu], b[iu]).statistic)
    assert got["spearman"] == pytest.approx(scipy_stats.spearmanr(a[iu], b[iu]).statistic)
    with pytest.raises(InputFormatError):
        g3d.compare([TAD, [[1, 1], [1, 1]]])


def test_compartments_recover_planted_labels_and_follow_reference() -> None:
    labels = [0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1]
    matrix = [
        [(1 + 4 * (labels[i] == labels[j])) * 10 / (1 + abs(i - j)) for j in range(12)]
        for i in range(12)
    ]
    called = ok(g3d.annotate(matrix, options={"entity": "compartment"}))
    letters = called["labels"]
    assert letters in (
        ["A" if x == 0 else "B" for x in labels],
        ["B" if x == 0 else "A" for x in labels],
    )
    gc = [0.6 if x == 1 else 0.4 for x in labels]
    oriented = ok(g3d.annotate(matrix, options={"entity": "compartment", "reference": gc}))
    assert (
        oriented["labels"] == ["A" if x == 1 else "B" for x in labels]
        and oriented["oriented_by_reference"] is True
    )
    eig = ok(g3d.calculate(matrix, options={"entity": "compartment", "reference": gc}))[
        "eigenvector"
    ]
    assert all((v >= 0) == (x == 1) for v, x in zip(eig, labels, strict=False))


def test_compartments_agree_with_independent_numpy_eigendecomposition() -> None:
    rng = np.random.default_rng(3)
    base = rng.random((10, 10)) * 5 + 1
    sym = base + base.T
    out = ok(g3d.calculate(sym.tolist(), options={"entity": "compartment"}))["eigenvector"]
    n = 10
    expected = np.array([np.mean(np.diagonal(sym, k)) for k in range(n)])
    oe = sym / expected[np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])]
    _values, vectors = np.linalg.eigh(np.corrcoef(oe))
    vec = vectors[:, -1]
    vec = vec if vec[np.argmax(np.abs(vec))] > 0 else -vec
    assert np.allclose(out, vec, atol=1e-9)


def test_heatmap_plot_cells() -> None:
    svg = ok(g3d.plot(TAD))["svg"]
    assert len(re.findall(r'class="cell"', svg)) == 64 and 'data-bins="8"' in svg
    cell = re.search(r'data-i="0" data-j="7" data-value="1"', svg)
    assert cell is not None
    big = [[1.0] * 30 for _ in range(30)]
    small = ok(g3d.plot(big, options={"max_cells": 10}))["svg"]
    assert len(re.findall(r'class="cell"', small)) == 100 and 'data-factor="3"' in small
    assert 'data-value="9"' in small  # 3x3 blocks of ones sum to 9


GENOME3D_BAD_INPUT = [
    [[1, 2], [3, 1]],  # asymmetric
    [[1, -2], [-2, 1]],  # negative contacts
    [[1, 2, 3], [2, 1, 1]],  # not square
    "a\tb\nc\td",  # non-numeric
    [],
]


@pytest.mark.parametrize("matrix", GENOME3D_BAD_INPUT)
def test_genome3d_rejects_invalid_matrices(matrix) -> None:
    for function in (g3d.read, g3d.normalize, g3d.plot, g3d.aggregate):
        with pytest.raises(DNAKitError):
            function(matrix, options={"factor": 2})


def test_genome3d_unknown_entity_and_method() -> None:
    with pytest.raises(ConfigurationError):
        g3d.calculate(TAD, options={"entity": "tad"})
    with pytest.raises(ConfigurationError):
        g3d.normalize(TAD, options={"method": "kr"})


# ============================================================================ honest leftovers
NOT_IMPLEMENTED = [
    (micro.assemble, "microbiome"),
    (micro.bin, "microbiome"),
    (micro.build_microbial_gene_catalog, "microbiome"),
    (micro.build_pathogen_transmission_tree, "microbiome"),
    (micro.calculate_digital_ddh, "microbiome"),
    (micro.classify, "microbiome"),
    (micro.infer_viral_quasispecies, "microbiome"),
    (micro.predict, "microbiome"),
    (micro.track, "microbiome"),
    (crispr.predict_editing_outcomes, "design_crispr"),
    (epi.align, "epigenome"),
    (epi.predict, "epigenome"),
    (reg.predict, "regulation"),
]


@pytest.mark.parametrize(
    ("function", "module"), NOT_IMPLEMENTED, ids=lambda v: getattr(v, "__name__", str(v))
)
def test_functions_without_an_honest_implementation_report_not_implemented(
    function, module: str
) -> None:
    envelope = function("ACGTACGT", options={"anything": 1})
    assert envelope["status"] == "not_implemented"
    assert envelope["module"] == module and envelope["operation"] == function.__name__
    assert envelope["required_dependency"] and envelope["message"]
    assert "result" not in envelope


# ============================================================================ envelope / signature contracts
def test_envelope_shape_of_implemented_functions() -> None:
    envelope = micro.qc(TABLE)
    assert (
        envelope["status"] == "ok"
        and envelope["module"] == "microbiome"
        and envelope["operation"] == "qc"
    )
    assert (
        envelope["provenance"]["implementation"] == "native"
        and "result" in envelope
        and "input_summary" in envelope
    )


@pytest.mark.parametrize("module", [micro, crispr, epi, reg, g3d], ids=lambda m: m.__name__)
def test_public_signatures_are_unchanged(module) -> None:
    import inspect

    for name in module.__all__:
        function = getattr(module, name)
        if getattr(function, "__module__", "") != module.__name__:
            continue
        parameters = list(inspect.signature(function).parameters)
        assert parameters[:2] == ["input", "options"], name
        assert inspect.signature(function).parameters["options"].default is None
        doc = (function.__doc__ or "").strip()
        if doc:  # the untouched pre-existing "specific" adapters have no docstring
            assert doc.splitlines()[0].endswith((".", ")")), name
