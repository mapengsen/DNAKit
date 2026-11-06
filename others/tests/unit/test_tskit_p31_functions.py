"""Tests for the tskit P31 functions (tree-sequence editing, statistics, exports)."""

from __future__ import annotations

import inspect

import pytest

pytest.importorskip("tskit")
pytest.importorskip("numpy")

from dnakit.exceptions import ConfigurationError
from dnakit.population import tskit_p31_functions as m

# 4 samples, two trees: [0,50) ((0,1)4,(2,3)5)6 ; [50,100) ((0,2)4',(1,3)5')6'
TABLES = {
    "sequence_length": 100,
    "nodes": [
        {"time": 0, "sample": True},
        {"time": 0, "sample": True},
        {"time": 0, "sample": True},
        {"time": 0, "sample": True},
        {"time": 1},
        {"time": 2},
        {"time": 4},
        {"time": 1},
        {"time": 3},
        {"time": 5},
    ],
    "edges": [
        {"left": 0, "right": 50, "parent": 4, "child": 0},
        {"left": 0, "right": 50, "parent": 4, "child": 1},
        {"left": 0, "right": 50, "parent": 5, "child": 2},
        {"left": 0, "right": 50, "parent": 5, "child": 3},
        {"left": 0, "right": 50, "parent": 6, "child": 4},
        {"left": 0, "right": 50, "parent": 6, "child": 5},
        {"left": 50, "right": 100, "parent": 7, "child": 0},
        {"left": 50, "right": 100, "parent": 7, "child": 2},
        {"left": 50, "right": 100, "parent": 8, "child": 1},
        {"left": 50, "right": 100, "parent": 8, "child": 3},
        {"left": 50, "right": 100, "parent": 9, "child": 7},
        {"left": 50, "right": 100, "parent": 9, "child": 8},
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "C"},
        {"position": 60, "ancestral_state": "G"},
        {"position": 80, "ancestral_state": "T"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "T"},
        {"site": 1, "node": 5, "derived_state": "G"},
        {"site": 2, "node": 7, "derived_state": "A"},
        {"site": 3, "node": 8, "derived_state": "C"},
    ],
}


@pytest.fixture(scope="module")
def ts():
    import tskit

    tskit_mod = tskit
    tc = tskit_mod.TableCollection(100)
    for row in TABLES["nodes"]:
        tc.nodes.add_row(flags=1 if row.get("sample") else 0, time=row["time"])
    for e in TABLES["edges"]:
        tc.edges.add_row(**e)
    for s in TABLES["sites"]:
        tc.sites.add_row(**s)
    for mu in TABLES["mutations"]:
        tc.mutations.add_row(**mu)
    tc.sort()
    return tc.tree_sequence()


def test_all_public_functions_are_keyword_only_with_verbose():
    names = [n for n in dir(m) if n.startswith("tskit_") and callable(getattr(m, n))]
    assert len(names) == 39
    for name in names:
        sig = inspect.signature(getattr(m, name))
        assert sig.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY


def test_listings_and_breakpoints(ts):
    assert m.tskit_list_sites(TABLES)["num_sites"] == 4
    assert m.tskit_list_sites(TABLES, tree_index=1)["num_sites"] == 2
    muts = m.tskit_list_mutations(TABLES)
    assert muts["num_mutations"] == 4 and muts["mutations"][0]["node"] == 4
    assert m.tskit_breakpoints(TABLES)["breakpoints"] == [0.0, 50.0, 100.0]


def test_genotypes_variants_alignments_haplotypes(ts):
    gm = m.tskit_genotype_matrix(TABLES)
    assert gm["matrix"] == ts.genotype_matrix().tolist()
    var = m.tskit_variants(TABLES)
    assert var["variants"][0]["states"] == ["T", "T", "A", "A"]
    assert var["variants"][0]["counts"] == {"T": 2, "A": 2}
    al = m.tskit_alignments(TABLES, reference_sequence="A" * 100)
    assert [len(a["sequence"]) for a in al["alignments"]] == [100] * 4
    assert al["alignments"][0]["sequence"][10] == "T"
    hp = m.tskit_haplotypes(TABLES)
    assert hp["haplotypes"][0]["haplotype"] == "TCAT"


def test_statistics_match_tskit(ts):
    sets = [[0, 1], [2, 3]]
    seg = m.tskit_segregating_sites(TABLES, sample_sets=sets, span_normalise=False)
    assert seg["values"] == ts.segregating_sites(sets, span_normalise=False).tolist()
    f4 = m.tskit_sample_set_statistic(
        TABLES, statistic="f4", sample_sets=[[0], [1], [2], [3]], indexes=[[0, 1, 2, 3]]
    )
    assert f4["values"] == ts.f4([[0], [1], [2], [3]], indexes=[(0, 1, 2, 3)]).tolist()
    y1 = m.tskit_sample_set_statistic(TABLES, statistic="Y1", sample_sets=[[0, 1, 2]])
    assert y1["values"] == ts.Y1([[0, 1, 2]]).tolist()
    y2 = m.tskit_sample_set_statistic(TABLES, statistic="Y2", sample_sets=sets, indexes=[[0, 1]])
    assert y2["values"] == ts.Y2(sets, indexes=[(0, 1)]).tolist()
    rel = m.tskit_genetic_relatedness(TABLES, sample_sets=[[0], [1], [2]], indexes=[[0, 1]])
    assert rel["values"] == ts.genetic_relatedness([[0], [1], [2]], indexes=[(0, 1)]).tolist()
    mat = m.tskit_genetic_relatedness(TABLES, sample_sets=[[0], [1]])
    assert mat["output"] == "matrix" and len(mat["values"]) == 2
    import numpy as np

    w = [[1.0], [2.0], [0.5], [3.0]]
    wa = np.array(w)
    for stat, ref in (("correlation", ts.trait_correlation), ("covariance", ts.trait_covariance)):
        assert (
            m.tskit_trait_statistic(TABLES, statistic=stat, weights=w)["values"] == ref(wa).tolist()
        )
    assert (
        m.tskit_trait_statistic(TABLES, statistic="linear_model", weights=w)["values"]
        == ts.trait_linear_model(wa).tolist()
    )
    assert (
        m.tskit_genetic_relatedness_vector(TABLES, weights=w)["values"]
        == ts.genetic_relatedness_vector(wa, mode="branch").tolist()
    )
    gs = m.tskit_general_stat(
        TABLES, weights=[[1.0], [-1.0], [2.0], [-2.0]], summary="square", polarised=True
    )
    assert len(gs["values"]) == 1
    sc = m.tskit_general_stat(
        TABLES, sample_sets=sets, summary="product", polarised=True, strict=False
    )
    assert sc["input"] == "sample_sets"
    md = m.tskit_mean_descendants(TABLES, sample_sets=sets)
    assert md["values"] == ts.mean_descendants(sets).tolist()
    gnn = m.tskit_genealogical_nearest_neighbours(TABLES, focal=[0], sample_sets=[[1], [2, 3]])
    assert gnn["values"] == ts.genealogical_nearest_neighbours([0], [[1], [2, 3]]).tolist()


def test_pair_coalescence_and_ibd(ts):
    import numpy as np

    sets = [[0, 1, 2, 3]]
    counts = m.tskit_pair_coalescence(TABLES, output="counts", sample_sets=sets)
    assert counts["values"] == ts.pair_coalescence_counts(sample_sets=sets).tolist()
    rates = m.tskit_pair_coalescence(TABLES, output="rates", time_windows=[0, 2], sample_sets=sets)
    assert (
        rates["values"]
        == ts.pair_coalescence_rates(np.array([0, 2, np.inf]), sample_sets=sets).tolist()
    )
    q = m.tskit_pair_coalescence(TABLES, output="quantiles", quantiles=[0.5, 1.0], sample_sets=sets)
    assert len(q["values"]) == 2
    with pytest.raises(ConfigurationError):
        m.tskit_pair_coalescence(TABLES, output="rates", time_windows="nodes")
    ibd = m.tskit_ibd_segments(TABLES)
    assert ibd["num_segments"] == len(
        list(ts.ibd_segments(store_segments=True, store_pairs=True).values())
        and [s for v in ts.ibd_segments(store_segments=True, store_pairs=True).values() for s in v]
    )


def test_ld_and_relatedness_inputs(ts):
    import tskit

    calc = tskit.LdCalculator(ts)
    assert m.tskit_ld_r2(TABLES)["values"] == calc.r2_matrix().tolist()
    assert m.tskit_ld_r2(TABLES, output="pair", site_a=0, site_b=1)["values"] == calc.r2(0, 1)
    assert (
        m.tskit_ld_r2(TABLES, output="array", focal_site=0)["values"] == calc.r2_array(0).tolist()
    )
    with pytest.raises(ConfigurationError):
        m.tskit_ld_r2(TABLES, output="pair")


def test_editing_functions(ts):
    out = m.tskit_select_intervals(TABLES, intervals=[[0, 50]], action="keep")
    assert out["result_summary"]["num_sites"] == 2
    gone = m.tskit_select_intervals(TABLES, intervals=[[0, 50]], action="delete")
    assert gone["result_summary"]["num_sites"] == 2
    assert m.tskit_delete_sites(TABLES, site_ids=[0, 1])["result_summary"]["num_sites"] == 2
    simp = m.tskit_simplify(TABLES, samples=[0, 1])
    assert simp["result_summary"]["num_samples"] == 2 and len(simp["node_map"]) == 10
    assert m.tskit_subset(TABLES, nodes=[0, 1, 4])["result_summary"]["num_nodes"] == 3
    dec = m.tskit_remove_older_genealogy(TABLES, time=3.5)
    assert all(dec["tables"]["nodes"][e["parent"]]["time"] <= 3.5 for e in dec["tables"]["edges"])
    old = m.tskit_remove_older_genealogy(TABLES, time=3.5, add_roots=False)
    assert all(old["tables"]["nodes"][e["parent"]]["time"] < 3.5 for e in old["tables"]["edges"])
    shifted = m.tskit_shift_coordinates(TABLES, value=10)
    assert shifted["tables"]["sequence_length"] == 110.0
    trimmed = m.tskit_trim_genome(
        m.tskit_shift_coordinates(TABLES, value=10)["tables"], side="left"
    )
    assert trimmed["tables"]["sequence_length"] == 100.0
    assert m.tskit_split_edges(TABLES, time=1.5)["result_summary"]["num_nodes"] > 10
    ext = m.tskit_extend_haplotypes(m.tskit_compute_mutation_times(TABLES)["tables"])
    assert ext["result_summary"]["num_samples"] == 4
    cat = m.tskit_concatenate(TABLES, [TABLES])
    assert cat["result_summary"]["sequence_length"] == 200.0
    uni = m.tskit_union(TABLES, TABLES, node_mapping=list(range(10)))
    assert uni["result_summary"]["num_nodes"] == 10
    la = m.tskit_link_ancestors(TABLES, samples=[0, 1], ancestors=[4])
    assert la["num_edges"] >= 1
    mp = m.tskit_compute_mutation_parents(TABLES)
    assert mp["mutation_parents"] == [-1, -1, -1, -1]
    mt = m.tskit_compute_mutation_times(TABLES, method="tables")
    assert all(t is not None for t in mt["mutation_times"])
    imp = m.tskit_compute_mutation_times(TABLES, method="min")
    assert imp["mutation_times"] == ts.impute_unknown_mutations_time().tolist()
    dup = m.tskit_deduplicate_sites(TABLES)
    assert dup["sites_removed"] == 0


def test_deduplicate_merges_duplicate_positions(tmp_path):
    import tskit

    tc = tskit.TableCollection(10)
    tc.nodes.add_row(flags=1, time=0)
    tc.nodes.add_row(flags=1, time=0)
    tc.nodes.add_row(time=1)
    tc.edges.add_row(0, 10, 2, 0)
    tc.edges.add_row(0, 10, 2, 1)
    tc.sites.add_row(5, "A")
    tc.sites.add_row(5, "A")
    tc.mutations.add_row(0, 0, "T")
    tc.mutations.add_row(1, 1, "G")
    tc.sort()
    path = tmp_path / "dup.trees"
    tc.dump(str(path))
    res = m.tskit_deduplicate_sites(str(path))
    assert res["sites_removed"] == 1 and res["result_summary"]["num_mutations"] == 2


def test_map_mutations_and_vcf_model(ts):
    res = m.tskit_map_mutations(TABLES, genotypes=[1, 1, 0, 0], tree_index=0)
    assert res["ancestral_state"] == "0" and res["num_mutations"] == 1
    assert res["mutations"][0]["node"] == 4
    vm = m.tskit_map_to_vcf_model(TABLES, ploidy=2)
    assert vm["individual_names"] == ["tsk_0", "tsk_1"]
    assert m.tskit_sample_nodes_by_ploidy(TABLES, ploidy=2)["groups"] == [[0, 1], [2, 3]]


def test_recombination_map(tmp_path):
    hap = tmp_path / "map.txt"
    hap.write_text(
        "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n"
        "chr1\t0\t1.0\t0.0\nchr1\t1000000\t2.0\t1.0\nchr1\t2000000\t0.0\t3.0\n"
    )
    rm = m.tskit_read_hapmap(hap)
    assert rm["sequence_length"] == 2000000.0 and len(rm["intervals"]) >= 2
    rate = m.tskit_recombination_map_lookup([500000, 1500000], hapmap_file=hap)
    assert rate["values"][0] == pytest.approx(1e-8) and rate["values"][1] == pytest.approx(2e-8)
    mass = m.tskit_recombination_map_lookup([1000000], quantity="cumulative_mass", hapmap_file=hap)
    assert mass["values"][0] == pytest.approx(1e-8 * 1e6)
    manual = m.tskit_recombination_map_lookup([5], map_positions=[0, 10, 20], map_rates=[0.1, 0.2])
    assert manual["values"] == [0.1]
    with pytest.raises(ConfigurationError):
        m.tskit_recombination_map_lookup([5])


def test_validation_errors():
    with pytest.raises(ConfigurationError):
        m.tskit_delete_sites(TABLES, site_ids=[9])
    with pytest.raises(ConfigurationError):
        m.tskit_select_intervals(TABLES, intervals=[[0, 1]], action="remove")
    with pytest.raises(ConfigurationError):
        m.tskit_sample_set_statistic(TABLES, statistic="D")
    with pytest.raises(ConfigurationError):
        m.tskit_general_stat(TABLES)
    with pytest.raises(ConfigurationError):
        m.tskit_simplify(TABLES, samples=[99])
