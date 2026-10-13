"""Contract tests for the tskit and dadi wrappers in dnakit.population (group G13).

Each wrapped upstream function is checked against a hand calculation or a direct call to the
upstream function. Inputs are written inside tmp_path, and no network is used.
"""

from __future__ import annotations

import io
import math
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.population import dadi_functions as dadi_api
from dnakit.population import tskit_functions as tskit_api

np = pytest.importorskip("numpy")

# Four samples (0-3) under two cherries: ((0,1),(2,3)); sites at 10, 20 and 30 on 0..100.
TABLES: dict[str, Any] = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}


def _star_tables() -> dict[str, Any]:
    """Four samples hanging from one root at time 1: a single polytomy."""

    return {
        "sequence_length": 100,
        "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}],
        "edges": [{"left": 0, "right": 100, "parent": 4, "child": child} for child in range(4)],
        "sites": [],
        "mutations": [],
    }


def _direct_ts(tskit: Any) -> Any:
    """The TABLES fixture built with tskit calls only, for direct comparisons."""

    tables = tskit.TableCollection(sequence_length=100)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    for time in (1, 2, 3):
        tables.nodes.add_row(flags=0, time=time)
    for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]:
        tables.edges.add_row(left=0, right=100, parent=parent, child=child)
    for position in (10, 20, 30):
        tables.sites.add_row(position=position, ancestral_state="A")
    for site, node in [(0, 0), (1, 4), (2, 5)]:
        tables.mutations.add_row(site=site, node=node, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


@pytest.fixture
def tskit_mod() -> Any:
    return pytest.importorskip("tskit")


@pytest.fixture
def dadi_mod() -> Any:
    pytest.importorskip("numpy")
    return pytest.importorskip("dadi")


def _same(left: Any, right: Any) -> None:
    """Compare nested lists of floats and None, with a relative tolerance on floats."""

    if isinstance(left, list):
        assert isinstance(right, list) and len(left) == len(right)
        for a, b in zip(left, right, strict=True):
            _same(a, b)
    elif left is None or right is None:
        assert left is None and right is None
    else:
        assert left == pytest.approx(right, rel=1e-9, abs=1e-12)


def _masked_to_lists(np: Any, array: Any) -> Any:
    data = np.ma.getdata(array).tolist()
    mask = np.ma.getmaskarray(array).tolist()

    def walk(value: Any, flag: Any) -> Any:
        if isinstance(value, list):
            return [walk(v, f) for v, f in zip(value, flag, strict=True)]
        return None if flag else float(value)

    return walk(data, mask)


# ----------------------------------------------------------------------------------------------
# tskit
# ----------------------------------------------------------------------------------------------


def test_export_sample_sequences_fasta_matches_hand_haplotypes(tmp_path, tskit_mod):
    target = tmp_path / "samples.fa"
    result = tskit_api.export_sample_sequences(TABLES, target, wrap_width=60)
    assert result["num_samples"] == 4 and result["num_sites"] == 3
    records = {}
    name = None
    for line in target.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        else:
            records[name] += line
    # Site positions are 0-based indices into the 100-character sequence.
    assert [len(records[f"n{i}"]) for i in range(4)] == [100] * 4
    expected = {
        "n0": ("G", "G", "A"),
        "n1": ("A", "G", "A"),
        "n2": ("A", "A", "G"),
        "n3": ("A", "A", "G"),
    }
    for sample, bases in expected.items():
        assert (records[sample][10], records[sample][20], records[sample][30]) == bases
    direct = tmp_path / "direct.fa"
    tskit_mod.write_fasta(
        _direct_ts(tskit_mod),
        str(direct),
        wrap_width=60,
        reference_sequence=None,
        missing_data_character="N",
    )
    assert target.read_text() == direct.read_text()


def test_export_sample_sequences_nexus_matches_direct_call(tmp_path, tskit_mod):
    target = tmp_path / "samples.nex"
    result = tskit_api.export_sample_sequences(TABLES, target, file_format="nexus")
    assert result["file_format"] == "nexus"
    direct = tmp_path / "direct.nex"
    tskit_mod.write_nexus(
        _direct_ts(tskit_mod),
        str(direct),
        precision=None,
        include_trees=False,
        include_alignments=None,
        reference_sequence=None,
        missing_data_character="N",
        node_labels=None,
    )
    assert target.read_text() == direct.read_text()
    assert "MISSING=N" in target.read_text()


def test_export_ms_format_matches_hand_genotypes(tmp_path, tskit_mod):
    target = tmp_path / "samples.ms"
    result = tskit_api.export_ms_format(TABLES, target, precision=4)
    assert result["num_segregating_sites"] == 3
    lines = target.read_text().splitlines()
    assert lines[4] == "segsites: 3"
    # Derived-allele rows: n0 carries sites 0 and 1, n1 carries site 1, n2 and n3 carry site 2.
    assert lines[6:10] == ["110", "010", "001", "001"]
    buffer = io.StringIO()
    tskit_mod.write_ms(_direct_ts(tskit_mod), buffer, precision=4, num_replicates=1)
    assert target.read_text() == buffer.getvalue()


def test_complement_intervals_matches_hand_calculation():
    result = tskit_api.complement_intervals([(10, 20), (30, 40)], end=100)
    assert result["intervals"] == [[0.0, 10.0], [20.0, 30.0], [40.0, 100.0]]
    assert result["complement_length"] == pytest.approx(80.0)
    assert result["covered_length"] == pytest.approx(20.0)
    assert tskit_api.complement_intervals([], end=50)["intervals"] == [[0.0, 50.0]]
    with pytest.raises(ConfigurationError):
        tskit_api.complement_intervals([(10, 40), (30, 50)], end=100)
    with pytest.raises(ConfigurationError):
        tskit_api.complement_intervals([(10, 120)], end=100)


def test_read_text_table_matches_direct_parser(tmp_path, tskit_mod):
    tree_dir = tmp_path / "text"
    tskit_api.write_text_tree_sequence(TABLES, tree_dir)
    nodes = tskit_api.read_text_table(tree_dir / "nodes.txt", table="nodes")
    assert nodes["num_rows"] == 7
    assert [row["time"] for row in nodes["rows"]] == [0, 0, 0, 0, 1, 2, 3]
    assert [row["flags"] for row in nodes["rows"][:4]] == [1, 1, 1, 1]
    direct = tskit_mod.parse_nodes(io.StringIO((tree_dir / "nodes.txt").read_text()), strict=True)
    assert direct.num_rows == nodes["num_rows"]
    sites = tskit_api.read_text_table(tree_dir / "sites.txt", table="sites")
    assert [row["position"] for row in sites["rows"]] == [10.0, 20.0, 30.0]
    assert {row["ancestral_state"] for row in sites["rows"]} == {"A"}
    with pytest.raises(ConfigurationError):
        tskit_api.read_text_table(tree_dir / "sites.txt", table="genomes")


def test_text_tree_sequence_round_trip_matches_load_text(tmp_path, tskit_mod):
    tree_dir = tmp_path / "text"
    written = tskit_api.write_text_tree_sequence(TABLES, tree_dir)
    assert written["files"] == ["nodes.txt", "edges.txt", "sites.txt", "mutations.txt"]
    read = tskit_api.read_text_tree_sequence(
        tree_dir / "nodes.txt",
        tree_dir / "edges.txt",
        sites=tree_dir / "sites.txt",
        mutations=tree_dir / "mutations.txt",
    )
    assert read["num_samples"] == 4 and read["num_trees"] == 1
    assert read["num_mutations"] == 3 and read["sequence_length"] == pytest.approx(100.0)
    with (tree_dir / "nodes.txt").open() as nodes, (tree_dir / "edges.txt").open() as edges:
        direct = tskit_mod.load_text(nodes=nodes, edges=edges, sequence_length=100)
    assert direct.num_nodes == read["num_nodes"] and direct.num_edges == read["num_edges"]
    with pytest.raises(ConfigurationError):
        tskit_api.write_text_tree_sequence(TABLES, tree_dir)


def test_count_sample_topologies_matches_direct_counter(tskit_mod):
    result = tskit_api.count_sample_topologies(TABLES, {"A": [0, 1], "B": [2, 3]})
    assert result["num_trees"] == 1 and result["sample_sets"] == {"A": 2, "B": 2}
    counter = next(
        iter(
            tskit_mod.combinatorics.treeseq_count_topologies(
                _direct_ts(tskit_mod), [[0, 1], [2, 3]]
            )
        )
    )
    direct = {
        "+".join("AB"[position] for position in key): sum(counter.topologies[key].values())
        for key in counter.topologies
    }
    observed = {
        label: sum(entry["count"] for entry in entries)
        for label, entries in result["trees"][0]["counts"].items()
    }
    assert observed == direct
    assert set(result["trees"][0]["counts"]) == {"A", "A+B", "B"}
    with pytest.raises(ConfigurationError):
        tskit_api.count_sample_topologies(TABLES, None)
    with pytest.raises(ConfigurationError):
        tskit_api.count_sample_topologies(TABLES, {"A": [0, 1], "B": [2, 6]})


def test_resolve_polytomies_turns_star_into_binary_tree(tskit_mod):
    result = tskit_api.resolve_polytomies(_star_tables(), random_seed=7)
    # A binary tree on 4 leaves has 3 internal nodes; the star had 1, so 2 nodes are added.
    assert result["added_nodes"] == 2
    assert result["newick"].count("(") == 3
    assert all(f"n{leaf}" in result["newick"] for leaf in range(4))
    again = tskit_api.resolve_polytomies(_star_tables(), random_seed=7)
    assert again["newick"] == result["newick"]
    binary = tskit_api.resolve_polytomies(TABLES)
    assert binary["added_nodes"] == 0


def test_tree_to_newick_matches_hand_branch_lengths(tskit_mod):
    result = tskit_api.tree_to_newick(TABLES, taxon_labels=["a", "b", "c", "d"])
    # Leaves at time 0 under nodes at times 1 and 2, root at 3: branch lengths 1, 1, 2, 2 and 2, 1.
    assert result["newick"] == (
        "((a:1.000000,b:1.000000):2.000000,(c:2.000000,d:2.000000):1.000000);"
    )
    default = tskit_api.tree_to_newick(TABLES, precision=0, include_branch_lengths=False)
    assert default["newick"] == "((n0,n1),(n2,n3));"
    with pytest.raises(ConfigurationError):
        tskit_api.tree_to_newick(TABLES, taxon_labels=["a", "a", "c", "d"])
    with pytest.raises(ConfigurationError):
        tskit_api.tree_to_newick(TABLES, taxon_labels=["a(", "b", "c", "d"])


def test_random_tree_newick_matches_direct_generator(tskit_mod):
    labelled = {f"n{node}": node for node in range(5)}
    result = tskit_api.random_tree_newick(5, seed=3)
    tree = tskit_mod.combinatorics.generate_random_binary(
        5, span=1, branch_length=1.0, random_seed=3, record_provenance=False
    )
    labels = {node: f"n{node}" for node in tree.tree_sequence.samples()}
    direct = tskit_mod.build_newick(
        tree, root=tree.root, precision=6, node_labels=labels, include_branch_lengths=True
    )
    assert result["newick"] == direct
    assert all(f"{label}:" in result["newick"] for label in labelled)
    star = tskit_api.random_tree_newick(4, shape="star")
    star_tree = tskit_mod.combinatorics.generate_star(
        4, span=1, branch_length=1.0, record_provenance=False
    )
    star_labels = {node: f"n{node}" for node in star_tree.tree_sequence.samples()}
    assert star["newick"] == tskit_mod.build_newick(
        star_tree,
        root=star_tree.root,
        precision=6,
        node_labels=star_labels,
        include_branch_lengths=True,
    )
    balanced = tskit_api.random_tree_newick(4, shape="balanced", arity=2)
    assert balanced["newick"].count("(") == 3
    with pytest.raises(ConfigurationError):
        tskit_api.random_tree_newick(4, shape="star", seed=1)


def test_map_allele_indices_matches_hand_mapping(tskit_mod):
    result = tskit_api.map_allele_indices(["A", "G", "T", "X"], ["A", "G"])
    assert result["indices"] == [0, 1, 2, 3]
    assert result["unmapped_alleles"] == ["T", "X"]
    other = tskit_api.map_allele_indices(["G", "A"], ["A", "C", "G"])
    assert other["indices"] == [2, 0]
    assert other["unmapped_alleles"] == []
    with pytest.raises(ConfigurationError):
        tskit_api.map_allele_indices(["A"], ["A", "A"])


def test_read_plink_fam_counts_and_direct_parse(tmp_path, tskit_mod):
    fam = tmp_path / "study.fam"
    fam.write_text("f1 i1 0 0 1\nf1 i2 0 0 2\nf2 i3 0 0 0\n")
    result = tskit_api.read_plink_fam(fam)
    assert result["num_individuals"] == 3
    assert (result["num_males"], result["num_females"], result["num_unknown_sex"]) == (1, 1, 1)
    assert [item["iid"] for item in result["individuals"]] == ["i1", "i2", "i3"]
    direct = tskit_mod.parse_fam(io.StringIO(fam.read_text()))
    assert [direct[index].metadata["sex"] for index in range(3)] == [1, 2, 0]


# ----------------------------------------------------------------------------------------------
# dadi
# ----------------------------------------------------------------------------------------------

SNP_TEXT = (
    "Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid\n"
    "ACA TCA A 2 1 G 0 1 snp1\n"
    "GCT GCT C 1 2 T 1 0 snp2\n"
    "TTG TAG T 2 2 G 0 0 snp3\n"
    "ACG ACG A 1 1 C 1 1 snp4\n"
)

VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
chr1\t100\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|0\t0|1\t1|1\t0|0
chr1\t200\t.\tC\tT\t.\tPASS\tAA=T\tGT\t0|1\t1|1\t0|0\t0|0
chr1\t300\t.\tG\tA\t.\tPASS\tAA=G\tGT\t0|0\t0|0\t0|1\t1|1
"""


def _chunked_snp_file(path, positions) -> None:
    """Haploid SNP rows with IDs chr1_<position>, two populations with two chromosomes each."""

    rows = ["Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid"]
    patterns = [(2, 1, 0, 1), (1, 0, 1, 2), (0, 2, 2, 0), (1, 1, 1, 1), (2, 0, 0, 2), (0, 1, 2, 1)]
    for position, (a, b, c, d) in zip(positions, patterns, strict=False):
        rows.append(f"ACG ACG A {a} {b} G {c} {d} chr1_{position}")
    path.write_text("\n".join(rows) + "\n")


def test_read_snp_data_matches_hand_polarization(tmp_path, dadi_mod):
    snp = tmp_path / "snps.txt"
    snp.write_text(SNP_TEXT)
    result = dadi_api.read_snp_data(snp)
    assert result["num_snps"] == 4 and result["num_biallelic"] == 4
    assert result["populations"] == ["pop1", "pop2"]
    # snp2 (outgroup C, alleles C and T) and snp4 (outgroup C, alleles A and C) are polarized.
    assert result["num_polarized"] == 2


def test_snp_configuration_counts_matches_upstream(tmp_path, dadi_mod):
    snp = tmp_path / "snps.txt"
    snp.write_text(SNP_TEXT)
    result = dadi_api.snp_configuration_counts(snp, pop_ids=["pop1", "pop2"])
    direct = dadi_mod.Misc.count_data_dict(dadi_mod.Misc.make_data_dict(str(snp)), ["pop1", "pop2"])
    assert result["num_configurations"] == len(direct)
    assert result["num_snps_counted"] == sum(direct.values()) == 4
    for row in result["configurations"]:
        key = (tuple(row["successful_calls"]), tuple(row["derived_calls"]), row["polarized"])
        assert direct[key] == row["count"]
    with pytest.raises(ConfigurationError):
        dadi_api.snp_configuration_counts(snp, pop_ids=["pop1", "pop9"])


def test_read_vcf_data_builtin_and_cyvcf2_agree(tmp_path, dadi_mod):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(VCF_TEXT)
    popinfo = tmp_path / "pop.txt"
    popinfo.write_text("s1 popA\ns2 popA\ns3 popB\ns4 popB\n")
    builtin = dadi_api.read_vcf_data(vcf, popinfo)
    assert builtin["num_snps"] == 3 and builtin["populations"] == ["popA", "popB"]
    direct = dadi_mod.Misc.make_data_dict_vcf(str(vcf), str(popinfo), filter=True)
    assert len(direct) == builtin["num_snps"]
    pytest.importorskip("cyvcf2")
    cyvcf = dadi_api.read_vcf_data(vcf, popinfo, parser="cyvcf2")
    assert cyvcf["num_snps"] == builtin["num_snps"]
    with pytest.raises(ConfigurationError):
        dadi_api.read_vcf_data(vcf, popinfo, parser="pysam")


@pytest.mark.filterwarnings("ignore:Mean of empty slice")
@pytest.mark.filterwarnings("ignore:Degrees of freedom")
def test_bootstrap_snp_spectra_matches_upstream_bootstraps(tmp_path, dadi_mod, monkeypatch):
    snp = tmp_path / "chunks.txt"
    _chunked_snp_file(snp, [100, 200, 1100, 1200, 2100, 2200])
    result = dadi_api.bootstrap_snp_spectra(
        snp, pop_ids=["pop1", "pop2"], projections=[2, 2], chunk_size=1000, num_boot=20, seed=5
    )
    assert result["num_chunks"] == 3 and result["num_boot"] == 20
    data = dadi_mod.Misc.make_data_dict(str(snp))
    fragments = dadi_mod.Misc.fragment_data_dict(data, 1000)
    import random

    state = random.getstate()
    random.seed(5)
    spectra = dadi_mod.Misc.bootstraps_from_dd_chunks(
        fragments, 20, ["pop1", "pop2"], [2, 2], mask_corners=True, polarized=True
    )
    random.setstate(state)
    mean = np.nanmean(
        np.array([np.ma.filled(np.ma.asarray(s, dtype=float), np.nan) for s in spectra]), axis=0
    )
    expected = np.ma.masked_invalid(mean)
    _same(result["mean_spectrum"], _masked_to_lists(np, expected))
    again = dadi_api.bootstrap_snp_spectra(
        snp, pop_ids=["pop1", "pop2"], projections=[2, 2], chunk_size=1000, num_boot=20, seed=5
    )
    assert again["mean_spectrum"] == result["mean_spectrum"]


def test_pool_sfs_populations_matches_hand_sums(dadi_mod):
    result = dadi_api.pool_sfs_populations([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # Pooled index i + j sums the cells on each anti-diagonal; the corners are masked.
    assert result["sample_sizes"] == [4]
    assert result["spectrum"] == [None, 6.0, 15.0, 14.0, None]
    three = [float(index) for index in range(27)]
    shaped = [[[three[i * 9 + j * 3 + k] for k in range(3)] for j in range(3)] for i in range(3)]
    pooled = dadi_api.pool_sfs_populations(shaped, populations=(0, 2))
    direct = dadi_mod.Misc.combine_pops(dadi_mod.Spectrum(shaped), idx=[0, 2])
    _same(pooled["spectrum"], _masked_to_lists(np, dadi_mod.Spectrum(direct)))
    with pytest.raises(ConfigurationError):
        dadi_api.pool_sfs_populations(shaped, populations=(0, 0))


def test_sfs_file_round_trip_matches_upstream_reader(tmp_path, dadi_mod):
    target = tmp_path / "spectrum.sfs"
    written = dadi_api.write_sfs_file([[1.5, 2.0], [3.0, 4.0]], target, comment_lines=["hello"])
    assert written["shape"] == [2, 2] and written["num_values"] == 4
    lines = target.read_text().splitlines()
    assert lines[0] == "# hello"
    assert lines[1] == "2 2 "
    read = dadi_api.read_sfs_file(target)
    assert read["spectrum"] == [[1.5, 2.0], [3.0, 4.0]]
    assert read["comments"] == ["hello"]
    with target.open() as handle:
        direct, comments = dadi_mod.Numerics.array_from_file(handle, return_comments=True)
    assert direct.tolist() == read["spectrum"] and comments == read["comments"]
    with pytest.raises(ConfigurationError):
        dadi_api.write_sfs_file([[1.0]], target)


def test_apply_ancestral_misidentification_matches_hand_formula(dadi_mod):
    result = dadi_api.apply_ancestral_misidentification([0, 1, 2, 3, 0], 0.1)
    # (1 - p) * fs + p * reversed(fs) with p = 0.1: 0.9 * 1 + 0.1 * 3 = 1.2, and so on.
    assert result["spectrum"][0] is None and result["spectrum"][4] is None
    assert result["spectrum"][1:4] == pytest.approx([1.2, 2.0, 2.8])
    with pytest.raises(ConfigurationError):
        dadi_api.apply_ancestral_misidentification([0, 1, 2, 3, 0], 1.5)


def test_expected_sfs_1d_matches_upstream_model(dadi_mod):
    result = dadi_api.expected_sfs_1d("two_epoch", [2.0, 0.5], 10, grid_points=50)
    direct = dadi_mod.Demographics1D.two_epoch((2.0, 0.5), (10,), 50)
    _same(result["spectrum"], _masked_to_lists(np, direct))
    assert result["spectrum"][0] is None and result["spectrum"][10] is None
    with pytest.raises(ConfigurationError):
        dadi_api.expected_sfs_1d("two_epoch", [2.0], 10, grid_points=50)
    with pytest.raises(ConfigurationError):
        dadi_api.expected_sfs_1d("no_such_model", [2.0, 0.5], 10)


def test_expected_sfs_2d_matches_upstream_model(dadi_mod):
    result = dadi_api.expected_sfs_2d("split_mig", [1.0, 1.0, 0.5, 1.0], [4, 4], grid_points=30)
    direct = dadi_mod.Demographics2D.split_mig((1.0, 1.0, 0.5, 1.0), (4, 4), 30)
    _same(result["spectrum"], _masked_to_lists(np, direct))
    with pytest.raises(ConfigurationError):
        dadi_api.expected_sfs_2d("split_mig", [1.0, 1.0, 0.5, 1.0], [4], grid_points=30)


def test_expected_sfs_3d_matches_upstream_model(dadi_mod):
    params = [1.0, 1.0, 0.5, 1.0, 0.5, 0.5]
    result = dadi_api.expected_sfs_3d("split_nomig", params, [2, 2, 2], grid_points=20)
    direct = dadi_mod.Demographics3D.split_nomig(tuple(params), (2, 2, 2), 20)
    _same(result["spectrum"], _masked_to_lists(np, direct))


def test_expected_sfs_selection_matches_upstream_model(dadi_mod):
    params = [1.0, 1.0, 0.5, 1.0, 1.0, 1.0]
    result = dadi_api.expected_sfs_selection("split_mig_sel", params, [4, 4], grid_points=20)
    from dadi.DFE import DemogSelModels

    direct = DemogSelModels.split_mig_sel(tuple(params), (4, 4), 20)
    _same(result["spectrum"], _masked_to_lists(np, direct))


def test_fit_demographic_model_matches_upstream_optimizer(dadi_mod):
    truth = dadi_mod.Demographics1D.two_epoch((2.0, 0.5), (10,), 50)
    observed = [0.0 if value is None else value for value in _masked_to_lists(np, truth)]
    data = [0.0] + [value * 100.0 for value in observed[1:-1]] + [0.0]
    result = dadi_api.fit_demographic_model(
        data,
        model="two_epoch",
        p0=[1.0, 0.3],
        lower_bound=[0.01, 0.001],
        upper_bound=[10.0, 5.0],
        grid_points=50,
        maxiter=3,
    )
    observed_spectrum = dadi_mod.Spectrum(np.array(data))
    direct = dadi_mod.Inference.optimize_log(
        [1.0, 0.3],
        observed_spectrum,
        dadi_mod.Demographics1D.two_epoch,
        50,
        lower_bound=[0.01, 0.001],
        upper_bound=[10.0, 5.0],
        verbose=0,
        multinom=True,
        maxiter=3,
    )
    assert result["params"] == pytest.approx(list(direct), rel=1e-9)
    model = dadi_mod.Demographics1D.two_epoch(tuple(result["params"]), (10,), 50)
    assert result["log_likelihood"] == pytest.approx(
        float(dadi_mod.Inference.ll_multinom(model, observed_spectrum)), rel=1e-9
    )
    with pytest.raises(ConfigurationError):
        dadi_api.fit_demographic_model(data, model="split_mig", p0=[1.0, 0.3], grid_points=50)


def test_sfs_residuals_match_hand_poisson_values():
    model = [1.0, 4.0, 9.0, 4.0, 1.0]
    data = [1.0, 2.0, 2.0, 2.0, 1.0]
    result = dadi_api.sfs_residuals(model, data, residual="poisson")
    # Interior residuals (model - data) / sqrt(model): 2/2, 7/3, 2/2.
    assert result["residuals"][1:4] == pytest.approx([1.0, 7.0 / 3.0, 1.0])
    assert result["residuals"][0] is None and result["residuals"][4] is None
    assert result["sum_of_squares"] == pytest.approx(2.0 + 49.0 / 9.0)
    assert result["num_entries"] == 3


def test_sfs_residuals_anscombe_matches_upstream(dadi_mod):
    model = [1.0, 4.0, 9.0, 4.0, 1.0]
    data = [1.0, 2.0, 2.0, 2.0, 1.0]
    result = dadi_api.sfs_residuals(model, data)
    direct = dadi_mod.Inference.Anscombe_Poisson_residual(
        dadi_mod.Spectrum(model), dadi_mod.Spectrum(data)
    )
    _same(result["residuals"], _masked_to_lists(np, direct))


def test_parameter_uncertainty_fim_matches_upstream(dadi_mod):
    truth = dadi_mod.Demographics1D.two_epoch((2.0, 0.5), (6,), 20)
    counts = (np.ma.filled(truth, 0.0) * 100.0).tolist()
    data = dadi_mod.Spectrum(np.array(counts))
    result = dadi_api.parameter_uncertainty_fim(
        counts, model="two_epoch", p0=[2.0, 0.5], grid_points=(20, 30, 40)
    )
    function = dadi_mod.Numerics.make_extrap_func(dadi_mod.Demographics1D.two_epoch)
    direct = dadi_mod.Godambe.FIM_uncert(
        function, [20, 30, 40], [2.0, 0.5], data, log=False, multinom=True, eps=0.01
    )
    assert result["uncertainties"] == pytest.approx([float(value) for value in direct], rel=1e-9)
    assert len(result["uncertainties"]) == 3 and all(
        math.isfinite(v) for v in result["uncertainties"]
    )


def test_annotate_snps_from_annovar_keeps_nonsynonymous(tmp_path, dadi_mod):
    snp = tmp_path / "chunks.txt"
    snp.write_text(
        "Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid\n"
        "ACA TCA A 2 1 G 0 1 chr1_100\n"
        "GCT GCT C 1 2 T 1 0 chr1_200\n"
    )
    annovar = tmp_path / "sample.exonic_variant_function"
    annovar.write_text(
        "1\tnonsynonymous SNV\tGENE1\tchr1 100\n2\tsynonymous SNV\tGENE2\tchr1 200\n"
    )
    result = dadi_api.annotate_snps_from_annovar(snp, annovar, variant_type="nonsynonymous")
    assert result["num_snps_total"] == 2 and result["num_snps_kept"] == 1
    assert result["kept_snp_ids"] == ["chr1_100"]


def test_trinucleotide_frequency_array_uses_dadi_order(dadi_mod):
    bases = "ACGT"
    trinucleotides = [a + b + c for a in bases for b in bases for c in bases]
    result = dadi_api.trinucleotide_frequency_array(
        {key: 0.0 for key in trinucleotides} | {"GTA": 1.0}
    )
    position = result["trinucleotides"].index("GTA")
    assert result["values"][position] == 1.0
    assert sum(result["values"]) == 1.0
    # dadi's own order is C, G, T, A for each base, so GTA is entry 16 + 8 + 3 = 27.
    assert position == 27
    direct = dadi_mod.Misc.tri_freq_dict_to_array(
        {key: 0.0 for key in trinucleotides} | {"GTA": 1.0}
    )
    assert list(direct) == result["values"]
    with pytest.raises(ConfigurationError):
        dadi_api.trinucleotide_frequency_array({"AAA": 1.0})


DEMES_TEXT = """description: split with migration
time_units: generations
demes:
  - name: ancestral
    epochs:
      - end_time: 2000
        start_size: 1000
  - name: A
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
  - name: B
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
migrations:
  - source: A
    dest: B
    start_time: 2000
    end_time: 0
    rate: 0.0001
"""


def test_expected_sfs_from_demes_matches_upstream(tmp_path, dadi_mod):
    demes = pytest.importorskip("demes")
    graph = tmp_path / "split.yaml"
    graph.write_text(DEMES_TEXT)
    result = dadi_api.expected_sfs_from_demes(graph, ["A", "B"], [4, 4], grid_points=40)
    direct = dadi_mod.Demes.SFS(demes.load(str(graph)), ["A", "B"], [4, 4], 40, theta=1.0)
    _same(result["spectrum"], _masked_to_lists(np, direct))
    assert result["spectrum"][0][0] is None and len(result["spectrum"]) == 5
    single = dadi_api.expected_sfs_from_demes(graph, ["A"], [4], grid_points=40)
    assert len(single["spectrum"]) == 5
    with pytest.raises(ConfigurationError):
        dadi_api.expected_sfs_from_demes(graph, ["Z"], [4])
