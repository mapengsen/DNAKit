"""Tests for dnakit.phylogeny.iqtree_p11_cli_functions (IQ-TREE 3 facades); real runs skip when iqtree3 is missing."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import iqtree_p11_cli_functions as m

requires_iqtree = pytest.mark.skipif(
    shutil.which("iqtree3") is None, reason="iqtree3 is not installed"
)
TREE = "((A:0.1,B:0.12):0.05,((C:0.2,D:0.1):0.07,(E:0.15,F:0.11):0.06):0.04,(G:0.3,H:0.2):0.05);\n"


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    if shutil.which("iqtree3") is None:
        pytest.skip("iqtree3 is not installed")
    root = tmp_path_factory.mktemp("iq")
    tree = root / "tree.nwk"
    tree.write_text(TREE)
    sim = m.simulate_alignment(
        root / "sim", tree=tree, model="HKY", length=400, seed=7, verbose=False
    )
    trees = root / "trees.nwk"
    trees.write_text(TREE + "((A,C),((B,D),(E,F)),(G,H));\n")
    return {"tree": tree, "aln": Path(sim["alignment_files"][0]), "trees": trees, "root": root}


def test_missing_program_and_validation(tmp_path: Path) -> None:
    aln = tmp_path / "a.fa"
    aln.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(BackendUnavailableError):
        m.select_substitution_model(
            aln, tmp_path / "o", iqtree_path=str(tmp_path / "none"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.select_substitution_model(aln, tmp_path / "o2", search="BAD", verbose=False)
    with pytest.raises(ConfigurationError):
        m.compute_branch_support(aln, tmp_path / "o3", methods=[], verbose=False)


@requires_iqtree
def test_simulation_and_site_statistics(data: dict[str, Path]) -> None:
    result = m.alignment_site_statistics(data["aln"], data["root"] / "stats", verbose=False)
    assert result["sites"] == 400 and result["version"] == "3.1.4"
    assert sum(result["counts"].values()) == 400


@requires_iqtree
def test_evaluate_tree_likelihood(data: dict[str, Path]) -> None:
    fixed = m.evaluate_tree_likelihood(
        data["aln"], data["tree"], data["root"] / "e1", model="HKY", verbose=False
    )
    shown = m.evaluate_tree_likelihood(
        data["aln"],
        data["tree"],
        data["root"] / "e2",
        model="HKY",
        optimise_branch_lengths=False,
        verbose=False,
    )
    assert fixed["log_likelihood"] >= shown["log_likelihood"]


@requires_iqtree
def test_model_selection_and_support(data: dict[str, Path]) -> None:
    sel = m.select_substitution_model(
        data["aln"], data["root"] / "mf", search="MF", seed=1, verbose=False
    )
    assert sel["best_model"] and sel["model_table"][0]["model"]
    sup = m.compute_branch_support(
        data["aln"],
        data["root"] / "sup",
        methods=["ufboot", "alrt", "abayes"],
        model="HKY",
        seed=1,
        verbose=False,
    )
    assert Path(sup["outputs"]["treefile"]).is_file() and sup["consensus_tree"]
    with pytest.raises(ConfigurationError):
        m.compute_branch_support(
            data["aln"], data["root"] / "bad", methods=["ufboot"], replicates=10, verbose=False
        )


@requires_iqtree
def test_topology_tests_mapping_concordance_ancestral(data: dict[str, Path]) -> None:
    tests = m.compare_tree_topologies(
        data["aln"],
        data["trees"],
        data["root"] / "tt",
        model="HKY",
        replicates=1000,
        seed=1,
        verbose=False,
    )
    assert len(tests["tests"]) == 2 and tests["tests"][0]["p_au"]["passed"] is True
    lmap = m.likelihood_mapping(
        data["aln"], data["root"] / "lm", quartets=200, model="HKY", seed=1, verbose=False
    )
    assert lmap["quartets"] == 200 and "fully_resolved" in lmap["resolution"]
    scf = m.compute_site_concordance_factors(
        data["aln"], data["tree"], data["root"] / "scf", model="HKY", seed=1, verbose=False
    )
    assert scf["branches"] and 0 <= scf["branches"][0]["scf"] <= 100
    anc = m.reconstruct_ancestral_sequences(
        data["aln"], data["tree"], data["root"] / "anc", model="HKY", verbose=False
    )
    assert anc["sites"] == 400 and all(len(s) == 400 for s in anc["sequences"].values())
    sym = m.run_symmetry_tests(data["aln"], data["root"] / "sym", seed=1, verbose=False)
    assert sym["tests"][0]["sym_non"] + sym["tests"][0]["sym_sig"] == 28


@requires_iqtree
def test_pathogen_dating_pomo_and_partitions(data: dict[str, Path]) -> None:
    path = m.infer_pathogen_tree(
        data["aln"], data["root"] / "pat", model="HKY", force_cmaple=True, seed=1, verbose=False
    )
    assert Path(path["outputs"]["treefile"]).is_file()
    dates = data["root"] / "dates.txt"
    dates.write_text("".join(f"{n} {2000 + i}\n" for i, n in enumerate("ABCDEFGH")))
    dated = m.date_phylogeny(
        data["aln"],
        dates,
        data["root"] / "date",
        tree=data["tree"],
        model="HKY",
        seed=1,
        verbose=False,
    )
    assert Path(dated["outputs"]["timetree.nwk"]).is_file()
    counts = data["root"] / "counts.txt"
    counts.write_text(
        "COUNTSFILE NPOP 3 NSITES 4\nCHROM POS s1 s2 s3\nc 1 3,0,0,0 3,0,0,0 0,0,3,0\n"
        "c 2 0,3,0,0 0,2,1,0 0,3,0,0\nc 3 3,0,0,0 2,1,0,0 0,0,0,3\nc 4 0,0,3,0 0,0,3,0 0,0,3,0\n"
    )
    pomo = m.infer_polymorphism_aware_tree(
        counts, data["root"] / "pomo", population_size=9, seed=1, verbose=False
    )
    assert pomo["pomo_model"] == "HKY+P+N9"
    parts = data["root"] / "part.nex"
    parts.write_text("#nexus\nbegin sets;\n charset p1 = 1-200;\n charset p2 = 201-400;\nend;\n")
    merged = m.select_substitution_model(
        data["aln"],
        data["root"] / "merge",
        search="TESTONLY",
        partitions=parts,
        merge=True,
        seed=1,
        verbose=False,
    )
    assert "best_scheme" in merged
