"""Hand-verified tests for the real assembly, haplotypes and phylogeny candidate functions."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

import dnakit.assembly as asm
import dnakit.haplotypes as hap
import dnakit.phylogeny as phylo
from dnakit.assembly import _impl_assembly
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.haplotypes import _impl_haplotypes


def result_of(envelope: object) -> dict:
    assert isinstance(envelope, dict)
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# --------------------------------------------------------------------------- assembly: pangenome


PAN = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}


def test_add_genome_to_pangenome() -> None:
    res = result_of(
        asm.add_genome_to_pangenome(
            {"A": ["g1", "g2"], "B": ["g2", "g3"]}, options={"genome": "C", "genes": ["g2", "g4"]}
        )
    )
    assert res["genomes"] == ["A", "B", "C"]
    assert res["genes"] == ["g1", "g2", "g3", "g4"]
    assert res["matrix"] == [[1, 0, 0], [1, 1, 1], [0, 1, 0], [0, 0, 1]]
    assert (res["new_gene_count"], res["pangenome_size_before"], res["pangenome_size_after"]) == (
        1,
        3,
        4,
    )


def test_build_pangenome_partition_counts() -> None:
    res = result_of(asm.build_pangenome(PAN, options={"cloud_threshold": 0.5}))
    # a in 3/3 (core); b in 2/3 (shell); c, d, e in 1/3 < 0.5 (cloud)
    assert (res["pangenome_size"], res["core_count"], res["shell_count"], res["cloud_count"]) == (
        5,
        1,
        1,
        3,
    )
    assert res["genes_per_genome"] == {"G1": 3, "G2": 3, "G3": 2}
    assert res["unique_genes_per_genome"] == {"G1": 1, "G2": 1, "G3": 1}


@pytest.mark.parametrize(
    ("threshold", "core"),
    [(1.0, ["a"]), (0.66, ["a", "b"]), (0.3, ["a", "b", "c", "d", "e"])],
)
def test_extract_core_pangenome(threshold: float, core: list[str]) -> None:
    res = result_of(asm.extract_core_pangenome(PAN, options={"core_threshold": threshold}))
    assert res["core_genes"] == core
    assert res["core_count"] == len(core)


def test_extract_accessory_pangenome() -> None:
    res = result_of(
        asm.extract_accessory_pangenome(
            PAN, options={"core_threshold": 0.66, "cloud_threshold": 0.5}
        )
    )
    assert res["accessory_genes"] == ["c", "d", "e"]
    assert res["cloud_genes"] == ["c", "d", "e"] and res["shell_genes"] == []


def test_pangenome_matrix_forms_agree() -> None:
    explicit = {
        "genomes": ["G1", "G2", "G3"],
        "genes": ["a", "b", "c", "d", "e"],
        "matrix": [[1, 1, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    }
    nested = {
        "a": {"G1": 1, "G2": 1, "G3": 1},
        "b": {"G1": 1, "G2": 1, "G3": 0},
        "c": {"G1": 1},
        "d": {"G2": 1},
        "e": {"G3": 1},
    }
    for form in (explicit, nested):
        assert result_of(asm.extract_core_pangenome(form, options={"core_threshold": 1.0}))[
            "core_genes"
        ] == ["a"]


def test_pangenome_openness_new_gene_curve_and_fit() -> None:
    sets = {
        "g1": list("abcdefgh"),
        "g2": list("abcd") + list("ijkl"),  # 4 new
        "g3": list("ab") + list("mn"),  # 2 new
        "g4": ["a", "o"],  # 1 new
        "g5": ["a", "p"],  # 1 new
    }
    res = result_of(asm.calculate_pangenome_openness(sets, options={"n_permutations": 0}))
    assert res["mean_new_genes_by_genome_number"] == [8, 4, 2, 1, 1]
    assert res["mean_pangenome_size_by_genome_number"] == [8, 12, 14, 15, 16]
    # OLS of ln(delta) on ln(N) for N = 2..5, delta = 4, 2, 1, 1 (computed by hand: slope -1.6325)
    assert res["alpha"] == pytest.approx(1.6325, abs=2e-3)
    assert res["kappa"] == pytest.approx(math.exp(0.52 + 1.6325 * 1.19685), rel=5e-3)
    assert res["open"] is False


def test_heaps_fit_recovers_exact_power_law() -> None:
    deltas = [50.0 * n**-0.5 for n in range(2, 8)]
    fit = _impl_assembly.fit_heaps(deltas)
    assert fit["alpha"] == pytest.approx(0.5)
    assert fit["kappa"] == pytest.approx(50.0)


def test_pangenome_openness_seeded_permutations_are_reproducible() -> None:
    sets = {
        f"g{i}": [f"core{j}" for j in range(5)] + [f"acc{i}_{j}" for j in range(4)]
        for i in range(6)
    }
    first = result_of(
        asm.calculate_pangenome_openness(sets, options={"n_permutations": 20, "seed": 3})
    )
    second = result_of(
        asm.calculate_pangenome_openness(sets, options={"n_permutations": 20, "seed": 3})
    )
    assert first == second
    # every genome adds exactly 4 private genes regardless of order: delta = 4 for N >= 2 -> alpha = 0 (open)
    assert first["alpha"] == pytest.approx(0.0, abs=1e-9)
    assert first["open"] is True


# --------------------------------------------------------------------------- assembly: statistics and k-mer models


def test_estimate_size_from_histogram() -> None:
    res = result_of(asm.estimate_size({1: 500, 10: 20000, 20: 90000}))
    # valley at coverage 2, peak at 20: (10*20000 + 20*90000) / 20 = 100000
    assert res["genome_size"] == pytest.approx(100000)
    assert (res["min_coverage"], res["peak_coverage"]) == (2, 20)


def test_estimate_size_from_jellyfish_text_and_explicit_options() -> None:
    text = "1 50\n2 3\n3 100\n4 400\n5 100\n6 5\n"
    # valley at 2 (3 <= 100), peak 4: (2*3 + 3*100 + 4*400 + 5*100 + 6*5) / 4 = 609
    assert result_of(asm.estimate_size(text))["genome_size"] == pytest.approx(609)
    explicit = result_of(asm.estimate_size(text, options={"min_coverage": 3, "peak_coverage": 4}))
    assert explicit["genome_size"] == pytest.approx((300 + 1600 + 500 + 30) / 4)


def test_estimate_size_counts_kmers_of_sequences() -> None:
    # "AAAAAAAA" has 6 AAA k-mers (canonical AAA); every k-mer instance has coverage 6 => one distinct k-mer at 6
    res = result_of(
        asm.estimate_size("AAAAAAAA", options={"k": 3, "min_coverage": 1, "peak_coverage": 6})
    )
    assert res["genome_size"] == pytest.approx(1.0)


def test_estimate_heterozygosity_recovers_simulated_rate() -> None:
    # k = 10, r = 0.01, L = 1e5: het k-mers 2*k*r*L = 20000 at half coverage, hom k-mers L*(1-k*r) = 90000
    res = result_of(asm.estimate_heterozygosity({1: 500, 10: 20000, 20: 90000}, options={"k": 10}))
    assert res["heterozygosity"] == pytest.approx(0.01)
    assert res["het_to_hom_ratio"] == pytest.approx(2 / 9)


def test_estimate_repeats_from_histogram() -> None:
    res = result_of(asm.estimate_repeats({1: 100, 10: 1000, 20: 50, 30: 10}))
    # T = 10*1000 + 20*50 + 30*10 = 11300, repeat part = 1300
    assert res["repeat_fraction"] == pytest.approx(1300 / 11300)
    assert res["copy_number_classes"] == {"2": 50, "3": 10}
    assert res["genome_size_estimate"] == pytest.approx(1130)


def test_calculate_assembly_completeness() -> None:
    by_length = result_of(
        asm.calculate_assembly_completeness([100, 50], options={"reference_length": 200})
    )
    assert by_length["length_ratio"] == pytest.approx(0.75)
    assert by_length["length_completeness"] == pytest.approx(0.75)
    capped = result_of(
        asm.calculate_assembly_completeness([300], options={"reference_length": 200})
    )
    assert capped["length_ratio"] == pytest.approx(1.5) and capped["length_completeness"] == 1.0
    # reference ACGTACGT (k=4): canonical k-mers {ACGT, CGTA, GTAC}; assembly ACGTA has ACGT, CGTA -> 2/3
    kmer = result_of(
        asm.calculate_assembly_completeness("ACGTA", options={"reference": "ACGTACGT", "k": 4})
    )
    assert kmer["kmer_completeness"] == pytest.approx(2 / 3)
    assert (kmer["reference_distinct_kmers"], kmer["shared_kmers"]) == (3, 2)


def test_calculate_busco_completeness() -> None:
    line = "C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614"
    res = result_of(asm.calculate_busco_completeness(f"# BUSCO\n\t{line}\n"))
    assert res["complete_percent"] == pytest.approx(100 * (1566 + 24) / 1614)
    assert (res["single"], res["duplicated"], res["fragmented"], res["missing"], res["n"]) == (
        1566,
        24,
        8,
        16,
        1614,
    )
    counts = result_of(
        asm.calculate_busco_completeness(
            {"single": 90, "duplicated": 5, "fragmented": 3, "missing": 2}
        )
    )
    assert counts["complete_percent"] == pytest.approx(95.0)
    assert counts["summary"] == "C:95.0%[S:90.0%,D:5.0%],F:3.0%,M:2.0%,n:100"


def test_check_input() -> None:
    good = result_of(asm.check_input(">c1\nACGT\n>c2\nGGNN\n"))
    assert good["valid"] is True and good["contig_count"] == 2 and good["total_length"] == 8
    bad = result_of(asm.check_input(">c1\nACGT\n>c1\nAXGT\n"))
    assert bad["valid"] is False
    assert bad["duplicate_names"] == ["c1"] and bad["invalid_symbols"] == {"c1": ["X"]}
    assert result_of(asm.check_input([]))["valid"] is False  # empty input is reported, not raised


def test_compare_assemblies() -> None:
    left = {"c1": "AAAA", "c2": "CCCCCC"}
    right = {"x": "TTTT", "y": "GGGGGGG"}
    res = result_of(asm.compare_assemblies([left, right]))
    assert res["left"]["total_length"] == 10 and res["right"]["total_length"] == 11
    assert res["left"]["n50"] == 6 and res["right"]["n50"] == 7
    assert res["n50_ratio_right_to_left"] == pytest.approx(7 / 6)
    assert res["difference_right_minus_left"]["total_length"] == 1
    # AAAA and TTTT are reverse complements: one shared contig sequence
    assert (res["shared_contig_sequences"], res["shared_bases"]) == (1, 4)
    assert res["left"]["gc_fraction"] == pytest.approx(0.6)  # 6 C of 10 bases
    via_option = result_of(asm.compare_assemblies(left, options={"other": right}))
    assert via_option["shared_contig_sequences"] == 1


def test_assembly_contiguity_stats_auN_and_gaps() -> None:
    res = result_of(asm.compare_assemblies([">a\nACGTNNNNAC\n>b\nNNAC\n", [10, 20]]))
    left = res["left"]
    assert left["gap_count"] == 2 and left["n_bases"] == 6  # runs NNNN and NN
    assert left["gc_fraction"] == pytest.approx(
        4 / 8
    )  # C,G,C (contig a) + C (contig b) of 8 ACGT bases
    assert left["aun"] == pytest.approx((100 + 16) / 14)
    assert res["right"]["aun"] == pytest.approx((100 + 400) / 30)


def test_identify_repeats_in_assembly() -> None:
    seq = "GATTACACCGATTACA"
    res = result_of(asm.identify_repeats_in_assembly({"c": seq}, options={"k": 5}))
    assert res["repeats"]["c"] == [
        {"start": 0, "end": 7, "length": 7},
        {"start": 9, "end": 16, "length": 7},
    ]
    assert res["repeat_bases"] == 14 and res["repeat_fraction"] == pytest.approx(14 / 16)


def test_identify_telomeres() -> None:
    contigs = {
        "left_tel": "TTAGGG" * 6 + "ACGTACGTACGTACGGATC",
        "right_tel": "GATCGATCGATCGATC" + "CCCTAA" * 5,
        "none": "ACGTACGTACGTACGTACGTACGT",
        "short": "TTAGGG" * 4 + "ACGTACGTACGT",
    }
    res = result_of(asm.identify_telomeres(contigs))
    assert res["contigs"]["left_tel"]["left"] == {"start": 0, "end": 36, "copies": 6}
    assert res["contigs"]["left_tel"]["right"] is None
    assert res["contigs"]["right_tel"]["right"] == {"start": 16, "end": 46, "copies": 5}
    assert res["contigs"]["none"]["telomeres_at_ends"] == 0
    assert res["contigs"]["short"]["telomeres_at_ends"] == 0  # 4 copies < min_repeats = 5
    assert res["contigs_with_telomere"] == 2
    relaxed = result_of(asm.identify_telomeres(contigs, options={"min_repeats": 4}))
    assert relaxed["contigs"]["short"]["left"]["copies"] == 4


# --------------------------------------------------------------------------- assembly: graphs and files


def test_build_graph_genome_unitigs() -> None:
    res = result_of(asm.build_graph_genome({"s1": "ACGTT", "s2": "ACGTA"}, options={"k": 3}))
    assert sorted(u["sequence"] for u in res["unitigs"]) == ["ACGT", "GTA", "GTT"]
    assert res["link_count"] == 2 and res["kmer_count"] == 4  # ACG CGT GTT GTA
    by_seq = {u["sequence"]: u for u in res["unitigs"]}
    assert by_seq["ACGT"]["sources"] == ["s1", "s2"] and by_seq["GTT"]["sources"] == ["s1"]
    assert "\tACGT\t" in res["gfa"] and "+\t2M" in res["gfa"]
    linear = result_of(asm.build_graph_genome("ACGTAC", options={"k": 3}))
    assert [u["sequence"] for u in linear["unitigs"]] == ["ACGTAC"] and linear["link_count"] == 0


GFA_A = "H\tVN:Z:1.0\nS\ta\tACGT\nS\tb\tAAAC\nL\ta\t+\tb\t+\t0M\n"
# q stores the reverse complement of b, so a+ -> b+ is the same link as p+ -> q-
GFA_B = "S\tp\tACGT\nS\tq\tGTTT\nS\tr\tCCCC\nL\tp\t+\tq\t-\t0M\nL\tq\t+\tr\t+\t0M\n"


def test_compare_assembly_graphs() -> None:
    res = result_of(asm.compare_assembly_graphs([GFA_A, GFA_B]))
    assert res["compared_by"] == "canonical segment sequence"
    assert (
        res["shared_segments"] == 2
        and res["right_only_segments"] == 1
        and res["left_only_segments"] == 0
    )
    assert res["segment_jaccard"] == pytest.approx(2 / 3)
    assert res["shared_links"] == 1 and res["link_jaccard"] == pytest.approx(1 / 2)
    assert res["identical"] is False
    same = result_of(asm.compare_assembly_graphs([GFA_A, GFA_A]))
    assert same["identical"] is True and same["link_jaccard"] == 1.0


def test_plot_graph_writes_svg(tmp_path: Path) -> None:
    target = tmp_path / "graph.svg"
    res = result_of(asm.plot_graph(GFA_B, options={"output": str(target)}))
    assert res["node_count"] == 3 and res["edge_count"] == 2
    assert res["svg"].count("<circle") == 3 and res["svg"].count("<line") == 2
    assert target.read_text(encoding="utf-8") == res["svg"]


def test_parse_output_fasta_and_gfa() -> None:
    fasta = ">NODE_1_length_8_cov_5.5\nACGTACGT\n>NODE_2_length_4_cov_2.0\nGGCC\n"
    res = result_of(asm.parse_output(fasta))
    assert res["format"] == "fasta" and res["header_convention"] == "spades"
    assert [r["length"] for r in res["records"]] == [8, 4]
    assert [r["coverage"] for r in res["records"]] == [5.5, 2.0]
    assert (
        res["records"][0]["gc_fraction"] == pytest.approx(0.5)
        and res["records"][1]["gc_fraction"] == 1.0
    )
    assert res["stats"]["n50"] == 8
    mega = result_of(asm.parse_output(">k141_0 flag=1 multi=2.5000 len=4\nACGT\n"))
    assert mega["header_convention"] == "megahit" and mega["records"][0]["multiplicity"] == 2.5
    gfa = result_of(
        asm.parse_output("S\ts1\tACGT\tDP:f:3.5\nS\ts2\t*\tLN:i:100\nL\ts1\t+\ts2\t+\t0M\n")
    )
    assert gfa["format"] == "gfa" and gfa["link_count"] == 1
    assert gfa["records"] == [
        {"name": "s1", "length": 4, "depth": 3.5},
        {"name": "s2", "length": 100, "depth": None},
    ]


def test_parse_log_counts_levels() -> None:
    log = "INFO started\nwarning: low coverage\nERROR failed to open file\nINFO done\n[debug] details\n"
    res = result_of(asm.parse_log(log))
    assert res["level_counts"] == {"CRITICAL": 0, "ERROR": 1, "WARNING": 1, "INFO": 2, "DEBUG": 1}
    assert res["has_errors"] is True and res["errors"][0]["line"] == 3
    assert res["line_count"] == 5


def test_remove_contamination() -> None:
    res = result_of(
        asm.remove_contamination(
            {"c1": "AAAA", "c2": "CCCCCC", "c3": "GG"}, options={"contaminants": ["c2", "zz"]}
        )
    )
    assert res["kept"] == ["c1", "c3"] and res["removed"] == ["c2"] and res["not_found"] == ["zz"]
    assert res["fasta"] == ">c1\nAAAA\n>c3\nGG\n"
    assert (res["stats_before"]["total_length"], res["stats_after"]["total_length"]) == (12, 6)


def test_scaffold_with_explicit_layout() -> None:
    contigs = {"c1": "AACC", "c2": "GGTT", "c3": "ACGT"}
    layout = [{"name": "S1", "contigs": ["c1", {"contig": "c2", "orientation": "-"}]}]
    res = result_of(asm.scaffold(contigs, options={"layout": layout, "gap": 5}))
    # c2 reverse-complemented: GGTT -> AACC
    assert res["scaffolds"] == [
        {"name": "S1", "sequence": "AACC" + "N" * 5 + "AACC", "length": 13, "contig_count": 2}
    ]
    assert [(row["type"], row["start"], row["end"]) for row in res["agp"]] == [
        ("contig", 0, 4),
        ("gap", 4, 9),
        ("contig", 9, 13),
    ]
    assert res["unplaced_contigs"] == ["c3"]


@pytest.mark.parametrize(
    ("tool", "settings", "argv"),
    [
        (
            "flye",
            {"reads": "r.fq", "genome_size": "5m", "threads": 4},
            [
                "flye",
                "--nano-hq",
                "r.fq",
                "--genome-size",
                "5m",
                "--out-dir",
                "out",
                "--threads",
                "4",
            ],
        ),
        (
            "hifiasm",
            {"reads": "r.fq", "threads": 2},
            ["hifiasm", "-o", "out/assembly", "-t", "2", "r.fq"],
        ),
        (
            "megahit",
            {"read1": "a.fq", "read2": "b.fq"},
            ["megahit", "-o", "out/megahit", "-t", "1", "-1", "a.fq", "-2", "b.fq"],
        ),
        (
            "unicycler",
            {"read1": "a.fq", "mode": "bold"},
            ["unicycler", "-1", "a.fq", "-o", "out/unicycler", "-m", "bold", "-t", "1"],
        ),
        (
            "shasta",
            {"reads": "r.fq"},
            [
                "shasta",
                "--input",
                "r.fq",
                "--config",
                "Nanopore-May2022",
                "--assemblyDirectory",
                "out",
            ],
        ),
    ],
)
def test_build_command(tool: str, settings: dict, argv: list[str]) -> None:
    res = result_of(asm.build_command({"tool": tool, "output_directory": "out", **settings}))
    assert res["argv"] == argv


def test_run_delegates_to_installed_assembler(tmp_path: Path) -> None:
    fake = tmp_path / "fake-shasta"
    fake.write_text(
        f"#!{sys.executable}\nfrom pathlib import Path\nimport sys\n"
        "d = Path(sys.argv[sys.argv.index('--assemblyDirectory') + 1])\n"
        "(d / 'Assembly.fasta').write_text('>c\\nACGT\\n')\n(d / 'Assembly.gfa').write_text('H\\tVN:Z:1.0\\n')\n",
        encoding="ascii",
    )
    fake.chmod(0o755)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n", encoding="ascii")
    out = asm.run(
        {
            "tool": "shasta",
            "reads": str(reads),
            "output_directory": str(tmp_path / "out"),
            "program": str(fake),
        }
    )
    assert out["status"] == "ok"
    assert (
        out["result"]["backend"] == "shasta"
        and "Assembly.fasta" in out["result"]["result"]["artifacts"]
    )


def test_run_without_program_is_not_implemented(tmp_path: Path) -> None:
    out = asm.run(
        {
            "tool": "flye",
            "reads": "r.fq",
            "genome_size": "5m",
            "output_directory": str(tmp_path),
            "program": "no-such-flye-binary",
        }
    )
    assert out["status"] == "not_implemented"
    assert "no-such-flye-binary" in out["required_dependency"]


# --------------------------------------------------------------------------- haplotypes


TRIO = {
    "child": ["0/1", "1/1", "0/0", "./."],
    "father": ["0/0", "0/1", "0/1", "0/0"],
    "mother": ["0/1", "0/1", "1/1", "0/0"],
}


def test_analyze_trio_inheritance() -> None:
    res = result_of(hap.analyze_trio_inheritance(TRIO))
    assert [s["status"] for s in res["sites"]] == [
        "consistent",
        "consistent",
        "inconsistent",
        "missing",
    ]
    assert (res["informative_sites"], res["consistent"], res["inconsistent"]) == (3, 2, 1)
    assert res["mendelian_error_rate"] == pytest.approx(1 / 3)
    assert res["inconsistent_indices"] == [2]
    assert (res["sites"][0]["paternal_allele"], res["sites"][0]["maternal_allele"]) == (0, 1)
    assert (res["sites"][1]["paternal_allele"], res["sites"][1]["maternal_allele"]) == (1, 1)


def test_calculate_pedigree_likelihood() -> None:
    ped = [{"id": "f"}, {"id": "m"}, {"id": "c", "father": "f", "mother": "m"}]
    het = result_of(
        hap.calculate_pedigree_likelihood(
            {"pedigree": ped, "genotypes": {"f": [1], "m": [1], "c": [1]}}
        )
    )
    assert het["likelihood"] == pytest.approx(0.125)  # 0.5 * 0.5 * P(child Aa | Aa x Aa) = 0.5
    rare = result_of(
        hap.calculate_pedigree_likelihood(
            {"pedigree": ped, "genotypes": {"f": ["0/0"], "m": ["1/1"], "c": ["0/1"]}},
            options={"allele_frequency": 0.2},
        )
    )
    assert rare["likelihood"] == pytest.approx(0.64 * 0.04 * 1.0)
    unobserved_father = result_of(
        hap.calculate_pedigree_likelihood({"pedigree": ped, "genotypes": {"m": [2], "c": [2]}})
    )
    assert unobserved_father["likelihood"] == pytest.approx(
        0.125
    )  # (0.5*0.5 + 0.25*1) * P(mother = 2) = 0.5 * 0.25
    impossible = result_of(
        hap.calculate_pedigree_likelihood(
            {"pedigree": ped, "genotypes": {"f": [0], "m": [0], "c": [2]}}
        )
    )
    assert impossible["impossible"] is True and impossible["log10_likelihood"] is None
    two_sites = result_of(
        hap.calculate_pedigree_likelihood(
            {"pedigree": ped, "genotypes": {"f": [1, 1], "m": [1, 1], "c": [1, 1]}}
        )
    )
    assert two_sites["log10_likelihood"] == pytest.approx(2 * math.log10(0.125))


def test_cluster_haplotypes() -> None:
    haps = ["0000", "0001", "1111", "1110"]
    res = result_of(hap.cluster(haps, options={"max_distance": 1}))
    assert [c["indices"] for c in res["clusters"]] == [[0, 1], [2, 3]]
    assert result_of(hap.cluster(haps, options={"max_distance": 0}))["cluster_count"] == 4
    assert result_of(hap.cluster(haps, options={"max_distance": 4}))["cluster_count"] == 1
    assert (
        result_of(hap.cluster(["00N0", "0000", "1111"], options={"max_distance": 0}))[
            "cluster_count"
        ]
        == 2
    )


def test_compare_haplotypes_hamming_and_switch_error() -> None:
    res = result_of(hap.compare(["0011", "0101"]))
    assert res["distance"] == 2 and res["differing_positions"] == [1, 2]
    assert res["distance_matrix"] == [[0, 2], [2, 0]]
    sw = result_of(
        hap.compare(
            {"truth": ["0000", "1111"], "inferred": ["0011", "1100"]},
            options={"mode": "switch_error"},
        )
    )
    # relative orientation along the 4 het sites: same, same, flipped, flipped -> 1 switch over 3 intervals
    assert (sw["heterozygous_sites"], sw["switches"]) == (4, 1)
    assert sw["switch_error_rate"] == pytest.approx(1 / 3)
    perfect = result_of(
        hap.compare(
            {"truth": ["0000", "1111"], "inferred": ["1111", "0000"]},
            options={"mode": "switch_error"},
        )
    )
    assert perfect["switches"] == 0


def test_construct_ld_blocks() -> None:
    haps = ["0000", "0011", "1100", "1111"]
    res = result_of(hap.construct(haps))
    assert [b["sites"] for b in res["blocks"]] == [[0, 1], [2, 3]]
    assert [
        b["sites"]
        for b in result_of(hap.construct(haps, options={"statistic": "r2", "threshold": 0.9}))[
            "blocks"
        ]
    ] == [[0, 1], [2, 3]]
    positions = result_of(hap.construct(haps, options={"positions": [100, 200, 300, 400]}))
    assert (positions["blocks"][1]["start_position"], positions["blocks"][1]["end_position"]) == (
        300,
        400,
    )
    # site 0 is monomorphic here and is excluded; sites 1-2 and 3 are analysed
    mono = result_of(hap.construct(["0010", "0100", "0010", "0100"]))
    assert mono["sites_excluded_monomorphic_or_low_maf"] == [0, 3]
    assert [b["sites"] for b in mono["blocks"]] == [[1, 2]]


def test_ld_pair_matches_textbook_values() -> None:
    # pA = pB = 0.5, pAB = 0.4: D = 0.15, r2 = 0.15^2 / (0.25 * 0.25) = 0.36, Dmax = 0.25, D' = 0.6
    x = [1] * 5 + [0] * 5
    y = [1] * 4 + [0] * 1 + [1] * 1 + [0] * 4
    r2, dprime = _impl_haplotypes.ld_pair(x, y)
    assert r2 == pytest.approx(0.36) and dprime == pytest.approx(0.6)


def test_tag_snps() -> None:
    res = result_of(hap.tag(["0000", "0011", "1100", "1111"]))
    assert res["tag_sites"] == [0, 2] and res["bins"] == {"0": [0, 1], "2": [2, 3]}
    redundant = result_of(hap.tag(["000", "111", "000", "111"]))
    assert redundant["tag_sites"] == [0] and redundant["bins"] == {"0": [0, 1, 2]}


def test_detect_compound_heterozygotes() -> None:
    trio = {
        "genes": ["G1", "G1", "G2"],
        "child": ["0/1", "0/1", "0/1"],
        "father": ["0/1", "0/0", "0/1"],
        "mother": ["0/0", "0/1", "0/1"],
    }
    res = result_of(hap.detect_compound_heterozygotes(trio))
    assert res["mode"] == "trio"
    assert res["compound_heterozygous_genes"] == [
        {"gene": "G1", "paternal": [0], "maternal": [1], "undetermined": []}
    ]
    assert res["undetermined_origin_sites"] == 1  # G2: both parents het
    phased = result_of(
        hap.detect_compound_heterozygotes(
            {"genes": ["G1", "G1", "G2", "G2"], "child": ["0|1", "1|0", "0|1", "0|1"]}
        )
    )
    assert [g["gene"] for g in phased["compound_heterozygous_genes"]] == ["G1"]
    with pytest.raises(InputFormatError):
        hap.detect_compound_heterozygotes({"genes": ["G1"], "child": ["0/1"]})


def test_detect_denovo_variants() -> None:
    trio = {
        "child": ["0/1", "0/1", "1/1", "0/0"],
        "father": ["0/0", "0/1", "0/0", "0/0"],
        "mother": ["0/0", "0/0", "0/1", "0/0"],
    }
    res = result_of(hap.detect_denovo_variants(trio))
    assert [d["index"] for d in res["denovo"]] == [0]
    assert res["denovo"][0]["novel_alleles"] == [1]
    assert res["other_mendelian_inconsistent_indices"] == [2]
    assert res["evaluated_sites"] == 4


def test_frequency() -> None:
    res = result_of(hap.frequency(["01", "01", "10", "11"]))
    assert [(r["haplotype"], r["count"]) for r in res["haplotypes"]] == [
        ("01", 2),
        ("10", 1),
        ("11", 1),
    ]
    assert res["haplotypes"][0]["frequency"] == pytest.approx(0.5)
    assert res["haplotype_diversity"] == pytest.approx(4 / 3 * (1 - (0.25 + 0.0625 + 0.0625)))
    by_sample = result_of(hap.frequency({"s1": ["01", "10"], "s2": ["01", "01"]}))
    assert by_sample["total"] == 4 and by_sample["haplotypes"][0] == {
        "haplotype": "01",
        "count": 3,
        "frequency": 0.75,
    }
    counts = result_of(hap.frequency({"AB": 6, "ab": 2}))
    assert counts["haplotypes"][1]["frequency"] == pytest.approx(0.25)


def test_infer_parentage() -> None:
    data = {
        "child": ["0/0", "1/1", "0/1"],
        "candidates": {"A": ["0/1", "1/1", "0/0"], "B": ["1/1", "0/0", "0/1"]},
    }
    res = result_of(hap.infer_parentage(data))
    assert res["candidates"]["A"]["mismatches"] == 0 and res["candidates"]["B"]["mismatches"] == 2
    assert res["not_excluded"] == ["A"] and res["best_candidates"] == ["A"]
    lenient = result_of(hap.infer_parentage(data, options={"max_mismatches": 2}))
    assert lenient["not_excluded"] == ["A", "B"] and lenient["best_candidates"] == ["A"]


def test_phase_trio() -> None:
    trio = {
        "child": ["0/1", "0/1", "0/1", "1/1", "./.", "0/1"],
        "father": ["0/0", "0/1", "0/0", "1/1", "0/0", "0/1"],
        "mother": ["0/1", "0/1", "1/1", "0/1", "0/0", "0/0"],
    }
    res = result_of(hap.phase(trio))
    assert res["phased_child"] == ["0|1", "0/1", "0|1", "1|1", None, "1|0"]
    assert res["status"] == ["phased", "ambiguous", "phased", "homozygous", "missing", "phased"]
    assert res["phased_het_sites"] == 3 and res["ambiguous_het_sites"] == 1
    assert res["fraction_het_sites_phased"] == pytest.approx(3 / 4)


def test_reconstruct_em_resolves_double_heterozygote() -> None:
    genotypes = {"a": [0, 0], "b": [0, 0], "c": [2, 2], "d": [2, 2], "e": [1, 1]}
    res = result_of(hap.reconstruct(genotypes, options={"tolerance": 1e-12}))
    # 4 chromosomes of each 00 and 11 are certain; the double heterozygote is 00/11 (not 01/10): f(00) = f(11) = 0.5
    freqs = {row["haplotype"]: row["frequency"] for row in res["frequencies"]}
    assert freqs["00"] == pytest.approx(0.5, abs=1e-6) and freqs["11"] == pytest.approx(
        0.5, abs=1e-6
    )
    assert res["converged"] is True
    assert res["individuals"]["e"]["haplotypes"] == ["00", "11"]
    assert res["individuals"]["a"]["haplotypes"] == ["00", "00"]


def test_reconstruct_em_first_iteration_matches_hand_calculation() -> None:
    genotypes = {"a": [0, 0], "c": [2, 2], "e": [1, 1]}
    res = result_of(hap.reconstruct(genotypes, options={"max_iterations": 1}))
    # uniform start over {00, 01, 10, 11}: posterior of 00/11 and 01/10 for e is 1/2 each ->
    # counts 00: 2 + 0.5, 11: 2 + 0.5, 01: 0.5, 10: 0.5 over 6 chromosomes
    freqs = {row["haplotype"]: row["frequency"] for row in res["frequencies"]}
    assert freqs == pytest.approx({"00": 2.5 / 6, "11": 2.5 / 6, "01": 0.5 / 6, "10": 0.5 / 6})


def test_select_haplotypes() -> None:
    sample = ["01"] * 5 + ["10"] * 4 + ["11"]
    res = result_of(hap.select(sample, options={"min_frequency": 0.15}))
    assert [r["haplotype"] for r in res["selected"]] == ["01", "10"]
    assert res["covered_frequency"] == pytest.approx(0.9)
    assert [
        r["haplotype"]
        for r in result_of(hap.select(sample, options={"min_frequency": 0.0, "top_n": 1}))[
            "selected"
        ]
    ] == ["01"]


def test_simulate_haplotypes_is_deterministic_and_respects_support() -> None:
    fixed = result_of(hap.simulate({"01": 1.0, "10": 0.0}, options={"n_individuals": 5}))
    assert all(ind["haplotypes"] == ["01", "01"] for ind in fixed["individuals"])
    assert fixed["observed_haplotype_counts"] == {"01": 10}
    a = result_of(hap.simulate({"01": 0.3, "10": 0.7}, options={"n_individuals": 50, "seed": 7}))
    b = result_of(hap.simulate({"01": 3, "10": 7}, options={"n_individuals": 50, "seed": 7}))
    assert a["individuals"] == b["individuals"]
    assert sum(a["observed_haplotype_counts"].values()) == 100


def test_visualize_haplotypes(tmp_path: Path) -> None:
    target = tmp_path / "haps.svg"
    res = result_of(
        hap.visualize({"h1": "01", "h2": "10", "h3": "1N"}, options={"output": str(target)})
    )
    assert res["svg"].count("<rect") == 6 and (res["haplotype_count"], res["site_count"]) == (3, 2)
    assert "#999999" in res["svg"] and target.read_text(encoding="utf-8") == res["svg"]


# --------------------------------------------------------------------------- phylogeny


WIKI_NJ = {
    "names": ["a", "b", "c", "d", "e"],
    "matrix": [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0],
    ],
}
WIKI_UPGMA = {
    "names": ["a", "b", "c", "d", "e"],
    "matrix": [
        [0, 17, 21, 31, 23],
        [17, 0, 30, 34, 21],
        [21, 30, 0, 28, 39],
        [31, 34, 28, 0, 43],
        [23, 21, 39, 43, 0],
    ],
}


def test_build_neighbor_joining_textbook_example() -> None:
    res = result_of(phylo.build(WIKI_NJ))
    assert res["rooted"] is False and res["leaf_count"] == 5
    # Saitou & Nei example: ((a:2,b:3):3, c:4, (d:2,e:1):2) as an unrooted tree
    expected = "((a:2,b:3):3,c:4,(d:2,e:1):2);"
    assert result_of(phylo.compare([res["newick"], expected]))["robinson_foulds"] == 0
    patristic = result_of(phylo.calculate(res["newick"], options={"what": "patristic_matrix"}))
    order = {name: i for i, name in enumerate(patristic["names"])}
    for i, a in enumerate(WIKI_NJ["names"]):
        for j, b in enumerate(WIKI_NJ["names"]):
            assert patristic["matrix"][order[a]][order[b]] == pytest.approx(WIKI_NJ["matrix"][i][j])
    assert res["total_branch_length"] == pytest.approx(2 + 3 + 3 + 4 + 2 + 2 + 1)


def test_build_upgma_textbook_example() -> None:
    res = result_of(phylo.build(WIKI_UPGMA, options={"method": "upgma"}))
    assert res["rooted"] is True
    assert res["newick"] == "((e:11,(a:8.5,b:8.5):2.5):5.5,(c:14,d:14):2.5);"
    stats = result_of(phylo.calculate(res["newick"]))
    assert stats["ultrametric"] is True and stats["max_root_to_tip"] == pytest.approx(16.5)


def test_build_two_taxa_and_nested_mapping() -> None:
    assert result_of(phylo.build({"x": {"y": 4.0}, "y": {"x": 4.0}}))["newick"] == "(x:2,y:2);"


def test_calculate_tree_statistics() -> None:
    res = result_of(phylo.calculate("((A:1,B:2)X:3,C:4);"))
    assert (res["leaf_count"], res["internal_node_count"], res["node_count"]) == (3, 2, 5)
    assert res["leaves"] == ["A", "B", "C"] and res["rooted_binary"] is True
    assert res["total_branch_length"] == pytest.approx(10.0)
    # root-to-tip: A 1 + 3 = 4, B 2 + 3 = 5, C 4
    assert (res["min_root_to_tip"], res["max_root_to_tip"]) == (
        pytest.approx(4.0),
        pytest.approx(5.0),
    )
    assert res["mean_root_to_tip"] == pytest.approx(13 / 3) and res["ultrametric"] is False
    assert res["treeness"] == pytest.approx(3 / 10)
    assert (
        res["sackin_index"] == 2 + 2 + 1 and res["colless_index"] == 1
    )  # |2-1| at the root, |1-1| below
    ultra = result_of(phylo.calculate("((A:1,B:1):3,C:4);"))
    assert ultra["ultrametric"] is True and ultra["treeness"] == pytest.approx(3 / 9)
    caterpillar = result_of(phylo.calculate("(((A,B),C),D);"))
    assert caterpillar["has_branch_lengths"] is False and caterpillar["max_depth_edges"] == 3
    star = result_of(phylo.calculate("(A:1,B:1,C:1,D:1);"))
    assert (
        star["root_degree"] == 4 and star["polytomy_count"] == 1 and star["rooted_binary"] is False
    )


def test_calculate_distance_matrix_p_and_jukes_cantor() -> None:
    aln = {"A": "ACGTACGTAC", "B": "ACGTACGTTT", "C": "-CGTACGTAC"}
    p = result_of(phylo.calculate(aln, options={"distance": "p"}))
    assert p["names"] == ["A", "B", "C"]
    assert p["matrix"][0][1] == pytest.approx(0.2)  # 2 mismatches / 10 sites
    assert p["matrix"][0][2] == 0.0 and p["compared_sites"][0][2] == 9  # gap skipped pairwise
    assert p["matrix"][1][2] == pytest.approx(2 / 9)
    jc = result_of(phylo.calculate(aln))
    assert jc["matrix"][0][1] == pytest.approx(-0.75 * math.log(1 - 4 * 0.2 / 3))
    assert jc["matrix"][0][1] == pytest.approx(0.232616, abs=1e-6)
    fasta = result_of(phylo.calculate(">A\nACGT\n>B\nACGA\n", options={"distance": "p"}))
    assert fasta["matrix"][0][1] == pytest.approx(0.25)


def test_calculate_jukes_cantor_saturation() -> None:
    aln = {"A": "AAAA", "B": "CCCC"}
    with pytest.raises(InputFormatError):
        phylo.calculate(aln)
    capped = result_of(phylo.calculate(aln, options={"saturation_distance": 5.0}))
    assert capped["matrix"][0][1] == 5.0


def test_compare_robinson_foulds() -> None:
    res = result_of(phylo.compare(["((A,B),(C,D),(E,F));", "(((A,B),C),(D,(E,F)));"]))
    # splits: {AB, CD, EF} vs {ABC|DEF, AB, EF}; shared AB, EF -> RF = 2
    assert res["robinson_foulds"] == 2 and res["shared_splits"] == 2
    assert res["normalized_rf"] == pytest.approx(2 / 6) and res["max_rf_binary"] == 6
    quartet = result_of(phylo.compare({"a": "((A,B),(C,D));", "b": "((A,C),(B,D));"}))
    assert quartet["robinson_foulds"] == 2 and quartet["normalized_rf"] == 1.0
    assert (
        result_of(phylo.compare(["((A,B),(C,D));", "(A,(B,(C,D)));"]))["robinson_foulds"] == 0
    )  # same unrooted split AB|CD
    rooted = result_of(
        phylo.compare(["((A,B),(C,D));", "(A,(B,(C,D)));"], options={"rooted": True})
    )
    assert (
        rooted["robinson_foulds"] == 2
    )  # rooted clades {AB, CD} vs {CD, BCD}: AB and BCD are unique
    ignoring_lengths = result_of(
        phylo.compare(["((A:1,B:2):3,(C:1,D:1):1);", "((B:9,A:9):1,(D:1,C:1):1);"])
    )
    assert ignoring_lengths["robinson_foulds"] == 0


def test_annotate_tree() -> None:
    res = result_of(
        phylo.annotate("((A:1,B:2):3,C:4);", options={"metadata": {"A": "x", "B": "x", "C": "y"}})
    )
    nodes = res["nodes"]
    assert [n["name"] for n in nodes] == [
        "N0",
        "N1",
        "A",
        "B",
        "C",
    ]  # unnamed internal nodes are labelled
    internal = nodes[1]
    assert (
        internal["leaf_count"] == 2
        and internal["clade"] == ["A", "B"]
        and internal["root_distance"] == 3.0
    )
    assert internal["metadata_values"] == ["x"] and internal["metadata_homogeneous"] is True
    assert nodes[0]["metadata_homogeneous"] is False and nodes[0]["leaf_count"] == 3
    assert nodes[3]["depth_edges"] == 2 and nodes[3]["root_distance"] == pytest.approx(5.0)
    assert res["annotated_newick"] == "((A:1,B:2)N1:3,C:4)N0;"
    with pytest.raises(ConfigurationError):
        phylo.annotate("((A,B),C);", options={"metadata": {"Z": 1}})


def test_export_formats(tmp_path: Path) -> None:
    tree = "((A:1,B:2):3,'C d':4);"
    newick = result_of(phylo.export(tree))
    assert newick["content"] == "((A:1,B:2):3,'C d':4);"
    nexus = result_of(phylo.export(tree, options={"format": "nexus"}))["content"]
    assert (
        "DIMENSIONS NTAX=3;" in nexus
        and "TAXLABELS A B 'C d';" in nexus
        and nexus.startswith("#NEXUS")
    )
    as_json = result_of(phylo.export(tree, options={"format": "json"}))["content"]
    assert as_json["children"][0]["children"][1] == {"name": "B", "branch_length": 2.0}
    edges = result_of(phylo.export(tree, options={"format": "edgelist"}))["content"]
    assert (edges[0]["parent"], edges[0]["child"], edges[0]["length"]) == (0, 1, 3.0) and len(
        edges
    ) == 4
    target = tmp_path / "t.nwk"
    result_of(phylo.export("(A,B);", options={"output": str(target)}))
    assert target.read_text(encoding="utf-8") == "(A,B);"


ALN = {"A": "ACGTACGTACGT", "B": "ACGTACGTACGA", "C": "TGCATGGTACGT", "D": "TGCATGGTACGA"}


def test_infer_distance_tree_groups_close_sequences() -> None:
    res = result_of(phylo.infer(ALN))
    assert result_of(phylo.compare([res["newick"], "((A,B),(C,D));"]))["robinson_foulds"] == 0
    assert res["distance"] == "jc69" and res["leaf_count"] == 4
    assert res["distance_matrix"][0][1] == pytest.approx(-0.75 * math.log(1 - 4 / 3 * (1 / 12)))
    upgma = result_of(phylo.infer(ALN, options={"method": "upgma", "distance": "p"}))
    assert upgma["rooted"] is True
    assert result_of(phylo.compare([upgma["newick"], "((A,B),(C,D));"]))["robinson_foulds"] == 0


def test_test_bootstrap_support_is_reproducible_and_hand_checkable() -> None:
    strong = result_of(phylo.test(ALN, options={"replicates": 50, "seed": 0}))
    assert [c["support"] for c in strong["clades"]] == [1.0]
    assert strong["clades"][0]["clade"] in (["A", "B"], ["C", "D"])
    assert strong["mean_support"] == 1.0 and "100" in strong["newick_with_support"]
    # 6 columns say AB|CD, 4 columns say AC|BD (the full-alignment NJ tree is AB|CD); with p distance
    # a replicate recovers AB|CD when it draws more AB|CD columns than AC|BD columns: P(>5 of 10) = 0.633
    # (+ ties), so the support is far from 1
    conflicting = {"A": "AAAAAAGGGG", "B": "AAAAAATTTT", "C": "CCCCCCGGGG", "D": "CCCCCCTTTT"}
    first = result_of(
        phylo.test(conflicting, options={"replicates": 400, "seed": 1, "distance": "p"})
    )
    again = result_of(
        phylo.test(conflicting, options={"replicates": 400, "seed": 1, "distance": "p"})
    )
    assert first == again
    assert 0.5 < first["clades"][0]["support"] < 0.9
    other_seed = result_of(
        phylo.test(conflicting, options={"replicates": 400, "seed": 2, "distance": "p"})
    )
    assert other_seed["clades"][0]["support"] != 1.0


def test_plot_tree_svg(tmp_path: Path) -> None:
    target = tmp_path / "tree.svg"
    res = result_of(phylo.plot("((A:1,B:2):3,C:4);", options={"output": str(target)}))
    assert res["leaf_count"] == 3 and res["svg"].count("<text") == 3
    assert (
        ">A<" in res["svg"]
        and ">C<" in res["svg"]
        and target.read_text(encoding="utf-8") == res["svg"]
    )
    clado = result_of(phylo.plot("((A,B),C);", options={"branch_lengths": False}))
    assert clado["svg"].count("<text") == 3
    with pytest.raises(InputFormatError):
        phylo.plot("((A,B),C);")  # phylogram without branch lengths


# --------------------------------------------------------------------------- not implemented and invalid input


NOT_IMPLEMENTED = [
    (asm.build_chromosome_level_assembly, "x"),
    (asm.build_haplotype_assembly, "x"),
    (asm.detect_assembly_contamination, "x"),
    (asm.detect_chimeric_contigs, "x"),
    (asm.detect_misassemblies, "x"),
    (asm.identify_centromeres, "x"),
    (asm.phase_assembly, "x"),
    (asm.polish, "x"),
    (asm.purge_haplotigs, "x"),
    (hap.impute, "x"),
]


@pytest.mark.parametrize(("function", "payload"), NOT_IMPLEMENTED)
def test_leftover_functions_report_not_implemented(function, payload) -> None:
    out = function(payload)
    assert out["status"] == "not_implemented"
    assert out["required_dependency"] and out["message"]
    assert "result" not in out


def test_infer_with_ml_method_is_not_implemented() -> None:
    out = phylo.infer(ALN, options={"method": "ml"})
    assert out["status"] == "not_implemented"


@pytest.mark.parametrize(
    ("call", "error"),
    [
        (lambda: phylo.calculate("((A,B);"), InputFormatError),
        (lambda: phylo.calculate("((A,B),(A,C));"), InputFormatError),  # duplicate leaf
        (lambda: phylo.compare(["((A,B),C);", "((A,B),D);"]), InputFormatError),
        (lambda: phylo.compare(["(A,B);"]), InputFormatError),
        (lambda: phylo.build({"names": ["a", "b"], "matrix": [[0, 1], [2, 0]]}), InputFormatError),
        (
            lambda: phylo.build(
                {"names": ["a", "b"], "matrix": [[0, 1], [1, 0]]}, options={"method": "ml"}
            ),
            ConfigurationError,
        ),
        (lambda: phylo.infer({"A": "ACGT", "B": "ACG"}), InputFormatError),
        (lambda: phylo.test(ALN, options={"replicates": 0}), ConfigurationError),
        (lambda: asm.calculate_assembly_completeness([100]), ConfigurationError),
        (lambda: asm.calculate_busco_completeness("no summary here"), InputFormatError),
        (
            lambda: asm.estimate_size({5: 100, 4: 200, 3: 300}),
            InputFormatError,
        ),  # monotone: no valley
        (lambda: asm.calculate_pangenome_openness(PAN), InputFormatError),  # only three genomes
        (
            lambda: asm.extract_core_pangenome(PAN, options={"core_threshold": 1.5}),
            ConfigurationError,
        ),
        (lambda: asm.extract_core_pangenome("not a matrix"), InputFormatError),
        (
            lambda: asm.remove_contamination({"c": "AC"}, options={"contaminants": ["c"]}),
            ConfigurationError,
        ),
        (
            lambda: asm.scaffold({"c": "ACGT"}, options={"layout": [{"contigs": ["zz"]}]}),
            ConfigurationError,
        ),
        (lambda: asm.build_command({"tool": "bwa", "output_directory": "o"}), ConfigurationError),
        (lambda: asm.identify_repeats_in_assembly([100, 200]), InputFormatError),  # lengths only
        (lambda: asm.compare_assemblies([100, 200]), InputFormatError),
        (
            lambda: hap.analyze_trio_inheritance(
                {"child": ["0/1"], "father": ["0/0", "0/1"], "mother": ["0/0"]}
            ),
            InputFormatError,
        ),
        (
            lambda: hap.analyze_trio_inheritance(
                {"child": ["0/x"], "father": ["0/0"], "mother": ["0/0"]}
            ),
            InputFormatError,
        ),
        (lambda: hap.phase({"child": ["0/1"]}), InputFormatError),
        (lambda: hap.construct(["01", "012"]), InputFormatError),
        (lambda: hap.compare(["01"]), InputFormatError),
        (lambda: hap.compare(["01", "10"], options={"mode": "bogus"}), ConfigurationError),
        (lambda: hap.reconstruct({"a": [0, None]}), InputFormatError),
        (lambda: hap.frequency("not-a-list-of-haplotypes"), InputFormatError),
        (lambda: hap.simulate({"01": -1}), InputFormatError),
        (
            lambda: hap.calculate_pedigree_likelihood(
                {"pedigree": [{"id": "c", "father": "f"}], "genotypes": {"c": [1]}}
            ),
            InputFormatError,
        ),
    ],
)
def test_invalid_inputs_raise_dnakit_errors(call, error) -> None:
    with pytest.raises(error) as info:
        call()
    assert getattr(info.value, "code", None)
    assert getattr(info.value, "hint", None) is not None or isinstance(
        info.value, (InputFormatError, ConfigurationError)
    )


def test_every_ok_envelope_has_module_operation_and_status() -> None:
    out = phylo.compare(["((A,B),C);", "((A,B),C);"])
    assert (out["module"], out["operation"], out["status"]) == ("phylogeny", "compare", "ok")
    assert json.dumps(out)  # JSON-serialisable
