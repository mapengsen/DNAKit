"""MAGeCK, DrugZ and Cas-OFFinder facades (each skips when its program or script is missing)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from _p16_data import need, rand_dna, mageck_data, write_fasta

import dnakit.design_crispr.p16_crispr_screen_cli_functions as cr
from dnakit.exceptions import ConfigurationError


def test_mageck_count_rank_mle_pathway(tmp_path: Path) -> None:
    program = need("mageck")
    data = mageck_data(tmp_path)
    labels = ["c1", "c2", "t1", "t2"]
    counted = cr.count_sgrnas_mageck(data["lib.txt"], str(tmp_path / "cnt"), fastq=[data[f"{n}.fq"] for n in labels], sample_labels=labels, trim_5=2,
                                     mageck_path=program, verbose=False)
    assert [row["Reads"] for row in counted["summary"]] == [row["Mapped"] for row in counted["summary"]] and len(counted["summary"]) == 4
    ranked = cr.rank_genes_mageck(counted["count_table"], str(tmp_path / "rra"), treatment=["t1", "t2"], control=["c1", "c2"], mageck_path=program, verbose=False)
    assert {row["id"] for row in ranked["genes"][:2]} == {"G1", "G2"}
    mle = cr.estimate_gene_essentiality_mageck(counted["count_table"], data["design.txt"], str(tmp_path / "mle"), permutation_rounds=2, mageck_path=program, verbose=False)
    assert {row["Gene"] for row in mle["genes"]}
    top = min(cr.cli.read_table(Path(mle["gene_summary"])), key=lambda r: r["treat|beta"])["Gene"]
    assert top in {"G1", "G2"}
    pathway = cr.enrich_pathways_mageck(ranked["gene_summary"], data["sets.gmt"], str(tmp_path / "pw"), method="rra", mageck_path=program, verbose=False)
    assert {"setA", "setB"} <= {row["id"] for row in pathway["pathways"]}
    with pytest.raises(ConfigurationError):
        cr.rank_genes_mageck(counted["count_table"], str(tmp_path / "bad"), treatment=["t1"], verbose=False)


def test_mageck_plot_and_cnv_inputs(tmp_path: Path) -> None:
    program = need("mageck")
    data = mageck_data(tmp_path)
    labels = ["c1", "c2", "t1", "t2"]
    counted = cr.count_sgrnas_mageck(data["lib.txt"], str(tmp_path / "cnt"), fastq=[data[f"{n}.fq"] for n in labels], sample_labels=labels, trim_5=2,
                                     mageck_path=program, verbose=False)
    ranked = cr.rank_genes_mageck(counted["count_table"], str(tmp_path / "rra"), treatment=["t1", "t2"], control=["c1", "c2"], mageck_path=program, verbose=False)
    try:
        plot = cr.plot_screen_results_mageck(counted["count_table"], ranked["gene_summary"], str(tmp_path / "pl"), genes=["G1"], mageck_path=program, verbose=False)
    except Exception as exc:  # MAGeCK plot needs R for some figures
        pytest.skip(f"mageck plot unavailable: {exc}")
    assert any(path.endswith(".pdf") for path in plot["outputs"])
    cnv = tmp_path / "cnv.txt"
    cnv.write_text("SYMBOL\tline1\n" + "".join(f"G{i}\t{1.0 + (i * 7 % 40) * 0.1}\n" for i in range(40)))
    # The copy-number regression of MAGeCK needs genome-scale input, so a 40-gene toy matrix is rejected by MAGeCK itself;
    # the test checks that the flags reach the program.
    from dnakit.exceptions import BackendExecutionError

    with pytest.raises(BackendExecutionError) as raised:
        cr.rank_genes_mageck(counted["count_table"], str(tmp_path / "cnvrra"), treatment=["t1", "t2"], control=["c1", "c2"],
                             cnv_normalization=str(cnv), cell_line="line1", mageck_path=program, verbose=False)
    command = raised.value.context["command"]
    assert "--cnv-norm" in command and command[command.index("--cell-line") + 1] == "line1"


def test_drugz_paired_and_unpaired(tmp_path: Path) -> None:
    pytest.importorskip("pandas")
    pytest.importorskip("scipy")
    pytest.importorskip("six")
    script = os.environ.get("DRUGZ_SCRIPT")
    if not script or not Path(script).is_file():
        pytest.skip("set DRUGZ_SCRIPT to the downloaded drugz.py (https://github.com/hart-lab/drugz)")
    data = mageck_data(tmp_path)
    paired = cr.score_drug_gene_interactions_drugz(data["drugz_counts.txt"], str(tmp_path / "p.txt"), control_samples=["c1", "c2"], drug_samples=["d1", "d2"],
                                                   drugz_script=script, half_window_size=20, verbose=False)
    unpaired = cr.score_drug_gene_interactions_drugz(data["drugz_counts.txt"], str(tmp_path / "u.txt"), control_samples=["c1", "c2"], drug_samples=["d1", "d2"],
                                                     drugz_script=script, half_window_size=20, unpaired=True, verbose=False)
    assert paired["genes"][0]["GENE"] == unpaired["genes"][0]["GENE"] == "G5" and paired["genes"][0]["normZ"] < -3
    assert paired["genes"][0]["numObs"] == 8 and unpaired["genes"][0]["numObs"] == 4


def test_cas_offinder_finds_the_guide_site(tmp_path: Path) -> None:
    program = need("cas-offinder")
    import random

    rng = random.Random(6)
    chrom = rand_dna(rng, 6000)
    site = chrom[1000:1020]
    chrom = chrom[:1020] + "TGG" + chrom[1023:]
    fasta = write_fasta(tmp_path / "g.fa", {"chrT": chrom})
    result = cr.search_offtargets_cas_offinder(fasta, {site + "NNN": 2}, str(tmp_path / "hits.txt"), pattern="N" * 20 + "NGG", cas_offinder_path=program, verbose=False)
    exact = [s for s in result["sites"] if s["mismatches"] == 0]
    assert exact and exact[0]["position"] == 1000 and exact[0]["strand"] == "+" and exact[0]["site"].upper().endswith("TGG")
    with pytest.raises(ConfigurationError, match="same length"):
        cr.search_offtargets_cas_offinder(fasta, {"ACGT": 1}, str(tmp_path / "h2.txt"), verbose=False)
