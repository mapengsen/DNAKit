"""Tests for dnakit.gwas.gcta_p08_cli_functions (GCTA command-line operations).

The tests run the real GCTA (1.95.4 tested) on a small synthetic PLINK fileset.
They skip when `gcta64` is not installed.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

import _gcta_p08_data as data

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.gwas import gcta_p08_cli_functions as g

requires_gcta = pytest.mark.skipif(not data.gcta_available(), reason="gcta64 is not installed")


@pytest.fixture(scope="module")
def toy(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return data.write_toy(tmp_path_factory.mktemp("gcta_toy"))


@pytest.fixture(scope="module")
def grm(toy: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> str:
    if not data.gcta_available():
        pytest.skip("gcta64 is not installed")
    out = tmp_path_factory.mktemp("gcta_grm")
    return g.gcta_make_grm(str(out / "g"), bfile=str(toy["bfile"]), verbose=False)["grm_prefix"]


@pytest.fixture(scope="module")
def sparse(grm: str, tmp_path_factory: pytest.TempPathFactory) -> str:
    out = tmp_path_factory.mktemp("gcta_sparse")
    return g.gcta_make_sparse_grm(grm, str(out / "sp"), verbose=False)["grm_prefix"]


@pytest.fixture(scope="module")
def summary(toy: dict[str, Path], grm: str, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("gcta_ma")
    result = g.gcta_mlma(
        str(toy["bfile"]), str(toy["toy.phen"]), str(out / "ml"), grm_prefix=grm, verbose=False
    )
    return data.ma_from_mlma(Path(result["outputs"]["mlma"]), out / "ml.ma")


def _ok(result: dict) -> dict:
    assert result["status"] == "completed", result["log_tail"][-1500:]
    assert result["backend"] == "gcta" and result["version"]
    return result


def _bfile(toy: dict[str, Path]) -> str:
    return str(toy["bfile"])


def _read_grm(prefix: str) -> list[float]:
    import struct

    raw = Path(f"{prefix}.grm.bin").read_bytes()
    return list(struct.unpack(f"<{len(raw) // 4}f", raw))


def test_missing_program_and_validation(toy: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        g.gcta_make_grm(
            str(tmp_path / "g"), bfile=_bfile(toy), gcta_path=str(tmp_path / "none"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        g.gcta_make_grm(str(tmp_path / "g"), verbose=False)  # neither bfile nor mbfile
    with pytest.raises(ConfigurationError):
        g.gcta_make_grm(str(tmp_path / "g"), bfile=_bfile(toy), grm_type="sparse", verbose=False)
    with pytest.raises(ConfigurationError):
        g.gcta_make_grm(str(tmp_path / "g"), bfile=str(tmp_path / "nofiles"), verbose=False)
    (tmp_path / "g.log").write_text("old")
    with pytest.raises(ConfigurationError):
        g.gcta_make_grm(str(tmp_path / "g"), bfile=_bfile(toy), verbose=False)
    with pytest.raises(ConfigurationError):
        g.gcta_fastgwa(
            str(toy["toy.phen"]), str(tmp_path / "f"), bfile=_bfile(toy), model="mlm", verbose=False
        )
    with pytest.raises(ConfigurationError):
        g.gcta_fastgwa(
            str(toy["toy.phen"]),
            str(tmp_path / "f"),
            bfile=_bfile(toy),
            model="lr",
            sparse_grm_prefix="x",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        g.gcta_cojo(
            _bfile(toy), str(toy["toy.phen"]), str(tmp_path / "c"), mode="cond", verbose=False
        )
    with pytest.raises(ConfigurationError):
        g.gcta_fastbat(_bfile(toy), str(toy["toy.phen"]), str(tmp_path / "b"), verbose=False)
    with pytest.raises(ConfigurationError):
        g.gcta_simulate_phenotypes(
            _bfile(toy), str(tmp_path / "s"), trait="case_control", verbose=False
        )
    with pytest.raises(ConfigurationError):
        g.gcta_unify_grm(["only_one"], str(tmp_path / "u"), verbose=False)


@requires_gcta
def test_recode_standardized_genotypes_matches_numpy(toy: dict[str, Path], tmp_path: Path) -> None:
    result = _ok(
        g.gcta_recode_standardized_genotypes(
            _bfile(toy), str(tmp_path / "rs"), chromosome=1, maf=0.05, verbose=False
        )
    )
    assert result["n_individuals"] == 200 and result["n_snps"] == 1300
    with gzip.open(result["outputs"]["xmat"], "rt") as handle:
        rows = [line.split() for line in handle]
    matrix = [[float(v) for v in row[2:]] for row in rows[2:]]
    assert len(matrix) == 200
    columns = list(zip(*matrix, strict=True))
    for column in columns[:20]:
        assert sum(column) / 200 == pytest.approx(0.0, abs=1e-3)  # standardised columns are centred


@requires_gcta
def test_make_grm_types(toy: dict[str, Path], tmp_path: Path) -> None:
    additive = _ok(g.gcta_make_grm(str(tmp_path / "a"), bfile=_bfile(toy), verbose=False))
    assert additive["n_individuals"] == 200
    values = _read_grm(additive["grm_prefix"])
    assert len(values) == 200 * 201 // 2
    diagonal = [values[i * (i + 1) // 2 + i] for i in range(200)]
    assert sum(diagonal) / 200 == pytest.approx(1.0, abs=0.05)  # random unrelated samples
    dominance = _ok(
        g.gcta_make_grm(str(tmp_path / "d"), bfile=_bfile(toy), grm_type="dominance", verbose=False)
    )
    assert (
        dominance["grm_prefix"].endswith("d.d")
        and Path(dominance["outputs"]["grm_bin"]).name == "d.d.grm.bin"
    )
    for kind in ("xchr", "inbred"):
        _ok(g.gcta_make_grm(str(tmp_path / kind), bfile=_bfile(toy), grm_type=kind, verbose=False))
    listing = tmp_path / "chromosomes.txt"
    listing.write_text(_bfile(toy) + "\n")
    multi = _ok(g.gcta_make_grm(str(tmp_path / "m"), mbfile=str(listing), verbose=False))
    assert _read_grm(multi["grm_prefix"]) == pytest.approx(values)


@requires_gcta
def test_sparse_grm_and_unify(grm: str, toy: dict[str, Path], tmp_path: Path) -> None:
    sparse = _ok(g.gcta_make_sparse_grm(grm, str(tmp_path / "sp"), cutoff=0.05, verbose=False))
    assert sparse["n_pairs"] >= 200  # at least the diagonal
    assert Path(sparse["outputs"]["grm_sp"]).name == "sp.grm.sp"
    dense = _ok(
        g.gcta_make_sparse_grm(grm, str(tmp_path / "bk"), output_format="dense", verbose=False)
    )
    thresholded = _read_grm(dense["grm_prefix"])
    assert sum(1 for v in thresholded if v != 0) <= len(thresholded)
    subset = _ok(
        g.gcta_make_grm(
            str(tmp_path / "s150"), bfile=_bfile(toy), keep=str(toy["keep150.txt"]), verbose=False
        )
    )
    assert subset["n_individuals"] == 150
    unified = _ok(
        g.gcta_unify_grm([grm, subset["grm_prefix"]], str(tmp_path / "un"), verbose=False)
    )
    assert unified["common_individuals"] == 150 and len(unified["unified"]) == 2
    for prefix in unified["unified"]:
        assert len(_read_grm(prefix)) == 150 * 151 // 2
    with pytest.raises(ConfigurationError):
        g.gcta_unify_grm([grm, subset["grm_prefix"]], str(tmp_path / "un"), verbose=False)


@requires_gcta
def test_inbreeding_coefficient(toy: dict[str, Path], tmp_path: Path) -> None:
    result = _ok(g.gcta_inbreeding_coefficient(_bfile(toy), str(tmp_path / "ibc"), verbose=False))
    assert result["summary"]["individuals"] == 200
    assert set(result["summary"]["mean"]) == {"Fhat1", "Fhat2", "Fhat3"}
    assert abs(result["summary"]["mean"]["Fhat3"]) < 0.05


@requires_gcta
def test_reml_heritability_and_predictions(toy: dict[str, Path], grm: str, tmp_path: Path) -> None:
    result = _ok(
        g.gcta_reml(
            str(toy["toy.phen"]), str(tmp_path / "r"), grm_prefix=grm, predict="blup", verbose=False
        )
    )
    components = {row["source"]: row for row in result["summary"]["variance_components"]}
    assert components["V(G)/Vp"]["variance"] == pytest.approx(
        components["V(G)"]["variance"] / components["Vp"]["variance"], abs=2e-3
    )
    assert result["summary"]["statistics"]["n"] == 200 and result["n_predictions"] == 200
    cv = _ok(
        g.gcta_reml(
            str(toy["toy.phen"]),
            str(tmp_path / "cv"),
            grm_prefix=grm,
            predict="cvblup",
            verbose=False,
        )
    )
    assert cv["n_predictions"] == 200 and "indi_cvblp" in cv["outputs"]
    gxe = _ok(
        g.gcta_reml(
            str(toy["toy.phen"]),
            str(tmp_path / "gx"),
            grm_prefix=grm,
            gxe_file=str(toy["toy.env_class"]),
            verbose=False,
        )
    )
    assert any(row["source"].startswith("V(GxE)") for row in gxe["summary"]["variance_components"])
    with pytest.raises(ConfigurationError):
        g.gcta_reml(str(toy["toy.phen"]), str(tmp_path / "bad"), verbose=False)
    with pytest.raises(ConfigurationError):
        g.gcta_reml(
            str(toy["toy.phen"]),
            str(tmp_path / "bad"),
            grm_prefix=grm,
            predict="mean",
            verbose=False,
        )


@requires_gcta
def test_blup_snp_effects_are_proportional_to_the_blup(
    toy: dict[str, Path], grm: str, tmp_path: Path
) -> None:
    pred = _ok(
        g.gcta_reml(
            str(toy["toy.phen"]), str(tmp_path / "r"), grm_prefix=grm, predict="blup", verbose=False
        )
    )
    result = _ok(
        g.gcta_blup_snp(
            _bfile(toy),
            pred["outputs"]["indi_blp"],
            str(tmp_path / "bs"),
            autosome=True,
            verbose=False,
        )
    )
    assert result["n_snps"] == 2500
    failing = g.gcta_blup_snp(
        _bfile(toy), pred["outputs"]["indi_blp"], str(tmp_path / "bs2"), verbose=False
    )
    assert failing["status"] == "failed"  # the X chromosome needs autosome=True


@requires_gcta
def test_bivariate_reml_and_haseman_elston(toy: dict[str, Path], grm: str, tmp_path: Path) -> None:
    biv = _ok(
        g.gcta_reml_bivar(
            str(toy["toy2.phen"]), str(tmp_path / "bv"), grm_prefix=grm, verbose=False
        )
    )
    sources = [row["source"] for row in biv["summary"]["variance_components"]]
    assert "rG" in sources and "V(G)_tr1" in sources
    he = _ok(
        g.gcta_he_regression(
            str(toy["toy.phen"]), str(tmp_path / "he"), grm_prefix=grm, verbose=False
        )
    )
    assert set(he["summary"]) == {"HE-CP", "HE-SD"}
    assert {row["coefficient"] for row in he["summary"]["HE-CP"]} == {"Intercept", "V(G)/Vp"}
    pair = _ok(
        g.gcta_he_regression(
            str(toy["toy2.phen"]),
            str(tmp_path / "hb"),
            grm_prefix=grm,
            bivariate=True,
            verbose=False,
        )
    )
    assert any(row["coefficient"] == "rG" for row in pair["summary"]["HE-CP"])


def _p_values(path: str, column: str) -> list[float]:
    lines = Path(path).read_text().splitlines()
    index = lines[0].split().index(column)
    return [float(line.split()[index]) for line in lines[1:]]


@requires_gcta
def test_fastgwa_models(toy: dict[str, Path], sparse: str, tmp_path: Path) -> None:
    linear = _ok(
        g.gcta_fastgwa(
            str(toy["toy.phen"]), str(tmp_path / "lr"), bfile=_bfile(toy), model="lr", verbose=False
        )
    )
    assert linear["summary"]["rows"] == 2600
    assert linear["summary"]["top_rows"][0]["P"] == min(
        _p_values(linear["outputs"]["fastgwa"], "P")
    )
    mixed = _ok(
        g.gcta_fastgwa(
            str(toy["toy.phen"]),
            str(tmp_path / "mlm"),
            bfile=_bfile(toy),
            sparse_grm_prefix=sparse,
            verbose=False,
        )
    )
    exact = _ok(
        g.gcta_fastgwa(
            str(toy["toy.phen"]),
            str(tmp_path / "ex"),
            bfile=_bfile(toy),
            model="mlm_exact",
            sparse_grm_prefix=sparse,
            verbose=False,
        )
    )
    assert mixed["summary"]["rows"] == exact["summary"]["rows"] == 2600
    interaction = _ok(
        g.gcta_fastgwa(
            str(toy["toy.phen"]),
            str(tmp_path / "ge"),
            bfile=_bfile(toy),
            model="lr",
            envir_file=str(toy["toy.envir"]),
            verbose=False,
        )
    )
    assert "BETA_G_by_E" in interaction["summary"]["columns"]
    binary = _ok(
        g.gcta_fastgwa(
            str(toy["toy.bin.phen"]),
            str(tmp_path / "bin"),
            bfile=_bfile(toy),
            model="mlm_binary",
            sparse_grm_prefix=sparse,
            quantitative_covariate_file=str(toy["toy.qcovar"]),
            timeout_seconds=300,
            verbose=False,
        )
    )
    assert "CONVERGE" in binary["summary"]["columns"]


@requires_gcta
def test_binary_model_burden_and_acat(toy: dict[str, Path], sparse: str, tmp_path: Path) -> None:
    kwargs = {
        "bfile": _bfile(toy),
        "model": "mlm_binary",
        "sparse_grm_prefix": sparse,
        "quantitative_covariate_file": str(toy["toy.qcovar"]),
    }
    saved = _ok(
        g.gcta_fastgwa(
            str(toy["toy.bin.phen"]),
            str(tmp_path / "mdl"),
            model_only=True,
            timeout_seconds=300,
            verbose=False,
            **kwargs,
        )
    )
    assert set(saved["outputs"]) >= {"model", "model_ids"}
    burden = _ok(
        g.gcta_fastgwa_burden(
            _bfile(toy),
            saved["model_prefix"],
            str(toy["genes_burden.txt"]),
            str(tmp_path / "bur"),
            max_maf=0.5,
            maf=0.0001,
            verbose=False,
        )
    )
    assert burden["summary"]["rows"] == 3 and "VAR_N" in burden["summary"]["columns"]
    single = _ok(
        g.gcta_fastgwa(
            str(toy["toy.bin.phen"]),
            str(tmp_path / "single"),
            timeout_seconds=300,
            verbose=False,
            **kwargs,
        )
    )
    genes = tmp_path / "genes.txt"
    genes.write_text("1 1000 100000 G1\n1 200000 400000 G2\n2 1000 300000 G3\n")
    acat = _ok(
        g.gcta_acat(
            single["outputs"]["fastgwa"],
            str(genes),
            str(tmp_path / "acat.res"),
            max_maf=0.5,
            min_mac=5,
            verbose=False,
        )
    )
    rows = {row["GENE"]: row for row in acat["summary"]["top_rows"]}
    assert (
        rows["G1"]["SNP_NUM"] == 100
        and rows["G2"]["SNP_NUM"] == 201
        and rows["G3"]["SNP_NUM"] == 300
    )
    assert all(0 < row["P_ACAT"] <= 1 for row in rows.values())
    with pytest.raises(ConfigurationError):
        g.gcta_fastgwa_burden(
            _bfile(toy),
            str(tmp_path / "nomodel"),
            str(toy["genes_burden.txt"]),
            str(tmp_path / "b2"),
            verbose=False,
        )


@requires_gcta
def test_mlma_and_loco(toy: dict[str, Path], grm: str, tmp_path: Path) -> None:
    plain = _ok(
        g.gcta_mlma(
            _bfile(toy), str(toy["toy.phen"]), str(tmp_path / "m"), grm_prefix=grm, verbose=False
        )
    )
    loco = _ok(
        g.gcta_mlma(
            _bfile(toy), str(toy["toy.phen"]), str(tmp_path / "l"), loco=True, verbose=False
        )
    )
    assert plain["summary"]["rows"] == 2600 and Path(loco["outputs"]["mlma"]).name == "l.loco.mlma"
    assert plain["summary"]["columns"] == ["Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p"]
    assert plain["summary"]["top_rows"][0]["p"] == min(_p_values(plain["outputs"]["mlma"], "p"))
    with pytest.raises(ConfigurationError):
        g.gcta_mlma(
            _bfile(toy),
            str(toy["toy.phen"]),
            str(tmp_path / "x"),
            grm_prefix=grm,
            loco=True,
            verbose=False,
        )


@requires_gcta
def test_cojo_modes_and_sblup(toy: dict[str, Path], summary: Path, tmp_path: Path) -> None:
    select = _ok(
        g.gcta_cojo(
            _bfile(toy),
            str(summary),
            str(tmp_path / "s"),
            chromosome=1,
            p_threshold=1e-2,
            verbose=False,
        )
    )
    assert select["summary"]["rows"] >= 1 and {"jma", "ldr"} <= set(select["outputs"])
    top = _ok(
        g.gcta_cojo(
            _bfile(toy), str(summary), str(tmp_path / "t"), chromosome=1, top_snps=3, verbose=False
        )
    )
    assert top["summary"]["rows"] == 3
    joint = _ok(
        g.gcta_cojo(
            _bfile(toy),
            str(summary),
            str(tmp_path / "j"),
            mode="joint",
            extract=str(toy["joint.snplist"]),
            chromosome=1,
            verbose=False,
        )
    )
    assert joint["summary"]["rows"] == 3
    cond = _ok(
        g.gcta_cojo(
            _bfile(toy),
            str(summary),
            str(tmp_path / "c"),
            mode="cond",
            snp_list_file=str(toy["cond.snplist"]),
            chromosome=1,
            verbose=False,
        )
    )
    assert cond["summary"]["rows"] == 1299 and "given" in cond["outputs"]
    blup = _ok(
        g.gcta_sblup(
            _bfile(toy), str(summary), str(tmp_path / "sb"), 2000, chromosome=1, verbose=False
        )
    )
    assert blup["n_snps"] == 1300
    with pytest.raises(ConfigurationError):
        g.gcta_cojo(
            _bfile(toy), str(summary), str(tmp_path / "z"), top_snps=2, mode="joint", verbose=False
        )


@requires_gcta
def test_gene_based_tests_and_snp_sets(toy: dict[str, Path], summary: Path, tmp_path: Path) -> None:
    gene = _ok(
        g.gcta_fastbat(
            _bfile(toy),
            str(summary),
            str(tmp_path / "g"),
            gene_list_file=str(toy["genes.txt"]),
            verbose=False,
        )
    )
    assert (
        gene["summary"]["rows"] == 3 and Path(gene["outputs"]["fastbat"]).name == "g.gene.fastbat"
    )
    segment = _ok(
        g.gcta_fastbat(
            _bfile(toy), str(summary), str(tmp_path / "s"), segment_kb=100, verbose=False
        )
    )
    assert Path(segment["outputs"]["fastbat"]).name == "s.seg.fastbat"
    sets = _ok(
        g.gcta_fastbat(
            _bfile(toy),
            str(summary),
            str(tmp_path / "t"),
            set_list_file=str(toy["sets.txt"]),
            verbose=False,
        )
    )
    assert sets["summary"]["rows"] == 2
    combo = _ok(
        g.gcta_mbat_combo(
            _bfile(toy),
            str(summary),
            str(toy["genes.txt"]),
            str(tmp_path / "mb"),
            print_all_p=True,
            verbose=False,
        )
    )
    assert {"P_mBATcombo", "P_mBAT", "P_fastBAT"} <= set(combo["summary"]["columns"])
    made = _ok(
        g.gcta_write_gene_snpsets(
            _bfile(toy), str(summary), str(toy["genes.txt"]), str(tmp_path / "ss"), verbose=False
        )
    )
    assert made["genes"] == 3
    reuse = _ok(
        g.gcta_fastbat(
            _bfile(toy),
            str(summary),
            str(tmp_path / "again"),
            set_list_file=made["snpset"],
            verbose=False,
        )
    )
    assert reuse["summary"]["rows"] == 3
    mbat = _ok(
        g.gcta_write_gene_snpsets(
            _bfile(toy),
            str(summary),
            str(toy["genes.txt"]),
            str(tmp_path / "ms"),
            method="mbat",
            verbose=False,
        )
    )
    assert mbat["genes"] == 3 and mbat["snpset"].endswith(".gene.snpset.mbat")


@requires_gcta
def test_gsmr_and_mtcojo(toy: dict[str, Path], summary: Path, tmp_path: Path) -> None:
    exposure = tmp_path / "exposure.ma"
    lines = summary.read_text().splitlines()
    scaled = [lines[0]]
    for line in lines[1:]:
        parts = line.split()
        parts[4] = str(float(parts[4]) * 0.5 + 0.01)
        scaled.append(" ".join(parts))
    exposure.write_text("\n".join(scaled) + "\n")
    (tmp_path / "expo.txt").write_text(f"bmi {exposure}\n")
    (tmp_path / "outc.txt").write_text(f"t2d {summary}\n")
    gsmr = _ok(
        g.gcta_gsmr(
            str(tmp_path / "expo.txt"),
            str(tmp_path / "outc.txt"),
            str(tmp_path / "gs"),
            bfile=_bfile(toy),
            snp_min=3,
            gwas_thresh=0.05,
            clump_r2=0.5,
            heidi_thresh=0.0,
            effect_plot=True,
            verbose=False,
        )
    )
    row = gsmr["summary"]["top_rows"][0]
    assert (row["Exposure"], row["Outcome"]) == ("bmi", "t2d") and row["nsnp"] >= 3
    assert row["bxy"] > 0  # the outcome effects are the exposure effects scaled by 2
    assert "effect_plot" in gsmr["outputs"]
    ld = data.write_ldsc(tmp_path / "ldsc", Path(_bfile(toy) + ".bim"))
    (tmp_path / "mt.list").write_text(f"t2d {summary}\nbmi {exposure}\n")
    (tmp_path / "bxy.txt").write_text("bmi 0.5\n")
    cojo = _ok(
        g.gcta_mtcojo(
            str(tmp_path / "mt.list"),
            str(tmp_path / "mt"),
            str(ld),
            str(ld),
            bfile=_bfile(toy),
            bxy_file=str(tmp_path / "bxy.txt"),
            verbose=False,
        )
    )
    assert cojo["summary"]["rows"] == 2600 and "bC" in cojo["summary"]["columns"]


@requires_gcta
def test_ld_scores_friends_and_population_statistics(
    toy: dict[str, Path], grm: str, tmp_path: Path
) -> None:
    friends = _ok(
        g.gcta_ld_friends(
            _bfile(toy),
            str(toy["ld.snplist"]),
            str(tmp_path / "ldf"),
            window_kb=5000,
            significance=0.05,
            verbose=False,
        )
    )
    assert [row["target_SNP"] for row in friends["summary"]["first_rows"]] == ["rs5", "rs100"]
    assert {"rsq", "r", "snp"} <= set(friends["outputs"])
    segment = _ok(g.gcta_segment_ld_score(_bfile(toy), str(tmp_path / "seg"), verbose=False))
    assert segment["summary"]["rows"] == 2500 and "ldscore_SNP" in segment["summary"]["columns"]
    multi = _ok(
        g.gcta_multi_set_ld_score(
            _bfile(toy),
            [str(toy["set1.snps"]), str(toy["set2.snps"])],
            str(tmp_path / "multi"),
            window_kb=1000,
            verbose=False,
        )
    )
    assert multi["n_sets"] == 2 and multi["summary"]["rows"] == 2500
    fst = _ok(g.gcta_fst(_bfile(toy), str(toy["toy.subpop"]), str(tmp_path / "fst"), verbose=False))
    assert fst["summary"]["rows"] == 2600 and abs(fst["summary"]["mean_fst"]) < 0.02
    subprocess.run(
        ["gcta64", "--grm", grm, "--pca", "5", "--out", str(tmp_path / "pc")],
        check=True,
        capture_output=True,
    )
    loading = _ok(
        g.gcta_pc_loading(_bfile(toy), str(tmp_path / "pc"), str(tmp_path / "pl"), verbose=False)
    )
    assert loading["n_snps"] == 2500
    projected = _ok(
        g.gcta_project_loading(
            _bfile(toy), str(tmp_path / "pl"), str(tmp_path / "proj"), n_components=5, verbose=False
        )
    )
    assert projected["n_individuals"] == 200
    reference = [
        line.split()[2:]
        for line in Path(str(tmp_path / "pc") + ".eigenvec").read_text().splitlines()
    ]
    again = [
        line.split()[2:] for line in Path(projected["outputs"]["eigenvec"]).read_text().splitlines()
    ]
    corr = _correlation([float(row[0]) for row in reference], [float(row[0]) for row in again])
    assert abs(corr) > 0.99  # projecting the reference onto its own loadings reproduces its PCs


def _correlation(x: list[float], y: list[float]) -> float:
    mx, my = sum(x) / len(x), sum(y) / len(y)
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True))
    return cov / (sum((a - mx) ** 2 for a in x) ** 0.5 * sum((b - my) ** 2 for b in y) ** 0.5)


@requires_gcta
def test_simulated_phenotypes(toy: dict[str, Path], tmp_path: Path) -> None:
    quantitative = _ok(
        g.gcta_simulate_phenotypes(
            _bfile(toy),
            str(tmp_path / "q"),
            causal_loci_file=str(toy["causal.snplist"]),
            heritability=0.5,
            replicates=2,
            verbose=False,
        )
    )
    assert quantitative["n_individuals"] == 200
    assert {row["QTL"]: row["Effect"] for row in quantitative["causal_variants"]} == {
        "rs5": 0.1,
        "rs100": -0.2,
    }
    assert (
        len(Path(quantitative["outputs"]["phen"]).read_text().splitlines()[0].split()) == 4
    )  # FID IID + 2 replicates
    case_control = _ok(
        g.gcta_simulate_phenotypes(
            _bfile(toy),
            str(tmp_path / "cc"),
            trait="case_control",
            n_cases=40,
            n_controls=60,
            prevalence=0.3,
            causal_loci_file=str(toy["causal.snplist"]),
            verbose=False,
        )
    )
    values = [
        line.split()[2] for line in Path(case_control["outputs"]["phen"]).read_text().splitlines()
    ]
    assert values.count("2") == 40 and values.count("1") == 60
