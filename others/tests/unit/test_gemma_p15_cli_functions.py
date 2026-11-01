from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.gwas import gemma_p15_cli_functions as gm

needs_gemma = pytest.mark.skipif(shutil.which("gemma") is None, reason="gemma is not installed")
N, P = 80, 60


def toy(folder: Path) -> dict[str, Path]:
    random.seed(7)
    freq = [random.uniform(0.2, 0.5) for _ in range(P)]
    geno = [[sum(random.random() < f for _ in range(2)) for _ in range(N)] for f in freq]
    (folder / "geno.txt").write_text("".join(f"rs{j + 1}, A, G, " + ", ".join(map(str, row)) + "\n" for j, row in enumerate(geno)))
    rows = []
    for i in range(N):
        y1 = 0.8 * geno[3][i] - 0.6 * geno[20][i] + random.gauss(0, 1)
        rows.append((y1, 0.5 * y1 + random.gauss(0, 1)))
    (folder / "pheno.txt").write_text("".join(f"{'NA' if i >= 70 else format(a, '.4f')}\t{b:.4f}\n" for i, (a, b) in enumerate(rows)))
    (folder / "env.txt").write_text("".join(f"{random.gauss(0, 1):.4f}\n" for _ in range(N)))
    code = {2: 0b00, 1: 0b10, 0: 0b11}
    with (folder / "toy.bed").open("wb") as bed, (folder / "toy.bim").open("w") as bim:
        bed.write(bytes([0x6C, 0x1B, 0x01]))
        for j, row in enumerate(geno):
            out = bytearray()
            for k in range(0, N, 4):
                out.append(sum(code[row[k + m]] << (2 * m) for m in range(4) if k + m < N))
            bed.write(out)
            bim.write(f"1\trs{j + 1}\t0\t{1000 * (j + 1)}\tA\tG\n")
    (folder / "toy.fam").write_text("".join(f"f{i} i{i} 0 0 1 -9\n" for i in range(N)))
    (folder / "beta.txt").write_text("SNP N Z INC_ALLELE DEC_ALLELE\n" + "".join(f"rs{j + 1} 70 {((j * 7) % 11 - 5) / 5.0} A G\n" for j in range(P)))
    return {name: folder / name for name in ("geno.txt", "pheno.txt", "env.txt", "beta.txt")}


def test_missing_program_and_validation(tmp_path: Path) -> None:
    files = toy(tmp_path)
    with pytest.raises(BackendUnavailableError):
        gm.gemma_eigen_decomposition(str(tmp_path / "o"), str(files["geno.txt"]), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), gemma_path=str(tmp_path / "missing"), verbose=False)
    with pytest.raises(ConfigurationError):
        gm.gemma_bslmm(str(tmp_path / "o1"), verbose=False)  # no genotypes
    with pytest.raises(ConfigurationError):
        gm.gemma_bslmm(str(tmp_path / "o2"), genotype_file=str(files["geno.txt"]), verbose=False)  # BIMBAM needs phenotypes
    with pytest.raises(ConfigurationError):
        gm.gemma_multivariate_association(str(tmp_path / "o3"), [1], genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), kinship_file=str(files["geno.txt"]), verbose=False)
    with pytest.raises(ConfigurationError):
        gm.gemma_predict_phenotypes(str(tmp_path / "o4"), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    with pytest.raises(ConfigurationError):
        gm.gemma_summary_variance_components(str(tmp_path / "o5"), str(files["beta.txt"]), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), method="mqs_ldw", verbose=False)


@needs_gemma
def test_kinship_workflows(tmp_path: Path) -> None:
    files = toy(tmp_path)
    import subprocess

    subprocess.run(["gemma", "-g", str(files["geno.txt"]), "-p", str(files["pheno.txt"]), "-gk", "1", "-outdir", str(tmp_path), "-o", "kin"], check=True, capture_output=True, cwd=tmp_path)
    kinship = tmp_path / "kin.cXX.txt"
    eigen = gm.gemma_eigen_decomposition(str(tmp_path / "eig"), str(kinship), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    assert eigen["individual_count"] == 70 and eigen["eigen_values"][0] >= 0
    mv = gm.gemma_multivariate_association(str(tmp_path / "mv"), [1, 2], genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), kinship_file=str(kinship), verbose=False)
    assert mv["association_count"] == P and "beta_2" in mv["associations"][0]
    bs = gm.gemma_bslmm(str(tmp_path / "bs"), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), burn_in=100, sampling_steps=500, seed=1, verbose=False)
    assert bs["snp_count"] == P and 0 <= bs["hyper_parameter_means"]["pve"] <= 1
    ridge = gm.gemma_bslmm(str(tmp_path / "ridge"), model="ridge", genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    assert ridge["snp_count"] == P
    predicted = gm.gemma_predict_phenotypes(
        str(tmp_path / "pr"), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]),
        parameter_file=str(tmp_path / "bs" / "result.param.txt"), log_file=str(tmp_path / "bs" / "result.log.txt"), verbose=False,
    )
    assert predicted["predicted_count"] == 10 and predicted["predictions"][0] is None
    blup = gm.gemma_predict_phenotypes(str(tmp_path / "blup"), source="kinship", genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), kinship_file=str(kinship), phenotype_columns=[1, 2], verbose=False)
    assert len(blup["predictions"]) == N
    gxe = gm.gemma_gxe_association(str(tmp_path / "gxe"), str(tmp_path / "toy"), str(files["env.txt"]), kinship_file=str(kinship), phenotype_file=str(files["pheno.txt"]), verbose=False)
    assert gxe["association_count"] == P


@needs_gemma
def test_summary_statistic_functions(tmp_path: Path) -> None:
    files = toy(tmp_path)
    s = gm.gemma_snp_covariance_matrix(str(tmp_path / "s"), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    assert s["sizes"][0] == P and s["status"] == "completed"
    q = gm.gemma_summary_statistic_q_vector(str(tmp_path / "q"), str(files["beta.txt"]), str(files["pheno.txt"]), verbose=False)
    assert len(q["q"]) == 1
    vc = gm.gemma_summary_variance_components(str(tmp_path / "vc"), str(files["beta.txt"]), genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    pve = vc["log_values"]["pve_estimates"]
    ci = gm.gemma_summary_variance_component_errors(str(tmp_path / "ci"), str(files["beta.txt"]), str(tmp_path / "vc" / "result"), pve, genotype_file=str(files["geno.txt"]), phenotype_file=str(files["pheno.txt"]), verbose=False)
    assert ci["log_values"]["se_pve"] > 0
