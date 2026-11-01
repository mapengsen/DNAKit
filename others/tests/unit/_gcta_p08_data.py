"""Small deterministic PLINK data and phenotypes for the GCTA P08 facade tests."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

SAMPLES = 200
CHR1, CHR2, CHRX = 1300, 1200, 100


def gcta_available() -> bool:
    return shutil.which("gcta64") is not None or shutil.which("gcta") is not None


def write_plink(
    prefix: Path, genotypes: np.ndarray, chromosomes: list[int], positions: list[int]
) -> None:
    """Write a PLINK .bed/.bim/.fam fileset from a samples x SNPs matrix of allele counts."""
    samples, snps = genotypes.shape
    code = {0: 0b00, 1: 0b10, 2: 0b11}
    with open(f"{prefix}.bed", "wb") as handle:
        handle.write(bytes([0x6C, 0x1B, 0x01]))
        for j in range(snps):
            packed = bytearray((samples + 3) // 4)
            for i, value in enumerate(genotypes[:, j]):
                packed[i // 4] |= code[int(value)] << (2 * (i % 4))
            handle.write(bytes(packed))
    Path(f"{prefix}.bim").write_text(
        "".join(f"{chromosomes[j]}\trs{j + 1}\t0\t{positions[j]}\tA\tG\n" for j in range(snps)),
        encoding="utf-8",
    )
    Path(f"{prefix}.fam").write_text(
        "".join(f"F{i + 1} I{i + 1} 0 0 {1 + i % 2} -9\n" for i in range(samples)), encoding="utf-8"
    )


def write_toy(directory: Path, seed: int = 1) -> dict[str, Path]:
    """Write ``toy`` (200 samples, 2600 SNPs on chromosomes 1, 2 and X) with its side files."""
    rng = np.random.default_rng(seed)
    snps = CHR1 + CHR2 + CHRX
    frequency = rng.uniform(0.1, 0.5, snps)
    genotypes = rng.binomial(2, frequency, size=(SAMPLES, snps))
    chromosomes = [1] * CHR1 + [2] * CHR2 + [23] * CHRX
    positions = (
        [1000 * (k + 1) for k in range(CHR1)]
        + [1000 * (k + 1) for k in range(CHR2)]
        + [1000 * (k + 1) for k in range(CHRX)]
    )
    prefix = directory / "toy"
    write_plink(prefix, genotypes, chromosomes, positions)
    beta = np.zeros(snps)
    causal = rng.choice(snps, 60, replace=False)
    beta[causal] = rng.normal(0, 0.3, 60)
    standardised = (genotypes - 2 * frequency) / np.sqrt(2 * frequency * (1 - frequency))
    trait = standardised @ beta + rng.normal(0, 1.0, SAMPLES)
    second = 0.5 * trait + rng.normal(0, 1.0, SAMPLES)
    binary = (trait > np.median(trait)).astype(int)
    covariate = rng.normal(0, 1.0, SAMPLES)
    lines = {
        "toy.phen": [f"F{i + 1} I{i + 1} {trait[i]:.5f}" for i in range(SAMPLES)],
        "toy2.phen": [f"F{i + 1} I{i + 1} {trait[i]:.5f} {second[i]:.5f}" for i in range(SAMPLES)],
        "toy.bin.phen": [f"F{i + 1} I{i + 1} {binary[i]}" for i in range(SAMPLES)],
        "toy.qcovar": [f"F{i + 1} I{i + 1} {covariate[i]:.4f}" for i in range(SAMPLES)],
        "toy.envir": [
            f"F{i + 1} I{i + 1} {covariate[(i + 7) % SAMPLES]:.4f}" for i in range(SAMPLES)
        ],
        "toy.env_class": [f"F{i + 1} I{i + 1} {i % 3}" for i in range(SAMPLES)],
        "toy.subpop": [
            f"F{i + 1} I{i + 1} {'PopA' if i % 2 == 0 else 'PopB'}" for i in range(SAMPLES)
        ],
        "keep150.txt": [f"F{i + 1} I{i + 1}" for i in range(150)],
        "causal.snplist": ["rs5 0.1", "rs100 -0.2"],
        "ld.snplist": ["rs5", "rs100"],
        "joint.snplist": ["rs5", "rs100", "rs300"],
        "cond.snplist": ["rs5"],
        "genes.txt": ["1 1 400000 Gene1", "1 400001 800000 Gene2", "2 1 400000 Gene3"],
        "genes_burden.txt": [
            "gene chr start end strand",
            "G1 1 1 100000 +",
            "G2 1 200000 400000 +",
            "G3 2 1 300000 -",
        ],
        "sets.txt": ["S1", "rs1", "rs2", "rs3", "rs4", "END", "S2", "rs10", "rs11", "rs12", "END"],
        "set1.snps": [f"rs{k}" for k in range(1, 301)],
        "set2.snps": [f"rs{k}" for k in range(301, 701)],
    }
    for name, content in lines.items():
        (directory / name).write_text("\n".join(content) + "\n", encoding="utf-8")
    return {"bfile": prefix, **{name: directory / name for name in lines}}


def write_ldsc(directory: Path, bim: Path) -> Path:
    """LDSC-format LD score files for chromosomes 1 to 22 (constant score) for ``bim``."""
    import gzip

    directory.mkdir()
    rows = [line.split() for line in bim.read_text(encoding="utf-8").splitlines()]
    for chromosome in range(1, 23):
        selected = [row for row in rows if row[0] == str(chromosome)]
        with gzip.open(directory / f"{chromosome}.l2.ldscore.gz", "wt") as handle:
            handle.write("CHR\tSNP\tBP\tL2\n")
            for row in selected:
                handle.write(f"{row[0]}\t{row[1]}\t{row[3]}\t1.5\n")
        (directory / f"{chromosome}.l2.M_5_50").write_text(f"{len(selected)}\n", encoding="utf-8")
        (directory / f"{chromosome}.l2.M").write_text(f"{len(selected)}\n", encoding="utf-8")
    return directory


def ma_from_mlma(mlma: Path, destination: Path, sample_size: int = SAMPLES) -> Path:
    """A COJO ``.ma`` summary file (SNP A1 A2 freq b se p N) from an MLMA result table."""
    lines = ["SNP A1 A2 freq b se p N"]
    for row in mlma.read_text(encoding="utf-8").splitlines()[1:]:
        parts = row.split()
        lines.append(
            " ".join(
                [
                    parts[1],
                    parts[3],
                    parts[4],
                    parts[5],
                    parts[6],
                    parts[7],
                    parts[8],
                    str(sample_size),
                ]
            )
        )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination
