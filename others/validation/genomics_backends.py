#!/usr/bin/env python3
"""Validate four PLINK adapters on synthetic data using an installed PLINK 2.

This version creates deterministic synthetic genotypes, runs PCA/KING/LD/GWAS,
checks numerical/file contracts, and records backend version and provenance.
It does not download software or establish biological prediction accuracy.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

from rich.progress import track

from dnakit.integrations import run_analysis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plink2", default="plink2", help="Installed PLINK 2 executable")
    parser.add_argument("--output", type=Path, required=True, help="New validation directory")
    args = parser.parse_args()
    binary = shutil.which(args.plink2)
    if binary is None:
        parser.error("PLINK 2 is not installed or executable.")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    version = subprocess.check_output([binary, "--version"], text=True, timeout=30).strip()
    log = subprocess.run(
        [
            binary,
            "--dummy",
            "100",
            "1000",
            "--seed",
            "42",
            "--make-bed",
            "--out",
            str(root / "toy"),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    (root / "fixture.log").write_text(log.stdout + log.stderr, encoding="utf-8")
    rows = [x.split() for x in (root / "toy.fam").read_text(encoding="utf-8").splitlines()]
    for name, column, values in (
        ("pheno.tsv", "trait", [(i % 11) + i / 100 for i in range(100)]),
        ("covar.tsv", "age", [20 + (i * 7) % 53 for i in range(100)]),
    ):
        body = "#FID\tIID\t" + column + "\n"
        body += "".join(
            f"{row[0]}\t{row[1]}\t{value}\n" for row, value in zip(rows, values, strict=True)
        )
        (root / name).write_text(body, encoding="utf-8")
    results = []
    for analysis in track(
        ["plink2-pca", "plink2-king", "plink2-ld-prune", "plink2-gwas"],
        description="Validate real PLINK adapters",
    ):
        inputs = {"genotypes": str(root / "toy")}
        if analysis == "plink2-gwas":
            inputs.update(phenotypes=str(root / "pheno.tsv"), covariates=str(root / "covar.tsv"))
        result = run_analysis(
            analysis,
            inputs,
            str(root / analysis),
            executable_paths={"plink2": binary},
            timeout_seconds=120,
            verbose=False,
        )
        results.append(
            {
                "analysis": analysis,
                "status": result["status"],
                "artifacts": sorted(result["artifacts"]),
            }
        )
    pca = (root / "plink2-pca/pca.eigenvec").read_text(encoding="utf-8").splitlines()
    if len(pca) != 101:
        raise AssertionError("PCA must report one row per synthetic sample.")
    if not all(math.isfinite(float(x)) for row in pca[1:] for x in row.split()[2:]):
        raise AssertionError("PCA values must be finite.")
    sample_ids = {(row[0], row[1]) for row in rows}
    if {(row.split()[0], row.split()[1]) for row in pca[1:]} != sample_ids:
        raise AssertionError("PCA sample identifiers differ from the input.")
    kinship_path = root / "plink2-king/kinship.kin0"
    pairs = set()
    with kinship_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pair = frozenset(((row["#FID1"], row["IID1"]), (row["FID2"], row["IID2"])))
            if len(pair) != 2 or pair in pairs or not pair <= sample_ids:
                raise AssertionError("KING must report unique pairs of input samples.")
            if not math.isfinite(float(row["KINSHIP"])):
                raise AssertionError("KING coefficients must be finite.")
            pairs.add(pair)
    if len(pairs) != 100 * 99 // 2:
        raise AssertionError("KING must report every pair of synthetic samples.")
    kept = set((root / "plink2-ld-prune/ld.prune.in").read_text().split())
    removed = set((root / "plink2-ld-prune/ld.prune.out").read_text().split())
    variant_ids = {r.split()[1] for r in (root / "toy.bim").read_text().splitlines()}
    if kept & removed or kept | removed != variant_ids:
        raise AssertionError("LD outputs must partition all input variants.")
    gwas = list((root / "plink2-gwas").glob("gwas.*.glm.*"))
    if [p.name for p in gwas] != ["gwas.trait.glm.linear"]:
        raise AssertionError("Only explicitly supplied phenotypes may be analyzed.")
    gwas_rows = gwas[0].read_text(encoding="utf-8").splitlines()
    if len(gwas_rows) != 1001:
        raise AssertionError("GWAS must report every synthetic variant.")
    tested_variants = 0
    constant_variants = 0
    for row in csv.DictReader(gwas_rows, delimiter="\t"):
        if row["ERRCODE"] == "CONST_OMITTED_ALLELE":
            if row["P"] != "NA" or float(row["A1_FREQ"]) not in {0.0, 1.0}:
                raise AssertionError("Constant variants must be explicitly untestable.")
            constant_variants += 1
            continue
        probability = float(row["P"])
        if not math.isfinite(probability) or not 0 <= probability <= 1 or row["ERRCODE"] != ".":
            raise AssertionError("Synthetic GWAS must return finite probabilities without errors.")
        tested_variants += 1
    if not tested_variants or tested_variants + constant_variants != len(variant_ids):
        raise AssertionError("GWAS must account for tested and constant variants.")
    report = {
        "snapshot_date": "2026-09-20",
        "backend_version": version,
        "binary_sha256": hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        "fixture": "PLINK --dummy 100 1000 --seed 42; synthetic quantitative trait and covariate",
        "scientific_accuracy": "not assessed; no real biological dataset",
        "checks": [
            "100 PCA rows",
            "Finite PCA coordinates with original sample IDs",
            "4950 unique KING sample pairs with finite coefficients",
            "LD outputs partition 1000 variants",
            "GWAS reports only supplied trait with 1000 rows",
            "Finite GWAS probabilities in [0, 1]; constant variants explicitly untestable",
        ],
        "gwas_tested_variants": tested_variants,
        "gwas_constant_variants": constant_variants,
        "results": results,
    }
    target = root / "validation.json"
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
