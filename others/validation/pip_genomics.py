"""Validate pip-installed genomics APIs with deterministic synthetic inputs.

This version verifies PCA, ANI, CDS prediction and BAM statistics without
launching external processes, downloading data or requiring user input files.
Run after installing the current checkout with ``pip install -e '.[genomics]'``.
"""

from __future__ import annotations

import argparse
import json
import random
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any
from unittest.mock import patch

from dnakit._progress import Progress
from dnakit.bioinformatics.reads import alignment_summary
from dnakit.genomics import (
    average_nucleotide_identity,
    genotype_pca,
    predict_prokaryotic_genes,
)


def validate(*, quiet: bool = False) -> dict[str, Any]:
    """Return evidence from real Python backends; failures raise an exception."""
    checks: dict[str, Any] = {}
    with (
        patch("subprocess.Popen", side_effect=RuntimeError("External process is forbidden")),
        Progress(disable=quiet) as progress,
        TemporaryDirectory(prefix="dnakit-pip-genomics-") as directory,
    ):
        task = progress.add_task("Validate pip genomics", total=4)
        pca = genotype_pca([[0, 0], [1, 2], [2, 1]], components=2, verbose=False)
        if any(
            abs(actual - expected) > 1e-10
            for actual, expected in zip(pca["eigenvalues"], [3.0, 1.0], strict=True)
        ):
            raise RuntimeError("PCA differs from the analytical covariance eigenvalues.")
        checks["pca"] = {"eigenvalues": pca["eigenvalues"], "status": "passed"}
        progress.advance(task)

        rng = random.Random(42)
        genome = "".join(rng.choices("ACGT", k=60_000))
        ani = average_nucleotide_identity(genome, genome, verbose=False)
        if ani["status"] != "matched" or abs(ani["ani_percent"] - 100.0) > 0.01:
            raise RuntimeError("Identical synthetic genomes did not produce near-100% ANI.")
        checks["ani"] = {
            "ani_percent": ani["ani_percent"],
            "query_fragments": ani["query_fragments"],
            "status": "passed",
        }
        progress.advance(task)

        cds = "ATG" + "GCTGACGAACTGAAAGGT" * 20 + "TAA"
        sequence = "TAA" * 10 + "AGGAGG" + "A" * 7 + cds + "TAA" * 10
        genes = predict_prokaryotic_genes(sequence, closed=True, verbose=False)
        if len(genes.genes) != 1 or genes.genes[0].cds != cds:
            raise RuntimeError("Pyrodigal did not recover the known synthetic CDS.")
        checks["cds_prediction"] = {"genes": len(genes.genes), "status": "passed"}
        progress.advance(task)

        import pysam

        path = Path(directory) / "synthetic.bam"
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "synthetic", "LN": 100}]}
        with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
            read = pysam.AlignedSegment()
            read.query_name = "synthetic_read"
            read.query_sequence = "ACGTACGT"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigarstring = "8M"
            read.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
            handle.write(read)
        summary = alignment_summary(path, verbose=False)
        if summary["records"] != 1 or summary["mapped"] != 1:
            raise RuntimeError("BAM statistics differ from the constructed alignment.")
        checks["bam_statistics"] = {"records": summary["records"], "status": "passed"}
        progress.advance(task)

    return {
        "status": "passed",
        "data": "deterministic synthetic fixtures; no biological accuracy claim",
        "subprocess_popen": "blocked during all four checks; in-process APIs inspected",
        "versions": {name: version(name) for name in ("numpy", "pyfastani", "pyrodigal", "pysam")},
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write evidence to a new JSON file.")
    parser.add_argument("--quiet", action="store_true", help="Disable the progress bar.")
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        parser.error("--output must be a new file")
    result = validate(quiet=args.quiet)
    text = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as handle:
            handle.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
