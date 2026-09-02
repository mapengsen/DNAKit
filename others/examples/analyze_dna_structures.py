"""Analyze the downloaded DNA structure fixtures and save auditable JSON results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskID, TextColumn

from dnakit.core import DNASequence
from dnakit.secondary_structure import analyze_dot_bracket, probe_nupack
from dnakit.structure3d import (
    analyze_ensemble_flexibility,
    analyze_structure,
    load_pdb_ensemble,
    read_dssr_json,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze local DNA PDB ensembles and documented secondary structures."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("temp/dna_structures"),
        help="Directory containing manifest.json and downloaded fixtures.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("temp/dna_structures/analysis_results.json"),
        help="JSON result path.",
    )
    parser.add_argument(
        "--sasa-points",
        type=int,
        default=96,
        help="Fibonacci-sphere points per atom for native SASA.",
    )
    parser.add_argument(
        "--volume-grid",
        type=float,
        default=0.75,
        help="Voxel spacing in angstrom for native molecular volume.",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _manifest(input_dir: Path) -> tuple[dict[str, Any], dict[str, str]]:
    manifest_path = input_dir / "manifest.json"
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or not isinstance(raw.get("files"), list):
        raise ValueError("manifest.json must contain a files array")
    expected = {
        str(item["filename"]): str(item["sha256"])
        for item in raw["files"]
        if isinstance(item, dict) and "filename" in item and "sha256" in item
    }
    return raw, expected


def _secondary_records(path: Path) -> list[dict[str, Any]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    records: list[dict[str, Any]] = []
    for index in range(0, len(lines), 3):
        if index + 2 >= len(lines) or not lines[index].startswith(">"):
            raise ValueError(f"invalid three-line DBN record at line {index + 1}")
        summary = analyze_dot_bracket(
            (DNASequence(lines[index + 1]),),
            lines[index + 2],
        )
        records.append({"name": lines[index][1:], "summary": summary.to_dict()})
    return records


def main() -> int:
    args = _arguments()
    input_dir = args.input_dir.resolve()
    output = args.output.resolve()
    manifest, expected_hashes = _manifest(input_dir)
    pdb_paths = sorted(input_dir.glob("*.pdb"))
    if not pdb_paths:
        raise FileNotFoundError(f"no PDB fixtures found under {input_dir}")

    structure_results: list[dict[str, Any]] = []
    progress = Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        console=Console(),
    )
    with progress:
        for path in pdb_paths:
            observed_hash = _sha256(path)
            expected_hash = expected_hashes.get(path.name)
            if expected_hash is not None and observed_hash != expected_hash:
                raise ValueError(f"SHA-256 mismatch for {path.name}")
            models = load_pdb_ensemble(path)
            total_atoms = len(models[0].atoms)
            task = progress.add_task(path.name, total=total_atoms * 2)

            def update(
                stage: str,
                completed: int,
                total: int,
                task_id: TaskID = task,
            ) -> None:
                offset = 0 if stage == "sasa" else total
                progress.update(task_id, completed=offset + completed)

            geometry = analyze_structure(
                models[0],
                sasa_points_per_atom=args.sasa_points,
                volume_grid_spacing_angstrom=args.volume_grid,
                progress=update,
            )
            flexibility = None
            if len(models) > 1:
                flexibility = analyze_ensemble_flexibility(models).to_dict()
            structure_results.append(
                {
                    "filename": path.name,
                    "sha256": observed_hash,
                    "checksum_verified": expected_hash == observed_hash,
                    "model_count": len(models),
                    "geometry": geometry.to_dict(),
                    "ensemble_flexibility": flexibility,
                }
            )

    dbn_path = input_dir / "1AC7_secondary_structure.dbn"
    dssr_path = input_dir / "1ehz-dssr-example.json"
    nupack = probe_nupack()
    result = {
        "manifest": manifest,
        "structures": structure_results,
        "secondary_structure_annotations": (
            _secondary_records(dbn_path) if dbn_path.is_file() else []
        ),
        "dssr_schema_example": read_dssr_json(dssr_path).to_dict() if dssr_path.is_file() else None,
        "nupack_backend": nupack.to_dict(),
        "methods": {
            "sasa_points_per_atom": args.sasa_points,
            "volume_grid_spacing_angstrom": args.volume_grid,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    Console().print(f"[green]已写入[/green] {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
