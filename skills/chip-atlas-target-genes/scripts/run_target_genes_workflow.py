"""
Orchestrate the ChIP-Atlas target-genes pipeline end to end.

Steps
-----
1. Validate antigen availability (HTTP HEAD)
2. Download and parse the wide-format TSV
3. Apply filters (score, cell types, STRING, binding rate, top-N)
4. Detect co-located genes (genes with identical score vectors)
5. Return a results dict for `generate_all_plots` and `export_all`

Inputs : protein name + parameters
Outputs: dict with `target_genes`, `experiment_data`, `cell_types`,
         `protein`, `parameters`, `metadata`.

Example
-------
    python run_target_genes_workflow.py --example tp53
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, List, Optional

try:
    from scripts.download_target_genes import (
        SUPPORTED_GENOMES, VALID_DISTANCES,
        check_antigen_available, download_target_genes, get_cell_types,
    )
    from scripts.filter_targets import detect_colocated_genes, filter_targets
    from scripts.load_example_query import list_examples, load_example_query
except ImportError:
    from download_target_genes import (
        SUPPORTED_GENOMES, VALID_DISTANCES,
        check_antigen_available, download_target_genes, get_cell_types,
    )
    from filter_targets import detect_colocated_genes, filter_targets
    from load_example_query import list_examples, load_example_query


def run_target_genes_workflow(
    protein: str,
    genome: str = "hg38",
    distance: int = 5,
    min_score: float = 0.0,
    top_n: Optional[int] = 500,
    cell_types: Optional[Iterable[str]] = None,
    min_string_score: float = 0.0,
    min_binding_rate: float = 0.0,
    output_dir: str = "target_genes_results",
) -> dict:
    if not protein.strip():
        raise ValueError("Protein name must be non-empty (e.g. 'TP53', 'CTCF')")
    if genome not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'")
    if distance not in VALID_DISTANCES:
        raise ValueError(f"Unsupported distance {distance}")

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 72)
    print("ChIP-Atlas target-genes workflow")
    print("=" * 72)
    print(f"  protein:  {protein}")
    print(f"  genome:   {genome}")
    print(f"  distance: ±{distance} kb")

    print("\nStep 1 — validate availability")
    if not check_antigen_available(protein, genome, distance):
        raise ValueError(
            f"No target-gene data available for '{protein}' ({genome}, ±{distance}kb). "
            "Check the spelling (case-sensitive) and that this is a non-histone antigen."
        )

    print("\nStep 2 — download and parse")
    summary_df, experiment_df = download_target_genes(protein, genome, distance)
    total_genes = len(summary_df)
    total_experiments = len(experiment_df.columns) - 1
    all_cell_types = get_cell_types(experiment_df)

    print("\nStep 3 — apply filters")
    has_filters = (
        min_score > 0
        or cell_types is not None
        or min_string_score > 0
        or min_binding_rate > 0
        or top_n is not None
    )
    if has_filters:
        summary_df, experiment_df = filter_targets(
            summary_df,
            experiment_df=experiment_df,
            min_avg_score=min_score,
            cell_types=cell_types,
            min_string_score=min_string_score,
            top_n=top_n,
            min_binding_rate=min_binding_rate,
        )
    else:
        print("   no filters applied")

    summary_df, colocated_info = detect_colocated_genes(summary_df, experiment_df)

    parameters = {
        "protein": protein,
        "genome": genome,
        "distance_kb": distance,
        "min_score": min_score,
        "top_n": top_n,
        "cell_types": list(cell_types) if cell_types else None,
        "min_string_score": min_string_score,
        "min_binding_rate": min_binding_rate,
        "output_dir": output_dir,
    }

    metadata = {
        "total_genes_before_filter": total_genes,
        "total_genes_after_filter": len(summary_df),
        "total_experiments": total_experiments,
        "genome": genome,
        "distance_kb": distance,
        "unique_cell_types": len(all_cell_types),
        "colocated_info": colocated_info,
    }

    print("\n" + "=" * 72)
    print(f"  ✓ {len(summary_df)} target genes retained")
    print(f"  ✓ {total_experiments} experiments parsed")
    print(f"  ✓ {len(all_cell_types)} unique cell types")
    print("=" * 72)

    return {
        "target_genes": summary_df,
        "experiment_data": experiment_df,
        "cell_types": all_cell_types,
        "protein": protein,
        "parameters": parameters,
        "metadata": metadata,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--example", choices=list_examples(), help="Run a bundled query")
    parser.add_argument("--protein", help="Protein/TF name (case-sensitive)")
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--distance", type=int, default=5, choices=VALID_DISTANCES)
    parser.add_argument("--min-score", type=float, default=0.0)
    parser.add_argument("--top-n", type=int, default=500)
    parser.add_argument("--cell-types", nargs="*", default=None)
    parser.add_argument("--min-string-score", type=float, default=0.0)
    parser.add_argument("--min-binding-rate", type=float, default=0.0)
    parser.add_argument("--output-dir", default="target_genes_results")
    args = parser.parse_args(argv)

    if args.example:
        query = load_example_query(args.example)
        protein, genome, distance = query["protein"], query["genome"], query["distance"]
    else:
        if not args.protein:
            parser.error("Provide --protein or --example")
        protein, genome, distance = args.protein, args.genome, args.distance

    results = run_target_genes_workflow(
        protein=protein,
        genome=genome,
        distance=distance,
        min_score=args.min_score,
        top_n=args.top_n,
        cell_types=args.cell_types,
        min_string_score=args.min_string_score,
        min_binding_rate=args.min_binding_rate,
        output_dir=args.output_dir,
    )
    print(f"✓ Workflow done; {len(results['target_genes'])} retained target genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
