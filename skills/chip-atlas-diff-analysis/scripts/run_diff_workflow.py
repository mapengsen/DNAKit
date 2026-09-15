"""
Orchestrate the full ChIP-Atlas diff-analysis pipeline.

Steps
-----
1. Submit groups A and B to the ChIP-Atlas WABI endpoint
2. Poll until the job completes (or fails / times out)
3. Download the result ZIP and extract
4. Parse the 8-column BED into a DataFrame with derived columns
5. Run automated QC checks on the unfiltered table
6. Apply default QC trims (non-standard contigs + <10 bp regions)
7. Apply optional user filters
8. Return a results dict ready for `generate_all_plots` and `export_all`

Inputs : two accession lists, parameters
Outputs: a dict with `diff_regions`, `diff_regions_unfiltered`,
         `qc_warnings`, `experiments_a/b`, `raw_files`, `log_content`,
         `parameters`.

Example
-------
    python run_diff_workflow.py --example tp53_k562_dmso_vs_daunorubicin
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

try:
    from scripts.query_chip_atlas_diff import (
        ANALYSIS_TYPES, poll_job_status, retrieve_zip_results, submit_diff_job,
    )
    from scripts.parse_diff_bed import parse_diff_bed, summarize_regions
    from scripts.filter_regions import (
        filter_min_region_size, filter_regions, filter_standard_chromosomes,
    )
    from scripts.qc_checks import run_qc_checks
    from scripts.load_example_data import load_example_data
except ImportError:
    from query_chip_atlas_diff import (
        ANALYSIS_TYPES, poll_job_status, retrieve_zip_results, submit_diff_job,
    )
    from parse_diff_bed import parse_diff_bed, summarize_regions
    from filter_regions import (
        filter_min_region_size, filter_regions, filter_standard_chromosomes,
    )
    from qc_checks import run_qc_checks
    from load_example_data import load_example_data


def run_diff_workflow(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str = "hg38",
    analysis_type: str = "diffbind",
    title: str = "chip_atlas_diff_run",
    description_a: str = "group_A",
    description_b: str = "group_B",
    min_score: int = 0,
    min_region_size: int = 0,
    direction_filter: Optional[str] = None,
    design_caveats: Optional[List[str]] = None,
    output_dir: str = "diff_results",
    poll_timeout: int = 900,
) -> dict:
    if len(experiments_a) < 2 or len(experiments_b) < 2:
        raise ValueError("Both groups need at least 2 experiment IDs")

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 72)
    print("ChIP-Atlas diff-analysis workflow")
    print("=" * 72)
    print(f"  analysis_type: {analysis_type} ({ANALYSIS_TYPES.get(analysis_type, '')})")
    print(f"  genome:        {genome}")
    print(f"  group A ({description_a}): {len(experiments_a)} experiments")
    print(f"  group B ({description_b}): {len(experiments_b)} experiments")

    request_id = submit_diff_job(
        experiments_a=experiments_a,
        experiments_b=experiments_b,
        genome=genome,
        analysis_type=analysis_type,
        title=title,
        description_a=description_a,
        description_b=description_b,
    )
    poll_job_status(request_id, timeout=poll_timeout)
    raw_files = retrieve_zip_results(request_id, output_dir)

    log_content = ""
    if "log" in raw_files:
        with open(raw_files["log"], "r", encoding="utf-8", errors="replace") as handle:
            log_content = handle.read()

    print("\nParsing BED…")
    parsed = parse_diff_bed(raw_files["bed"])
    summarize_regions(parsed)

    qc_warnings = run_qc_checks(
        parsed,
        n_samples_a=len(experiments_a),
        n_samples_b=len(experiments_b),
        analysis_type=analysis_type,
        genome=genome,
    )

    unfiltered = parsed.copy()
    print("\nApplying default QC filters…")
    filtered = filter_standard_chromosomes(parsed, genome=genome)
    filtered = filter_min_region_size(filtered, min_size=10)

    if min_score > 0 or min_region_size > 0 or direction_filter:
        print("Applying user filters…")
        filtered = filter_regions(
            filtered,
            min_score=min_score,
            min_region_size=min_region_size,
            direction=direction_filter,
        )

    summarize_regions(filtered)
    print("\n✓ Diff analysis completed successfully")

    return {
        "diff_regions": filtered,
        "diff_regions_unfiltered": unfiltered,
        "qc_warnings": qc_warnings,
        "experiments_a": experiments_a,
        "experiments_b": experiments_b,
        "raw_files": raw_files,
        "log_content": log_content,
        "parameters": {
            "genome": genome,
            "analysis_type": analysis_type,
            "title": title,
            "description_a": description_a,
            "description_b": description_b,
            "min_score": min_score,
            "min_region_size": min_region_size,
            "direction_filter": direction_filter,
            "design_caveats": list(design_caveats) if design_caveats else [],
        },
    }


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--example", help="Run with a bundled example dataset")
    parser.add_argument("--experiments-a", nargs="+")
    parser.add_argument("--experiments-b", nargs="+")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--analysis-type", default="diffbind", choices=list(ANALYSIS_TYPES))
    parser.add_argument("--description-a", default="group_A")
    parser.add_argument("--description-b", default="group_B")
    parser.add_argument("--output-dir", default="diff_results")
    parser.add_argument("--poll-timeout", type=int, default=900)
    args = parser.parse_args(argv)

    if args.example:
        payload = load_example_data(args.example)
        results = run_diff_workflow(
            experiments_a=payload["experiments_a"],
            experiments_b=payload["experiments_b"],
            genome=payload["genome"],
            analysis_type=payload["analysis_type"],
            description_a=payload["description_a"],
            description_b=payload["description_b"],
            design_caveats=payload.get("design_caveats"),
            output_dir=args.output_dir,
            poll_timeout=args.poll_timeout,
        )
    else:
        if not args.experiments_a or not args.experiments_b:
            parser.error("Provide --experiments-a and --experiments-b, or use --example")
        results = run_diff_workflow(
            experiments_a=args.experiments_a,
            experiments_b=args.experiments_b,
            genome=args.genome,
            analysis_type=args.analysis_type,
            description_a=args.description_a,
            description_b=args.description_b,
            output_dir=args.output_dir,
            poll_timeout=args.poll_timeout,
        )
    print(f"✓ Workflow done. {len(results['diff_regions'])} regions retained.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
