"""
Orchestrate the ChIP-Atlas peak enrichment pipeline end to end.

Steps
-----
1. Submit the gene list to the WABI endpoint
2. Poll until completion
3. Retrieve and parse the TSV result
4. Apply optional client-side filters
5. (Optional) Ensembl coordinate verification — non-blocking
6. Return a results dict for `generate_all_plots` and `export_all`

Inputs : gene_list + API parameters + optional filters
Outputs: dict with `enrichment_results`, `input_genes`, `input_regions`,
         `failed_genes`, `metadata`, `api_total_regions`, `parameters`.

Example
-------
    python run_enrichment_workflow.py --example tp53_targets
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Iterable, List, Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")

try:
    from scripts.query_chip_atlas_enrichment import (
        SUPPORTED_GENOMES, VALID_ANTIGEN_CLASSES, VALID_THRESHOLDS,
        parse_api_results, poll_job_status, retrieve_results, submit_enrichment_job,
    )
    from scripts.filter_experiments import filter_experiments
    from scripts.convert_genes_to_regions import convert_genes_to_regions
    from scripts.load_example_data import load_example_data, list_examples
except ImportError:
    from query_chip_atlas_enrichment import (
        SUPPORTED_GENOMES, VALID_ANTIGEN_CLASSES, VALID_THRESHOLDS,
        parse_api_results, poll_job_status, retrieve_results, submit_enrichment_job,
    )
    from filter_experiments import filter_experiments
    from convert_genes_to_regions import convert_genes_to_regions
    from load_example_data import load_example_data, list_examples


def run_enrichment_workflow(
    gene_list: List[str],
    genome: str = "hg38",
    antigen_class: str = "TFs and others",
    cell_class: str = "All cell types",
    threshold: int = 50,
    distance_up: int = 5000,
    distance_down: int = 5000,
    antigen_filter: Optional[Iterable[str]] = None,
    cell_type_filter: Optional[Iterable[str]] = None,
    min_fold_enrichment: float = 1.0,
    max_qvalue: Optional[float] = None,
    output_dir: str = "enrichment_results",
    poll_timeout: int = 600,
    ensembl_retry: bool = True,
) -> dict:
    if len(gene_list) < 3:
        raise ValueError("gene_list must contain at least 3 symbols")
    if genome not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'")
    if antigen_class not in VALID_ANTIGEN_CLASSES:
        raise ValueError(
            f"antigen_class '{antigen_class}' not in {VALID_ANTIGEN_CLASSES}"
        )
    if threshold not in VALID_THRESHOLDS:
        raise ValueError(f"threshold must be in {sorted(VALID_THRESHOLDS)}")

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 72)
    print("ChIP-Atlas enrichment workflow")
    print("=" * 72)
    print(f"  genes:         {len(gene_list)}")
    print(f"  genome:        {genome}")
    print(f"  antigen_class: {antigen_class}")
    print(f"  cell_class:    {cell_class}")
    print(f"  threshold:     {threshold}")
    print(f"  TSS window:    -{distance_up} / +{distance_down}")

    request_id = submit_enrichment_job(
        gene_list=gene_list,
        genome=genome,
        antigen_class=antigen_class,
        cell_class=cell_class,
        threshold=threshold,
        distance_up=distance_up,
        distance_down=distance_down,
    )
    poll_job_status(request_id, timeout=poll_timeout)
    raw_frame = retrieve_results(request_id)
    results = parse_api_results(raw_frame)

    api_total_regions = int(results["total_regions"].iloc[0]) if len(results) else len(gene_list)
    if api_total_regions != len(gene_list):
        print(
            f"   ⚠ API analysed {api_total_regions} regions for {len(gene_list)} submitted genes."
        )

    if antigen_filter or cell_type_filter or min_fold_enrichment > 1.0 or max_qvalue is not None:
        print("Applying client-side filters…")
        results = filter_experiments(
            results,
            antigen_filter=antigen_filter,
            cell_type_filter=cell_type_filter,
            min_fold_enrichment=min_fold_enrichment,
            max_qvalue=max_qvalue,
        )

    print("Looking up gene coordinates via Ensembl (optional)…")
    input_regions: list = []
    failed_genes: List[str] = []
    try:
        input_regions, failed_genes = convert_genes_to_regions(
            gene_list=gene_list,
            genome=genome,
            upstream=distance_up,
            downstream=distance_down,
        )
    except Exception as exc:
        print(f"   ! Ensembl lookup failed ({exc}); enrichment results are unaffected")

    if ensembl_retry and not input_regions and len(gene_list) > 0:
        print("   no regions resolved — retrying once after 60 s")
        time.sleep(60)
        try:
            input_regions, failed_genes = convert_genes_to_regions(
                gene_list=gene_list,
                genome=genome,
                upstream=distance_up,
                downstream=distance_down,
            )
        except Exception as exc:
            print(f"   ! Ensembl retry failed ({exc})")

    if failed_genes:
        print(f"   ! {len(failed_genes)} gene(s) unresolved by Ensembl: {failed_genes}")

    n_significant = int((results["q_value"] < 0.05).sum()) if len(results) else 0
    print(f"\n✓ {len(results)} experiments returned; {n_significant} significant (q<0.05)")
    if len(results):
        top_row = results.iloc[0]
        print(
            f"   top: {top_row['antigen']} in {top_row['cell_type']} "
            f"(q={top_row['q_value']:.2e}, fold={top_row['fold_enrichment']:.1f}×)"
        )

    return {
        "enrichment_results": results,
        "input_genes": list(gene_list),
        "input_regions": input_regions,
        "failed_genes": failed_genes,
        "metadata": results.copy(),
        "api_total_regions": api_total_regions,
        "parameters": {
            "genome": genome,
            "antigen_class": antigen_class,
            "cell_class": cell_class,
            "threshold": threshold,
            "distance_up": distance_up,
            "distance_down": distance_down,
            "antigen_filter": list(antigen_filter) if antigen_filter else None,
            "cell_type_filter": list(cell_type_filter) if cell_type_filter else None,
        },
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--example", choices=list_examples(), help="Use a bundled gene set")
    parser.add_argument("--genes", nargs="+", help="Symbols to analyze")
    parser.add_argument("--genes-file", help="Plain text file with one gene per line")
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--antigen-class", default="TFs and others", choices=list(VALID_ANTIGEN_CLASSES))
    parser.add_argument("--cell-class", default="All cell types")
    parser.add_argument("--threshold", type=int, default=50, choices=sorted(VALID_THRESHOLDS))
    parser.add_argument("--distance-up", type=int, default=5000)
    parser.add_argument("--distance-down", type=int, default=5000)
    parser.add_argument("--output-dir", default="enrichment_results")
    parser.add_argument("--poll-timeout", type=int, default=600)
    args = parser.parse_args(argv)

    if args.example:
        payload = load_example_data(args.example)
        genes = payload["genes"]
    elif args.genes:
        genes = args.genes
    elif args.genes_file:
        with open(args.genes_file, "r", encoding="utf-8") as handle:
            genes = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    else:
        parser.error("Provide --example, --genes, or --genes-file")
        return 1

    results = run_enrichment_workflow(
        gene_list=genes,
        genome=args.genome,
        antigen_class=args.antigen_class,
        cell_class=args.cell_class,
        threshold=args.threshold,
        distance_up=args.distance_up,
        distance_down=args.distance_down,
        output_dir=args.output_dir,
        poll_timeout=args.poll_timeout,
    )
    print(f"✓ Workflow done; {len(results['enrichment_results'])} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
