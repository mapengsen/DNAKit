"""
Convert gene symbols to promoter regions via the Ensembl REST API.

Output is a list of (chrom, start, end, gene, strand) tuples. Genes that
cannot be resolved are reported in the second return value so the caller
can surface them in the summary.

Behaviour
---------
* Tries a single batch POST (`/lookup/symbol/{species}`) first.
* Falls back to per-symbol GETs for genes the batch missed.
* Retries each per-symbol GET with exponential backoff (max 3 attempts).
* Applies a small alias map for known retired HGNC symbols.

Inputs : gene list, genome assembly, TSS window
Outputs: (list of regions, list of failed symbols)

Example
-------
    python convert_genes_to_regions.py --genes TP53 MYC CDKN1A INVALID
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import List, Optional, Tuple

try:
    import requests
except ImportError:
    sys.exit("Install with: pip install requests")


ENSEMBL_REST = "https://rest.ensembl.org"

GENOME_TO_SPECIES = {
    "hg": "human",
    "mm": "mouse",
    "rn": "rat",
    "dm": "drosophila_melanogaster",
    "ce": "caenorhabditis_elegans",
    "sacCer": "saccharomyces_cerevisiae",
}

# Retired -> current HGNC symbols. Add to this dict when you spot a new
# pair in the wild; the enrichment results are still correct without
# Ensembl lookup, but the input_regions table is nicer with it.
GENE_ALIASES = {
    "IL8": "CXCL8", "IL8A": "CXCL8",
    "GRO1": "CXCL1", "GRO2": "CXCL2", "GRO3": "CXCL3",
    "SDF1": "CXCL12", "MIG": "CXCL9", "IP10": "CXCL10",
    "MCP1": "CCL2", "MIP1A": "CCL3", "MIP1B": "CCL4", "RANTES": "CCL5",
    "TNFSF6": "FASLG", "APO1": "FAS", "TRAIL": "TNFSF10",
    "MLL": "KMT2A", "MLL2": "KMT2D", "MLL3": "KMT2C",
    "P53": "TP53", "ERK1": "MAPK3", "ERK2": "MAPK1", "P38": "MAPK14",
}

Region = Tuple[str, int, int, str, str]


def convert_genes_to_regions(
    gene_list: List[str],
    genome: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
) -> Tuple[List[Region], List[str]]:
    species = _resolve_species(genome)
    regions, missing = _batch_lookup(gene_list, species, upstream, downstream)

    if not missing:
        return regions, missing

    failed: List[str] = []
    for gene in missing:
        region = _sequential_lookup(gene, species, upstream, downstream)
        if region is not None:
            regions.append(region)
            time.sleep(0.5)
            continue
        alias = GENE_ALIASES.get(gene.upper())
        if alias:
            print(f"   retrying {gene} as alias {alias}")
            region = _sequential_lookup(alias, species, upstream, downstream)
            if region is not None:
                chrom, start, end, _, strand = region
                regions.append((chrom, start, end, gene, strand))
                time.sleep(0.5)
                continue
        failed.append(gene)
        time.sleep(0.5)
    return regions, failed


def _resolve_species(genome: str) -> str:
    for prefix, species in GENOME_TO_SPECIES.items():
        if genome.startswith(prefix):
            return species
    return "human"


def _batch_lookup(
    gene_list: List[str], species: str, upstream: int, downstream: int,
) -> Tuple[List[Region], List[str]]:
    url = f"{ENSEMBL_REST}/lookup/symbol/{species}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    regions: List[Region] = []
    missing: List[str] = []
    try:
        response = requests.post(url, json={"symbols": gene_list}, headers=headers, timeout=30)
        response.raise_for_status()
        payload = response.json()
    except Exception as exc:
        print(f"   ! batch Ensembl lookup failed ({exc}); falling back to per-symbol")
        return [], list(gene_list)

    for gene in gene_list:
        record = payload.get(gene)
        if record is None:
            missing.append(gene)
            continue
        region = _record_to_region(record, gene, upstream, downstream)
        if region is None:
            missing.append(gene)
            continue
        regions.append(region)
    print(f"   batch lookup: {len(regions)}/{len(gene_list)} resolved")
    return regions, missing


def _sequential_lookup(
    gene: str, species: str, upstream: int, downstream: int,
) -> Optional[Region]:
    url = f"{ENSEMBL_REST}/lookup/symbol/{species}/{gene}"
    headers = {"Content-Type": "application/json"}
    for attempt in range(3):
        try:
            response = requests.get(url, headers=headers, params={"expand": "0"}, timeout=30)
            response.raise_for_status()
            return _record_to_region(response.json(), gene, upstream, downstream)
        except requests.exceptions.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else 0
            if status == 404:
                return None
            if status == 429:
                print("   ! Ensembl rate limit; sleeping 60s")
                time.sleep(60)
                continue
            return None
        except requests.exceptions.Timeout:
            time.sleep(2 ** attempt)
        except Exception:
            return None
    return None


def _record_to_region(record: dict, gene: str, upstream: int, downstream: int) -> Optional[Region]:
    try:
        chrom = record["seq_region_name"]
        start = int(record["start"])
        end = int(record["end"])
        strand = int(record["strand"])
    except (KeyError, TypeError, ValueError):
        return None

    if strand == 1:
        tss = start
        promoter_start = max(1, tss - upstream)
        promoter_end = tss + downstream
        strand_symbol = "+"
    else:
        tss = end
        promoter_start = max(1, tss - downstream)
        promoter_end = tss + upstream
        strand_symbol = "-"

    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    return chrom, promoter_start, promoter_end, gene, strand_symbol


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", nargs="+", required=True)
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--upstream", type=int, default=2000)
    parser.add_argument("--downstream", type=int, default=500)
    args = parser.parse_args(argv)
    regions, failed = convert_genes_to_regions(args.genes, args.genome, args.upstream, args.downstream)
    print(f"✓ Resolved {len(regions)} of {len(args.genes)} genes")
    for region in regions:
        print(f"   {region[3]}: {region[0]}:{region[1]}-{region[2]} ({region[4]} strand)")
    if failed:
        print(f"   failed: {failed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
