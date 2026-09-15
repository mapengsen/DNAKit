"""
Validate a user-supplied target-genes query.

Inputs : protein name, genome assembly, distance
Outputs: dict (protein, genome, distance, description)

The validator catches the easiest mistakes (case, unsupported genome,
unsupported distance) before the workflow hits the network.

Example
-------
    python load_user_query.py --protein TP53 --genome hg38 --distance 5
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict


SUPPORTED_GENOMES = (
    "hg38", "hg19", "mm10", "mm9", "rn6",
    "dm6", "dm3", "ce11", "ce10", "sacCer3",
)
VALID_DISTANCES = (1, 5, 10)


def load_user_query(protein: str, genome: str = "hg38", distance: int = 5) -> Dict:
    if not isinstance(protein, str) or not protein.strip():
        raise ValueError("Protein name must be a non-empty string (e.g. 'TP53', 'CTCF', 'MYC')")
    protein = protein.strip()

    if genome not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'. Supported: {SUPPORTED_GENOMES}")
    if distance not in VALID_DISTANCES:
        raise ValueError(f"Unsupported distance {distance}. Supported: {VALID_DISTANCES}")

    query = {
        "protein": protein,
        "genome": genome,
        "distance": distance,
        "description": f"User query: {protein} target genes ({genome}, ±{distance}kb from TSS)",
    }
    print(f"✓ Query validated: {protein} ({genome}, ±{distance}kb)")
    return query


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--distance", type=int, default=5, choices=VALID_DISTANCES)
    args = parser.parse_args(argv)
    load_user_query(args.protein, args.genome, args.distance)
    return 0


if __name__ == "__main__":
    sys.exit(main())
