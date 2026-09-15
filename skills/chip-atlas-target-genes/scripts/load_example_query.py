"""
Bundled example queries for the ChIP-Atlas target-genes workflow.

Each query specifies a protein, genome assembly, and TSS distance — the
three required inputs to `run_target_genes_workflow`. Use them when you
need a quick demo or smoke test.

Inputs : query name string
Outputs: dict with protein/genome/distance/description

Example
-------
    python load_example_query.py --name tp53
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List


EXAMPLES: Dict[str, Dict] = {
    "tp53": {
        "protein": "TP53",
        "genome": "hg38",
        "distance": 5,
        "description": (
            "Tumour suppressor and master transcription factor. Large dataset "
            "(~395 experiments, ~16K candidate target genes); ~13 MB TSV. "
            "Canonical targets include CDKN1A, BAX, and MDM2; runtime ~10–30 s."
        ),
    },
    "e2f1": {
        "protein": "E2F1",
        "genome": "hg38",
        "distance": 5,
        "description": (
            "Cell-cycle TF controlling G1/S transition. Moderate dataset; "
            "expected to enrich for replication and DNA-damage genes. Runtime ~5–15 s."
        ),
    },
    "myc": {
        "protein": "MYC",
        "genome": "hg38",
        "distance": 5,
        "description": (
            "Proto-oncogene with broad genome-wide binding. Many target genes; "
            "useful for demonstrating that strong universal binders look very "
            "different from cell-type-specific ones. Runtime ~5–15 s."
        ),
    },
}


def load_example_query(name: str = "tp53") -> Dict:
    key = name.lower().strip()
    if key not in EXAMPLES:
        available = sorted(EXAMPLES)
        raise ValueError(f"Unknown example '{name}'. Available: {available}")
    payload = {**EXAMPLES[key]}
    payload["name"] = key
    print(f"✓ Query loaded: {payload['protein']} target genes "
          f"({payload['genome']}, ±{payload['distance']}kb)")
    print(f"   {payload['description']}")
    return payload


def list_examples() -> List[str]:
    return sorted(EXAMPLES)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="tp53", choices=list_examples())
    args = parser.parse_args(argv)
    load_example_query(args.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
