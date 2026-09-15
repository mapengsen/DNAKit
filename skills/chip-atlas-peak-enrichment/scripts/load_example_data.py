"""
Return bundled example gene sets for quick demos.

Bundled sets
------------
* `tp53_targets`        — 5 canonical TP53 targets (~2 min runtime)
* `immune_response`     — 20 cytokine/TF genes (~5 min runtime)
* `oncogenic_signature` — 50-gene panel covering proliferation, apoptosis,
                          and immune signaling

Inputs : example name string
Outputs: dict with `genes` and `description`

Example
-------
    python load_example_data.py --name tp53_targets
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List


EXAMPLES: Dict[str, Dict] = {
    "tp53_targets": {
        "genes": ["CDKN1A", "BAX", "BBC3", "GADD45A", "MDM2"],
        "description": (
            "Five canonical TP53 transcriptional targets covering cell-cycle arrest "
            "(CDKN1A/p21, GADD45A), apoptosis (BAX, BBC3/PUMA), and the negative "
            "feedback loop (MDM2). A fast demonstration set; results should be "
            "interpreted as exploratory at N=5."
        ),
    },
    "immune_response": {
        "genes": [
            "IL1A", "IL1B", "IL6", "CXCL8", "TNF", "IFNG", "IL10",
            "CCL2", "CCL5", "CXCL10", "CD40", "CD80", "CD86",
            "STAT1", "STAT3", "NFKB1", "RELA", "JUN", "FOS", "MYD88",
        ],
        "description": (
            "Twenty inflammation-associated genes (cytokines, chemokines, and "
            "transcription factors). Useful for surfacing NF-kB, STAT, and AP-1 family "
            "regulators in the enrichment output."
        ),
    },
    "oncogenic_signature": {
        "genes": [
            "MYC", "CCND1", "CDK4", "CDK6", "E2F1", "E2F2", "MCM2", "MCM4",
            "BIRC5", "CDC20", "AURKA", "AURKB", "FOXM1", "TOP2A", "MKI67",
            "BAX", "BAK1", "BBC3", "BIRC3", "TP53", "TP63", "TP73",
            "CDKN1A", "CDKN2A", "RB1", "MDM2", "MDM4",
            "STAT1", "STAT3", "IRF1", "IRF3", "IRF7",
            "NFKB1", "RELA", "RELB", "NFKBIA", "JUN", "JUNB", "FOS",
            "EGFR", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MTOR",
            "PIK3CA", "AKT1", "PTEN", "GSK3B", "CTNNB1",
        ],
        "description": (
            "Fifty-gene oncogenic panel spanning cell cycle, apoptosis, p53 family, "
            "MAPK/PI3K signaling, and immune-response factors. Useful for stress-testing "
            "the workflow with a more biologically heterogeneous list."
        ),
    },
}


def load_example_data(name: str = "tp53_targets") -> Dict:
    if name not in EXAMPLES:
        raise ValueError(
            f"Unknown example '{name}'. Available: {sorted(EXAMPLES)}"
        )
    payload = {**EXAMPLES[name]}
    payload["name"] = name
    print(f"✓ Example loaded: {name} ({len(payload['genes'])} genes)")
    print(f"   {payload['description']}")
    return payload


def list_examples() -> List[str]:
    return sorted(EXAMPLES)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="tp53_targets", choices=list_examples())
    args = parser.parse_args(argv)
    payload = load_example_data(args.name)
    print(f"✓ Genes: {payload['genes'][:10]}{'...' if len(payload['genes']) > 10 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
