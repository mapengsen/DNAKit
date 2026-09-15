"""
Return pre-baked experiment-group payloads for quick demos.

The bundled examples both use TP53 ChIP-seq from GSE131484. They are
intentionally non-ideal designs (each group mixes multiple TP53 genotypes)
so the workflow's QC layer has something to flag.

Inputs : a dataset key string
Outputs: a dict ready to splat into `run_diff_workflow`

Example
-------
    python load_example_data.py --dataset tp53_k562_dmso_vs_daunorubicin
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List

EXAMPLES: Dict[str, Dict] = {
    "tp53_k562_dmso_vs_daunorubicin": {
        "experiments_a": [
            "SRX5865959", "SRX5865960", "SRX5865961", "SRX5865962",
            "SRX5865963", "SRX5865964", "SRX5865965", "SRX5865966",
        ],
        "experiments_b": [
            "SRX5865967", "SRX5865968", "SRX5865969", "SRX5865970",
            "SRX5865971", "SRX5865972", "SRX5865973", "SRX5865974",
        ],
        "genome": "hg38",
        "analysis_type": "diffbind",
        "description_a": "K562_DMSO",
        "description_b": "K562_Daunorubicin",
        "summary": (
            "TP53 ChIP-seq in K-562, vehicle (DMSO) vs daunorubicin treatment. "
            "8 vs 8 experiments spanning the WT and 6 mutant TP53 backgrounds from "
            "GSE131484. Identifies how daunorubicin-induced DNA damage reshapes "
            "TP53 occupancy."
        ),
        "design_caveats": [
            "Each group bundles eight different TP53 genotypes (WT + 6 hotspot mutants + KO), "
            "not eight clean biological replicates. edgeR treats within-group experiments "
            "as replicates, so the result is the average DMSO-vs-daunorubicin effect across "
            "all genotypes. Genotype-specific responses are averaged out and cannot be "
            "disentangled from noise."
        ],
    },
    "tp53_molm13_vs_k562": {
        "experiments_a": [
            "SRX5865975", "SRX5865976", "SRX5865977", "SRX5865982",
        ],
        "experiments_b": [
            "SRX5865959", "SRX5865960", "SRX5865961", "SRX5865966",
        ],
        "genome": "hg38",
        "analysis_type": "diffbind",
        "description_a": "MOLM-13",
        "description_b": "K-562",
        "summary": (
            "TP53 ChIP-seq comparing two leukemia lines under DMSO. "
            "MOLM-13 (AML, male) vs K-562 (CML, female) with matched TP53 backgrounds "
            "(R175H, Y220C, M237I, WT) from GSE131484. Triggers the sex-chromosome "
            "QC warning because the two lines differ in sex."
        ),
        "design_caveats": [
            "Both groups span four TP53 genotypes (R175H, Y220C, M237I, WT). edgeR "
            "treats them as replicates, so the analysis recovers the average MOLM-13-vs-"
            "K-562 difference across those genotypes.",
            "MOLM-13 is male and K-562 is female, so any chrY regions and many chrX "
            "regions will look strongly differential. The QC system flags this; do not "
            "interpret sex-chromosome hits as real differential binding.",
        ],
    },
}


def load_example_data(dataset: str = "tp53_k562_dmso_vs_daunorubicin") -> Dict:
    if dataset not in EXAMPLES:
        raise ValueError(
            f"Unknown example '{dataset}'. Available: {sorted(EXAMPLES)}"
        )
    payload = dict(EXAMPLES[dataset])
    n_a = len(payload["experiments_a"])
    n_b = len(payload["experiments_b"])
    print(f"✓ Example loaded: {dataset} ({n_a} vs {n_b} experiments, genome={payload['genome']})")
    print(f"   {payload['summary']}")
    return payload


def list_examples() -> List[str]:
    return sorted(EXAMPLES)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", default="tp53_k562_dmso_vs_daunorubicin", choices=list_examples())
    args = parser.parse_args(argv)
    data = load_example_data(args.dataset)
    print(f"✓ Group A: {len(data['experiments_a'])} IDs; Group B: {len(data['experiments_b'])} IDs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
