"""Offline Python API example: native features, metrics, optimization and a DAG.

This version demonstrates the new bioinformatics APIs without files, network
requests or specialist executables. Install DNAKit's bio-analysis extra first.
"""

from __future__ import annotations

import json

from dnakit.bioinformatics.alignment import MultipleAlignment
from dnakit.bioinformatics.cloning import CloningPlan, Duplex
from dnakit.bioinformatics.encoding import encode_dna
from dnakit.bioinformatics.metrics import kernel_distance, property_hit_rate, sequence_properties
from dnakit.bioinformatics.optimization import optimize_sequence, specification
from dnakit.bioinformatics.phylogeny import build_tree
from dnakit.bioinformatics.workflows import Step, run_workflow


def gc_fraction(sequence: str) -> float:
    return sum(base in "GC" for base in sequence) / len(sequence)


def main() -> None:
    sequences = {"a": "AACGTG", "b": "CCGTAT", "c": "ATGCAA"}
    encoded = encode_dna(sequences, "RCKmer", parameters={"k": 2})
    properties = sequence_properties(sequences, {"gc": gc_fraction})
    optimized = optimize_sequence(
        "GCTGCT",
        constraints=[specification("avoid_pattern", pattern="GCT")],
        preserve_translation=True,
        seed=0,
    )
    msa = MultipleAlignment.from_sequences(sequences)
    tree = build_tree(msa.distances(), msa.ids)
    cloning = CloningPlan()
    cloning.add("template", Duplex.from_sequence("AACCGGTT"))
    cut = cloning.apply("cut", lambda molecule: molecule.cut(3, 5), parents=["template"])
    joined = cut.product[0].ligate(cut.product[1])
    workflow = run_workflow(
        [
            Step("input", lambda: sequences),
            Step(
                "encoding",
                encode_dna,
                dependencies={"sequences": "input"},
                parameters={"descriptor": "NAC"},
            ),
        ]
    )
    print(
        json.dumps(
            {
                "encoding": {"schema": encoded.schema_version, "shape": list(encoded.values.shape)},
                "self_biased_mmd2": kernel_distance(encoded.values, encoded.values, unbiased=False),
                "gc_hit_rate": property_hit_rate(properties["values"], lower=[0.4], upper=[0.6])[
                    "joint"
                ],
                "optimization": {
                    "sequence": optimized.sequence,
                    "satisfied": optimized.constraints_satisfied,
                },
                "tree_tips": tree.tips,
                "ligation_roundtrip": joined.watson == "AACCGGTT",
                "workflow_status": {step.name: step.status for step in workflow},
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
