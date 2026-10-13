"""Regression tests for the formal ``dnakit.<category>`` candidate surface."""

from __future__ import annotations

import csv
import importlib
from pathlib import Path

import dnakit
from dnakit.tools.catalog import build_tool_catalog
from dnakit.tools.registry import default_tool_registry


def test_candidate_categories_are_available_from_the_formal_namespace() -> None:
    assert dnakit.sequence.reverse_complement_dna("ATGC") == "GCAT"
    assert dnakit.sequence.mutate_sequence_at_positions("AAAA", {1: "C"}) == "ACAA"
    n50 = dnakit.assembly.calculate_n50([100, 60, 40], options={"allow_local_fallback": True})
    assert n50["status"] == "completed"
    assert n50["result"]["n50"] == 100
    figure = dnakit.visualization.plot({"records": ["ATGC"]})
    assert figure["status"] == "ok"
    assert figure["result"].startswith("<svg")


def test_candidate_functions_are_discoverable_by_the_agent_registry() -> None:
    catalog = build_tool_catalog(verbose=False)
    names = {tool.name for tool in catalog}
    assert "dnakit.sequence.reverse_complement_dna" in names
    assert "dnakit.assembly.calculate_n50" in names
    assert "dnakit.visualization.plot" in names

    registry = default_tool_registry()
    stats = registry.stats()
    assert stats["total"] == len(catalog)
    assert stats["agent_compatible"] == stats["total"] - 3
    assert (
        registry.execute(
            "dnakit.sequence.reverse_complement_dna",
            {"sequence": "ATGC"},
        )
        == "GCAT"
    )


def test_all_in_scope_catalog_candidates_are_formal_dnakit_exports() -> None:
    root = Path(__file__).parents[3]
    catalog = list(csv.DictReader((root / "new_vision/DNAKit_unified_function_catalog.csv").open(encoding="utf-8")))
    comparison = {
        row["capability_id"]: row
        for row in csv.DictReader((root / "new_vision/DNAKit_unified_comparison.csv").open(encoding="utf-8"))
    }
    missing: list[str] = []
    for row in catalog:
        if comparison[row["capability_id"]]["decision"] not in {"new_candidate", "review_needed"}:
            continue
        module = importlib.import_module(f"dnakit.{row['module']}")
        if not callable(getattr(module, row["operation"], None)):
            missing.append(f"dnakit.{row['module']}.{row['operation']}")
    assert not missing


def test_external_qc_facades_share_the_agent_mcp_effect_contract() -> None:
    registry = default_tool_registry()
    operations = tuple(
        ("dnakit.sequencing_qc", name)
        for name in ("run_cutadapt", "run_fastp", "run_fastqc", "run_multiqc")
    ) + tuple(
        ("dnakit.ngs", name)
        for name in ("map_reads", "prepare_bam", "call_variants", "assemble_reads")
    )
    for module, operation in operations:
        tool = registry.resolve(f"{module}.{operation}")
        assert tool.effect == "external_process"
        assert tool.agent_compatible
        assert tool.requires_confirmation
