"""Keep the fixed, read-only documentation demo synchronized with live APIs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dnakit.datasets import deduplicate
from dnakit.descriptors import (
    exact_repeat_fraction,
    gc_at_content,
    kmer_statistics,
    linguistic_complexity,
)
from dnakit.evaluation import evaluate_synthesis_risk
from dnakit.io import read_set
from dnakit.ops import reverse_complement
from dnakit.thermodynamics import melting_temperature


def test_fixed_demo_matches_current_deterministic_apis() -> None:
    root = Path(__file__).resolve().parents[1]
    source_payload = json.loads(
        (root / "examples/fixed_demo_expected.json").read_text(encoding="utf-8")
    )
    docs_payload = json.loads((root / "docs/demo/data/fixed_demo.json").read_text(encoding="utf-8"))
    records = read_set(root / "examples/fixed_demo.fasta")
    record = records[0]
    expected = source_payload["expected"]["seq-a"]

    assert docs_payload["expected"] == source_payload["expected"]
    assert docs_payload["records"] == [
        {"id": item.id, "sequence": item.sequence.symbols} for item in records
    ]
    assert source_payload["source"] == "fixed_demo.fasta"
    assert expected["symbol_length"] == record.sequence.symbol_length
    assert expected["gc_fraction"] == gc_at_content(record).gc_fraction
    assert expected["kmer_counts_k2"] == dict(kmer_statistics(record, 2).counts)
    assert expected["reverse_complement"] == reverse_complement(record.sequence).symbols
    assert expected["linguistic_complexity_k1_to_k3"] == pytest.approx(
        linguistic_complexity(record, max_word_size=3).score
    )
    assert expected["exact_tandem_repeat_fraction"] == exact_repeat_fraction(record).repeat_fraction
    assert (
        expected["wallace_tm_celsius"]
        == melting_temperature(record.sequence, method="wallace").tm_celsius
    )
    risk = evaluate_synthesis_risk(record)
    assert expected["synthesis_risk_level"] == risk.entries[0].metrics["risk_level"]
    dataset_expected = source_payload["expected"]["dataset"]
    deduplicated = deduplicate(records)
    assert dataset_expected["record_count"] == len(records)
    assert dataset_expected["unique_raw_sequence_count"] == deduplicated.output_count
    assert dataset_expected["exact_duplicate_groups"] == [
        list(group.member_ids) for group in deduplicated.groups if len(group.member_ids) > 1
    ]
