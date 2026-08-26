"""Human-verifiable tests for the fixed 240-field descriptor schema."""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap, Topology
from dnakit.core._json import FrozenDict
from dnakit.descriptors import (
    DESCRIPTOR_NAMES_V1,
    DESCRIPTOR_SCHEMA_V1,
    DINUCLEOTIDE_PROPERTY_SPECS,
    DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
    DINUCLEOTIDES,
    DinucleotideProperty,
    DinucleotidePropertyTable,
    all_descriptors,
    descriptor_schema_v1,
    load_dinucleotide_property_table,
)
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.thermodynamics import ThermodynamicConditions


def _user_table(tmp_path: Path) -> DinucleotidePropertyTable:
    payload = {
        "schema_version": DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
        "name": "synthetic-test-table",
        "version": "1",
        "source": "generated only inside this unit test",
        "properties": {
            spec.key: {
                "unit": spec.unit,
                "values": {
                    word: property_index * 100.0 + word_index
                    for word_index, word in enumerate(DINUCLEOTIDES)
                },
            }
            for property_index, spec in enumerate(DINUCLEOTIDE_PROPERTY_SPECS)
        },
    }
    path = tmp_path / "user-table.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return load_dinucleotide_property_table(path)


def test_descriptor_schema_has_exact_stable_240_field_layout() -> None:
    schema = descriptor_schema_v1()

    assert schema is DESCRIPTOR_SCHEMA_V1
    assert len(schema) == 240
    assert tuple(field.index for field in schema) == tuple(range(1, 241))
    assert tuple(field.name for field in schema) == DESCRIPTOR_NAMES_V1
    assert len(set(DESCRIPTOR_NAMES_V1)) == 240
    assert Counter(field.category for field in schema) == {
        "basic": 12,
        "composition": 16,
        "kmer": 84,
        "skew_cpg": 16,
        "complexity": 20,
        "coding": 16,
        "physicochemical": 16,
        "dinucleotide_property": 60,
    }
    assert all(field.unit and field.formula and field.source for field in schema)
    assert schema[164].name == "mw_ss_oh_da"
    assert schema[180].name == "diprodb_twist_mean"


def test_complete_descriptor_document_lists_schema_in_exact_order() -> None:
    document = (Path(__file__).parents[2] / "docs/api/features/05_all_descriptors.md").read_text(
        encoding="utf-8"
    )
    documented_rows = tuple(
        tuple(cell.strip().strip("`") for cell in line.split("|")[1:-1])
        for line in document.splitlines()
        if line.startswith("| ") and line.split("|")[1].strip().isdigit()
    )

    assert documented_rows == tuple(
        (
            str(field.index),
            field.name,
            field.category,
            field.unit,
            field.formula,
            field.source,
        )
        for field in DESCRIPTOR_SCHEMA_V1
    )


def test_all_descriptors_returns_ordered_values_and_auditable_context() -> None:
    record = DNARecord(DNASequence("ACGT"), "example")
    result = all_descriptors(record)

    assert result.schema_version == "descriptor_schema_v1"
    assert result.sequence_id == "example"
    assert tuple(result.values) == DESCRIPTOR_NAMES_V1
    assert len(result.values) == 240
    assert set(result.unavailable_reasons) == {
        name for name, value in result.values.items() if value is None
    }
    assert result.conditions["ambiguity_policy"] == "ignore"
    assert result.conditions["cross_gaps"] is False
    assert result.conditions["dinucleotide_property_table"] is None
    assert result.provenance.reference is None
    assert all(result.values[name] is None for name in DESCRIPTOR_NAMES_V1[180:])
    assert len(result.to_dict()["values"]) == 240


def test_basic_composition_kmer_skew_and_entropy_are_hand_checkable() -> None:
    values = all_descriptors(DNASequence("ACGT")).values

    assert values["symbol_length"] == 4
    assert values["coordinate_span"] == 4
    assert values["canonical_run_count"] == 1
    assert values["canonical_symbol_fraction"] == 1.0
    assert values["purine_count"] == 2
    assert values["purine_fraction"] == 0.5
    assert values["gc_at_ratio"] == 1.0
    assert values["k1_A_frequency"] == 0.25
    assert values["k2_AC_frequency"] == pytest.approx(1 / 3)
    assert values["k2_AA_frequency"] == 0.0
    assert values["k3_ACG_frequency"] == 0.5
    assert values["gc_skew"] == 0.0
    assert values["at_skew"] == 0.0
    assert values["cpg_count"] == 1
    assert values["cpg_density"] == pytest.approx(1 / 3)
    assert values["cpg_observed_expected"] == 4.0
    assert values["cumulative_gc_skew_range"] == 1
    assert values["cumulative_at_skew_range"] == 1
    assert values["dinucleotide_rc_total_variation"] == 0.0
    assert values["mono_chargaff_l1_distance"] == 0.0
    assert values["shannon_entropy_k1_bits"] == 2.0
    assert values["shannon_entropy_k2_bits"] == pytest.approx(math.log2(3))
    assert values["normalized_entropy_k1"] == 1.0


def test_complexity_coding_physchem_and_user_table_values_are_calculated(tmp_path: Path) -> None:
    table = _user_table(tmp_path)
    acgt = all_descriptors(DNASequence("ACGT")).values
    coding = all_descriptors(DNASequence("ATGAAATAA")).values
    aa = all_descriptors(DNASequence("AA"), dinucleotide_property_table=table).values

    assert acgt["linguistic_complexity_k2"] == 1.0
    assert acgt["linguistic_complexity_product_k1_k6"] == 1.0
    assert acgt["lz76_complexity"] == 4
    assert acgt["longest_homopolymer_nt"] == 1
    assert acgt["exact_tandem_repeat_coverage_fraction"] == 0.0
    assert coding["frame0_codon_count"] == 3
    assert coding["frame0_unique_codon_count"] == 3
    assert coding["frame0_start_codon_fraction"] == pytest.approx(1 / 3)
    assert coding["frame0_stop_codon_fraction"] == pytest.approx(1 / 3)
    assert coding["frame0_effective_number_of_codons"] == pytest.approx(3.0)
    assert coding["six_frame_complete_orf_count"] == 1
    assert coding["six_frame_forward_complete_orf_count"] == 1
    assert coding["six_frame_reverse_complete_orf_count"] == 0
    assert coding["six_frame_longest_complete_orf_nt"] == 9
    assert coding["six_frame_complete_orf_coverage_fraction"] == 1.0
    assert acgt["mw_ss_oh_da"] == pytest.approx(1173.84)
    assert acgt["epsilon260_ss_m_inverse_cm_inverse"] == 40300.0
    assert acgt["nmol_per_a260_1ml_1cm"] == pytest.approx(1_000_000 / 40300)
    assert acgt["ug_per_a260_1ml_1cm"] == pytest.approx(1000 * 1173.84 / 40300)
    assert acgt["tm_wallace_c"] == 12.0
    assert acgt["self_complementary"] is True
    assert aa["diprodb_twist_mean"] == 0.0
    assert aa["diprodb_twist_sd"] == 0.0
    assert aa["diprodb_stacking_energy_mean"] == 1300.0


def test_gap_and_ambiguity_are_boundaries_and_unavailable_values_have_reasons(
    tmp_path: Path,
) -> None:
    table = _user_table(tmp_path)
    gapped = DNASequence(["CG", Gap(None), "GC"])
    result = all_descriptors(gapped, dinucleotide_property_table=table)
    values = result.values

    assert values["coordinate_span"] is None
    assert values["canonical_run_count"] == 2
    assert values["k2_CG_frequency"] == 0.5
    assert values["k2_GC_frequency"] == 0.5
    assert values["k2_GG_frequency"] == 0.0
    assert values["cpg_gpc_ratio"] == 1.0
    assert values["diprodb_twist_mean"] == 7.5
    assert values["mw_ss_oh_da"] is None
    assert "ungapped" in str(result.unavailable_reasons["mw_ss_oh_da"])
    assert values["lz76_complexity"] is None

    ambiguous = DNASequence("CNG", alphabet=DNAAlphabet.IUPAC)
    ignored = all_descriptors(ambiguous)
    assert ignored.values["canonical_base_count"] == 2
    assert ignored.values["cpg_count"] == 0
    assert ignored.values["cpg_density"] is None
    assert ignored.values["diprodb_twist_mean"] is None
    with pytest.raises(InvalidAlphabetError):
        all_descriptors(ambiguous, ambiguity_policy="error")


def test_user_table_is_checksummed_and_invalid_tables_are_rejected(tmp_path: Path) -> None:
    table = _user_table(tmp_path)
    result = all_descriptors(DNASequence("AA"), dinucleotide_property_table=table)

    assert result.provenance.reference is not None
    assert result.provenance.reference.name == "synthetic-test-table"
    assert result.provenance.reference.checksum == table.sha256
    table_conditions = result.conditions["dinucleotide_property_table"]
    assert isinstance(table_conditions, FrozenDict)
    assert table_conditions["sha256"] == table.sha256

    invalid = tmp_path / "invalid.json"
    invalid.write_text("{}", encoding="utf-8")
    with pytest.raises(ConfigurationError) as error:
        load_dinucleotide_property_table(invalid)
    assert error.value.code == "INVALID_DINUCLEOTIDE_TABLE_SCHEMA"

    with pytest.raises(ConfigurationError):
        DinucleotideProperty(
            key="twist",
            display_name="Twist",
            unit="degree",
            values={"AA": 1.0},  # type: ignore[arg-type]
        )

    oversized = tmp_path / "oversized.json"
    oversized.write_text(" " * 129, encoding="utf-8")
    with pytest.raises(ConfigurationError, match="exceeds max_bytes"):
        load_dinucleotide_property_table(oversized, max_bytes=128)


def test_empty_short_long_and_circular_domains_use_none_not_fake_zero() -> None:
    empty = all_descriptors(DNASequence(""))
    long_result = all_descriptors(DNASequence("A" * 61))
    circular = all_descriptors(DNASequence("ATGAAATAA", topology=Topology.CIRCULAR))

    assert empty.values["canonical_base_count"] == 0
    assert empty.values["canonical_symbol_fraction"] is None
    assert empty.values["k1_A_frequency"] is None
    assert empty.values["frame0_codon_count"] == 0
    assert empty.values["frame0_codon_entropy_bits"] is None
    assert empty.values["mw_ss_oh_da"] is None
    assert len(empty.unavailable_reasons) == sum(value is None for value in empty.values.values())
    assert long_result.values["mw_ss_oh_da"] is not None
    assert long_result.values["tm_wallace_c"] is None
    assert long_result.values["nn_tm_c"] is None
    assert circular.values["purine_count"] == 7
    assert circular.values["six_frame_complete_orf_count"] is None
    assert circular.values["exact_tandem_repeat_coverage_fraction"] is None
    assert circular.values["mw_ss_oh_da"] is None


def test_conditions_are_recorded_and_nn_delta_g_remains_at_37_celsius() -> None:
    result = all_descriptors(
        DNASequence("ACGTACGT"),
        conditions=ThermodynamicConditions(
            temperature_celsius=25.0,
            sodium_molar=0.1,
            strand_concentration_molar=500e-9,
        ),
    )

    assert result.conditions["temperature_celsius"] == 25.0
    assert result.conditions["free_energy_reference_temperature_celsius"] == 37.0
    assert result.conditions["sodium_molar"] == 0.1
    assert result.values["nn_delta_g37_kcal_per_mol"] is not None


def test_all_descriptor_result_and_schema_are_immutable() -> None:
    result = all_descriptors(DNASequence("ACGT"))

    with pytest.raises(TypeError):
        result.values["symbol_length"] = 0  # type: ignore[index]
    with pytest.raises(FrozenInstanceError):
        result.schema_version = "changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        DESCRIPTOR_SCHEMA_V1[0].name = "changed"  # type: ignore[misc]
