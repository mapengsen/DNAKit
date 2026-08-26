"""Tests for interpretable advanced fingerprints and preprocessing."""

import json
from collections.abc import Iterator

import pytest

from dnakit.core import (
    BackendInfo,
    DNARecord,
    DNASequence,
    ExecutionMode,
    Gap,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
)
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.fingerprints import (
    coding_fingerprint,
    fit_preprocessor,
    gc_spatial_fingerprint,
    hybrid_fingerprint,
    motif_fingerprint,
    multiscale_fingerprint,
    repeat_fingerprint,
    restriction_fingerprint,
    thermodynamic_fingerprint,
)
from dnakit.thermodynamics import Primer3ThermodynamicResult, ThermodynamicConditions


class _FakeFingerprintStructureAdapter:
    def __init__(self) -> None:
        self.info = BackendInfo(
            "fake-primer3",
            version="1.0",
            capabilities=("hairpin", "self_dimer", "heterodimer"),
        )
        self.calls: list[str] = []

    def _result(
        self,
        capability: str,
        sequences: tuple[DNASequence, ...],
        *,
        conditions: ThermodynamicConditions | None,
        max_loop: int,
    ) -> Primer3ThermodynamicResult:
        self.calls.append(capability)
        return Primer3ThermodynamicResult(
            capability=capability,
            sequences_5to3=tuple(sequence.symbols for sequence in sequences),
            structure_found=capability != "self_dimer",
            tm_celsius=41.0,
            delta_g_kcal_per_mol=-3.0,
            delta_h_kcal_per_mol=-11.0,
            delta_s_cal_per_k_mol=-21.0,
            ascii_structure=None,
            conditions=ThermodynamicConditions() if conditions is None else conditions,
            max_loop=max_loop,
            method="fake-primer3",
            algorithm_version="1",
            backend=self.info,
            provenance=Provenance(
                implementation=ImplementationInfo(
                    label=ImplementationLabel.ADAPTER,
                    execution_mode=ExecutionMode.EXTERNAL,
                    origin_class=OriginClass.INTEGRATION,
                ),
                backend=self.info,
            ),
            issues=(),
        )

    def hairpin(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._result("hairpin", (sequence,), conditions=conditions, max_loop=max_loop)

    def self_dimer(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._result("self_dimer", (sequence,), conditions=conditions, max_loop=max_loop)

    def heterodimer(
        self,
        sequence_a: DNASequence,
        sequence_b: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._result(
            "heterodimer",
            (sequence_a, sequence_b),
            conditions=conditions,
            max_loop=max_loop,
        )


def test_interpretable_fingerprints_have_fixed_schemas() -> None:
    value = DNARecord(DNASequence("ATGATATATGAATTC"), "x")
    motif = motif_fingerprint(value, {"start": "ATG", "eco": "GAATTC"})
    restriction = restriction_fingerprint(value, ["EcoRI", "HaeIII"])
    gc = gc_spatial_fingerprint(value, bins=3)
    repeat = repeat_fingerprint(value)
    coding = coding_fingerprint(value)

    assert motif.feature_names == ("motif:eco", "motif:start")
    assert motif.values == (1.0, 2.0)
    assert restriction.feature_names == ("restriction:EcoRI", "restriction:HaeIII")
    assert len(gc.values) == 3
    assert repeat.values[0] > 0
    assert coding.values[0] == 5
    for result in (motif, restriction, gc, repeat, coding):
        json.dumps(result.to_dict(), sort_keys=True)


def test_thermodynamic_multiscale_and_hybrid_fingerprints() -> None:
    value = DNARecord(DNASequence("CGTTCCAAAGATGTGGGCATGAGCTTAC"), "x")
    thermo = thermodynamic_fingerprint(value)
    multiscale = multiscale_fingerprint(value, k_values=(1, 2))
    hybrid = hybrid_fingerprint({"thermo": thermo, "multi": multiscale}, weights={"thermo": 0.5})

    assert thermo.feature_names[0] == "tm_celsius"
    assert thermo.schema_version == "dnakit.thermodynamic_fingerprint.v2"
    assert len(thermo.feature_names) == 16
    assert thermo.parameters["hairpin_dimer_included"] is False
    assert thermo.values[4:8] == (0.0, 0.0, 0.0, 0.0)
    assert multiscale.feature_names[-1] == "global_gc"
    assert hybrid.sequence_id == "x"
    assert len(hybrid.values) == len(thermo.values) + len(multiscale.values)


def test_thermodynamic_fingerprint_explicit_adapter_and_missing_strategies() -> None:
    value = DNARecord(DNASequence("CGTTCCAAAGATGTGGGCATGAGCTTAC"), "x")
    adapter = _FakeFingerprintStructureAdapter()
    result = thermodynamic_fingerprint(
        value,
        structure_adapter=adapter,
        paired_value=DNASequence("GTAAGCTCATGCCCACATCTTTGGAACG"),
    )
    sentinel = thermodynamic_fingerprint(
        value,
        missing_strategy="sentinel",
        missing_value=-123.0,
    )

    assert adapter.calls == ["hairpin", "self_dimer", "heterodimer"]
    assert result.feature_names == sentinel.feature_names
    assert result.values[4:8] == (1.0, 1.0, 41.0, -3.0)
    assert result.values[8:12] == (1.0, 0.0, 41.0, -3.0)
    assert result.values[12:16] == (1.0, 1.0, 41.0, -3.0)
    assert result.parameters["automatic_backend_probe"] is False
    assert result.parameters["heterodimer_included"] is True
    assert result.provenance.backend == adapter.info
    assert sentinel.values[4:8] == (0.0, -123.0, -123.0, -123.0)
    with pytest.raises(ConfigurationError, match="missing_strategy=error"):
        thermodynamic_fingerprint(value, missing_strategy="error")
    with pytest.raises(ConfigurationError, match="structure_adapter"):
        thermodynamic_fingerprint(value, structure_adapter=object())  # type: ignore[arg-type]


def test_feature_preprocessor_fit_transform_and_variance_filter() -> None:
    processor = fit_preprocessor(
        [[1.0, None, 3.0], [3.0, 6.0, 3.0], [5.0, 4.0, 3.0]],
        feature_names=("a", "b", "constant"),
        mode="standard",
        missing_strategy="mean",
        variance_threshold=0.1,
    )
    transformed = processor.transform([[1.0, None, 3.0]])

    assert processor.kept_indices == (0, 1)
    assert transformed[0][0] == pytest.approx(-1.224744871)
    assert transformed[0][1] == pytest.approx(0)
    json.dumps(processor.to_dict(), sort_keys=True)


def test_preprocessor_none_minmax_and_norm_modes_do_not_accidentally_center() -> None:
    none = fit_preprocessor(
        [[1.0], [3.0]],
        feature_names=("x",),
        mode="none",
        missing_strategy="mean",
    )
    minmax = fit_preprocessor(
        [[1.0], [3.0]],
        feature_names=("x",),
        mode="minmax",
        missing_strategy="mean",
    )
    l1 = fit_preprocessor(
        [[1.0, 3.0]],
        feature_names=("x", "y"),
        mode="l1",
        missing_strategy="mean",
    )

    assert none.transform([[None], [1.0]]) == ((2.0,), (1.0,))
    assert minmax.transform([[None]]) == ((0.5,),)
    assert l1.transform([[1.0, 3.0]]) == ((0.25, 0.75),)


def test_spatial_fingerprint_rejects_silent_gap_omission() -> None:
    with pytest.raises(UnsupportedGapOperationError) as error:
        gc_spatial_fingerprint(DNASequence(["AC", Gap(2), "GT"]))
    assert error.value.code == "SPATIAL_FINGERPRINT_GAP_NOT_ALLOWED"


def test_advanced_fingerprint_parameter_errors_are_clear() -> None:
    value = DNASequence("ACGT")
    with pytest.raises(ConfigurationError):
        motif_fingerprint(value, {})
    with pytest.raises(ConfigurationError):
        gc_spatial_fingerprint(value, bins=True)
    with pytest.raises(ConfigurationError):
        multiscale_fingerprint(value, k_values=(1, 1))
    with pytest.raises(ConfigurationError):
        fit_preprocessor([[1, None]], feature_names=("a", "b"))
    with pytest.raises(ConfigurationError):
        motif_fingerprint(value, {1: "AC"})  # type: ignore[dict-item]
    with pytest.raises(ConfigurationError):
        restriction_fingerprint(value, "EcoRI")
    with pytest.raises(ConfigurationError):
        fit_preprocessor([[1]], feature_names="a")


def test_restriction_fingerprint_bounds_enzyme_iterable() -> None:
    consumed = 0

    def enzymes() -> Iterator[str]:
        nonlocal consumed
        while True:
            consumed += 1
            yield "EcoRI"

    with pytest.raises(ConfigurationError, match="max_enzymes"):
        restriction_fingerprint(DNASequence("GAATTC"), enzymes(), max_enzymes=2)
    assert consumed == 3


def test_preprocessor_bounds_each_feature_row_to_schema_width() -> None:
    consumed = 0

    def row() -> Iterator[float]:
        nonlocal consumed
        while True:
            consumed += 1
            yield 1.0

    with pytest.raises(ConfigurationError, match="width"):
        fit_preprocessor([row()], feature_names=("x",))
    assert consumed == 2
