from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

import pytest

from dnakit.core import (
    BackendInfo,
    DNASequence,
    ExecutionMode,
    Gap,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    Strand,
    Topology,
)
from dnakit.core._json import FrozenDict
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    InvalidAlphabetError,
    UnsupportedGapOperationError,
)
from dnakit.molbio import (
    Primer3CLIDesignAdapter,
    match_primer,
    prepare_primer_design,
    primer_properties,
    scan_crispr_candidates,
    simulate_pcr,
)
from dnakit.thermodynamics import Primer3ThermodynamicResult, ThermodynamicConditions


class _FakeStructureAdapter:
    def __init__(self) -> None:
        self.info = BackendInfo(
            "fake-primer3",
            version="1.0",
            capabilities=("hairpin", "self_dimer", "heterodimer", "tm"),
        )
        self.calls: list[tuple[str, tuple[str, ...]]] = []

    def _result(
        self,
        capability: str,
        sequences: tuple[DNASequence, ...],
        *,
        conditions: ThermodynamicConditions | None,
        max_loop: int | None,
        output_structure: bool,
    ) -> Primer3ThermodynamicResult:
        symbols = tuple(sequence.symbols for sequence in sequences)
        self.calls.append((capability, symbols))
        resolved_conditions = ThermodynamicConditions() if conditions is None else conditions
        tm_only = capability == "tm"
        return Primer3ThermodynamicResult(
            capability=capability,
            sequences_5to3=symbols,
            structure_found=None if tm_only else True,
            tm_celsius=42.0,
            delta_g_kcal_per_mol=None if tm_only else -2.5,
            delta_h_kcal_per_mol=None if tm_only else -10.0,
            delta_s_cal_per_k_mol=None if tm_only else -20.0,
            ascii_structure="structure" if output_structure and not tm_only else None,
            conditions=resolved_conditions,
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

    def tm(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
    ) -> Primer3ThermodynamicResult:
        return self._result(
            "tm",
            (sequence,),
            conditions=conditions,
            max_loop=None,
            output_structure=False,
        )

    def hairpin(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._result(
            "hairpin",
            (sequence,),
            conditions=conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )

    def self_dimer(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._result(
            "self_dimer",
            (sequence,),
            conditions=conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )

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
            output_structure=output_structure,
        )


def _fake_primer3_core(
    path: Path,
    output: str,
    *,
    capture_path: Path | None = None,
    exit_code: int = 0,
) -> Path:
    capture = (
        ""
        if capture_path is None
        else f"Path({str(capture_path)!r}).write_text(input_text, encoding='utf-8')\n"
    )
    body = (
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        "input_text = Path(sys.argv[-1]).read_text(encoding='utf-8')\n"
        f"{capture}"
        "output_name = next(item.split('=', 1)[1] for item in sys.argv "
        "if item.startswith('--output='))\n"
        "error_name = next(item.split('=', 1)[1] for item in sys.argv "
        "if item.startswith('--error='))\n"
        f"Path(output_name).write_text({output!r}, encoding='utf-8')\n"
        "Path(error_name).write_text('', encoding='utf-8')\n"
        f"raise SystemExit({exit_code})\n"
    )
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_primer_matching_reports_orientation_mismatches_and_three_prime_policy() -> None:
    template = DNASequence("AAAACCCCTTTT")
    exact = match_primer(DNASequence("AAAA"), template)
    mismatch = match_primer(
        DNASequence("AAAT"),
        template,
        strand="forward",
        max_mismatches=1,
    )

    assert {(hit.strand, hit.start) for hit in exact.hits} == {
        (Strand.FORWARD, 0),
        (Strand.REVERSE, 8),
    }
    assert mismatch.hits[0].mismatch_positions_5to3 == (3,)
    assert (
        match_primer(
            DNASequence("AAAT"),
            template,
            strand="forward",
            max_mismatches=1,
            strict_three_prime_bases=1,
        ).hits
        == ()
    )


def test_circular_primer_full_length_hit_uses_unambiguous_end_coordinate() -> None:
    result = match_primer(
        DNASequence("ACGT"),
        DNASequence("ACGT", topology="circular"),
        strand="forward",
    )

    assert result.hits[0].start == 0
    assert result.hits[0].end == 4
    assert not result.hits[0].wraps_origin


def test_linear_and_circular_pcr_construct_primer_defined_products() -> None:
    linear = simulate_pcr(
        DNASequence("AAAACCCCGGGGTTTT"),
        DNASequence("AAAA"),
        DNASequence("AAAA"),
        min_product_length=16,
    )
    circular = simulate_pcr(
        DNASequence("TTTTCCCCGGGGAAAA", topology="circular"),
        DNASequence("AAAA"),
        DNASequence("AAAA"),
        min_product_length=8,
        max_product_length=8,
    )

    assert [amplicon.sequence.symbols for amplicon in linear.amplicons] == ["AAAACCCCGGGGTTTT"]
    assert circular.amplicons[0].sequence.symbols == "AAAATTTT"
    assert circular.amplicons[0].wraps_origin
    assert circular.amplicons[0].sequence.topology is Topology.LINEAR


def test_primer_resources_and_sequence_boundaries_are_explicit() -> None:
    with pytest.raises(UnsupportedGapOperationError):
        match_primer(DNASequence("AAAA"), DNASequence(["AA", Gap(1), "AA"]))
    with pytest.raises(InvalidAlphabetError):
        match_primer(
            DNASequence("AANA", alphabet="iupac"),
            DNASequence("AAAA"),
        )
    with pytest.raises(ConfigurationError, match="max_comparison_cells"):
        match_primer(
            DNASequence("AAAA"),
            DNASequence("A" * 100),
            max_comparison_cells=1,
        )
    with pytest.raises(ConfigurationError, match="max_pair_checks"):
        simulate_pcr(
            DNASequence("A" * 20 + "T" * 20),
            DNASequence("AA"),
            DNASequence("AA"),
            max_pair_checks=1,
        )


def test_primer_properties_compute_gc_tm_but_do_not_claim_structures() -> None:
    result = primer_properties(
        DNASequence("GTGCAT"),
        paired_primer=DNASequence("ATGCAC"),
    )

    assert result.gc_fraction == pytest.approx(0.5)
    assert isinstance(result.tm_celsius, float)
    assert not result.hairpin.available
    assert result.hairpin.automatic_probe is False
    assert result.heterodimer is not None
    assert result.paired_gc_fraction == pytest.approx(0.5)
    assert isinstance(result.paired_tm_celsius, float)
    assert result.parameters["hairpin_backend_executed"] is False
    assert json.loads(json.dumps(result.to_dict()))["issues"][0]["code"] == (
        "PRIMER_STRUCTURE_BACKEND_REQUIRED"
    )
    with pytest.raises(ConfigurationError, match="structure_adapter"):
        primer_properties(DNASequence("GTGCAT"), structure_adapter=object())  # type: ignore[arg-type]


def test_primer_properties_executes_only_an_explicit_structure_adapter() -> None:
    adapter = _FakeStructureAdapter()
    result = primer_properties(
        DNASequence("GTGCAT"),
        paired_primer=DNASequence("ATGCAC"),
        structure_adapter=adapter,
        output_structure=True,
    )

    assert adapter.calls == [
        ("hairpin", ("GTGCAT",)),
        ("self_dimer", ("GTGCAT",)),
        ("hairpin", ("ATGCAC",)),
        ("self_dimer", ("ATGCAC",)),
        ("heterodimer", ("GTGCAT", "ATGCAC")),
    ]
    assert result.hairpin.available
    assert result.hairpin.execution_performed
    assert not result.hairpin.automatic_execution
    assert result.hairpin.delta_g_kcal_per_mol == pytest.approx(-2.5)
    assert result.paired_hairpin is not None and result.paired_hairpin.available
    assert result.paired_self_dimer is not None and result.paired_self_dimer.available
    assert result.heterodimer is not None and result.heterodimer.available
    assert result.parameters["automatic_backend_probe"] is False
    assert result.parameters["structure_adapter_supplied"] is True
    assert result.issues == ()
    json.dumps(result.to_dict(), sort_keys=True)


def test_primer_properties_supports_explicit_primer3_tm_with_divalent_conditions() -> None:
    adapter = _FakeStructureAdapter()
    conditions = ThermodynamicConditions(
        magnesium_molar=0.0015,
        dntp_molar=0.0006,
    )
    result = primer_properties(
        DNASequence("GTGCAT"),
        paired_primer=DNASequence("ATGCAC"),
        tm_method="primer3-cli",
        conditions=conditions,
        structure_adapter=adapter,
    )

    assert result.tm_method == "fake-primer3"
    assert result.tm_celsius == pytest.approx(42.0)
    assert result.paired_tm_celsius == pytest.approx(42.0)
    audited_conditions = result.parameters["conditions"]
    assert isinstance(audited_conditions, FrozenDict)
    assert audited_conditions["magnesium_molar"] == pytest.approx(0.0015)
    assert result.provenance.implementation.label is ImplementationLabel.ADAPTER
    assert adapter.calls[:2] == [("tm", ("GTGCAT",)), ("tm", ("ATGCAC",))]


@pytest.mark.parametrize(
    ("keyword", "value"),
    (
        ("primer_length_range", None),
        ("primer_length_range", (18,)),
        ("tm_range_celsius", (50.0, 60.0, 70.0)),
        ("gc_range", (0.4,)),
        ("product_length_range", "100,200"),
        ("excluded_regions", None),
    ),
)
def test_primer_design_request_rejects_malformed_range_inputs(keyword: str, value: object) -> None:
    options: dict[str, object] = {keyword: value}
    with pytest.raises(ConfigurationError):
        prepare_primer_design(
            DNASequence("ACGT" * 100),
            target_start=50,
            target_end=100,
            **options,  # type: ignore[arg-type]
        )


def test_primer_design_interface_validates_request_without_execution() -> None:
    result = prepare_primer_design(
        DNASequence("ACGT" * 100),
        target_start=50,
        target_end=100,
        excluded_regions=((10, 20),),
    )

    assert result.status == "request-ready-execution-not-performed"
    assert not result.execution_performed
    assert result.candidates == ()
    assert result.parameters["automatic_probe"] is False
    assert result.request.thermodynamic_conditions["sodium_molar"] == 0.05
    with pytest.raises(ConfigurationError, match="ThermodynamicConditions"):
        prepare_primer_design(
            DNASequence("ACGT"),
            target_start=1,
            target_end=3,
            conditions={},  # type: ignore[arg-type]
        )
    with pytest.raises(ConfigurationError, match="valid 0-based"):
        prepare_primer_design(
            DNASequence("ACGT"),
            target_start=2,
            target_end=5,
        )


def test_primer3_design_adapter_executes_prepared_request_with_bounded_output(
    tmp_path: Path,
) -> None:
    prepared = prepare_primer_design(
        DNASequence("ACGT" * 100),
        target_start=50,
        target_end=100,
        excluded_regions=((10, 20),),
        candidate_count=2,
    )
    output = "\n".join(
        (
            "PRIMER_PAIR_NUM_RETURNED=1",
            "PRIMER_LEFT_0=20,20",
            "PRIMER_RIGHT_0=139,20",
            f"PRIMER_LEFT_0_SEQUENCE={'ACGT' * 5}",
            f"PRIMER_RIGHT_0_SEQUENCE={'ACGT' * 5}",
            "PRIMER_LEFT_0_TM=60.1",
            "PRIMER_RIGHT_0_TM=60.2",
            "PRIMER_LEFT_0_GC_PERCENT=50.0",
            "PRIMER_RIGHT_0_GC_PERCENT=50.0",
            "PRIMER_PAIR_0_PRODUCT_SIZE=120",
            "PRIMER_PAIR_0_PENALTY=0.25",
            "PRIMER_WARNING=bounded warning",
            "=",
            "",
        )
    )
    capture = tmp_path / "captured-input.boulder"
    executable = _fake_primer3_core(
        tmp_path / "primer3_core",
        output,
        capture_path=capture,
    )
    adapter = Primer3CLIDesignAdapter(executable)

    result = adapter.design(prepared)

    assert result.execution_performed
    assert result.status == "execution-complete-with-candidates"
    assert result.candidates[0].left_start == 20
    assert result.candidates[0].right_start == 120
    assert result.candidates[0].right_end == 140
    assert result.candidates[0].left_gc_fraction == pytest.approx(0.5)
    assert result.provenance.backend == adapter.info
    assert result.issues[0].code == "PRIMER3_DESIGN_WARNING"
    captured = capture.read_text(encoding="utf-8")
    assert "SEQUENCE_TARGET=50,50" in captured
    assert "SEQUENCE_EXCLUDED_REGION=10,10" in captured
    assert "PRIMER_NUM_RETURN=2" in captured
    assert "PRIMER_MIN_GC=40.0" in captured
    json.dumps(result.to_dict(), sort_keys=True)


def test_primer3_design_rejects_sequence_coordinate_mismatch(
    tmp_path: Path,
) -> None:
    prepared = prepare_primer_design(DNASequence("ACGT" * 100), target_start=50, target_end=100)
    output = "\n".join(
        (
            "PRIMER_PAIR_NUM_RETURNED=1",
            "PRIMER_LEFT_0=20,20",
            "PRIMER_RIGHT_0=139,20",
            f"PRIMER_LEFT_0_SEQUENCE={'ACGT' * 5}",
            f"PRIMER_RIGHT_0_SEQUENCE={'TGCA' * 5}",
            "PRIMER_LEFT_0_TM=60.1",
            "PRIMER_RIGHT_0_TM=60.2",
            "PRIMER_LEFT_0_GC_PERCENT=50.0",
            "PRIMER_RIGHT_0_GC_PERCENT=50.0",
            "PRIMER_PAIR_0_PRODUCT_SIZE=120",
            "=",
            "",
        )
    )
    adapter = Primer3CLIDesignAdapter(_fake_primer3_core(tmp_path / "primer3_core", output))

    with pytest.raises(BackendExecutionError) as error:
        adapter.design(prepared)
    assert error.value.code == "INVALID_PRIMER3_DESIGN_OUTPUT"


def test_primer3_design_adapter_rejects_backend_error_and_output_overflow(
    tmp_path: Path,
) -> None:
    prepared = prepare_primer_design(
        DNASequence("ACGT" * 100),
        target_start=50,
        target_end=100,
    )
    rejected = _fake_primer3_core(
        tmp_path / "primer3_core_rejected",
        "PRIMER_ERROR=bad\n=\n",
    )
    adapter = Primer3CLIDesignAdapter(rejected)
    with pytest.raises(BackendExecutionError, match="rejected"):
        adapter.design(prepared)

    overflow = _fake_primer3_core(
        tmp_path / "primer3_core_overflow",
        "PRIMER_PAIR_NUM_RETURNED=2\n=\n",
    )
    adapter = Primer3CLIDesignAdapter(overflow)
    with pytest.raises(BackendExecutionError, match="more candidates"):
        adapter.design(prepared, max_returned_candidates=1)
    with pytest.raises(ConfigurationError, match="max_template_length"):
        adapter.design(prepared, max_template_length=100)


def test_primer3_design_counts_actual_mapping_keys(
    tmp_path: Path,
) -> None:
    prepared = prepare_primer_design(
        DNASequence("ACGT" * 100),
        target_start=50,
        target_end=100,
    )
    output = "\n".join((*[f"KEY_{index}=0" for index in range(11)], "=", ""))
    adapter = Primer3CLIDesignAdapter(_fake_primer3_core(tmp_path / "primer3_core", output))

    with pytest.raises(BackendExecutionError) as error:
        adapter.design(prepared, max_result_keys=10)
    assert error.value.code == "PRIMER3_DESIGN_OUTPUT_LIMIT_EXCEEDED"

    trailing = _fake_primer3_core(
        tmp_path / "primer3_core_trailing",
        "PRIMER_PAIR_NUM_RETURNED=0\n=\nUNEXPECTED=1\n",
    )
    with pytest.raises(BackendExecutionError) as trailing_error:
        Primer3CLIDesignAdapter(trailing).design(prepared)
    assert trailing_error.value.code == "INVALID_PRIMER3_DESIGN_OUTPUT"


def test_crispr_reuses_pam_candidates_and_enumerates_sequence_only_matches() -> None:
    target = DNASequence("A" * 20 + "TGG")
    result = scan_crispr_candidates(
        target,
        "SpCas9",
        references={"ref": DNASequence("A" * 20 + "CCCC")},
    )

    assert result.candidates[0].guide_sequence == "A" * 20
    assert any(hit.reference_id == "ref" and hit.start == 0 for hit in result.off_targets)
    assert not result.efficiency_prediction_performed
    assert result.parameters["allow_indels_or_bulges"] is False


def test_crispr_without_reference_and_comparison_limit_are_audited() -> None:
    no_reference = scan_crispr_candidates(DNASequence("A" * 20 + "TGG"), "SpCas9")
    assert no_reference.off_targets == ()
    assert {issue.code for issue in no_reference.issues} >= {
        "CRISPR_REFERENCE_NOT_PROVIDED",
        "CRISPR_BIOLOGICAL_RISK_NOT_PREDICTED",
    }
    with pytest.raises(ConfigurationError, match="max_total_comparison_cells"):
        scan_crispr_candidates(
            DNASequence("A" * 20 + "TGG"),
            "SpCas9",
            references={"ref": DNASequence("A" * 100)},
            max_total_comparison_cells=1,
        )
