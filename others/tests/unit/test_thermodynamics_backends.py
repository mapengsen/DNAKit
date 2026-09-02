"""Tests for explicit, bounded Primer3 CLI execution."""

from __future__ import annotations

import stat
import sys
from pathlib import Path

import pytest

from dnakit.core import DNASequence
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.thermodynamics import (
    Primer3CLIAdapter,
    ThermodynamicConditions,
    ThermodynamicsBackend,
    conditional_capabilities,
    probe_primer3,
)


def _python_executable(path: Path, body: str) -> Path:
    path.write_text(f"#!{sys.executable}\n{body}\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_missing_primer3_probe_is_passive_and_has_actionable_error() -> None:
    adapter = Primer3CLIAdapter(ntthal_path="/definitely/missing/ntthal")
    info = adapter.info

    assert isinstance(adapter, ThermodynamicsBackend)
    assert adapter.supports("hairpin")
    assert not adapter.supports("secondary_structure")
    assert not info.available
    assert info.version is None
    assert info.metadata["version_probe_executed"] is False
    assert info.metadata["automatic_path_search"] is False
    with pytest.raises(BackendUnavailableError) as error:
        adapter.ensure_available("hairpin")
    assert error.value.code == "PRIMER3_CLI_UNAVAILABLE"


def test_existing_executable_is_not_run_during_construction_or_probe(tmp_path: Path) -> None:
    marker = tmp_path / "executed"
    executable = _python_executable(
        tmp_path / "ntthal",
        f"from pathlib import Path\nPath({str(marker)!r}).touch()\nprint('invalid')",
    )
    adapter = Primer3CLIAdapter(ntthal_path=executable)

    assert adapter.info.available
    assert adapter.info.executable_path == str(executable.resolve())
    assert not marker.exists()
    assert probe_primer3(ntthal_path=executable).available
    assert not marker.exists()
    with pytest.raises(BackendExecutionError) as error:
        adapter.hairpin(DNASequence("ACGT"))
    assert error.value.code == "INVALID_PRIMER3_CLI_OUTPUT"
    assert marker.exists()


def test_cli_adapter_parses_oligotm_and_ntthal_outputs(tmp_path: Path) -> None:
    oligotm = _python_executable(
        tmp_path / "oligotm",
        (
            "import sys\n"
            "assert sys.argv[sys.argv.index('-tp') + 1] == '1'\n"
            "assert sys.argv[sys.argv.index('-sc') + 1] == '1'\n"
            "print('42.500000')"
        ),
    )
    ntthal = _python_executable(
        tmp_path / "ntthal",
        (
            "print('Calculated thermodynamical parameters for dimer:\\t'"
            "      'dS = -100\\tdH = -30000\\tdG = -2500\\tt = 55')\n"
            "print('SEQ\\tACGT')\nprint('STR\\tTGCA')"
        ),
    )
    adapter = Primer3CLIAdapter(oligotm_path=oligotm, ntthal_path=ntthal)
    conditions = ThermodynamicConditions(strand_concentration_molar=50e-9)

    tm = adapter.tm(DNASequence("ACGT"), conditions=conditions)
    hairpin = adapter.hairpin(DNASequence("ACGT"), conditions=conditions, output_structure=True)
    dimer = adapter.self_dimer(DNASequence("ACGT"), conditions=conditions)
    heterodimer = adapter.heterodimer(
        DNASequence("ACGT"), DNASequence("TGCA"), conditions=conditions
    )

    assert tm.tm_celsius == 42.5
    assert tm.method == "primer3-cli-oligotm-santalucia"
    assert hairpin.structure_found is True
    assert hairpin.delta_h_kcal_per_mol == -30.0
    assert hairpin.delta_g_kcal_per_mol == -2.5
    assert hairpin.delta_s_cal_per_k_mol == -100.0
    assert hairpin.ascii_structure == "SEQ\tACGT\nSTR\tTGCA"
    assert dimer.capability == "self_dimer"
    assert heterodimer.sequences_5to3 == ("ACGT", "TGCA")
    assert tm.backend.name == "primer3-cli"
    assert tm.backend.license_expression == "GPL-2.0-or-later"
    assert tm.provenance.backend == tm.backend


def test_cli_adapter_parses_no_structure_output(tmp_path: Path) -> None:
    ntthal = _python_executable(
        tmp_path / "ntthal",
        "print('0\\tdS = inf\\tdH = inf\\tinf\\tinf')",
    )
    result = Primer3CLIAdapter(ntthal_path=ntthal).hairpin(DNASequence("AAAA"))

    assert result.structure_found is False
    assert result.tm_celsius == 0.0
    assert result.delta_g_kcal_per_mol == 0.0
    assert result.issues[0].code == "PRIMER3_NO_STRUCTURE"


def test_conditional_metadata_uses_only_explicit_ntthal_path(tmp_path: Path) -> None:
    ntthal = _python_executable(tmp_path / "ntthal", "raise SystemExit(99)")
    capabilities = conditional_capabilities(primer3_ntthal_path=ntthal)

    assert [item.requirement_id for item in capabilities] == [
        "THERMO-008",
        "THERMO-009",
        "THERMO-010",
        "THERMO-011",
    ]
    assert all(item.execution_supported for item in capabilities[:3])
    assert all(not item.automatic_install for item in capabilities)
    assert all(not item.automatic_download for item in capabilities)
    assert capabilities[0].backend_info is not None
    assert capabilities[0].backend_info.name == "primer3-cli"
    assert capabilities[-1].backend_info is not None
    assert capabilities[-1].backend_info.name == "nupack"
    assert "separately licensed" in capabilities[-1].reason


def test_cli_adapter_validates_paths_conditions_and_boundaries(tmp_path: Path) -> None:
    ntthal = _python_executable(tmp_path / "ntthal", "print('invalid')")
    oligotm = _python_executable(tmp_path / "oligotm", "print('1')")

    with pytest.raises(ConfigurationError):
        Primer3CLIAdapter(ntthal_path="ntthal")
    with pytest.raises(BackendExecutionError) as unsupported:
        Primer3CLIAdapter().ensure_available("secondary_structure")
    assert unsupported.value.code == "BACKEND_CAPABILITY_UNSUPPORTED"
    with pytest.raises(ConfigurationError):
        Primer3CLIAdapter(ntthal_path=ntthal).hairpin(DNASequence("A"), max_loop=0)
    with pytest.raises(ConfigurationError):
        Primer3CLIAdapter(oligotm_path=oligotm).tm(DNASequence("A"))
    with pytest.raises(ConfigurationError) as cosolvent:
        Primer3CLIAdapter(ntthal_path=ntthal).hairpin(
            DNASequence("ACGT"),
            conditions=ThermodynamicConditions(dmso_percent=1.0),
        )
    assert cosolvent.value.code == "PRIMER3_NTTHAL_COSOLVENT_UNSUPPORTED"

    invalid_parameters = tmp_path / "missing-parameters"
    tm_only = Primer3CLIAdapter(
        oligotm_path=oligotm,
        thermodynamic_parameters_path=invalid_parameters,
    )
    assert tm_only.info.capabilities == frozenset({"tm"})
    assert tm_only.tm(DNASequence("ACGT")).tm_celsius == 1.0


def test_cli_adapter_wraps_nonzero_exit(tmp_path: Path) -> None:
    executable = _python_executable(
        tmp_path / "ntthal",
        "print('backend failed')\nraise SystemExit(7)",
    )
    with pytest.raises(BackendExecutionError) as error:
        Primer3CLIAdapter(ntthal_path=executable).hairpin(DNASequence("ACGT"))
    assert error.value.code == "PRIMER3_CLI_EXECUTION_FAILED"
