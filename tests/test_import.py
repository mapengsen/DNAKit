"""Smoke tests for package import and the minimal CLI."""

from typer.testing import CliRunner

import dnakit
from dnakit.cli.app import app

runner = CliRunner()


def test_import_exposes_version() -> None:
    assert dnakit.__version__ == "0.1.0.dev0"


def test_import_exposes_core_and_standardization_api() -> None:
    result = dnakit.normalize(" acgu ")

    assert result.sequence == dnakit.DNASequence("ACG")
    assert dnakit.validate(result.sequence).is_valid


def test_import_exposes_stable_io_api() -> None:
    """The architecture's lightweight I/O entry points live at package root."""

    assert callable(dnakit.read)
    assert callable(dnakit.read_one)
    assert callable(dnakit.read_set)
    assert callable(dnakit.write)
    assert dnakit.ReadConfig.__module__ == "dnakit.io.config"
    assert dnakit.WriteResult.__module__ == "dnakit.io.results"


def test_cli_info() -> None:
    result = runner.invoke(app, ["info"])

    assert result.exit_code == 0, result.output
    assert "DNAKit runtime" in result.output
    assert dnakit.__version__ in result.output


def test_cli_backends_does_not_require_external_programs() -> None:
    result = runner.invoke(app, ["backends"])

    assert result.exit_code == 0, result.output
    assert "primer3-cli" in result.output
    assert "primer3-py" not in result.output
    assert "nupack" in result.output
    assert result.output.count("not probed") == 6
