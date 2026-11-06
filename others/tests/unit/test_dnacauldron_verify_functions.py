"""Tests for dnacauldron_verify_functions.simulate_dnacauldron_plan.

A circular Gibson plan of three parts is written into tmp_path; the simulation must be valid and the
folder report must contain the plan summary. The same plan run through the upstream dnacauldron command
line gives an identical report (checked once by hand: ten files, no differing bytes). Ordinary tests
never reach the network.
"""

from __future__ import annotations

import importlib.util
import random
import socket

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import dnacauldron_verify_functions as module

requires_dnacauldron = pytest.mark.skipif(
    importlib.util.find_spec("dnacauldron") is None, reason="dnacauldron is not installed"
)


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)


def _write_plan(root):
    rng = random.Random(3)

    def dna(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    o0, o1, o2 = dna(40), dna(40), dna(40)
    parts = root / "parts"
    parts.mkdir()
    records = {
        "partA": o2 + dna(160) + o0,
        "partB": o0 + dna(160) + o1,
        "partC": o1 + dna(160) + o2,
    }
    for name, sequence in records.items():
        (parts / f"{name}.fa").write_text(f">{name}\n{sequence}\n")
    plan = root / "plan.csv"
    plan.write_text("class: gibson_assembly,construct_1,partA,partB,partC\n")
    connectors = root / "connectors"
    connectors.mkdir()
    return plan, parts, connectors


@requires_dnacauldron
def test_valid_gibson_plan_writes_the_folder_report(tmp_path):
    plan, parts, connectors = _write_plan(tmp_path)
    target = tmp_path / "report"
    result = module.simulate_dnacauldron_plan(
        plan, parts, target, connectors_path=connectors, verbose=False
    )
    assert result["outcomes"] == {
        "cancelled_assemblies": 0,
        "errored_assemblies": 0,
        "valid_assemblies": 1,
    }
    assert result["assemblies"] == 1
    assert (target / "plan_simulation" / "plan_summary.csv").is_file()
    assert (target / "plan_simulation" / "construct_1" / "construct_1.gb").is_file()


@requires_dnacauldron
def test_existing_report_is_never_overwritten(tmp_path):
    plan, parts, _ = _write_plan(tmp_path)
    target = tmp_path / "taken"
    target.mkdir()
    with pytest.raises(ConfigurationError):
        module.simulate_dnacauldron_plan(plan, parts, target, verbose=False)


@requires_dnacauldron
def test_missing_plan_is_rejected(tmp_path):
    _, parts, _ = _write_plan(tmp_path)
    with pytest.raises(ConfigurationError):
        module.simulate_dnacauldron_plan(
            tmp_path / "none.csv", parts, tmp_path / "r", verbose=False
        )
