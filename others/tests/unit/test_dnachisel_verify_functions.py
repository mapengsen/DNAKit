"""Tests for dnachisel_verify_functions.optimize_dnachisel_source_record.

The record carries its constraints as GenBank labels (``@AvoidPattern(GGTCTC)`` and an GC window), so the
DnaChisel command line can read them. The DNAKit function must give the same sequence as the command
line when both are seeded after their imports (random and numpy seeds). The comparison runs the command
line in a subprocess; no network is used.
"""

from __future__ import annotations

import importlib.util
import random
import socket
import subprocess
import sys

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import dnachisel_verify_functions as module

requires_dnachisel = pytest.mark.skipif(
    importlib.util.find_spec("dnachisel") is None or importlib.util.find_spec("Bio") is None,
    reason="dnachisel and biopython are not installed",
)

SEEDED_CLI = """
import random, sys
from dnachisel.cli import main
import numpy as np
random.seed(0); np.random.seed(0)
sys.argv = ['dnachisel', {source!r}, {target!r}]
main()
"""


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)


def _write_record(path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    rng = random.Random(7)
    sequence = "".join(rng.choice("ACGT") for _ in range(300))
    sequence = sequence[:50] + "GGTCTC" + sequence[56:150] + "GGTCTC" + sequence[156:]
    record = SeqRecord(
        Seq(sequence),
        id="fixture",
        name="fixture",
        description="fixture",
        annotations={"molecule_type": "DNA"},
    )
    for label in ("@AvoidPattern(GGTCTC)", "@EnforceGCContent(mini=0.4, maxi=0.6, window=50)"):
        record.features.append(
            SeqFeature(
                FeatureLocation(0, 300, 1), type="misc_feature", qualifiers={"label": [label]}
            )
        )
    SeqIO.write(record, str(path), "genbank")


@requires_dnachisel
def test_record_is_optimised_and_matches_the_command_line(tmp_path):
    from Bio import SeqIO

    source = tmp_path / "fixture.gb"
    _write_record(source)
    ours = tmp_path / "dnakit.gb"
    cli = tmp_path / "cli.gb"
    import dnachisel  # noqa: F401  - imported before seeding, as the command line does
    import numpy as np

    random.seed(0)
    np.random.seed(0)
    result = module.optimize_dnachisel_source_record(str(source), str(ours), verbose=False)
    assert result["output_kind"] == "genbank"
    subprocess.run(
        [sys.executable, "-c", SEEDED_CLI.format(source=str(source), target=str(cli))],
        check=True,
        capture_output=True,
        timeout=280,
    )
    ours_seq = str(SeqIO.read(str(ours), "genbank").seq)
    cli_seq = str(SeqIO.read(str(cli), "genbank").seq)
    assert "GGTCTC" not in ours_seq
    assert ours_seq == cli_seq


@requires_dnachisel
def test_existing_output_is_never_overwritten(tmp_path):
    source = tmp_path / "fixture.gb"
    _write_record(source)
    target = tmp_path / "taken.gb"
    target.write_text("occupied")
    with pytest.raises(ConfigurationError):
        module.optimize_dnachisel_source_record(str(source), str(target), verbose=False)
