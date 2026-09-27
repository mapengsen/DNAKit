"""Initial ANI tests: real scientific sanity checks and strict input/output contracts."""

from __future__ import annotations

import gzip
import json
import random
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.genomics import ani
from dnakit.genomics.ani import average_nucleotide_identity


@pytest.fixture(scope="module")
def genome() -> str:
    # Independent random DNA has no organism origin or external download.
    return "".join(random.Random(1729).choices("ACGT", k=300_000))


@pytest.fixture
def real_pyfastani() -> None:
    pytest.importorskip("pyfastani", minversion="0.6.2")


@pytest.mark.usefixtures("real_pyfastani")
def test_identical_genomes_and_reverse_complement(genome: str, capsys: Any) -> None:
    reverse = genome.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for query in (genome.lower(), reverse):
        result = average_nucleotide_identity(query, genome, threads=2, verbose=False)
        assert result["status"] == "matched"
        assert result["ani_percent"] == pytest.approx(100.0, abs=0.01)
        assert result["query_fragments"] == 100
        assert result["matching_fragments"] >= 95
        assert 0.95 <= result["query_aligned_fraction"] <= 1.0
        assert result["backend"] == "pyfastani"
        assert result["backend_version"]
        assert result["direction"] == "query_to_reference"
        assert result["query"]["contigs"] == 1
        assert result["query"]["length"] == len(genome)
        assert result["parameters"]["fragment_length"] == 3000
        assert json.loads(json.dumps(result, allow_nan=False)) == result
    captured = capsys.readouterr()
    assert not captured.out and not captured.err


@pytest.mark.usefixtures("real_pyfastani")
def test_unrelated_genomes_have_no_reported_identity(genome: str) -> None:
    other = "".join(random.Random(2718).choices("ACGT", k=len(genome)))
    result = average_nucleotide_identity(genome, other, verbose=False)
    assert result["status"] == "no_match"
    assert result["reason"] == "no_reportable_ani"
    assert result["ani_percent"] is None
    assert result["matching_fragments"] is None
    assert result["query_aligned_fraction"] is None
    assert result["query_fragments"] == 100
    json.dumps(result, allow_nan=False)


@pytest.mark.usefixtures("real_pyfastani")
def test_multicontig_gzip_fasta_preserves_boundaries(tmp_path: Path, genome: str) -> None:
    source = tmp_path / "draft.fasta.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(f">first\n{genome[:150_099]}\n>second\n{genome[150_099:]}\n")
    result = average_nucleotide_identity(source, str(source), verbose=False)
    assert result["status"] == "matched"
    assert result["ani_percent"] == pytest.approx(100.0, abs=0.01)
    assert result["query"]["contigs"] == result["reference"]["contigs"] == 2
    # Joining these contigs would incorrectly produce 100 complete fragments.
    assert result["query_fragments"] == 99
    assert result["query"]["source"] == str(source)
    unbroken = average_nucleotide_identity(genome, genome, verbose=False)
    assert result["query"]["input_sha256"] != unbroken["query"]["input_sha256"]


@pytest.mark.usefixtures("real_pyfastani")
def test_swapping_genomes_changes_query_fragment_denominator(genome: str) -> None:
    forward = average_nucleotide_identity(genome[:150_000], genome, verbose=False)
    reverse = average_nucleotide_identity(genome, genome[:150_000], verbose=False)
    assert forward["status"] == reverse["status"] == "matched"
    assert forward["query_fragments"] == 50
    assert reverse["query_fragments"] == 100
    assert forward["ani_percent"] == pytest.approx(100.0, abs=0.01)
    assert reverse["ani_percent"] == pytest.approx(100.0, abs=0.01)
    assert forward["query_aligned_fraction"] > reverse["query_aligned_fraction"]


@pytest.mark.usefixtures("real_pyfastani")
def test_short_query_and_explicit_alphabet_only_path(tmp_path: Path, genome: str) -> None:
    source = tmp_path / "ACGT"
    source.write_text(f">contig\n{genome}\n")
    result = average_nucleotide_identity("ACGTRYN", source, verbose=False)
    assert result["status"] == "no_match"
    assert result["reason"] == "no_complete_query_fragments"
    assert result["query_fragments"] == 0
    assert result["query"]["ambiguous_bases_to_n"] == 3
    assert result["reference"]["source"] == str(source)
    assert result["reference"]["length"] == len(genome)


@pytest.mark.usefixtures("real_pyfastani")
def test_uninformative_reference_is_no_match(genome: str) -> None:
    result = average_nucleotide_identity(genome, "N" * len(genome), verbose=False)
    assert result["status"] == "no_match"
    assert result["reason"] == "no_reference_minimizers"
    assert result["ani_percent"] is None
    assert result["reference"]["ambiguous_bases_to_n"] == len(genome)


@pytest.mark.usefixtures("real_pyfastani")
def test_all_short_contigs_are_not_joined(tmp_path: Path, genome: str) -> None:
    source = tmp_path / "short.fa"
    source.write_text(">one\n" + genome[:2000] + "\n>two\n" + genome[2000:4000] + "\n")
    result = average_nucleotide_identity(source, genome, verbose=False)
    assert result["query"]["length"] == 4000
    assert result["query_fragments"] == 0
    assert result["status"] == "no_match"
    assert result["reason"] == "no_complete_query_fragments"


@pytest.mark.parametrize("value", ["", "AC-GT", "AC GT", "ACUGT", ">seq\nACGT", 42, None])
def test_invalid_input_does_not_reach_backend(value: Any, monkeypatch: Any) -> None:
    def unexpected_import(*args: Any) -> None:
        pytest.fail("Invalid input reached the optional backend import.")

    monkeypatch.setattr(ani, "require", unexpected_import)
    with pytest.raises(ConfigurationError):
        average_nucleotide_identity(value, "ACGT", verbose=False)
    with pytest.raises(ConfigurationError):
        average_nucleotide_identity("ACGT", value, verbose=False)


@pytest.mark.parametrize("contents", ["", "ACGT\n", ">seq\nACXGT\n", ">seq\n"])
def test_invalid_fasta_rejected(tmp_path: Path, contents: str) -> None:
    source = tmp_path / "invalid.fa"
    source.write_text(contents)
    with pytest.raises((ConfigurationError, InputFormatError)):
        average_nucleotide_identity(source, "ACGT", verbose=False)


@pytest.mark.parametrize("threads", [True, 0, -1, 1.5, "2", 1025, 2**31])
def test_invalid_threads_rejected(threads: Any) -> None:
    with pytest.raises(ConfigurationError, match="threads"):
        average_nucleotide_identity("ACGT", "ACGT", threads=threads, verbose=False)


@pytest.mark.parametrize("identity", [float("nan"), float("inf"), -1.0, 100.1])
def test_invalid_backend_estimate_raises(identity: float, monkeypatch: Any) -> None:
    # A scientific backend must not leak NaN, infinity, or out-of-range ANI
    # into a successful JSON response. This failure is injected deliberately.
    hit = SimpleNamespace(name="reference", identity=identity, matches=1, fragments=1)
    mapper = SimpleNamespace(
        query_draft=lambda *args, **kwargs: [hit],
        minimizers=[object()],
        fragment_length=3000,
        k=16,
        minimum_fraction=0.2,
        percentage_identity=80,
    )
    sketch = SimpleNamespace(add_draft=lambda *args: None, index=lambda: mapper)
    backend = SimpleNamespace(Sketch=lambda: sketch, __version__="test")
    monkeypatch.setattr(ani, "require", lambda *args: backend)
    with pytest.raises(BackendExecutionError, match="invalid ANI"):
        average_nucleotide_identity("ACGT" * 750, "ACGT" * 750, verbose=False)
