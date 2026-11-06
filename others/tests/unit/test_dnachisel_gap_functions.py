"""Tests for dnachisel_gap_functions: constraint breach plots and PDFs, optimization and no-solution
reports, optimization change plots, BLAST+ hits and Bowtie 1 index and matches.

Only the functions behind the batch D worklist rows are tested here. DnaChisel's reports need the
optional pdf_reports, sequenticon and geneblocks packages; their happy paths are skipped with the
reason when those are missing, while the guard path always runs. The BLAST+ and Bowtie 1 binaries
are not installed here, so their real-binary happy paths are skipped. Bowtie's index and match
pipeline is exercised with small stand-in executables written under tmp_path, which return the
output format the wrapper parses. A socket guard blocks any connection.
"""

from __future__ import annotations

import importlib.util
import shutil
import socket
import stat
import sys
import types
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import dnachisel_gap_functions as module

requires_pkg = pytest.mark.skipif(
    importlib.util.find_spec("dnachisel") is None,
    reason="dnachisel is not installed",
)


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _is_png(path: Path) -> bool:
    return path.is_file() and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def _is_pdf(path: Path) -> bool:
    return path.is_file() and path.read_bytes()[:5] == b"%PDF-"


def _installed(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


# 100 bases: the 50-base unit repeated, with GGATCC at 31 and 81.
UNIT = "ATGCGTACGTTAGCATCGATCGGCTAGCTAGGGATCCAAAGCTTGCATGC"
SEQUENCE = UNIT * 2
GC_CONSTRAINTS = [
    {"name": "EnforceGCContent", "mini": 0.45, "maxi": 0.55, "window": 20},
    {"name": "AvoidPattern", "pattern": "GGATCC"},
]
IMPOSSIBLE_GC = [
    {"name": "EnforceGCContent", "mini": 0.9, "maxi": 1.0, "window": 20},
    {"name": "AvoidChanges", "location": [0, 100]},
]


def _write_stand_in(path: Path, body: str) -> str:
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return str(path)


# ---------------------------------------------------------------------------
# Constraint breach plots and PDFs
# ---------------------------------------------------------------------------


@requires_pkg
def test_breach_plot_marks_the_avoided_pattern(tmp_path):
    target = tmp_path / "sub" / "breaches.png"
    result = module.plot_constraint_breaches_png(
        SEQUENCE, output_path=str(target), constraints=GC_CONSTRAINTS, verbose=False
    )
    assert _is_png(target)
    assert result["image_path"] == str(target.resolve())
    assert result["sequence_length"] == 100
    # The GGATCC pattern occurs at 31 and 81 on the plus strand (DnaChisel's breach location).
    assert SEQUENCE.find("GGATCC") == 31
    assert result["breaches"]["GGATCC"] == [[31, 37, 1], [81, 87, 1]]
    assert len(result["breaches"]) == 2


@requires_pkg
def test_breach_plot_refuses_existing_output_and_non_png(tmp_path):
    existing = tmp_path / "exists.png"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError, match="already exists"):
        module.plot_constraint_breaches_png(
            SEQUENCE, output_path=str(existing), constraints=GC_CONSTRAINTS, verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"\.png"):
        module.plot_constraint_breaches_png(
            SEQUENCE, output_path=str(tmp_path / "b.pdf"), constraints=GC_CONSTRAINTS, verbose=False
        )
    assert existing.read_bytes() == b"x"


@requires_pkg
def test_breach_pdf_has_one_page_per_sequence(tmp_path):
    target = tmp_path / "breaches.pdf"
    result = module.write_constraint_breaches_pdf(
        {"first": SEQUENCE, "second": "ATGC" * 10},
        output_path=str(target),
        constraints=GC_CONSTRAINTS,
        verbose=False,
    )
    assert _is_pdf(target)
    assert result["pdf_path"] == str(target.resolve())
    assert result["page_count"] == 2
    assert result["sequence_count"] == 2
    assert result["constraint_labels"] == ["GC outside 45-55%/20bp", "GGATCC"]
    names = [row["sequence"] for row in result["rows"]]
    assert names == ["first", "second"]


@requires_pkg
def test_breach_pdf_rejects_bad_specs_before_writing(tmp_path):
    with pytest.raises(ConfigurationError, match="Unsupported DnaChisel specification"):
        module.write_constraint_breaches_pdf(
            {"s": SEQUENCE},
            output_path=str(tmp_path / "x.pdf"),
            constraints=[{"name": "NotASpec"}],
            verbose=False,
        )
    assert not (tmp_path / "x.pdf").exists()


# ---------------------------------------------------------------------------
# Optimization and no-solution reports
# ---------------------------------------------------------------------------


@requires_pkg
def test_no_solution_report_writes_plots_for_impossible_constraints(tmp_path):
    target = tmp_path / "no_solution"
    result = module.write_dnachisel_no_solution_report(
        SEQUENCE, output_dir=str(target), constraints=IMPOSSIBLE_GC, verbose=False
    )
    assert result["status"] == "no_solution"
    assert result["report_path"] == str(target.resolve())
    assert "plots.pdf" in result["report_files"]
    assert _is_pdf(target / "plots.pdf")
    start, end, strand = result["location"]
    assert start == 0 and end > start and strand == 0
    assert "EnforceGCContent" in result["message"]


@requires_pkg
def test_no_solution_report_reports_satisfiable_constraints_without_writing(tmp_path):
    target = tmp_path / "unused"
    result = module.write_dnachisel_no_solution_report(
        SEQUENCE, output_dir=str(target), constraints=[GC_CONSTRAINTS[0]], verbose=False
    )
    assert result["status"] == "solved"
    assert result["report_path"] is None and result["report_files"] == []
    assert not target.exists()


@requires_pkg
def test_no_solution_report_needs_a_constraint_and_a_new_directory(tmp_path):
    with pytest.raises(ConfigurationError, match="at least one constraint"):
        module.write_dnachisel_no_solution_report(
            SEQUENCE, output_dir=str(tmp_path / "a"), verbose=False
        )
    (tmp_path / "taken").mkdir()
    with pytest.raises(ConfigurationError, match="already exists"):
        module.write_dnachisel_no_solution_report(
            SEQUENCE, output_dir=str(tmp_path / "taken"), constraints=GC_CONSTRAINTS, verbose=False
        )


@requires_pkg
def test_optimization_report_guard_without_report_packages(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "pdf_reports", None)
    with pytest.raises(BackendUnavailableError) as info:
        module.write_dnachisel_optimization_report(
            SEQUENCE, output_dir=str(tmp_path / "report"), constraints=GC_CONSTRAINTS, verbose=False
        )
    assert "pdf_reports" in str(info.value)
    assert not (tmp_path / "report").exists()


@pytest.mark.skipif(
    not (_installed("pdf_reports") and _installed("sequenticon")),
    reason="pdf_reports and sequenticon (extra bio-dnachisel-reports) are not installed here",
)
@requires_pkg
def test_optimization_report_happy_path(tmp_path):
    result = module.write_dnachisel_optimization_report(
        SEQUENCE, output_dir=str(tmp_path / "report"), constraints=GC_CONSTRAINTS, verbose=False
    )
    assert result["all_constraints_pass"] is True
    assert result["sequence_length"] == 100


@requires_pkg
def test_optimization_report_needs_constraints_or_objectives(tmp_path, monkeypatch):
    # The report packages are required before the specs are read; stand-ins let the check run.
    monkeypatch.setitem(sys.modules, "pdf_reports", types.ModuleType("pdf_reports"))
    monkeypatch.setitem(sys.modules, "sequenticon", types.ModuleType("sequenticon"))
    with pytest.raises(ConfigurationError, match="constraint or objective"):
        module.write_dnachisel_optimization_report(
            SEQUENCE, output_dir=str(tmp_path / "r"), verbose=False
        )
    assert not (tmp_path / "r").exists()


# ---------------------------------------------------------------------------
# Optimization change plot
# ---------------------------------------------------------------------------


@requires_pkg
def test_changes_plot_guard_without_geneblocks(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "geneblocks", None)
    with pytest.raises(BackendUnavailableError) as info:
        module.plot_optimization_changes_png(
            SEQUENCE,
            output_path=str(tmp_path / "changes.png"),
            constraints=GC_CONSTRAINTS,
            verbose=False,
        )
    assert "bio-geneblocks" in (info.value.hint or "")
    assert not (tmp_path / "changes.png").exists()


@pytest.mark.skipif(
    not _installed("geneblocks"),
    reason="geneblocks (extra bio-geneblocks) is not installed here",
)
@requires_pkg
def test_changes_plot_happy_path(tmp_path):
    result = module.plot_optimization_changes_png(
        SEQUENCE,
        output_path=str(tmp_path / "changes.png"),
        constraints=GC_CONSTRAINTS,
        verbose=False,
    )
    assert _is_png(tmp_path / "changes.png")
    assert result["all_constraints_pass"] is True


# ---------------------------------------------------------------------------
# BLAST+ hits (blastn)
# ---------------------------------------------------------------------------


def test_blast_argument_checks_run_before_the_binary_lookup(tmp_path):
    with pytest.raises(ConfigurationError, match="exactly one"):
        module.blast_sequence_hits(
            SEQUENCE, blast_db="db", subject_fasta=str(tmp_path / "s.fa"), verbose=False
        )
    with pytest.raises(ConfigurationError, match="task"):
        module.blast_sequence_hits(SEQUENCE, blast_db="db", task="tblastx", verbose=False)
    with pytest.raises(ConfigurationError, match="dust"):
        module.blast_sequence_hits(SEQUENCE, blast_db="db", dust="maybe", verbose=False)
    with pytest.raises(ConfigurationError, match="existing FASTA"):
        module.blast_sequence_hits(
            SEQUENCE, subject_fasta=str(tmp_path / "missing.fa"), verbose=False
        )


@requires_pkg
def test_blast_missing_binary_raises_backend_error(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(BackendUnavailableError) as info:
        module.blast_sequence_hits(SEQUENCE, blast_db=str(tmp_path / "db"), verbose=False)
    assert info.value.code == "DNACHISEL_EXTERNAL_UNAVAILABLE"
    assert "BLAST+" in (info.value.hint or "")
    explicit = tmp_path / "no-such-blastn"
    with pytest.raises(BackendUnavailableError, match="missing or not executable"):
        module.blast_sequence_hits(
            SEQUENCE, blast_db=str(tmp_path / "db"), blastn_path=str(explicit), verbose=False
        )


@pytest.mark.skipif(
    shutil.which("blastn") is None,
    reason="blastn (BLAST+ command-line tools) is not installed here; no real BLAST run",
)
def test_blast_happy_path_with_real_blastn(tmp_path):
    subject = tmp_path / "subject.fa"
    subject.write_text(f">subject\n{SEQUENCE}\n", encoding="utf-8")
    result = module.blast_sequence_hits(SEQUENCE, subject_fasta=str(subject), verbose=False)
    assert result["hit_count"] >= 1
    first = result["queries"][0]["alignments"][0]["hsps"][0]
    assert first["identities"] == 100 and first["query_start"] == 1


# ---------------------------------------------------------------------------
# Bowtie 1 index and matches
# ---------------------------------------------------------------------------


def test_bowtie_index_argument_checks(tmp_path):
    with pytest.raises(ConfigurationError, match="non-empty list"):
        module.create_bowtie_index([], index_dir=str(tmp_path / "idx"), verbose=False)
    with pytest.raises(ConfigurationError, match="sequences\\[0\\]"):
        module.create_bowtie_index(["ACGT!"], index_dir=str(tmp_path / "idx"), verbose=False)
    assert not (tmp_path / "idx").exists()


def test_bowtie_index_missing_binary_raises_backend_error(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(BackendUnavailableError, match="bowtie-build"):
        module.create_bowtie_index([SEQUENCE], index_dir=str(tmp_path / "idx"), verbose=False)


def test_bowtie_index_with_stand_in_builder(tmp_path):
    builder = _write_stand_in(
        tmp_path / "bowtie-build",
        'if [ "$1" = "--version" ]; then echo "bowtie-build version 1.3.1"; exit 0; fi\n'
        'touch "$3.1.ebwt"\n',
    )
    target = tmp_path / "index"
    result = module.create_bowtie_index(
        [SEQUENCE, "ACGTACGT"], index_dir=str(target), bowtie_build_path=builder, verbose=False
    )
    assert result["version"] == "1.3.1"
    assert result["sequence_count"] == 2
    assert result["index_prefix"] == str(target.resolve() / "bowtie")
    assert "bowtie.1.ebwt" in result["index_files"]
    fasta = (target / "sequences.fa").read_text(encoding="utf-8")
    assert fasta == f">0\n{SEQUENCE}\n>1\nACGTACGT\n"


def test_bowtie_match_argument_checks(tmp_path):
    with pytest.raises(ConfigurationError, match="No Bowtie index"):
        module.find_bowtie_matches(
            SEQUENCE, index_prefix=str(tmp_path / "none"), match_length=10, verbose=False
        )
    (tmp_path / "idx.1.ebwt").write_bytes(b"x")
    with pytest.raises(ConfigurationError, match=r"max_mismatches|between"):
        module.find_bowtie_matches(
            SEQUENCE,
            index_prefix=str(tmp_path / "idx"),
            match_length=10,
            max_mismatches=4,
            verbose=False,
        )


def test_bowtie_match_missing_binary_raises_backend_error(tmp_path, monkeypatch):
    (tmp_path / "idx.1.ebwt").write_bytes(b"x")
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(BackendUnavailableError, match="bowtie"):
        module.find_bowtie_matches(
            SEQUENCE, index_prefix=str(tmp_path / "idx"), match_length=10, verbose=False
        )


def test_bowtie_matches_with_stand_in_aligner(tmp_path):
    (tmp_path / "idx.1.ebwt").write_bytes(b"x")
    aligner = _write_stand_in(
        tmp_path / "bowtie",
        'if [ "$1" = "--version" ]; then echo "bowtie version 1.3.1"; exit 0; fi\n'
        # One line per k-mer: index, then the mismatch field (empty for an exact hit).
        "printf '0\\t\\n2\\t1:A>C\\n'\n",
    )
    result = module.find_bowtie_matches(
        "ACGTACGT",
        index_prefix=str(tmp_path / "idx"),
        match_length=4,
        max_mismatches=1,
        bowtie_path=aligner,
        verbose=False,
    )
    assert result["version"] == "1.3.1"
    assert result["kmer_length"] == 4
    assert result["matches"] == [
        {"start": 0, "end": 4, "mismatches": 0},
        {"start": 2, "end": 6, "mismatches": 1},
    ]
    assert result["match_count"] == 2


@pytest.mark.skipif(
    shutil.which("bowtie") is None or shutil.which("bowtie-build") is None,
    reason="bowtie and bowtie-build (Bowtie 1) are not installed here; no real index or alignment",
)
def test_bowtie_happy_path_with_real_binaries(tmp_path):
    index = module.create_bowtie_index([SEQUENCE], index_dir=str(tmp_path / "idx"), verbose=False)
    result = module.find_bowtie_matches(
        SEQUENCE, index_prefix=index["index_prefix"], match_length=20, verbose=False
    )
    # Every 20-mer of the indexed sequence is present exactly, so 81 matches with no mismatch.
    assert result["match_count"] == 81
    assert all(match["mismatches"] == 0 for match in result["matches"])
