"""Tests for dnakit.regulation.fimo_p17_cli_functions (FIMO with position-specific priors).

Command assembly and result parsing are tested with a small stand-in ``fimo`` script written into tmp_path (it records its
arguments and writes the output files FIMO writes). Tests marked ``requires_fimo`` run the real program and skip when
the MEME suite is not installed. No test uses the network.
"""

from __future__ import annotations

import inspect
import json
import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

import dnakit.regulation.fimo_p17_cli_functions as fimo_module
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.regulation.fimo_p17_cli_functions import scan_motifs_with_position_priors

requires_fimo = pytest.mark.skipif(
    shutil.which("fimo") is None or shutil.which("compute-prior-dist") is None,
    reason="the MEME suite (fimo, compute-prior-dist) is not on PATH",
)

MOTIF = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF TESTMOTIF
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.05 0.05 0.85 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.05 0.85
0.05 0.85 0.05 0.05
0.85 0.05 0.05 0.05
0.05 0.05 0.85 0.05
"""

FAKE_FIMO = """#!{python}
import json, sys
from pathlib import Path

args = sys.argv[1:]
if args == ["--version"]:
    print("5.5.9")
    sys.exit(0)
out = Path(args[args.index("--oc") + 1])
(out / "args.json").write_text(json.dumps(args))
(out / "fimo.tsv").write_text(
    "motif_id\\tmotif_alt_id\\tsequence_name\\tstart\\tstop\\tstrand\\tscore\\tp-value\\tq-value\\tmatched_sequence\\n"
    "TESTMOTIF\\t\\tseq1\\t11\\t16\\t+\\t7.5\\t1.2e-05\\t0.03\\tGATCAG\\n"
    "TESTMOTIF\\t\\tseq1\\t40\\t45\\t-\\t6.1\\t4e-05\\t0.04\\tGATCAG\\n"
    "\\n# FIMO (Find Individual Motif Occurrences): Version 5.5.9\\n"
)
(out / "fimo.gff").write_text("##gff-version 3\\n")
"""


def _fake_fimo(tmp_path: Path) -> str:
    path = tmp_path / "fake_bin" / "fimo"
    path.parent.mkdir()
    path.write_text(FAKE_FIMO.replace("{python}", sys.executable), encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return str(path)


def _inputs(tmp_path: Path) -> dict[str, str]:
    files = {
        "motif_file": ("motif.meme", MOTIF),
        "sequence_file": ("seq.fa", ">seq1\nACGTACGTACGATCAGACGTACGT\n"),
        "prior_file": ("prior.wig", "variableStep chrom=seq1\n1 0.5\n"),
        "prior_distribution_file": ("prior.dist", "#min 0.5\n#max 0.5\n"),
    }
    for filename, text in files.values():
        (tmp_path / filename).write_text(text, encoding="utf-8")
    return {key: str(tmp_path / filename) for key, (filename, _) in files.items()}


def test_public_function_has_keyword_only_verbose() -> None:
    assert fimo_module.__all__ == ["scan_motifs_with_position_priors"]
    parameter = inspect.signature(scan_motifs_with_position_priors).parameters["verbose"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is True


def test_command_assembly_and_parsed_matches(tmp_path: Path) -> None:
    files = _inputs(tmp_path)
    out = tmp_path / "out"
    result = scan_motifs_with_position_priors(
        **files,
        output_dir=str(out),
        prior_alpha=0.5,
        p_value_threshold=1e-3,
        scan_reverse_complement=False,
        motif_ids=["TESTMOTIF"],
        background_source="uniform",
        fimo_path=_fake_fimo(tmp_path),
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["version"] == "5.5.9"
    arguments = json.loads((out / "args.json").read_text())
    assert arguments[arguments.index("--psp") + 1] == files["prior_file"]
    assert arguments[arguments.index("--prior-dist") + 1] == files["prior_distribution_file"]
    assert arguments[arguments.index("--alpha") + 1] == "0.5"
    assert arguments[arguments.index("--thresh") + 1] == "0.001"
    assert arguments[arguments.index("--bfile") + 1] == "--uniform--"
    assert arguments[arguments.index("--motif") + 1] == "TESTMOTIF"
    assert "--norc" in arguments
    assert arguments[-2:] == [files["motif_file"], files["sequence_file"]]
    assert result["match_count"] == 2
    assert result["matches_truncated"] is False
    first = result["matches"][0]
    assert (first["motif_id"], first["start"], first["stop"], first["strand"]) == ("TESTMOTIF", 11, 16, "+")
    assert first["score"] == 7.5 and first["p_value"] == 1.2e-05 and first["q_value"] == 0.03
    assert first["matched_sequence"] == "GATCAG"
    assert result["tsv_file"] == str(out / "fimo.tsv")
    json.dumps(result)


def test_max_rows_truncates_the_parsed_rows_only(tmp_path: Path) -> None:
    result = scan_motifs_with_position_priors(
        **_inputs(tmp_path),
        output_dir=str(tmp_path / "out"),
        max_rows=1,
        fimo_path=_fake_fimo(tmp_path),
        verbose=False,
    )
    assert result["match_count"] == 2
    assert len(result["matches"]) == 1
    assert result["matches_truncated"] is True


def test_q_value_options_are_passed_or_refused(tmp_path: Path) -> None:
    files = _inputs(tmp_path)
    fake = _fake_fimo(tmp_path)
    result = scan_motifs_with_position_priors(
        **files, output_dir=str(tmp_path / "q"), use_q_value_threshold=True, max_strand=True,
        parse_genomic_coordinates=False, fimo_path=fake, verbose=False,
    )
    arguments = json.loads((tmp_path / "q" / "args.json").read_text())
    assert result["status"] == "completed"
    assert {"--qv-thresh", "--max-strand", "--no-pgc"} <= set(arguments)
    with pytest.raises(ConfigurationError, match="compute_q_values"):
        scan_motifs_with_position_priors(
            **files, output_dir=str(tmp_path / "bad"), use_q_value_threshold=True, compute_q_values=False,
            fimo_path=fake, verbose=False,
        )
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"prior_alpha": 0.0}, "prior_alpha"),
        ({"prior_alpha": 1.5}, "prior_alpha"),
        ({"p_value_threshold": -1.0}, "p_value_threshold"),
        ({"background_source": "genome"}, "background_source"),
        ({"background_source": "uniform", "background_file": "x"}, "either"),
        ({"motif_ids": []}, "motif_ids"),
        ({"motif_ids": ["--x"]}, "motif_ids"),
        ({"max_stored_scores": 0}, "max_stored_scores"),
        ({"prior_file": "/nonexistent/prior.wig"}, "prior_file"),
    ],
)
def test_invalid_arguments_are_refused_before_anything_runs(tmp_path: Path, override: dict, message: str) -> None:
    files = _inputs(tmp_path) | override
    with pytest.raises(ConfigurationError, match=message):
        scan_motifs_with_position_priors(
            **files, output_dir=str(tmp_path / "out"), fimo_path=_fake_fimo(tmp_path), verbose=False
        )
    assert not (tmp_path / "out").exists()


def test_existing_output_directory_is_never_overwritten(tmp_path: Path) -> None:
    (tmp_path / "out").mkdir()
    with pytest.raises(ConfigurationError, match="must not exist"):
        scan_motifs_with_position_priors(
            **_inputs(tmp_path), output_dir=str(tmp_path / "out"), fimo_path=_fake_fimo(tmp_path), verbose=False
        )


def test_missing_fimo_raises_install_hint_and_creates_nothing(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python"))
    with pytest.raises(BackendUnavailableError) as raised:
        scan_motifs_with_position_priors(**_inputs(tmp_path), output_dir=str(tmp_path / "out"), verbose=False)
    assert "meme" in (raised.value.hint or "").lower()
    assert not (tmp_path / "out").exists()
    with pytest.raises(BackendUnavailableError):
        scan_motifs_with_position_priors(
            **_inputs(tmp_path), output_dir=str(tmp_path / "out2"), fimo_path=str(tmp_path / "nope"), verbose=False
        )


# ---------------------------------------------------------------------------
# Real FIMO (skipped without the MEME suite)
# ---------------------------------------------------------------------------


def _real_inputs(tmp_path: Path) -> dict[str, str]:
    """A 600 bp sequence with GATCAG planted at 101 and 401 (1-based); a wiggle prior that favours only site 101.

    ``compute-prior-dist`` reads MEME PSP format, so the same priors are also written as a PSP file to make the
    binned distribution that FIMO needs next to the wiggle prior.
    """
    import subprocess

    rng = random.Random(11)
    sequence = [rng.choice("ACGT") for _ in range(600)]
    for start in (100, 400):  # 0-based
        sequence[start : start + 6] = "GATCAG"
    (tmp_path / "seq.fa").write_text(">seq1\n" + "".join(sequence) + "\n", encoding="utf-8")
    (tmp_path / "motif.meme").write_text(MOTIF, encoding="utf-8")
    priors = [0.9 if abs(position - 101) <= 5 else 0.1 for position in range(1, 601)]
    total = sum(priors)
    priors = [value / total for value in priors]
    priors[-5:] = [0.0] * 5  # a width-6 motif cannot start in the last five positions
    wiggle = ["variableStep chrom=seq1 span=1"] + [f"{i} {v:.8f}" for i, v in enumerate(priors, start=1) if v > 0]
    (tmp_path / "prior.wig").write_text("\n".join(wiggle) + "\n", encoding="utf-8")
    (tmp_path / "prior.psp").write_text(">seq1 6\n" + "\n".join(f"{v:.8f}" for v in priors) + "\n", encoding="utf-8")
    distribution = subprocess.run(
        ["compute-prior-dist", "10", str(tmp_path / "prior.psp")], capture_output=True, text=True, check=True
    )
    (tmp_path / "prior.dist").write_text(distribution.stdout, encoding="utf-8")
    return {
        "motif_file": str(tmp_path / "motif.meme"),
        "sequence_file": str(tmp_path / "seq.fa"),
        "prior_file": str(tmp_path / "prior.wig"),
        "prior_distribution_file": str(tmp_path / "prior.dist"),
    }


@requires_fimo
def test_real_fimo_scores_with_priors_and_matches_the_program(tmp_path: Path) -> None:
    import subprocess

    files = _real_inputs(tmp_path)
    result = scan_motifs_with_position_priors(
        **files, output_dir=str(tmp_path / "out"), p_value_threshold=1e-3, verbose=False
    )
    assert result["status"] == "completed", result["log_tail"]
    assert (tmp_path / "out" / "fimo.gff").is_file()
    by_start = {row["start"]: row for row in result["matches"] if row["strand"] == "+"}
    assert set(by_start) == {101, 401}
    # the prior raises the posterior score and lowers the p-value of the favoured site only
    assert by_start[101]["score"] > by_start[401]["score"]
    assert by_start[101]["p_value"] < by_start[401]["p_value"]
    direct = subprocess.run(
        ["fimo", "--text", "--thresh", "0.001", "--psp", files["prior_file"], "--prior-dist",
         files["prior_distribution_file"], files["motif_file"], files["sequence_file"]],
        capture_output=True, text=True, check=True,
    ).stdout
    rows = [line.split("\t") for line in direct.splitlines() if line.startswith("TESTMOTIF")]
    assert sorted((int(r[3]), float(r[6])) for r in rows) == sorted(
        (row["start"], row["score"]) for row in result["matches"]
    )


@requires_fimo
def test_real_fimo_default_threshold_keeps_only_the_favoured_site(tmp_path: Path) -> None:
    result = scan_motifs_with_position_priors(**_real_inputs(tmp_path), output_dir=str(tmp_path / "out"), verbose=False)
    assert [(row["start"], row["strand"]) for row in result["matches"]] == [(101, "+")]
    assert result["matches"][0]["matched_sequence"] == "GATCAG"
    assert result["matches"][0]["q_value"] is not None


@requires_fimo
def test_real_fimo_rejects_prior_without_distribution_at_the_program(tmp_path: Path) -> None:
    files = _real_inputs(tmp_path)
    files["prior_distribution_file"] = str(tmp_path / "empty.dist")
    (tmp_path / "empty.dist").write_text("", encoding="utf-8")
    result = scan_motifs_with_position_priors(**files, output_dir=str(tmp_path / "out"), verbose=False)
    assert result["status"] == "failed" and result["matches"] == []
