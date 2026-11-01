"""Tests for dnakit.annotation.repeatmasker_p11_cli_functions; real runs skip when RepeatMasker (configured with RMBlast
and TRF) or its util scripts are missing. Programs that need FamDB or other unavailable data use a stand-in script."""

from __future__ import annotations

import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.annotation import repeatmasker_p11_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

have_rm = shutil.which("RepeatMasker") is not None and shutil.which("rmOutToGFF3.pl") is not None
requires_rm = pytest.mark.skipif(not have_rm, reason="RepeatMasker with its util scripts is not installed")

OUT = """   SW   perc perc perc  query     position in query    matching repeat       position in repeat
score   div. del. ins.  sequence  begin end   (left)   repeat   class/family begin  end    (left)  ID

  734   27.6  9.9  4.0  seq1       1483  1806  (326) + L1M5     LINE/L1      5368   5710  (414)   1
 2134   10.7  3.3  1.3  seq1       1809  2111   (21) + AluSx    SINE/Alu        3    311    (1)   2
"""


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    rng = random.Random(5)
    root = tmp_path_factory.mktemp("rm")
    repeat = "".join(rng.choice("ACGT") for _ in range(300))
    seq = ""
    for _ in range(4):
        seq += "".join(rng.choice("ACGT") for _ in range(900)) + repeat
    seq += "".join(rng.choice("ACGT") for _ in range(900))
    fa = root / "genome.fa"
    fa.write_text(f">seq1\n{seq}\n")
    lib = root / "lib.fa"
    lib.write_text(f">Rep1#SINE/Alu\n{repeat}\n")
    out = root / "x.out"
    out.write_text(OUT)
    seqfile = root / "seq1.fa"
    seqfile.write_text(">seq1\n" + "ACGT" * 700 + "\n")
    return {"root": root, "fa": fa, "lib": lib, "out": out, "seq": seqfile}


def test_missing_program(tmp_path: Path, work: dict[str, Path]) -> None:
    with pytest.raises(BackendUnavailableError):
        m.convert_repeats_to_gff3(work["out"], tmp_path / "o", repeatmasker_dir=str(tmp_path))
    with pytest.raises(ConfigurationError):
        m.mask_repeats(work["fa"], tmp_path / "o2", species="x", library=work["lib"], repeatmasker_dir=str(tmp_path), verbose=False)


@requires_rm
def test_mask_process_and_tools(work: dict[str, Path]) -> None:
    root = work["root"]
    run = m.mask_repeats(work["fa"], root / "run", library=work["lib"], alignments=True, mask_style="soft", verbose=False)
    assert run["annotations"] == 4 and run["summary"]["masked_bases"] >= 1100 and run["version"]
    out, cat, aln = (run["outputs"][f"genome.fa.{k}"] for k in ("out", "cat", "align"))
    again = m.process_repeat_hits(cat, root / "pr", library=work["lib"], verbose=False)
    assert again["annotations"] == 4
    gff = m.convert_repeats_to_gff3(out, root / "gff", verbose=False)
    assert gff["features"] == 4
    bed = m.convert_repeats_to_bed(out, root / "bed", verbose=False)
    assert bed["intervals"] == 4
    assert m.export_repeats_ucsc_tables(out, root / "ucsc", verbose=False)["tables"]["genome.fa.out.tsv"] == 4
    assert m.extract_repeat_sequences(work["fa"], out, root / "ex", verbose=False)["sequences"] == 4
    masked = m.mask_sequences_from_annotation(work["fa"], out, root / "mf", soft_mask=True, verbose=False)
    assert masked["masked_bases"] >= 1100
    div = m.recalculate_repeat_divergence(aln, root / "div", verbose=False)
    assert div["families"] and div["families"][0]["repeat"] == "Rep1"
    land = m.plot_repeat_landscape(div["outputs"]["repeats.divsum"], root / "land", genome_size=5200, verbose=False)
    assert land["classes"] == ["SINE/Alu"] and land["landscape"]
    comb = m.combine_repeat_annotations(out, aln, out, aln, root / "comb", verbose=False)
    assert comb["annotations"] >= 4
    assert m.mask_tandem_repeats(work["fa"], root / "trf", verbose=False)["summary"]["sequences"] == 1


@requires_rm
def test_summary_needs_two_bit_info(work: dict[str, Path], tmp_path: Path) -> None:
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "twoBitInfo").write_text("#!/bin/sh\nexit 0\n")
    (stub / "twoBitInfo").chmod(0o755)
    sizes = tmp_path / "g.tsv"
    sizes.write_text("seq1\t2132\n")
    import os

    old = os.environ["PATH"]
    os.environ["PATH"] = f"{stub}:{old}"
    try:
        report = m.summarize_repeat_annotation(work["out"], tmp_path / "s", genome=sizes, verbose=False)
    finally:
        os.environ["PATH"] = old
    assert "Alu" in report["report"]


@pytest.mark.skipif(sys.platform == "win32", reason="stand-in executables are POSIX scripts")
def test_standin_programs_receive_expected_arguments(tmp_path: Path, work: dict[str, Path]) -> None:
    home = tmp_path / "home"
    (home / "util").mkdir(parents=True)
    body = f"#!{sys.executable}\nimport sys, pathlib\na = sys.argv[1:]\n"
    body += "if a == ['-v']:\n    print('RepeatMasker version 4.2.6'); sys.exit(0)\n"
    body += "pathlib.Path('args.txt').write_text(' '.join(a))\n"
    for name in ("RepeatMasker", "DateRepeats", "DupMasker", "RepeatProteinMask", "ProcessRepeats"):
        (home / name).write_text(body + "n = a[-1]\nfor suffix in ('.withoutIS', '.duplicons', '.annot'):\n    pathlib.Path(n + suffix).write_text('x')\n")
        (home / name).chmod(0o755)
    clip = m.clip_bacterial_insertion_sequences(work["seq"], tmp_path / "c", repeatmasker_dir=home, verbose=False)
    assert "-is_only" in (Path(clip["output_directory"]) / "args.txt").read_text().split() and clip["version"] == "4.2.6"
    dup = m.annotate_segmental_duplications(work["seq"], tmp_path / "d", max_divergence=10.0, gff_output=True, repeatmasker_dir=home, verbose=False)
    assert dup["duplicons"] and "-maxDiv" in (Path(dup["output_directory"]) / "args.txt").read_text()
    prot = m.mask_repeat_proteins(work["seq"], tmp_path / "p", p_value=0.01, repeatmasker_dir=home, verbose=False)
    assert "-pvalue" in (Path(prot["output_directory"]) / "args.txt").read_text()
    dated = m.date_repeats(
        work["out"], tmp_path / "dt", query_species="human", comparison_species=["mouse"], sequence_file=work["seq"],
        mask_species="mouse", repeatmasker_dir=home, verbose=False,
    )
    text = (Path(dated["output_directory"]) / "args.txt").read_text()
    assert "-query human" in text and "-comp mouse" in text and "-mask mouse" in text
    with pytest.raises(ConfigurationError):
        m.date_repeats(work["out"], tmp_path / "dt2", query_species="human", comparison_species=["mouse"], mask_species="mouse", repeatmasker_dir=home, verbose=False)
