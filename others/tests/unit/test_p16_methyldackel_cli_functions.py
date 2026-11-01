"""MethylDackel facades on a synthetic directional bisulfite BAM (skip without MethylDackel)."""

from __future__ import annotations

from pathlib import Path

import pytest
from _p16_data import bisulfite_data, need

import dnakit.epigenome.p16_methyldackel_cli_functions as md
from dnakit.exceptions import ConfigurationError


def test_report_merge_bias_per_read_and_bbm(tmp_path: Path) -> None:
    program = need("MethylDackel")
    data = bisulfite_data(tmp_path)
    ref, bam = data["ref.fa"], data["bs.bam"]
    cytosine = md.extract_methylation_report_methyldackel(ref, bam, str(tmp_path / "cyt"), methyldackel_path=program, verbose=False)
    rows = [line.split("\t") for line in (tmp_path / "cyt.cytosine_report.txt").read_text().splitlines()]
    assert all(len(r) == 7 for r in rows) and {r[2] for r in rows} == {"+", "-"} and cytosine["report"] == "cytosine_report"
    assert any(r[5] == "CG" and int(r[3]) + int(r[4]) > 0 for r in rows)
    merged_direct = md.extract_methylation_report_methyldackel(ref, bam, str(tmp_path / "mc"), report="merge_context", methyldackel_path=program, verbose=False)
    plain = md.extract_methylation_report_methyldackel(ref, bam, str(tmp_path / "pl"), report="per_cytosine", methyldackel_path=program, verbose=False)
    two_step = md.merge_cytosine_context_methyldackel(ref, str(tmp_path / "pl_CpG.bedGraph"), str(tmp_path / "m2.bedGraph"), methyldackel_path=program, verbose=False)
    body = lambda p: [l for l in Path(p).read_text().splitlines() if not l.startswith("track")]  # noqa: E731
    assert body(tmp_path / "mc_CpG.bedGraph") == body(tmp_path / "m2.bedGraph") and two_step["sites"] > 0 and merged_direct["files"]
    bias = md.methylation_bias_methyldackel(ref, bam, str(tmp_path / "mb"), methyldackel_path=program, verbose=False)
    assert bias["suggested_bounds"].get("OT") and all(p.endswith(".svg") for p in bias["outputs"])
    per_read = md.summarize_per_read_methylation_methyldackel(ref, bam, str(tmp_path / "pr.tsv"), methyldackel_path=program, verbose=False)
    assert per_read["rows"][0]["read"] == "r0" and 0 <= per_read["rows"][0]["methylation_percent"] <= 100
    bbm = md.convert_mappability_bigwig_to_bbm_methyldackel(ref, bam, data["map.bw"], str(tmp_path / "map.bbm"), str(tmp_path / "bb"), methyldackel_path=program, verbose=False)
    assert Path(bbm["bbm"]).read_bytes()[:3] == b"BBM" or bbm["bbm_bytes"] > 0
    with pytest.raises(ConfigurationError):
        md.extract_methylation_report_methyldackel(ref, bam, str(tmp_path / "x"), report="bogus", verbose=False)
    with pytest.raises(ConfigurationError):
        md.extract_methylation_report_methyldackel(ref, bam, str(tmp_path / "y"), inclusion_bounds_ot="1,2", verbose=False)
