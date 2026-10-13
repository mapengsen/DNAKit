"""Contract tests for the G18m review wrappers: methylpy, MACS3 and deepTools (other rows)."""

from __future__ import annotations

import gzip
import importlib
import inspect
import json
import math
import os
import pkgutil
import random
import shutil
import struct
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

import dnakit.epigenome as macs3_epigenome_package
import dnakit.epigenome as methylpy_epigenome_package
from dnakit.epigenome import (
    MACS3_review_functions as macs3_review,
)
from dnakit.epigenome import (
    deeptools_review_functions as deeptools_review,
)
from dnakit.epigenome import (
    methylpy_review_functions as methylpy_review,
)
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _methylpy_write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _methylpy_residuals(mc: list[int], total: list[int]) -> list[list[float]]:
    """Residuals of the RMS test, as computed in methylpy's rms.cpp, per sample [mc, uc]."""

    grand = float(sum(total))
    rows = [[m, h - m] for m, h in zip(mc, total, strict=True)]
    column = [sum(row[j] for row in rows) for j in range(2)]
    result = []
    for row in rows:
        row_total = sum(row)
        line = []
        for j in range(2):
            expected = row_total * column[j] / grand
            denominator = (expected * (1.0 - row_total / grand) * (1 - column[j] / grand)) ** 0.5
            line.append(0.0 if denominator == 0 else (row[j] - expected) / denominator)
        result.append(line)
    return result


def _methylpy_rms_table(
    path: Path, rows: list[tuple[Any, ...]], samples: list[str], residuals: dict[int, list[float]]
) -> Path:
    """Write a gzip RMS results table.

    Each row is (chrom, pos, context, pvalue, mc list, h list, sims_sig, sims_run). ``residuals``
    maps a row index to an explicit [mc_res per sample..., uc_res per sample...] list.
    """

    header = ["chr", "pos", "strand", "mc_class", "pvalue"]
    for prefix in ("mc_", "h_", "frac_", "mc_residual_", "uc_residual_"):
        header += [prefix + name for name in samples]
    header += ["num_simulations_sig", "num_simulations_run"]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for index, (chrom, pos, context, pvalue, mc, h, sig, run) in enumerate(rows):
            if index in residuals:
                values = residuals[index]
            else:
                pairs = _methylpy_residuals(mc, h)
                values = [pair[0] for pair in pairs] + [pair[1] for pair in pairs]
            fractions = [m / n for m, n in zip(mc, h, strict=True)]
            fields = [chrom, str(pos), "+", context, pvalue, *map(str, mc), *map(str, h)]
            fields += [f"{value:.6f}" for value in fractions]
            fields += [f"{value:.6f}" for value in values]
            fields += [str(sig), str(run)]
            handle.write("\t".join(fields) + "\n")
    return path


# Four tests, three of them significant: used for a hand-calculated FDR cutoff.
_METHYLPY_SMALL_ROWS = [
    ("1", 100, "CGN", "0.01", [9, 1], [10, 10], 1, 100),
    ("1", 120, "CGN", "0.01", [10, 0], [10, 10], 1, 100),
    ("1", 140, "CGN", "0.02", [8, 1], [10, 10], 2, 100),
    ("1", 500, "CGN", "0.50", [5, 5], [10, 10], 5, 10),
]

# Ten tests: two DMR-like clusters on chromosome 1, one single site on chromosome 2.
_METHYLPY_MERGE_ROWS = [
    ("1", 100, "CGN", "0.01", [9, 1], [10, 10], 1, 100),
    ("1", 120, "CGN", "0.01", [10, 0], [10, 10], 1, 100),
    ("1", 140, "CGN", "0.02", [8, 1], [10, 10], 2, 100),
    ("1", 300, "CGN", "0.50", [5, 6], [10, 10], 5, 10),
    ("1", 320, "CGN", "0.42", [4, 5], [10, 10], 5, 12),
    ("1", 340, "CGN", "0.62", [3, 4], [10, 10], 5, 8),
    ("1", 600, "CGN", "0.01", [1, 9], [10, 10], 1, 100),
    ("1", 620, "CGN", "0.03", [0, 10], [10, 10], 3, 100),
    ("2", 50, "CGN", "0.50", [5, 5], [10, 10], 5, 10),
    ("2", 80, "CGN", "0.01", [9, 2], [10, 10], 1, 100),
]

_METHYLPY_ALLC_A = (
    "chr1\t100\t+\tCGA\t9\t10\t1\n"
    "chr1\t120\t+\tCGT\t10\t10\t1\n"
    "chr1\t125\t+\tCHH\t1\t4\t0\n"
    "chr1\t140\t+\tCGC\t8\t10\t1\n"
    "chr1\t500\t+\tCGG\t5\t10\t1\n"
    "chr2\t30\t+\tCGA\t0\t9\t0\n"
)
_METHYLPY_ALLC_B = (
    "chr1\t100\t+\tCGA\t1\t10\t0\n"
    "chr1\t120\t+\tCGT\t0\t10\t0\n"
    "chr1\t140\t+\tCGC\t1\t10\t0\n"
    "chr1\t500\t+\tCGG\t5\t10\t1\n"
    "chr2\t30\t+\tCGA\t3\t6\t1\n"
)
_METHYLPY_DMRS = (
    "#chr\tstart\tend\tnumber_of_dms\thypermethylated_samples\thypomethylated_samples\n"
    "chr1\t100\t140\t3\tA\tB\n"
    "chr1\t120\t500\t2\tA\tB\n"
    "chr1\t600\t700\t1\tA\tB\n"
    "chr2\t1\t40\t1\tA\tB\n"
)


def _methylpy_allc_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    return (
        _methylpy_write(tmp_path / "allc_A.tsv", _METHYLPY_ALLC_A),
        _methylpy_write(tmp_path / "allc_B.tsv", _METHYLPY_ALLC_B),
        _methylpy_write(tmp_path / "dmrs.tsv", _METHYLPY_DMRS),
    )


def _methylpy_read_gz(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return handle.read()


def _methylpy_tree(directory: Path) -> list[str]:
    return sorted(item.name for item in directory.iterdir())


# --- module contract ----------------------------------------------------------------------


def test_methylpy_review_public_functions_take_keyword_verbose() -> None:
    assert len(methylpy_review.__all__) == 12
    for name in methylpy_review.__all__:
        parameter = inspect.signature(getattr(methylpy_review, name)).parameters.get("verbose")
        assert parameter is not None, name
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_methylpy_review_names_do_not_collide_with_other_epigenome_modules() -> None:
    ours = set(methylpy_review.__all__)
    for info in pkgutil.iter_modules(methylpy_epigenome_package.__path__):
        if not info.name.endswith("_functions") or info.name == "methylpy_review_functions":
            continue
        other = importlib.import_module(f"dnakit.epigenome.{info.name}")
        assert not ours & set(getattr(other, "__all__", ())), info.name


def test_methylpy_missing_package_reports_install_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mpileup = _methylpy_write(tmp_path / "x.mpileup", "chr1\t1\tC\t1\t.\n" * 4)
    for name in ("methylpy", "methylpy.utilities"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError, match="bio-methylpy"):
        methylpy_review.split_mpileup_chunks(mpileup, tmp_path / "chunk_", 2, verbose=False)


# --- estimate_dmr_fdr_cutoff --------------------------------------------------------------


def test_methylpy_estimate_dmr_fdr_cutoff_hand_calculation(tmp_path: Path) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    table = _methylpy_rms_table(
        tmp_path / "rms.tsv.gz", _METHYLPY_SMALL_ROWS, ["A", "B"], residuals={}
    )
    result = methylpy_review.estimate_dmr_fdr_cutoff(
        table, num_sims=100, num_sig_tests=5, target_fdr=0.05, convergence_diff=1, verbose=False
    )
    # Hand calculation with N = 4 tests. Iteration 0 starts with m0 = 4. Observed counts are 2 at
    # (1, 100) and 1 at (2, 100) against expected 0.04 and 0.08, so 1.96 + 0.92 are excess and
    # m0 = 4 - 2.88 = 1.12. With m0 = 1.12 the expected counts are 0.0112 per step, the excess is
    # 1.9888 + 0.9888 and m0 = 4 - 2.9776 = 1.0224; the change 0.0976 is below 1, so it stops.
    assert result["total_tests"] == 4
    assert result["null_tests_estimate"] == pytest.approx(1.0224)
    # Three tests are significant at (5, d) for d near 37: the first denominator whose rounded
    # P-value is 0.14 is 37, so the cutoff is 5/37 and its FDR is 0.14 * 1.0224 / 3.
    assert result["p_value_cutoff"] == pytest.approx(5 / 37)
    assert result["estimated_fdr"] == pytest.approx(0.14 * 1.0224 / 3)
    assert result["outcome"] == "converged"
    assert result["p_value_cutoff"] == dmrfind.histogram_correction_DMRfind(
        str(table), 100, 5, target_fdr=0.05, convergence_diff=1
    )


def test_methylpy_estimate_dmr_fdr_cutoff_fallback_is_smallest_pvalue(tmp_path: Path) -> None:
    rows = [("1", 100 + 10 * i, "CGN", "0.50", [5, 5], [10, 10], 5, 10) for i in range(4)]
    table = _methylpy_rms_table(tmp_path / "null.tsv.gz", rows, ["A", "B"], residuals={})
    result = methylpy_review.estimate_dmr_fdr_cutoff(
        table, num_sims=100, num_sig_tests=5, convergence_diff=1, verbose=False
    )
    # Nothing is below the expected null count, so m0 stays at 4 after one pass and methylpy
    # falls back to 1 / num_sims.
    assert result["outcome"] == "not_converged"
    assert result["p_value_cutoff"] == pytest.approx(0.01)


def test_methylpy_estimate_dmr_fdr_cutoff_rejects_bad_tables(tmp_path: Path) -> None:
    good = _methylpy_rms_table(tmp_path / "ok.tsv.gz", _METHYLPY_SMALL_ROWS, ["A", "B"], {})
    with pytest.raises(ConfigurationError, match="impossible"):
        methylpy_review.estimate_dmr_fdr_cutoff(good, num_sims=100, num_sig_tests=4, verbose=False)
    with pytest.raises(ConfigurationError, match="num_sig_tests"):
        methylpy_review.estimate_dmr_fdr_cutoff(good, num_sims=10, num_sig_tests=11, verbose=False)
    with pytest.raises(ConfigurationError, match="target_fdr"):
        methylpy_review.estimate_dmr_fdr_cutoff(
            good, num_sims=100, num_sig_tests=5, target_fdr=0.0, verbose=False
        )
    plain = _methylpy_write(tmp_path / "plain.tsv", "chr\tpos\n")
    with pytest.raises(ConfigurationError, match="gzip"):
        methylpy_review.estimate_dmr_fdr_cutoff(plain, num_sims=100, num_sig_tests=5, verbose=False)
    bad_header = tmp_path / "header.tsv.gz"
    with gzip.open(bad_header, "wt", encoding="utf-8") as handle:
        handle.write("chr\tpos\tvalue\n")
    with pytest.raises(ConfigurationError, match="header"):
        methylpy_review.estimate_dmr_fdr_cutoff(
            bad_header, num_sims=100, num_sig_tests=5, verbose=False
        )
    empty = _methylpy_rms_table(tmp_path / "empty.tsv.gz", [], ["A", "B"], {})
    with pytest.raises(ConfigurationError, match="no tests"):
        methylpy_review.estimate_dmr_fdr_cutoff(empty, num_sims=100, num_sig_tests=5, verbose=False)


# --- dmr_residual_cutoff ------------------------------------------------------------------


def test_methylpy_dmr_residual_cutoff_hand_percentile(tmp_path: Path) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    rows = [
        ("1", 100, "CGN", "0.01", [9, 1], [10, 10], 1, 100),
        ("1", 200, "CGN", "0.90", [5, 5], [10, 10], 5, 10),
        ("1", 300, "CGN", "0.80", [5, 5], [10, 10], 5, 12),
    ]
    explicit = {1: [1.0, -1.0, -1.0, 1.0], 2: [2.0, -2.0, -2.0, 2.0]}
    table = _methylpy_rms_table(tmp_path / "rms.tsv.gz", rows, ["A", "B"], explicit)
    result = methylpy_review.dmr_residual_cutoff(table, 0.5, resid_cutoff=0.25, verbose=False)
    # Two sites are above P = 0.5; their eight residuals sorted are -2 -2 -1 -1 1 1 2 2. The 75th
    # percentile sits at index 0.75 * 7 = 5.25: 1 + 0.25 * (2 - 1) = 1.25.
    assert result["null_sites"] == 2
    assert result["default_used"] is False
    assert result["residual_cutoff"] == pytest.approx(1.25)
    assert result["residual_cutoff"] == pytest.approx(
        float(dmrfind.get_resid_cutoff(0.25, 0.5, str(table)))
    )
    # No site is above P = 0.95, so methylpy's documented default of 2 is used.
    fallback = methylpy_review.dmr_residual_cutoff(table, 0.95, verbose=False)
    assert fallback["default_used"] is True
    assert fallback["residual_cutoff"] == 2.0


def test_methylpy_dmr_residual_cutoff_from_counts(tmp_path: Path) -> None:
    table = _methylpy_rms_table(tmp_path / "rms.tsv.gz", _METHYLPY_MERGE_ROWS, ["A", "B"], {})
    result = methylpy_review.dmr_residual_cutoff(table, 0.135135, resid_cutoff=0.2, verbose=False)
    # Sites above P = 0.135135: chr1:300, 320, 340 and chr2:50. For mc counts 5 and 6 of 10 reads
    # the residual is 0.5 / sqrt(5.5 * 0.5 * 0.45) = 0.4495 in absolute value, for 4 and 5 it is the
    # same, for 3 and 4 it is 0.5 / sqrt(3.5 * 0.5 * 0.65) = 0.4688, and chr2:50 has residual 0.
    # Of the 16 residuals 6 are negative, 4 are zero and 6 are positive; the 80th percentile
    # (index 0.8 * 15 = 12 of the sorted values) is the third largest, 0.4495.
    assert result["null_sites"] == 4
    assert result["residual_cutoff"] == pytest.approx(0.5 / (5.5 * 0.5 * 0.45) ** 0.5, abs=1e-5)


def test_methylpy_dmr_residual_cutoff_rejects_less_than_pvalues(tmp_path: Path) -> None:
    rows = [("1", 100, "CGN", "<1.0e-02", [9, 1], [10, 10], 1, 100)]
    table = _methylpy_rms_table(tmp_path / "lt.tsv.gz", rows, ["A", "B"], {})
    with pytest.raises(ConfigurationError, match="starting with '<'"):
        methylpy_review.dmr_residual_cutoff(table, 0.5, verbose=False)
    with pytest.raises(ConfigurationError, match="pvalue_cutoff"):
        methylpy_review.dmr_residual_cutoff(table, 1.5, verbose=False)


# --- collapse_dmr_windows -----------------------------------------------------------------


def _methylpy_collapse_table(tmp_path: Path) -> Path:
    rows = [
        ("1", 100, "CGN", "0.001", [9, 1], [10, 10], 1, 1000),
        ("1", 120, "CGN", "0.001", [10, 0], [10, 10], 1, 1000),
        ("1", 140, "CGN", "0.001", [8, 1], [10, 10], 1, 1000),
        ("1", 400, "CGN", "0.001", [1, 9], [10, 10], 1, 1000),
        ("1", 410, "CGN", "0.001", [0, 10], [10, 10], 1, 1000),
        ("1", 700, "CGN", "0.5", [5, 5], [10, 10], 5, 10),
    ]
    return _methylpy_rms_table(tmp_path / "rms.tsv.gz", rows, ["A", "B"], {})


def test_methylpy_collapse_dmr_windows_hand_regions(tmp_path: Path) -> None:
    table = _methylpy_collapse_table(tmp_path)
    result = methylpy_review.collapse_dmr_windows(
        table, tmp_path / "dmrs.tsv", max_dist=50, sig_cutoff=0.01, resid_cutoff=2.0, verbose=False
    )
    # Sites 100, 120 and 140 are 20 apart and sample A is above expectation in all of them, so
    # they form one region; 400 is 260 away and starts a new one that 410 joins. Sample A is below
    # expectation there, so the hyper- and hypo-methylated samples swap.
    assert result["samples"] == ["A", "B"]
    assert result["dmr_count"] == 2
    assert result["dms_count"] == 5
    assert result["regions"] == [
        {
            "chrom": "1",
            "start": 100,
            "end": 140,
            "dms": 3,
            "hypermethylated_samples": ["A"],
            "hypomethylated_samples": ["B"],
        },
        {
            "chrom": "1",
            "start": 400,
            "end": 410,
            "dms": 2,
            "hypermethylated_samples": ["B"],
            "hypomethylated_samples": ["A"],
        },
    ]
    dms = Path(result["files"]["dms_bed"]).read_text(encoding="utf-8").splitlines()
    assert dms[0] == "chr1\t99\t100\tCGN;+\t0.001"
    assert len(dms) == 5
    dmr = Path(result["files"]["dmr_bed"]).read_text(encoding="utf-8").splitlines()
    assert dmr == ["chr1\t99\t140\t3", "chr1\t399\t410\t2"]


def test_methylpy_collapse_dmr_windows_matches_upstream_and_filters(tmp_path: Path) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    table = _methylpy_collapse_table(tmp_path)
    result = methylpy_review.collapse_dmr_windows(
        table,
        tmp_path / "wrapped.tsv",
        max_dist=50,
        sig_cutoff=0.01,
        resid_cutoff=2.0,
        min_num_dms=3,
        verbose=False,
    )
    direct = tmp_path / "direct.tsv"
    count = dmrfind.collapse_dmr_windows(
        str(table), str(direct), max_dist=50, sig_cutoff=0.01, resid_cutoff=2.0, min_num_dms=3
    )
    assert count == result["dmr_count"] == 1
    assert (tmp_path / "wrapped.tsv").read_text() == direct.read_text()
    assert Path(result["files"]["dmr_bed"]).read_text() == (
        (tmp_path / "direct.tsv.DMR.bed").read_text()
    )


def test_methylpy_collapse_dmr_windows_without_significant_sites(tmp_path: Path) -> None:
    rows = [("1", 100, "CGN", "0.5", [5, 5], [10, 10], 5, 10)]
    table = _methylpy_rms_table(tmp_path / "none.tsv.gz", rows, ["A", "B"], {})
    result = methylpy_review.collapse_dmr_windows(table, tmp_path / "out.tsv", verbose=False)
    assert result["dmr_count"] == 0
    assert result["regions"] == []
    assert result["files"]["dms_bed"] is None
    assert result["files"]["dmr_bed"] is None
    assert (tmp_path / "out.tsv").read_text().count("\n") == 1


def test_methylpy_collapse_dmr_windows_validation(tmp_path: Path) -> None:
    table = _methylpy_collapse_table(tmp_path)
    taken = _methylpy_write(tmp_path / "taken.tsv", "x")
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.collapse_dmr_windows(table, taken, verbose=False)
    _methylpy_write(tmp_path / "bed.tsv.DMS.bed", "x")
    with pytest.raises(ConfigurationError, match=r"DMS\.bed"):
        methylpy_review.collapse_dmr_windows(table, tmp_path / "bed.tsv", verbose=False)
    with pytest.raises(ConfigurationError, match="collapse_samples"):
        methylpy_review.collapse_dmr_windows(
            table, tmp_path / "a.tsv", collapse_samples=["A", "Z"], verbose=False
        )
    with pytest.raises(ConfigurationError, match="needs collapse_samples"):
        methylpy_review.collapse_dmr_windows(
            table, tmp_path / "b.tsv", sample_category=[0, 1], verbose=False
        )
    with pytest.raises(ConfigurationError, match="one entry per"):
        methylpy_review.collapse_dmr_windows(
            table,
            tmp_path / "c.tsv",
            collapse_samples=["A", "B"],
            sample_category=[0],
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="sig_cutoff"):
        methylpy_review.collapse_dmr_windows(table, tmp_path / "d.tsv", sig_cutoff=2, verbose=False)
    assert not (tmp_path / "a.tsv").exists()


def test_methylpy_collapse_dmr_windows_sample_categories(tmp_path: Path) -> None:
    table = _methylpy_collapse_table(tmp_path)
    kept = methylpy_review.collapse_dmr_windows(
        table,
        tmp_path / "kept.tsv",
        max_dist=50,
        sig_cutoff=0.01,
        resid_cutoff=2.0,
        collapse_samples=["A", "B"],
        sample_category=[0, 1],
        min_cluster=1,
        verbose=False,
    )
    # A and B sit in different categories and move in opposite directions, so both regions stay.
    assert kept["dmr_count"] == 2


# --- merge_dms_to_dmrs --------------------------------------------------------------------


def test_methylpy_merge_dms_to_dmrs_matches_upstream_and_hand_regions(tmp_path: Path) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    table = _methylpy_rms_table(tmp_path / "rms.tsv.gz", _METHYLPY_MERGE_ROWS, ["A", "B"], {})
    result = methylpy_review.merge_dms_to_dmrs(
        table,
        tmp_path / "wrapped.tsv",
        num_sims=100,
        num_sig_tests=5,
        sig_cutoff=0.1,
        dmr_max_dist=50,
        resid_cutoff=0.2,
        verbose=False,
    )
    direct = tmp_path / "direct.tsv"
    dmrfind.merge_DMS_to_DMR(
        str(table),
        str(direct),
        collapse_samples=False,
        sample_category=False,
        min_cluster=0,
        sig_cutoff=0.1,
        dmr_max_dist=50,
        min_num_dms=0,
        resid_cutoff=0.2,
        num_sims=100,
        num_sig_tests=5,
        max_iterations=1000,
        convergence_diff=1,
    )
    assert (tmp_path / "wrapped.tsv").read_text() == direct.read_text()
    assert (tmp_path / "wrapped.tsv.DMS.bed").read_text() == (
        (tmp_path / "direct.tsv.DMS.bed").read_text()
    )
    assert result["fdr_outcome"] == "converged"
    assert result["p_value_cutoff"] == pytest.approx(5 / 37)
    # Significant sites (P <= 5/37): 100, 120, 140 with A above expectation, 600 and 620 with A
    # below it, and 80 on chromosome 2. The five other tests are above the cutoff.
    regions = [
        (r["chrom"], r["start"], r["end"], r["dms"], r["hypermethylated_samples"])
        for r in result["regions"]
    ]
    assert regions == [
        ("1", 100, 140, 3, ["A"]),
        ("1", 600, 620, 2, ["B"]),
        ("2", 80, 80, 1, ["A"]),
    ]
    assert result["dms_count"] == 6
    assert result["residual_cutoff"] == pytest.approx(
        float(dmrfind.get_resid_cutoff(0.2, result["p_value_cutoff"], str(table)))
    )


def test_methylpy_merge_dms_to_dmrs_rejects_wrong_simulation_settings(tmp_path: Path) -> None:
    table = _methylpy_rms_table(tmp_path / "rms.tsv.gz", _METHYLPY_MERGE_ROWS, ["A", "B"], {})
    with pytest.raises(ConfigurationError, match="impossible"):
        methylpy_review.merge_dms_to_dmrs(table, tmp_path / "out.tsv", verbose=False)
    assert _methylpy_tree(tmp_path) == ["rms.tsv.gz"]


# --- add_methylation_levels_to_dmrs and add_cytosine_info_to_dmrs ------------------------------


def test_methylpy_add_methylation_levels_hand_values(tmp_path: Path) -> None:
    a, b, dmrs = _methylpy_allc_inputs(tmp_path)
    result = methylpy_review.add_methylation_levels_to_dmrs(
        dmrs, [a, b], ["A", "B"], tmp_path / "levels.tsv", verbose=False
    )
    levels = [(r["chrom"], r["start"], r["methylation_level"]) for r in result["rows"]]
    # CGN sites of region chr1:100-140 in A: 9 + 10 + 8 methylated of 30 reads (the CHH site at 125
    # is skipped) = 0.9; in B: 1 + 0 + 1 of 30. Region chr1:120-500 overlaps the first one:
    # A (10 + 8 + 5) / 30, B (0 + 1 + 5) / 30. Region chr1:600-700 holds no site. chr2:1-40: A 0/9,
    # B 3/6.
    assert levels[0] == ("chr1", 100, {"A": pytest.approx(0.9), "B": pytest.approx(2 / 30)})
    assert levels[1] == ("chr1", 120, {"A": pytest.approx(23 / 30), "B": pytest.approx(0.2)})
    assert levels[2] == ("chr1", 600, {"A": None, "B": None})
    assert levels[3] == ("chr2", 1, {"A": 0.0, "B": 0.5})
    text = (tmp_path / "levels.tsv").read_text().splitlines()
    assert text[0].endswith("methylation_level_A\tmethylation_level_B")
    assert result["dmr_count"] == 4
    assert not [name for name in _methylpy_tree(tmp_path) if "temp" in name]


def test_methylpy_add_methylation_levels_matches_upstream_and_context_filter(
    tmp_path: Path,
) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    a, b, dmrs = _methylpy_allc_inputs(tmp_path)
    methylpy_review.add_methylation_levels_to_dmrs(
        dmrs, [a, b], ["A", "B"], tmp_path / "wrapped.tsv", min_cov=5, verbose=False
    )
    dmrfind.get_methylation_levels_DMRfind(
        str(dmrs), str(tmp_path / "direct.tsv"), [str(a), str(b)], samples=["A", "B"], min_cov=5
    )
    assert (tmp_path / "wrapped.tsv").read_text() == (tmp_path / "direct.tsv").read_text()
    chh = methylpy_review.add_methylation_levels_to_dmrs(
        dmrs, [a, b], ["A", "B"], tmp_path / "chh.tsv", mc_type="CHH", verbose=False
    )
    # Only the CHH site at chr1:125 (1 of 4 reads in A) lies in region chr1:100-140; B has none.
    assert chh["rows"][0]["methylation_level"] == {"A": 0.25, "B": None}


def test_methylpy_add_methylation_levels_without_header(tmp_path: Path) -> None:
    a, b, _ = _methylpy_allc_inputs(tmp_path)
    plain = _methylpy_write(tmp_path / "plain.tsv", "chr1\t100\t140\nchr2\t1\t40\n")
    result = methylpy_review.add_methylation_levels_to_dmrs(
        plain, [a, b], ["A", "B"], tmp_path / "out.tsv", input_no_header=True, verbose=False
    )
    assert [r["methylation_level"]["A"] for r in result["rows"]] == [pytest.approx(0.9), 0.0]
    assert (tmp_path / "out.tsv").read_text().splitlines()[0].startswith("chr\tstart\tend")


def test_methylpy_add_cytosine_info_hand_values(tmp_path: Path) -> None:
    a, b, dmrs = _methylpy_allc_inputs(tmp_path)
    result = methylpy_review.add_cytosine_info_to_dmrs(
        dmrs, [a, b], ["A", "B"], tmp_path / "info.tsv", verbose=False
    )
    first = result["rows"][0]["samples"]
    assert first["A"] == {
        "methylation_level": pytest.approx(0.9),
        "total_basecalls": 30,
        "covered_sites": 3,
    }
    assert first["B"]["total_basecalls"] == 30
    assert result["rows"][2]["samples"]["A"] == {
        "methylation_level": None,
        "total_basecalls": None,
        "covered_sites": None,
    }
    last = result["rows"][3]["samples"]
    assert (last["A"]["total_basecalls"], last["A"]["covered_sites"]) == (9, 1)
    assert (last["B"]["total_basecalls"], last["B"]["covered_sites"]) == (6, 1)
    assert Path(result["total_basecalls_output"]).name == "info.tsv.total_basecalls.tsv"
    assert Path(result["covered_sites_output"]).is_file()
    assert not [name for name in _methylpy_tree(tmp_path) if "temp" in name]


def test_methylpy_add_cytosine_info_matches_upstream(tmp_path: Path) -> None:
    dmrfind = pytest.importorskip("methylpy.DMRfind")
    a, b, dmrs = _methylpy_allc_inputs(tmp_path)
    methylpy_review.add_cytosine_info_to_dmrs(
        dmrs, [a, b], ["A", "B"], tmp_path / "wrapped.tsv", max_cov=9, verbose=False
    )
    dmrfind.get_c_info_DMRfind(
        str(dmrs), [str(a), str(b)], ["A", "B"], str(tmp_path / "direct.tsv"), max_cov=9
    )
    for suffix in ("", ".total_basecalls.tsv", ".covered_sites.tsv"):
        wrapped = (tmp_path / f"wrapped.tsv{suffix}").read_text()
        direct = (tmp_path / f"direct.tsv{suffix}").read_text()
        # Only the header differs between the two files if the sample columns agree.
        assert wrapped.splitlines()[1:] == direct.splitlines()[1:]


def test_methylpy_dmr_joins_reject_bad_inputs_and_leave_no_files(tmp_path: Path) -> None:
    a, b, dmrs = _methylpy_allc_inputs(tmp_path)
    before = _methylpy_tree(tmp_path)
    missing = _methylpy_write(tmp_path / "missing.tsv", "#c\ts\te\nchr1\t100\t140\nchr3\t1\t40\n")
    with pytest.raises(ConfigurationError, match="chr3"):
        methylpy_review.add_methylation_levels_to_dmrs(
            missing, [a, b], ["A", "B"], tmp_path / "x.tsv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="chr3"):
        methylpy_review.add_cytosine_info_to_dmrs(
            missing, [a, b], ["A", "B"], tmp_path / "x.tsv", verbose=False
        )
    unsorted = _methylpy_write(
        tmp_path / "unsorted.tsv", "#c\ts\te\nchr1\t500\t600\nchr1\t100\t140\n"
    )
    with pytest.raises(ConfigurationError, match="sorted"):
        methylpy_review.add_methylation_levels_to_dmrs(
            unsorted, [a, b], ["A", "B"], tmp_path / "x.tsv", verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"\.tsv"):
        methylpy_review.add_methylation_levels_to_dmrs(
            dmrs, [a, b], ["A", "B"], tmp_path / "levels.txt", verbose=False
        )
    with pytest.raises(ConfigurationError, match="samples"):
        methylpy_review.add_methylation_levels_to_dmrs(
            dmrs, [a, b], ["A"], tmp_path / "x.tsv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="samples"):
        methylpy_review.add_methylation_levels_to_dmrs(
            dmrs, [a, b], ["A", "A"], tmp_path / "x.tsv", verbose=False
        )
    spaced = tmp_path / "with space"
    spaced.mkdir()
    with pytest.raises(ConfigurationError, match="shell"):
        methylpy_review.add_methylation_levels_to_dmrs(
            dmrs, [a, b], ["A", "B"], spaced / "x.tsv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="min_cov"):
        methylpy_review.add_cytosine_info_to_dmrs(
            dmrs, [a, b], ["A", "B"], tmp_path / "x.tsv", min_cov=-1, verbose=False
        )
    with pytest.raises(ConfigurationError, match="mc_type"):
        methylpy_review.add_cytosine_info_to_dmrs(
            dmrs, [a, b], ["A", "B"], tmp_path / "x.tsv", mc_type="XYZ!", verbose=False
        )
    created = set(_methylpy_tree(tmp_path)) - set(before)
    assert created <= {
        "missing.tsv",
        "unsorted.tsv",
        "with space",
        "allc_A.tsv.idx",
        "allc_B.tsv.idx",
    }


# --- split_allc_by_chromosome and split_allc_into_chunks ----------------------------------------


def test_methylpy_split_allc_by_chromosome_reports_real_files(tmp_path: Path) -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    a, _, _ = _methylpy_allc_inputs(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    result = methylpy_review.split_allc_by_chromosome(a, "A", out, verbose=False)
    assert [item["chrom"] for item in result["files"]] == ["chr1", "chr2"]
    assert [item["records"] for item in result["files"]] == [5, 1]
    assert _methylpy_tree(out) == ["allc_A_chr1.tsv.gz", "allc_A_chr2.tsv.gz"]
    assert _methylpy_read_gz(Path(result["files"][1]["path"])) == "chr2\t30\t+\tCGA\t0\t9\t0\n"
    first_lines = _methylpy_read_gz(Path(result["files"][0]["path"]))
    assert first_lines == _METHYLPY_ALLC_A.split("chr2")[0]
    # methylpy lists the names without .gz even though it wrote .tsv.gz files.
    direct_dir = tmp_path / "direct"
    direct_dir.mkdir()
    names = call_mc.do_split_allc_file(str(a), "A", str(direct_dir), compress_output=True)
    assert [Path(name).name for name in names] == ["allc_A_chr1.tsv", "allc_A_chr2.tsv"]
    assert _methylpy_tree(direct_dir) == _methylpy_tree(out)
    plain = tmp_path / "plain"
    plain.mkdir()
    flat = methylpy_review.split_allc_by_chromosome(
        a, "A", plain, compress_output=False, verbose=False
    )
    assert _methylpy_tree(plain) == ["allc_A_chr1.tsv", "allc_A_chr2.tsv"]
    assert flat["compressed"] is False


def test_methylpy_split_allc_by_chromosome_rejects_unsafe_tables(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    empty = _methylpy_write(tmp_path / "empty.tsv", "")
    with pytest.raises(ConfigurationError, match="empty"):
        methylpy_review.split_allc_by_chromosome(empty, "E", out, verbose=False)
    interleaved = _methylpy_write(
        tmp_path / "interleaved.tsv",
        "chr1\t1\t+\tCGA\t1\t2\t0\nchr2\t1\t+\tCGA\t1\t2\t0\nchr1\t5\t+\tCGA\t1\t2\t0\n",
    )
    with pytest.raises(ConfigurationError, match="more than one block"):
        methylpy_review.split_allc_by_chromosome(interleaved, "U", out, verbose=False)
    short = _methylpy_write(tmp_path / "short.tsv", "chr1\t1\t+\n")
    with pytest.raises(ConfigurationError, match="fewer than 6"):
        methylpy_review.split_allc_by_chromosome(short, "S", out, verbose=False)
    odd = _methylpy_write(tmp_path / "odd.tsv", "chr/1\t1\t+\tCGA\t1\t2\t0\n")
    with pytest.raises(ConfigurationError, match="file name"):
        methylpy_review.split_allc_by_chromosome(odd, "S", out, verbose=False)
    a, _, _ = _methylpy_allc_inputs(tmp_path)
    methylpy_review.split_allc_by_chromosome(a, "A", out, verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.split_allc_by_chromosome(a, "A", out, verbose=False)
    with pytest.raises(ConfigurationError, match="sample"):
        methylpy_review.split_allc_by_chromosome(a, "bad name", out, verbose=False)
    with pytest.raises(ConfigurationError, match="output_dir"):
        methylpy_review.split_allc_by_chromosome(a, "A", tmp_path / "nope", verbose=False)
    assert _methylpy_tree(out) == ["allc_A_chr1.tsv.gz", "allc_A_chr2.tsv.gz"]


def test_methylpy_split_allc_into_chunks_hand_packing(tmp_path: Path) -> None:
    lines = (
        "chr1\t1\t+\tCGA\t1\t2\t0\nchr1\t2\t+\tCGA\t1\t2\t0\nchr1\t3\t+\tCGA\t1\t2\t0\n"
        "chr2\t1\t+\tCGA\t1\t2\t0\n"
        "chr3\t1\t+\tCGA\t1\t2\t0\nchr3\t2\t+\tCGA\t1\t2\t0\n"
    )
    source = _methylpy_write(tmp_path / "allc.tsv", lines)
    out = tmp_path / "out"
    out.mkdir()
    result = methylpy_review.split_allc_into_chunks(source, "S", 2, out, verbose=False)
    # Six records in two chunks give a chunk size of 3. Chromosomes are packed smallest first:
    # chr2 (1) and chr3 (2) reach 3 records in chunk 0, then chr1 (3) opens chunk 1.
    assert [(c["index"], c["records"], c["chromosomes"]) for c in result["chunks"]] == [
        (0, 3, ["chr2", "chr3"]),
        (1, 3, ["chr1"]),
    ]
    assert _methylpy_tree(out) == ["allc_S_0.tsv.gz", "allc_S_1.tsv.gz"]
    first = _methylpy_read_gz(out / "allc_S_0.tsv.gz").splitlines()
    assert first[0].startswith("chr2") and first[-1].startswith("chr3")
    assert result["requested_chunks"] == 2


def test_methylpy_split_allc_into_chunks_fewer_files_than_requested(tmp_path: Path) -> None:
    a, _, _ = _methylpy_allc_inputs(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    result = methylpy_review.split_allc_into_chunks(
        a, "A", 8, out, compress_output=False, verbose=False
    )
    assert 1 <= len(result["chunks"]) <= 2
    assert sum(chunk["records"] for chunk in result["chunks"]) == 6
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.split_allc_into_chunks(a, "A", 8, out, compress_output=False, verbose=False)
    with pytest.raises(ConfigurationError, match="num_chunks"):
        methylpy_review.split_allc_into_chunks(a, "A", 0, out, verbose=False)
    empty = _methylpy_write(tmp_path / "empty.tsv", "")
    with pytest.raises(ConfigurationError, match="empty"):
        methylpy_review.split_allc_into_chunks(empty, "E", 2, out, verbose=False)


# --- run_allc_binomial_test ---------------------------------------------------------------


_METHYLPY_BINOMIAL_ALLC = (
    "chr1\t10\t+\tCGA\t8\t10\t1\n"
    "chr1\t20\t+\tCGT\t1\t10\t1\n"
    "chr1\t30\t+\tCHH\t0\t9\t0\n"
    "chr1\t40\t+\tCHH\t5\t9\t1\n"
    "chr1\t50\t+\tCHH\t1\t50\t0\n"
    "chr1\t60\t+\tCGC\t0\t0\t0\n"
    "chrL\t10\t+\tCHH\t1\t100\t0\n"
    "chrL\t20\t+\tCHH\t1\t100\t0\n"
)


def test_methylpy_run_allc_binomial_test_hand_calculation(tmp_path: Path) -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    source = _methylpy_write(tmp_path / "allc_S.tsv", _METHYLPY_BINOMIAL_ALLC)
    out = tmp_path / "out"
    out.mkdir()
    result = methylpy_review.run_allc_binomial_test(source, "S", out, "chrL", verbose=False)
    # The control chrL has 2 methylated of 200 reads: non-conversion rate 0.01. P(X >= 8 of 10) is
    # about 4.5e-15 and P(X >= 5 of 9) about 1.2e-8; with one test (CGA) or four tests (CHH) the
    # Benjamini-Hochberg values are far below 0.01. P(X >= 1 of 10) = 1 - 0.99**10 = 0.0956 is not
    # significant for the only CGT test, and the CHH sites with 1 of 50 or 1 of 100 reads have
    # P = 0.395 and 0.634. Sites without methylated reads or without any reads are never called.
    assert result["non_conversion_rate"] == pytest.approx(0.01)
    final = out / "allc_S.tsv.gz"
    assert result["output"] == str(final)
    lines = [line.split("\t") for line in _methylpy_read_gz(final).splitlines()]
    flags = {(row[0], int(row[1])): row[6] for row in lines}
    assert flags == {
        ("chr1", 10): "1",
        ("chr1", 20): "0",
        ("chr1", 30): "0",
        ("chr1", 40): "1",
        ("chr1", 50): "0",
        ("chrL", 10): "0",
        ("chrL", 20): "0",
    }
    assert result["sites"] == 7
    assert result["significant_sites"] == 2
    assert result["contexts"] == [
        {"context": "CGA", "sites": 1, "significant_sites": 1},
        {"context": "CGT", "sites": 1, "significant_sites": 0},
        {"context": "CHH", "sites": 5, "significant_sites": 1},
    ]
    assert Path(result["index"]).is_file()
    assert _methylpy_tree(out) == ["allc_S.tsv.gz"]
    # The same call with the numeric rate gives the same file as methylpy itself.
    direct = tmp_path / "direct"
    direct.mkdir()
    call_mc.perform_binomial_test(
        str(source), "S", str(direct), 0.01, min_cov=2, sig_cutoff=0.01, compress_output=False
    )
    plain = tmp_path / "plain"
    plain.mkdir()
    methylpy_review.run_allc_binomial_test(
        source, "S", plain, 0.01, compress_output=False, verbose=False
    )
    assert (plain / "allc_S.tsv").read_text() == (direct / "allc_S.tsv").read_text()


def test_methylpy_run_allc_binomial_test_validation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _methylpy_write(tmp_path / "allc_S.tsv", _METHYLPY_BINOMIAL_ALLC)
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(ConfigurationError, match="unmethylated_control"):
        methylpy_review.run_allc_binomial_test(source, "S", out, "1", verbose=False)
    with pytest.raises(ConfigurationError, match="unmethylated_control"):
        methylpy_review.run_allc_binomial_test(source, "S", out, 1.5, verbose=False)
    with pytest.raises(ConfigurationError, match="not in the output allc"):
        methylpy_review.run_allc_binomial_test(source, "S", out, "chrZ", verbose=False)
    with pytest.raises(ConfigurationError, match="sort_mem"):
        methylpy_review.run_allc_binomial_test(
            source, "S", out, 0.01, sort_mem="lots", verbose=False
        )
    with pytest.raises(ConfigurationError, match="min_cov"):
        methylpy_review.run_allc_binomial_test(source, "S", out, 0.01, min_cov=0, verbose=False)
    broken = _methylpy_write(tmp_path / "broken.tsv", "chr1\t10\t+\tCGA\tx\t10\t1\n")
    with pytest.raises(ConfigurationError, match="integer"):
        methylpy_review.run_allc_binomial_test(broken, "B", out, 0.01, verbose=False)
    methylpy_review.run_allc_binomial_test(source, "S", out, 0.01, verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.run_allc_binomial_test(source, "S", out, 0.01, verbose=False)
    assert _methylpy_tree(out) == ["allc_S.tsv.gz"]
    # sort runs with LC_ALL=C and the caller's setting is put back.
    other = tmp_path / "other"
    other.mkdir()
    monkeypatch.setenv("LC_ALL", "en_US.UTF-8")
    methylpy_review.run_allc_binomial_test(source, "S", other, 0.01, verbose=False)
    assert os.environ["LC_ALL"] == "en_US.UTF-8"
    third = tmp_path / "third"
    third.mkdir()
    monkeypatch.delenv("LC_ALL")
    methylpy_review.run_allc_binomial_test(source, "S", third, 0.01, verbose=False)
    assert "LC_ALL" not in os.environ


# --- filter_allc_by_pvalue_cutoffs ----------------------------------------------------------


def test_methylpy_filter_allc_by_pvalue_cutoffs_hand_values_and_inputs_untouched(
    tmp_path: Path,
) -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    first_text = (
        "chr1\t20\t+\tCGN\t1\t10\t0.2\n"
        "chr1\t10\t+\tCGN\t8\t10\t0.0001\n"
        "chr2\t5\t+\tCHH\t0\t9\t2.0\n"
    )
    second_text = "chr1\t40\t+\tCHH\t5\t9\t0.001\nchr1\t30\t+\tCHH\t0\t9\t0.6\n"
    first = _methylpy_write(tmp_path / "b1.tsv", first_text)
    second = _methylpy_write(tmp_path / "b2.tsv", second_text)
    result = methylpy_review.filter_allc_by_pvalue_cutoffs(
        [first, second], tmp_path / "filtered.tsv", {"CGN": 0.05, "CHH": 0.01}, verbose=False
    )
    # Each input is sorted by position, then P <= cutoff gives 1: chr1:10 (0.0001 <= 0.05) and
    # chr1:40 (0.001 <= 0.01) are significant; 0.2 > 0.05, 0.6 > 0.01 and the untested 2.0 are not.
    final = tmp_path / "filtered.tsv.gz"
    assert result["output"] == str(final)
    assert _methylpy_read_gz(final) == (
        "chr1\t10\t+\tCGN\t8\t10\t1\n"
        "chr1\t20\t+\tCGN\t1\t10\t0\n"
        "chr2\t5\t+\tCHH\t0\t9\t0\n"
        "chr1\t30\t+\tCHH\t0\t9\t0\n"
        "chr1\t40\t+\tCHH\t5\t9\t1\n"
    )
    assert result["records"] == 5
    assert result["significant_sites"] == 2
    # methylpy sorts its input in place; the wrapper keeps the originals as they were.
    assert first.read_text() == first_text
    assert second.read_text() == second_text
    assert not [name for name in _methylpy_tree(tmp_path) if name.startswith("methylpy_")]
    direct = tmp_path / "direct"
    direct.mkdir()
    copies = []
    for name, text in (("c1.tsv", first_text), ("c2.tsv", second_text)):
        copies.append(str(_methylpy_write(direct / name, text)))
    call_mc.filter_files_by_pvalue_combined(
        copies, str(direct / "out.tsv"), {"CGN": 0.05, "CHH": 0.01}, 1
    )
    assert _methylpy_read_gz(direct / "out.tsv.gz") == _methylpy_read_gz(final)
    assert Path(copies[0]).read_text() != first_text


def test_methylpy_filter_allc_by_pvalue_cutoffs_validation(tmp_path: Path) -> None:
    source = _methylpy_write(tmp_path / "b.tsv", "chr1\t10\t+\tCGN\t8\t10\t0.01\n")
    with pytest.raises(ConfigurationError, match="no cutoff"):
        methylpy_review.filter_allc_by_pvalue_cutoffs(
            [source], tmp_path / "out.tsv", {"CHH": 0.01}, verbose=False
        )
    with pytest.raises(ConfigurationError, match="cutoff for CGN"):
        methylpy_review.filter_allc_by_pvalue_cutoffs(
            [source], tmp_path / "out.tsv", {"CGN": 2.0}, verbose=False
        )
    short = _methylpy_write(tmp_path / "short.tsv", "chr1\t10\t+\tCGN\t8\t10\n")
    with pytest.raises(ConfigurationError, match="fewer than 7"):
        methylpy_review.filter_allc_by_pvalue_cutoffs(
            [short], tmp_path / "out.tsv", {"CGN": 0.01}, verbose=False
        )
    with pytest.raises(ConfigurationError, match="gz"):
        methylpy_review.filter_allc_by_pvalue_cutoffs(
            [source], tmp_path / "out.gz", {"CGN": 0.01}, compress_output=False, verbose=False
        )
    methylpy_review.filter_allc_by_pvalue_cutoffs(
        [source], tmp_path / "ok.tsv", {"CGN": 0.01}, compress_output=False, verbose=False
    )
    assert (tmp_path / "ok.tsv").read_text() == "chr1\t10\t+\tCGN\t8\t10\t1\n"
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.filter_allc_by_pvalue_cutoffs(
            [source], tmp_path / "ok.tsv", {"CGN": 0.01}, compress_output=False, verbose=False
        )


# --- merge_allc_tables ----------------------------------------------------------------------


def test_methylpy_merge_allc_tables_hand_sums_and_real_output_path(tmp_path: Path) -> None:
    a, b, _ = _methylpy_allc_inputs(tmp_path)
    result = methylpy_review.merge_allc_tables([a, b], tmp_path / "merged.tsv", verbose=False)
    # methylpy appends .gz when it compresses; the wrapper reports the file that exists.
    final = tmp_path / "merged.tsv.gz"
    assert result["output"] == str(final)
    assert final.is_file() and not (tmp_path / "merged.tsv").exists()
    # Position by position: methylated 9 + 1, 10 + 0, 1 (A only), 8 + 1, 5 + 5, and chr2 0 + 3;
    # totals 10 + 10, 10 + 10, 4, 10 + 10, 10 + 10 and 9 + 6.
    assert _methylpy_read_gz(final) == (
        "chr1\t100\t+\tCGA\t10\t20\t1\n"
        "chr1\t120\t+\tCGT\t10\t20\t1\n"
        "chr1\t125\t+\tCHH\t1\t4\t1\n"
        "chr1\t140\t+\tCGC\t9\t20\t1\n"
        "chr1\t500\t+\tCGG\t10\t20\t1\n"
        "chr2\t30\t+\tCGA\t3\t15\t1\n"
    )
    assert result["records"] == 6
    assert result["chromosomes"] == [
        {"chrom": "chr1", "records": 5},
        {"chrom": "chr2", "records": 1},
    ]
    offsets = methylpy_review.index_allc_table(final, verbose=False)["chromosomes"]
    chr1_bytes = len(_methylpy_read_gz(final).split("chr2")[0].encode())
    assert offsets == [{"chrom": "chr1", "offset": 0}, {"chrom": "chr2", "offset": chr1_bytes}]
    assert Path(result["index"]).is_file()


def test_methylpy_merge_allc_tables_plain_output_and_snp_columns(tmp_path: Path) -> None:
    utilities = pytest.importorskip("methylpy.utilities")
    first = _methylpy_write(
        tmp_path / "s1.tsv",
        "chr1\t10\t+\tCGA\t3\t5\t1\t2,1\t0,1\nchr1\t20\t+\tCGT\t1\t4\t0\t1,1\t0,0\n",
    )
    second = _methylpy_write(tmp_path / "s2.tsv", "chr1\t10\t+\tCGA\t2\t6\t1\t3,0\t1,2\n")
    result = methylpy_review.merge_allc_tables(
        [first, second],
        tmp_path / "snp.tsv",
        compress_output=False,
        skip_snp_info=False,
        verbose=False,
    )
    assert result["output"] == str(tmp_path / "snp.tsv")
    # Counts add up (3 + 2 of 5 + 6); the match and mismatch lists add element by element.
    assert (tmp_path / "snp.tsv").read_text() == (
        "chr1\t10\t+\tCGA\t5\t11\t1\t5,1\t1,3\nchr1\t20\t+\tCGT\t1\t4\t1\t1,1\t0,0\n"
    )
    direct = tmp_path / "direct.tsv"
    utilities.merge_allc_files(
        [str(first), str(second)], str(direct), compress_output=False, skip_snp_info=False
    )
    assert direct.read_text() == (tmp_path / "snp.tsv").read_text()
    plain = _methylpy_write(tmp_path / "plain.tsv", "chr1\t10\t+\tCGA\t1\t2\t1\n")
    with pytest.raises(ConfigurationError, match="fewer than 9"):
        methylpy_review.merge_allc_tables(
            [first, plain], tmp_path / "bad.tsv", skip_snp_info=False, verbose=False
        )


def test_methylpy_merge_allc_tables_order_does_not_depend_on_hash_seed(tmp_path: Path) -> None:
    a, b, _ = _methylpy_allc_inputs(tmp_path)
    code = (
        "import sys\n"
        "from dnakit.epigenome import methylpy_review_functions as m\n"
        "m.merge_allc_tables([sys.argv[1], sys.argv[2]], sys.argv[3], verbose=False)\n"
    )
    texts = []
    for seed in (1, 2, 3):
        output = tmp_path / f"merged_{seed}.tsv"
        for index_file in tmp_path.glob("*.idx"):
            index_file.unlink()
        subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code, str(a), str(b), str(output)],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            check=True,
            capture_output=True,
            text=True,
        )
        texts.append(_methylpy_read_gz(Path(str(output) + ".gz")))
    assert texts[0] == texts[1] == texts[2]
    assert [line.split("\t")[0] for line in texts[0].splitlines()] == ["chr1"] * 5 + ["chr2"]


def test_methylpy_merge_allc_tables_validation(tmp_path: Path) -> None:
    a, b, _ = _methylpy_allc_inputs(tmp_path)
    with pytest.raises(ConfigurationError, match="2 to"):
        methylpy_review.merge_allc_tables([a], tmp_path / "x.tsv", verbose=False)
    with pytest.raises(ConfigurationError, match="same file twice"):
        methylpy_review.merge_allc_tables([a, a], tmp_path / "x.tsv", verbose=False)
    with pytest.raises(ConfigurationError, match="mini_batch"):
        methylpy_review.merge_allc_tables([a, b], tmp_path / "x.tsv", mini_batch=1, verbose=False)
    with pytest.raises(ConfigurationError, match="gz"):
        methylpy_review.merge_allc_tables(
            [a, b], tmp_path / "x.gz", compress_output=False, verbose=False
        )
    interleaved = _methylpy_write(
        tmp_path / "bad.tsv",
        "chr1\t1\t+\tCGA\t1\t2\t0\nchr2\t1\t+\tCGA\t1\t2\t0\nchr1\t5\t+\tCGA\t1\t2\t0\n",
    )
    with pytest.raises(ConfigurationError, match="more than one block"):
        methylpy_review.merge_allc_tables([a, interleaved], tmp_path / "x.tsv", verbose=False)
    methylpy_review.merge_allc_tables([a, b], tmp_path / "x.tsv", verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.merge_allc_tables([a, b], tmp_path / "x.tsv", verbose=False)


# --- split_mpileup_chunks -------------------------------------------------------------------


def test_methylpy_split_mpileup_chunks_hand_overlap(tmp_path: Path) -> None:
    utilities = pytest.importorskip("methylpy.utilities")
    lines = [f"chr1\t{i}\t{base}\t5\t.....\n" for i, base in enumerate("CGACCTC", start=1)]
    source = _methylpy_write(tmp_path / "reads.mpileup", "".join(lines))
    result = methylpy_review.split_mpileup_chunks(source, tmp_path / "chunk_", 3, verbose=False)
    # Seven lines in three chunks give chunk_size = ceil(7 / 3) = 3. Lines 1-3 form chunk 0; chunk
    # 1 starts with lines 2-3 again (two overlap lines) and adds lines 4-6; chunk 2 repeats lines
    # 5-6 and adds line 7.
    assert result["records"] == 7
    assert result["chunk_size"] == 3
    assert [(c["records"], c["overlap_records"]) for c in result["chunks"]] == [
        (3, 0),
        (5, 2),
        (3, 2),
    ]
    chunk1 = (tmp_path / "chunk_1").read_text().splitlines()
    assert chunk1 == [line.rstrip("\n") for line in lines[1:6]]
    chunk2 = (tmp_path / "chunk_2").read_text().splitlines()
    assert chunk2 == [line.rstrip("\n") for line in lines[4:7]]
    direct = tmp_path / "direct_"
    utilities.split_mpileup_file(3, str(source), str(direct))
    assert (tmp_path / "direct_1").read_text() == (tmp_path / "chunk_1").read_text()


def test_methylpy_split_mpileup_chunks_rejects_bad_inputs(tmp_path: Path) -> None:
    empty = _methylpy_write(tmp_path / "empty.mpileup", "")
    with pytest.raises(ConfigurationError, match="methylpy rejected"):
        methylpy_review.split_mpileup_chunks(empty, tmp_path / "e_", 2, verbose=False)
    zipped = tmp_path / "reads.mpileup.gz"
    with gzip.open(zipped, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t1\tC\t5\t.....\n" * 5)
    with pytest.raises(ConfigurationError, match="plain text"):
        methylpy_review.split_mpileup_chunks(zipped, tmp_path / "z_", 2, verbose=False)
    source = _methylpy_write(tmp_path / "reads.mpileup", "chr1\t1\tC\t5\t.....\n" * 4)
    methylpy_review.split_mpileup_chunks(source, tmp_path / "c_", 2, verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        methylpy_review.split_mpileup_chunks(source, tmp_path / "c_", 2, verbose=False)
    with pytest.raises(ConfigurationError, match="num_chunks"):
        methylpy_review.split_mpileup_chunks(source, tmp_path / "d_", 0, verbose=False)
    with pytest.raises(ConfigurationError, match="prefix"):
        methylpy_review.split_mpileup_chunks(source, str(tmp_path) + "/", 2, verbose=False)
    assert not [name for name in _methylpy_tree(tmp_path) if name.startswith(("e_", "z_", "d_"))]


# --- upstream behaviour that justifies not wrapping fasta_iter ---------------------------------


def test_methylpy_fasta_iter_fails_on_python3(tmp_path: Path) -> None:
    call_mc = pytest.importorskip("methylpy.call_mc_se")
    fasta = _methylpy_write(tmp_path / "tiny.fa", ">chr1 note\nACGT\nGG\n>chr2\nTT\n")
    with pytest.raises(AttributeError, match="next"):
        call_mc.fasta_iter(str(fasta), "chr1")


_MACS3_CLASSES = ("short", "mono", "di", "tri")
_MACS3_MEANS = [50.0, 200.0, 400.0, 600.0]
_MACS3_STDDEVS = [20.0, 20.0, 20.0, 20.0]

# Probability file with columns chrom,end,state0,state1,state2. Here state 2 is open, state 0 is
# nucleosomal and state 1 is background; the bin size is 10.
_MACS3_PROBA = (
    "chr1,20,0.1,0.8,0.1\n"
    "chr1,30,0.8,0.1,0.1\n"
    "chr1,40,0.1,0.1,0.8\n"
    "chr1,50,0.1,0.1,0.8\n"
    "chr1,60,0.8,0.1,0.1\n"
    "chr1,70,0.1,0.8,0.1\n"
    "chr1,100,0.8,0.1,0.1\n"
    "chr1,110,0.1,0.1,0.8\n"
    "chr1,120,0.8,0.1,0.1\n"
    "chr2,40,0.1,0.8,0.1\n"
    "chr2,50,0.8,0.1,0.1\n"
    "chr2,60,0.1,0.1,0.8\n"
    "chr2,70,0.8,0.1,0.1\n"
)
_MACS3_INDEXES = {"open_state": 2, "nucleosomal_state": 0, "background_state": 1}

# The states path that hand-reading the file above gives: a bin without a neighbour on its left
# starts after a background gap, and a run of equal labels is one entry.
_MACS3_STATES = [
    ("chr1", 0, 10, "bg"),
    ("chr1", 10, 20, "bg"),
    ("chr1", 20, 30, "nuc"),
    ("chr1", 30, 50, "open"),
    ("chr1", 50, 60, "nuc"),
    ("chr1", 60, 70, "bg"),
    ("chr1", 70, 90, "bg"),
    ("chr1", 90, 100, "nuc"),
    ("chr1", 100, 110, "open"),
    ("chr1", 110, 120, "nuc"),
    ("chr2", 0, 30, "bg"),
    ("chr2", 30, 40, "bg"),
    ("chr2", 40, 50, "nuc"),
    ("chr2", 50, 60, "open"),
    ("chr2", 60, 70, "nuc"),
]


def _macs3_write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _macs3_fragments() -> list[tuple[str, int, int]]:
    """Fragments with short ones inside four open sites and mono-nucleosomal ones around them."""

    rng = random.Random(5)
    fragments: list[tuple[str, int, int]] = []

    def add(chrom: str, low: int, high: int, count: int, shortest: int, longest: int) -> None:
        for _ in range(count):
            length = rng.randint(shortest, longest)
            start = rng.randint(low, high - length)
            fragments.append((chrom, start, start + length))

    for chrom in ("chr1", "chr2"):
        for site in (2000, 6000):
            add(chrom, site - 600, site, 40, 150, 200)
            add(chrom, site - 600, site, 12, 370, 430)
            add(chrom, site, site + 400, 100, 40, 100)
            add(chrom, site + 400, site + 1000, 40, 150, 200)
            add(chrom, site + 400, site + 1100, 8, 570, 630)
    return fragments


def _macs3_flat(values: Any) -> list[float]:
    """Flatten nested lists of numbers so that pytest.approx can compare them."""

    if isinstance(values, (list, tuple)):
        return [number for item in values for number in _macs3_flat(item)]
    return [float(values)]


def _macs3_regions() -> list[tuple[str, int, int]]:
    return [(chrom, site - 800, site + 1200) for chrom in ("chr1", "chr2") for site in (2000, 6000)]


def _macs3_training(seed: int = 3) -> tuple[list[list[float]], list[int]]:
    """Three groups of 15 bins that cycle through background, nucleosomal and open levels."""

    rng = random.Random(seed)
    levels = {
        "bg": [0.3, 0.6, 0.3, 0.1],
        "nuc": [1.0, 6.0, 2.0, 0.5],
        "open": [8.0, 2.0, 0.8, 0.3],
    }
    pattern = ["bg"] * 4 + ["nuc"] * 3 + ["open"] * 4 + ["nuc"] * 2 + ["bg"] * 2
    rows = []
    for _ in range(3):
        for state in pattern:
            rows.append(
                [round(max(0.0001, level * (1 + rng.gauss(0, 0.12))), 4) for level in levels[state]]
            )
    return rows, [len(pattern)] * 3


def _macs3_upstream_digest(
    fragments: list[tuple[str, int, int]], means: list[float], stddevs: list[float]
) -> list[Any]:
    paired = pytest.importorskip("MACS3.Signal.PairedEndTrack")
    processing = pytest.importorskip("MACS3.Signal.HMMR_Signal_Processing")
    track = paired.PETrackI()
    for chrom, left, right in fragments:
        track.add_loc(chrom.encode(), left, right)
    track.finalize()
    mapping = processing.generate_weight_mapping(
        sorted(track.count_fraglengths()), means, stddevs, 0.001
    )
    return list(processing.generate_digested_signals(track, mapping))


def _macs3_model_json(path: Path, kind: str = "gaussian") -> Path:
    """A hand-written model: state 0 is background, state 1 nucleosomal and state 2 open."""

    model: dict[str, Any] = {
        "startprob": [0.2, 0.3, 0.5],
        "transmat": [[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]],
        "n_features": 4,
        "i_open_region": 2,
        "i_background_region": 0,
        "i_nucleosomal_region": 1,
        "hmm_binsize": 10,
        "hmm_type": kind,
    }
    if kind == "gaussian":
        model.update(
            {
                "means": [[1.0] * 4, [5.0] * 4, [10.0] * 4],
                "covars": [[1.0] * 4, [4.0] * 4, [9.0] * 4],
                "covariance_type": "diag",
            }
        )
    else:
        model["lambdas"] = [[1.0] * 4, [4.0] * 4, [9.0] * 4]
    return _macs3_write(path, json.dumps(model))


def _macs3_proba(tmp_path: Path) -> Path:
    return _macs3_write(tmp_path / "proba.csv", _MACS3_PROBA)


def _macs3_state_dicts(states: list[tuple[str, int, int, str]]) -> list[dict[str, Any]]:
    return [{"chrom": c, "start": s, "end": e, "state": t} for c, s, e, t in states]


# --- module contract ---------------------------------------------------------------------------


def test_macs3_review_public_functions_take_keyword_verbose() -> None:
    assert len(macs3_review.__all__) == 10
    for name in macs3_review.__all__:
        parameter = inspect.signature(getattr(macs3_review, name)).parameters.get("verbose")
        assert parameter is not None, name
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_macs3_review_names_do_not_collide_with_other_epigenome_modules() -> None:
    ours = set(macs3_review.__all__)
    for info in pkgutil.iter_modules(macs3_epigenome_package.__path__):
        if not info.name.endswith("_functions") or info.name == "MACS3_review_functions":
            continue
        other = importlib.import_module(f"dnakit.epigenome.{info.name}")
        assert not ours & set(getattr(other, "__all__", ())), info.name


def test_macs3_default_verbose_still_returns_the_result(tmp_path: Path) -> None:
    result = macs3_review.hmmratac_states_path(_macs3_proba(tmp_path), **_MACS3_INDEXES)
    assert result["status"] == "completed"
    assert result["operation"] == "hmmratac_states_path"
    assert result["provenance"]["package"].startswith("MACS3 ")


def test_macs3_missing_package_reports_install_extra(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    proba = _macs3_proba(tmp_path)
    for name in ("MACS3", "MACS3.Commands.hmmratac_cmd"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError, match="bio-macs3"):
        macs3_review.hmmratac_states_path(proba, **_MACS3_INDEXES, verbose=False)


def test_macs3_upstream_errors_become_configuration_errors() -> None:
    def failing_assert() -> None:
        raise AssertionError

    def failing_exit() -> None:
        sys.exit("bad input")

    with pytest.raises(ConfigurationError, match="rejected the input"):
        macs3_review._upstream_call(failing_assert)
    with pytest.raises(ConfigurationError, match="bad input"):
        macs3_review._upstream_call(failing_exit)


# --- digest_hmmratac_fragments -----------------------------------------------------------------


def test_macs3_digest_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.HMMR_Signal_Processing")
    stddevs = [40.0] * 4
    result = macs3_review.digest_hmmratac_fragments(
        [("chr1", 10, 135)], _MACS3_MEANS, stddevs, verbose=False
    )

    def density(x: float, mean: float, sd: float) -> float:
        return math.exp(-((x - mean) ** 2) / (2 * sd**2)) / (sd * math.sqrt(2 * math.pi))

    densities = [density(125, mean, 40.0) for mean in _MACS3_MEANS]
    for label, value in zip(_MACS3_CLASSES, densities, strict=True):
        track = result["tracks"][label]["chr1"]
        # Segments are [0, 10) at the baseline, then [10, 135) at the fragment's weight.
        assert track["positions"] == [10, 135]
        assert track["values"][0] == 0.0
        assert track["values"][1] == pytest.approx(value / sum(densities), rel=1e-6)
    # A 125 bp fragment is as far from the short mean as from the mono mean.
    assert result["tracks"]["short"]["chr1"]["values"][1] == pytest.approx(0.5, abs=1e-6)
    assert result["tracks"]["mono"]["chr1"]["values"][1] == pytest.approx(0.5, abs=1e-6)
    assert result["intervals"] == {"short": 2, "mono": 2, "di": 2, "tri": 2}
    assert result["fragments"] == 1
    assert result["chromosomes"] == ["chr1"]


def test_macs3_digest_matches_upstream() -> None:
    fragments = _macs3_fragments()
    result = macs3_review.digest_hmmratac_fragments(
        fragments, _MACS3_MEANS, _MACS3_STDDEVS, verbose=False
    )
    signals = _macs3_upstream_digest(fragments, _MACS3_MEANS, _MACS3_STDDEVS)
    assert len(signals) == 4
    for label, signal in zip(_MACS3_CLASSES, signals, strict=True):
        names = sorted(signal.get_chr_names())
        assert sorted(result["tracks"][label]) == [name.decode() for name in names]
        for name in names:
            positions, values = signal.get_data_by_chr(name)
            track = result["tracks"][label][name.decode()]
            assert track["positions"] == list(positions)
            assert track["values"] == list(values)
    assert result["fragments"] == len(fragments)
    assert result["chromosomes"] == ["chr1", "chr2"]
    assert result["intervals"]["short"] == sum(
        len(track["positions"]) for track in result["tracks"]["short"].values()
    )


def test_macs3_digest_writes_four_bedgraph_files(tmp_path: Path) -> None:
    bedgraph = pytest.importorskip("MACS3.IO.BedGraphIO")
    fragments = [("chr1", 10, 60), ("chr1", 30, 230)]
    prefix = tmp_path / "sample"
    result = macs3_review.digest_hmmratac_fragments(
        fragments,
        _MACS3_MEANS,
        _MACS3_STDDEVS,
        output_prefix=prefix,
        include_tracks=False,
        verbose=False,
    )
    assert "tracks" not in result
    assert result["files"] == {label: f"{prefix}_{label}.bdg" for label in _MACS3_CLASSES}
    short = Path(result["files"]["short"]).read_text(encoding="utf-8").splitlines()
    assert short[0].startswith('track type=bedGraph name="short"')
    assert short[1:] == ["chr1\t0\t10\t0.00000", "chr1\t10\t60\t1.00000", "chr1\t60\t230\t0.00000"]
    # The mono weight of the short fragment is 6e-13, which prints as 0.00000. Reading the file
    # back merges the neighbouring equal values.
    mono = Path(result["files"]["mono"]).read_text(encoding="utf-8").splitlines()
    assert mono[1:] == ["chr1\t0\t10\t0.00000", "chr1\t10\t30\t0.00000", "chr1\t30\t230\t1.00000"]
    track = bedgraph.bedGraphIO(result["files"]["mono"]).read_bedGraph(0)
    positions, values = track.get_data_by_chr(b"chr1")
    assert list(positions) == [30, 230]
    assert list(values) == pytest.approx([0.0, 1.0], abs=1e-5)
    with pytest.raises(ConfigurationError, match="already exists"):
        macs3_review.digest_hmmratac_fragments(
            fragments, _MACS3_MEANS, _MACS3_STDDEVS, output_prefix=prefix, verbose=False
        )


def test_macs3_digest_validation() -> None:
    good = [("chr1", 10, 60)]

    def call(
        fragments: Any = good, means: Any = _MACS3_MEANS, stddevs: Any = _MACS3_STDDEVS
    ) -> Any:
        return macs3_review.digest_hmmratac_fragments(fragments, means, stddevs, verbose=False)

    with pytest.raises(ConfigurationError, match="at least one"):
        call(fragments=[])
    with pytest.raises(ConfigurationError, match="must be a list"):
        call(fragments="chr1")
    with pytest.raises(ConfigurationError, match="triple"):
        call(fragments=[("chr1", 10)])
    with pytest.raises(ConfigurationError, match="chromosome name"):
        call(fragments=[("chr 1", 10, 60)])
    with pytest.raises(ConfigurationError, match="chromosome name"):
        call(fragments=[("chr1,x", 10, 60)])
    with pytest.raises(ConfigurationError, match="right coordinate"):
        call(fragments=[("chr1", 60, 60)])
    with pytest.raises(ConfigurationError, match="at most 10000"):
        call(fragments=[("chr1", 0, 10_001)])
    with pytest.raises(ConfigurationError, match="integer"):
        call(fragments=[("chr1", 1.5, 60)])
    with pytest.raises(ConfigurationError, match="four values"):
        call(means=[50.0, 200.0])
    with pytest.raises(ConfigurationError, match="stddevs entries"):
        call(stddevs=[20.0, 20.0, 20.0, 0.0])
    with pytest.raises(ConfigurationError, match="include_tracks"):
        macs3_review.digest_hmmratac_fragments(
            good, _MACS3_MEANS, _MACS3_STDDEVS, include_tracks="yes", verbose=False
        )


# --- extract_hmmratac_training_data ------------------------------------------------------------


def test_macs3_extract_hand_calculation() -> None:
    pytest.importorskip("MACS3.Signal.HMMR_Signal_Processing")
    # Two 50 bp fragments overlap on 120-150, so the short track is 1, 2 and 1 over 100-170.
    result = macs3_review.extract_hmmratac_training_data(
        [("chr1", 100, 150), ("chr1", 120, 170)],
        [("chr1", 100, 170)],
        _MACS3_MEANS,
        _MACS3_STDDEVS,
        binsize=10,
        verbose=False,
    )
    assert result["bins"] == [["chr1", end] for end in range(110, 171, 10)]
    assert result["lengths"] == [7]
    assert result["regions"] == 1
    short = [row[0] for row in result["data"]]
    assert short == pytest.approx([1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0], abs=1e-5)
    # The other three tracks are 0 here and are raised to 0.0001 like the pipeline does.
    assert all(row[1:] == [0.0001, 0.0001, 0.0001] for row in result["data"])


def test_macs3_extract_matches_upstream() -> None:
    region_module = pytest.importorskip("MACS3.Signal.Region")
    processing = pytest.importorskip("MACS3.Signal.HMMR_Signal_Processing")
    fragments = _macs3_fragments()
    regions = _macs3_regions()
    result = macs3_review.extract_hmmratac_training_data(
        fragments, regions, _MACS3_MEANS, _MACS3_STDDEVS, binsize=20, flank=100, verbose=False
    )
    signals = _macs3_upstream_digest(fragments, _MACS3_MEANS, _MACS3_STDDEVS)
    expected = region_module.Regions()
    for chrom, start, end in regions:
        expected.add_loc(chrom.encode(), start, end)
    expected.expand(100)
    expected.merge_overlap()
    bins, data, lengths = processing.extract_signals_from_regions(
        signals, expected, binsize=20, hmm_type="gaussian"
    )
    assert result["bins"] == [[chrom.decode(), end] for chrom, end in bins]
    assert result["data"] == [[float(item) for item in row] for row in data]
    assert result["lengths"] == list(lengths)
    assert sum(result["lengths"]) == len(result["bins"])
    assert len(result["lengths"]) == result["regions"] == 4
    assert next(name for name, _ in result["bins"]) == "chr2"


def test_macs3_extract_merges_overlapping_regions_and_flank() -> None:
    fragments = [("chr1", 100, 150), ("chr1", 120, 170), ("chr1", 400, 450)]
    split = macs3_review.extract_hmmratac_training_data(
        fragments,
        [("chr1", 100, 130), ("chr1", 120, 170)],
        _MACS3_MEANS,
        _MACS3_STDDEVS,
        verbose=False,
    )
    # The overlapping regions are one region of the bins 110 to 170.
    assert split["regions"] == 1
    assert split["lengths"] == [7]
    wide = macs3_review.extract_hmmratac_training_data(
        fragments,
        [("chr1", 150, 170), ("chr1", 400, 420)],
        _MACS3_MEANS,
        _MACS3_STDDEVS,
        flank=50,
        verbose=False,
    )
    # With 50 bp of flank the regions become 100-220 and 350-470. The signal stops at 450, so the
    # second region keeps the bins ending 360 to 450.
    assert wide["regions"] == 2
    assert wide["lengths"] == [12, 10]
    assert wide["bins"][0] == ["chr1", 110]
    assert wide["bins"][11] == ["chr1", 220]
    assert wide["bins"][12] == ["chr1", 360]


def test_macs3_extract_poisson_returns_whole_numbers() -> None:
    result = macs3_review.extract_hmmratac_training_data(
        [("chr1", 100, 150), ("chr1", 120, 170)],
        [("chr1", 100, 170)],
        _MACS3_MEANS,
        _MACS3_STDDEVS,
        hmm_type="poisson",
        verbose=False,
    )
    # The Poisson model needs whole numbers, so MACS3 truncates 0.0001 to 0 and keeps 1 and 2.
    assert all(isinstance(item, int) for row in result["data"] for item in row)
    assert result["data"] == [[value, 0, 0, 0] for value in (1, 1, 2, 2, 2, 1, 1)]
    assert result["hmm_type"] == "poisson"


def test_macs3_extract_validation() -> None:
    fragments = [("chr1", 100, 150)]
    with pytest.raises(ConfigurationError, match="No signal bin"):
        macs3_review.extract_hmmratac_training_data(
            fragments, [("chr9", 0, 500)], _MACS3_MEANS, _MACS3_STDDEVS, verbose=False
        )
    with pytest.raises(ConfigurationError, match="No signal bin"):
        macs3_review.extract_hmmratac_training_data(
            fragments, [("chr1", 5000, 6000)], _MACS3_MEANS, _MACS3_STDDEVS, verbose=False
        )
    for bad, message in (
        ({"regions": []}, "at least one"),
        ({"regions": [("chr1", 50, 40)]}, "end after it starts"),
        ({"regions": [("chr1", 0)]}, "triple"),
        ({"binsize": 0}, "binsize"),
        ({"flank": -1}, "flank"),
        ({"hmm_type": "beta"}, "hmm_type"),
        ({"min_fragment_probability": 2.0}, "min_fragment_probability"),
    ):
        options: dict[str, Any] = {"regions": [("chr1", 100, 150)], **bad}
        with pytest.raises(ConfigurationError, match=message):
            macs3_review.extract_hmmratac_training_data(
                fragments,
                options.pop("regions"),
                _MACS3_MEANS,
                _MACS3_STDDEVS,
                **options,
                verbose=False,
            )


def test_macs3_digest_and_extract_do_not_depend_on_hash_seed() -> None:
    code = (
        "import json\n"
        "from dnakit.epigenome import MACS3_review_functions as m\n"
        "frags = [('chr1', 100, 150), ('chr2', 120, 170), ('chr10', 90, 140), ('chr1', 110, 160)]\n"
        "regions = [('chr1', 100, 170), ('chr2', 100, 170), ('chr10', 90, 170)]\n"
        "d = m.digest_hmmratac_fragments(frags, [50, 200, 400, 600], [20] * 4, verbose=False)\n"
        "e = m.extract_hmmratac_training_data(\n"
        "    frags, regions, [50, 200, 400, 600], [20] * 4, verbose=False\n"
        ")\n"
        "print(json.dumps([d, e], sort_keys=True))\n"
    )
    outputs = []
    for seed in (1, 2, 3):
        done = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            check=True,
            capture_output=True,
            text=True,
        )
        outputs.append(done.stdout.strip().splitlines()[-1])
    assert outputs[0] == outputs[1] == outputs[2]
    chromosomes = [name for name, _ in json.loads(outputs[0])[1]["bins"]]
    assert sorted(set(chromosomes), reverse=True) == list(dict.fromkeys(chromosomes))


# --- train_hmmratac_model ----------------------------------------------------------------------


def test_macs3_train_matches_upstream_and_labels_states() -> None:
    model_module = pytest.importorskip("MACS3.Signal.HMMR_HMM")
    rows, lengths = _macs3_training()
    result = macs3_review.train_hmmratac_model(rows, lengths, verbose=False)
    model = model_module.hmm_training(rows, lengths, random_seed=12345)
    assert result["converged"] is bool(model.monitor_.converged)
    assert result["iterations"] == model.monitor_.iter
    assert result["covariance_type"] == "full"
    assert _macs3_flat(result["means"]) == pytest.approx(
        _macs3_flat(model.means_.tolist()), rel=1e-6, abs=1e-8
    )
    assert _macs3_flat(result["transition_matrix"]) == pytest.approx(
        _macs3_flat(model.transmat_.tolist()), rel=1e-6, abs=1e-8
    )
    assert result["start_probabilities"] == pytest.approx(model.startprob_.tolist(), abs=1e-8)
    assert len(result["covariances"]) == 3
    assert len(result["covariances"][0]) == 4
    # The open state has the largest summed emission, the background state the smallest.
    totals = [sum(row) for row in result["means"]]
    assert totals.index(max(totals)) == result["open_state"]
    assert totals.index(min(totals)) == result["background_state"]
    assert {
        result["open_state"],
        result["nucleosomal_state"],
        result["background_state"],
    } == {0, 1, 2}
    for row in result["transition_matrix"]:
        assert sum(row) == pytest.approx(1.0)
    assert sum(result["start_probabilities"]) == pytest.approx(1.0)


def test_macs3_train_is_seeded_and_leaves_global_random_state_alone() -> None:
    numpy = pytest.importorskip("numpy")
    rows, lengths = _macs3_training()
    random.seed(11)
    numpy.random.seed(11)
    python_state = random.getstate()
    numpy_state = numpy.random.get_state()
    first = macs3_review.train_hmmratac_model(rows, lengths, random_seed=7, verbose=False)
    second = macs3_review.train_hmmratac_model(rows, lengths, random_seed=7, verbose=False)
    assert random.getstate() == python_state
    assert numpy.random.get_state()[1].tolist() == numpy_state[1].tolist()
    assert numpy.random.get_state()[2:] == numpy_state[2:]
    assert _macs3_flat(first["means"]) == pytest.approx(
        _macs3_flat(second["means"]), rel=1e-6, abs=1e-8
    )
    assert first["open_state"] == second["open_state"]


def test_macs3_train_diag_and_poisson_variants() -> None:
    rows, lengths = _macs3_training()
    diag = macs3_review.train_hmmratac_model(rows, lengths, covariance_type="diag", verbose=False)
    assert diag["covariance_type"] == "diag"
    off_diagonal = diag["covariances"][0][0][1]
    assert off_diagonal == 0.0
    counts = [[round(item) for item in row] for row in rows]
    poisson = macs3_review.train_hmmratac_model(counts, lengths, hmm_type="poisson", verbose=False)
    assert poisson["hmm_type"] == "poisson"
    assert "means" not in poisson
    assert len(poisson["lambdas"]) == 3
    totals = [sum(row) for row in poisson["lambdas"]]
    assert totals.index(max(totals)) == poisson["open_state"]
    assert totals.index(min(totals)) == poisson["background_state"]


def test_macs3_train_validation() -> None:
    rows, lengths = _macs3_training()
    with pytest.raises(ConfigurationError, match="four values"):
        macs3_review.train_hmmratac_model([[1.0, 2.0, 3.0]] * 6, [6], verbose=False)
    with pytest.raises(ConfigurationError, match="add up"):
        macs3_review.train_hmmratac_model(rows, [10, 10], verbose=False)
    with pytest.raises(ConfigurationError, match="at least one row"):
        macs3_review.train_hmmratac_model([], [], verbose=False)
    with pytest.raises(ConfigurationError, match="data entries"):
        macs3_review.train_hmmratac_model([[1.0, 2.0, 3.0, -1.0]] * 6, [6], verbose=False)
    with pytest.raises(ConfigurationError, match="data entries"):
        macs3_review.train_hmmratac_model([[1.0, 2.0, 3.0, math.nan]] * 6, [6], verbose=False)
    with pytest.raises(ConfigurationError, match="integer"):
        macs3_review.train_hmmratac_model(rows, lengths, hmm_type="poisson", verbose=False)
    with pytest.raises(ConfigurationError, match="hmm_type"):
        macs3_review.train_hmmratac_model(rows, lengths, hmm_type="beta", verbose=False)
    with pytest.raises(ConfigurationError, match="covariance_type"):
        macs3_review.train_hmmratac_model(rows, lengths, covariance_type="tied", verbose=False)
    with pytest.raises(ConfigurationError, match="random_seed"):
        macs3_review.train_hmmratac_model(rows, lengths, random_seed=-1, verbose=False)
    # hmmlearn cannot fit three states to two rows, or to identical rows.
    with pytest.raises(ConfigurationError, match="rejected the input"):
        macs3_review.train_hmmratac_model([[1.0, 1.0, 1.0, 1.0]] * 2, [2], verbose=False)
    with pytest.raises(ConfigurationError, match=r"rejected the input|equal emissions"):
        macs3_review.train_hmmratac_model([[1.0, 1.0, 1.0, 1.0]] * 30, [30], verbose=False)


# --- save_hmmratac_model and load_hmmratac_model ------------------------------------------------


def test_macs3_save_and_load_round_trip(tmp_path: Path) -> None:
    model_module = pytest.importorskip("MACS3.Signal.HMMR_HMM")
    rows, lengths = _macs3_training()
    path = tmp_path / "model.json"
    saved = macs3_review.save_hmmratac_model(rows, lengths, path, binsize=25, verbose=False)
    trained = macs3_review.train_hmmratac_model(rows, lengths, verbose=False)
    assert saved["output"] == str(path)
    assert saved["binsize"] == 25
    assert saved["converged"] is trained["converged"]
    file_content = json.loads(path.read_text(encoding="utf-8"))
    assert file_content["hmm_binsize"] == 25
    assert file_content["hmm_type"] == "gaussian"
    assert file_content["n_features"] == 4
    assert file_content["i_open_region"] == saved["open_state"] == trained["open_state"]
    assert file_content["i_nucleosomal_region"] == saved["nucleosomal_state"]
    assert file_content["i_background_region"] == saved["background_state"]
    loaded = macs3_review.load_hmmratac_model(path, verbose=False)
    upstream, i_open, i_background, i_nucleosomal, binsize, kind = model_module.hmm_model_init(
        str(path)
    )
    assert (loaded["open_state"], loaded["nucleosomal_state"], loaded["background_state"]) == (
        i_open,
        i_nucleosomal,
        i_background,
    )
    assert loaded["binsize"] == binsize == 25
    assert loaded["hmm_type"] == kind == "gaussian"
    assert loaded["means"] == upstream.means_.tolist()
    assert loaded["covariances"] == upstream.covars_.tolist()
    assert loaded["transition_matrix"] == upstream.transmat_.tolist()
    assert loaded["start_probabilities"] == upstream.startprob_.tolist()
    assert _macs3_flat(loaded["means"]) == pytest.approx(
        _macs3_flat(trained["means"]), rel=1e-6, abs=1e-8
    )
    assert sorted(loaded["state_labels"]) == ["bg", "nuc", "open"]
    assert loaded["state_labels"][loaded["open_state"]] == "open"
    assert loaded["features"] == list(_MACS3_CLASSES)
    with pytest.raises(ConfigurationError, match="already exists"):
        macs3_review.save_hmmratac_model(rows, lengths, path, verbose=False)


def test_macs3_save_matches_the_upstream_model_file(tmp_path: Path) -> None:
    model_module = pytest.importorskip("MACS3.Signal.HMMR_HMM")
    rows, lengths = _macs3_training()
    ours = tmp_path / "ours.json"
    macs3_review.save_hmmratac_model(rows, lengths, ours, binsize=15, random_seed=99, verbose=False)
    model = model_module.hmm_training(rows, lengths, random_seed=99)
    sums = model.means_.sum(axis=1)
    i_open = int(sums.argmax())
    i_background = int(sums.argmin())
    i_nucleosomal = ({0, 1, 2} - {i_open, i_background}).pop()
    theirs = tmp_path / "theirs.json"
    model_module.hmm_model_save(
        str(theirs), model, 15, i_open, i_nucleosomal, i_background, "gaussian"
    )
    mine = json.loads(ours.read_text(encoding="utf-8"))
    expected = json.loads(theirs.read_text(encoding="utf-8"))
    assert sorted(mine) == sorted(expected)
    for key in ("i_open_region", "i_nucleosomal_region", "i_background_region", "hmm_binsize"):
        assert mine[key] == expected[key]
    for key in ("n_features", "hmm_type", "covariance_type"):
        assert mine[key] == expected[key]
    for key in ("startprob", "transmat", "means", "covars"):
        assert _macs3_flat(mine[key]) == pytest.approx(
            _macs3_flat(expected[key]), rel=1e-6, abs=1e-8
        )


def test_macs3_save_and_load_diag_and_poisson(tmp_path: Path) -> None:
    rows, lengths = _macs3_training()
    diag_path = tmp_path / "diag.json"
    macs3_review.save_hmmratac_model(
        rows, lengths, diag_path, covariance_type="diag", verbose=False
    )
    diag = macs3_review.load_hmmratac_model(diag_path, verbose=False)
    assert diag["covariance_type"] == "diag"
    assert json.loads(diag_path.read_text(encoding="utf-8"))["covars"][0][0] == pytest.approx(
        diag["covariances"][0][0][0]
    )
    counts = [[round(item) for item in row] for row in rows]
    poisson_path = tmp_path / "poisson.json"
    macs3_review.save_hmmratac_model(
        counts, lengths, poisson_path, hmm_type="poisson", binsize=5, verbose=False
    )
    poisson = macs3_review.load_hmmratac_model(poisson_path, verbose=False)
    assert poisson["hmm_type"] == "poisson"
    assert poisson["binsize"] == 5
    assert len(poisson["lambdas"]) == 3
    assert "means" not in poisson


def test_macs3_load_reads_a_hand_written_model_and_old_files(tmp_path: Path) -> None:
    path = _macs3_model_json(tmp_path / "hand.json")
    loaded = macs3_review.load_hmmratac_model(path, verbose=False)
    assert loaded["open_state"] == 2
    assert loaded["nucleosomal_state"] == 1
    assert loaded["background_state"] == 0
    assert loaded["state_labels"] == ["bg", "nuc", "open"]
    assert loaded["start_probabilities"] == [0.2, 0.3, 0.5]
    assert loaded["means"][2] == [10.0] * 4
    assert loaded["covariances"][1][0] == [4.0, 0.0, 0.0, 0.0]
    # A file without hmm_type was written by an older MACS3 and is read as Gaussian.
    content = json.loads(path.read_text(encoding="utf-8"))
    del content["hmm_type"]
    old = _macs3_write(tmp_path / "old.json", json.dumps(content))
    assert macs3_review.load_hmmratac_model(old, verbose=False)["hmm_type"] == "gaussian"


def test_macs3_load_and_save_validation(tmp_path: Path) -> None:
    good = json.loads(_macs3_model_json(tmp_path / "good.json").read_text(encoding="utf-8"))

    def broken(name: str, **changes: Any) -> Path:
        content = {**good, **changes}
        content = {key: value for key, value in content.items() if value is not None}
        return _macs3_write(tmp_path / f"{name}.json", json.dumps(content))

    with pytest.raises(ConfigurationError, match="does not exist"):
        macs3_review.load_hmmratac_model(tmp_path / "missing.json", verbose=False)
    with pytest.raises(ConfigurationError, match="not a JSON model file"):
        macs3_review.load_hmmratac_model(_macs3_write(tmp_path / "x.json", "{"), verbose=False)
    with pytest.raises(ConfigurationError, match="JSON object"):
        macs3_review.load_hmmratac_model(_macs3_write(tmp_path / "y.json", "[1]"), verbose=False)
    cases = (
        ("kind", {"hmm_type": "beta"}, "hmm_type"),
        ("start", {"startprob": [1.0, 0.0]}, "startprob"),
        ("trans", {"transmat": [[1.0, 0.0, 0.0]] * 2}, "transmat"),
        ("means", {"means": None}, "means"),
        ("features", {"n_features": 3}, "n_features"),
        ("cov", {"covariance_type": "tied"}, "covariance_type"),
        ("index", {"i_open_region": 1}, "state indexes"),
        ("flag", {"i_open_region": True}, "state indexes"),
        ("bin", {"hmm_binsize": 0}, "hmm_binsize"),
        ("nan", {"means": [[1.0, 1.0, 1.0, "x"]] * 3}, "means"),
    )
    for name, changes, message in cases:
        with pytest.raises(ConfigurationError, match=message):
            macs3_review.load_hmmratac_model(broken(name, **changes), verbose=False)
    # A covariance that is not positive definite is refused by hmmlearn, not by the file check.
    flat = broken("flat", covars=[[0.0] * 4] * 3)
    with pytest.raises(ConfigurationError, match="rejected the input"):
        macs3_review.load_hmmratac_model(flat, verbose=False)
    rows, lengths = _macs3_training()
    with pytest.raises(ConfigurationError, match="binsize"):
        macs3_review.save_hmmratac_model(
            rows, lengths, tmp_path / "m.json", binsize=0, verbose=False
        )
    with pytest.raises(ConfigurationError, match="directory"):
        macs3_review.save_hmmratac_model(
            rows, lengths, tmp_path / "nodir" / "m.json", verbose=False
        )
    assert not (tmp_path / "m.json").exists()


# --- predict_hmmratac_states -------------------------------------------------------------------


def test_macs3_predict_gaussian_hand_calculation(tmp_path: Path) -> None:
    path = _macs3_model_json(tmp_path / "hand.json")
    observed = [2.0, 2.0, 2.0, 2.0]
    result = macs3_review.predict_hmmratac_states(path, [observed], [1], verbose=False)
    # With one bin the posterior is start probability times the diagonal Gaussian density.
    weights = []
    for start, mean, variance in ((0.2, 1.0, 1.0), (0.3, 5.0, 4.0), (0.5, 10.0, 9.0)):
        density = 1.0
        for value in observed:
            density *= math.exp(-((value - mean) ** 2) / (2 * variance)) / math.sqrt(
                2 * math.pi * variance
            )
        weights.append(start * density)
    posterior = [weight / sum(weights) for weight in weights]
    # Columns are open (state 2), nucleosomal (state 1) and background (state 0).
    assert result["probabilities"][0] == pytest.approx(
        [posterior[2], posterior[1], posterior[0]], rel=1e-9
    )
    assert result["labels"] == ["bg"]
    assert result["label_counts"] == {"open": 0, "nuc": 0, "bg": 1}
    assert result["rows"] == 1
    assert result["columns"] == ["open", "nucleosomal", "background"]
    assert (result["open_state"], result["nucleosomal_state"], result["background_state"]) == (
        2,
        1,
        0,
    )


def test_macs3_predict_poisson_hand_calculation(tmp_path: Path) -> None:
    path = _macs3_model_json(tmp_path / "pois.json", kind="poisson")
    result = macs3_review.predict_hmmratac_states(path, [[3, 3, 3, 3]], [1], verbose=False)
    weights = []
    for start, rate in ((0.2, 1.0), (0.3, 4.0), (0.5, 9.0)):
        density = (math.exp(-rate) * rate**3 / math.factorial(3)) ** 4
        weights.append(start * density)
    posterior = [weight / sum(weights) for weight in weights]
    assert result["hmm_type"] == "poisson"
    assert result["probabilities"][0] == pytest.approx(
        [posterior[2], posterior[1], posterior[0]], rel=1e-9
    )
    assert result["labels"] == ["nuc"]
    with pytest.raises(ConfigurationError, match="integer"):
        macs3_review.predict_hmmratac_states(path, [[3.5, 3, 3, 3]], [1], verbose=False)


def test_macs3_predict_matches_upstream_and_writes_probability_file(tmp_path: Path) -> None:
    model_module = pytest.importorskip("MACS3.Signal.HMMR_HMM")
    rows, lengths = _macs3_training()
    model_path = tmp_path / "model.json"
    macs3_review.save_hmmratac_model(rows, lengths, model_path, verbose=False)
    bins = [["chr1", 10 * (index + 1)] for index in range(30)] + [
        ["chr2", 10 * (index + 1)] for index in range(15)
    ]
    output = tmp_path / "proba.csv"
    result = macs3_review.predict_hmmratac_states(
        model_path, rows, lengths, bins=bins, output_path=output, verbose=False
    )
    upstream, i_open, _i_background, i_nucleosomal, _binsize, _kind = model_module.hmm_model_init(
        str(model_path)
    )
    raw = model_module.hmm_predict(rows, lengths, upstream)
    i_background = ({0, 1, 2} - {i_open, i_nucleosomal}).pop()
    assert result["rows"] == 45
    assert len(result["probabilities"]) == len(result["labels"]) == 45
    for row, expected in zip(result["probabilities"], raw, strict=True):
        assert row == [
            float(expected[i_open]),
            float(expected[i_nucleosomal]),
            float(expected[i_background]),
        ]
        assert sum(row) == pytest.approx(1.0)
    assert sum(result["label_counts"].values()) == 45
    # The file keeps the model's own state order, with six decimals, in the MACS3 pipeline format.
    lines = output.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 45
    assert lines[0].startswith("chr1,10,")
    assert lines[30].startswith("chr2,10,")
    first = lines[0].split(",")
    assert [float(item) for item in first[2:]] == pytest.approx(raw[0].tolist(), abs=5e-7)
    assert all(len(item.split(".")[1]) == 6 for item in first[2:])
    assert result["output"] == str(output)


def test_macs3_predict_probability_file_feeds_state_calls(tmp_path: Path) -> None:
    rows, lengths = _macs3_training()
    model_path = tmp_path / "model.json"
    saved = macs3_review.save_hmmratac_model(rows, lengths, model_path, verbose=False)
    bins = [["chr1", 10 * (index + 1)] for index in range(45)]
    output = tmp_path / "proba.csv"
    predicted = macs3_review.predict_hmmratac_states(
        model_path, rows, [45], bins=bins, output_path=output, verbose=False
    )
    indexes = {
        "open_state": saved["open_state"],
        "nucleosomal_state": saved["nucleosomal_state"],
        "background_state": saved["background_state"],
    }
    states = macs3_review.hmmratac_states_path(output, **indexes, verbose=False)
    # Neighbouring bins with the same most-probable label are merged into one entry.
    assert states["states"][0]["start"] == 0
    assert states["states"][-1]["end"] == 450
    assert sum(item["end"] - item["start"] for item in states["states"]) == 450
    for label in ("open", "nuc", "bg"):
        assert states["state_bases"][label] == 10 * predicted["label_counts"][label]


def test_macs3_predict_validation(tmp_path: Path) -> None:
    path = _macs3_model_json(tmp_path / "hand.json")
    rows = [[1.0, 1.0, 1.0, 1.0]] * 3
    with pytest.raises(ConfigurationError, match="go together"):
        macs3_review.predict_hmmratac_states(
            path, rows, [3], output_path=tmp_path / "p.csv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="go together"):
        macs3_review.predict_hmmratac_states(
            path, rows, [3], bins=[["chr1", 10]] * 3, verbose=False
        )
    with pytest.raises(ConfigurationError, match="one \\(chrom, end\\) pair"):
        macs3_review.predict_hmmratac_states(
            path, rows, [3], bins=[["chr1", 10]], output_path=tmp_path / "p.csv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="pair"):
        macs3_review.predict_hmmratac_states(
            path, rows, [3], bins=[["chr1"]] * 3, output_path=tmp_path / "p.csv", verbose=False
        )
    with pytest.raises(ConfigurationError, match="add up"):
        macs3_review.predict_hmmratac_states(path, rows, [2], verbose=False)
    with pytest.raises(ConfigurationError, match="four values"):
        macs3_review.predict_hmmratac_states(path, [[1.0, 2.0]], [1], verbose=False)
    existing = _macs3_write(tmp_path / "exists.csv", "")
    with pytest.raises(ConfigurationError, match="already exists"):
        macs3_review.predict_hmmratac_states(
            path, rows, [3], bins=[["chr1", 10]] * 3, output_path=existing, verbose=False
        )
    assert not (tmp_path / "p.csv").exists()


# --- hmmratac_proba_to_bedgraph ----------------------------------------------------------------


def test_macs3_proba_to_bedgraph_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    proba = _macs3_proba(tmp_path)
    result = macs3_review.hmmratac_proba_to_bedgraph(
        proba, tmp_path / "s", **_MACS3_INDEXES, verbose=False
    )
    assert result["files"] == {
        "open": str(tmp_path / "s_open.bdg"),
        "nucleosomal": str(tmp_path / "s_nuc.bdg"),
        "background": str(tmp_path / "s_bg.bdg"),
    }
    assert result["bin_count"] == 13
    assert result["chromosomes"] == ["chr1", "chr2"]
    # Open is column 2. The gap before the first bin is 0.0, equal neighbours merge.
    open_lines = (tmp_path / "s_open.bdg").read_text(encoding="utf-8").splitlines()
    assert open_lines == [
        "chr1\t0\t10\t0.00000",
        "chr1\t10\t30\t0.10000",
        "chr1\t30\t50\t0.80000",
        "chr1\t50\t70\t0.10000",
        "chr1\t70\t90\t0.00000",
        "chr1\t90\t100\t0.10000",
        "chr1\t100\t110\t0.80000",
        "chr1\t110\t120\t0.10000",
        "chr2\t0\t30\t0.00000",
        "chr2\t30\t50\t0.10000",
        "chr2\t50\t60\t0.80000",
        "chr2\t60\t70\t0.10000",
    ]
    # Background is column 1 and is 1.0 in the gaps.
    background = (tmp_path / "s_bg.bdg").read_text(encoding="utf-8").splitlines()
    assert background == [
        "chr1\t0\t10\t1.00000",
        "chr1\t10\t20\t0.80000",
        "chr1\t20\t60\t0.10000",
        "chr1\t60\t70\t0.80000",
        "chr1\t70\t90\t1.00000",
        "chr1\t90\t120\t0.10000",
        "chr2\t0\t30\t1.00000",
        "chr2\t30\t40\t0.80000",
        "chr2\t40\t70\t0.10000",
    ]
    # Nucleosomal is column 0.
    nucleosomal = (tmp_path / "s_nuc.bdg").read_text(encoding="utf-8").splitlines()
    assert nucleosomal[:3] == [
        "chr1\t0\t10\t0.00000",
        "chr1\t10\t20\t0.10000",
        "chr1\t20\t30\t0.80000",
    ]
    assert result["intervals"] == {
        "open": len(open_lines),
        "nucleosomal": len(nucleosomal),
        "background": len(background),
    }


def test_macs3_proba_to_bedgraph_matches_upstream(tmp_path: Path) -> None:
    command = pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    proba = _macs3_proba(tmp_path)
    macs3_review.hmmratac_proba_to_bedgraph(
        proba,
        tmp_path / "ours",
        open_state=0,
        nucleosomal_state=2,
        background_state=1,
        verbose=False,
    )
    with proba.open("rb") as handle:
        command.save_proba_to_bedGraph(
            handle,
            10,
            str(tmp_path / "up_open.bdg"),
            str(tmp_path / "up_nuc.bdg"),
            str(tmp_path / "up_bg.bdg"),
            0,
            2,
            1,
        )
    for ours, theirs in (("open", "open"), ("nuc", "nuc"), ("bg", "bg")):
        assert (tmp_path / f"ours_{ours}.bdg").read_bytes() == (
            tmp_path / f"up_{theirs}.bdg"
        ).read_bytes()


def test_macs3_proba_files_are_checked(tmp_path: Path) -> None:
    def check(text: str, message: str) -> None:
        path = _macs3_write(tmp_path / "bad.csv", text)
        with pytest.raises(ConfigurationError, match=message):
            macs3_review.hmmratac_states_path(path, **_MACS3_INDEXES, verbose=False)

    check("", "empty")
    check("chr1,20,0.1,0.8\n", "three probabilities")
    check("chr1,x,0.1,0.8,0.1\n", "malformed")
    check("chr1,20,0.1,0.8,abc\n", "malformed")
    check("chr1,5,0.1,0.8,0.1\n", "invalid chromosome or bin end")
    check("chr1,20,0.1,-0.8,0.1\n", "invalid probability")
    check("chr1,20,0.1,nan,0.1\n", "invalid probability")
    check("chr1,20,0.1,0.8,0.1\nchr1,25,0.1,0.8,0.1\n", "must not overlap")
    check("chr1,20,0.1,0.8,0.1\nchr2,20,0.1,0.8,0.1\nchr1,40,0.1,0.8,0.1\n", "stay together")
    proba = _macs3_proba(tmp_path)
    prefix = tmp_path / "s"
    for options, message in (
        ({"open_state": 0, "nucleosomal_state": 0, "background_state": 1}, "0, 1 and 2"),
        ({"open_state": 3, "nucleosomal_state": 0, "background_state": 1}, "open_state"),
        ({**_MACS3_INDEXES, "binsize": 0}, "binsize"),
    ):
        with pytest.raises(ConfigurationError, match=message):
            macs3_review.hmmratac_proba_to_bedgraph(proba, prefix, **options, verbose=False)
    (tmp_path / "s_nuc.bdg").write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        macs3_review.hmmratac_proba_to_bedgraph(proba, prefix, **_MACS3_INDEXES, verbose=False)
    with pytest.raises(ConfigurationError, match="not a directory"):
        macs3_review.hmmratac_proba_to_bedgraph(proba, tmp_path, **_MACS3_INDEXES, verbose=False)
    assert not (tmp_path / "s_open.bdg").exists()


# --- hmmratac_states_path ----------------------------------------------------------------------


def test_macs3_states_path_hand_calculation(tmp_path: Path) -> None:
    result = macs3_review.hmmratac_states_path(
        _macs3_proba(tmp_path), **_MACS3_INDEXES, verbose=False
    )
    assert result["states"] == _macs3_state_dicts(_MACS3_STATES)
    assert result["bin_count"] == 13
    assert result["chromosomes"] == ["chr1", "chr2"]
    assert result["state_counts"] == {"open": 3, "nuc": 6, "bg": 6}
    # chr1 is covered over 0-120 and chr2 over 0-70, which is 190 bases in all.
    assert result["state_bases"] == {"open": 40, "nuc": 60, "bg": 90}


def test_macs3_states_path_matches_upstream(tmp_path: Path) -> None:
    command = pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    proba = _macs3_proba(tmp_path)
    result = macs3_review.hmmratac_states_path(proba, **_MACS3_INDEXES, verbose=False)
    with proba.open("rb") as handle:
        expected = command.generate_states_path(handle, 10, 2, 0, 1)
    assert [(s["chrom"], s["start"], s["end"], s["state"]) for s in result["states"]] == [
        (chrom.decode(), start, end, label) for chrom, start, end, label in expected
    ]


def test_macs3_states_path_ties_go_to_open_then_nucleosomal(tmp_path: Path) -> None:
    proba = _macs3_write(
        tmp_path / "ties.csv",
        "chr1,10,0.3,0.3,0.3\nchr1,20,0.4,0.4,0.1\nchr1,30,0.1,0.45,0.45\n",
    )
    # Columns: state 0 nuc, state 1 bg, state 2 open. Row 1 is a three-way tie, row 2 a
    # nuc and bg tie, row 3 a bg and open tie.
    result = macs3_review.hmmratac_states_path(proba, **_MACS3_INDEXES, verbose=False)
    assert [item["state"] for item in result["states"]] == ["open", "nuc", "open"]
    assert [(item["start"], item["end"]) for item in result["states"]] == [
        (0, 10),
        (10, 20),
        (20, 30),
    ]
    merged = macs3_review.hmmratac_states_path(
        _macs3_write(tmp_path / "same.csv", "chr1,10,0.1,0.1,0.8\nchr1,20,0.1,0.1,0.8\n"),
        **_MACS3_INDEXES,
        binsize=10,
        verbose=False,
    )
    assert merged["states"] == [{"chrom": "chr1", "start": 0, "end": 20, "state": "open"}]


# --- write_hmmratac_states_bed -----------------------------------------------------------------


def test_macs3_states_bed_skips_background(tmp_path: Path) -> None:
    command = pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    states = _macs3_state_dicts(_MACS3_STATES)
    output = tmp_path / "states.bed"
    result = macs3_review.write_hmmratac_states_bed(states, output, verbose=False)
    assert output.read_text(encoding="utf-8") == (
        "chr1\t20\t30\tnuc\n"
        "chr1\t30\t50\topen\n"
        "chr1\t50\t60\tnuc\n"
        "chr1\t90\t100\tnuc\n"
        "chr1\t100\t110\topen\n"
        "chr1\t110\t120\tnuc\n"
        "chr2\t40\t50\tnuc\n"
        "chr2\t50\t60\topen\n"
        "chr2\t60\t70\tnuc\n"
    )
    assert result["records"] == 9
    assert result["state_counts"] == {"open": 3, "nuc": 6}
    assert result["skipped_background"] == 6
    with (tmp_path / "up.bed").open("w", encoding="utf-8") as handle:
        command.save_states_bed([(c.encode(), s, e, t) for c, s, e, t in _MACS3_STATES], handle)
    assert (tmp_path / "up.bed").read_text(encoding="utf-8") == output.read_text(encoding="utf-8")
    # A list of four items per state works the same way.
    listed = macs3_review.write_hmmratac_states_bed(
        [list(item) for item in _MACS3_STATES], tmp_path / "listed.bed", verbose=False
    )
    assert listed["records"] == 9


def test_macs3_states_validation(tmp_path: Path) -> None:
    output = tmp_path / "x.bed"
    for states, message in (
        ([{"chrom": "chr1", "start": 0, "end": 10}], "state label"),
        ([{"chrom": "chr1", "start": 0, "end": 10, "state": "closed"}], "state label"),
        ([{"chrom": "chr1", "start": 10, "end": 10, "state": "bg"}], "end after it starts"),
        ([("chr1", 0, 10)], "mapping with chrom"),
        ("chr1", "must be a list"),
        (
            [("chr1", 0, 20, "bg"), ("chr1", 10, 30, "open")],
            "without overlap",
        ),
    ):
        with pytest.raises(ConfigurationError, match=message):
            macs3_review.write_hmmratac_states_bed(states, output, verbose=False)
    assert not output.exists()
    # An empty list is a valid states path and gives an empty file.
    empty = macs3_review.write_hmmratac_states_bed([], output, verbose=False)
    assert empty["records"] == 0
    assert output.read_text(encoding="utf-8") == ""
    with pytest.raises(ConfigurationError, match="already exists"):
        macs3_review.write_hmmratac_states_bed([], output, verbose=False)


# --- write_hmmratac_accessible_regions ---------------------------------------------------------


def _macs3_open_bedgraph(tmp_path: Path) -> Path:
    macs3_review.hmmratac_proba_to_bedgraph(
        _macs3_proba(tmp_path), tmp_path / "s", **_MACS3_INDEXES, verbose=False
    )
    return tmp_path / "s_open.bdg"


def test_macs3_accessible_regions_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    scores = _macs3_open_bedgraph(tmp_path)
    output = tmp_path / "regions.narrowPeak"
    result = macs3_review.write_hmmratac_accessible_regions(
        _macs3_state_dicts(_MACS3_STATES), scores, output, min_length=25, verbose=False
    )
    # Three nuc-open-nuc runs span 40, 30 and 30 bases: chr1 30-50, chr1 100-110 and chr2 50-60.
    # The score is 10 times the open probability 0.8 there and the summit is the segment middle.
    expected = [
        {
            "chrom": "chr1",
            "start": 30,
            "end": 50,
            "name": "MACS_peak_1",
            "score": 8,
            "summit": 40,
        },
        {
            "chrom": "chr1",
            "start": 100,
            "end": 110,
            "name": "MACS_peak_2",
            "score": 8,
            "summit": 105,
        },
        {
            "chrom": "chr2",
            "start": 50,
            "end": 60,
            "name": "MACS_peak_3",
            "score": 8,
            "summit": 55,
        },
    ]
    # MACS3 3.0.5 drops the last region of the list; a release that fixes this writes all three.
    assert result["regions"] in (expected, expected[:2])
    assert result["count"] == len(result["regions"])
    lines = output.read_text(encoding="utf-8").splitlines()
    assert len(lines) == result["count"]
    assert lines[0] == "chr1\t30\t50\tMACS_peak_1\t8\t.\t0\t0\t0\t10"
    assert result["min_length"] == 25


def test_macs3_accessible_regions_min_length_filters_runs(tmp_path: Path) -> None:
    scores = _macs3_open_bedgraph(tmp_path)
    states = _macs3_state_dicts(_MACS3_STATES)
    # Runs span 40, 30 and 30 bases; a run must be strictly longer than min_length.
    none = macs3_review.write_hmmratac_accessible_regions(
        states, scores, tmp_path / "none.np", min_length=40, verbose=False
    )
    assert none["regions"] == []
    assert (tmp_path / "none.np").read_text(encoding="utf-8") == ""
    one = macs3_review.write_hmmratac_accessible_regions(
        states, scores, tmp_path / "one.np", min_length=35, verbose=False
    )
    # Only the 40 base run is long enough, and it is the last region, which MACS3 3.0.5 drops.
    assert one["regions"] == [] or (one["regions"][0]["start"], one["regions"][0]["end"]) == (
        30,
        50,
    )


def test_macs3_accessible_regions_matches_upstream(tmp_path: Path) -> None:
    command = pytest.importorskip("MACS3.Commands.hmmratac_cmd")
    bedgraph = pytest.importorskip("MACS3.IO.BedGraphIO")
    scores = _macs3_open_bedgraph(tmp_path)
    output = tmp_path / "ours.np"
    macs3_review.write_hmmratac_accessible_regions(
        _macs3_state_dicts(_MACS3_STATES), scores, output, min_length=25, verbose=False
    )
    track = bedgraph.bedGraphIO(str(scores)).read_bedGraph(0)
    with (tmp_path / "up.np").open("w", encoding="utf-8") as handle:
        command.save_accessible_regions(
            [(c.encode(), s, e, t) for c, s, e, t in _MACS3_STATES], handle, 25, track
        )
    assert output.read_text(encoding="utf-8") == (tmp_path / "up.np").read_text(encoding="utf-8")


def test_macs3_accessible_regions_use_the_bedgraph_scores(tmp_path: Path) -> None:
    states = [
        ("chr1", 0, 100, "nuc"),
        ("chr1", 100, 300, "open"),
        ("chr1", 300, 400, "nuc"),
        ("chr1", 400, 500, "bg"),
        ("chr1", 500, 600, "nuc"),
        ("chr1", 600, 800, "open"),
        ("chr1", 800, 900, "nuc"),
    ]
    scores = _macs3_write(
        tmp_path / "fold.bdg",
        "track type=bedGraph name=fc\n"
        "chr1\t0\t100\t1.5\n"
        "chr1\t100\t200\t4.0\n"
        "chr1\t200\t300\t12.5\n"
        "chr1\t300\t600\t2.0\n"
        "chr1\t600\t700\t7.0\n"
        "chr1\t700\t900\t3.0\n",
    )
    result = macs3_review.write_hmmratac_accessible_regions(
        _macs3_state_dicts(states), scores, tmp_path / "fc.np", min_length=100, verbose=False
    )
    # The first region scores 10 times its highest value 12.5 with the summit in that segment.
    first = result["regions"][0]
    assert (first["chrom"], first["start"], first["end"]) == ("chr1", 100, 300)
    assert first["score"] == 125
    assert first["summit"] == 250


def test_macs3_accessible_regions_validation(tmp_path: Path) -> None:
    scores = _macs3_open_bedgraph(tmp_path)
    states = _macs3_state_dicts(_MACS3_STATES)

    def call(score_file: Path = scores, **options: Any) -> Any:
        return macs3_review.write_hmmratac_accessible_regions(
            states, score_file, tmp_path / "out.np", verbose=False, **options
        )

    with pytest.raises(ConfigurationError, match="min_length"):
        call(min_length=-1)
    with pytest.raises(ConfigurationError, match="does not exist"):
        call(tmp_path / "missing.bdg")
    with pytest.raises(ConfigurationError, match="four fields"):
        call(_macs3_write(tmp_path / "a.bdg", "chr1\t0\t10\n"))
    with pytest.raises(ConfigurationError, match="malformed"):
        call(_macs3_write(tmp_path / "b.bdg", "chr1\t0\t10\tx\n"))
    with pytest.raises(ConfigurationError, match="invalid interval"):
        call(_macs3_write(tmp_path / "c.bdg", "chr1\t10\t10\t1\n"))
    with pytest.raises(ConfigurationError, match="sorted and contiguous"):
        call(_macs3_write(tmp_path / "d.bdg", "chr1\t0\t10\t1\nchr1\t20\t30\t1\n"))
    with pytest.raises(ConfigurationError, match="sorted and contiguous"):
        call(_macs3_write(tmp_path / "e.bdg", "chr1\t0\t10\t1\nchr2\t0\t10\t1\nchr1\t10\t20\t1\n"))
    with pytest.raises(ConfigurationError, match="no intervals"):
        call(_macs3_write(tmp_path / "f.bdg", "track type=bedGraph\n"))
    assert not (tmp_path / "out.np").exists()
    call(min_length=25)
    with pytest.raises(ConfigurationError, match="already exists"):
        call(min_length=25)


def test_macs3_save_never_leaves_a_model_file_that_cannot_be_read(tmp_path: Path) -> None:
    # Only short and mono fragments: the di and tri tracks are constant, so the covariances are
    # singular. Either the fit is refused, or the file that was written loads again.
    fragments = [item for item in _macs3_fragments() if item[2] - item[1] < 300]
    extracted = macs3_review.extract_hmmratac_training_data(
        fragments, _macs3_regions(), _MACS3_MEANS, _MACS3_STDDEVS, verbose=False
    )
    path = tmp_path / "flat.json"
    try:
        macs3_review.save_hmmratac_model(
            extracted["data"], extracted["lengths"], path, verbose=False
        )
    except ConfigurationError as error:
        assert "positive" in str(error) or "rejected the input" in str(error)
        assert not path.exists()
    else:
        assert macs3_review.load_hmmratac_model(path, verbose=False)["binsize"] == 10


# --- whole chain -------------------------------------------------------------------------------


def test_macs3_chain_from_fragments_to_regions(tmp_path: Path) -> None:
    fragments = _macs3_fragments()
    extracted = macs3_review.extract_hmmratac_training_data(
        fragments, _macs3_regions(), _MACS3_MEANS, _MACS3_STDDEVS, binsize=10, verbose=False
    )
    model_path = tmp_path / "model.json"
    saved = macs3_review.save_hmmratac_model(
        extracted["data"], extracted["lengths"], model_path, binsize=10, verbose=False
    )
    proba = tmp_path / "proba.csv"
    predicted = macs3_review.predict_hmmratac_states(
        model_path,
        extracted["data"],
        extracted["lengths"],
        bins=extracted["bins"],
        output_path=proba,
        verbose=False,
    )
    indexes = {
        "open_state": saved["open_state"],
        "nucleosomal_state": saved["nucleosomal_state"],
        "background_state": saved["background_state"],
    }
    bedgraphs = macs3_review.hmmratac_proba_to_bedgraph(
        proba, tmp_path / "chain", binsize=saved["binsize"], **indexes, verbose=False
    )
    path = macs3_review.hmmratac_states_path(
        proba, binsize=saved["binsize"], **indexes, verbose=False
    )
    assert path["bin_count"] == predicted["rows"] == len(extracted["bins"])
    assert bedgraphs["bin_count"] == predicted["rows"]
    bed = macs3_review.write_hmmratac_states_bed(
        path["states"], tmp_path / "chain.bed", verbose=False
    )
    assert bed["records"] == path["state_counts"]["open"] + path["state_counts"]["nuc"]
    regions = macs3_review.write_hmmratac_accessible_regions(
        path["states"], bedgraphs["files"]["open"], tmp_path / "chain.np", verbose=False
    )
    assert regions["count"] <= path["state_counts"]["open"]
    # Every labelled bin is covered exactly once by the merged states.
    assert sum(path["state_bases"].values()) >= 10 * predicted["rows"]


_DEEPTOOLS_CODES = {"T": 0, "C": 1, "A": 2, "G": 3}

# A 40 bp chromosome of ten 4-mers. Their GC counts are 2, 0, 4, 0, 2, 0, 3, 1, 4, 1.
_DEEPTOOLS_BLOCKS = ["ACGT", "AAAA", "GGCC", "ATAT", "GCAT", "TTTT", "CCCA", "AGAT", "GCGC", "ACTA"]
_DEEPTOOLS_SEQUENCE = "".join(_DEEPTOOLS_BLOCKS)
# Starts (0-based) of 5 bp reads. The first eleven sit on sampled positions (multiples of 4).
_DEEPTOOLS_FORWARD = [0, 4, 4, 8, 16, 24, 24, 24, 28, 32, 32, 2, 5, 17]
_DEEPTOOLS_REVERSE = [6, 20, 21]
_DEEPTOOLS_READ = 5


def _deeptools_write_twobit(path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a little-endian UCSC .2bit file; N is hard-masked, lowercase soft-masked."""

    def runs(flags: list[bool]) -> tuple[list[int], list[int]]:
        starts: list[int] = []
        sizes: list[int] = []
        index = 0
        while index < len(flags):
            if flags[index]:
                stop = index
                while stop < len(flags) and flags[stop]:
                    stop += 1
                starts.append(index)
                sizes.append(stop - index)
                index = stop
            else:
                index += 1
        return starts, sizes

    header = struct.pack("<IIII", 0x1A412743, 0, len(records), 0)
    index_bytes = b""
    body = b""
    offset = len(header) + sum(1 + len(name.encode()) + 4 for name, _ in records)
    for name, sequence in records:
        upper = sequence.upper()
        n_starts, n_sizes = runs([symbol == "N" for symbol in upper])
        mask_starts, mask_sizes = runs([symbol.islower() for symbol in sequence])
        packed = bytearray((len(sequence) + 3) // 4)
        for position, symbol in enumerate(upper):
            packed[position // 4] |= _DEEPTOOLS_CODES.get(symbol, 0) << (6 - 2 * (position % 4))
        record = struct.pack("<I", len(sequence))
        record += struct.pack("<I", len(n_starts)) + struct.pack(f"<{len(n_starts)}I", *n_starts)
        record += struct.pack(f"<{len(n_sizes)}I", *n_sizes)
        record += struct.pack("<I", len(mask_starts))
        record += struct.pack(f"<{len(mask_starts)}I", *mask_starts)
        record += struct.pack(f"<{len(mask_sizes)}I", *mask_sizes)
        record += struct.pack("<I", 0) + bytes(packed)
        raw = name.encode()
        index_bytes += bytes([len(raw)]) + raw + struct.pack("<I", offset)
        body += record
        offset += len(record)
    path.write_bytes(header + index_bytes + body)
    return path


def _deeptools_write_bam(
    directory: Path,
    name: str,
    forward: list[int],
    reverse: list[int],
    *,
    reference: str = "chr1",
    length: int = 40,
) -> Path:
    """Sorted, indexed BAM of 5 bp reads at the given 0-based starts."""

    pysam = pytest.importorskip("pysam")
    rows = [
        f"f{index}\t0\t{reference}\t{start + 1}\t60\t5M\t*\t0\t0\tAAAAA\tIIIII"
        for index, start in enumerate(forward)
    ]
    rows += [
        f"r{index}\t16\t{reference}\t{start + 1}\t60\t5M\t*\t0\t0\tAAAAA\tIIIII"
        for index, start in enumerate(reverse)
    ]
    header = f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{reference}\tLN:{length}\n"
    sam = directory / f"{name}.sam"
    sam.write_text(header + "".join(row + "\n" for row in rows), encoding="utf-8")
    unsorted = directory / f"{name}.unsorted.bam"
    bam = directory / f"{name}.bam"
    pysam.view("-b", "-o", str(unsorted), str(sam), catch_stdout=False)
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam


def _deeptools_gc_fixture(
    directory: Path, sequence: str = _DEEPTOOLS_SEQUENCE
) -> tuple[Path, Path]:
    pytest.importorskip("py2bit")
    pytest.importorskip("deeptools.computeGCBias")
    twobit = _deeptools_write_twobit(directory / "genome.2bit", [("chr1", sequence)])
    bam = _deeptools_write_bam(directory, "reads", _DEEPTOOLS_FORWARD, _DEEPTOOLS_REVERSE)
    return twobit, bam


def _deeptools_expected_gc_table(
    sequence: str, fragment: int, step: int, starts: list[int], max_reads: float = float("inf")
) -> tuple[list[int], list[int]]:
    """Positions and forward read starts per GC count, counted straight from the sequence text."""

    positions = [0] * (fragment + 1)
    reads = [0] * (fragment + 1)
    for position in range(0, len(sequence), step):
        if position + fragment > len(sequence):
            break
        piece = sequence[position : position + fragment]
        if piece.count("N") > 0.05 * fragment:
            continue
        gc = sum(1 for symbol in piece if symbol in "GC")
        positions[gc] += 1
        count = starts.count(position)
        if count >= max_reads:
            continue
        reads[gc] += count
    return positions, reads


def _deeptools_expected_windows(
    sequence: str, size: int, step: int, forward: list[int], reverse: list[int]
) -> list[tuple[int, float]]:
    """(overlapping reads, GC fraction) for windows that start every ``step`` bases."""

    starts = forward + reverse
    windows: list[tuple[int, float]] = []
    for left in range(0, len(sequence), step):
        right = left + size
        if right > len(sequence):
            break
        piece = sequence[left:right]
        if piece.count("N") > 0.05 * size:
            continue
        overlap = sum(1 for s in starts if s < right and s + _DEEPTOOLS_READ > left)
        windows.append((overlap, sum(1 for symbol in piece if symbol in "GC") / size))
    return windows


def _deeptools_write_bigwig(path: Path, sizes: dict[str, int], entries: list[tuple]) -> Path:
    pybigwig = pytest.importorskip("pyBigWig")
    handle = pybigwig.open(str(path), "w")
    handle.addHeader(list(sizes.items()))
    for chrom, start, end, value in entries:
        handle.addEntries([chrom], [start], ends=[end], values=[float(value)])
    handle.close()
    return path


def _deeptools_bigwig_pair(directory: Path) -> tuple[Path, Path]:
    pytest.importorskip("deeptools.writeBedGraph_bam_and_bw")
    sizes = {"chr1": 100, "chr2": 60}
    first = _deeptools_write_bigwig(
        directory / "a.bw",
        sizes,
        [("chr1", 0, 50, 2), ("chr1", 50, 100, 4), ("chr2", 0, 30, 6)],
    )
    second = _deeptools_write_bigwig(
        directory / "b.bw",
        sizes,
        [("chr1", 0, 25, 1), ("chr1", 25, 100, 3), ("chr2", 0, 60, 1)],
    )
    return first, second


def _deeptools_bedgraph_lines(path: str | Path) -> list[str]:
    return Path(path).read_text(encoding="utf-8").splitlines()


def _deeptools_reference_jsd(values: list[float], maxlen: int) -> float:
    """Independent calculation of the synthetic JSD with numpy and scipy only."""

    numpy = pytest.importorskip("numpy")
    poisson = pytest.importorskip("scipy.stats").poisson
    vec = numpy.asarray(values, dtype=float)
    observed = numpy.zeros(maxlen)
    for value in vec:
        value = min(value, maxlen - 1)
        if value > 0:
            observed[int(value)] += 1
    expected = vec.sum() * poisson.pmf(numpy.arange(1, maxlen), vec.mean())

    def lorenz(counts: Any) -> Any:
        index = numpy.arange(len(counts))
        bins = numpy.cumsum(counts) / counts.sum()
        signal = numpy.cumsum(counts * index) / (counts * index).sum()
        grid = numpy.arange(0, 1.00001, 0.00001)
        return numpy.interp(grid, bins, signal, left=0, right=1)

    p = numpy.diff(lorenz(expected))
    q = numpy.diff(lorenz(observed))
    mix = (p + q) / 2
    with numpy.errstate(all="ignore"):
        divergence = 0.5 * numpy.nansum(p * numpy.log2(p / mix))
        divergence += 0.5 * numpy.nansum(q * numpy.log2(q / mix))
    return round(float(numpy.sqrt(divergence)), 15)


# contract ------------------------------------------------------------------------------


def test_deeptools_public_functions_take_keyword_verbose() -> None:
    assert deeptools_review.__all__ == sorted(deeptools_review.__all__)
    for name in deeptools_review.__all__:
        parameter = inspect.signature(getattr(deeptools_review, name)).parameters.get("verbose")
        assert parameter is not None, name
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_deeptools_names_do_not_collide_with_existing_epigenome_modules() -> None:
    existing: set[str] = set()
    for module_name in (
        "MACS3_functions",
        "cooler_functions",
        "cooltools_functions",
        "deeptools_functions",
        "methylpy_functions",
        "pairtools_functions",
    ):
        module = importlib.import_module(f"dnakit.epigenome.{module_name}")
        existing.update(module.__all__)
    assert existing.isdisjoint(deeptools_review.__all__)


# tabulate_gc_fragment_counts -----------------------------------------------------------


def test_deeptools_tabulate_gc_matches_hand_count(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    result = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=10, verbose=False
    )
    json.dumps(result)
    assert result["status"] == "completed"
    assert result["operation"] == "tabulate_gc_fragment_counts"
    assert result["provenance"]["package"].startswith("deeptools ")
    # Sampled positions 0, 4, ..., 36 hold the 4-mers listed above, so the positions per GC
    # count 0..4 are 3, 2, 2, 1, 2. Forward read starts on those positions: 2 at GC 0
    # (position 4 twice), 1 at GC 1 (position 28), 1 + 1 at GC 2 (positions 0 and 16), 3 at GC 3
    # (position 24) and 1 + 2 at GC 4 (positions 8 and 32). Reverse reads and starts off the
    # grid (2, 5, 17) are not counted.
    table = result["gc_table"]
    assert [row["gc_bases"] for row in table] == [0, 1, 2, 3, 4]
    assert [row["genome_positions"] for row in table] == [3, 2, 2, 1, 2]
    assert [row["read_starts"] for row in table] == [2, 1, 2, 3, 3]
    # The ratio is (reads / positions) over the overall density 11 / 10.
    scaling = 10 / 11
    expected = [
        reads / positions * scaling
        for reads, positions in zip([2, 1, 2, 3, 3], [3, 2, 2, 1, 2], strict=True)
    ]
    assert [row["read_ratio"] for row in table] == pytest.approx(expected)
    assert result["positions_sampled"] == 10
    assert result["read_starts_counted"] == 11
    assert result["step_size"] == 4
    assert result["mapped_reads"] == 17
    assert result["reads_per_bp"] == pytest.approx(17 / 40)
    positions, reads = _deeptools_expected_gc_table(_DEEPTOOLS_SEQUENCE, 4, 4, _DEEPTOOLS_FORWARD)
    assert [row["genome_positions"] for row in table] == positions
    assert [row["read_starts"] for row in table] == reads


def test_deeptools_tabulate_gc_skips_peak_positions(tmp_path: Path) -> None:
    stats = pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    # reads_per_bp = 17 / 400 and the mean in a 4 bp fragment is 4 * 0.0425 * 4 = 0.68.
    result = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 400, sample_size=10, verbose=False
    )
    assert result["max_reads"] == pytest.approx(float(stats.poisson(0.68).isf(0.1)))
    assert result["max_reads"] == 2.0
    # Positions 4, 24 and 32 carry two or more read starts and are skipped as peaks, but they
    # are still counted as sampled positions.
    table = result["gc_table"]
    assert [row["genome_positions"] for row in table] == [3, 2, 2, 1, 2]
    assert [row["read_starts"] for row in table] == [0, 1, 2, 0, 1]
    positions, reads = _deeptools_expected_gc_table(
        _DEEPTOOLS_SEQUENCE, 4, 4, _DEEPTOOLS_FORWARD, max_reads=2
    )
    assert [row["read_starts"] for row in table] == reads
    assert [row["genome_positions"] for row in table] == positions


def test_deeptools_tabulate_gc_leaves_out_n_fragments(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    sequence = "".join(
        "NNNN" if index == 3 else block for index, block in enumerate(_DEEPTOOLS_BLOCKS)
    )
    twobit, bam = _deeptools_gc_fixture(tmp_path, sequence)
    result = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=10, verbose=False
    )
    # The N block at positions 12 to 15 is the 4-mer with GC 0 and no read starts.
    assert result["positions_sampled"] == 9
    assert [row["genome_positions"] for row in result["gc_table"]] == [2, 2, 2, 1, 2]
    assert [row["read_starts"] for row in result["gc_table"]] == [2, 1, 2, 3, 3]


def test_deeptools_tabulate_gc_region_and_default_step(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    region = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=10, region="chr1:0:20", verbose=False
    )
    # Positions 0, 4, 8, 12, 16 with GC 2, 0, 4, 0, 2 and read starts 1, 2, 1, 0, 1.
    assert region["region"] == "chr1:0:20"
    assert [row["genome_positions"] for row in region["gc_table"]] == [2, 0, 2, 0, 1]
    assert [row["read_starts"] for row in region["gc_table"]] == [2, 0, 2, 0, 1]
    default = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=20, verbose=False
    )
    # The default step is max(int(genome size / sample size), 1) = 2.
    assert default["step_size"] == 2
    assert default["positions_sampled"] == 19
    explicit = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=20, step_size=8, verbose=False
    )
    assert explicit["step_size"] == 8
    assert explicit["positions_sampled"] == 5


def test_deeptools_tabulate_gc_matches_bam_without_chr_prefix(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    plain = _deeptools_write_bam(
        tmp_path, "plain", _DEEPTOOLS_FORWARD, _DEEPTOOLS_REVERSE, reference="1"
    )
    named = deeptools_review.tabulate_gc_fragment_counts(
        twobit, bam, 4, 40, sample_size=10, verbose=False
    )
    mapped = deeptools_review.tabulate_gc_fragment_counts(
        twobit, plain, 4, 40, sample_size=10, verbose=False
    )
    # computeGCBias itself fails for this pair; the wrapper maps 1 to chr1.
    assert mapped["gc_table"] == named["gc_table"]


def test_deeptools_gc_functions_restore_module_state(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("scipy.stats")
    numpy = pytest.importorskip("numpy")
    module = importlib.import_module("deeptools.computeGCBias")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    monkeypatch.delattr(module, "global_vars", raising=False)
    saved = numpy.geterr()
    try:
        numpy.seterr(all="warn")
        deeptools_review.tabulate_gc_fragment_counts(
            twobit, bam, 4, 40, sample_size=10, verbose=False
        )
        deeptools_review.count_reads_per_gc_bins(
            twobit, bam, 40, region_size=8, sample_size=40, step_size=4, verbose=False
        )
        # A failing upstream call (fragments longer than the chromosome) restores the state too.
        with pytest.raises(ConfigurationError, match="No fragments"):
            deeptools_review.tabulate_gc_fragment_counts(
                twobit, bam, 400, 40, sample_size=10, verbose=False
            )
        assert not hasattr(module, "global_vars")
        assert numpy.geterr()["divide"] == "warn"
        sentinel = {"marker": 1}
        monkeypatch.setattr(module, "global_vars", sentinel, raising=False)
        deeptools_review.count_reads_per_gc_bins(
            twobit, bam, 40, region_size=8, sample_size=40, step_size=4, verbose=False
        )
        assert module.global_vars is sentinel
        assert sentinel == {"marker": 1}
    finally:
        numpy.seterr(**saved)


def test_deeptools_wrappers_keep_numpy_error_settings_after_first_import(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    # deepTools sets numpy.seterr(all="ignore") when its modules are imported; in a fresh
    # interpreter the first call must leave the caller's settings as they were.
    code = (
        "import sys, numpy\n"
        "numpy.seterr(all='warn')\n"
        "from dnakit.epigenome import deeptools_review_functions as d\n"
        "d.tabulate_gc_fragment_counts(sys.argv[1], sys.argv[2], 4, 40, sample_size=10)\n"
        "import deeptools.computeGCBias as c\n"
        "print(numpy.geterr()['divide'], hasattr(c, 'global_vars'))\n"
    )
    done = subprocess.run(
        [sys.executable, "-c", code, str(twobit), str(bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode == 0, done.stderr[-500:]
    assert done.stdout.split() == ["warn", "False"]


def test_deeptools_tabulate_gc_rejects_bad_input(tmp_path: Path) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    call = deeptools_review.tabulate_gc_fragment_counts
    for bad in (0, -4, True, 4.0, "4", None):
        with pytest.raises(ConfigurationError, match="fragment_length"):
            call(twobit, bam, bad, 40, verbose=False)
    for bad in (0, True, 1.5):
        with pytest.raises(ConfigurationError, match="effective_genome_size"):
            call(twobit, bam, 4, bad, verbose=False)
    with pytest.raises(ConfigurationError, match="sample_size"):
        call(twobit, bam, 4, 40, sample_size=1, verbose=False)
    with pytest.raises(ConfigurationError, match="step_size"):
        call(twobit, bam, 4, 40, step_size=0, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        call(tmp_path / "missing.2bit", bam, 4, 40, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        call(twobit, tmp_path / "missing.bam", 4, 40, verbose=False)
    with pytest.raises(ConfigurationError, match="region"):
        call(twobit, bam, 4, 40, region="chr1", verbose=False)
    with pytest.raises(ConfigurationError, match="region"):
        call(twobit, bam, 4, 40, region="chr1:20:10", verbose=False)
    with pytest.raises(ConfigurationError, match="Unknown chromosome"):
        call(twobit, bam, 4, 40, region="chrZ:0:10", verbose=False)
    not_twobit = tmp_path / "not.2bit"
    not_twobit.write_text("text", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="2bit"):
        call(not_twobit, bam, 4, 40, verbose=False)
    unindexed = tmp_path / "unindexed.bam"
    shutil.copy(bam, unindexed)
    with pytest.raises(ConfigurationError, match="index"):
        call(twobit, unindexed, 4, 40, verbose=False)
    empty = _deeptools_write_bam(tmp_path, "empty", [], [])
    with pytest.raises(ConfigurationError, match="no mapped reads"):
        call(twobit, empty, 4, 40, verbose=False)
    other = _deeptools_write_twobit(tmp_path / "other.2bit", [("chrX", "ACGT" * 10)])
    with pytest.raises(ConfigurationError, match="cannot be matched"):
        call(other, bam, 4, 40, verbose=False)


def test_deeptools_tabulate_gc_limits_positions(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("scipy.stats")
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    monkeypatch.setattr(deeptools_review, "_MAX_POSITIONS", 5)
    with pytest.raises(ConfigurationError, match="limit"):
        deeptools_review.tabulate_gc_fragment_counts(
            twobit, bam, 4, 40, sample_size=40, verbose=False
        )


# count_reads_per_gc_bins ---------------------------------------------------------------


def test_deeptools_count_reads_per_gc_matches_hand_count(tmp_path: Path) -> None:
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    result = deeptools_review.count_reads_per_gc_bins(
        twobit, bam, 40, region_size=8, sample_size=40, step_size=4, verbose=False
    )
    json.dumps(result)
    assert result["operation"] == "count_reads_per_gc_bins"
    assert result["window_count"] == 9
    # Window [0, 8) is ACGTAAAA: GC 2 / 8 = 0.25, and the reads starting at 0, 2, 4, 4, 5 (forward)
    # and 6 (reverse) overlap it, so 6 reads. Window [8, 16) is GGCCATAT: 0.5, with reads
    # starting at 4, 4, 5, 6 (they end after 8), 8 and 12 -> 5 plus none after.
    assert result["windows"][0] == {"reads": 6, "gc_fraction": 0.25}
    assert result["windows"][2] == {"reads": 5, "gc_fraction": 0.5}
    expected = _deeptools_expected_windows(
        _DEEPTOOLS_SEQUENCE, 8, 4, _DEEPTOOLS_FORWARD, _DEEPTOOLS_REVERSE
    )
    assert [(w["reads"], w["gc_fraction"]) for w in result["windows"]] == expected
    assert result["mean_reads"] == pytest.approx(sum(reads for reads, _ in expected) / 9)


def test_deeptools_count_reads_per_gc_default_step_region_and_n(tmp_path: Path) -> None:
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    default = deeptools_review.count_reads_per_gc_bins(
        twobit, bam, 40, region_size=8, sample_size=4, verbose=False
    )
    # The default step is ten times int(40 / 4) = 100, which leaves the window at 0 only.
    assert default["step_size"] == 100
    assert default["windows"] == [{"reads": 6, "gc_fraction": 0.25}]
    region = deeptools_review.count_reads_per_gc_bins(
        twobit, bam, 40, region_size=8, step_size=4, region="chr1:0:20", verbose=False
    )
    assert region["window_count"] == 5
    assert region["region"] == "chr1:0:20"
    # A window longer than the chromosome gives no window at all.
    long = deeptools_review.count_reads_per_gc_bins(
        twobit, bam, 40, region_size=100, step_size=4, verbose=False
    )
    assert long["windows"] == []
    assert long["mean_reads"] is None


def test_deeptools_count_reads_per_gc_leaves_out_n_windows(tmp_path: Path) -> None:
    pytest.importorskip("py2bit")
    pytest.importorskip("deeptools.computeGCBias")
    sequence = "".join(
        "NNNN" if index == 3 else block for index, block in enumerate(_DEEPTOOLS_BLOCKS)
    )
    twobit, bam = _deeptools_gc_fixture(tmp_path, sequence)
    result = deeptools_review.count_reads_per_gc_bins(
        twobit, bam, 40, region_size=8, step_size=4, verbose=False
    )
    # Windows starting at 8 and 12 contain four N of eight bases and are left out.
    assert result["window_count"] == 7
    expected = _deeptools_expected_windows(sequence, 8, 4, _DEEPTOOLS_FORWARD, _DEEPTOOLS_REVERSE)
    assert len(expected) == 7
    assert [(w["reads"], w["gc_fraction"]) for w in result["windows"]] == expected


def test_deeptools_count_reads_per_gc_rejects_bad_input(tmp_path: Path, monkeypatch) -> None:
    twobit, bam = _deeptools_gc_fixture(tmp_path)
    call = deeptools_review.count_reads_per_gc_bins
    for bad in (0, True, 8.0, "8"):
        with pytest.raises(ConfigurationError, match="region_size"):
            call(twobit, bam, 40, region_size=bad, verbose=False)
    with pytest.raises(ConfigurationError, match="effective_genome_size"):
        call(twobit, bam, 0, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        call(tmp_path / "none.2bit", bam, 40, verbose=False)
    monkeypatch.setattr(deeptools_review, "_MAX_WINDOWS", 5)
    with pytest.raises(ConfigurationError, match="windows"):
        call(twobit, bam, 40, region_size=8, step_size=4, verbose=False)


# fingerprint_synthetic_jsd -------------------------------------------------------------


def test_deeptools_synthetic_jsd_matches_independent_calculation(monkeypatch) -> None:
    plot = pytest.importorskip("deeptools.plotFingerprint")
    # A shorter coverage axis gives the same value for these vectors and keeps the test fast.
    monkeypatch.setattr(plot, "MAXLEN", 200)
    vectors = [
        [0, 1, 1, 2, 3, 5, 8, 2, 1, 0],
        [4] * 20,
        [0, 0, 0, 10, 10, 10, 0, 0],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [3.7, 2.2, 0.4, 5.9],
    ]
    for vector in vectors:
        result = deeptools_review.fingerprint_synthetic_jsd(vector, verbose=False)
        json.dumps(result)
        assert result["jsd"] == pytest.approx(_deeptools_reference_jsd(vector, 200), abs=1e-12)
        assert 0 < result["jsd"] < 1
        assert result["bins"] == len(vector)
        assert result["mean_coverage"] == pytest.approx(sum(vector) / len(vector))
    # The wrapper is deterministic.
    first = deeptools_review.fingerprint_synthetic_jsd(vectors[0], verbose=False)
    second = deeptools_review.fingerprint_synthetic_jsd(vectors[0], verbose=False)
    assert first == second
    assert first["provenance"]["package"].startswith("deeptools ")


def test_deeptools_synthetic_jsd_full_length_values() -> None:
    pytest.importorskip("deeptools.plotFingerprint")
    # Values of the unmodified upstream call (coverage axis of 10 million), which also agree
    # with the independent calculation above.
    result = deeptools_review.fingerprint_synthetic_jsd(
        [0, 1, 1, 2, 3, 5, 8, 2, 1, 0], verbose=False
    )
    assert result["jsd"] == pytest.approx(0.242168914850192, abs=1e-12)
    assert result["bins_above_maximum"] == 0
    uniform = deeptools_review.fingerprint_synthetic_jsd([4] * 20, verbose=False)
    assert uniform["jsd"] == pytest.approx(0.299153600431795, abs=1e-12)


def test_deeptools_synthetic_jsd_counts_bins_above_maximum(monkeypatch) -> None:
    plot = pytest.importorskip("deeptools.plotFingerprint")
    monkeypatch.setattr(plot, "MAXLEN", 50)
    result = deeptools_review.fingerprint_synthetic_jsd([60, 1, 3, 50, 2], verbose=False)
    assert result["bins_above_maximum"] == 2
    assert 0 < result["jsd"] < 1


def test_deeptools_synthetic_jsd_rejects_bad_input() -> None:
    call = deeptools_review.fingerprint_synthetic_jsd
    for bad in ([], "abc", {"a": 1}, None, 5):
        with pytest.raises(ConfigurationError, match="coverage"):
            call(bad, verbose=False)
    for bad in ([1, -1], [1, True], [1, "2"], [1, float("nan")], [1, float("inf")], [2, None]):
        with pytest.raises(ConfigurationError, match="coverage entry"):
            call(bad, verbose=False)
    with pytest.raises(ConfigurationError, match="1 or more"):
        call([0, 0.5, 0.2], verbose=False)


# bigwig_bedgraph_tiles -----------------------------------------------------------------


def test_deeptools_bigwig_tiles_mean_matches_hand_calculation(tmp_path: Path) -> None:
    first, second = _deeptools_bigwig_pair(tmp_path)
    output = tmp_path / "mean.bdg"
    result = deeptools_review.bigwig_bedgraph_tiles(
        [first, second], output, tile_size=25, verbose=False
    )
    json.dumps(result)
    assert result["operation"] == "bigwig_bedgraph_tiles"
    # chr1 tiles of 25: a = 2, 2, 4, 4 and b = 1, 3, 3, 3 -> means 1.5, 2.5, 3.5, 3.5 (the last
    # two are merged). chr2: a = 6 on [0, 30), b = 1 everywhere: tile [0, 25) is (6 + 1) / 2,
    # tile [25, 50) is (1.2 + 1) / 2 because only 5 of its 25 bases are covered in a, and tile
    # [50, 60) is (0 + 1) / 2.
    assert _deeptools_bedgraph_lines(output) == [
        "chr1\t0\t25\t1.5",
        "chr1\t25\t50\t2.5",
        "chr1\t50\t100\t3.5",
        "chr2\t0\t25\t3.5",
        "chr2\t25\t50\t1.1",
        "chr2\t50\t60\t0.5",
    ]
    assert result["intervals"] == 6
    assert result["chromosomes"] == ["chr1", "chr2"]
    assert result["first_intervals"][0] == {"chrom": "chr1", "start": 0, "end": 25, "value": 1.5}
    assert result["scale_factors"] == [1.0, 1.0]
    assert result["output"] == str(output.resolve())


def test_deeptools_bigwig_tiles_sum_with_scale_factors(tmp_path: Path) -> None:
    first, second = _deeptools_bigwig_pair(tmp_path)
    output = tmp_path / "difference.bdg"
    result = deeptools_review.bigwig_bedgraph_tiles(
        [first, second],
        output,
        tile_size=25,
        reduction="sum",
        scale_factors=[1, -1],
        verbose=False,
    )
    # a - b per tile: chr1 1, -1, 1, 1 and chr2 5, 0.2, -1.
    assert _deeptools_bedgraph_lines(output) == [
        "chr1\t0\t25\t1",
        "chr1\t25\t50\t-1",
        "chr1\t50\t100\t1",
        "chr2\t0\t25\t5",
        "chr2\t25\t50\t0.2",
        "chr2\t50\t60\t-1",
    ]
    assert result["reduction"] == "sum"
    assert result["scale_factors"] == [1.0, -1.0]
    single = tmp_path / "doubled.bdg"
    deeptools_review.bigwig_bedgraph_tiles(
        [first], single, tile_size=50, scale_factors=[2], verbose=False
    )
    # Tiles of 50: chr1 2 and 4 -> 4 and 8; chr2 tile [0, 50) is 6 * 30 / 50 = 3.6 -> 7.2 and
    # the last 10 bases are empty (zero, which upstream does not write at the end of a chromosome).
    assert _deeptools_bedgraph_lines(single) == [
        "chr1\t0\t50\t4",
        "chr1\t50\t100\t8",
        "chr2\t0\t50\t7.2",
    ]


def test_deeptools_bigwig_tiles_region_missing_data_and_fixed_step(tmp_path: Path) -> None:
    first, _second = _deeptools_bigwig_pair(tmp_path)
    region_output = tmp_path / "region.bdg"
    region = deeptools_review.bigwig_bedgraph_tiles(
        [first], region_output, tile_size=25, region="chr1:10:90", verbose=False
    )
    # The region start moves down to the tile boundary 0 and the end stays at 90.
    assert _deeptools_bedgraph_lines(region_output) == ["chr1\t0\t50\t2", "chr1\t50\t90\t4"]
    assert region["region"] == "chr1:10:90"
    nan_output = tmp_path / "nan.bdg"
    deeptools_review.bigwig_bedgraph_tiles(
        [first], nan_output, tile_size=25, missing_data_as_zero=False, verbose=False
    )
    # Tiles that include uncovered bases are NaN and are not written.
    assert _deeptools_bedgraph_lines(nan_output) == [
        "chr1\t0\t50\t2",
        "chr1\t50\t100\t4",
        "chr2\t0\t25\t6",
    ]
    fixed_output = tmp_path / "fixed.bdg"
    fixed = deeptools_review.bigwig_bedgraph_tiles(
        [first],
        fixed_output,
        tile_size=25,
        fixed_step=True,
        missing_data_as_zero=False,
        verbose=False,
    )
    assert _deeptools_bedgraph_lines(fixed_output) == [
        "chr1\t0\t25\t2",
        "chr1\t25\t50\t2",
        "chr1\t50\t75\t4",
        "chr1\t75\t100\t4",
        "chr2\t0\t25\t6",
        "chr2\t25\t50\tnan",
        "chr2\t50\t60\tnan",
    ]
    assert fixed["intervals"] == 7
    assert fixed["first_intervals"][-1]["value"] is None


def test_deeptools_bigwig_tiles_rejects_bad_input(tmp_path: Path) -> None:
    first, second = _deeptools_bigwig_pair(tmp_path)
    call = deeptools_review.bigwig_bedgraph_tiles
    target = tmp_path / "out.bdg"
    for bad in ([], "a.bw", None):
        with pytest.raises(ConfigurationError, match="bigwig_paths"):
            call(bad, target, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        call([tmp_path / "missing.bw"], target, verbose=False)
    existing = tmp_path / "existing.bdg"
    existing.write_text("old\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        call([first], existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "old\n"
    with pytest.raises(ConfigurationError, match="directory"):
        call([first], tmp_path / "no_dir" / "out.bdg", verbose=False)
    for bad in (0, -5, True, 2.5, "25"):
        with pytest.raises(ConfigurationError, match="tile_size"):
            call([first], target, tile_size=bad, verbose=False)
    with pytest.raises(ConfigurationError, match="reduction"):
        call([first], target, reduction="median", verbose=False)
    with pytest.raises(ConfigurationError, match="one value per"):
        call([first, second], target, scale_factors=[1], verbose=False)
    with pytest.raises(ConfigurationError, match="scale_factors entry"):
        call([first], target, scale_factors=[float("nan")], verbose=False)
    for bad in ("chr1", "chr1:5", "chr1:50:10", 5):
        with pytest.raises(ConfigurationError, match="region"):
            call([first], target, region=bad, verbose=False)
    with pytest.raises(ConfigurationError, match="Unknown chromosome"):
        call([first], target, region="chrZ:0:10", verbose=False)
    with pytest.raises(ConfigurationError, match="missing_data_as_zero"):
        call([first], target, missing_data_as_zero="yes", verbose=False)
    with pytest.raises(ConfigurationError, match="fixed_step"):
        call([first], target, fixed_step=1, verbose=False)
    # Neither a BAM file (the BAM branch of upstream writeBedGraph fails) nor text is a bigWig.
    pytest.importorskip("pysam")
    bam = _deeptools_write_bam(tmp_path, "reads", [0], [])
    with pytest.raises(ConfigurationError, match="not a bigWig"):
        call([bam], target, verbose=False)
    text = tmp_path / "text.bw"
    text.write_text("chr1 0 10 1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="not a bigWig"):
        call([text], target, verbose=False)
    assert not target.exists()


def test_deeptools_bigwig_tiles_rejects_incompatible_files(tmp_path: Path) -> None:
    first, _second = _deeptools_bigwig_pair(tmp_path)
    shorter = _deeptools_write_bigwig(
        tmp_path / "shorter.bw", {"chr1": 80, "chr2": 60}, [("chr1", 0, 80, 1)]
    )
    other = _deeptools_write_bigwig(tmp_path / "other.bw", {"chrX": 50}, [("chrX", 0, 50, 1)])
    call = deeptools_review.bigwig_bedgraph_tiles
    with pytest.raises(ConfigurationError, match="different lengths"):
        call([first, shorter], tmp_path / "x.bdg", verbose=False)
    with pytest.raises(ConfigurationError, match="share no chromosome"):
        call([first, other], tmp_path / "y.bdg", verbose=False)
    assert not (tmp_path / "x.bdg").exists()
    assert not (tmp_path / "y.bdg").exists()
