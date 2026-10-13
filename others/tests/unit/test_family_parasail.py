"""parasail kernels chosen by name parts and SSW with a reusable query profile."""

from __future__ import annotations

import inspect
import itertools

import pytest

pytest.importorskip("parasail")
import parasail

from dnakit.alignment.parasail_family_functions import (
    _ALGORITHMS,
    _function_name,
    parasail_dp_variant,
    parasail_ssw_profile_align,
)
from dnakit.alignment.parasail_functions import parasail_ssw_align
from dnakit.exceptions import ConfigurationError

QUERY = "ACGTTACGGATC"
TARGET = "ACGATACGATCC"
MATRIX = parasail.matrix_create("ACGTN", 5, -4)


def _direct(name: str, profile_width: str | None) -> object:
    function = getattr(parasail, name)
    if profile_width is None:
        return function(QUERY, TARGET, 10, 1, MATRIX)
    profile = getattr(parasail, f"profile_create_{profile_width}")(QUERY, MATRIX)
    return function(profile, TARGET, 10, 1)


@pytest.mark.parametrize("algorithm", _ALGORITHMS)
@pytest.mark.parametrize("output", ["score", "stats", "table", "rowcol", "trace"])
def test_every_algorithm_and_output_matches_the_direct_call(algorithm: str, output: str) -> None:
    result = parasail_dp_variant(
        QUERY,
        TARGET,
        algorithm=algorithm,
        output=output,
        method="striped",
        width="sat",
        verbose=False,
    )
    direct = _direct(result["parasail_function"], None)
    assert result["parasail_function"] == _function_name(algorithm, output, "striped", False, "sat")
    assert result["score"] == int(direct.score)
    assert (result["end_query"], result["end_target"]) == (
        int(direct.end_query),
        int(direct.end_ref),
    )


@pytest.mark.parametrize(
    ("method", "width", "profile"),
    [
        *itertools.product(["scan", "striped", "diag"], ["8", "16", "32", "64", "sat"], [False]),
        *itertools.product(["scan", "striped"], ["8", "16", "32", "64", "sat"], [True]),
        ("reference", None, False),
        ("scan", None, False),
    ],
)
def test_method_width_and_profile_reach_the_named_kernel(
    method: str, width: str | None, profile: bool
) -> None:
    result = parasail_dp_variant(
        QUERY,
        TARGET,
        algorithm="sg_qe_de",
        output="rowcol",
        method=method,
        width=width,
        profile=profile,
        verbose=False,
    )
    name = result["parasail_function"]
    assert name.startswith("sg_qe_de_rowcol")
    assert name.endswith("_profile" + (f"_{width}" if width else "")) if profile else True
    direct = _direct(name, width if profile else None)
    assert result["score_row"] == [int(value) for value in direct.score_row]
    assert result["score_col"] == [int(value) for value in direct.score_col]
    assert result["profile"] is profile


def test_stats_equal_hand_count_for_identical_sequences() -> None:
    result = parasail_dp_variant(
        "ACGTACGT", "ACGTACGT", algorithm="nw", output="stats", method="reference", verbose=False
    )
    assert result["score"] == 8 * 5
    assert (result["matches"], result["similar"], result["length"]) == (8, 8, 8)
    assert result["parasail_function"] == "nw_stats"
    assert result["width"] is None


def test_table_has_one_row_per_query_base_and_corner_equals_global_score() -> None:
    result = parasail_dp_variant(
        "ACGT", "ACGA", algorithm="nw", output="table", method="scan", width="32", verbose=False
    )
    table = result["score_table"]
    assert len(table) == 4 and all(len(row) == 4 for row in table)
    assert table[-1][-1] == result["score"] == 3 * 5 - 4


def test_stats_table_and_rowcol_expose_every_stat_grid() -> None:
    both = parasail_dp_variant(
        "ACGT",
        "ACGA",
        algorithm="sw",
        output="stats_table",
        method="diag",
        width="16",
        verbose=False,
    )
    assert {"score_table", "matches_table", "similar_table", "length_table"} <= set(both)
    cells = parasail_dp_variant(
        "ACGT",
        "ACGA",
        algorithm="sw",
        output="stats_rowcol",
        method="diag",
        width="16",
        verbose=False,
    )
    assert {"matches_row", "matches_col", "similar_row", "length_col"} <= set(cells)
    assert cells["score_row"] == both["score_table"][-1]


def test_trace_output_returns_aligned_strings_and_cigar() -> None:
    result = parasail_dp_variant(
        QUERY, TARGET, algorithm="nw", output="trace", method="scan", width="sat", verbose=False
    )
    assert len(result["query_aligned"]) == len(result["target_aligned"])
    assert result["query_aligned"].replace("-", "") == QUERY
    assert result["target_aligned"].replace("-", "") == TARGET
    assert result["cigar"]


def test_results_are_json_serializable_and_repeatable() -> None:
    import json

    first = parasail_dp_variant(QUERY, TARGET, output="stats_table", verbose=False)
    second = parasail_dp_variant(QUERY, TARGET, output="stats_table", verbose=False)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"algorithm": "xx"},
        {"output": "everything"},
        {"method": "fast"},
        {"width": 12},
        {"width": "sat", "profile": True, "method": "reference"},
        {"profile": True, "output": "stats"},
        {"profile": True, "output": "stats_rowcol", "method": "scan"},
        {"profile": True, "width": None},
        {"profile": "yes"},
        {"max_cells": 5},
        {"gap_open": -1},
    ],
)
def test_bad_arguments_are_refused(kwargs: dict[str, object]) -> None:
    with pytest.raises(ConfigurationError):
        parasail_dp_variant(QUERY, TARGET, verbose=False, **kwargs)


def test_non_dna_and_oversized_table_are_refused() -> None:
    with pytest.raises(ConfigurationError):
        parasail_dp_variant("ACGU", TARGET, verbose=False)
    long = "A" * 600
    with pytest.raises(ConfigurationError, match="table"):
        parasail_dp_variant(long, long, output="table", verbose=False)


def test_narrow_width_saturation_is_reported() -> None:
    from dnakit.exceptions import BackendExecutionError

    long = "ACGT" * 60
    with pytest.raises(BackendExecutionError, match="saturated"):
        parasail_dp_variant(long, long, method="striped", width="8", verbose=False)


def test_ledger_names_rebuild_exactly() -> None:
    names = {
        "sg_qb_db_stats_rowcol_striped_profile_8": (
            "sg_qb_db",
            "stats_rowcol",
            "striped",
            True,
            "8",
        ),
        "nw_trace_scan_sat": ("nw", "trace", "scan", False, "sat"),
        "sw_table_diag_64": ("sw", "table", "diag", False, "64"),
        "nw_stats": ("nw", "stats", "reference", False, "sat"),
        "sg_rowcol_scan": ("sg", "rowcol", "scan", False, ""),
    }
    for name, parts in names.items():
        assert _function_name(*parts) == name
        assert callable(getattr(parasail, name))


# --- parasail_ssw_profile_align ---------------------------------------------------------------


def test_profile_align_matches_one_shot_ssw_for_each_target() -> None:
    targets = [TARGET, "TTTTTTTT", "GGACGTTACGGATCAA"]
    result = parasail_ssw_profile_align(
        QUERY, targets, match=2, mismatch=-3, gap_open=5, gap_extend=2, verbose=False
    )
    assert result["target_count"] == 3
    for hit, target in zip(result["hits"], targets, strict=True):
        single = parasail_ssw_align(
            QUERY, target, match=2, mismatch=-3, gap_open=5, gap_extend=2, verbose=False
        )
        for key in (
            "score",
            "found",
            "query_start",
            "query_end",
            "target_start",
            "target_end",
            "cigar",
            "matches",
            "mismatches",
            "identity",
        ):
            assert hit[key] == single[key], key


@pytest.mark.parametrize("score_size", [0, 1, 2])
def test_profile_align_accepts_every_supported_score_size(score_size: int) -> None:
    result = parasail_ssw_profile_align(QUERY, [TARGET], score_size=score_size, verbose=False)
    assert result["hits"][0]["score"] > 0 and result["score_size"] == score_size


@pytest.mark.parametrize("score_size", [-1, 3, 16, True, "1"])
def test_profile_align_refuses_score_sizes_that_crash_parasail(score_size: object) -> None:
    with pytest.raises(ConfigurationError):
        parasail_ssw_profile_align(QUERY, [TARGET], score_size=score_size, verbose=False)


@pytest.mark.parametrize("targets", [TARGET, [], [TARGET, "ACGX"], b"ACGT"])
def test_profile_align_refuses_bad_targets(targets: object) -> None:
    with pytest.raises(ConfigurationError):
        parasail_ssw_profile_align(QUERY, targets, verbose=False)


# --- exports ----------------------------------------------------------------------------------


def test_exported_from_alignment_and_signatures_are_agent_friendly() -> None:
    import dnakit.alignment as alignment
    import dnakit.alignment.parasail_family_functions as module

    for name in module.__all__:
        assert name in alignment.__all__
        assert inspect.unwrap(getattr(alignment, name)) is inspect.unwrap(getattr(module, name))
        assert inspect.signature(getattr(module, name)).parameters["verbose"].default is True
