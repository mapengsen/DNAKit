"""Tests for dnakit.variants.bgen_reader_gap_functions (bgen_reader example file fetch).

Only the refusal of unknown names runs by default; the download runs with DNAKIT_LIVE_TESTS=1.
"""

from __future__ import annotations

import inspect
import os

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("bgen_reader")

from dnakit.variants import bgen_reader_gap_functions as facades  # noqa: E402


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_unknown_example_name_is_refused_before_any_download():
    with pytest.raises(ConfigurationError) as info:
        facades.fetch_bgen_example("not-a-sample.bgen", verbose=False)
    assert "example.bgen" in str(info.value)


@pytest.mark.skipif(os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="network test; set DNAKIT_LIVE_TESTS=1")
def test_live_fetch_example_bgen():
    pytest.importorskip("requests")
    result = facades.fetch_bgen_example("example.bgen", verbose=False)
    assert result["size_bytes"] > 0
    assert result["path"].endswith("example.bgen")
