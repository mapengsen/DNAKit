"""Upstream checks behind the G18c review decisions for sbol3, screed and pyfaidx.

The review wraps no new function in this batch. Each test calls the real upstream API, or
works a hand calculation, to confirm why a reviewed row is not a DNA-level wrapper.
"""

from __future__ import annotations

import os
from io import BytesIO
from pathlib import Path

import pytest

from dnakit.molecular_lab import sbol3_review_functions
from dnakit.sequence import pyfaidx_review_functions, screed_review_functions
from dnakit.sequence.screed_functions import format_fastx_records


def test_review_modules_export_no_new_wrappers() -> None:
    assert list(sbol3_review_functions.__all__) == []
    assert list(screed_review_functions.__all__) == []
    assert list(pyfaidx_review_functions.__all__) == []


def test_sbol3_namespace_helpers_share_one_global_state() -> None:
    sbol3 = pytest.importorskip("sbol3")
    original = sbol3.get_namespace()
    try:
        # Class-module re-exports are the same objects as the top-level configuration helpers.
        assert sbol3.component.set_namespace is sbol3.set_namespace
        assert sbol3.attachment.get_namespace is sbol3.get_namespace
        assert sbol3.document.set_defaults is sbol3.set_defaults

        sbol3.component.set_namespace("https://example.org/review")
        assert sbol3.get_namespace() == "https://example.org/review"
        sbol3.attachment.set_defaults()
        assert sbol3.get_namespace() is None
    finally:
        sbol3.set_namespace(original)


def test_sbol3_namespace_setter_rejects_non_urls_without_changing_state() -> None:
    sbol3 = pytest.importorskip("sbol3")
    original = sbol3.get_namespace()
    try:
        with pytest.raises(ValueError, match="Expected URL"):
            sbol3.set_namespace("not a url")
        assert sbol3.get_namespace() == original
    finally:
        sbol3.set_namespace(original)


def test_sbol3_data_path_resolves_inside_package_directory() -> None:
    sbol3 = pytest.importorskip("sbol3")
    package_dir = os.path.dirname(os.path.realpath(sbol3.document.__file__))
    assert sbol3.document.data_path("resources") == os.path.join(package_dir, "resources")


def test_hava_iter_reads_six_placeholder_fields_per_record() -> None:
    hava = pytest.importorskip("screed.hava")
    assert [name for name, _ in hava.FieldTypes] == [
        "hava",
        "quarzk",
        "muchalo",
        "fakours",
        "selimizicka",
        "marshoon",
    ]
    # Two hand-written records of six lines each; the iterator reuses one dict, so copy it.
    handle = BytesIO(b"r1\na\nb\nc\nd\ne\nr2\nf\ng\nh\ni\nj\n")
    records = [dict(item) for item in hava.hava_iter(handle)]
    assert records == [
        {
            "hava": "r1",
            "quarzk": "a",
            "muchalo": "b",
            "fakours": "c",
            "selimizicka": "d",
            "marshoon": "e",
        },
        {
            "hava": "r2",
            "quarzk": "f",
            "muchalo": "g",
            "fakours": "h",
            "selimizicka": "i",
            "marshoon": "j",
        },
    ]


def test_write_fastx_pair_fails_upstream_and_format_fastx_records_writes_the_pair() -> None:
    screed = pytest.importorskip("screed")
    read1 = screed.Record(name="read1", sequence="ACGT")
    read2 = screed.Record(name="read2", sequence="GGCC")
    with pytest.raises(NameError, match="write_record"):
        screed.screedRecord.write_fastx_pair(read1, read2, BytesIO())

    # The pair is the same as two single-record writes, which is what format_fastx_records does.
    manual = BytesIO()
    screed.screedRecord.write_fastx(read1, manual)
    screed.screedRecord.write_fastx(read2, manual)
    pair = format_fastx_records(
        [{"name": "read1", "sequence": "ACGT"}, {"name": "read2", "sequence": "GGCC"}],
        output_format="fasta",
        verbose=False,
    )
    assert pair["text"] == ">read1\nACGT\n>read2\nGGCC\n"
    assert pair["text"] == manual.getvalue().decode("utf-8")
    assert pair["record_count"] == 2


def test_get_valid_filename_replaces_spaces_and_drops_other_characters() -> None:
    pyfaidx = pytest.importorskip("pyfaidx")
    assert pyfaidx.get_valid_filename("chromosome 6.fa") == "chromosome_6.fa"
    # Strip the ends, turn the space into '_', then drop ':' and '?'.
    assert pyfaidx.get_valid_filename("  my seq:1?.fa ") == "my_seq1.fa"


def test_getmtime_fsspec_for_local_paths_matches_os_getmtime(tmp_path: Path) -> None:
    pyfaidx = pytest.importorskip("pyfaidx")
    path = tmp_path / "reference.fa"
    path.write_text(">chr1\nACGT\n", encoding="utf-8")
    assert pyfaidx.getmtime_fsspec(str(path), None) == os.path.getmtime(path)
