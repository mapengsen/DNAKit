"""Wave 1d sbol3 rows: the DNA design-record operations are covered by existing molecular_lab wrappers.

No new module is needed; these tests check the covered wrappers in dnakit.molecular_lab.sbol3_functions
against pySBOL3 directly. They skip when pySBOL3 is not installed.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

sbol3 = pytest.importorskip("sbol3")
sbol3_functions = importlib.import_module("dnakit.molecular_lab.sbol3_functions")

DNA_CLASSES = [
    "component", "sequence", "feature", "location", "seqfeat", "subcomponent", "toplevel",
    "document", "collection", "combderiv", "compref", "varcomp", "localsub",
]


def test_display_id_checks_match_sbol3():
    for text in ("gene_1", "1gene", "", "promoter-A"):
        assert sbol3_functions.check_sbol3_display_id(text)["valid"] == bool(sbol3.is_valid_display_id(text))


def test_extract_display_id_matches_sbol3_for_non_urls():
    for text in ("promoter", ""):
        assert sbol3_functions.extract_sbol3_display_id(text)["display_id"] is None
        assert sbol3.extract_display_id(text) is None


def test_extract_display_id_refuses_invalid_url_tail():
    with pytest.raises(Exception):
        sbol3_functions.extract_sbol3_display_id("https://example.org/dnakit/w1d/promoter/1")


def test_to_display_id_matches_sbol3_string_to_display_id():
    for text in ("my part", "123abc", "gfp-1.2"):
        assert sbol3_functions.to_sbol3_display_id(text)["display_id"] == str(
            sbol3.util.string_to_display_id(text)
        )


def test_copy_design_matches_sbol3_copy(tmp_path: Path):
    sbol3.set_namespace("https://example.org/dnakit/w1d")
    doc = sbol3.Document()
    seq = sbol3.Sequence("promoter_seq", elements="ttgacaaaaagcaatt", encoding=sbol3.IUPAC_DNA_ENCODING)
    doc.add(seq)
    comp = sbol3.Component("promoter", sbol3.SBO_DNA)
    comp.sequences = [seq]
    doc.add(comp)
    source = tmp_path / "design.nt"
    doc.write(str(source), sbol3.SORTED_NTRIPLES)
    result = sbol3_functions.copy_sbol3_design(str(source), into_namespace="https://example.org/dnakit/copy")
    assert result["copied_count"] == 2
    assert sorted(row["display_id"] for row in result["copied"]) == ["promoter", "promoter_seq"]


@pytest.mark.parametrize("cls", DNA_CLASSES)
def test_dna_class_methods_are_the_top_level_sbol3_functions(cls: str):
    module = importlib.import_module(f"sbol3.{cls}")
    assert module.extract_display_id is sbol3.extract_display_id
    assert module.is_valid_display_id is sbol3.is_valid_display_id
