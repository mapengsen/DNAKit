"""Tests for manuscript citation evidence records."""

from dnakit.integrations import (
    build_software_citation,
    describe_publication_record,
    publication_readiness,
    publication_registry_summary,
)
from dnakit.exceptions import ConfigurationError


def test_registry_contains_recoverable_software_dois() -> None:
    assert publication_registry_summary(verbose=False)["records"] == 19
    assert describe_publication_record("pyfaidx", verbose=False)["doi"] == "10.7287/peerj.preprints.970v1"
    assert describe_publication_record("dnaio", verbose=False)["doi"] == "10.5281/zenodo.10548864"


def test_readiness_and_citation_include_exact_version() -> None:
    result = publication_readiness(
        "SeqKit",
        version="2.10.0",
        accessed="2026-10-12",
        verbose=False,
    )
    assert result["status"] == "ready_for_manuscript_record"
    citation = build_software_citation(
        "scikit-bio",
        "0.7.1",
        accessed="2026-10-12",
        verbose=False,
    )
    assert "0.7.1" in citation["citation"]
    assert citation["doi"] == "10.1038/s41592-025-02981-z"


def test_readiness_rejects_invalid_access_date() -> None:
    try:
        publication_readiness("Biopython", version="1.86", accessed="2026/10/12", verbose=False)
    except ConfigurationError:
        return
    raise AssertionError("invalid ISO access dates must be rejected")
