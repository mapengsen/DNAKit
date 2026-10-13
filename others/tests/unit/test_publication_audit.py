"""Coverage for the complete 320-record publication audit."""

from dnakit.integrations import (
    audit_publication_evidence,
    describe_publication_audit,
    publication_audit_summary,
)


def test_audit_covers_every_survey_record_and_five_categories() -> None:
    summary = publication_audit_summary(verbose=False)
    assert summary["records"] == 320
    assert summary["survey_complete"] is True
    assert set(summary["manuscript_categories"]) == {
        "manuscript_ready",
        "needs_version_or_access_date",
        "project_only",
        "missing_reliable_citation",
        "license_restricted",
    }
    assert summary["journal_citation_status"]["journal_article"] > 0
    assert summary["journal_citation_status"]["preprint"] > 0
    assert summary["journal_citation_status"]["software_archive"] > 0


def test_journal_status_does_not_promote_preprint_or_archive() -> None:
    pyfaidx = describe_publication_audit("pyfaidx", verbose=False)
    dnaio = describe_publication_audit("dnaio", verbose=False)
    assert pyfaidx["journal_citation_status"] == "preprint"
    assert dnaio["journal_citation_status"] == "software_archive"
    assert describe_publication_audit("Centrifuge", verbose=False)["journal_citation_status"] == "journal_article"
    assert pyfaidx["manuscript_category"] != "manuscript_ready"
    assert dnaio["manuscript_category"] != "manuscript_ready"


def test_verified_journal_record_exposes_independent_validation_and_license() -> None:
    row = describe_publication_audit("BWA", verbose=False)
    assert row["journal_citation_status"] == "journal_article"
    assert row["manuscript_category"] == "manuscript_ready"
    assert row["citation_evidence"]["doi"] == "10.1093/bioinformatics/btp324"
    assert row["scientific_validation"] == "not_assessed_by_publication_audit"
    assert "license_status" in row


def test_filters_are_reproducible() -> None:
    rows = audit_publication_evidence(journal_citation_status="journal_article", verbose=False)
    assert rows
    assert all(row["journal_citation_status"] == "journal_article" for row in rows)


def test_high_priority_dna_tools_have_confirmed_primary_sources() -> None:
    expected = {
        "HUMAnN": "10.7554/eLife.65088",
        "Galaxy（2026 年 2 月 25 日后的 MIT 版本）": "10.1093/nar/gkag469",
        "SURVIVOR": "10.1038/ncomms14061",
        "gget": "10.1093/bioinformatics/btac836",
        "DnaFeaturesViewer": "10.1093/bioinformatics/btaa213",
    }
    for name, doi in expected.items():
        row = describe_publication_audit(name, verbose=False)
        assert row["journal_citation_status"] == "journal_article"
        assert row["citation_evidence"]["doi"] == doi
