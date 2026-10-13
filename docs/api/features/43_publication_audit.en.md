# 43) Publication audit for 320 DNA and genomics tools

This page audits the 320 deduplicated projects in `public_dna_tools.md` for journal use. Peer-reviewed articles, preprints, software archives, project documentation, and missing evidence are kept separate. Redistribution rights and scientific validation are independent fields; a DOI is not automatically a journal article.

## 1) Generate the complete audit

- **Purpose:** Produce one reviewable row per surveyed project with citation, license, and missing metadata.
- **Method:** Read the frozen survey, verified DOI records, and project sources without guessing a DOI online.
- **API:** `dnakit.integrations.audit_publication_evidence()`.
- **Input:** No arguments audits all projects; use `query`, `manuscript_category`, or `journal_citation_status` to filter.
- **Output:** Each row contains `journal_citation_status`, `manuscript_category`, `citation_evidence`, `license_status`, `required_metadata`, and `scientific_validation`.
- **Example code:**

```python
from dnakit.integrations import audit_publication_evidence

rows = audit_publication_evidence(verbose=False)
print(len(rows))
print(rows[0]["journal_citation_status"])
```

- **Example output:**

```text
320
journal_article
```

## 2) Review the five manuscript categories

- **Purpose:** Separate manuscript-ready citations, missing version/access metadata, project-only sources, missing reliable papers/DOIs, and license restrictions.
- **Method:** License restrictions take precedence in the manuscript category; citation and license evidence remain independent.
- **API:** `dnakit.integrations.publication_audit_summary()`.
- **Input:** None.
- **Output:** Counts by manuscript category, strict journal status, license status, and `ready_for_journal`.
- **Example code:**

```python
from dnakit.integrations import publication_audit_summary

summary = publication_audit_summary(verbose=False)
print(summary["manuscript_categories"])
print(summary["journal_citation_status"])
```

- **Example output:**

```text
{'manuscript_ready': 84, 'needs_version_or_access_date': 4,
 'project_only': 47, 'missing_reliable_citation': 150,
 'license_restricted': 35}
```

## 3) Check whether one project has a journal article

- **Purpose:** Verify that a tool can be cited with a peer-reviewed article.
- **Method:** Only `journal_article` is treated as a journal source. `preprint`, `software_archive`, and `project_only` remain explicitly non-journal.
- **API:** `dnakit.integrations.describe_publication_audit(name[required])`.
- **Input:** A surveyed project name such as `"BWA"`, `"pyfaidx"`, or `"dnaio"`.
- **Output:** DOI, publication type, required version/access metadata, license state, and scientific-validation boundary.
- **Example code:**

```python
from dnakit.integrations import describe_publication_audit

for name in ("BWA", "pyfaidx", "dnaio"):
    row = describe_publication_audit(name, verbose=False)
    print(name, row["journal_citation_status"], row["manuscript_category"])
```

- **Example output:**

```text
BWA journal_article manuscript_ready
pyfaidx preprint needs_version_or_access_date
dnaio software_archive needs_version_or_access_date
```

## 4) Metadata required in a paper

- **Purpose:** Keep the Methods section reproducible.
- **Method:** Journal articles require the exact software version; preprints and archives also require source-type labeling and access date; project-only sources require the official URL.
- **API:** The `required_metadata` field in each audit row.
- **Input:** One audit row.
- **Output:** Values such as `exact_software_version`, `access_date`, and `label_non_journal_source`.

## 5) License and scientific-validation boundary

- A DOI or paper identifies a citation source; it does not grant redistribution rights for code, models, databases, or services.
- `license_status` reports only the license evidence in the frozen survey. `restricted_or_unverified` requires a separate license review.
- `scientific_validation` is always `not_assessed_by_publication_audit`; citation auditing does not replace benchmarks, controls, statistical tests, or independent replication.
