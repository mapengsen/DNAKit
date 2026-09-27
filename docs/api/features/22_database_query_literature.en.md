# Literature, preprint and clinical study queries

Search publications, preprints and clinical study records.

<!-- BEGIN GENERATED API: literature -->

<span id="literature"></span>

## 1) Search arXiv papers by keyword {#dnakit.databases.ar_xiv_search_papers}

- **Purpose:** Search arXiv papers by keyword.
- **Method:** Use `ArXiv` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ar_xiv_search_papers(*, query, limit=None, sort_by=None, sort_order=None, date_from=None, date_to=None, category=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `title`, `abstract`, `authors`, `published` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ar_xiv_search_papers(
    query='machine learning',
    limit=2,
    sort_by='relevance',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "title": "<string>",
    "abstract": "<string>",
    "authors": [
      "<string>"
    ],
    "published": "<string>",
    "updated": "<string>",
    "category": "<string>",
    "url": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.ar_xiv_search_papers-parameters)

## 2) Retrieve bioRxiv preprint metadata by DOI {#dnakit.databases.bio_rxiv_get_preprint}

- **Purpose:** Retrieve bioRxiv preprint metadata by DOI.
- **Method:** Use `BioRxiv` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.bio_rxiv_get_preprint(*, doi, server=None, verbose=True)`
- **Input:** Provide `doi`.
- **Output:** Returns a mapping with `doi`, `title`, `authors`, `author_corresponding` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.bio_rxiv_get_preprint(
    doi='10.1101/2023.12.01.569554',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "doi": "<string>",
  "title": "<string>",
  "authors": [
    "<string>"
  ],
  "author_corresponding": "<string>",
  "author_corresponding_institution": "<string>",
  "abstract": "<string>",
  "date": "<string>",
  "version": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.bio_rxiv_get_preprint-parameters)

## 3) List recent preprints within a date range {#dnakit.databases.bio_rxiv_list_recent_preprints}

- **Purpose:** List recent preprints within a date range.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.bio_rxiv_list_recent_preprints(*, start_date=None, end_date=None, server=None, cursor=None, query=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `messages`, `collection` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.bio_rxiv_list_recent_preprints(
    start_date='2024-01-01',
    end_date='2024-01-03',
    server='biorxiv',
    cursor=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "messages": [
    {
      "status": "<string>",
      "interval": "<string>",
      "cursor": "<integer>",
      "count": "<integer>",
      "count_new_papers": "<string>",
      "total": "<string>"
    }
  ],
  "collection": [
    {
      "doi": "<string>",
      "title": "<string>",
      "authors": "<string>",
      "author_corresponding": "<string>",
      "author_corresponding_institution": "<string>",
      "date": "<string>",
      "version": "<string>",
      "category": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.bio_rxiv_list_recent_preprints-parameters)

## 4) Search clinical studies by condition, intervention or sponsor {#dnakit.databases.clinical_trials_search_studies}

- **Purpose:** Search clinical studies by condition, intervention or sponsor.
- **Method:** Use `ClinicalTrials` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.clinical_trials_search_studies(*, operation=None, query_cond=None, query_intr=None, query_term=None, filter_status=None, filter_phase=None, filter_study_type=None, page_size=None, next_page_token=None, query=None, condition=None, status=None, max_results=None, intervention=None, sponsor=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `studies`, `total_count`, `total_count_note`, `next_page_token` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.clinical_trials_search_studies(
    query_cond='breast cancer',
    query_intr='pembrolizumab',
    filter_status='RECRUITING',
    page_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "studies": [
    {
      "nct_id": "<string>",
      "brief_title": "<string>",
      "status": "<string>",
      "study_type": "<string>",
      "phases": "<array>",
      "enrollment": "<integer>",
      "conditions": "<array>",
      "interventions": "<array>"
    }
  ],
  "total_count": "<integer>",
  "total_count_note": "<string>",
  "next_page_token": "<string>",
  "executed_query": "<object>",
  "query_rewrites": [
    {
      "parameter": "<string>",
      "api_field": "<string>",
      "submitted": "<string>",
      "executed": "<string>",
      "restricted_to_fields": "<array>"
    }
  ],
  "relaxed_match_check": {
    "relaxed_query": "<object>",
    "relaxed_total_count": "<integer>",
    "note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.clinical_trials_search_studies-parameters)

## 5) Search Crossref publications and DOI metadata {#dnakit.databases.crossref_search_works}

- **Purpose:** Search Crossref publications and DOI metadata.
- **Method:** Use `Crossref` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.crossref_search_works(*, query, limit=None, offset=None, filter=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `DOI`, `title`, `author`, `container-title` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.crossref_search_works(
    query='artificial intelligence',
    limit=2,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "DOI": "<string>",
    "title": [
      "<string>"
    ],
    "author": [
      "<value>"
    ],
    "container-title": [
      "<value>"
    ],
    "published": "<object>",
    "publisher": "<string>",
    "type": "<string>",
    "is-referenced-by-count": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.crossref_search_works-parameters)

## 6) Search Europe PMC articles and abstracts {#dnakit.databases.europe_pmc_search_articles}

- **Purpose:** Search Europe PMC articles and abstracts.
- **Method:** Use `EuropePMC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.europe_pmc_search_articles(*, query, limit=None, require_has_ft=None, fulltext_terms=None, enrich_missing_abstract=None, extract_terms_from_fulltext=None, page_size=None, pageSize=None, source=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `title`, `abstract`, `abstract_retrieval_trace`, `authors` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.europe_pmc_search_articles(
    query='cancer therapy',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "title": "<string>",
    "abstract": "<string>",
    "abstract_retrieval_trace": [
      "<object>"
    ],
    "authors": [
      "<string>"
    ],
    "journal": "<string>",
    "year": "<integer>",
    "doi": "<string>",
    "doi_url": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.europe_pmc_search_articles-parameters)

## 7) Retrieve medRxiv preprint metadata by DOI {#dnakit.databases.med_rxiv_get_preprint}

- **Purpose:** Retrieve medRxiv preprint metadata by DOI.
- **Method:** Use `BioRxiv` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.med_rxiv_get_preprint(*, doi=None, server=None, query=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `doi`, `title`, `authors`, `author_corresponding` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.med_rxiv_get_preprint(
    doi='10.1101/2021.04.29.21256344',
    server='medrxiv',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "doi": "<string>",
  "title": "<string>",
  "authors": [
    "<string>"
  ],
  "author_corresponding": "<string>",
  "author_corresponding_institution": "<string>",
  "abstract": "<string>",
  "date": "<string>",
  "version": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.med_rxiv_get_preprint-parameters)

## 8) Retrieve OpenAlex author details {#dnakit.databases.openalex_get_author}

- **Purpose:** Retrieve OpenAlex author details.
- **Method:** Use `OpenAlex` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.openalex_get_author(*, author_id, mailto=None, verbose=True)`
- **Input:** Provide `author_id`.
- **Output:** Returns an OpenAlex author record; the declared response also allows null.
- **Example code:**

```python
import dnakit

result = dnakit.databases.openalex_get_author(
    author_id='A5001226970',
    verbose=False,
)
print(result)
```

- **Example output:** Returns an OpenAlex author record; the declared response also allows null.

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_get_author-parameters)

## 9) Search scholarly literature with OpenAlex {#dnakit.databases.openalex_literature_search}

- **Purpose:** Search scholarly literature with OpenAlex.
- **Method:** Use `OpenAlex` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.openalex_literature_search(*, search_keywords=None, query=None, max_results=None, limit=None, num_results=None, year_from=None, start_year=None, year_to=None, end_year=None, open_access=None, require_has_fulltext=None, fulltext_terms=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select at least one argument group: `search_keywords` / `query`.
- **Output:** Returns a list of records; each record includes `title`, `abstract`, `authors`, `year` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.openalex_literature_search(
    search_keywords='deep learning',
    max_results=10,
    year_from=2020,
    year_to=2023,
    open_access=True,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "title": "<string>",
    "abstract": "<string>",
    "authors": [
      "<string>"
    ],
    "year": "<integer>",
    "organizations": [
      "<string>"
    ],
    "venue": "<string>",
    "doi": "<string>",
    "citation_count": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_literature_search-parameters)

## 10) Search OpenAlex works with filters {#dnakit.databases.openalex_search_works}

- **Purpose:** Search OpenAlex works with filters.
- **Method:** Use `OpenAlex` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.openalex_search_works(*, search=None, query=None, filter=None, require_has_fulltext=None, fulltext_terms=None, per_page=None, limit=None, page=None, sort=None, mailto=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select at least one argument group: `search` / `query`.
- **Output:** Returns a mapping with `meta`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.openalex_search_works(
    search='attention is all you need',
    per_page=3,
    page=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": "<object>",
  "results": [
    {
      "id": "<string>",
      "title": "<string>",
      "doi": "<string>",
      "publication_year": "<integer>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_search_works-parameters)

## 11) Retrieve full PubMed article metadata by PMID {#dnakit.databases.pub_med_get_article}

- **Purpose:** Retrieve full PubMed article metadata by PMID.
- **Method:** Use `PubMed` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_med_get_article(*, pmid, verbose=True)`
- **Input:** Provide `pmid`.
- **Output:** Returns article XML or parsed metadata for the PMID, including title, abstract, authors and journal information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_med_get_article(
    pmid='17284678',
    verbose=False,
)
print(result)
```

- **Example output:** Returns article XML or parsed metadata for the PMID, including title, abstract, authors and journal information.

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_article-parameters)

## 12) Find articles citing a PubMed paper {#dnakit.databases.pub_med_get_cited_by}

- **Purpose:** Find articles citing a PubMed paper.
- **Method:** Use `PubMed` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_med_get_cited_by(*, pmid, limit=None, verbose=True)`
- **Input:** Provide `pmid`.
- **Output:** Returns a list for “Find articles citing a PubMed paper”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_med_get_cited_by(
    pmid='20210808',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_cited_by-parameters)

## 13) Find articles related to a PubMed paper {#dnakit.databases.pub_med_get_related}

- **Purpose:** Find articles related to a PubMed paper.
- **Method:** Use `PubMed` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_med_get_related(*, pmid, limit=None, verbose=True)`
- **Input:** Provide `pmid`.
- **Output:** Returns a list for “Find articles related to a PubMed paper”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_med_get_related(
    pmid='20210808',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_related-parameters)

## 14) Search PubMed biomedical literature {#dnakit.databases.pub_med_search_articles}

- **Purpose:** Search PubMed biomedical literature.
- **Method:** Use `PubMed` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_med_search_articles(*, query, limit=None, mindate=None, maxdate=None, datetype=None, include_abstract=None, sort=None, max_results=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `pmid`, `title`, `authors`, `author_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_med_search_articles(
    query='cancer research',
    limit=2,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "pmid": "<string>",
    "title": "<string>",
    "authors": [
      "<string>"
    ],
    "author_count": "<integer>",
    "authors_truncated": "<boolean>",
    "journal": "<string>",
    "pub_date": "<string>",
    "pub_year": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_search_articles-parameters)

## 15) Find clinical trials by disease, drug or keyword {#dnakit.databases.search_clinical_trials}

- **Purpose:** Find clinical trials by disease, drug or keyword.
- **Method:** Use `ClinicalTrialsSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.search_clinical_trials(*, condition=None, intervention=None, query_term=None, pageSize=None, pageToken=None, overall_status=None, max_results=None, limit=None, status=None, keyword=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a string.
- **Example code:**

```python
import dnakit

result = dnakit.databases.search_clinical_trials(
    condition='diabetes mellitus',
    intervention='metformin',
    pageSize=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
"<string>"
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.search_clinical_trials-parameters)

## 16) Search papers, abstracts and summaries with Semantic Scholar {#dnakit.databases.semantic_scholar_search_papers}

- **Purpose:** Search papers, abstracts and summaries with Semantic Scholar.
- **Method:** Use `SemanticScholar` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.semantic_scholar_search_papers(*, query, limit=None, year=None, sort=None, include_abstract=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `title`, `abstract`, `tldr`, `journal` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.semantic_scholar_search_papers(
    query='machine learning optimization',
    limit=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "title": "<string>",
    "abstract": "<string>",
    "tldr": "<string>",
    "journal": "<string>",
    "year": "<integer>",
    "url": "<string>",
    "paper_id": "<string>",
    "doi": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.semantic_scholar_search_papers-parameters)

## 17) Search the web using search engines {#dnakit.databases.web_search}

- **Purpose:** Search the web using search engines.
- **Method:** Use `WebSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.web_search(*, query, max_results=None, search_type=None, backend=None, region=None, safesearch=None, verbose=True)`
- **Input:** Provide `query`. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns a mapping with `status`, `query`, `search_type`, `total_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.web_search(
    query='NCBI genome assembly documentation',
    max_results=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "query": "<string>",
  "search_type": "<string>",
  "total_results": "<integer>",
  "results": [
    {
      "title": "<string>",
      "url": "<string>",
      "snippet": "<string>",
      "rank": "<integer>"
    }
  ],
  "backend_used": "<string>",
  "attempted_backends": [
    "<string>"
  ],
  "provider_errors": {
    "duckduckgo": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_literature.md#dnakit.databases.web_search-parameters)

<!-- END GENERATED API: literature -->
