# Literature, preprints and study discovery

This page documents **17 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

## Function index

- [`dnakit.databases.ar_xiv_search_papers`](#dnakit.databases.ar_xiv_search_papers)
- [`dnakit.databases.bio_rxiv_get_preprint`](#dnakit.databases.bio_rxiv_get_preprint)
- [`dnakit.databases.bio_rxiv_list_recent_preprints`](#dnakit.databases.bio_rxiv_list_recent_preprints)
- [`dnakit.databases.clinical_trials_search_studies`](#dnakit.databases.clinical_trials_search_studies)
- [`dnakit.databases.crossref_search_works`](#dnakit.databases.crossref_search_works)
- [`dnakit.databases.europe_pmc_search_articles`](#dnakit.databases.europe_pmc_search_articles)
- [`dnakit.databases.med_rxiv_get_preprint`](#dnakit.databases.med_rxiv_get_preprint)
- [`dnakit.databases.openalex_get_author`](#dnakit.databases.openalex_get_author)
- [`dnakit.databases.openalex_literature_search`](#dnakit.databases.openalex_literature_search)
- [`dnakit.databases.openalex_search_works`](#dnakit.databases.openalex_search_works)
- [`dnakit.databases.pub_med_get_article`](#dnakit.databases.pub_med_get_article)
- [`dnakit.databases.pub_med_get_cited_by`](#dnakit.databases.pub_med_get_cited_by)
- [`dnakit.databases.pub_med_get_related`](#dnakit.databases.pub_med_get_related)
- [`dnakit.databases.pub_med_search_articles`](#dnakit.databases.pub_med_search_articles)
- [`dnakit.databases.search_clinical_trials`](#dnakit.databases.search_clinical_trials)
- [`dnakit.databases.semantic_scholar_search_papers`](#dnakit.databases.semantic_scholar_search_papers)
- [`dnakit.databases.web_search`](#dnakit.databases.web_search)

## `dnakit.databases.ar_xiv_search_papers` {#dnakit.databases.ar_xiv_search_papers}

Search arXiv for papers by keyword using the public arXiv API. Returns papers with title, abstract, authors, publication date, category, and URL.

```text
dnakit.databases.ar_xiv_search_papers(*, query, limit=None, sort_by=None, sort_order=None, date_from=None, date_to=None, category=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query for arXiv papers. Unquoted, unprefixed space-separated keywords are ANDed together (every keyword must match) -- adding more keywords NARROWS results, it does not broaden them. To search multiple synonyms/phrases and broaden recall, join them with explicit OR, e.g. &#x27;&quot;machine learning&quot; OR &quot;deep learning&quot; OR &quot;neural networks&quot;&#x27;. Supports field prefixes (au:, ti:, abs:, cat:, all:) and AND/OR/ANDNOT boolean operators, which are passed through as-is. |
| `limit` | integer | None | Number of papers to return. This sets the maximum number of papers retrieved from arXiv. {&quot;default&quot;: 10} |
| `sort_by` | string | None | Sort order for results. Options: &#x27;relevance&#x27;, &#x27;lastUpdatedDate&#x27;, &#x27;submittedDate&#x27; {&quot;default&quot;: &quot;relevance&quot;} |
| `sort_order` | string | None | Sort direction. Options: &#x27;ascending&#x27;, &#x27;descending&#x27; {&quot;default&quot;: &quot;descending&quot;} |
| `date_from` | string | None | Filter results from this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper&#x27;s FIRST version -- a paper revised/updated after this date but originally submitted earlier is still excluded. |
| `date_to` | string | None | Filter results up to this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper&#x27;s FIRST version -- a paper revised/updated after this date but originally submitted earlier is still included. |
| `category` | string | None | Category or result class defined by this operation&#x27;s provider/schema. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query for arXiv papers. Unquoted, unprefixed space-separated keywords are ANDed together (every keyword must match) -- adding more keywords NARROWS results, it does not broaden them. To search multiple synonyms/phrases and broaden recall, join them with explicit OR, e.g. '\"machine learning\" OR \"deep learning\" OR \"neural networks\"'. Supports field prefixes (au:, ti:, abs:, cat:, all:) and AND/OR/ANDNOT boolean operators, which are passed through as-is."
        },
        "limit": {
          "type": "integer",
          "description": "Number of papers to return. This sets the maximum number of papers retrieved from arXiv.",
          "default": 10
        },
        "sort_by": {
          "type": "string",
          "description": "Sort order for results. Options: 'relevance', 'lastUpdatedDate', 'submittedDate'",
          "default": "relevance"
        },
        "sort_order": {
          "type": "string",
          "description": "Sort direction. Options: 'ascending', 'descending'",
          "default": "descending"
        },
        "date_from": {
          "type": "string",
          "description": "Filter results from this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper's FIRST version -- a paper revised/updated after this date but originally submitted earlier is still excluded."
        },
        "date_to": {
          "type": "string",
          "description": "Filter results up to this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper's FIRST version -- a paper revised/updated after this date but originally submitted earlier is still included."
        },
        "category": {
          "type": "string"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "List of arXiv papers matching the search query",
          "items": {
            "type": "object",
            "properties": {
              "title": {
                "type": "string",
                "description": "Title of the paper"
              },
              "abstract": {
                "type": "string",
                "description": "Abstract of the paper"
              },
              "authors": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of author names"
              },
              "published": {
                "type": "string",
                "description": "Publication date in ISO format"
              },
              "updated": {
                "type": "string",
                "description": "Last updated date in ISO format"
              },
              "category": {
                "type": "string",
                "description": "arXiv category (e.g., cs.LG, math.CO)"
              },
              "url": {
                "type": "string",
                "description": "URL to the paper on arXiv"
              }
            }
          }
        },
        {
          "type": "object",
          "description": "Error response",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "parameters": {
        "category": {
          "type": "string"
        }
      },
      "handler": "arxiv_category"
    }
    ```

**Execution and provenance:** `network`；upstream name `ArXiv_search_papers`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ar_xiv_search_papers", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.bio_rxiv_get_preprint` {#dnakit.databases.bio_rxiv_get_preprint}

Get full metadata for a specific bioRxiv or medRxiv preprint by its DOI. Returns comprehensive metadata including title, authors, abstract, publication date, version, category, corresponding author, funding information, and URLs (HTML, PDF, XML). Response time: <1 second. Use when you have a bioRxiv DOI (format: 10.1101/YYYY.MM.DD.XXXXXX). To SEARCH for preprints by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter (fast, relevance-ranked), or web_search with 'site:biorxiv.org'. bioRxiv API does not support text search.

```text
dnakit.databases.bio_rxiv_get_preprint(*, doi, server=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `doi` | string | Required by input schema | bioRxiv or medRxiv DOI. Can be full DOI (e.g., &#x27;10.1101/2023.12.01.569554&#x27;) or just the numeric part after &#x27;10.1101/&#x27; (e.g., &#x27;2023.12.01.569554&#x27;). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations. |
| `server` | string | None | Server to query: &#x27;biorxiv&#x27; for bioRxiv preprints or &#x27;medrxiv&#x27; for medRxiv preprints. Default is &#x27;biorxiv&#x27;. {&quot;default&quot;: &quot;biorxiv&quot;, &quot;enum&quot;: &#91;&quot;biorxiv&quot;, &quot;medrxiv&quot;&#93;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "doi": {
          "type": "string",
          "description": "bioRxiv or medRxiv DOI. Can be full DOI (e.g., '10.1101/2023.12.01.569554') or just the numeric part after '10.1101/' (e.g., '2023.12.01.569554'). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations."
        },
        "server": {
          "type": "string",
          "description": "Server to query: 'biorxiv' for bioRxiv preprints or 'medrxiv' for medRxiv preprints. Default is 'biorxiv'.",
          "default": "biorxiv",
          "enum": [
            "biorxiv",
            "medrxiv"
          ]
        }
      },
      "required": [
        "doi"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Full preprint metadata from bioRxiv/medRxiv",
      "properties": {
        "doi": {
          "type": "string",
          "description": "Digital Object Identifier"
        },
        "title": {
          "type": "string",
          "description": "Title of the preprint"
        },
        "authors": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of author names"
        },
        "author_corresponding": {
          "type": [
            "string",
            "null"
          ],
          "description": "Corresponding author name"
        },
        "author_corresponding_institution": {
          "type": [
            "string",
            "null"
          ],
          "description": "Corresponding author's institution"
        },
        "abstract": {
          "type": "string",
          "description": "Full abstract text"
        },
        "date": {
          "type": "string",
          "description": "Publication date (YYYY-MM-DD format)"
        },
        "version": {
          "type": "string",
          "description": "Version number of the preprint"
        },
        "type": {
          "type": "string",
          "description": "Type of submission (e.g., 'new results', 'confirmatory')"
        },
        "license": {
          "type": "string",
          "description": "License type (e.g., 'cc_by', 'cc0')"
        },
        "category": {
          "type": "string",
          "description": "Subject category (e.g., 'microbiology', 'neuroscience')"
        },
        "published": {
          "type": [
            "string",
            "null"
          ],
          "description": "Published DOI if the preprint has been published in a journal"
        },
        "url": {
          "type": [
            "string",
            "null"
          ],
          "description": "URL to the preprint on bioRxiv/medRxiv"
        },
        "pdf_url": {
          "type": [
            "string",
            "null"
          ],
          "description": "Direct PDF download URL"
        },
        "xml_url": {
          "type": [
            "string",
            "null"
          ],
          "description": "JATS XML source URL"
        },
        "server": {
          "type": "string",
          "description": "Server the preprint is hosted on ('biorxiv' or 'medrxiv')"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `BioRxiv_get_preprint`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("bio_rxiv_get_preprint", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.bio_rxiv_list_recent_preprints` {#dnakit.databases.bio_rxiv_list_recent_preprints}

List recent bioRxiv or medRxiv preprints by date range. Returns up to 100 preprints posted between the given start and end dates, with title, authors, abstract, DOI, category, and publication status. The bioRxiv API returns results in batches of up to 100; use the cursor parameter for pagination (0 for first page, 100 for second, etc.). Use this to monitor new preprints in a date range or discover recent work in a field. For searching by DOI, use BioRxiv_get_preprint instead. For searching by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter.

```text
dnakit.databases.bio_rxiv_list_recent_preprints(*, start_date=None, end_date=None, server=None, cursor=None, query=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `start_date` | string | None | Start date in YYYY-MM-DD format (e.g., &#x27;2024-01-01&#x27;). Date range must not exceed 60 days. |
| `end_date` | string | None | End date in YYYY-MM-DD format (e.g., &#x27;2024-01-03&#x27;). Date range must not exceed 60 days. |
| `server` | string | None | Server: &#x27;biorxiv&#x27; for biology preprints, &#x27;medrxiv&#x27; for health sciences preprints {&quot;default&quot;: &quot;biorxiv&quot;, &quot;enum&quot;: &#91;&quot;biorxiv&quot;, &quot;medrxiv&quot;&#93;} |
| `cursor` | integer | None | Pagination cursor (0 for first 100 results, 100 for next 100, etc.) {&quot;default&quot;: 0} |
| `query` | string | None | Search term or expression in the current provider&#x27;s query syntax. |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "start_date": {
          "type": "string",
          "description": "Start date in YYYY-MM-DD format (e.g., '2024-01-01'). Date range must not exceed 60 days."
        },
        "end_date": {
          "type": "string",
          "description": "End date in YYYY-MM-DD format (e.g., '2024-01-03'). Date range must not exceed 60 days."
        },
        "server": {
          "type": "string",
          "description": "Server: 'biorxiv' for biology preprints, 'medrxiv' for health sciences preprints",
          "default": "biorxiv",
          "enum": [
            "biorxiv",
            "medrxiv"
          ]
        },
        "cursor": {
          "type": "integer",
          "description": "Pagination cursor (0 for first 100 results, 100 for next 100, etc.)",
          "default": 0
        },
        "query": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "description": "bioRxiv/medRxiv recent preprints response",
          "properties": {
            "messages": {
              "type": "array",
              "description": "Status messages with pagination info",
              "items": {
                "type": "object",
                "properties": {
                  "status": {
                    "type": "string"
                  },
                  "interval": {
                    "type": "string",
                    "description": "Date range queried"
                  },
                  "cursor": {
                    "type": "integer",
                    "description": "Current cursor position"
                  },
                  "count": {
                    "type": "integer",
                    "description": "Number of results in this batch (max 100)"
                  },
                  "count_new_papers": {
                    "type": "string",
                    "description": "Total unique new papers in date range"
                  },
                  "total": {
                    "type": "string",
                    "description": "Total records including revisions"
                  }
                }
              }
            },
            "collection": {
              "type": "array",
              "description": "Array of preprint metadata",
              "items": {
                "type": "object",
                "properties": {
                  "doi": {
                    "type": "string",
                    "description": "bioRxiv/medRxiv DOI"
                  },
                  "title": {
                    "type": "string",
                    "description": "Preprint title"
                  },
                  "authors": {
                    "type": "string",
                    "description": "Author list as semicolon-separated string"
                  },
                  "author_corresponding": {
                    "type": "string",
                    "description": "Corresponding author name"
                  },
                  "author_corresponding_institution": {
                    "type": "string",
                    "description": "Corresponding author institution"
                  },
                  "date": {
                    "type": "string",
                    "description": "Publication date"
                  },
                  "version": {
                    "type": "string",
                    "description": "Version number"
                  },
                  "category": {
                    "type": "string",
                    "description": "Subject category (e.g., neuroscience, bioinformatics, genomics)"
                  },
                  "abstract": {
                    "type": "string",
                    "description": "Full abstract text"
                  },
                  "published": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Published journal DOI if accepted"
                  },
                  "license": {
                    "type": "string",
                    "description": "License type"
                  }
                }
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "parameters": {
        "query": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
        }
      },
      "optional": [
        "start_date",
        "end_date"
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `BioRxiv_list_recent_preprints`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("bio_rxiv_list_recent_preprints", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.clinical_trials_search_studies` {#dnakit.databases.clinical_trials_search_studies}

Search ClinicalTrials.gov for clinical trial studies by condition, intervention, sponsor, or other criteria. Returns NCT IDs, titles, status, phase, enrollment, and key trial metadata. Example: search for 'breast cancer' trials with 'pembrolizumab' that are RECRUITING in Phase 3. Query rewriting: multi-word intervention/condition/sponsor values are sent to ClinicalTrials.gov as an EXACT PHRASE, and intervention/sponsor values are further restricted to specific registry fields, so a phrase no registrant recorded verbatim can match 0 studies even when many trials mention it. The query actually executed is always returned in `executed_query` (with per-parameter detail in `query_rewrites`); when a rewritten query matches nothing, `relaxed_match_check` reports how many studies your original wording would have matched and how to retrieve them. To opt out of the rewrite, wrap a value in parentheses (e.g. '(donor derived cell therapy)') and it is sent verbatim as a loose keyword match.

```text
dnakit.databases.clinical_trials_search_studies(*, operation=None, query_cond=None, query_intr=None, query_term=None, filter_status=None, filter_phase=None, filter_study_type=None, page_size=None, next_page_token=None, query=None, condition=None, status=None, max_results=None, intervention=None, sponsor=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;search&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;search&quot;&#93;, &quot;default&quot;: &quot;search&quot;} |
| `query_cond` | string &#124; null | None | Disease or condition to search for (e.g., &#x27;breast cancer&#x27;, &#x27;type 2 diabetes&#x27;, &#x27;COVID-19&#x27;). Supports Boolean operators (AND, OR, NOT). Multi-word values are matched as an exact phrase (e.g. &#x27;kidney transplantation&#x27; is sent as &quot;kidney transplantation&quot;), which is narrower than a loose keyword match. For looser matching, wrap the value in parentheses (e.g. &#x27;(kidney transplantation)&#x27;) to send it verbatim. Values already containing Essie syntax (AND/OR/NOT/AREA&#91;&#93;/quotes) are passed through unchanged. See `executed_query` in the response for what actually ran. |
| `query_intr` | string &#124; null | None | Intervention or drug to search for (e.g., &#x27;pembrolizumab&#x27;, &#x27;metformin&#x27;, &#x27;chemotherapy&#x27;). Matched against each study&#x27;s registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so &#x27;donor derived cell therapy&#x27; matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. &#x27;(donor derived cell therapy)&#x27;) to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response. |
| `query_term` | string &#124; null | None | Free-text search across all study fields. Use for general keyword search. |
| `filter_status` | string &#124; null | None | Filter by recruitment status. Options: &#x27;RECRUITING&#x27;, &#x27;NOT_YET_RECRUITING&#x27;, &#x27;ACTIVE_NOT_RECRUITING&#x27;, &#x27;COMPLETED&#x27;, &#x27;WITHDRAWN&#x27;, &#x27;SUSPENDED&#x27;, &#x27;TERMINATED&#x27;, &#x27;ENROLLING_BY_INVITATION&#x27;. Comma-separate multiple values. |
| `filter_phase` | string &#124; null | None | Filter by trial phase. Options: &#x27;EARLY_PHASE1&#x27;, &#x27;PHASE1&#x27;, &#x27;PHASE2&#x27;, &#x27;PHASE3&#x27;, &#x27;PHASE4&#x27;, &#x27;NA&#x27;. Comma-separate multiple (e.g., &#x27;PHASE2,PHASE3&#x27;). |
| `filter_study_type` | string &#124; null | None | Filter by study type: &#x27;INTERVENTIONAL&#x27;, &#x27;OBSERVATIONAL&#x27;, &#x27;EXPANDED_ACCESS&#x27;. |
| `page_size` | integer | None | Number of results per page (default 10, max 1000). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `next_page_token` | string &#124; null | None | Token for retrieving the next page of results (obtained from previous response). |
| `query` | string &#124; null | None | General keyword search across all fields. Alias for query_term. E.g., &quot;FLT3 AML&quot;, &quot;KRAS NSCLC&quot;. |
| `condition` | string &#124; null | None | Disease or condition to search for. Alias for query_cond. E.g., &quot;acute myeloid leukemia&quot;. Multi-word values are matched as an exact phrase (e.g. &#x27;kidney transplantation&#x27; is sent as &quot;kidney transplantation&quot;), which is narrower than a loose keyword match. For looser matching, wrap the value in parentheses (e.g. &#x27;(kidney transplantation)&#x27;) to send it verbatim. Values already containing Essie syntax (AND/OR/NOT/AREA&#91;&#93;/quotes) are passed through unchanged. See `executed_query` in the response for what actually ran. |
| `status` | string &#124; null | None | Recruitment status filter. Alias for filter_status. E.g., &quot;RECRUITING&quot;, &quot;COMPLETED&quot;. |
| `max_results` | integer &#124; null | None | Maximum number of results to return. Alias for page_size. |
| `intervention` | string &#124; null | None | Intervention or drug to search for. Alias for query_intr. E.g., &quot;pembrolizumab&quot;, &quot;metformin&quot;. Matched against each study&#x27;s registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so &#x27;donor derived cell therapy&#x27; matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. &#x27;(donor derived cell therapy)&#x27;) to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response. |
| `sponsor` | string &#124; null | None | Lead sponsor or collaborator name to search for. E.g., &quot;National Cancer Institute&quot;. Matched against each study&#x27;s registered lead-sponsor name (AREA&#91;LeadSponsorName&#93;) only, and multi-word values are matched as an EXACT PHRASE. For looser matching, wrap the value in parentheses (e.g. &#x27;(National Cancer Institute)&#x27;) to send it verbatim. See `executed_query` and `relaxed_match_check` in the response. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "operation": {
          "type": "string",
          "enum": [
            "search"
          ],
          "description": "Operation (optional; defaults to 'search' for this tool).",
          "default": "search"
        },
        "query_cond": {
          "type": [
            "string",
            "null"
          ],
          "description": "Disease or condition to search for (e.g., 'breast cancer', 'type 2 diabetes', 'COVID-19'). Supports Boolean operators (AND, OR, NOT). Multi-word values are matched as an exact phrase (e.g. 'kidney transplantation' is sent as \"kidney transplantation\"), which is narrower than a loose keyword match. For looser matching, wrap the value in parentheses (e.g. '(kidney transplantation)') to send it verbatim. Values already containing Essie syntax (AND/OR/NOT/AREA[]/quotes) are passed through unchanged. See `executed_query` in the response for what actually ran."
        },
        "query_intr": {
          "type": [
            "string",
            "null"
          ],
          "description": "Intervention or drug to search for (e.g., 'pembrolizumab', 'metformin', 'chemotherapy'). Matched against each study's registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so 'donor derived cell therapy' matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. '(donor derived cell therapy)') to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response."
        },
        "query_term": {
          "type": [
            "string",
            "null"
          ],
          "description": "Free-text search across all study fields. Use for general keyword search."
        },
        "filter_status": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by recruitment status. Options: 'RECRUITING', 'NOT_YET_RECRUITING', 'ACTIVE_NOT_RECRUITING', 'COMPLETED', 'WITHDRAWN', 'SUSPENDED', 'TERMINATED', 'ENROLLING_BY_INVITATION'. Comma-separate multiple values."
        },
        "filter_phase": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by trial phase. Options: 'EARLY_PHASE1', 'PHASE1', 'PHASE2', 'PHASE3', 'PHASE4', 'NA'. Comma-separate multiple (e.g., 'PHASE2,PHASE3')."
        },
        "filter_study_type": {
          "type": [
            "string",
            "null"
          ],
          "description": "Filter by study type: 'INTERVENTIONAL', 'OBSERVATIONAL', 'EXPANDED_ACCESS'."
        },
        "page_size": {
          "type": "integer",
          "description": "Number of results per page (default 10, max 1000).",
          "default": 10,
          "minimum": 1,
          "maximum": 1000
        },
        "next_page_token": {
          "type": [
            "string",
            "null"
          ],
          "description": "Token for retrieving the next page of results (obtained from previous response)."
        },
        "query": {
          "type": [
            "string",
            "null"
          ],
          "description": "General keyword search across all fields. Alias for query_term. E.g., \"FLT3 AML\", \"KRAS NSCLC\"."
        },
        "condition": {
          "type": [
            "string",
            "null"
          ],
          "description": "Disease or condition to search for. Alias for query_cond. E.g., \"acute myeloid leukemia\". Multi-word values are matched as an exact phrase (e.g. 'kidney transplantation' is sent as \"kidney transplantation\"), which is narrower than a loose keyword match. For looser matching, wrap the value in parentheses (e.g. '(kidney transplantation)') to send it verbatim. Values already containing Essie syntax (AND/OR/NOT/AREA[]/quotes) are passed through unchanged. See `executed_query` in the response for what actually ran."
        },
        "status": {
          "type": [
            "string",
            "null"
          ],
          "description": "Recruitment status filter. Alias for filter_status. E.g., \"RECRUITING\", \"COMPLETED\"."
        },
        "max_results": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of results to return. Alias for page_size."
        },
        "intervention": {
          "type": [
            "string",
            "null"
          ],
          "description": "Intervention or drug to search for. Alias for query_intr. E.g., \"pembrolizumab\", \"metformin\". Matched against each study's registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so 'donor derived cell therapy' matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. '(donor derived cell therapy)') to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response."
        },
        "sponsor": {
          "type": [
            "string",
            "null"
          ],
          "description": "Lead sponsor or collaborator name to search for. E.g., \"National Cancer Institute\". Matched against each study's registered lead-sponsor name (AREA[LeadSponsorName]) only, and multi-word values are matched as an EXACT PHRASE. For looser matching, wrap the value in parentheses (e.g. '(National Cancer Institute)') to send it verbatim. See `executed_query` and `relaxed_match_check` in the response."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "studies": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "nct_id": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "brief_title": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "status": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "study_type": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "phases": {
                    "type": "array"
                  },
                  "enrollment": {
                    "type": [
                      "integer",
                      "null"
                    ]
                  },
                  "conditions": {
                    "type": "array"
                  },
                  "interventions": {
                    "type": "array"
                  },
                  "sponsor": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "start_date": {
                    "type": [
                      "string",
                      "null"
                    ]
                  },
                  "completion_date": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "total_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Studies matching the query, independent of page_size. NULL on any page fetched with next_page_token: ClinicalTrials.gov reports this figure only on a first page, and a page's row count is not a substitute for it. Read the total from the first page of the same query. When null, total_count_note says so."
            },
            "total_count_note": {
              "type": "string",
              "description": "Present only when total_count is null, explaining that the registry does not report a matching-study count on a paged response and where to read the real total."
            },
            "next_page_token": {
              "type": [
                "string",
                "null"
              ]
            },
            "executed_query": {
              "type": "object",
              "description": "The query parameters actually sent to ClinicalTrials.gov, after phrase-quoting and AREA[] field restriction. Compare against what you passed: they differ whenever `query_rewrites` is non-empty."
            },
            "query_rewrites": {
              "type": "array",
              "description": "One entry per parameter whose value was rewritten before being sent. Empty when the query ran exactly as supplied.",
              "items": {
                "type": "object",
                "properties": {
                  "parameter": {
                    "type": "string",
                    "description": "The parameter name you passed."
                  },
                  "api_field": {
                    "type": "string",
                    "description": "The ClinicalTrials.gov query field it was sent as."
                  },
                  "submitted": {
                    "type": "string",
                    "description": "The value you supplied."
                  },
                  "executed": {
                    "type": "string",
                    "description": "The rewritten value that was actually sent."
                  },
                  "restricted_to_fields": {
                    "type": "array",
                    "description": "Registry fields the value was restricted to, if any.",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "relaxed_match_check": {
              "type": "object",
              "description": "Present ONLY when a rewritten query returned zero studies. Reports whether your original wording would have matched, so an empty `studies` list can be told apart from a genuine absence of trials. The relaxed results are never substituted for the ones you asked for -- follow `note` to retrieve them.",
              "properties": {
                "relaxed_query": {
                  "type": "object",
                  "description": "The un-rewritten query used for the check."
                },
                "relaxed_total_count": {
                  "type": [
                    "integer",
                    "null"
                  ],
                  "description": "Studies matching your original wording; null if the check request failed."
                },
                "note": {
                  "type": "string",
                  "description": "Plain-language reading of the check, naming the escape hatch when the zero is an artifact of the rewrite."
                }
              }
            }
          },
          "required": [
            "studies"
          ]
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `ClinicalTrials_search_studies`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("clinical_trials_search_studies", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.crossref_search_works` {#dnakit.databases.crossref_search_works}

Search Crossref Works API for scholarly articles, publications, and research outputs by keyword. Returns metadata including title, abstract, journal, authors, publication year, DOI, URL, citations, and publisher. Supports filtering by publication type (journal article, conference paper, etc.) and date range. Use this to discover academic literature, verify citations, or build bibliographies. Find DOIs to use with Crossref_get_work for complete metadata. The response reports 'count' (works in this page) alongside 'total_count' (works matching the query in all of Crossref) -- these differ on almost every real search, since 'limit' defaults to 10 and caps at 100, so treat 'count' as a page size and never as the size of the literature. When they differ, 'truncated' is true and 'truncation_note' gives the 'offset' for the next page. Note: Crossref's own API weakens keyword relevance-matching whenever a non-default 'sort' query param is combined with 'query' (e.g. sort=is-referenced-by-count), returning globally top-cited works largely unrelated to the search term -- omit 'sort' (the default relevance ranking) unless you have separately verified the combination for your query.

```text
dnakit.databases.crossref_search_works(*, query, limit=None, offset=None, filter=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query for Crossref works. Use keywords separated by spaces to refine your search across titles, abstracts, authors, and other bibliographic fields. |
| `limit` | integer | None | Number of articles to return in this page. This sets the maximum number of articles retrieved from Crossref. Max 100 per request. It does not change &#x27;total_count&#x27;, which always reports how many works match the query. {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Skip this many matching works before returning results. Use with &#x27;limit&#x27; to page through a result set larger than 100; the response &#x27;truncation_note&#x27; gives the offset for the next page. {&quot;minimum&quot;: 0} |
| `filter` | string | None | Optional filter string for Crossref API. Examples: &#x27;type:journal-article&#x27; (only journal articles), &#x27;from-pub-date:2020-01-01&#x27; (published after date), &#x27;has-abstract:true&#x27; (only with abstracts). Multiple filters comma-separated. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query for Crossref works. Use keywords separated by spaces to refine your search across titles, abstracts, authors, and other bibliographic fields."
        },
        "limit": {
          "type": "integer",
          "description": "Number of articles to return in this page. This sets the maximum number of articles retrieved from Crossref. Max 100 per request. It does not change 'total_count', which always reports how many works match the query.",
          "default": 10,
          "minimum": 1,
          "maximum": 100
        },
        "offset": {
          "type": "integer",
          "minimum": 0,
          "description": "Skip this many matching works before returning results. Use with 'limit' to page through a result set larger than 100; the response 'truncation_note' gives the offset for the next page."
        },
        "filter": {
          "type": "string",
          "description": "Optional filter string for Crossref API. Examples: 'type:journal-article' (only journal articles), 'from-pub-date:2020-01-01' (published after date), 'has-abstract:true' (only with abstracts). Multiple filters comma-separated."
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array",
      "description": "Array of Crossref work objects (articles/publications)",
      "items": {
        "type": "object",
        "properties": {
          "DOI": {
            "type": "string"
          },
          "title": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "author": {
            "type": "array"
          },
          "container-title": {
            "type": "array",
            "description": "Journal name"
          },
          "published": {
            "type": "object"
          },
          "publisher": {
            "type": "string"
          },
          "type": {
            "type": "string",
            "description": "Publication type (e.g., journal-article)"
          },
          "is-referenced-by-count": {
            "type": "integer",
            "description": "Citation count"
          },
          "abstract": {
            "type": "string"
          },
          "reference": {
            "type": "array",
            "description": "List of references"
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Crossref_search_works`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("crossref_search_works", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.europe_pmc_search_articles` {#dnakit.databases.europe_pmc_search_articles}

Search for articles on Europe PMC including abstracts and metadata. Returns the top `limit` results ranked by relevance, which is normally a SLICE of the matches: `metadata.total_results` carries Europe PMC's own hitCount (every article matching the query) while `metadata.count` is only how many were returned, and a top-level `truncated: true` plus `truncation_note` appear whenever matches were left behind. Do not read a small `count` as 'few articles exist' — check `total_results`. Europe PMC supports fielded queries for indexed full text (e.g., BODY:"term") but only for records where Europe PMC has full text indexed (HAS_FT:Y). Use require_has_ft/fulltext_terms to avoid confusing 'has a full-text link somewhere' with 'full text searchable in Europe PMC'.

```text
dnakit.databases.europe_pmc_search_articles(*, query, limit=None, require_has_ft=None, fulltext_terms=None, enrich_missing_abstract=None, extract_terms_from_fulltext=None, page_size=None, pageSize=None, source=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query for Europe PMC. Supports Lucene-like fielded syntax (e.g., BODY:&quot;term&quot;, INTRO:&quot;term&quot;, METHODS:&quot;term&quot;). |
| `limit` | integer | None | Number of articles to return. This sets the maximum number of articles retrieved from Europe PMC; anything beyond it is silently left behind, so compare it against `metadata.total_results` in the response (Europe PMC&#x27;s hitCount) before treating the returned set as complete. Europe PMC accepts up to 1000 per request. {&quot;default&quot;: 5} |
| `require_has_ft` | boolean | None | If true, appends `HAS_FT:Y` to the query to restrict results to records where Europe PMC has full text indexed. This is NOT the same as having a link to free full text elsewhere (e.g., PMC) and is often the key reason body-only keywords cannot be found via Europe PMC search. {&quot;default&quot;: false} |
| `fulltext_terms` | array | None | Optional list of terms that must occur in the indexed full text (Europe PMC BODY field). When provided, the tool adds a BODY:&quot;...&quot; clause and automatically enables `require_has_ft`. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enrich_missing_abstract` | boolean | None | If true, best-effort fills missing abstracts by fetching Europe PMC fullTextXML (bounded to a few results). {&quot;default&quot;: false} |
| `extract_terms_from_fulltext` | array | None | Optional list of terms to extract from full text (open access only). When provided, automatically fetches fullTextXML for open-access articles and returns snippets around these terms. Terms are processed in batches of 5 internally, so any number of terms is accepted. Bounded to first 3 OA articles to avoid latency. Returns snippets in a &#x27;fulltext_snippets&#x27; field for each article. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `page_size` | integer | None | Alias for limit — number of results to return |
| `pageSize` | integer | None | DNAKit optional alias for limit |
| `source` | string | None | Source filter using values accepted by this operation. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query for Europe PMC. Supports Lucene-like fielded syntax (e.g., BODY:\"term\", INTRO:\"term\", METHODS:\"term\")."
        },
        "limit": {
          "type": "integer",
          "description": "Number of articles to return. This sets the maximum number of articles retrieved from Europe PMC; anything beyond it is silently left behind, so compare it against `metadata.total_results` in the response (Europe PMC's hitCount) before treating the returned set as complete. Europe PMC accepts up to 1000 per request.",
          "default": 5
        },
        "require_has_ft": {
          "type": "boolean",
          "description": "If true, appends `HAS_FT:Y` to the query to restrict results to records where Europe PMC has full text indexed. This is NOT the same as having a link to free full text elsewhere (e.g., PMC) and is often the key reason body-only keywords cannot be found via Europe PMC search.",
          "default": false
        },
        "fulltext_terms": {
          "type": "array",
          "description": "Optional list of terms that must occur in the indexed full text (Europe PMC BODY field). When provided, the tool adds a BODY:\"...\" clause and automatically enables `require_has_ft`.",
          "items": {
            "type": "string"
          }
        },
        "enrich_missing_abstract": {
          "type": "boolean",
          "description": "If true, best-effort fills missing abstracts by fetching Europe PMC fullTextXML (bounded to a few results).",
          "default": false
        },
        "extract_terms_from_fulltext": {
          "type": "array",
          "description": "Optional list of terms to extract from full text (open access only). When provided, automatically fetches fullTextXML for open-access articles and returns snippets around these terms. Terms are processed in batches of 5 internally, so any number of terms is accepted. Bounded to first 3 OA articles to avoid latency. Returns snippets in a 'fulltext_snippets' field for each article.",
          "items": {
            "type": "string"
          }
        },
        "page_size": {
          "type": "integer",
          "description": "Alias for limit — number of results to return"
        },
        "pageSize": {
          "type": "integer",
          "description": "DNAKit optional alias for limit"
        },
        "source": {
          "type": "string"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "List of articles from Europe PMC matching the search query",
          "items": {
            "type": "object",
            "properties": {
              "title": {
                "type": "string",
                "description": "Title of the article"
              },
              "abstract": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Abstract of the article"
              },
              "abstract_retrieval_trace": {
                "type": "array",
                "description": "Machine-readable trace of attempts used to retrieve/derive the abstract when enrichment is enabled.",
                "items": {
                  "type": "object",
                  "properties": {
                    "attempt": {
                      "type": "string"
                    },
                    "url": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "status_code": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "content_type": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "note": {
                      "type": [
                        "string",
                        "null"
                      ]
                    }
                  }
                }
              },
              "authors": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of author names"
              },
              "journal": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Journal name where the article was published"
              },
              "year": {
                "type": [
                  "integer",
                  "string",
                  "null"
                ],
                "description": "Publication year"
              },
              "doi": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Digital Object Identifier"
              },
              "doi_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "DOI resolver URL"
              },
              "url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "URL to the article on Europe PMC"
              },
              "source_db": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Europe PMC source database for the record (e.g., MED, PMC)."
              },
              "article_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Record identifier within the Europe PMC source database (e.g., PMID for MED)."
              },
              "pmid": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "PubMed ID (when available)"
              },
              "pmcid": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "PubMed Central ID (when available)"
              },
              "fulltext_xml_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Europe PMC fullTextXML URL (when derivable)"
              },
              "abstract_source": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Source of abstract text (e.g., Europe PMC core or fullTextXML)"
              },
              "fulltext_snippets": {
                "type": "array",
                "description": "Full-text snippets around search terms (only present when extract_terms_from_fulltext is provided and article is OA)",
                "items": {
                  "type": "object",
                  "properties": {
                    "term": {
                      "type": "string",
                      "description": "The search term that was matched"
                    },
                    "snippet": {
                      "type": "string",
                      "description": "Text snippet with context around the matched term"
                    }
                  }
                }
              },
              "fulltext_snippets_count": {
                "type": "integer",
                "description": "Number of snippets extracted (only present when fulltext_snippets exists)"
              },
              "fulltext_snippets_source": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Provenance for snippet extraction (e.g., Europe PMC fullTextXML, NCBI PMC OAI (JATS), NCBI PMC HTML)."
              },
              "fulltext_snippets_retrieval_trace": {
                "type": "array",
                "description": "Machine-readable trace of attempts used to retrieve full text for snippet extraction.",
                "items": {
                  "type": "object",
                  "properties": {
                    "attempt": {
                      "type": "string"
                    },
                    "url": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "status_code": {
                      "type": [
                        "integer",
                        "null"
                      ]
                    },
                    "content_type": {
                      "type": [
                        "string",
                        "null"
                      ]
                    },
                    "note": {
                      "type": [
                        "string",
                        "null"
                      ]
                    }
                  }
                }
              },
              "citations": {
                "type": "integer",
                "description": "Number of citations"
              },
              "open_access": {
                "type": "boolean",
                "description": "Open access status"
              },
              "keywords": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of keywords"
              },
              "source": {
                "type": "string",
                "description": "Source identifier (always 'Europe PMC')"
              },
              "data_quality": {
                "type": "object",
                "description": "Data quality indicators",
                "properties": {
                  "has_abstract": {
                    "type": "boolean"
                  },
                  "has_authors": {
                    "type": "boolean"
                  },
                  "has_journal": {
                    "type": "boolean"
                  },
                  "has_year": {
                    "type": "boolean"
                  },
                  "has_doi": {
                    "type": "boolean"
                  },
                  "has_citations": {
                    "type": "boolean"
                  },
                  "has_keywords": {
                    "type": "boolean"
                  },
                  "has_url": {
                    "type": "boolean"
                  }
                }
              }
            }
          },
          "required": [
            "data"
          ]
        },
        {
          "type": "object",
          "description": "Error response",
          "properties": {
            "error": {
              "type": "string",
              "description": "Error message"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "pageSize": "limit"
      },
      "parameters": {
        "source": {
          "type": "string"
        }
      },
      "handler": "literature_source"
    }
    ```

**Execution and provenance:** `network`；upstream name `EuropePMC_search_articles`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("europe_pmc_search_articles", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.med_rxiv_get_preprint` {#dnakit.databases.med_rxiv_get_preprint}

Get full metadata for a specific medRxiv preprint by its DOI. Returns comprehensive metadata including title, authors, abstract, publication date, version, category, and URLs. Response time: <1 second. Use when you have a medRxiv DOI (format: 10.1101/YYYY.MM.DD.XXXXXXX). To SEARCH for medRxiv preprints by keywords, use EuropePMC_search_articles with 'SRC:PPR' filter (fast, relevance-ranked), or web_search with 'site:medrxiv.org'. medRxiv API does not support text search.

```text
dnakit.databases.med_rxiv_get_preprint(*, doi=None, server=None, query=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `doi` | string | None | medRxiv DOI. Can be full DOI (e.g., &#x27;10.1101/2021.04.29.21256344&#x27;) or just the numeric part after &#x27;10.1101/&#x27; (e.g., &#x27;2021.04.29.21256344&#x27;). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations. |
| `server` | string | None | Server to query - always &#x27;medrxiv&#x27; for this tool. {&quot;default&quot;: &quot;medrxiv&quot;, &quot;enum&quot;: &#91;&quot;medrxiv&quot;&#93;} |
| `query` | string | None | Search term or expression in the current provider&#x27;s query syntax. |
| `limit` | integer | None | Maximum results for the current request, not a database-wide total. {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "doi": {
          "type": "string",
          "description": "medRxiv DOI. Can be full DOI (e.g., '10.1101/2021.04.29.21256344') or just the numeric part after '10.1101/' (e.g., '2021.04.29.21256344'). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations."
        },
        "server": {
          "type": "string",
          "description": "Server to query - always 'medrxiv' for this tool.",
          "default": "medrxiv",
          "enum": [
            "medrxiv"
          ]
        },
        "query": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "description": "Full preprint metadata from medRxiv",
      "properties": {
        "doi": {
          "type": "string",
          "description": "Digital Object Identifier"
        },
        "title": {
          "type": "string",
          "description": "Title of the preprint"
        },
        "authors": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of author names"
        },
        "author_corresponding": {
          "type": [
            "string",
            "null"
          ],
          "description": "Corresponding author name"
        },
        "author_corresponding_institution": {
          "type": [
            "string",
            "null"
          ],
          "description": "Corresponding author's institution"
        },
        "abstract": {
          "type": "string",
          "description": "Full abstract text"
        },
        "date": {
          "type": "string",
          "description": "Publication date (YYYY-MM-DD format)"
        },
        "version": {
          "type": "string",
          "description": "Version number of the preprint"
        },
        "type": {
          "type": "string",
          "description": "Type of submission"
        },
        "license": {
          "type": "string",
          "description": "License type (e.g., 'cc_by', 'cc0')"
        },
        "category": {
          "type": "string",
          "description": "Subject category"
        },
        "published": {
          "type": [
            "string",
            "null"
          ],
          "description": "Published DOI if the preprint has been published in a journal"
        },
        "url": {
          "type": [
            "string",
            "null"
          ],
          "description": "URL to the preprint on medRxiv"
        },
        "pdf_url": {
          "type": [
            "string",
            "null"
          ],
          "description": "Direct PDF download URL"
        },
        "xml_url": {
          "type": [
            "string",
            "null"
          ],
          "description": "JATS XML source URL"
        },
        "server": {
          "type": "string",
          "description": "Server the preprint is hosted on (will be 'medrxiv')"
        }
      }
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "parameters": {
        "query": {
          "type": "string"
        },
        "limit": {
          "type": "integer",
          "minimum": 1,
          "maximum": 100
        }
      },
      "optional": [
        "doi"
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `MedRxiv_get_preprint`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("med_rxiv_get_preprint", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.openalex_get_author` {#dnakit.databases.openalex_get_author}

Get a single OpenAlex author by Author ID (A...). You can pass either the short ID (e.g., "A5001226970") or the full URL (e.g., "https://openalex.org/A5001226970").

```text
dnakit.databases.openalex_get_author(*, author_id, mailto=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `author_id` | string | Required by input schema | OpenAlex Author ID (A...) or full OpenAlex URL. |
| `mailto` | string | None | Optional contact email for OpenAlex polite pool. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "author_id": {
          "type": "string",
          "description": "OpenAlex Author ID (A...) or full OpenAlex URL."
        },
        "mailto": {
          "type": "string",
          "description": "Optional contact email for OpenAlex polite pool."
        }
      },
      "required": [
        "author_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "object",
        "null"
      ],
      "description": "OpenAlex author data.",
      "additionalProperties": true
    }
    ```

**Execution and provenance:** `network`；upstream name `openalex_get_author`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("openalex_get_author", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.openalex_literature_search` {#dnakit.databases.openalex_literature_search}

Search for academic literature using OpenAlex. Supports optional full-text-index filtering (has_fulltext + fulltext.search) when you need body-only term discovery, but note: OpenAlex full-text indexing is incomplete, and a work can be open-access elsewhere (e.g., PMC) while still having has_fulltext=false in OpenAlex. Returns paper metadata and (when available) OA links/content URLs.

```text
dnakit.databases.openalex_literature_search(*, search_keywords=None, query=None, max_results=None, limit=None, num_results=None, year_from=None, start_year=None, year_to=None, end_year=None, open_access=None, require_has_fulltext=None, fulltext_terms=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `search_keywords` | string | None | Keywords to search for in paper titles/abstracts/etc. For full-text-index-only matching, also use require_has_fulltext/fulltext_terms. |
| `query` | string | None | Alias for `search_keywords` (recommended when you standardize on `query` across multiple paper-search tools). |
| `max_results` | integer | None | Maximum number of papers to retrieve (default: 10, maximum: 200). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `limit` | integer | None | Alias for `max_results` (OpenAlex max 200). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `num_results` | integer | None | Alias for `max_results` (OpenAlex max 200). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `year_from` | integer | None | Start year for publication date filter (e.g., 2020). Optional parameter to limit search to papers published from this year onwards. |
| `start_year` | integer | None | Alias for `year_from`. |
| `year_to` | integer | None | End year for publication date filter (e.g., 2023). Optional parameter to limit search to papers published up to this year. |
| `end_year` | integer | None | Alias for `year_to`. |
| `open_access` | boolean | None | Filter for open access papers only. Set to true for open access papers, false for non-open access, or omit for all papers. |
| `require_has_fulltext` | boolean | None | If true, filters to works where OpenAlex indicates a full-text index is available (has_fulltext:true). {&quot;default&quot;: false} |
| `fulltext_terms` | array | None | Optional list of terms that must occur in OpenAlex full-text index (adds fulltext.search:&lt;term&gt; filters and implicitly enables require_has_fulltext). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "search_keywords": {
          "type": "string",
          "description": "Keywords to search for in paper titles/abstracts/etc. For full-text-index-only matching, also use require_has_fulltext/fulltext_terms."
        },
        "query": {
          "type": "string",
          "description": "Alias for `search_keywords` (recommended when you standardize on `query` across multiple paper-search tools)."
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of papers to retrieve (default: 10, maximum: 200).",
          "default": 10,
          "minimum": 1,
          "maximum": 200
        },
        "limit": {
          "type": "integer",
          "description": "Alias for `max_results` (OpenAlex max 200).",
          "minimum": 1,
          "maximum": 200
        },
        "num_results": {
          "type": "integer",
          "description": "Alias for `max_results` (OpenAlex max 200).",
          "minimum": 1,
          "maximum": 200
        },
        "year_from": {
          "type": "integer",
          "description": "Start year for publication date filter (e.g., 2020). Optional parameter to limit search to papers published from this year onwards."
        },
        "start_year": {
          "type": "integer",
          "description": "Alias for `year_from`."
        },
        "year_to": {
          "type": "integer",
          "description": "End year for publication date filter (e.g., 2023). Optional parameter to limit search to papers published up to this year."
        },
        "end_year": {
          "type": "integer",
          "description": "Alias for `year_to`."
        },
        "open_access": {
          "type": "boolean",
          "description": "Filter for open access papers only. Set to true for open access papers, false for non-open access, or omit for all papers."
        },
        "require_has_fulltext": {
          "type": "boolean",
          "description": "If true, filters to works where OpenAlex indicates a full-text index is available (has_fulltext:true).",
          "default": false
        },
        "fulltext_terms": {
          "type": "array",
          "description": "Optional list of terms that must occur in OpenAlex full-text index (adds fulltext.search:<term> filters and implicitly enables require_has_fulltext).",
          "items": {
            "type": "string"
          }
        }
      },
      "anyOf": [
        {
          "required": [
            "search_keywords"
          ]
        },
        {
          "required": [
            "query"
          ]
        }
      ],
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array",
      "description": "List of academic papers from OpenAlex matching the search query",
      "items": {
        "type": "object",
        "properties": {
          "title": {
            "type": "string",
            "description": "Title of the paper"
          },
          "abstract": {
            "type": "string",
            "description": "Abstract of the paper"
          },
          "authors": {
            "type": "array",
            "items": {
              "type": "string"
            },
            "description": "List of author names"
          },
          "year": {
            "type": [
              "integer",
              "string"
            ],
            "description": "Publication year (integer) or 'Year not available' (string)"
          },
          "organizations": {
            "type": "array",
            "items": {
              "type": "string"
            },
            "description": "List of organizational affiliations"
          },
          "venue": {
            "type": "string",
            "description": "Journal or venue name"
          },
          "doi": {
            "type": "string",
            "description": "Digital Object Identifier (DOI) of the paper"
          },
          "citation_count": {
            "type": "integer",
            "description": "Number of citations"
          },
          "open_access": {
            "type": "boolean",
            "description": "Open access status"
          },
          "pdf_url": {
            "type": [
              "string",
              "null"
            ],
            "description": "URL to PDF if available, or null if not available"
          },
          "keywords": {
            "type": [
              "array",
              "string"
            ],
            "description": "List of keywords/concepts associated with the paper, or 'Keywords not available' if none"
          },
          "article_type": {
            "type": "string",
            "description": "Type of article (e.g., 'article', 'review', etc.)"
          },
          "publisher": {
            "type": "string",
            "description": "Publisher name"
          },
          "openalex_id": {
            "type": "string",
            "description": "OpenAlex identifier (URL format)"
          },
          "url": {
            "type": "string",
            "description": "URL to the paper (DOI if available, otherwise OpenAlex ID)"
          },
          "data_quality": {
            "type": "object",
            "description": "Data quality indicators for the paper",
            "properties": {
              "has_abstract": {
                "type": "boolean",
                "description": "Whether the paper has an abstract"
              },
              "has_authors": {
                "type": "boolean",
                "description": "Whether author information is available"
              },
              "has_venue": {
                "type": "boolean",
                "description": "Whether venue/journal information is available"
              },
              "has_year": {
                "type": "boolean",
                "description": "Whether publication year is available"
              },
              "has_doi": {
                "type": "boolean",
                "description": "Whether DOI is available"
              },
              "has_citation_count": {
                "type": "boolean",
                "description": "Whether citation count is available and greater than 0"
              },
              "has_keywords": {
                "type": "boolean",
                "description": "Whether keywords/concepts are available"
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `openalex_literature_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("openalex_literature_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.openalex_search_works` {#dnakit.databases.openalex_search_works}

Search OpenAlex works (papers) via the /works endpoint. Supports both general search and optional full-text-index filtering: set require_has_fulltext=true and/or pass fulltext_terms to add OpenAlex filters (has_fulltext:true and fulltext.search:<term>). Important: has_fulltext/fulltext.search only match works where OpenAlex has indexed full text; this can miss papers whose full text is available elsewhere (e.g., PMC) but not indexed by OpenAlex. Use this to discover Work IDs (W...) and DOIs, then follow up with openalex_get_work/openalex_get_work_by_doi to retrieve OA links and (sometimes) content_urls.

```text
dnakit.databases.openalex_search_works(*, search=None, query=None, filter=None, require_has_fulltext=None, fulltext_terms=None, per_page=None, limit=None, page=None, sort=None, mailto=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `search` | string | None | Search query for works. Use filter + fulltext_terms/require_has_fulltext when you need full-text-index-only matching. |
| `query` | string | None | Alias for `search` (recommended when you standardize on `query` across multiple paper-search tools). |
| `filter` | string | None | OpenAlex filter string (comma-separated). Example: &quot;from_publication_date:2020-01-01,is_oa:true&quot;. |
| `require_has_fulltext` | boolean | None | If true, appends OpenAlex filter has_fulltext:true (keeps only works with full-text index available). {&quot;default&quot;: false} |
| `fulltext_terms` | array | None | Optional list of terms to match in OpenAlex full-text index. Adds one or more fulltext.search:&lt;term&gt; filters and implicitly enables require_has_fulltext. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `per_page` | integer | None | Results per page (OpenAlex max 200). {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `limit` | integer | None | Alias for `per_page` (OpenAlex max 200). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 200} |
| `page` | integer | None | Page number (1-indexed). {&quot;default&quot;: 1, &quot;minimum&quot;: 1} |
| `sort` | string | None | Sort order string, e.g. &quot;cited_by_count:desc&quot;. |
| `mailto` | string | None | Optional contact email for OpenAlex polite pool. If omitted, ToolUniverse uses a default. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "search": {
          "type": "string",
          "description": "Search query for works. Use filter + fulltext_terms/require_has_fulltext when you need full-text-index-only matching."
        },
        "query": {
          "type": "string",
          "description": "Alias for `search` (recommended when you standardize on `query` across multiple paper-search tools)."
        },
        "filter": {
          "type": "string",
          "description": "OpenAlex filter string (comma-separated). Example: \"from_publication_date:2020-01-01,is_oa:true\"."
        },
        "require_has_fulltext": {
          "type": "boolean",
          "description": "If true, appends OpenAlex filter has_fulltext:true (keeps only works with full-text index available).",
          "default": false
        },
        "fulltext_terms": {
          "type": "array",
          "description": "Optional list of terms to match in OpenAlex full-text index. Adds one or more fulltext.search:<term> filters and implicitly enables require_has_fulltext.",
          "items": {
            "type": "string"
          }
        },
        "per_page": {
          "type": "integer",
          "description": "Results per page (OpenAlex max 200).",
          "default": 10,
          "minimum": 1,
          "maximum": 200
        },
        "limit": {
          "type": "integer",
          "description": "Alias for `per_page` (OpenAlex max 200).",
          "minimum": 1,
          "maximum": 200
        },
        "page": {
          "type": "integer",
          "description": "Page number (1-indexed).",
          "default": 1,
          "minimum": 1
        },
        "sort": {
          "type": "string",
          "description": "Sort order string, e.g. \"cited_by_count:desc\"."
        },
        "mailto": {
          "type": "string",
          "description": "Optional contact email for OpenAlex polite pool. If omitted, ToolUniverse uses a default."
        }
      },
      "anyOf": [
        {
          "required": [
            "search"
          ]
        },
        {
          "required": [
            "query"
          ]
        }
      ],
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "object",
        "null"
      ],
      "description": "OpenAlex /works search response data.",
      "additionalProperties": true,
      "properties": {
        "meta": {
          "type": "object",
          "additionalProperties": true
        },
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "additionalProperties": true,
            "properties": {
              "id": {
                "type": "string"
              },
              "title": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "doi": {
                "type": [
                  "string",
                  "null"
                ]
              },
              "publication_year": {
                "type": [
                  "integer",
                  "null"
                ]
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `openalex_search_works`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("openalex_search_works", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_med_get_article` {#dnakit.databases.pub_med_get_article}

Get complete metadata for a specific PubMed article by its PMID (PubMed ID) using efetch. Returns comprehensive article information including title, full abstract, complete author list with affiliations, journal details (name, volume, issue, pages), publication dates, DOI, MeSH terms, article type, grant information, and references. Essential for accessing full abstracts, detailed metadata, and citation information. Find PMIDs using PubMed_search_articles (format: integer like 12345678).

```text
dnakit.databases.pub_med_get_article(*, pmid, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `pmid` | string &#124; integer | Required by input schema | PubMed ID (PMID) as string or integer (e.g., &#x27;12345678&#x27; or 12345678, &#x27;33250470&#x27; or 33250470). Find PMIDs using PubMed_search_articles. Can also accept comma-separated PMIDs as string for multiple articles. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "pmid": {
          "type": [
            "string",
            "integer"
          ],
          "description": "PubMed ID (PMID) as string or integer (e.g., '12345678' or 12345678, '33250470' or 33250470). Find PMIDs using PubMed_search_articles. Can also accept comma-separated PMIDs as string for multiple articles."
        }
      },
      "required": [
        "pmid"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "object",
        "string"
      ],
      "description": "Full PubMed article metadata in XML format (as string) or parsed structure"
    }
    ```

**Execution and provenance:** `network`；upstream name `PubMed_get_article`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_med_get_article", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_med_get_cited_by` {#dnakit.databases.pub_med_get_cited_by}

Get a list of PubMed articles that cite a specific PMID using elink. Returns PMIDs of articles that have cited the input article in their reference list. Use this to track citation impact, find more recent work building on a paper, or explore how research has been applied. Essential for literature reviews, impact analysis, and staying current with research developments. Find PMIDs using PubMed_search_articles.

```text
dnakit.databases.pub_med_get_cited_by(*, pmid, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `pmid` | string | Required by input schema | PubMed ID (PMID) for which to find citing articles (e.g., &#x27;12345678&#x27;). Find PMIDs using PubMed_search_articles. |
| `limit` | integer | None | Maximum number of citing articles to return (default: 20, max: 100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "pmid": {
          "type": "string",
          "description": "PubMed ID (PMID) for which to find citing articles (e.g., '12345678'). Find PMIDs using PubMed_search_articles."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of citing articles to return (default: 20, max: 100).",
          "default": 20,
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": [
        "pmid"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "array",
        "object",
        "string"
      ],
      "description": "List of PMIDs that cite the input article, linkset data structure, or error message"
    }
    ```

**Execution and provenance:** `network`；upstream name `PubMed_get_cited_by`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_med_get_cited_by", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_med_get_related` {#dnakit.databases.pub_med_get_related}

Get related PubMed articles for a specific PMID using elink. Returns a list of computationally similar articles based on title, abstract, and MeSH terms. Use this to discover additional relevant literature, find similar studies, or explore a research topic more comprehensively. The algorithm identifies articles with similar content, methodology, or subject matter. Find PMIDs using PubMed_search_articles.

```text
dnakit.databases.pub_med_get_related(*, pmid, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `pmid` | string | Required by input schema | PubMed ID (PMID) for which to find related articles (e.g., &#x27;20210808&#x27;, &#x27;19879512&#x27;). Find PMIDs using PubMed_search_articles. |
| `limit` | integer | None | Maximum number of related articles to return (default: 20, max: 100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "pmid": {
          "type": "string",
          "description": "PubMed ID (PMID) for which to find related articles (e.g., '20210808', '19879512'). Find PMIDs using PubMed_search_articles."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of related articles to return (default: 20, max: 100).",
          "default": 20,
          "minimum": 1,
          "maximum": 100
        }
      },
      "required": [
        "pmid"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": [
        "array",
        "object",
        "string"
      ],
      "description": "List of related PMIDs, linkset data structure, or error message"
    }
    ```

**Execution and provenance:** `network`；upstream name `PubMed_get_related`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_med_get_related", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.pub_med_search_articles` {#dnakit.databases.pub_med_search_articles}

Search PubMed biomedical literature database using NCBI E-utilities (esearch + esummary). Returns articles with PMID, title, authors, journal, publication year, DOI, article type, and PubMed URL. Automatically fetches rich metadata for each article. Use this to discover medical/life sciences literature, find articles by topic, or search for publications by author or journal. Rate limiting is automatically handled (3 req/sec without API key, 10 req/sec with key set via NCBI_API_KEY environment variable).

```text
dnakit.databases.pub_med_search_articles(*, query, limit=None, mindate=None, maxdate=None, datetype=None, include_abstract=None, sort=None, max_results=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query for PubMed articles. Use keywords, author names, journal names, or MeSH terms. Examples: &#x27;cancer immunotherapy&#x27;, &#x27;Smith J&#91;Author&#93;&#x27;, &#x27;Nature&#91;Journal&#93;&#x27;, &#x27;diabetes&#91;MeSH&#93;&#x27;. Use AND/OR/NOT for complex queries. Note: PubMed ANDs every space-separated keyword by default, so a long, highly specific free-text query (5+ concepts) often returns 0 results even when relevant literature exists. Prefer 2-4 key concepts, add explicit OR for synonyms, or try EuropePMC_search_articles for broader relevance-ranked full-text search if this returns 0 results. |
| `limit` | integer | None | Number of articles to return. This sets the maximum number of articles retrieved from PubMed (max: 200). {&quot;default&quot;: 10} |
| `mindate` | string | None | Minimum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set. |
| `maxdate` | string | None | Maximum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set. |
| `datetype` | string | None | Type of date to filter on: &#x27;pdat&#x27; (publication date), &#x27;edat&#x27; (Entrez date), or &#x27;mdat&#x27; (modification date). Required when using mindate/maxdate. {&quot;default&quot;: &quot;pdat&quot;} |
| `include_abstract` | boolean | None | If true, best-effort fetches abstracts via efetch (adds abstract/abstract_source fields). {&quot;default&quot;: false} |
| `sort` | string | None | Sort order for results. Valid values: &#x27;pub_date&#x27; (newest first), &#x27;Author&#x27; (alphabetical by first author), &#x27;JournalName&#x27; (alphabetical by journal), &#x27;relevance&#x27; (default PubMed ranking). |
| `max_results` | integer | None | Alias for limit — maximum number of results to return |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query for PubMed articles. Use keywords, author names, journal names, or MeSH terms. Examples: 'cancer immunotherapy', 'Smith J[Author]', 'Nature[Journal]', 'diabetes[MeSH]'. Use AND/OR/NOT for complex queries. Note: PubMed ANDs every space-separated keyword by default, so a long, highly specific free-text query (5+ concepts) often returns 0 results even when relevant literature exists. Prefer 2-4 key concepts, add explicit OR for synonyms, or try EuropePMC_search_articles for broader relevance-ranked full-text search if this returns 0 results."
        },
        "limit": {
          "type": "integer",
          "description": "Number of articles to return. This sets the maximum number of articles retrieved from PubMed (max: 200).",
          "default": 10
        },
        "mindate": {
          "type": "string",
          "description": "Minimum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set."
        },
        "maxdate": {
          "type": "string",
          "description": "Maximum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set."
        },
        "datetype": {
          "type": "string",
          "description": "Type of date to filter on: 'pdat' (publication date), 'edat' (Entrez date), or 'mdat' (modification date). Required when using mindate/maxdate.",
          "default": "pdat"
        },
        "include_abstract": {
          "type": "boolean",
          "description": "If true, best-effort fetches abstracts via efetch (adds abstract/abstract_source fields).",
          "default": false
        },
        "sort": {
          "type": "string",
          "description": "Sort order for results. Valid values: 'pub_date' (newest first), 'Author' (alphabetical by first author), 'JournalName' (alphabetical by journal), 'relevance' (default PubMed ranking)."
        },
        "max_results": {
          "type": "integer",
          "description": "Alias for limit — maximum number of results to return"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "Array of PubMed articles with metadata",
          "items": {
            "type": "object",
            "properties": {
              "pmid": {
                "type": "string",
                "description": "PubMed ID (PMID)"
              },
              "title": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Article title"
              },
              "authors": {
                "type": "array",
                "description": "List of author names, truncated to the first 5. Check authors_truncated/author_count before citing this list — an article with exactly 5 authors and an article whose 20 authors were cut to 5 look identical here. Call PubMed_get_article with the pmid for the complete author list with affiliations.",
                "items": {
                  "type": "string"
                }
              },
              "author_count": {
                "type": "integer",
                "description": "Total number of authors on the article, as reported by PubMed — not the length of the truncated `authors` array."
              },
              "authors_truncated": {
                "type": "boolean",
                "description": "True when the article has more than 5 authors and `authors` therefore lists only the first 5. False means `authors` is the complete author list."
              },
              "journal": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Full journal name"
              },
              "pub_date": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Publication date (e.g., '2023 Jan 15')"
              },
              "pub_year": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Publication year"
              },
              "doi": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Digital Object Identifier (DOI)"
              },
              "pmcid": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "PubMed Central ID (if available)"
              },
              "article_type": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Type of article (e.g., 'Journal Article, Review')"
              },
              "url": {
                "type": "string",
                "description": "Direct URL to the article on PubMed"
              },
              "abstract": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Abstract text (only when include_abstract=true and available)"
              },
              "abstract_source": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Source for abstract text (e.g., 'PubMed')"
              },
              "doi_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "DOI resolver URL"
              },
              "pmc_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "PMC full text URL when PMCID is available"
              }
            }
          }
        },
        {
          "type": "array",
          "description": "Stub response when summaries fail (contains PMIDs and URLs only)",
          "items": {
            "type": "object",
            "properties": {
              "pmid": {
                "type": "string"
              },
              "partial": {
                "type": "boolean"
              },
              "warning": {
                "type": "string"
              }
            },
            "required": [
              "partial"
            ]
          }
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `PubMed_search_articles`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("pub_med_search_articles", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.search_clinical_trials` {#dnakit.databases.search_clinical_trials}

Search ClinicalTrials.gov for clinical trials by disease/condition, drug/intervention, or keywords. This is the PRIMARY tool for finding clinical trials — use it when you need to discover trials for a drug (e.g., 'olaparib', 'pembrolizumab'), a disease (e.g., 'lung cancer', 'diabetes'), or any study criteria. Accepts condition names, intervention/treatment names, and free-text queries. Returns NCT IDs and brief descriptions with pagination. Use the returned NCT IDs with 'get_clinical_trial_*' tools for detailed protocol info, or 'extract_clinical_trial_*' tools for efficacy/safety results. Limited to trials beyond phase 1. Query rewriting: multi-word `intervention` values are matched as an EXACT PHRASE against each study's registered intervention names/synonyms, so a phrase no registrant recorded verbatim can match 0 studies even when many trials mention it. The query actually executed is returned in `executed_query` (detail in `query_rewrites`), and a zero-hit rewritten query reports in `relaxed_match_check` how many studies your original wording would have matched. Wrap a value in parentheses (e.g. '(cell therapy)') to send it verbatim.

```text
dnakit.databases.search_clinical_trials(*, condition=None, intervention=None, query_term=None, pageSize=None, pageToken=None, overall_status=None, max_results=None, limit=None, status=None, keyword=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `condition` | string | None | Query for condition or disease using Essie expression syntax (e.g., &#x27;lung cancer&#x27;, &#x27;(head OR neck) AND pain AND NOT &quot;back pain&quot;&#x27;). Sent to ClinicalTrials.gov verbatim (this tool does not rewrite condition values); see `executed_query` in the response for what actually ran. |
| `intervention` | string | None | Query for intervention/treatment using Essie expression syntax (e.g., &#x27;chemotherapy&#x27;, &#x27;immunotherapy&#x27;, &#x27;olaparib&#x27;, &#x27;combination therapy&#x27;). Matched against each study&#x27;s registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, so prefer the generic/INN name: a brand name that no registrant recorded in either field can return zero results. Values that already contain Essie syntax (AREA&#91;...&#93;/AND/OR/NOT/quotes) are passed through unchanged. Matched against each study&#x27;s registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so &#x27;donor derived cell therapy&#x27; matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. &#x27;(donor derived cell therapy)&#x27;) to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response. |
| `query_term` | string | None | Query for &#x27;other terms&#x27; with Essie expression syntax (e.g., &#x27;combination&#x27;, &#x27;AREA&#91;LastUpdatePostDate&#93;RANGE&#91;2023-01-15,MAX&#93;&#x27;, &#x27;Phase II&#x27;). Can be used to search for all other protocol fields, including but not limited to title, outcome measures, status, phase, location, etc. |
| `pageSize` | integer | None | Maximum number of studies to return per page (default 10, max 1000). |
| `pageToken` | string | None | Token to retrieve the next page of results, obtained from the &#x27;nextPageToken&#x27; field of the previous response. Do not specify it for first page. When you make an initial request to the API which supports pagination, the response will include a nextPageToken. This token can then be used as a parameter in the subsequent API request to retrieve the next set of data. |
| `overall_status` | array | None | Filter by overall study status (e.g., &#91;&#x27;RECRUITING&#x27;&#93;, &#91;&#x27;COMPLETED&#x27;&#93;, &#91;&#x27;RECRUITING&#x27;, &#x27;NOT_YET_RECRUITING&#x27;&#93;). Valid values: RECRUITING, NOT_YET_RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, ENROLLING_BY_INVITATION, SUSPENDED, TERMINATED, WITHDRAWN. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `max_results` | integer | None | Maximum number of studies to return (alias for pageSize, default 10, max 1000). |
| `limit` | integer | None | Alias for max_results: maximum number of studies to return (default 10, max 1000). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `status` | object &#124; None | None | Alias for overall_status. Filter by trial status, e.g. &quot;RECRUITING&quot;, &quot;COMPLETED&quot;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;string&quot;}, {&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}&#93;} |
| `keyword` | string | None | Alias for query_term. Free-text keyword search across all trial fields (e.g., drug name, condition, investigator). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "condition": {
          "type": "string",
          "description": "Query for condition or disease using Essie expression syntax (e.g., 'lung cancer', '(head OR neck) AND pain AND NOT \"back pain\"'). Sent to ClinicalTrials.gov verbatim (this tool does not rewrite condition values); see `executed_query` in the response for what actually ran."
        },
        "intervention": {
          "type": "string",
          "description": "Query for intervention/treatment using Essie expression syntax (e.g., 'chemotherapy', 'immunotherapy', 'olaparib', 'combination therapy'). Matched against each study's registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, so prefer the generic/INN name: a brand name that no registrant recorded in either field can return zero results. Values that already contain Essie syntax (AREA[...]/AND/OR/NOT/quotes) are passed through unchanged. Matched against each study's registered intervention names and synonyms (InterventionName/InterventionOtherName), not free text in study summaries, and multi-word values are matched as an EXACT PHRASE -- so 'donor derived cell therapy' matches 0 studies even though 251 trials match the same words loosely. For looser matching, wrap the value in parentheses (e.g. '(donor derived cell therapy)') to send it verbatim, or use ClinicalTrials_search_studies with `query_term` for free text across all study fields. Prefer the generic/INN name over a brand name. See `executed_query` and `relaxed_match_check` in the response."
        },
        "query_term": {
          "type": "string",
          "description": "Query for 'other terms' with Essie expression syntax (e.g., 'combination', 'AREA[LastUpdatePostDate]RANGE[2023-01-15,MAX]', 'Phase II'). Can be used to search for all other protocol fields, including but not limited to title, outcome measures, status, phase, location, etc."
        },
        "pageSize": {
          "type": "integer",
          "description": "Maximum number of studies to return per page (default 10, max 1000)."
        },
        "pageToken": {
          "type": "string",
          "description": "Token to retrieve the next page of results, obtained from the 'nextPageToken' field of the previous response. Do not specify it for first page. When you make an initial request to the API which supports pagination, the response will include a nextPageToken. This token can then be used as a parameter in the subsequent API request to retrieve the next set of data."
        },
        "overall_status": {
          "type": "array",
          "description": "Filter by overall study status (e.g., ['RECRUITING'], ['COMPLETED'], ['RECRUITING', 'NOT_YET_RECRUITING']). Valid values: RECRUITING, NOT_YET_RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, ENROLLING_BY_INVITATION, SUSPENDED, TERMINATED, WITHDRAWN.",
          "items": {
            "type": "string"
          }
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of studies to return (alias for pageSize, default 10, max 1000)."
        },
        "limit": {
          "type": "integer",
          "description": "Alias for max_results: maximum number of studies to return (default 10, max 1000).",
          "minimum": 1,
          "maximum": 1000
        },
        "status": {
          "oneOf": [
            {
              "type": "string"
            },
            {
              "type": "array",
              "items": {
                "type": "string"
              }
            }
          ],
          "description": "Alias for overall_status. Filter by trial status, e.g. \"RECRUITING\", \"COMPLETED\"."
        },
        "keyword": {
          "type": "string",
          "description": "Alias for query_term. Free-text keyword search across all trial fields (e.g., drug name, condition, investigator)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "string",
          "description": "Message when no studies are found"
        },
        {
          "type": "object",
          "properties": {
            "studies": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "NCT ID": {
                    "type": "string"
                  },
                  "brief_title": {
                    "type": "string"
                  },
                  "brief_summary": {
                    "type": "string"
                  },
                  "overall_status": {
                    "type": "string"
                  },
                  "condition": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  },
                  "phase": {
                    "type": "array",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "nextPageToken": {
              "type": "string"
            },
            "total_count": {
              "type": "integer"
            },
            "executed_query": {
              "type": "object",
              "description": "The query parameters actually sent to ClinicalTrials.gov, after phrase-quoting and AREA[] field restriction. Compare against what you passed: they differ whenever `query_rewrites` is non-empty."
            },
            "query_rewrites": {
              "type": "array",
              "description": "One entry per parameter whose value was rewritten before being sent. Empty when the query ran exactly as supplied.",
              "items": {
                "type": "object",
                "properties": {
                  "parameter": {
                    "type": "string",
                    "description": "The parameter name you passed."
                  },
                  "api_field": {
                    "type": "string",
                    "description": "The ClinicalTrials.gov query field it was sent as."
                  },
                  "submitted": {
                    "type": "string",
                    "description": "The value you supplied."
                  },
                  "executed": {
                    "type": "string",
                    "description": "The rewritten value that was actually sent."
                  },
                  "restricted_to_fields": {
                    "type": "array",
                    "description": "Registry fields the value was restricted to, if any.",
                    "items": {
                      "type": "string"
                    }
                  }
                }
              }
            },
            "relaxed_match_check": {
              "type": "object",
              "description": "Present ONLY when a rewritten query returned zero studies. Reports whether your original wording would have matched, so an empty `studies` list can be told apart from a genuine absence of trials. The relaxed results are never substituted for the ones you asked for -- follow `note` to retrieve them.",
              "properties": {
                "relaxed_query": {
                  "type": "object",
                  "description": "The un-rewritten query used for the check."
                },
                "relaxed_total_count": {
                  "type": [
                    "integer",
                    "null"
                  ],
                  "description": "Studies matching your original wording; null if the check request failed."
                },
                "note": {
                  "type": "string",
                  "description": "Plain-language reading of the check, naming the escape hatch when the zero is an artifact of the rewrite."
                }
              }
            }
          },
          "required": [
            "studies"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `search_clinical_trials`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("search_clinical_trials", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.semantic_scholar_search_papers` {#dnakit.databases.semantic_scholar_search_papers}

Search for papers on Semantic Scholar including abstracts and AI-generated TLDR summaries. This tool queries the Semantic Scholar API using natural language keywords and returns papers with details such as title, abstract, TLDR summary, fields of study, publication year, journal (venue), and URL. Returns only one page of results: `metadata.returned` is how many papers are in `data`, `metadata.total` is Semantic Scholar's approximate count of ALL papers matching the query in the corpus (null when the API does not report it, see `metadata.total_source`), and the top-level `truncated` flag plus `truncation_note` say when more matches exist than were returned. Do not treat `metadata.returned` as the size of the literature. Rate limits: 1 req/sec without API key, 100 req/sec with SEMANTIC_SCHOLAR_API_KEY environment variable set.

```text
dnakit.databases.semantic_scholar_search_papers(*, query, limit=None, year=None, sort=None, include_abstract=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query for Semantic Scholar. Use keywords separated by spaces to refine the search. |
| `limit` | integer | None | Maximum number of papers to return from Semantic Scholar. {&quot;default&quot;: 5} |
| `year` | string | None | Filter results by publication year. Use a single year (e.g., &#x27;2024&#x27;) or a range (e.g., &#x27;2020-2024&#x27;). |
| `sort` | string | None | Sort results. Options: &#x27;citationCount:desc&#x27;, &#x27;citationCount:asc&#x27;, &#x27;publicationDate:desc&#x27;, &#x27;publicationDate:asc&#x27;. |
| `include_abstract` | boolean | None | If true, best-effort fetches missing abstracts via the paper detail endpoint (only when abstract is missing in search results). {&quot;default&quot;: false} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query for Semantic Scholar. Use keywords separated by spaces to refine the search."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of papers to return from Semantic Scholar.",
          "default": 5
        },
        "year": {
          "type": "string",
          "description": "Filter results by publication year. Use a single year (e.g., '2024') or a range (e.g., '2020-2024')."
        },
        "sort": {
          "type": "string",
          "description": "Sort results. Options: 'citationCount:desc', 'citationCount:asc', 'publicationDate:desc', 'publicationDate:asc'."
        },
        "include_abstract": {
          "type": "boolean",
          "description": "If true, best-effort fetches missing abstracts via the paper detail endpoint (only when abstract is missing in search results).",
          "default": false
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "description": "List of papers from Semantic Scholar matching the search query",
          "items": {
            "type": "object",
            "properties": {
              "title": {
                "type": "string",
                "description": "Title of the paper"
              },
              "abstract": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Abstract of the paper"
              },
              "tldr": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "AI-generated one-sentence TLDR summary of the paper"
              },
              "journal": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Journal or venue name"
              },
              "year": {
                "type": [
                  "integer",
                  "null"
                ],
                "description": "Publication year"
              },
              "url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "URL to the paper on Semantic Scholar"
              },
              "paper_id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Semantic Scholar paperId"
              },
              "doi": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "DOI (if available in externalIds)"
              },
              "doi_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "DOI resolver URL"
              },
              "authors": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "List of author names"
              },
              "fields_of_study": {
                "type": [
                  "array",
                  "null"
                ],
                "description": "Fields of study (e.g., Biology, Medicine, Computer Science)",
                "items": {
                  "type": "string"
                }
              },
              "citation_count": {
                "type": [
                  "integer",
                  "null"
                ],
                "description": "Citation count (if provided)"
              },
              "reference_count": {
                "type": [
                  "integer",
                  "null"
                ],
                "description": "Reference count (if provided)"
              },
              "open_access": {
                "type": [
                  "boolean",
                  "null"
                ],
                "description": "Open access flag (if provided)"
              },
              "open_access_pdf_url": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Open access PDF URL (if provided)"
              }
            }
          }
        },
        {
          "type": "array",
          "description": "Error response (single-element array with error details)",
          "items": {
            "type": "object",
            "properties": {
              "title": {
                "type": "string",
                "description": "Always \"Error\" for error responses"
              },
              "error": {
                "type": "string",
                "description": "Error message"
              },
              "retryable": {
                "type": "boolean",
                "description": "Whether the request can be retried"
              },
              "suggestion": {
                "type": "string",
                "description": "Suggested next step"
              }
            },
            "required": [
              "error"
            ]
          }
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `SemanticScholar_search_papers`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("semantic_scholar_search_papers", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.web_search` {#dnakit.databases.web_search}

General web search using DDGS (Dux Distributed Global Search) supporting multiple search engines including Google, Bing, Brave, Yahoo, and DuckDuckGo. No API keys required. Parallel Search MCP and SerpBase are opt-in; selecting either sends the query to that hosted third-party service, so do not submit sensitive or patient-identifying information. SerpBase additionally requires SERPBASE_API_KEY and returns real Google SERP results (billed per request), useful when DDGS's Google backend is unavailable or rate-limited.

```text
dnakit.databases.web_search(*, query, max_results=None, search_type=None, backend=None, region=None, safesearch=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query string |
| `max_results` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `search_type` | string | None | Type of search to perform {&quot;enum&quot;: &#91;&quot;general&quot;, &quot;api_documentation&quot;, &quot;python_packages&quot;, &quot;github&quot;&#93;, &quot;default&quot;: &quot;general&quot;} |
| `backend` | string | None | Search engine backend to use. Selecting &#x27;parallel&#x27; sends the query to https://search.parallel.ai/mcp; that backend does not apply ToolUniverse region or safesearch controls. Selecting &#x27;serpbase&#x27; sends the query to https://api.serpbase.dev (requires SERPBASE_API_KEY) for real Google SERP results; that backend applies region but not safesearch. {&quot;enum&quot;: &#91;&quot;auto&quot;, &quot;parallel&quot;, &quot;serpbase&quot;, &quot;google&quot;, &quot;bing&quot;, &quot;brave&quot;, &quot;yahoo&quot;, &quot;duckduckgo&quot;&#93;, &quot;default&quot;: &quot;auto&quot;} |
| `region` | string | None | Search region/locale for DDGS backends, and country/language for the SerpBase backend (split on the hyphen into gl/hl). Parallel does not support this control; only the default value is accepted with backend=&#x27;parallel&#x27;. {&quot;enum&quot;: &#91;&quot;us-en&quot;, &quot;cn-zh&quot;, &quot;uk-en&quot;, &quot;de-de&quot;, &quot;fr-fr&quot;, &quot;ja-jp&quot;&#93;, &quot;default&quot;: &quot;us-en&quot;} |
| `safesearch` | string | None | Safe search level for DDGS backends. Parallel and SerpBase do not support this control; only the default value is accepted with those backends. {&quot;enum&quot;: &#91;&quot;on&quot;, &quot;moderate&quot;, &quot;off&quot;&#93;, &quot;default&quot;: &quot;moderate&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query string"
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of results to return",
          "default": 10,
          "minimum": 1,
          "maximum": 50
        },
        "search_type": {
          "type": "string",
          "description": "Type of search to perform",
          "enum": [
            "general",
            "api_documentation",
            "python_packages",
            "github"
          ],
          "default": "general"
        },
        "backend": {
          "type": "string",
          "description": "Search engine backend to use. Selecting 'parallel' sends the query to https://search.parallel.ai/mcp; that backend does not apply ToolUniverse region or safesearch controls. Selecting 'serpbase' sends the query to https://api.serpbase.dev (requires SERPBASE_API_KEY) for real Google SERP results; that backend applies region but not safesearch.",
          "enum": [
            "auto",
            "parallel",
            "serpbase",
            "google",
            "bing",
            "brave",
            "yahoo",
            "duckduckgo"
          ],
          "default": "auto"
        },
        "region": {
          "type": "string",
          "description": "Search region/locale for DDGS backends, and country/language for the SerpBase backend (split on the hyphen into gl/hl). Parallel does not support this control; only the default value is accepted with backend='parallel'.",
          "enum": [
            "us-en",
            "cn-zh",
            "uk-en",
            "de-de",
            "fr-fr",
            "ja-jp"
          ],
          "default": "us-en"
        },
        "safesearch": {
          "type": "string",
          "description": "Safe search level for DDGS backends. Parallel and SerpBase do not support this control; only the default value is accepted with those backends.",
          "enum": [
            "on",
            "moderate",
            "off"
          ],
          "default": "moderate"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "status": {
              "type": "string"
            },
            "query": {
              "type": "string"
            },
            "search_type": {
              "type": "string"
            },
            "total_results": {
              "type": [
                "integer",
                "number"
              ]
            },
            "results": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "title": {
                    "type": "string"
                  },
                  "url": {
                    "type": "string"
                  },
                  "snippet": {
                    "type": "string"
                  },
                  "rank": {
                    "type": [
                      "integer",
                      "number"
                    ]
                  }
                }
              }
            },
            "backend_used": {
              "type": "string"
            },
            "attempted_backends": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "provider_errors": {
              "type": "object",
              "properties": {
                "duckduckgo": {
                  "type": "string"
                }
              }
            },
            "provider_notice": {
              "type": "string",
              "description": "Disclosure included when the query was sent to a hosted third-party backend."
            }
          },
          "required": [
            "status"
          ]
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `external_process`；upstream name `web_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("web_search", verbose=False)
print(info["input_schema"])
```

