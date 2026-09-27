# 文献、预印本与临床研究查询

查询文献、预印本和临床研究记录。

<!-- BEGIN GENERATED API: literature -->

<span id="literature"></span>

## 1) 按关键词检索 arXiv 论文 {#dnakit.databases.ar_xiv_search_papers}

- **作用：** 按关键词检索 arXiv 论文。
- **处理方法：** 通过 `ArXiv` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ar_xiv_search_papers(*, query, limit=None, sort_by=None, sort_order=None, date_from=None, date_to=None, category=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `title`, `abstract`, `authors`, `published` 等字段。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.ar_xiv_search_papers-parameters)

## 2) 按 DOI 查询 bioRxiv 预印本元数据 {#dnakit.databases.bio_rxiv_get_preprint}

- **作用：** 按 DOI 查询 bioRxiv 预印本元数据。
- **处理方法：** 通过 `BioRxiv` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.bio_rxiv_get_preprint(*, doi, server=None, verbose=True)`
- **输入：** 必须提供 `doi`。
- **输出：** 返回字典，包含 `doi`, `title`, `authors`, `author_corresponding` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.bio_rxiv_get_preprint(
    doi='10.1101/2023.12.01.569554',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.bio_rxiv_get_preprint-parameters)

## 3) 按日期范围列出近期预印本 {#dnakit.databases.bio_rxiv_list_recent_preprints}

- **作用：** 按日期范围列出近期预印本。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.bio_rxiv_list_recent_preprints(*, start_date=None, end_date=None, server=None, cursor=None, query=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `messages`, `collection` 等字段。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.bio_rxiv_list_recent_preprints-parameters)

## 4) 按疾病、干预或申办方检索临床研究 {#dnakit.databases.clinical_trials_search_studies}

- **作用：** 按疾病、干预或申办方检索临床研究。
- **处理方法：** 通过 `ClinicalTrials` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clinical_trials_search_studies(*, operation=None, query_cond=None, query_intr=None, query_term=None, filter_status=None, filter_phase=None, filter_study_type=None, page_size=None, next_page_token=None, query=None, condition=None, status=None, max_results=None, intervention=None, sponsor=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `studies`, `total_count`, `total_count_note`, `next_page_token` 等字段。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.clinical_trials_search_studies-parameters)

## 5) 检索 Crossref 学术文献与 DOI 元数据 {#dnakit.databases.crossref_search_works}

- **作用：** 检索 Crossref 学术文献与 DOI 元数据。
- **处理方法：** 通过 `Crossref` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.crossref_search_works(*, query, limit=None, offset=None, filter=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `DOI`, `title`, `author`, `container-title` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.crossref_search_works(
    query='artificial intelligence',
    limit=2,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.crossref_search_works-parameters)

## 6) 检索 Europe PMC 文献与摘要 {#dnakit.databases.europe_pmc_search_articles}

- **作用：** 检索 Europe PMC 文献与摘要。
- **处理方法：** 通过 `EuropePMC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.europe_pmc_search_articles(*, query, limit=None, require_has_ft=None, fulltext_terms=None, enrich_missing_abstract=None, extract_terms_from_fulltext=None, page_size=None, pageSize=None, source=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `title`, `abstract`, `abstract_retrieval_trace`, `authors` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.europe_pmc_search_articles(
    query='cancer therapy',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.europe_pmc_search_articles-parameters)

## 7) 按 DOI 查询 medRxiv 预印本元数据 {#dnakit.databases.med_rxiv_get_preprint}

- **作用：** 按 DOI 查询 medRxiv 预印本元数据。
- **处理方法：** 通过 `BioRxiv` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.med_rxiv_get_preprint(*, doi=None, server=None, query=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `doi`, `title`, `authors`, `author_corresponding` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.med_rxiv_get_preprint(
    doi='10.1101/2021.04.29.21256344',
    server='medrxiv',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.med_rxiv_get_preprint-parameters)

## 8) 按 OpenAlex 编号查询作者详情 {#dnakit.databases.openalex_get_author}

- **作用：** 按 OpenAlex 编号查询作者详情。
- **处理方法：** 通过 `OpenAlex` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.openalex_get_author(*, author_id, mailto=None, verbose=True)`
- **输入：** 必须提供 `author_id`。
- **输出：** 返回 OpenAlex 作者记录；接口允许返回 null。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.openalex_get_author(
    author_id='A5001226970',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 OpenAlex 作者记录；接口允许返回 null。

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_get_author-parameters)

## 9) 使用 OpenAlex 检索学术文献 {#dnakit.databases.openalex_literature_search}

- **作用：** 使用 OpenAlex 检索学术文献。
- **处理方法：** 通过 `OpenAlex` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.openalex_literature_search(*, search_keywords=None, query=None, max_results=None, limit=None, num_results=None, year_from=None, start_year=None, year_to=None, end_year=None, open_access=None, require_has_fulltext=None, fulltext_terms=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组至少选一组：`search_keywords` / `query`。
- **输出：** 返回记录列表；每条记录包含 `title`, `abstract`, `authors`, `year` 等字段。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_literature_search-parameters)

## 10) 按筛选条件检索 OpenAlex 论文记录 {#dnakit.databases.openalex_search_works}

- **作用：** 按筛选条件检索 OpenAlex 论文记录。
- **处理方法：** 通过 `OpenAlex` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.openalex_search_works(*, search=None, query=None, filter=None, require_has_fulltext=None, fulltext_terms=None, per_page=None, limit=None, page=None, sort=None, mailto=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组至少选一组：`search` / `query`。
- **输出：** 返回字典，包含 `meta`, `results` 等字段。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.openalex_search_works-parameters)

## 11) 按 PMID 获取 PubMed 文献完整元数据 {#dnakit.databases.pub_med_get_article}

- **作用：** 按 PMID 获取 PubMed 文献完整元数据。
- **处理方法：** 通过 `PubMed` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_med_get_article(*, pmid, verbose=True)`
- **输入：** 必须提供 `pmid`。
- **输出：** 返回 PMID 对应文章的 XML 文本或已解析元数据，包含标题、摘要、作者和期刊信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_med_get_article(
    pmid='17284678',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 PMID 对应文章的 XML 文本或已解析元数据，包含标题、摘要、作者和期刊信息。

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_article-parameters)

## 12) 查询引用指定 PubMed 文献的论文 {#dnakit.databases.pub_med_get_cited_by}

- **作用：** 查询引用指定 PubMed 文献的论文。
- **处理方法：** 通过 `PubMed` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_med_get_cited_by(*, pmid, limit=None, verbose=True)`
- **输入：** 必须提供 `pmid`。
- **输出：** 返回“查询引用指定 PubMed 文献的论文”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_med_get_cited_by(
    pmid='20210808',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_cited_by-parameters)

## 13) 查找与指定 PubMed 文献相似的论文 {#dnakit.databases.pub_med_get_related}

- **作用：** 查找与指定 PubMed 文献相似的论文。
- **处理方法：** 通过 `PubMed` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_med_get_related(*, pmid, limit=None, verbose=True)`
- **输入：** 必须提供 `pmid`。
- **输出：** 返回“查找与指定 PubMed 文献相似的论文”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_med_get_related(
    pmid='20210808',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_get_related-parameters)

## 14) 按检索式查询 PubMed 生物医学文献 {#dnakit.databases.pub_med_search_articles}

- **作用：** 按检索式查询 PubMed 生物医学文献。
- **处理方法：** 通过 `PubMed` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_med_search_articles(*, query, limit=None, mindate=None, maxdate=None, datetype=None, include_abstract=None, sort=None, max_results=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `pmid`, `title`, `authors`, `author_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_med_search_articles(
    query='cancer research',
    limit=2,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.pub_med_search_articles-parameters)

## 15) 按疾病、药物或关键词查找临床试验 {#dnakit.databases.search_clinical_trials}

- **作用：** 按疾病、药物或关键词查找临床试验。
- **处理方法：** 通过 `ClinicalTrialsSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.search_clinical_trials(*, condition=None, intervention=None, query_term=None, pageSize=None, pageToken=None, overall_status=None, max_results=None, limit=None, status=None, keyword=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字符串。
- **示例代码：**

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

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
"<string>"
```

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.search_clinical_trials-parameters)

## 16) 检索论文、摘要与简短总结 {#dnakit.databases.semantic_scholar_search_papers}

- **作用：** 检索论文、摘要与简短总结。
- **处理方法：** 通过 `SemanticScholar` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.semantic_scholar_search_papers(*, query, limit=None, year=None, sort=None, include_abstract=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `title`, `abstract`, `tldr`, `journal` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.semantic_scholar_search_papers(
    query='machine learning optimization',
    limit=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.semantic_scholar_search_papers-parameters)

## 17) 通过搜索引擎检索网页 {#dnakit.databases.web_search}

- **作用：** 通过搜索引擎检索网页。
- **处理方法：** 通过 `WebSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.web_search(*, query, max_results=None, search_type=None, backend=None, region=None, safesearch=None, verbose=True)`
- **输入：** 必须提供 `query`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 返回字典，包含 `status`, `query`, `search_type`, `total_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.web_search(
    query='NCBI genome assembly documentation',
    max_results=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

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

[完整参数与使用说明](../../faq/api/22_database_query_literature.md#dnakit.databases.web_search-parameters)

<!-- END GENERATED API: literature -->
