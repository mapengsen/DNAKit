# 文献、预印本与临床研究查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_literature.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 按关键词检索 arXiv 论文 {#dnakit.databases.ar_xiv_search_papers-parameters}

`dnakit.databases.ar_xiv_search_papers(*, query, limit=None, sort_by=None, sort_order=None, date_from=None, date_to=None, category=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for arXiv papers. Unquoted, unprefixed space-separated keywords are ANDed together (every keyword must match) -- adding more keywords NARROWS results, it does not broaden them. To search multiple synonyms/phrases and broaden recall, join them with explicit OR, e.g. &#x27;&quot;machine learning&quot; OR &quot;deep learning&quot; OR &quot;neural networks&quot;&#x27;. Supports field prefixes (au:, ti:, abs:, cat:, all:) and AND/OR/ANDNOT boolean operators, which are passed through as-is. |
| `limit` | integer | None | Number of papers to return. This sets the maximum number of papers retrieved from arXiv. {&quot;default&quot;: 10} |
| `sort_by` | string | None | Sort order for results. Options: &#x27;relevance&#x27;, &#x27;lastUpdatedDate&#x27;, &#x27;submittedDate&#x27; {&quot;default&quot;: &quot;relevance&quot;} |
| `sort_order` | string | None | Sort direction. Options: &#x27;ascending&#x27;, &#x27;descending&#x27; {&quot;default&quot;: &quot;descending&quot;} |
| `date_from` | string | None | Filter results from this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper&#x27;s FIRST version -- a paper revised/updated after this date but originally submitted earlier is still excluded. |
| `date_to` | string | None | Filter results up to this date (format: YYYY-MM-DD). Uses submittedDate, the date of the paper&#x27;s FIRST version -- a paper revised/updated after this date but originally submitted earlier is still included. |
| `category` | string | None | 分类或结果类别，按当前接口的枚举/来源定义填写。 |
| `verbose` | bool | True | 运行提示开关。 |


## 按 DOI 查询 bioRxiv 预印本元数据 {#dnakit.databases.bio_rxiv_get_preprint-parameters}

`dnakit.databases.bio_rxiv_get_preprint(*, doi, server=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `doi` | string | 必填 | bioRxiv or medRxiv DOI. Can be full DOI (e.g., &#x27;10.1101/2023.12.01.569554&#x27;) or just the numeric part after &#x27;10.1101/&#x27; (e.g., &#x27;2023.12.01.569554&#x27;). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations. |
| `server` | string | None | Server to query: &#x27;biorxiv&#x27; for bioRxiv preprints or &#x27;medrxiv&#x27; for medRxiv preprints. Default is &#x27;biorxiv&#x27;. {&quot;default&quot;: &quot;biorxiv&quot;, &quot;enum&quot;: &#91;&quot;biorxiv&quot;, &quot;medrxiv&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按日期范围列出近期预印本 {#dnakit.databases.bio_rxiv_list_recent_preprints-parameters}

`dnakit.databases.bio_rxiv_list_recent_preprints(*, start_date=None, end_date=None, server=None, cursor=None, query=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `start_date` | string | None | Start date in YYYY-MM-DD format (e.g., &#x27;2024-01-01&#x27;). Date range must not exceed 60 days. |
| `end_date` | string | None | End date in YYYY-MM-DD format (e.g., &#x27;2024-01-03&#x27;). Date range must not exceed 60 days. |
| `server` | string | None | Server: &#x27;biorxiv&#x27; for biology preprints, &#x27;medrxiv&#x27; for health sciences preprints {&quot;default&quot;: &quot;biorxiv&quot;, &quot;enum&quot;: &#91;&quot;biorxiv&quot;, &quot;medrxiv&quot;&#93;} |
| `cursor` | integer | None | Pagination cursor (0 for first 100 results, 100 for next 100, etc.) {&quot;default&quot;: 0} |
| `query` | string | None | 检索词或查询表达式；语法由当前数据源定义。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按疾病、干预或申办方检索临床研究 {#dnakit.databases.clinical_trials_search_studies-parameters}

`dnakit.databases.clinical_trials_search_studies(*, operation=None, query_cond=None, query_intr=None, query_term=None, filter_status=None, filter_phase=None, filter_study_type=None, page_size=None, next_page_token=None, query=None, condition=None, status=None, max_results=None, intervention=None, sponsor=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
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
| `verbose` | bool | True | 运行提示开关。 |


## 检索 Crossref 学术文献与 DOI 元数据 {#dnakit.databases.crossref_search_works-parameters}

`dnakit.databases.crossref_search_works(*, query, limit=None, offset=None, filter=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for Crossref works. Use keywords separated by spaces to refine your search across titles, abstracts, authors, and other bibliographic fields. |
| `limit` | integer | None | Number of articles to return in this page. This sets the maximum number of articles retrieved from Crossref. Max 100 per request. It does not change &#x27;total_count&#x27;, which always reports how many works match the query. {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `offset` | integer | None | Skip this many matching works before returning results. Use with &#x27;limit&#x27; to page through a result set larger than 100; the response &#x27;truncation_note&#x27; gives the offset for the next page. {&quot;minimum&quot;: 0} |
| `filter` | string | None | Optional filter string for Crossref API. Examples: &#x27;type:journal-article&#x27; (only journal articles), &#x27;from-pub-date:2020-01-01&#x27; (published after date), &#x27;has-abstract:true&#x27; (only with abstracts). Multiple filters comma-separated. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 Europe PMC 文献与摘要 {#dnakit.databases.europe_pmc_search_articles-parameters}

`dnakit.databases.europe_pmc_search_articles(*, query, limit=None, require_has_ft=None, fulltext_terms=None, enrich_missing_abstract=None, extract_terms_from_fulltext=None, page_size=None, pageSize=None, source=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for Europe PMC. Supports Lucene-like fielded syntax (e.g., BODY:&quot;term&quot;, INTRO:&quot;term&quot;, METHODS:&quot;term&quot;). |
| `limit` | integer | None | Number of articles to return. This sets the maximum number of articles retrieved from Europe PMC; anything beyond it is silently left behind, so compare it against `metadata.total_results` in the response (Europe PMC&#x27;s hitCount) before treating the returned set as complete. Europe PMC accepts up to 1000 per request. {&quot;default&quot;: 5} |
| `require_has_ft` | boolean | None | If true, appends `HAS_FT:Y` to the query to restrict results to records where Europe PMC has full text indexed. This is NOT the same as having a link to free full text elsewhere (e.g., PMC) and is often the key reason body-only keywords cannot be found via Europe PMC search. {&quot;default&quot;: false} |
| `fulltext_terms` | array | None | Optional list of terms that must occur in the indexed full text (Europe PMC BODY field). When provided, the tool adds a BODY:&quot;...&quot; clause and automatically enables `require_has_ft`. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enrich_missing_abstract` | boolean | None | If true, best-effort fills missing abstracts by fetching Europe PMC fullTextXML (bounded to a few results). {&quot;default&quot;: false} |
| `extract_terms_from_fulltext` | array | None | Optional list of terms to extract from full text (open access only). When provided, automatically fetches fullTextXML for open-access articles and returns snippets around these terms. Terms are processed in batches of 5 internally, so any number of terms is accepted. Bounded to first 3 OA articles to avoid latency. Returns snippets in a &#x27;fulltext_snippets&#x27; field for each article. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `page_size` | integer | None | Alias for limit — number of results to return |
| `pageSize` | integer | None | DNAKit optional alias for limit |
| `source` | string | None | 数据来源筛选值；按当前接口允许的来源填写。 |
| `verbose` | bool | True | 运行提示开关。 |


## 按 DOI 查询 medRxiv 预印本元数据 {#dnakit.databases.med_rxiv_get_preprint-parameters}

`dnakit.databases.med_rxiv_get_preprint(*, doi=None, server=None, query=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `doi` | string | None | medRxiv DOI. Can be full DOI (e.g., &#x27;10.1101/2021.04.29.21256344&#x27;) or just the numeric part after &#x27;10.1101/&#x27; (e.g., &#x27;2021.04.29.21256344&#x27;). Find DOIs using EuropePMC_search_articles, web_search, or from paper citations. |
| `server` | string | None | Server to query - always &#x27;medrxiv&#x27; for this tool. {&quot;default&quot;: &quot;medrxiv&quot;, &quot;enum&quot;: &#91;&quot;medrxiv&quot;&#93;} |
| `query` | string | None | 检索词或查询表达式；语法由当前数据源定义。 |
| `limit` | integer | None | 当前请求的返回数量上限；不是全库命中总数。 {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 OpenAlex 编号查询作者详情 {#dnakit.databases.openalex_get_author-parameters}

`dnakit.databases.openalex_get_author(*, author_id, mailto=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `author_id` | string | 必填 | OpenAlex Author ID (A...) or full OpenAlex URL. |
| `mailto` | string | None | Optional contact email for OpenAlex polite pool. |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 OpenAlex 检索学术文献 {#dnakit.databases.openalex_literature_search-parameters}

`dnakit.databases.openalex_literature_search(*, search_keywords=None, query=None, max_results=None, limit=None, num_results=None, year_from=None, start_year=None, year_to=None, end_year=None, open_access=None, require_has_fulltext=None, fulltext_terms=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
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
| `verbose` | bool | True | 运行提示开关。 |

以下参数组至少提供一组：`search_keywords`; `query`。


## 按筛选条件检索 OpenAlex 论文记录 {#dnakit.databases.openalex_search_works-parameters}

`dnakit.databases.openalex_search_works(*, search=None, query=None, filter=None, require_has_fulltext=None, fulltext_terms=None, per_page=None, limit=None, page=None, sort=None, mailto=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
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
| `verbose` | bool | True | 运行提示开关。 |

以下参数组至少提供一组：`search`; `query`。


## 按 PMID 获取 PubMed 文献完整元数据 {#dnakit.databases.pub_med_get_article-parameters}

`dnakit.databases.pub_med_get_article(*, pmid, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pmid` | string &#124; integer | 必填 | PubMed ID (PMID) as string or integer (e.g., &#x27;12345678&#x27; or 12345678, &#x27;33250470&#x27; or 33250470). Find PMIDs using PubMed_search_articles. Can also accept comma-separated PMIDs as string for multiple articles. |
| `verbose` | bool | True | 运行提示开关。 |


## 查询引用指定 PubMed 文献的论文 {#dnakit.databases.pub_med_get_cited_by-parameters}

`dnakit.databases.pub_med_get_cited_by(*, pmid, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pmid` | string | 必填 | PubMed ID (PMID) for which to find citing articles (e.g., &#x27;12345678&#x27;). Find PMIDs using PubMed_search_articles. |
| `limit` | integer | None | Maximum number of citing articles to return (default: 20, max: 100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查找与指定 PubMed 文献相似的论文 {#dnakit.databases.pub_med_get_related-parameters}

`dnakit.databases.pub_med_get_related(*, pmid, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `pmid` | string | 必填 | PubMed ID (PMID) for which to find related articles (e.g., &#x27;20210808&#x27;, &#x27;19879512&#x27;). Find PMIDs using PubMed_search_articles. |
| `limit` | integer | None | Maximum number of related articles to return (default: 20, max: 100). {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按检索式查询 PubMed 生物医学文献 {#dnakit.databases.pub_med_search_articles-parameters}

`dnakit.databases.pub_med_search_articles(*, query, limit=None, mindate=None, maxdate=None, datetype=None, include_abstract=None, sort=None, max_results=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for PubMed articles. Use keywords, author names, journal names, or MeSH terms. Examples: &#x27;cancer immunotherapy&#x27;, &#x27;Smith J&#91;Author&#93;&#x27;, &#x27;Nature&#91;Journal&#93;&#x27;, &#x27;diabetes&#91;MeSH&#93;&#x27;. Use AND/OR/NOT for complex queries. Note: PubMed ANDs every space-separated keyword by default, so a long, highly specific free-text query (5+ concepts) often returns 0 results even when relevant literature exists. Prefer 2-4 key concepts, add explicit OR for synonyms, or try EuropePMC_search_articles for broader relevance-ranked full-text search if this returns 0 results. |
| `limit` | integer | None | Number of articles to return. This sets the maximum number of articles retrieved from PubMed (max: 200). {&quot;default&quot;: 10} |
| `mindate` | string | None | Minimum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set. |
| `maxdate` | string | None | Maximum date for results (format: YYYY/MM/DD or YYYY). Requires datetype to be set. |
| `datetype` | string | None | Type of date to filter on: &#x27;pdat&#x27; (publication date), &#x27;edat&#x27; (Entrez date), or &#x27;mdat&#x27; (modification date). Required when using mindate/maxdate. {&quot;default&quot;: &quot;pdat&quot;} |
| `include_abstract` | boolean | None | If true, best-effort fetches abstracts via efetch (adds abstract/abstract_source fields). {&quot;default&quot;: false} |
| `sort` | string | None | Sort order for results. Valid values: &#x27;pub_date&#x27; (newest first), &#x27;Author&#x27; (alphabetical by first author), &#x27;JournalName&#x27; (alphabetical by journal), &#x27;relevance&#x27; (default PubMed ranking). |
| `max_results` | integer | None | Alias for limit — maximum number of results to return |
| `verbose` | bool | True | 运行提示开关。 |


## 按疾病、药物或关键词查找临床试验 {#dnakit.databases.search_clinical_trials-parameters}

`dnakit.databases.search_clinical_trials(*, condition=None, intervention=None, query_term=None, pageSize=None, pageToken=None, overall_status=None, max_results=None, limit=None, status=None, keyword=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
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
| `verbose` | bool | True | 运行提示开关。 |


## 检索论文、摘要与简短总结 {#dnakit.databases.semantic_scholar_search_papers-parameters}

`dnakit.databases.semantic_scholar_search_papers(*, query, limit=None, year=None, sort=None, include_abstract=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query for Semantic Scholar. Use keywords separated by spaces to refine the search. |
| `limit` | integer | None | Maximum number of papers to return from Semantic Scholar. {&quot;default&quot;: 5} |
| `year` | string | None | Filter results by publication year. Use a single year (e.g., &#x27;2024&#x27;) or a range (e.g., &#x27;2020-2024&#x27;). |
| `sort` | string | None | Sort results. Options: &#x27;citationCount:desc&#x27;, &#x27;citationCount:asc&#x27;, &#x27;publicationDate:desc&#x27;, &#x27;publicationDate:asc&#x27;. |
| `include_abstract` | boolean | None | If true, best-effort fetches missing abstracts via the paper detail endpoint (only when abstract is missing in search results). {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |


## 通过搜索引擎检索网页 {#dnakit.databases.web_search-parameters}

`dnakit.databases.web_search(*, query, max_results=None, search_type=None, backend=None, region=None, safesearch=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Search query string |
| `max_results` | integer | None | Maximum number of results to return {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 50} |
| `search_type` | string | None | Type of search to perform {&quot;enum&quot;: &#91;&quot;general&quot;, &quot;api_documentation&quot;, &quot;python_packages&quot;, &quot;github&quot;&#93;, &quot;default&quot;: &quot;general&quot;} |
| `backend` | string | None | Search engine backend to use. Selecting &#x27;parallel&#x27; sends the query to https://search.parallel.ai/mcp; that backend does not apply ToolUniverse region or safesearch controls. Selecting &#x27;serpbase&#x27; sends the query to https://api.serpbase.dev (requires SERPBASE_API_KEY) for real Google SERP results; that backend applies region but not safesearch. {&quot;enum&quot;: &#91;&quot;auto&quot;, &quot;parallel&quot;, &quot;serpbase&quot;, &quot;google&quot;, &quot;bing&quot;, &quot;brave&quot;, &quot;yahoo&quot;, &quot;duckduckgo&quot;&#93;, &quot;default&quot;: &quot;auto&quot;} |
| `region` | string | None | Search region/locale for DDGS backends, and country/language for the SerpBase backend (split on the hyphen into gl/hl). Parallel does not support this control; only the default value is accepted with backend=&#x27;parallel&#x27;. {&quot;enum&quot;: &#91;&quot;us-en&quot;, &quot;cn-zh&quot;, &quot;uk-en&quot;, &quot;de-de&quot;, &quot;fr-fr&quot;, &quot;ja-jp&quot;&#93;, &quot;default&quot;: &quot;us-en&quot;} |
| `safesearch` | string | None | Safe search level for DDGS backends. Parallel and SerpBase do not support this control; only the default value is accepted with those backends. {&quot;enum&quot;: &#91;&quot;on&quot;, &quot;moderate&quot;, &quot;off&quot;&#93;, &quot;default&quot;: &quot;moderate&quot;} |
| `verbose` | bool | True | 运行提示开关。 |
