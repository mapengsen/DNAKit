# Universal search: parameters and notes

See the [tutorial](../../api/features/09_search.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Search similar protein sequences with remote BLAST {#dnakit.databases.blast_protein_search-parameters}

`dnakit.databases.blast_protein_search(*, sequence, database=None, expect=None, hitlist_size=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `sequence` | string | Required | Protein sequence to search (minimum 10 amino acids) |
| `database` | string | None | Database to search. Options: &#x27;nr&#x27; (comprehensive but slow), &#x27;swissprot&#x27; (curated, faster), &#x27;pdb&#x27; (protein structures only) {&quot;default&quot;: &quot;nr&quot;} |
| `expect` | number | None | E-value threshold for reporting hits {&quot;default&quot;: 10.0} |
| `hitlist_size` | integer | None | Maximum number of hits to return {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |
