# DNA descriptor: parameters and notes

See the [tutorial](../../api/features/05_all_descriptors.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Count a nucleotide or amino acid residue {#dnakit.databases.sequence_count_residues-parameters}

`dnakit.databases.sequence_count_residues(*, operation=None, sequence=None, uniprot_id=None, residue, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;count_residues&quot;&#93;, &quot;default&quot;: &quot;count_residues&quot;} |
| `sequence` | string | None | Sequence string (DNA/RNA/protein). Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch sequence from (e.g., &#x27;P00533&#x27; for EGFR). Optional if sequence is provided. |
| `residue` | string | Required | Single character residue to count (e.g., &#x27;C&#x27; for cysteine) |
| `verbose` | bool | True | Progress flag. |


## Calculate GC content of a nucleotide sequence {#dnakit.databases.sequence_gc_content-parameters}

`dnakit.databases.sequence_gc_content(*, operation=None, sequence, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;gc_content&quot;&#93;, &quot;default&quot;: &quot;gc_content&quot;} |
| `sequence` | string | Required | DNA or RNA sequence string |
| `verbose` | bool | True | Progress flag. |


## Summarize sequence length, composition and type-specific properties {#dnakit.databases.sequence_stats-parameters}

`dnakit.databases.sequence_stats(*, operation=None, sequence=None, uniprot_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;stats&quot;&#93;, &quot;default&quot;: &quot;stats&quot;} |
| `sequence` | string | None | Sequence string. Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch protein sequence from. Optional if sequence is provided. |
| `verbose` | bool | True | Progress flag. |
