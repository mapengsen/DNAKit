# Sequence Editing: parameters and notes

See the [tutorial](../../api/features/04_ops_004_sequence_editing.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Simulate restriction digestion and resulting fragments {#dnakit.databases.dna_virtual_digest-parameters}

`dnakit.databases.dna_virtual_digest(*, operation=None, sequence, enzymes=None, circular=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;virtual_digest&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;virtual_digest&quot;&#93;, &quot;default&quot;: &quot;virtual_digest&quot;} |
| `sequence` | string | Required | DNA sequence (A, T, G, C, N only). Case insensitive. |
| `enzymes` | array &#124; null | None | List of enzyme names to use for digestion (e.g., &#91;&#x27;EcoRI&#x27;, &#x27;BamHI&#x27;&#93;). If null, uses all 24 NEB enzymes. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `circular` | boolean &#124; null | None | Treat sequence as circular DNA (e.g., plasmid). Default: false (linear). {&quot;default&quot;: false} |
| `verbose` | bool | True | Progress flag. |
