# sequence completion: parameters and notes

See the [tutorial](../../api/features/04_ops_001_sequence_direction.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Compute a DNA reverse complement {#dnakit.databases.sequence_reverse_complement-parameters}

`dnakit.databases.sequence_reverse_complement(*, operation=None, sequence, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;reverse_complement&quot;&#93;, &quot;default&quot;: &quot;reverse_complement&quot;} |
| `sequence` | string | Required | DNA sequence (A/T/C/G/N only) |
| `verbose` | bool | True | Progress flag. |
