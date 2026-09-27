# Expression, Regulatory Tracks and Contact Maps: parameters and notes

See the [tutorial](../../api/features/23_prediction_regulatory_tracks.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

These functions call existing models or prediction services. Weights, SDKs, credentials and service terms remain provider requirements; DNAKit did not retrain these models.

## Predict multimodal genomic tracks with AlphaGenome {#dnakit.databases.alpha_genome_predict_interval-parameters}

`dnakit.databases.alpha_genome_predict_interval(*, chromosome, start, end, output_types=None, ontology_terms=None, organism=None, sequence_length=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `chromosome` | string | Required | Chromosome, e.g. &#x27;chr19&#x27;. |
| `start` | integer | Required | Interval start (0-based). |
| `end` | integer | Required | Interval end. |
| `output_types` | array &#124; null | None | Modalities to predict, e.g. &#91;&#x27;RNA_SEQ&#x27;,&#x27;ATAC&#x27;&#93; (default &#91;&#x27;RNA_SEQ&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `ontology_terms` | array &#124; null | None | Optional tissue/cell ontology terms (e.g. &#91;&#x27;UBERON:0001114&#x27;&#93; = liver) to restrict tracks. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `organism` | string &#124; null | None | &#x27;human&#x27; (default) or &#x27;mouse&#x27;. |
| `sequence_length` | string &#124; null | None | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `verbose` | bool | True | Progress flag. |

Configure before running: `ALPHA_GENOME_API_KEY`.
