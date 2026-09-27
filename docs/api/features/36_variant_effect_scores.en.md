# Variant Effect Scores and Experimental Data

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: variant-scores -->

<span id="variant-scores"></span>

## 1) Retrieve AlphaMissense classifications across a protein {#dnakit.databases.alpha_missense_get_protein_scores}

- **Purpose:** Aggregate substitution classes, residue mean scores and thresholds across a UniProt protein; max_residues limits the queried positions.
- **Method:** Use `AlphaMissense` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alpha_missense_get_protein_scores(*, uniprot_id, max_residues=None, verbose=True)`
- **Input:** Provide `uniprot_id`.
- **Output:** Returns a mapping with `uniprot_id`, `protein_length`, `scores`, `n_positions_returned` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_protein_scores(
    uniprot_id='P01116',
    max_residues=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "uniprot_id": "<string>",
  "protein_length": "<integer>",
  "scores": [
    "<object>"
  ],
  "n_positions_returned": "<integer>",
  "n_positions_attempted": "<integer>",
  "max_residues_cap": "<integer>",
  "thresholds": "<object>",
  "pdb_download": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_protein_scores-parameters)

## 2) Retrieve substitution classes and mean scores at a residue {#dnakit.databases.alpha_missense_get_residue_scores}

- **Purpose:** Retrieve benign, ambiguous and pathogenic substitution classes plus residue mean scores; these are not numeric scores for each individual substitution.
- **Method:** Use `AlphaMissense` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alpha_missense_get_residue_scores(*, uniprot_id, position, verbose=True)`
- **Input:** Provide `uniprot_id`, `position`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_residue_scores(
    uniprot_id='P15056',
    position=600,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "uniprot_id": "<string>",
    "position": "<integer>",
    "scores": "<object>",
    "thresholds": "<object>"
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_residue_scores-parameters)

## 3) Classify a missense variant with AlphaMissense {#dnakit.databases.alpha_missense_get_variant_score}

- **Purpose:** Retrieve an AlphaMissense classification for one missense variant. The upstream operation supplies no numeric substitution score: pathogenicity_score is null and score_available is false.
- **Method:** Use `AlphaMissense` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alpha_missense_get_variant_score(*, uniprot_id, variant, verbose=True)`
- **Input:** Provide `uniprot_id`, `variant`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_variant_score(
    uniprot_id='P15056',
    variant='V600E',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "uniprot_id": "<string>",
    "variant": "<string>",
    "position": "<integer>",
    "reference_aa": "<string>",
    "variant_aa": "<string>",
    "pathogenicity_score": "<number>",
    "score_available": "<boolean>",
    "classification": "<string>"
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_variant_score-parameters)

## 4) Retrieve CADD scores for substitutions at a position {#dnakit.databases.cadd_get_position_scores}

- **Purpose:** Retrieve CADD scores for substitutions at a position.
- **Method:** Use `CADD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cadd_get_position_scores(*, chrom, pos, version=None, verbose=True)`
- **Input:** Provide `chrom`, `pos`.
- **Output:** Returns a mapping with `chrom`, `pos`, `variants`, `version` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cadd_get_position_scores(
    chrom='7',
    pos=140753336,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "chrom": "<string>",
  "pos": "<integer>",
  "variants": [
    "<value>"
  ],
  "version": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_position_scores-parameters)

## 5) Retrieve CADD scores across a genomic interval {#dnakit.databases.cadd_get_range_scores}

- **Purpose:** Retrieve CADD scores across a genomic interval.
- **Method:** Use `CADD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cadd_get_range_scores(*, chrom, start, end, version=None, verbose=True)`
- **Input:** Provide `chrom`, `start`, `end`.
- **Output:** Returns a mapping with `chrom`, `start`, `end`, `variants` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cadd_get_range_scores(
    chrom='7',
    start=140753330,
    end=140753340,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "chrom": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "variants": [
    "<value>"
  ],
  "count": "<integer>",
  "version": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_range_scores-parameters)

## 6) Retrieve indexed CADD annotations through MyVariant {#dnakit.databases.cadd_get_scores}

- **Purpose:** Retrieve indexed CADD annotations through MyVariant.
- **Method:** Use `API` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cadd_get_scores(*, variant, assembly='hg38', verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns MyVariant query results; mapping results include dnakit_score_source identifying indexed CADD annotations. Missing scores remain missing; no new model prediction is performed.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cadd_get_scores(
    variant='chr7:g.140753336A>T',
    assembly='hg38',
    verbose=False,
)
print(result)
```

- **Example output:** Returns MyVariant query results; mapping results include dnakit_score_source identifying indexed CADD annotations. Missing scores remain missing; no new model prediction is performed.

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_scores-parameters)

## 7) Retrieve a variant's CADD deleteriousness score {#dnakit.databases.cadd_get_variant_score}

- **Purpose:** Retrieve a variant's CADD deleteriousness score.
- **Method:** Use `CADD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cadd_get_variant_score(*, chrom, pos, ref, alt, version=None, include_annotations=None, verbose=True)`
- **Input:** Provide `chrom`, `pos`, `ref`, `alt`.
- **Output:** Returns a mapping with `chrom`, `pos`, `ref`, `alt` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cadd_get_variant_score(
    chrom='7',
    pos=140753336,
    ref='A',
    alt='T',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "chrom": "<string>",
  "pos": "<integer>",
  "ref": "<string>",
  "alt": "<string>",
  "phred_score": "<number>",
  "raw_score": "<number>",
  "interpretation": "<string>",
  "version": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_variant_score-parameters)

## 8) Check EVE score availability for a gene {#dnakit.databases.eve_get_gene_info}

- **Purpose:** Check EVE score availability for a gene.
- **Method:** Use `EVE` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.eve_get_gene_info(*, gene_symbol, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.eve_get_gene_info(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "gene_symbol": "<string>",
    "ensembl_id": "<string>",
    "eve_website": "<string>"
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.eve_get_gene_info-parameters)

## 9) Retrieve an EVE variant pathogenicity score {#dnakit.databases.eve_get_variant_score}

- **Purpose:** Retrieve an EVE variant pathogenicity score.
- **Method:** Use `EVE` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.eve_get_variant_score(*, variant=None, chrom=None, pos=None, ref=None, alt=None, species=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.eve_get_variant_score(
    chrom='17',
    pos=7674220,
    ref='G',
    alt='A',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "variant": "<string>",
    "eve_scores": [
      "<value>"
    ],
    "threshold": "<object>"
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.eve_get_variant_score-parameters)

## 10) Build an amino acid effect matrix from MaveDB scores {#dnakit.databases.mave_db_get_effect_matrix}

- **Purpose:** Build an amino acid effect matrix from MaveDB scores.
- **Method:** Use `MaveDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.mave_db_get_effect_matrix(*, urn, uniprot_id=None, score_field=None, verbose=True)`
- **Input:** Provide `urn`.
- **Output:** Returns a mapping with `urn`, `matrix`, `positions`, `amino_acid_order` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.mave_db_get_effect_matrix(
    urn='urn:mavedb:00000115-a-7',
    uniprot_id='P01116',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "urn": "<string>",
  "matrix": [
    [
      "<number>"
    ]
  ],
  "positions": [
    "<integer>"
  ],
  "amino_acid_order": "<string>",
  "shape": [
    "<integer>"
  ],
  "n_parsed_single_missense": "<integer>",
  "n_dropped": "<integer>",
  "n_total_rows": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.mave_db_get_effect_matrix-parameters)

## 11) Retrieve multiple variant pathogenicity prediction scores {#dnakit.databases.my_variant_get_pathogenicity_scores}

- **Purpose:** Retrieve multiple variant pathogenicity prediction scores.
- **Method:** Use `MyVariant` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_variant_get_pathogenicity_scores(*, variant_id, fields=None, assembly=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns existing dbNSFP variant scores such as REVEL, CADD and SIFT where available; missing scores are not calculated.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_variant_get_pathogenicity_scores(
    variant_id='rs45478192',
    verbose=False,
)
print(result)
```

- **Example output:** Returns existing dbNSFP variant scores such as REVEL, CADD and SIFT where available; missing scores are not calculated.

[Full parameters and usage notes](../../faq/api/36_variant_effect_scores.md#dnakit.databases.my_variant_get_pathogenicity_scores-parameters)

<!-- END GENERATED API: variant-scores -->
