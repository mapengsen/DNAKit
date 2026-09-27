# Variant Annotation and Identifier Conversion

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

## 1) VEP annotation from HGVS {#1-vep-annotation-from-hgvs}

- **Purpose:** Obtain transcript consequences for one HGVS expression from Ensembl VEP.
- **Method:** Submit HGVS notation to Ensembl VEP and retain transcript consequence records.
- **API:** `dnakit.annotation.annotate_variant_vep(hgvs_notation[required], species[optional])`
- **Input:** One HGVS expression; species defaults to `human`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import annotate_variant_vep

result = annotate_variant_vep("NM_000546.6:c.215C>G")
print(result.data)
```

- **Example output:** Structured Ensembl VEP consequence records.

## 2) VEP annotation from rsID {#2-vep-annotation-from-rsid}

- **Purpose:** Annotate one dbSNP rsID with Ensembl VEP.
- **Method:** Submit the rsID to Ensembl VEP and return consequence annotations.
- **API:** `dnakit.annotation.annotate_rsid_vep(rsid[required], species[optional])`
- **Input:** An rsID such as `rs28934578`.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import annotate_rsid_vep

result = annotate_rsid_vep("rs28934578")
print(result.data)
```

- **Example output:** VEP consequences associated with the rsID.

## 3) Variant identifier recoding {#3-variant-identifier-recoding}

- **Purpose:** Resolve an rsID or HGVS-like identifier to equivalent Ensembl representations.
- **Method:** Query Ensembl variant recoding for equivalent representations.
- **API:** `dnakit.annotation.recode_variant(variant_id[required], species[optional])`
- **Input:** One variant identifier and an optional species.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit.annotation import recode_variant

result = recode_variant("rs28934578")
print(result.data)
```

- **Example output:** Equivalent HGVS, genomic-location, and database identifiers.

<!-- BEGIN GENERATED API: variant-annotation -->

<span id="variant-annotation"></span>

## 4) Retrieve Ensembl variant consequences and annotations {#dnakit.databases.ensembl_var_get_variant_consequences}

- **Purpose:** Retrieve Ensembl variant consequences and annotations.
- **Method:** Use `EnsemblVariationExt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_var_get_variant_consequences(*, variant_id, species=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `name`, `source`, `most_severe_consequence`, `ancestral_allele` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_var_get_variant_consequences(
    variant_id='rs429358',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "name": "<string>",
  "source": "<string>",
  "most_severe_consequence": "<string>",
  "ancestral_allele": "<string>",
  "minor_allele": "<string>",
  "MAF": "<number>",
  "synonyms": [
    "<string>"
  ],
  "mappings": [
    {
      "seq_region_name": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "allele_string": "<string>",
      "strand": "<integer>",
      "assembly_name": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_var_get_variant_consequences-parameters)

## 5) Annotate HGVS variants with Ensembl VEP {#dnakit.databases.ensembl_vep_annotate_hgvs}

- **Purpose:** Annotate HGVS variants with Ensembl VEP.
- **Method:** Use `EnsemblVEP` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_vep_annotate_hgvs(*, hgvs_notation, species=None, verbose=True)`
- **Input:** Provide `hgvs_notation`.
- **Output:** Returns a mapping with `input`, `assembly_name`, `seq_region_name`, `start` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_vep_annotate_hgvs(
    hgvs_notation='BRAF:p.Val600Glu',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "input": "<string>",
  "assembly_name": "<string>",
  "seq_region_name": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "strand": "<integer>",
  "allele_string": "<string>",
  "most_severe_consequence": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_annotate_hgvs-parameters)

## 6) Annotate rsID variants with Ensembl VEP {#dnakit.databases.ensembl_vep_annotate_rsid}

- **Purpose:** Annotate rsID variants with Ensembl VEP.
- **Method:** Use `EnsemblVEP` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_vep_annotate_rsid(*, variant_id=None, rsid=None, rs_id=None, species=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `input`, `assembly_name`, `seq_region_name`, `start` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_vep_annotate_rsid(
    variant_id='rs7903146',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "input": "<string>",
  "assembly_name": "<string>",
  "seq_region_name": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "allele_string": "<string>",
  "most_severe_consequence": "<string>",
  "transcript_consequences": [
    {
      "gene_symbol": "<string>",
      "variant_allele": "<string>",
      "consequence_terms": "<array>",
      "impact": "<string>",
      "amino_acids": "<string>",
      "sift_prediction": "<string>",
      "polyphen_prediction": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_annotate_rsid-parameters)

## 7) Annotate genomic variants with Ensembl VEP {#dnakit.databases.ensembl_vep_region}

- **Purpose:** Annotate genomic variants with Ensembl VEP.
- **Method:** Use `Ensembl` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_vep_region(*, species, region, allele=None, AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, GeneSplicer=None, GO=None, LoF=None, MaxEntScan=None, Phenotypes=None, SIFT=None, PolyPhen=None, verbose=True)`
- **Input:** Provide `species`, `region`.
- **Output:** Returns a list of records; each record includes `allele_string`, `consequence_terms`, `most_severe_consequence`, `transcript_consequences` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_vep_region(
    species='human',
    region='21:25891796-25891796',
    allele='T',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "allele_string": "<string>",
    "consequence_terms": [
      "<value>"
    ],
    "most_severe_consequence": "<string>",
    "transcript_consequences": [
      "<value>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_region-parameters)

## 8) Convert variant identifiers and HGVS representations {#dnakit.databases.ensembl_vep_variant_recoder}

- **Purpose:** Convert variant identifiers and HGVS representations.
- **Method:** Use `EnsemblVEP` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ensembl_vep_variant_recoder(*, variant_id, species=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a list of records; each record includes `allele`, `input`, `id`, `hgvsg` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ensembl_vep_variant_recoder(
    variant_id='rs429358',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "allele": "<string>",
    "input": "<string>",
    "id": [
      "<string>"
    ],
    "hgvsg": [
      "<string>"
    ],
    "hgvsc": [
      "<string>"
    ],
    "hgvsp": [
      "<string>"
    ],
    "spdi": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_variant_recoder-parameters)

## 9) Annotate GRCh38 variants with FAVOR {#dnakit.databases.favor_annotate_variant}

- **Purpose:** Annotate GRCh38 variants with FAVOR.
- **Method:** Use `FAVORVariantAnnotation` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.favor_annotate_variant(*, variant, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a mapping with `found`, `variant`, `gene_consequence`, `allele_frequency` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.favor_annotate_variant(
    variant='19-44908822-C-T',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "found": "<boolean>",
  "variant": "<object>",
  "gene_consequence": "<object>",
  "allele_frequency": "<object>",
  "deleteriousness": "<object>",
  "conservation": "<object>",
  "clinical": "<object>",
  "regulatory": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.favor_annotate_variant-parameters)

## 10) Retrieve variant annotations by HGVS or rsID {#dnakit.databases.my_variant_get_variant_annotation}

- **Purpose:** Retrieve variant annotations by HGVS or rsID.
- **Method:** Use `MyVariant` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_variant_get_variant_annotation(*, variant_id, fields=None, assembly=None, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `_id`, `dbsnp`, `clinvar`, `cadd` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_variant_get_variant_annotation(
    variant_id='chr1:g.218631822G>A',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "_id": "<string>",
  "dbsnp": {
    "rsid": "<string>",
    "gene": "<object>",
    "alleles": [
      "<value>"
    ]
  },
  "clinvar": "<object>",
  "cadd": {
    "phred": "<number>",
    "consequence": "<string>"
  },
  "gnomad_genome": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.my_variant_get_variant_annotation-parameters)

## 11) Search variants by gene, position or clinical significance {#dnakit.databases.my_variant_query_variants}

- **Purpose:** Search variants by gene, position or clinical significance.
- **Method:** Use `MyVariant` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.my_variant_query_variants(*, query=None, fields=None, size=None, assembly=None, variant_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `took`, `total`, `hits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.my_variant_query_variants(
    query='rs58991260',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "took": "<integer>",
  "total": "<integer>",
  "hits": [
    {
      "_id": "<string>",
      "dbsnp": "<object>",
      "clinvar": "<object>",
      "cadd": "<object>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.my_variant_query_variants-parameters)

## 12) Validate and convert HGVS variant descriptions {#dnakit.databases.variant_validator_validate_variant}

- **Purpose:** Validate and convert HGVS variant descriptions.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.variant_validator_validate_variant(*, genome_build, variant_description, select_transcripts, verbose=True)`
- **Input:** Provide `genome_build`, `variant_description`, `select_transcripts`.
- **Output:** Returns a mapping with `flag`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.variant_validator_validate_variant(
    genome_build='GRCh38',
    variant_description='NM_007294.4:c.5266dup',
    select_transcripts='NM_007294.4',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "flag": "<string>",
  "metadata": {
    "variantvalidator_version": "<string>",
    "variantvalidator_hgvs_version": "<string>",
    "vvdb_version": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/35_variant_annotation.md#dnakit.databases.variant_validator_validate_variant-parameters)

<!-- END GENERATED API: variant-annotation -->
