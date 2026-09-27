# Variant Annotation and Identifier Conversion: parameters and notes

See the [tutorial](../../api/features/35_variant_annotation.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Retrieve Ensembl variant consequences and annotations {#dnakit.databases.ensembl_var_get_variant_consequences-parameters}

`dnakit.databases.ensembl_var_get_variant_consequences(*, variant_id, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | dbSNP rsID of the variant. Examples: &#x27;rs429358&#x27; (APOE), &#x27;rs7903146&#x27; (TCF7L2), &#x27;rs1042779&#x27;. |
| `species` | string | None | Species name. Default &#x27;human&#x27;. Examples: &#x27;human&#x27;, &#x27;homo_sapiens&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Annotate HGVS variants with Ensembl VEP {#dnakit.databases.ensembl_vep_annotate_hgvs-parameters}

`dnakit.databases.ensembl_vep_annotate_hgvs(*, hgvs_notation, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `hgvs_notation` | string | Required | Variant in HGVS notation. Supports protein (p.), coding DNA (c.), and genomic (g.) notation. Examples: &#x27;BRAF:p.Val600Glu&#x27;, &#x27;ENST00000366667:c.803C&gt;T&#x27;, &#x27;NC_000007.14:g.140753336A&gt;T&#x27;. |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. Other options: &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;zebrafish&#x27;, etc. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | Progress flag. |


## Annotate rsID variants with Ensembl VEP {#dnakit.databases.ensembl_vep_annotate_rsid-parameters}

`dnakit.databases.ensembl_vep_annotate_rsid(*, variant_id=None, rsid=None, rs_id=None, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | None | dbSNP rs identifier. Examples: &#x27;rs7903146&#x27; (TCF7L2, type 2 diabetes), &#x27;rs429358&#x27; (APOE, Alzheimer&#x27;s), &#x27;rs1801133&#x27; (MTHFR C677T). Aliases: rsid, rs_id. |
| `rsid` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `rs_id` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | Progress flag. |


## Annotate genomic variants with Ensembl VEP {#dnakit.databases.ensembl_vep_region-parameters}

`dnakit.databases.ensembl_vep_region(*, species, region, allele=None, AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, GeneSplicer=None, GO=None, LoF=None, MaxEntScan=None, Phenotypes=None, SIFT=None, PolyPhen=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `species` | string | Required | Species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `region` | string | Required | Genomic region with variant (e.g., &#x27;21:25891796-25891796&#x27;) |
| `allele` | string | None | Allele argument for a region variant in Ensembl VEP format. |
| `AncestralAllele` | integer | None | Include ancestral allele information {&quot;default&quot;: 0} |
| `Blosum62` | integer | None | Include BLOSUM62 protein conservation score {&quot;default&quot;: 0} |
| `CADD` | integer | None | Include CADD scores {&quot;default&quot;: 0} |
| `Conservation` | integer | None | Include conservation scores {&quot;default&quot;: 0} |
| `DisGeNET` | integer | None | Include DisGeNET disease associations {&quot;default&quot;: 0} |
| `GeneSplicer` | integer | None | Include GeneSplicer predictions {&quot;default&quot;: 0} |
| `GO` | integer | None | Include Gene Ontology terms {&quot;default&quot;: 0} |
| `LoF` | integer | None | Include Loss-of-Function predictions {&quot;default&quot;: 0} |
| `MaxEntScan` | integer | None | Include MaxEntScan splice predictions {&quot;default&quot;: 0} |
| `Phenotypes` | integer | None | Include phenotype annotations {&quot;default&quot;: 0} |
| `SIFT` | string | None | SIFT prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `PolyPhen` | string | None | PolyPhen prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `verbose` | bool | True | Progress flag. |


## Convert variant identifiers and HGVS representations {#dnakit.databases.ensembl_vep_variant_recoder-parameters}

`dnakit.databases.ensembl_vep_variant_recoder(*, variant_id, species=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | Variant identifier to recode. Accepts: rsID (e.g., &#x27;rs429358&#x27;), HGVS notation (e.g., &#x27;AGT:c.803T&gt;C&#x27;), SPDI (e.g., &#x27;NC_000019.10:44908683:T:C&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | Progress flag. |


## Annotate GRCh38 variants with FAVOR {#dnakit.databases.favor_annotate_variant-parameters}

`dnakit.databases.favor_annotate_variant(*, variant, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant` | string | Required | GRCh38 SNV/indel as chr-pos-ref-alt. Accepts &#x27;19-44908822-C-T&#x27;, &#x27;chr19:44908822:C:T&#x27;, etc. |
| `verbose` | bool | True | Progress flag. |


## Retrieve variant annotations by HGVS or rsID {#dnakit.databases.my_variant_get_variant_annotation-parameters}

`dnakit.databases.my_variant_get_variant_annotation(*, variant_id, fields=None, assembly=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `variant_id` | string | Required | HGVS variant ID (e.g., &#x27;chr1:g.2186318G&gt;A&#x27;) or rsID (e.g., &#x27;rs123&#x27;). An HGVS id must be expressed in the assembly given by the &#x27;assembly&#x27; parameter. |
| `fields` | string | None | Fields to return. Useful fields: dbsnp, clinvar, cadd, gnomad_genome, gnomad_exome, dbnsfp (functional predictions). {&quot;default&quot;: &quot;dbsnp,clinvar,cadd,gnomad_genome,dbnsfp&quot;} |
| `assembly` | string | None | Reference genome assembly for both the HGVS variant_id and the returned coordinates: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `verbose` | bool | True | Progress flag. |


## Search variants by gene, position or clinical significance {#dnakit.databases.my_variant_query_variants-parameters}

`dnakit.databases.my_variant_query_variants(*, query=None, fields=None, size=None, assembly=None, variant_id=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Examples: &#x27;rs58991260&#x27; (rsID), &#x27;chr7:g.55249071G&gt;A&#x27; (HGVS), &#x27;clinvar.gene.symbol:BRCA1&#x27; (gene in ClinVar), &#x27;cadd.phred:&gt;30&#x27; (high CADD score). |
| `fields` | string | None | Comma-separated fields to return. Common: dbsnp, clinvar, cadd, gnomad_genome, dbnsfp. {&quot;default&quot;: &quot;dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,gnomad_genome.af.af&quot;} |
| `size` | integer | None | Maximum number of results (1-100). {&quot;default&quot;: 10} |
| `assembly` | string | None | Reference genome assembly for the coordinates in the query and in the results: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `variant_id` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | Progress flag. |


## Validate and convert HGVS variant descriptions {#dnakit.databases.variant_validator_validate_variant-parameters}

`dnakit.databases.variant_validator_validate_variant(*, genome_build, variant_description, select_transcripts, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `genome_build` | string | Required | Reference genome assembly: &#x27;GRCh37&#x27; (hg19) or &#x27;GRCh38&#x27; (hg38) |
| `variant_description` | string | Required | HGVS variant description (e.g., &#x27;NM_007294.4:c.5266dup&#x27; for BRCA1 c.5266dupC, &#x27;NM_000179.3:c.943C&gt;T&#x27;, &#x27;chr17:g.43092919del&#x27;) |
| `select_transcripts` | string | Required | Transcript to validate against (e.g., &#x27;NM_007294.4&#x27;). Use &#x27;all&#x27; to get all transcripts for a genomic variant. |
| `verbose` | bool | True | Progress flag. |
