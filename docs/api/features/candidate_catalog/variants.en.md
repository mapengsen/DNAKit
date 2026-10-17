# Variant analysis

This category contains **63** documented functions. 51 of them compute their result from the data passed in (VCF text, VCF paths or variant records, genotype counts, read counts, DNA sequences) with the methods described in each section; 12 always return `status: "not_implemented"` because they need annotation databases, trained predictors, external callers or aligned reads that DNAKit does not bundle (`annotate`, `annotate_conservation_score`, `annotate_mt`, `annotate_pathogenicity_score`, `annotate_population_frequency`, `annotate_regulatory_effect`, `annotate_splice_effect`, `annotate_variant_consequence`, `call`, `call_mt`, `detect_breakpoints`, `genotype_structural_variants`). Five of the computing functions need an option and raise `ConfigurationError` without it: `left_align`, `left_align_mt` and `left_align_variants` need `options.reference`, `detect_repeat_expansions` needs `options.motif` and `phase_genotypes` needs `options.trio`. Several functions are deliberate approximations and say so in their Method section: `fill_missing_genotypes_with_mode` is naive per-site modal-genotype filling, not Beagle-style imputation; `calculate_kinship_king` is the KING-robust kinship coefficient, not IBD-segment detection; the heuristics `detect_low_frequency_variants`, `detect_mosaic_variants`, `calculate_runs_of_homozygosity` and `prioritize_candidate_variants` are not validated callers or classifiers. The `_mt` functions take VCF text, paths or records like their plain counterparts, except `count_mt`, which only counts non-blank lines or list items; no Hail MatrixTable object is read. Malformed input or a missing required option raises `InputFormatError` or `ConfigurationError` (with a `code` and a `hint`) instead of returning an envelope. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Annotate variants with functional, population and clinical information. Not implemented: it only reports what is missing (transcript models, population frequencies and clinical databases).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Use `annotate_coding_effect` (works on a supplied CDS) or the external-program wrappers `annotate_vcf_ensembl_vep` and `annotate_vcf_vcfanno` instead.

- **API:** `dnakit.variants.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("annotation databases (Ensembl/VEP cache, gnomAD, ClinVar) or the vep / vcfanno programs") and `message` ("annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "annotation databases (Ensembl/VEP cache, gnomAD, ClinVar) or the vep / vcfanno programs",
  "message": "annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle."
}
```

## 2) annotate_coding_effect

- **Purpose:** Predict the protein-level effect (synonymous, missense, nonsense, start/stop loss, frameshift, in-frame indel) of variants on a user-supplied coding sequence.

- **Method:** The reference and the mutated CDS are translated with the standard genetic code (no transcript model and no genome are used) and compared. Positions are 1-based CDS coordinates and REF/ALT are on the coding strand; shared flanking bases of REF/ALT are trimmed first and REF must match the CDS. Same-length changes are classified codon by codon as `start_lost`, `stop_gained`, `stop_lost`, `missense_variant` or `synonymous_variant` (the most severe one wins for multi-codon changes); length changes that are not a multiple of 3 are `frameshift_variant` (or `stop_gained`/`start_lost` when applicable); other length changes are `inframe_insertion`, `inframe_deletion` or `protein_altering_variant`. HGVS protein strings (`p.Ala2Val`, `p.Lys3Ter`, `p.Ala2=`, `p.Met1?`, `fs`, `ins`, `del`, `delins`) are built from the first changed residue and are not 3' normalised. Not a substitute for VEP or a transcript-aware annotator: genomic variants, UTRs, introns and splicing are out of scope.

- **API:** `dnakit.variants.annotate_coding_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'cds': 'ATG...', 'variants': [...]}` where each variant is `{'cds_position': 5, 'ref': 'C', 'alt': 'T'}` or an `'c.5C>T'` string (a single mapping with `cds_position` is also accepted; `'-'` or `''` is an empty allele). `options.cds`: coding sequence (ACGT only, at most 1,000,000 bases) when it is not part of `input`. At most 10,000 variants.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `cds_length`, `protein_length` and `effects`, one entry per variant with `cds_position`, `ref`, `alt` (after trimming), `consequence`, `protein_position` and `hgvs_p`; substitutions also have `codons` (`codon_index`, `ref_codon`, `alt_codon`, `ref_aa`, `alt_aa`). `warnings` notes a CDS whose length is not a multiple of 3 or that does not start with ATG. A REF that does not match the CDS or a position outside the CDS raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_coding_effect

cds = "ATGGCCAAATGGTAG"  # M A K W stop

result = annotate_coding_effect({
    "cds": cds,
    "variants": [
        "c.5C>T",                                   # GCC -> GTC  Ala2Val
        {"cds_position": 7, "ref": "A", "alt": "T"},  # AAA -> TAA  Lys3Ter
        {"cds_position": 6, "ref": "C", "alt": "A"},  # GCC -> GCA  Ala2Ala
    ],
})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "annotate_coding_effect",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 140,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8a61ab07dc9daab0ac0e6b74d5f6a734627313a7cf8477a2dc29a261291ad048"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "translation of the reference and mutated CDS with the standard genetic code; REF/ALT are o...",
    "cds_length": 15,
    "protein_length": 4,
    "effects": [
      {
        "cds_position": 5,
        "ref": "C",
        "alt": "T",
        "consequence": "missense_variant",
        "codons": [{"codon_index": 2, "ref_codon": "GCC", "alt_codon": "GTC", "ref_aa": "A", "alt_aa": "V"}],
        "protein_position": 2,
        "hgvs_p": "p.Ala2Val"
      },
      {
        "cds_position": 7,
        "ref": "A",
        "alt": "T",
        "consequence": "stop_gained",
        "codons": [{"codon_index": 3, "ref_codon": "AAA", "alt_codon": "TAA", "ref_aa": "K", "alt_aa": "*"}],
        "protein_position": 3,
        "hgvs_p": "p.Lys3Ter"
      },
      {
        "cds_position": 6,
        "ref": "C",
        "alt": "A",
        "consequence": "synonymous_variant",
        "codons": [{"codon_index": 2, "ref_codon": "GCC", "alt_codon": "GCA", "ref_aa": "A", "alt_aa": "A"}],
        "protein_position": 2,
        "hgvs_p": "p.Ala2="
      }
    ]
  }
}
```

## 3) annotate_conservation_score

- **Purpose:** Annotate variants with phyloP, phastCons or GERP conservation scores. Not implemented: it only reports what is missing (a genome-wide per-base conservation track (bigWig/bedGraph)).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Nothing is read from a conservation track.

- **API:** `dnakit.variants.annotate_conservation_score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("per-base conservation track (phyloP, phastCons or GERP)") and `message` ("conservation scoring needs a genome-wide conservation track that is not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_conservation_score

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_conservation_score(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_conservation_score",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "per-base conservation track (phyloP, phastCons or GERP)",
  "message": "conservation scoring needs a genome-wide conservation track that is not available."
}
```

## 4) annotate_mt

- **Purpose:** Annotate mitochondrial variants. Not implemented: it only reports what is missing (mitochondrial annotation resources such as MITOMAP, gnomAD-mt or HmtVar).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Nothing is looked up in a mitochondrial database.

- **API:** `dnakit.variants.annotate_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("mitochondrial annotation resources (MITOMAP, gnomAD-mt, HmtVar)") and `message` ("mitochondrial annotation needs reference annotation databases that are not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_mt(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_mt",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "mitochondrial annotation resources (MITOMAP, gnomAD-mt, HmtVar)",
  "message": "mitochondrial annotation needs reference annotation databases that are not available."
}
```

## 5) annotate_pathogenicity_score

- **Purpose:** Annotate variants with a pathogenicity prediction score (CADD, REVEL, AlphaMissense, ...). Not implemented: it only reports what is missing (a trained predictor or precomputed score tables).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No score is estimated; `prioritize_candidate_variants` can combine scores that are already on the records.

- **API:** `dnakit.variants.annotate_pathogenicity_score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("trained pathogenicity predictor or precomputed score tables (CADD/REVEL/AlphaMissense)") and `message` ("pathogenicity scores come from trained models or precomputed score files that are not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_pathogenicity_score

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_pathogenicity_score(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_pathogenicity_score",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained pathogenicity predictor or precomputed score tables (CADD/REVEL/AlphaMissense)",
  "message": "pathogenicity scores come from trained models or precomputed score files that are not available."
}
```

## 6) annotate_population_frequency

- **Purpose:** Annotate variants with population allele frequencies (gnomAD, 1000 Genomes, dbSNP). Not implemented: it only reports what is missing (a frequency database or network access, neither of which is used).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No query is sent and no frequency file is read.

- **API:** `dnakit.variants.annotate_population_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("population frequency database (gnomAD, 1000 Genomes, dbSNP)") and `message` ("population frequencies need a frequency database or network access, neither of which is used here."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_population_frequency

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_population_frequency(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_population_frequency",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "population frequency database (gnomAD, 1000 Genomes, dbSNP)",
  "message": "population frequencies need a frequency database or network access, neither of which is used here."
}
```

## 7) annotate_regulatory_effect

- **Purpose:** Annotate variants with regulatory-element effects (Ensembl Regulatory Build, ENCODE cCREs, sequence-to-activity models). Not implemented: it only reports what is missing (a regulatory-region track or a trained model).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No regulatory annotation is attempted.

- **API:** `dnakit.variants.annotate_regulatory_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("regulatory annotation track or trained regulatory model") and `message` ("regulatory annotation needs regulatory-region tracks or trained models that are not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_regulatory_effect

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_regulatory_effect(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_regulatory_effect",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "regulatory annotation track or trained regulatory model",
  "message": "regulatory annotation needs regulatory-region tracks or trained models that are not available."
}
```

## 8) annotate_splice_effect

- **Purpose:** Annotate variants with predicted splicing effects (SpliceAI, MaxEntScan). Not implemented: it only reports what is missing (transcript models and a trained splice predictor).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No splice prediction is attempted.

- **API:** `dnakit.variants.annotate_splice_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("transcript annotation and a splice-effect predictor (SpliceAI)") and `message` ("splice-effect prediction needs transcript models and a trained predictor that are not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_splice_effect

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_splice_effect(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "annotate_splice_effect",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "transcript annotation and a splice-effect predictor (SpliceAI)",
  "message": "splice-effect prediction needs transcript models and a trained predictor that are not available."
}
```

## 9) annotate_variant_consequence

- **Purpose:** Annotate Sequence Ontology consequences of variants on transcripts. Not implemented: it only reports what is missing (transcript models (GTF/GFF) and a reference genome).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Only `annotate_coding_effect` works, and only from a user-supplied CDS; for genomic variants use `annotate_vcf_ensembl_vep`.

- **API:** `dnakit.variants.annotate_variant_consequence(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("transcript annotation (GTF/GFF) and genome FASTA, e.g. via Ensembl VEP") and `message` ("consequence calling needs transcript models and a reference genome; only annotate_coding_effect works from a user-supplied CDS."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import annotate_variant_consequence

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_variant_consequence(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "annotate_variant_consequence",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "transcript annotation (GTF/GFF) and genome FASTA, e.g. via Ensembl VEP",
  "message": "consequence calling needs transcript models and a reference genome; only annotate_coding_e..."
}
```

## 10) calculate_allele_frequency

- **Purpose:** Document a reusable calculate allele frequency capability for Variant analysis workflows. Candidate CAP-0768 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0768",
  "operation": "calculate_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'alt_count': 3, 'called_alleles': 6, 'frequency': 0.5}], 'total_alt_alleles': 3, 'total_called_alleles': 6, 'aggregate_frequency': 0.5}, 'provenance': ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 11) calculate_copy_number_profile

- **Purpose:** Convert per-bin sequencing depth into a copy-number profile: depth ratio to the median, copy number and merged segments.

- **Method:** ratio = depth / median depth (divided by the same ratio of `control_depth` when a matched control is given); copy number = (ratio * (purity * ploidy + 2 * (1 - purity)) - 2 * (1 - purity)) / purity, floored at 0. The median is taken over bins with positive depth outside `median_exclude_chroms`. Adjacent bins on the same chromosome with the same rounded copy number are merged into one segment and labelled `homozygous_deletion`, `loss`, `neutral`, `gain` (up to ploidy + 2) or `amplification`. This is a depth-ratio conversion only: no statistical segmentation (CBS, HMM), no GC or mappability correction and no B-allele information are used.

- **API:** `dnakit.variants.calculate_copy_number_profile(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of bins `{'chrom', 'start', 'end', 'depth'[, 'control_depth']}` (or `[chrom, start, end, depth]`), or `{'bins': [...]}`; `control_depth` must be given for every bin or for none. `options`: `purity` (default 1.0, 0.01-1), `ploidy` (default 2.0), `median_exclude_chroms` (default X, Y, M, MT, with or without `chr`).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `purity`, `ploidy`, `median_depth`, `control_median_depth`, `bins` (per bin: `chrom`, `start`, `end`, `depth`, `ratio`, `log2_ratio`, `copy_number`, `rounded_copy_number`) and `segments` (`chrom`, `start`, `end`, `copy_number`, `n_bins`, `mean_ratio`, `mean_log2_ratio`, `state`).

- **Example code:**

```python
import json

from dnakit.variants import calculate_copy_number_profile

bins = [
    ("chr1", 1, 1000, 30), ("chr1", 1001, 2000, 31), ("chr1", 2001, 3000, 60),
    ("chr1", 3001, 4000, 62), ("chr2", 1, 1000, 29), ("chr2", 1001, 2000, 15),
]

result = calculate_copy_number_profile(bins)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_copy_number_profile",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 150,
    "sequence_bases": 18,
    "lengths": [3, 3, 3, 3, 3, 3],
    "sha256": "711b9ccb7ed0139eca19642183f362a0f51ae054fb02c72c67491cc8e2648b6f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-bin depth ratio to the median (or to the control depth ratio), copy number = (ratio*(p...",
    "purity": 1.0,
    "ploidy": 2.0,
    "median_depth": 30.5,
    "control_median_depth": null,
    "bins": [
      {
        "chrom": "chr1",
        "start": 1,
        "end": 1000,
        "depth": 30.0,
        "ratio": 0.983607,
        "log2_ratio": -0.0238467,
        "copy_number": 1.96721,
        "rounded_copy_number": 2
      },
      {
        "chrom": "chr1",
        "start": 1001,
        "end": 2000,
        "depth": 31.0,
        "ratio": 1.01639,
        "log2_ratio": 0.023459,
        "copy_number": 2.03279,
        "rounded_copy_number": 2
      },
      "..."
    ],
    "segments": [
      {
        "chrom": "chr1",
        "start": 1,
        "end": 2000,
        "copy_number": 2,
        "n_bins": 2,
        "mean_ratio": 1.0,
        "mean_log2_ratio": 0.0,
        "state": "neutral"
      },
      {
        "chrom": "chr1",
        "start": 2001,
        "end": 4000,
        "copy_number": 4,
        "n_bins": 2,
        "mean_ratio": 2.0,
        "mean_log2_ratio": 1.0,
        "state": "gain"
      },
      "..."
    ]
  }
}
```

## 12) calculate_genotype_likelihood

- **Purpose:** Compute genotype likelihoods (log10), Phred-scaled PL values and GQ for a biallelic site from read evidence under a binomial error model.

- **Method:** For a genotype with k ALT copies out of m (ploidy), each read independently supports ALT with probability (k/m)(1-e) + (1-k/m)e, where e is the per-read error probability (from the base quality as 10^(-q/10), capped at 0.5, or `error_rate` when no quality is given). log10 L is the sum over reads; PL = -10 * (log10 L - max log10 L), rounded; GQ is the second-smallest PL capped at 99. Reads are treated as independent and there is no prior, mapping-quality or strand model, so this is the textbook model and not the GATK or bcftools implementation.

- **API:** `dnakit.variants.calculate_genotype_likelihood(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'ref_count': int, 'alt_count': int}`, `{'reads': [{'allele': 'ref'|'alt', 'base_quality': q}, ...]}`, a list of such mappings, or a VCF/records with `FORMAT/AD` per sample (multi-allelic sites are split first). `options`: `error_rate` (default 0.01, used when no base quality is given), `ploidy` (default 2, 1-6).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `model`, `error_rate`, `ploidy` and either `result` (single count or read input) or `results` (list input or VCF input, where each entry also has `site` and `sample`). Each likelihood entry has `genotypes`, `log10_likelihoods`, `phred_scaled`, `best_genotype`, `genotype_quality`, `ref_reads` and `alt_reads`. Calls without a usable AD are skipped with a warning.

- **Example code:**

```python
import json

from dnakit.variants import calculate_genotype_likelihood

result = calculate_genotype_likelihood({"ref_count": 10, "alt_count": 0})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_genotype_likelihood",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b06a2868736086281799ef5b0bab69005c6fccac827322fd56bb70f49b9e82b6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "each read independently supports ALT with probability (k/m)(1-e)+(1-k/m)e for a genotype w...",
    "error_rate": 0.01,
    "ploidy": 2,
    "result": {
      "genotypes": ["0/0", "0/1", "1/1"],
      "log10_likelihoods": [-0.0436481, -3.0103, -20.0],
      "phred_scaled": [0, 30, 200],
      "best_genotype": "0/0",
      "genotype_quality": 30,
      "ref_reads": 10,
      "alt_reads": 0
    }
  }
}
```

## 13) calculate_haplotype_frequency

- **Purpose:** Estimate haplotype frequencies at two or more loci: direct counting of phased haplotypes, or EM estimation from unphased genotypes.

- **Method:** Phased haplotypes are counted directly. Unphased diploid genotypes are resolved with the EM algorithm of Excoffier & Slatkin (1995): all haplotype pairs compatible with each individual's genotypes are enumerated, frequencies start uniform over the candidate haplotypes, and the E/M steps assume Hardy-Weinberg proportions until the largest change is below `tolerance`. Individuals with a missing genotype are excluded; EM supports 2-12 loci and at most 10 heterozygous loci per individual. Allele codes are integers (0 = REF, 1 = first ALT).

- **API:** `dnakit.variants.calculate_haplotype_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'haplotypes': ['010', '011', ...]}` (phased), `{'genotypes': [['0/1', '1/1'], ...]}` or a list of per-individual GT lists (one GT per locus), or a VCF/records with sample genotypes (loci are the biallelic records). `options`: `use_phase` (default True; for VCF input haplotypes are counted directly when every genotype is phased), `max_iterations` (default 1000), `tolerance` (default 1e-10).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `n_chromosomes`, `loci` and `haplotypes` (per haplotype: `haplotype` such as `0-1`, `alleles`, `count`, `frequency`, sorted by decreasing frequency); the EM branch also has `n_individuals` and `iterations`. `warnings` reports excluded individuals or skipped multi-allelic sites. EM counts and frequencies are floating-point estimates.

- **Example code:**

```python
import json

from dnakit.variants import calculate_haplotype_frequency

# unphased genotypes of 4 individuals at 2 loci (EM estimation)
result = calculate_haplotype_frequency({"genotypes": [["0/0", "0/0"], ["0/1", "0/1"], ["1/1", "1/1"], ["0/0", "0/1"]]})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_haplotype_frequency",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 79,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "17383568bffa3f8dde9185cc6b71f19f7de028428802b7bc3f05c24d779147b1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "EM algorithm for haplotype frequencies from unphased genotypes (Excoffier & Slatkin 1995;...",
    "n_individuals": 4,
    "n_chromosomes": 8,
    "iterations": 11,
    "loci": ["locus_1", "locus_2"],
    "haplotypes": [
      {"haplotype": "0-0", "alleles": [0, 0], "count": 4.0, "frequency": 0.5},
      {"haplotype": "1-1", "alleles": [1, 1], "count": 3.0, "frequency": 0.375},
      "..."
    ]
  }
}
```

## 14) calculate_heterozygosity

- **Purpose:** Document a reusable calculate heterozygosity capability for Variant analysis workflows. Candidate CAP-0772 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate heterozygosity operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_heterozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_heterozygosity` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_heterozygosity

result = calculate_heterozygosity(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0772",
  "operation": "calculate_heterozygosity",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_heterozygosity', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'heterozygous_samples': 1, 'called_samples': 3, 'observed_heterozygosity': 0.3333333333333333}], 'heterozygous_samples': 1, 'called_samples': 3, 'observe...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 15) calculate_hwe_pvalue

- **Purpose:** Test Hardy-Weinberg equilibrium with the exact test of Wigginton, Cutler & Abecasis (2005).

- **Method:** Given the genotype counts, the probabilities of all heterozygote counts compatible with the observed allele counts are computed with the exact recursion; the two-sided p-value is the sum of the probabilities of all tables that are not more likely than the observed one. With `mid_p` half the probability of the observed table is subtracted. For VCF input one p-value is computed per biallelic diploid site (multi-allelic and non-diploid calls are skipped with a warning).

- **API:** `dnakit.variants.calculate_hwe_pvalue(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: genotype counts `[hom_ref, het, hom_alt]` (or `{'counts': [...]}` or `{'hom_ref', 'het', 'hom_alt'}`), or a VCF/records with sample genotypes. `options`: `mid_p` (default False).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. For counts, `result` has `method`, `mid_p`, `n_hom_ref`, `n_het`, `n_hom_alt`, `n_samples`, `observed_het_fraction`, `expected_het_fraction`, `ref_allele_frequency` and `p_value`. For VCF input, `result.sites` holds the same fields plus `site` for every site.

- **Example code:**

```python
import json

from dnakit.variants import calculate_hwe_pvalue

result = calculate_hwe_pvalue([10, 20, 30])  # hom_ref, het, hom_alt counts

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_hwe_pvalue",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 12,
    "sequence_bases": 60,
    "lengths": [10, 20, 30],
    "sha256": "15d98dad2ba94449f66d415d63b1e633ed896a72c15c6b91407410b8e06f5e20"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "exact test of Hardy-Weinberg equilibrium (Wigginton, Cutler & Abecasis 2005), two-sided",
    "mid_p": false,
    "n_hom_ref": 10,
    "n_het": 20,
    "n_hom_alt": 30,
    "n_samples": 60,
    "observed_het_fraction": 0.333333,
    "expected_het_fraction": 0.444444,
    "ref_allele_frequency": 0.333333,
    "p_value": 0.0776747
  }
}
```

## 16) calculate_kinship_king

- **Purpose:** Estimate pairwise kinship between samples with the KING-robust estimator. This is a genotype-based kinship coefficient, not an IBD-segment detector.

- **Method:** For every sample pair, over biallelic sites called in both samples: kinship = (N_AaAa - 2 N_AAaa) / (N_Aa(i) + N_Aa(j)) (Manichaikul et al. 2010), where N_AaAa is the number of sites heterozygous in both, N_AAaa the number of opposite homozygotes and N_Aa(i) the heterozygote count of sample i. `pi_hat_approx` is simply 2 * kinship, `ibs0_fraction` is the opposite-homozygote fraction, and `relationship` is a label from the kinship thresholds of 0.3535 (duplicate or monozygotic twin), 0.1768 (first degree), 0.0884 (second degree) and 0.0442 (third degree), otherwise `unrelated` (a value exactly equal to a threshold falls in the lower class; the thresholds are 2^-1.5, 2^-2.5, 2^-3.5 and 2^-4.5 written to four decimals). No IBD segments are inferred, no linkage or recombination model is used, it uses genotype counts only (no allele frequencies), and the estimate needs many independent, common SNPs: with few sites (as in the example) the values are noisy and can be negative. At most 150 samples.

- **API:** `dnakit.variants.calculate_kinship_king(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a multi-sample VCF (text or path) or records with diploid genotypes; only biallelic sites are used. `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method` and `pairs`, one entry per sample pair with `sample_1`, `sample_2`, `n_sites`, `n_het_het`, `n_opposite_homozygous`, `ibs0_fraction`, `kinship`, `pi_hat_approx` and `relationship` (`kinship` is `null` when neither sample is heterozygous). Fewer than two samples, more than 150 samples or no biallelic site raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import calculate_kinship_king

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t0/1",
    "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/0\t0/1",
    "chr1\t300\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/0",
    "chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t1/1",
    "chr1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/1\t0/1",
    "chr1\t600\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t1/1\t0/0",
    "chr1\t700\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1",
    "chr1\t800\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1",
]) + "\n"

result = calculate_kinship_king(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_kinship_king",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 399,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7258ad9577e33ae92458e95ff141dfec5403d4cf8329ffb6e15f1b20c0d8b391"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "KING-robust kinship estimator (Manichaikul et al. 2010): (N_AaAa - 2 N_AAaa) / (N_Aa(i) +...",
    "pairs": [
      {
        "sample_1": "S1",
        "sample_2": "S2",
        "n_sites": 8,
        "n_het_het": 4,
        "n_opposite_homozygous": 0,
        "ibs0_fraction": 0.0,
        "kinship": 0.5,
        "pi_hat_approx": 1.0,
        "relationship": "duplicate_or_monozygotic_twin"
      },
      {
        "sample_1": "S1",
        "sample_2": "S3",
        "n_sites": 8,
        "n_het_het": 2,
        "n_opposite_homozygous": 2,
        "ibs0_fraction": 0.25,
        "kinship": -0.25,
        "pi_hat_approx": -0.5,
        "relationship": "unrelated"
      },
      "..."
    ]
  }
}
```

## 17) calculate_ld_matrix

- **Purpose:** Compute a pairwise linkage-disequilibrium matrix (r2, r or D') between biallelic variants.

- **Method:** With phased genotypes the haplotype frequencies of each variant pair give D = p_AB - p_A p_B, r = D / sqrt(p_A(1-p_A) p_B(1-p_B)), r2 = r^2 and D' = |D| / D_max. With unphased genotypes the Pearson correlation of ALT dosages (0/1/2) is used (Rogers-Huff r); D' is not available then and raises an error. Samples missing either genotype are dropped per pair; a pair with fewer than 2 usable samples or a monomorphic variant gives `null`. `phase: 'auto'` uses haplotypes only when every called heterozygous genotype is phased. Multi-allelic sites are skipped with a warning. This is the plain sample statistic with no small-sample correction.

- **API:** `dnakit.variants.calculate_ld_matrix(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text/path or records with sample genotypes (at least two samples and two biallelic variants). `options`: `statistic` (`'r2'` default, `'r'`, `'d_prime'`), `phase` (`'auto'` default, `'phased'`, `'unphased'`), `max_variants` (default 200, at most 1000).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `statistic`, `method` (`haplotype-based (phased)` or the dosage-correlation description), `variants` (site names such as `chr1:100:A:G`) and `matrix` (a symmetric square matrix, `null` where undefined).

- **Example code:**

```python
import json

from dnakit.variants import calculate_ld_matrix

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "chr1\t110\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "chr1\t120\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t0/1",
    "chr1\t130\t.\tT\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/1\t0/1",
]) + "\n"

result = calculate_ld_matrix(vcf)  # unphased genotypes: r^2 of ALT dosages

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_ld_matrix",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 239,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ce2023f5d41f80b38428eb4e237b65db249da7073fb2981324bde0150924df8c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "r2",
    "method": "squared/plain Pearson correlation of ALT-allele dosages (unphased, Rogers-Huff style r)",
    "variants": ["chr1:100:A:G", "chr1:110:C:T", "..."],
    "matrix": [[1.0, 1.0, "..."], [1.0, 1.0, "..."], "..."]
  }
}
```

## 18) calculate_minor_allele_frequency

- **Purpose:** Document a reusable calculate minor allele frequency capability for Variant analysis workflows. Candidate CAP-0776 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate minor allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_minor_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. Explicit `genotypes` are required; no missing-genotype imputation is performed.

- **Output:** The local `calculate_minor_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_minor_allele_frequency

result = calculate_minor_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0776",
  "operation": "calculate_minor_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_minor_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'called_alleles': 6, 'allele_frequencies': {'0': 0.5, '1': 0.5}, 'minor_allele_index': 0, 'minor_allele': 'A', 'maf': 0.5}], 'variant_count': 1, ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 19) calculate_missingness

- **Purpose:** Document a reusable calculate missingness capability for Variant analysis workflows. Candidate CAP-0777 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate missingness operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_missingness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_missingness` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_missingness

result = calculate_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0777",
  "operation": "calculate_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}, 'provenance': {'implementation'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 20) calculate_runs_of_homozygosity

- **Purpose:** Find runs of homozygosity in one sample from consecutive homozygous genotypes.

- **Method:** A greedy left-to-right scan per chromosome: a run starts at a homozygous call and extends while heterozygous or missing calls stay within `max_heterozygotes` and `max_missing` and the gap between consecutive variants stays within `max_gap_bp`; trailing heterozygous/missing calls are trimmed. Runs with fewer than `min_snps` calls or shorter than `min_length_bp` are discarded. It is a consecutive-SNP heuristic in the style of PLINK `--homozyg`, not an HMM, and uses no genetic map, allele frequencies or genotype-error model. The defaults (25 SNPs, 1 Mb) are meant for dense genome-scale data.

- **API:** `dnakit.variants.calculate_runs_of_homozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text/path or records with genotypes. `options`: `sample` (default the first sample), `min_snps` (25), `min_length_bp` (1000000), `max_heterozygotes` (1), `max_missing` (2), `max_gap_bp` (1000000).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `sample`, `method`, `parameters`, `runs` (each with `chrom`, `start`, `end`, `length_bp`, `n_snps`, `n_heterozygous`, `n_missing`), `n_runs` and `total_roh_bp`.

- **Example code:**

```python
import json

from dnakit.variants import calculate_runs_of_homozygosity

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t1000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t2000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t3000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t4000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t5000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t6000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t7000\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t8000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t9000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t10000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
]) + "\n"

# the defaults (25 SNPs, 1 Mb) are for genome-scale data; relaxed here for the toy example
result = calculate_runs_of_homozygosity(vcf, options={"min_snps": 5, "min_length_bp": 3000})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "calculate_runs_of_homozygosity",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 12,
    "text_length": 361,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "13481682ccc52cdb8c6279ba8202488902c4a2975beda2550094f35953e2d450"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sample": "S1",
    "method": "greedy left-to-right runs of homozygous genotypes tolerating a bounded number of heterozyg...",
    "parameters": {"min_snps": 5, "min_length_bp": 3000, "max_heterozygotes": 1, "max_missing": 2, "max_gap_bp": 1000000},
    "runs": [
      {
        "chrom": "chr1",
        "start": 1000,
        "end": 10000,
        "length_bp": 9001,
        "n_snps": 10,
        "n_heterozygous": 1,
        "n_missing": 0
      }
    ],
    "n_runs": 1,
    "total_roh_bp": 9001
  }
}
```

## 21) calculate_transition_transversion

- **Purpose:** Document a reusable calculate transition transversion capability for Variant analysis workflows. Candidate CAP-0779 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate transition transversion operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_transition_transversion(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_transition_transversion` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0779",
  "operation": "calculate_transition_transversion",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_transition_transversion', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'classification': 'transition'}], 'transitions': 1, 'transversions': 0, 'excluded': 0, 'ti_tv_ratio': None}, 'provenance': {'implementation': 'dna...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 22) calculate_variant_density

- **Purpose:** Document a reusable calculate variant density capability for Variant analysis workflows. Candidate CAP-0780 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate variant density operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_variant_density(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. `options.region_length` is required for variant density.

- **Output:** The local `calculate_variant_density` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_variant_density

result = calculate_variant_density(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"region_length": 100}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0780",
  "operation": "calculate_variant_density",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_variant_density', 'result': {'variant_count': 1, 'region_length': 100, 'density_per_base': 0.01, 'density_per_megabase': 10000.0, 'by_chromosome': {'chr1': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_d...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 23) call

- **Purpose:** Call variants from sequencing reads. Not implemented: it only reports what is missing (an external variant caller and aligned reads).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Use `call_freebayes`, `call_gatk_haplotypecaller`, `call_deepvariant`, `call_delly`, `call_cutesv` or `call_sniffles2`, which run the external program when it is installed.

- **API:** `dnakit.variants.call(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("variant caller binary (bcftools/FreeBayes/GATK/DeepVariant) and an aligned BAM/CRAM with a reference") and `message` ("variant calling needs an external caller and aligned reads; use the call_<caller> functions."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import call

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = call(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "call",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "variant caller binary (bcftools/FreeBayes/GATK/DeepVariant) and an aligned BAM/CRAM with a reference",
  "message": "variant calling needs an external caller and aligned reads; use the call_<caller> functions."
}
```

## 24) call_mt

- **Purpose:** Call mitochondrial variants from reads. Not implemented: it only reports what is missing (a mitochondria-aware caller (Mutect2 mitochondria mode, mutserve) and aligned reads).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No caller is run.

- **API:** `dnakit.variants.call_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("mitochondrial variant caller (Mutect2 --mitochondria-mode, mutserve) and aligned reads") and `message` ("mitochondrial variant calling needs an external caller and aligned reads."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import call_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = call_mt(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "call_mt",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "mitochondrial variant caller (Mutect2 --mitochondria-mode, mutserve) and aligned reads",
  "message": "mitochondrial variant calling needs an external caller and aligned reads."
}
```

## 25) compare

- **Purpose:** Document a reusable compare capability for Variant analysis workflows. Candidate CAP-0783 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `generic_operation_name`. Existing-name evidence: `compare|dnakit.similarity.compare`.

- **API:** `dnakit.variants.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping with `left` and `right` variant collections; matching uses exact CHROM/POS/REF/ALT identity without normalization.

- **Output:** The local `compare` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import compare

result = compare(
    input={"left": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]}, "right": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"}]}}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0783",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'compare', 'result': {'left_record_count': 1, 'right_record_count': 2, 'left_allele_count': 1, 'right_allele_count': 2, 'shared_allele_count': 1, 'only_left': [], 'only_right': [{'chrom': 'chr1', 'pos': 20, 'ref': 'C', 'alt': 'T'}], 'shared': [{'chrom': 'chr1'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 26) compare_mt

- **Purpose:** Compare two variant call sets: shared and private alleles, precision, recall, F1, Jaccard and genotype concordance.

- **Method:** `compare_mt` is the variant-by-sample (matrix-table style) name of `compare_vcfs`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. The left set is the truth/reference set and the right set is the query. Alleles are matched per ALT allele on (CHROM, POS, REF, ALT) (`match: 'exact'`) or on (CHROM, POS) (`'position'`). precision = shared / |right|, recall = shared / |left|, F1 is their harmonic mean and Jaccard = shared / |union|. For samples present in both sets, genotype concordance is computed over the shared sites (`exact` matching only), with alleles compared by sequence and order-insensitively; `non_reference_concordance` restricts this to sites where either genotype carries a non-reference allele. Without a reference, equivalent indel representations at different positions are different alleles; with `options.reference` both sets are normalised (trimmed and left-aligned) first. It is a set comparison, not a haplotype-aware tool such as vcfeval or hap.py.

- **API:** `dnakit.variants.compare_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'left': <VCF|records>, 'right': <VCF|records>}`, or one call set with the other in `options.other`. `options`: `other`, `match` (`'exact'` default or `'position'`), `reference` (plain sequence, FASTA text/path or `{contig: sequence}`) with `reference_start` (default 1) to normalise both sets first, `check_ref` (`'warn'` default or `'error'`).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `match`, `left_count`, `right_count`, `shared_count`, `only_left_count`, `only_right_count`, `precision_right_vs_left`, `recall_right_vs_left`, `f1`, `jaccard`, the lists `shared`, `only_left`, `only_right`, `genotype_concordance` (per shared sample: `compared_genotypes`, `concordant`, `discordant`, `genotype_concordance`, `non_reference_concordance`) and `normalized_with_reference`. A single call set without `options.other` raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import compare_mt

left = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"
right = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t1/1",
    "chr1\t400\t.\tG\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"

result = compare_mt({"left": left, "right": right})  # left = truth set

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "compare_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 417,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d82b05108a1b36c57a3f1642b1b8a7f45cb9b0ba3b49730adc578adc27802a35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "match": "exact",
    "left_count": 3,
    "right_count": 3,
    "shared_count": 2,
    "only_left_count": 1,
    "only_right_count": 1,
    "precision_right_vs_left": 0.666667,
    "recall_right_vs_left": 0.666667,
    "f1": 0.666667,
    "jaccard": 0.5,
    "shared": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 300, "ref": "AT", "alt": "A"}],
    "only_left": [{"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"}],
    "only_right": [{"chrom": "chr1", "pos": 400, "ref": "G", "alt": "A"}],
    "genotype_concordance": {
      "S1": {
        "compared_genotypes": 2,
        "concordant": 1,
        "discordant": 1,
        "genotype_concordance": 0.5,
        "non_reference_concordance": 0.5
      }
    },
    "normalized_with_reference": false
  }
}
```

## 27) compare_vcfs

- **Purpose:** Compare two variant call sets: shared and private alleles, precision, recall, F1, Jaccard and genotype concordance.

- **Method:** The left set is the truth/reference set and the right set is the query. Alleles are matched per ALT allele on (CHROM, POS, REF, ALT) (`match: 'exact'`) or on (CHROM, POS) (`'position'`). precision = shared / |right|, recall = shared / |left|, F1 is their harmonic mean and Jaccard = shared / |union|. For samples present in both sets, genotype concordance is computed over the shared sites (`exact` matching only), with alleles compared by sequence and order-insensitively; `non_reference_concordance` restricts this to sites where either genotype carries a non-reference allele. Without a reference, equivalent indel representations at different positions are different alleles; with `options.reference` both sets are normalised (trimmed and left-aligned) first. It is a set comparison, not a haplotype-aware tool such as vcfeval or hap.py.

- **API:** `dnakit.variants.compare_vcfs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'left': <VCF|records>, 'right': <VCF|records>}`, or one call set with the other in `options.other`. `options`: `other`, `match` (`'exact'` default or `'position'`), `reference` (plain sequence, FASTA text/path or `{contig: sequence}`) with `reference_start` (default 1) to normalise both sets first, `check_ref` (`'warn'` default or `'error'`).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `match`, `left_count`, `right_count`, `shared_count`, `only_left_count`, `only_right_count`, `precision_right_vs_left`, `recall_right_vs_left`, `f1`, `jaccard`, the lists `shared`, `only_left`, `only_right`, `genotype_concordance` (per shared sample: `compared_genotypes`, `concordant`, `discordant`, `genotype_concordance`, `non_reference_concordance`) and `normalized_with_reference`. A single call set without `options.other` raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import compare_vcfs

left = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"
right = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t1/1",
    "chr1\t400\t.\tG\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"

result = compare_vcfs({"left": left, "right": right})  # left = truth set

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "compare_vcfs",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 417,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d82b05108a1b36c57a3f1642b1b8a7f45cb9b0ba3b49730adc578adc27802a35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "match": "exact",
    "left_count": 3,
    "right_count": 3,
    "shared_count": 2,
    "only_left_count": 1,
    "only_right_count": 1,
    "precision_right_vs_left": 0.666667,
    "recall_right_vs_left": 0.666667,
    "f1": 0.666667,
    "jaccard": 0.5,
    "shared": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 300, "ref": "AT", "alt": "A"}],
    "only_left": [{"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"}],
    "only_right": [{"chrom": "chr1", "pos": 400, "ref": "G", "alt": "A"}],
    "genotype_concordance": {
      "S1": {
        "compared_genotypes": 2,
        "concordant": 1,
        "discordant": 1,
        "genotype_concordance": 0.5,
        "non_reference_concordance": 0.5
      }
    },
    "normalized_with_reference": false
  }
}
```

## 28) construct_haplotypes

- **Purpose:** Construct per-sample haplotypes (allele lists and, with a reference, sequences) from phased genotypes.

- **Method:** The first allele of each phased GT (`0|1`) is haplotype 1 and the second is haplotype 2. For every sample and chromosome the alleles are listed in position order. With `options.reference` the REF bases over the window are replaced by each haplotype's alleles to give a haplotype sequence (REF must match the reference, overlapping variants are an error, symbolic alleles are ignored with a warning). This function does not phase anything (see `phase_genotypes`): an unphased heterozygous genotype is an error, or skipped with `on_unphased: 'skip'`; missing and non-diploid calls are skipped with a warning.

- **API:** `dnakit.variants.construct_haplotypes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`); the genotypes should be phased. `options`: `sample`/`samples` (default all samples), `on_unphased` (`'error'` default or `'skip'`), `reference` (plain sequence, FASTA text/path or `{contig: sequence}`) with `reference_start` and `window_start`/`window_end` (default: the span of the variants) to also build sequences.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `haplotypes` (`{sample: {chrom: {'1': {...}, '2': {...}}}}`, each haplotype with `alleles` as a list of `{pos, ref, allele}` and, with a reference, `sequence`) and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import construct_haplotypes

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1",
    "chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t1|0\t0|0",
    "chr1\t8\t.\tT\tC\t.\tPASS\t.\tGT\t0|0\t0|1",
]) + "\n"

# reference positions: A1 C2 G3 T4 A5 C6 G7 T8 A9 C10; sequences are built over 2..8
result = construct_haplotypes(vcf, options={"reference": "ACGTACGTAC", "sample": "S1"})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "construct_haplotypes",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 172,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "94af7e21244e9ed9a82bd69757fe1112631c1f2ae81aadf937d83f5c4e4f8377"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-sample haplotype alleles taken from the phased GT (first allele = haplotype 1); option...",
    "haplotypes": {
      "S1": {
        "chr1": {
          "1": {
            "alleles": [
              {"pos": 2, "ref": "C", "allele": "C"},
              {"pos": 5, "ref": "A", "allele": "G"},
              {"pos": 8, "ref": "T", "allele": "T"}
            ],
            "sequence": "CGTGCGT"
          },
          "2": {
            "alleles": [
              {"pos": 2, "ref": "C", "allele": "T"},
              {"pos": 5, "ref": "A", "allele": "A"},
              {"pos": 8, "ref": "T", "allele": "T"}
            ],
            "sequence": "TGTACGT"
          }
        }
      }
    }
  }
}
```

## 29) count

- **Purpose:** Document a reusable count capability for Variant analysis workflows. Candidate CAP-0787 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.count(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a sequence or mapping containing VCF-like records; counts are based on explicit CHROM/POS/REF/ALT fields.

- **Output:** The local `count` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import count

result = count(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0787",
  "operation": "count",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'count', 'result': {'record_count': 1, 'alternate_allele_count': 1, 'by_chromosome': {'chr1': 1}, 'allele_type_counts': {'snp': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_dependency_free', 'coordinate_system': '...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 30) count_mt

- **Purpose:** Count the records of the input (non-blank lines of text, or items of a list). Despite the name it does not restrict the count to mitochondrial variants and does not read VCF structure.

- **Method:** Counts records in the shared local runtime (`run_candidate_operation`, helper `_local_result` in `src/dnakit/_candidate_runtime.py`): a text (or bytes, or the content of a `Path`) is split into non-blank lines and each line is one record; a list counts its items; a mapping that has one of the keys `records`, `items`, `values`, `data`, `sequences` or `reads` counts that value in the same way, any other mapping is one record. Nothing is parsed or validated: VCF header lines (`##...`, `#CHROM`) and records on other contigs are counted like any other line, so a VCF with a header of two lines and one variant gives 3. `None` counts as one record and an empty string as zero. For a VCF-aware count with per-chromosome and variant-class totals from explicit CHROM/POS/REF/ALT records use `dnakit.variants.count`.

- **API:** `dnakit.variants.count_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text, `bytes`, a `Path` of an existing file (first 1,000,000 bytes), a list of records, or a mapping (see Method). `options` is accepted but unused.

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"variants"`), `operation` (`"count_mt"`), `provenance` (`implementation` `"local_candidate_adapter"`, `backend`), `input_summary` (generic text summary: `record_count` = non-blank lines, `text_length`, `sequence_bases` = number of letters, `lengths`, `sha256`), `warnings` (empty) and `result`. `result` is an integer, the number of records.

- **Example code:**

```python
import json
from dnakit.variants import count_mt

result = count_mt("chrM\t73\t.\tA\tG\nchrM\t263\t.\tA\tG\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "count_mt",
  "provenance": {
    "implementation": "local_candidate_adapter",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "af190e6e8fb98b032e563306294f3551e543d377ca9a2b055fa3af98f76693e4"
  },
  "warnings": [],
  "status": "ok",
  "result": 2
}
```

## 31) detect_breakpoints

- **Purpose:** Detect structural-variant breakpoints from alignments. Not implemented: it only reports what is missing (split/discordant read alignments and an SV caller (DELLY, Manta, GRIDSS)).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. `call_delly`, `call_cutesv` and `call_sniffles2` run such tools when they are installed.

- **API:** `dnakit.variants.detect_breakpoints(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("SV caller (DELLY/Manta/GRIDSS) and aligned reads") and `message` ("breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and call_sniffles2 run such tools."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import detect_breakpoints

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = detect_breakpoints(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "detect_breakpoints",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "SV caller (DELLY/Manta/GRIDSS) and aligned reads",
  "message": "breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and c..."
}
```

## 32) detect_low_frequency_variants

- **Purpose:** Detect low-frequency variants from per-sample allelic depths with a binomial error-rate test.

- **Method:** For every sample and ALT allele (multi-allelic sites are split) with `FORMAT/AD`, a call is reported when `min_vaf <= alt/depth <= max_vaf`, `depth >= min_depth`, `alt >= min_alt_count` and the one-sided binomial P(X >= alt | depth, error_rate) <= alpha. The error rate is one value for the whole run and no strand-bias, sequence-context, base-quality or UMI model is used, so it is a heuristic screen and not a validated low-VAF caller; set `error_rate` from the assay's measured background.

- **API:** `dnakit.variants.detect_low_frequency_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text/path or records with `FORMAT/AD` (ref,alt) per sample. `options`: `min_vaf` (0.01), `max_vaf` (0.2), `min_depth` (100), `min_alt_count` (3), `error_rate` (0.001), `alpha` (0.01), `samples`/`sample` (default all).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `n_tested` (sample-site pairs with enough depth), `n_detected` and `variants` (each with `site`, `chrom`, `pos`, `ref`, `alt`, `sample`, `alt_count`, `depth`, `vaf`, `p_value_vs_error`). `warnings` notes when no call has AD.

- **Example code:**

```python
import json

from dnakit.variants import detect_low_frequency_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/0:995,5\t0/1:50,50",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/0:1000,0\t0/0:998,2",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/0:300,12\t0/0:40,1",
]) + "\n"

result = detect_low_frequency_variants(vcf)  # defaults: depth >= 100, 1% <= VAF <= 20%, >= 3 ALT reads

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "detect_low_frequency_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 224,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4fcd1dcadbc0bb575ea2cfc614e6d321fe5cadfa4d7763786d68656459193fe9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "calls where min_vaf <= alt/depth <= max_vaf, alt count >= min_alt_count, depth >= min_dept...",
    "parameters": {"min_vaf": 0.01, "max_vaf": 0.2, "min_depth": 100, "min_alt_count": 3, "error_rate": 0.001, "alpha": 0.01},
    "n_tested": 5,
    "n_detected": 1,
    "variants": [
      {
        "site": "chr1:300:G:A",
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": "A",
        "sample": "S1",
        "alt_count": 12,
        "depth": 312,
        "vaf": 0.0384615,
        "p_value_vs_error": 1.08712e-15
      }
    ]
  }
}
```

## 33) detect_mobile_element_insertions

- **Purpose:** Select mobile-element insertion records (Alu, LINE-1, SVA, HERV) from a structural-variant VCF.

- **Method:** A record is a mobile-element insertion when its ALT is a symbolic allele such as `<INS:ME:ALU>`, or `INFO/SVTYPE` is a known element (`ALU`, `L1`, `LINE1`, `SVA`, `HERV`, `HERVK`), or `INFO/MEINFO` is present; the family is normalised (for example `L1` becomes `LINE-1`). It only selects calls already made by an SV/MEI caller (MELT, Scramble and similar); it does not look at alignments or search for insertions.

- **API:** `dnakit.variants.detect_mobile_element_insertions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: SV VCF text/path or records with symbolic ALT or `INFO/SVTYPE`/`MEINFO`. `options`: `families` (optional list such as `['Alu', 'LINE-1']`; default all).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `n_records`, `n_insertions` and `insertions` (each with `chrom`, `pos`, `id`, `alt`, `family`, `svlen`, `filter`, `qual`).

- **Example code:**

```python
import json

from dnakit.variants import detect_mobile_element_insertions

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tsv1\tN\t<INS:ME:ALU>\t60\tPASS\tSVTYPE=INS;SVLEN=300",
    "chr1\t2000\tsv2\tN\t<INS:ME:L1>\t50\tPASS\tSVTYPE=INS;SVLEN=1500",
    "chr1\t3000\tsv3\tN\t<DEL>\t40\tPASS\tSVTYPE=DEL;SVLEN=-500",
]) + "\n"

result = detect_mobile_element_insertions(vcf, options={"families": ["Alu"]})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "detect_mobile_element_insertions",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 228,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "23103667ee2bb74d5020ced8f2e18261c7bf4e483df3da110656b25803022793"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "selection of mobile-element-insertion records from an SV VCF: ALT <INS:ME:FAMILY>, INFO/SV...",
    "n_records": 3,
    "n_insertions": 1,
    "insertions": [
      {
        "chrom": "chr1",
        "pos": 1000,
        "id": "sv1",
        "alt": "<INS:ME:ALU>",
        "family": "Alu",
        "svlen": "300",
        "filter": "PASS",
        "qual": 60.0
      }
    ]
  }
}
```

## 34) detect_mosaic_variants

- **Purpose:** Detect candidate mosaic variants: ALT fractions significantly above sequencing error and significantly below a germline 0.5.

- **Method:** For every sample and ALT allele with `FORMAT/AD`, a call is reported when `min_vaf <= alt/depth <= max_vaf`, `depth >= min_depth`, `alt >= min_alt_count`, the one-sided upper-tail binomial P(X >= alt | depth, error_rate) <= alpha and the lower-tail P(X <= alt | depth, 0.5) <= alpha. It is a heuristic: no population-frequency, strand-bias, panel-of-normals or parental-genotype filters are applied, and the single `error_rate` is not measured from the data.

- **API:** `dnakit.variants.detect_mosaic_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text/path or records with `FORMAT/AD` (ref,alt) per sample. `options`: `min_vaf` (0.03), `max_vaf` (0.4), `min_depth` (20), `min_alt_count` (3), `error_rate` (0.002), `alpha` (0.01), `samples`/`sample` (default all).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `n_tested`, `n_detected` and `variants` (each with `site`, `chrom`, `pos`, `ref`, `alt`, `sample`, `alt_count`, `depth`, `vaf`, `p_value_vs_error`, `p_value_vs_heterozygous`).

- **Example code:**

```python
import json

from dnakit.variants import detect_mosaic_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:80,20\t0/1:50,50",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/0:99,1\t0/0:98,2",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/1:45,55\t0/1:60,40",
]) + "\n"

result = detect_mosaic_variants(vcf)  # defaults: depth >= 20, 3% <= VAF <= 40%, >= 3 ALT reads

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "detect_mosaic_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 221,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "62da7281e20b1bd52f80c4f9a982472fe23027aa997553fe9c0cc010bfbff1e7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "mosaic candidates: ALT fraction in [min_vaf, max_vaf], enough ALT reads, and the ALT count...",
    "parameters": {"min_vaf": 0.03, "max_vaf": 0.4, "min_depth": 20, "min_alt_count": 3, "error_rate": 0.002, "alpha": 0.01},
    "n_tested": 6,
    "n_detected": 1,
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "sample": "S1",
        "alt_count": 20,
        "depth": 100,
        "vaf": 0.2,
        "p_value_vs_error": 4.82527e-34,
        "p_value_vs_heterozygous": 5.57954e-10
      }
    ]
  }
}
```

## 35) detect_repeat_expansions

- **Purpose:** Find the longest uninterrupted tandem repeat of a given unit in DNA sequences and flag expansions above a threshold.

- **Method:** For each sequence and motif, the longest run of exact, uninterrupted copies of the motif is found with a regular expression (any phase; one interruption ends a run). `expansion` is true when the run has at least `expansion_threshold` units. With `both_strands` the reverse complement of the motif is also searched. It reads a locus sequence you supply; it does not genotype repeats from sequencing reads like ExpansionHunter and does not tolerate impure repeats.

- **API:** `dnakit.variants.detect_repeat_expansions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a DNA string, FASTA text, a list of strings or `{name: sequence}` (each sequence at most 50 Mb). `options`: `motif` (required; 1-12 bases of ACGT, or a list), `expansion_threshold` (units; optional), `both_strands` (default False).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `expansion_threshold` and `results`, one entry per sequence, motif and strand with `sequence`, `motif`, `strand`, `longest_pure_run_units`, `start`, `end` (1-based; `null` when the motif never occurs) and `expansion` (`null` when no threshold is given). A missing or invalid `motif` raises `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.variants import detect_repeat_expansions

seq = "GGC" + "CAG" * 12 + "CAA" + "CAG" * 5 + "TTT"

result = detect_repeat_expansions(seq, options={"motif": "CAG", "expansion_threshold": 10})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "detect_repeat_expansions",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 60,
    "sequence_bases": 60,
    "lengths": [60],
    "sha256": "ffe9dfbc3a65a2a5108399d186ed70cdd7013fe8582ec07b76f6992253373167"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "longest uninterrupted tandem run of the repeat unit in each sequence (exact matches, every...",
    "expansion_threshold": 10,
    "results": [
      {
        "sequence": "sequence_1",
        "motif": "CAG",
        "strand": "+",
        "longest_pure_run_units": 12,
        "start": 4,
        "end": 39,
        "expansion": true
      }
    ]
  }
}
```

## 36) filter

- **Purpose:** Filter variant records (site and genotype level) with quality, depth, missingness, allele-frequency, HWE, variant-type and region criteria.

- **Method:** Every record is tested against each criterion that is set; a record is removed (or flagged) when any test fails, and the failing reasons are counted (`not_pass`, `low_qual`, `low_depth`, `high_depth`, `high_missing`, `af_out_of_range`, `hwe`, `variant_type`, `chrom`, `region`). A missing `QUAL` or depth fails the corresponding threshold, and `pass_only` accepts `PASS` or `.`. Genotype-level thresholds do not remove records but set the failing calls to missing (a call without `GQ`/`DP` fails the corresponding threshold too), and the missingness, allele-frequency and HWE tests see the genotypes after that masking. HWE uses the exact test of `calculate_hwe_pvalue`. The result can be written as a VCF with `output`.

- **API:** `dnakit.variants.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options` (all optional): `pass_only`, `min_qual`, `min_depth`, `max_depth` (site `DP`), `min_genotype_quality` and `min_sample_depth` (calls below the threshold are set to `./.`), `max_missing` (fraction), `min_af`, `max_af` (ALT frequency computed from the genotypes, else `INFO/AF`), `min_hwe_pvalue`, `variant_types` (`SNV`, `MNV`, `indel`, `SV`), `chroms`, `regions` (`'chr:start-end'`), `mode` (`'remove'` default, or `'flag'` to write the failing reasons into `FILTER` instead of dropping the record), `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf` (the filtered VCF text), `output_path`, `records` (unless `include_records` is False; per record `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`, `info`, `genotypes`), `input_records`, `output_records`, `failed_records`, `reason_counts`, `genotypes_masked` and `mode`.

- **Example code:**

```python
import json

from dnakit.variants import filter

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = filter(vcf, options={"pass_only": True, "min_qual": 30, "max_missing": 0.2})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "filter",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2", "S3"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {"DP": "30"},
        "genotypes": {
          "S1": {"GT": "0/0", "DP": "10", "GQ": "50"},
          "S2": {"GT": "0/1", "DP": "10", "GQ": "40"},
          "S3": {"GT": "1/1", "DP": "10", "GQ": "30"}
        }
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "AT",
        "alt": ["A"],
        "qual": 60.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A", "C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/2"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "input_records": 4,
    "output_records": 3,
    "failed_records": 1,
    "reason_counts": {"high_missing": 1, "low_qual": 1, "not_pass": 1},
    "genotypes_masked": 0,
    "mode": "remove"
  }
}
```

## 37) filter_mt

- **Purpose:** Filter a variant-by-sample (matrix-table style) call set with quality, depth, missingness, allele-frequency, HWE, type and region criteria.

- **Method:** `filter_mt` is the variant-by-sample (matrix-table style) name of `filter`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Every record is tested against each criterion that is set; a record is removed (or flagged) when any test fails, and the failing reasons are counted (`not_pass`, `low_qual`, `low_depth`, `high_depth`, `high_missing`, `af_out_of_range`, `hwe`, `variant_type`, `chrom`, `region`). A missing `QUAL` or depth fails the corresponding threshold, and `pass_only` accepts `PASS` or `.`. Genotype-level thresholds do not remove records but set the failing calls to missing (a call without `GQ`/`DP` fails the corresponding threshold too), and the missingness, allele-frequency and HWE tests see the genotypes after that masking. HWE uses the exact test of `calculate_hwe_pvalue`. The result can be written as a VCF with `output`.

- **API:** `dnakit.variants.filter_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options` (all optional): `pass_only`, `min_qual`, `min_depth`, `max_depth` (site `DP`), `min_genotype_quality` and `min_sample_depth` (calls below the threshold are set to `./.`), `max_missing` (fraction), `min_af`, `max_af` (ALT frequency computed from the genotypes, else `INFO/AF`), `min_hwe_pvalue`, `variant_types` (`SNV`, `MNV`, `indel`, `SV`), `chroms`, `regions` (`'chr:start-end'`), `mode` (`'remove'` default, or `'flag'` to write the failing reasons into `FILTER` instead of dropping the record), `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf` (the filtered VCF text), `output_path`, `records` (unless `include_records` is False; per record `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`, `info`, `genotypes`), `input_records`, `output_records`, `failed_records`, `reason_counts`, `genotypes_masked` and `mode`.

- **Example code:**

```python
import json

from dnakit.variants import filter_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = filter_mt(vcf, options={"pass_only": True, "min_qual": 30, "max_missing": 0.2})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "filter_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2", "S3"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {"DP": "30"},
        "genotypes": {
          "S1": {"GT": "0/0", "DP": "10", "GQ": "50"},
          "S2": {"GT": "0/1", "DP": "10", "GQ": "40"},
          "S3": {"GT": "1/1", "DP": "10", "GQ": "30"}
        }
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "AT",
        "alt": ["A"],
        "qual": 60.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A", "C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/2"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "input_records": 4,
    "output_records": 3,
    "failed_records": 1,
    "reason_counts": {"high_missing": 1, "low_qual": 1, "not_pass": 1},
    "genotypes_masked": 0,
    "mode": "remove"
  }
}
```

## 38) filter_variants_by_allele_frequency

- **Purpose:** Document a reusable filter variants by allele frequency capability for Variant analysis workflows. Candidate CAP-0796 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. Optional `options.min_frequency`, `max_frequency`, and `allele_index` select the explicit frequency range.

- **Output:** The local `filter_variants_by_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_allele_frequency

result = filter_variants_by_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"min_frequency": 0.1, "max_frequency": 0.9}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0796",
  "operation": "filter_variants_by_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/1', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'allele_index': 1, 'frequency': 0.5, 'kept': True}], 'input_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 39) filter_variants_by_depth

- **Purpose:** Document a reusable filter variants by depth capability for Variant analysis workflows. Candidate CAP-0797 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by depth operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_depth(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and an explicit `DP`/`depth` value. `options.min_depth` is required; optional `max_depth` and `keep_missing` bound explicit DP values.

- **Output:** The local `filter_variants_by_depth` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_depth

result = filter_variants_by_depth(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "DP": 30}]},
    options={"min_depth": 10, "max_depth": 100}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0797",
  "operation": "filter_variants_by_depth",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_depth', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'DP': 30}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'depth': 30.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_record_coun...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 40) filter_variants_by_genotype_quality

- **Purpose:** Document a reusable filter variants by genotype quality capability for Variant analysis workflows. Candidate CAP-0798 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by genotype quality operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_genotype_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires explicit `GQ`, `genotype_quality`, `genotype_qualities`, or sample GQ values. `options.min_genotype_quality` or `min_gq` is required; GQ values are taken from explicit records or samples.

- **Output:** The local `filter_variants_by_genotype_quality` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_genotype_quality

result = filter_variants_by_genotype_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "GQ": 45}]},
    options={"min_genotype_quality": 20}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0798",
  "operation": "filter_variants_by_genotype_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_genotype_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'GQ': 45}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'genotype_qualities': [45.0], 'minimum_observed': 45.0, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 41) filter_variants_by_missingness

- **Purpose:** Document a reusable filter variants by missingness capability for Variant analysis workflows. Candidate CAP-0799 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by missingness operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_missingness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and either explicit `genotypes` or `missing_count/sample_count`. `options.max_missingness` is required and must be between 0 and 1.

- **Output:** The local `filter_variants_by_missingness` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_missingness

result = filter_variants_by_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "./.", "0/0"]}]},
    options={"max_missingness": 0.5}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0799",
  "operation": "filter_variants_by_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', './.', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'missingness': 0.3333333333333333, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 42) filter_variants_by_quality

- **Purpose:** Document a reusable filter variants by quality capability for Variant analysis workflows. Candidate CAP-0800 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by quality operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and a finite `qual`/`QUAL` value. `options.min_quality` is required; missing QUAL values are excluded unless `keep_missing=True`.

- **Output:** The local `filter_variants_by_quality` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_quality

result = filter_variants_by_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50}]},
    options={"min_quality": 20}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0800",
  "operation": "filter_variants_by_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'qual': 50}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'quality': 50.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 43) genotype_structural_variants

- **Purpose:** Genotype known structural variants in a sample from its reads. Not implemented: it only reports what is missing (read alignments and an SV genotyper (SVTyper, Paragraph, GraphTyper)).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No genotyper is run.

- **API:** `dnakit.variants.genotype_structural_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("SV genotyper (SVTyper/Paragraph) and aligned reads") and `message` ("SV genotyping needs read alignments and an external genotyper."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.variants import genotype_structural_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = genotype_structural_variants(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "genotype_structural_variants",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "SV genotyper (SVTyper/Paragraph) and aligned reads",
  "message": "SV genotyping needs read alignments and an external genotyper."
}
```

## 44) fill_missing_genotypes_with_mode

- **Purpose:** Fill missing genotype calls with the most frequent called genotype at the same site (naive per-site modal-genotype filling).

- **Method:** For each record the called genotypes are counted (allele order ignored, so `1/0` counts as `0/1`); every missing call (`./.`) receives the most frequent genotype at that site, ties broken by the genotype string, and is marked with `FORMAT/IMPUTED=1`. Sites with no called genotype stay missing and are counted in `sites_without_any_call`. This is NOT genotype imputation in the sense of Beagle or IMPUTE5: there is no reference panel, no linkage disequilibrium, no haplotype model and no relatedness, so rare-variant carriers are systematically replaced by the majority genotype and the filled calls carry no uncertainty. It is only a quick way to remove missing values for tools that cannot handle them; use a real imputation program for analysis.

- **API:** `dnakit.variants.fill_missing_genotypes_with_mode(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`); sample genotypes are required. `options`: `output` (optional path to write the VCF), `overwrite` (default False).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `imputed_genotypes` (number of calls filled), `sites_without_any_call`, `records`, `vcf` (the filled VCF text) and `output_path`. Input without sample columns raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import fill_missing_genotypes_with_mode

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t./.\t0/1",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t./.\t./.\t0/0",
]) + "\n"

result = fill_missing_genotypes_with_mode(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "fill_missing_genotypes_with_mode",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 196,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "09e61722a5752bcc121f14eeac1a1019fd4d1f2abe1186409b50bdb06f7a06fb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "naive per-site modal-genotype imputation: each missing call receives the most frequent cal...",
    "imputed_genotypes": 3,
    "sites_without_any_call": 0,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "0/1", "IMPUTED": "1"}, "S3": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0", "IMPUTED": "1"}, "S2": {"GT": "0/0", "IMPUTED": "1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null
  }
}
```

## 45) left_align

- **Purpose:** Left-align and trim indels in variant records against a reference sequence.

- **Method:** Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. A reference is required: `options.reference` is looked up for every record, and a record whose REF disagrees with the reference is reported as a warning and kept as it is (or raises with `check_ref: 'error'`). The shift stops at the start of the supplied reference sequence. `hgvs_accession` additionally writes HGVS genomic strings, which are 3' shifted and use `dup` for tandem insertions as the HGVS standard requires.

- **API:** `dnakit.variants.left_align(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; required for left alignment), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. A missing `options.reference` raises `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.variants import left_align

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "left_align",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 46) left_align_mt

- **Purpose:** Left-align and trim indels in a variant-by-sample (matrix-table style) call set against a reference sequence.

- **Method:** `left_align_mt` is the variant-by-sample (matrix-table style) name of `left_align`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. A reference is required: `options.reference` is looked up for every record, and a record whose REF disagrees with the reference is reported as a warning and kept as it is (or raises with `check_ref: 'error'`). The shift stops at the start of the supplied reference sequence. `hgvs_accession` additionally writes HGVS genomic strings, which are 3' shifted and use `dup` for tandem insertions as the HGVS standard requires.

- **API:** `dnakit.variants.left_align_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; required for left alignment), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. A missing `options.reference` raises `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.variants import left_align_mt

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align_mt(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "left_align_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 47) left_align_variants

- **Purpose:** Left-align and trim indels in variant records against a reference sequence.

- **Method:** Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. A reference is required: `options.reference` is looked up for every record, and a record whose REF disagrees with the reference is reported as a warning and kept as it is (or raises with `check_ref: 'error'`). The shift stops at the start of the supplied reference sequence. `hgvs_accession` additionally writes HGVS genomic strings, which are 3' shifted and use `dup` for tandem insertions as the HGVS standard requires.

- **API:** `dnakit.variants.left_align_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; required for left alignment), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. A missing `options.reference` raises `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.variants import left_align_variants

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align_variants(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "left_align_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 48) merge

- **Purpose:** Merge several VCFs / call sets, combining samples (default) or concatenating records.

- **Method:** `mode: 'samples'` (the default, like `bcftools merge -m none`) combines records with identical CHROM/POS/REF/ALT into one record with the samples of all inputs side by side; sample names must be unique and samples absent from an input get `./.`. `mode: 'records'` (like `bcftools concat`) requires identical sample lists and keeps duplicate sites once. The merged QUAL is the maximum, FILTER is `PASS` only if every contributing record passes (otherwise the union of the failing labels), INFO comes from the first record that carries the key, and the output is sorted. FORMAT fields other than GT are carried per sample; no allele harmonisation or re-normalisation is performed, so differently written alleles stay separate records.

- **API:** `dnakit.variants.merge(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of at least two call sets (VCF text, VCF paths or `{'records': [...]}`) or `{'inputs': [...]}`. `options`: `mode` (`'samples'` default or `'records'`), `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf` (merged VCF text), `output_path`, `records` and `mode`. Fewer than two call sets, duplicate sample names (in `samples` mode) or different sample lists (in `records` mode) raise `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import merge

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "merge",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 49) merge_mt

- **Purpose:** Merge several variant-by-sample (matrix-table style) call sets, combining samples or concatenating records.

- **Method:** `merge_mt` is the variant-by-sample (matrix-table style) name of `merge`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. `mode: 'samples'` (the default, like `bcftools merge -m none`) combines records with identical CHROM/POS/REF/ALT into one record with the samples of all inputs side by side; sample names must be unique and samples absent from an input get `./.`. `mode: 'records'` (like `bcftools concat`) requires identical sample lists and keeps duplicate sites once. The merged QUAL is the maximum, FILTER is `PASS` only if every contributing record passes (otherwise the union of the failing labels), INFO comes from the first record that carries the key, and the output is sorted. FORMAT fields other than GT are carried per sample; no allele harmonisation or re-normalisation is performed, so differently written alleles stay separate records.

- **API:** `dnakit.variants.merge_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of at least two call sets (VCF text, VCF paths or `{'records': [...]}`) or `{'inputs': [...]}`. `options`: `mode` (`'samples'` default or `'records'`), `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf` (merged VCF text), `output_path`, `records` and `mode`. Fewer than two call sets, duplicate sample names (in `samples` mode) or different sample lists (in `records` mode) raise `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import merge_mt

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge_mt([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "merge_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 50) merge_vcf_files

- **Purpose:** Merge several VCF files, combining samples (default) or concatenating records.

- **Method:** `mode: 'samples'` (the default, like `bcftools merge -m none`) combines records with identical CHROM/POS/REF/ALT into one record with the samples of all inputs side by side; sample names must be unique and samples absent from an input get `./.`. `mode: 'records'` (like `bcftools concat`) requires identical sample lists and keeps duplicate sites once. The merged QUAL is the maximum, FILTER is `PASS` only if every contributing record passes (otherwise the union of the failing labels), INFO comes from the first record that carries the key, and the output is sorted. FORMAT fields other than GT are carried per sample; no allele harmonisation or re-normalisation is performed, so differently written alleles stay separate records.

- **API:** `dnakit.variants.merge_vcf_files(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of at least two call sets (VCF text, VCF paths or `{'records': [...]}`) or `{'inputs': [...]}`. `options`: `mode` (`'samples'` default or `'records'`), `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf` (merged VCF text), `output_path`, `records` and `mode`. Fewer than two call sets, duplicate sample names (in `samples` mode) or different sample lists (in `records` mode) raise `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.variants import merge_vcf_files

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge_vcf_files([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "merge_vcf_files",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 51) normalize

- **Purpose:** Normalise variant records: trim shared bases and, when a reference is given, left-align indels (bcftools norm style).

- **Method:** Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. Without `options.reference` only trimming is done (one anchor base is kept for indels) and `left_aligned` is reported as false. With a reference, a record whose REF does not match is warned about and left unchanged (`check_ref`), `split_multiallelic` additionally writes one record per ALT after normalisation, and `hgvs_accession` adds 3'-shifted HGVS genomic strings. The output is sorted unless `sort` is False. Only the reference you supply is used; no genome is downloaded.

- **API:** `dnakit.variants.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; optional here), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. Without a reference only trimming is done.

- **Example code:**

```python
import json

from dnakit.variants import normalize

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "normalize",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 52) normalize_mt

- **Purpose:** Normalise a variant-by-sample (matrix-table style) call set: trim shared bases and, with a reference, left-align indels.

- **Method:** `normalize_mt` is the variant-by-sample (matrix-table style) name of `normalize`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. Without `options.reference` only trimming is done (one anchor base is kept for indels) and `left_aligned` is reported as false. With a reference, a record whose REF does not match is warned about and left unchanged (`check_ref`), `split_multiallelic` additionally writes one record per ALT after normalisation, and `hgvs_accession` adds 3'-shifted HGVS genomic strings. The output is sorted unless `sort` is False. Only the reference you supply is used; no genome is downloaded.

- **API:** `dnakit.variants.normalize_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; optional here), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. Without a reference only trimming is done.

- **Example code:**

```python
import json

from dnakit.variants import normalize_mt

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize_mt(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "normalize_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 53) normalize_vcf

- **Purpose:** Normalise variant records: trim shared bases and, when a reference is given, left-align indels (bcftools norm style).

- **Method:** Algorithm of Tan, Abecasis & Kang (2015): while the alleles share their last base, drop it; when an allele becomes empty, extend all alleles to the left by the reference base; repeat until no more change, then drop shared leading bases. Multi-allelic records are processed jointly, so all alleles move together. Without `options.reference` only trimming is done (one anchor base is kept for indels) and `left_aligned` is reported as false. With a reference, a record whose REF does not match is warned about and left unchanged (`check_ref`), `split_multiallelic` additionally writes one record per ALT after normalisation, and `hgvs_accession` adds 3'-shifted HGVS genomic strings. The output is sorted unless `sort` is False. Only the reference you supply is used; no genome is downloaded.

- **API:** `dnakit.variants.normalize_vcf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `reference` (plain sequence, FASTA text/path or `{contig: sequence}`; optional here), `reference_start` (1-based position of a plain sequence's first base, default 1), `check_ref` (`'warn'` default: a record whose REF does not match is left unchanged; or `'error'`), `split_multiallelic` (default False), `sort` (default True), `hgvs_accession` (accession string or `{contig: accession}`) to add an HGVS `g.` string per ALT, `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `samples`, `vcf`, `output_path`, `records` (with `hgvs_g` when `hgvs_accession` is given), `changed_records`, `unchanged_or_skipped`, `skipped_records`, `left_aligned` (False when no reference was used) and `changes` (`chrom`, `from`, `to`). `warnings` lists reference mismatches. Without a reference only trimming is done.

- **Example code:**

```python
import json

from dnakit.variants import normalize_vcf

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize_vcf(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "variants",
  "operation": "normalize_vcf",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 54) parse_vcf_records

- **Purpose:** Document a reusable parse vcf records capability for Variant analysis workflows. Candidate CAP-0812 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse vcf records operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.parse_vcf_records(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be VCF text, UTF-8 bytes, a sequence of VCF lines, or a mapping containing `text`; mandatory VCF columns are validated.

- **Output:** The local `parse_vcf_records` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import parse_vcf_records

result = parse_vcf_records(
    input='##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\nchr1\t10\trs1\tA\tG,T\t50\tPASS\tDP=12\tGT\t0/1',
    options={"strict": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0812",
  "operation": "parse_vcf_records",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'parse_vcf_records', 'result': {'headers': ['##fileformat=VCFv4.3', '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\ts1'], 'samples': ['s1'], 'records': [{'chrom': 'chr1', 'pos': 10, 'id': 'rs1', 'ref': 'A', 'alt': ['G', 'T'], 'qual': 50.0, 'filter': 'P...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 55) phase_genotypes

- **Purpose:** Phase a child's heterozygous genotypes by transmission from the parents' genotypes (trio phasing).

- **Method:** For each heterozygous child genotype x/y, the two assignments (x from the father and y from the mother, or the reverse) are checked against the parents' alleles; when exactly one assignment is compatible the child's GT is rewritten as paternal|maternal. Homozygous child genotypes are simply marked phased; sites where both assignments fit (for example both parents heterozygous), where neither fits (Mendelian error) and sites with a missing child genotype stay unchanged. A missing parent genotype does not constrain its allele. This is deterministic phasing by transmission only; statistical phasing (Beagle, SHAPEIT) is not performed, so most sites in unrelated samples stay unphased.

- **API:** `dnakit.variants.phase_genotypes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`); it must contain the child, father and mother samples. `options`: `trio` = `{'child': ..., 'father': ..., 'mother': ...}` (or the keys `child`/`father`/`mother` directly; required), `output` (optional VCF path), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `trio`, `counts` (numbers of sites that were `phased`, `ambiguous`, `homozygous`, `mendelian_error` and `child_missing`; categories that did not occur are absent), `vcf` (the VCF text with the child's phased genotypes) and `output_path`. A missing `trio` or an unknown sample name raises `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.variants import phase_genotypes

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1",
    "chr1\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1",
]) + "\n"

result = phase_genotypes(vcf, options={"trio": {"child": "child", "father": "father", "mother": "mother"}})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "phase_genotypes",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 243,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "267f2e692be1c3ef4bfe87ba4f02cbf34e14b5201bf1f1c67af2e98903ab95bb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "trio phasing by transmission: the child's heterozygous genotype is written paternal|matern...",
    "trio": {"child": "child", "father": "father", "mother": "mother"},
    "counts": {"phased": 2, "ambiguous": 1, "homozygous": 1},
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother\nchr...",
    "output_path": null
  }
}
```

## 56) plot_density

- **Purpose:** Plot variant density along chromosomes as a binned SVG bar chart and return the bin counts.

- **Method:** Variants are counted in fixed-size bins of `bin_size` bases per chromosome (bin index = (POS - 1) // bin_size); chromosomes are ordered naturally and, when `chrom_lengths` is given, empty trailing bins up to the chromosome length are included. Density per Mb = count * 1,000,000 / bin_size. The SVG is a plain bar chart (bar height proportional to the count, scaled to the maximum bin) with chromosome labels, built as text; no plotting library is needed. More than 20,000 bins is an error (raise `bin_size`).

- **API:** `dnakit.variants.plot_density(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `bin_size` (default 1000000), `chroms` (list of contigs to include), `chrom_lengths` (`{contig: length}`), `svg_output` (optional path to write the SVG), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `bin_size`, `bins` (each with `chrom`, `start`, `end`, `count`, `density_per_mb`), `svg` (the SVG text), `output_path` and `max_count`.

- **Example code:**

```python
import json

from dnakit.variants import plot_density

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t1500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2200\t.\tA\tG\t50\tPASS\t.",
    "chr2\t300\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = plot_density(vcf, options={"bin_size": 1000})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "plot_density",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 8,
    "text_length": 213,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b43dee4028171f8b561a43f9b1b12264f39c2cfb625d601dd8e188346be3938"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_size": 1000,
    "bins": [
      {"chrom": "chr1", "start": 1, "end": 1000, "count": 2, "density_per_mb": 2000.0},
      {"chrom": "chr1", "start": 1001, "end": 2000, "count": 1, "density_per_mb": 1000.0},
      "..."
    ],
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"320\" height=\"220\" viewBox=\"0 0 320 220\"><te...",
    "output_path": null,
    "max_count": 2
  }
}
```

## 57) plot_density_mt

- **Purpose:** Plot variant density of a variant-by-sample (matrix-table style) call set as a binned SVG bar chart.

- **Method:** `plot_density_mt` is the variant-by-sample (matrix-table style) name of `plot_density`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Variants are counted in fixed-size bins of `bin_size` bases per chromosome (bin index = (POS - 1) // bin_size); chromosomes are ordered naturally and, when `chrom_lengths` is given, empty trailing bins up to the chromosome length are included. Density per Mb = count * 1,000,000 / bin_size. The SVG is a plain bar chart (bar height proportional to the count, scaled to the maximum bin) with chromosome labels, built as text; no plotting library is needed. More than 20,000 bins is an error (raise `bin_size`).

- **API:** `dnakit.variants.plot_density_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `bin_size` (default 1000000), `chroms` (list of contigs to include), `chrom_lengths` (`{contig: length}`), `svg_output` (optional path to write the SVG), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `bin_size`, `bins` (each with `chrom`, `start`, `end`, `count`, `density_per_mb`), `svg` (the SVG text), `output_path` and `max_count`.

- **Example code:**

```python
import json

from dnakit.variants import plot_density_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t1500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2200\t.\tA\tG\t50\tPASS\t.",
    "chr2\t300\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = plot_density_mt(vcf, options={"bin_size": 1000})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "plot_density_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 8,
    "text_length": 213,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b43dee4028171f8b561a43f9b1b12264f39c2cfb625d601dd8e188346be3938"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_size": 1000,
    "bins": [
      {"chrom": "chr1", "start": 1, "end": 1000, "count": 2, "density_per_mb": 2000.0},
      {"chrom": "chr1", "start": 1001, "end": 2000, "count": 1, "density_per_mb": 1000.0},
      "..."
    ],
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"320\" height=\"220\" viewBox=\"0 0 320 220\"><te...",
    "output_path": null,
    "max_count": 2
  }
}
```

## 58) prioritize_candidate_variants

- **Purpose:** Rank candidate variants with a transparent additive score computed from annotations that are already on the records.

- **Method:** Annotations are read from record keys or INFO (`consequence`/`most_severe_consequence`/`effect` or INFO `MC`, `clinvar`/`CLNSIG`, `gnomad_af`, `cadd`, `revel`, `gene`/`GENEINFO`); nothing is looked up. The score is the sum of: consequence severity from a fixed Sequence Ontology table (0-1: loss-of-function terms 1.0, `stop_lost` 0.8, in-frame 0.6, missense 0.5, splice region 0.3, synonymous 0.05, ...; the most severe term counts); ClinVar (pathogenic 1.0, likely pathogenic 0.8, uncertain 0.1, benign -1.0); rarity (0.5 when the population AF is absent or <= 1e-4, 0.3 up to 1e-3, 0.1 up to 1e-2, otherwise 0); predictors, 0.5 * max(CADD / 40, REVEL); and +0.5 for genes in `options.genes`. Variants with population AF above `max_af`, and ClinVar benign/likely benign variants when `exclude_benign` is true, are listed separately as excluded. This is a heuristic ranking with arbitrary weights, not a validated classifier; a variant without any population-frequency annotation keeps the full rarity score (0.5) but is flagged `frequency_unknown: true` and counted in one summary warning in `warnings`, so it can be told apart from a variant that is known to be rare. Every score component is returned so the ranking can be audited.

- **API:** `dnakit.variants.prioritize_candidate_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `max_af` (0.05; variants above it are excluded), `exclude_benign` (True), `genes` (panel list, +0.5), `top_n` (keep only the best n ranked variants).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `ranked` (each with `chrom`, `pos`, `ref`, `alt`, `id`, `score`, `components`, `excluded_reason` = `null`, `gene`, `rank`), `excluded` (same fields with the `excluded_reason`), `n_ranked` (before `top_n`), `n_excluded` and `n_frequency_unknown` (number of variants, ranked or excluded, without a population-frequency annotation; the envelope `warnings` then holds one message with that count). Every ranked and excluded row has the boolean `frequency_unknown`.

- **Example code:**

```python
import json

from dnakit.variants import prioritize_candidate_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=BRCA1;consequence=stop_gained;clinvar=Pathogenic;gnomad_af=0.00001;cadd=35",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=TTN;consequence=missense_variant;gnomad_af=0.002;revel=0.6",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=OBSCN;consequence=synonymous_variant;gnomad_af=0.2",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=APOE;consequence=missense_variant;clinvar=Benign",
]) + "\n"

result = prioritize_candidate_variants(vcf, options={"genes": ["BRCA1"]})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "prioritize_candidate_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 406,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1462d42c295ee82b6444284c676a34284ef9d2be57066750c3aee1651621c7f2"
  },
  "warnings": ["1 of 4 variants have no population-frequency annotation and received the full rarity score..."],
  "status": "ok",
  "result": {
    "method": "additive heuristic score: consequence severity (Sequence Ontology table, 0-1) + ClinVar (p...",
    "parameters": {"max_af": 0.05, "exclude_benign": true, "genes": ["BRCA1"]},
    "ranked": [
      {
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": ["G"],
        "id": null,
        "score": 3.4375,
        "components": {"consequence": 1.0, "clinvar": 1.0, "rarity": 0.5, "predictors": 0.4375, "panel_gene": 0.5},
        "excluded_reason": null,
        "gene": "BRCA1",
        "frequency_unknown": false,
        "rank": 1
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": ["T"],
        "id": null,
        "score": 0.9,
        "components": {"consequence": 0.5, "rarity": 0.1, "predictors": 0.3},
        "excluded_reason": null,
        "gene": "TTN",
        "frequency_unknown": false,
        "rank": 2
      }
    ],
    "excluded": [
      {
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": ["A"],
        "id": null,
        "score": 0.05,
        "components": {"consequence": 0.05, "rarity": 0.0},
        "excluded_reason": "population AF 0.2 > 0.05",
        "gene": "OBSCN",
        "frequency_unknown": false
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "ref": "T",
        "alt": ["C"],
        "id": null,
        "score": 0.0,
        "components": {"consequence": 0.5, "clinvar": -1.0, "rarity": 0.5},
        "excluded_reason": "ClinVar benign / likely benign",
        "gene": "APOE",
        "frequency_unknown": true
      }
    ],
    "n_ranked": 2,
    "n_excluded": 2,
    "n_frequency_unknown": 1
  }
}
```

## 59) split_multiallelic

- **Purpose:** Document a reusable split multiallelic capability for Variant analysis workflows. Candidate CAP-0817 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split multiallelic operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.split_multiallelic(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a sequence or mapping containing VCF-like records; each record requires `chrom`, positive 1-based `pos`, `ref`, and `alt`.

- **Output:** The local `split_multiallelic` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import split_multiallelic

result = split_multiallelic(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": ["G", "T"], "genotypes": ["0/1", "1/2"]}]},
    options={"remap_genotypes": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0817",
  "operation": "split_multiallelic",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'split_multiallelic', 'result': {'records': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/0'], 'original_alt_index': 1}, {'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'T', 'genotypes': ['0/0', '0/1'], 'original_alt_index': 2}]...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 60) split_multiallelic_mt

- **Purpose:** Split multi-allelic sites of a variant-by-sample (matrix-table style) call set into one record per ALT allele.

- **Method:** `split_multiallelic_mt` is the variant-by-sample (matrix-table style) name of `split_multiallelic_variants`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Each ALT of a multi-allelic record becomes its own record with the same CHROM/POS/REF. GT is remapped (the selected ALT becomes 1 and every other allele, including the REF, becomes 0, so a `1/2` genotype becomes `1/0` in the first record and `0/1` in the second); `FORMAT/AD` and `INFO` `AC`/`AF` are subset to the selected allele, and `FORMAT/PL` is dropped because it cannot be subset. The alleles are not trimmed or left-aligned (use `normalize` with `split_multiallelic: True` for that). Biallelic records are copied unchanged.

- **API:** `dnakit.variants.split_multiallelic_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records` (after splitting), `samples`, `vcf`, `output_path`, `records`, `input_records` and `split_records` (the number of records added).

- **Example code:**

```python
import json

from dnakit.variants import split_multiallelic_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\tAC=2,1\tGT:AD\t1/2:5,3,2\t0/1:6,4,0",
]) + "\n"

result = split_multiallelic_mt(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "split_multiallelic_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 131,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2cdf597390051f5fe42bb240739a02d0f21405e197041be6028c4f9d58996f97"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t400\t.\tG\tA\t30...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "2"},
        "genotypes": {"S1": {"GT": "1/0", "AD": "5,3"}, "S2": {"GT": "0/1", "AD": "6,4"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "1"},
        "genotypes": {"S1": {"GT": "0/1", "AD": "5,2"}, "S2": {"GT": "0/0", "AD": "6,0"}}
      }
    ],
    "input_records": 1,
    "split_records": 1
  }
}
```

## 61) split_multiallelic_variants

- **Purpose:** Split multi-allelic variant records into one record per ALT allele.

- **Method:** Each ALT of a multi-allelic record becomes its own record with the same CHROM/POS/REF. GT is remapped (the selected ALT becomes 1 and every other allele, including the REF, becomes 0, so a `1/2` genotype becomes `1/0` in the first record and `0/1` in the second); `FORMAT/AD` and `INFO` `AC`/`AF` are subset to the selected allele, and `FORMAT/PL` is dropped because it cannot be subset. The alleles are not trimmed or left-aligned (use `normalize` with `split_multiallelic: True` for that). Biallelic records are copied unchanged.

- **API:** `dnakit.variants.split_multiallelic_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `output` (VCF path), `include_records` (default True), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records` (after splitting), `samples`, `vcf`, `output_path`, `records`, `input_records` and `split_records` (the number of records added).

- **Example code:**

```python
import json

from dnakit.variants import split_multiallelic_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\tAC=2,1\tGT:AD\t1/2:5,3,2\t0/1:6,4,0",
]) + "\n"

result = split_multiallelic_variants(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "split_multiallelic_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 131,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2cdf597390051f5fe42bb240739a02d0f21405e197041be6028c4f9d58996f97"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t400\t.\tG\tA\t30...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "2"},
        "genotypes": {"S1": {"GT": "1/0", "AD": "5,3"}, "S2": {"GT": "0/1", "AD": "6,4"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "1"},
        "genotypes": {"S1": {"GT": "0/1", "AD": "5,2"}, "S2": {"GT": "0/0", "AD": "6,0"}}
      }
    ],
    "input_records": 1,
    "split_records": 1
  }
}
```

## 62) summarize

- **Purpose:** Summarise a call set: counts by variant type and chromosome, Ti/Tv, FILTER counts, multi-allelic sites and per-sample QC.

- **Method:** Variants are classified per ALT allele by the allele lengths after trimming shared bases: `SNV`, `MNV`, `indel` (subtype insertion, deletion or complex) and `SV` (symbolic alleles or a length change of at least 50 bases). Ti/Tv counts SNVs only (transitions A<->G and C<->T). The per-sample QC reports call rate, hom-ref/het/hom-alt counts, het/hom-alt ratio, SNV and indel counts, per-sample Ti/Tv and mean `DP` and `GQ`. These are plain counts of what is in the file; no filtering or reference is used.

- **API:** `dnakit.variants.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `n_alt_alleles`, `n_multiallelic_sites`, `n_samples`, `samples`, `type_counts`, `subtype_counts`, `by_chromosome`, `filter_counts`, `transitions`, `transversions`, `ti_tv`, `mean_qual`, `indel_length_counts`, `mean_sample_call_rate` and `sample_qc` (per sample: `call_rate`, `n_hom_ref`, `n_het`, `n_hom_alt`, `het_hom_alt_ratio`, `n_snv`, `n_indel`, `ti_tv`, `mean_depth`, `mean_genotype_quality`, ...). Ratios are `null` when the denominator is 0.

- **Example code:**

```python
import json

from dnakit.variants import summarize

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = summarize(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "summarize",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 4,
    "n_alt_alleles": 5,
    "n_multiallelic_sites": 1,
    "n_samples": 3,
    "samples": ["S1", "S2", "S3"],
    "type_counts": {"SNV": 4, "indel": 1},
    "subtype_counts": {"indel:deletion": 1},
    "by_chromosome": {"chr1": 4},
    "filter_counts": {"PASS": 3, "q10": 1},
    "transitions": 3,
    "transversions": 1,
    "ti_tv": 3.0,
    "mean_qual": 40.0,
    "indel_length_counts": {"-1": 1},
    "mean_sample_call_rate": 0.916667,
    "sample_qc": [
      {
        "sample": "S1",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 2,
        "n_het": 2,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 1,
        "ti_tv": 2.0,
        "mean_depth": 10.0,
        "mean_genotype_quality": 50.0
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 0,
        "n_het": 4,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 1,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 40.0
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 2,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 30.0
      }
    ]
  }
}
```

## 63) summarize_mt

- **Purpose:** Summarise a variant-by-sample (matrix-table style) call set: type counts, Ti/Tv, FILTER counts and per-sample QC.

- **Method:** `summarize_mt` is the variant-by-sample (matrix-table style) name of `summarize`: it runs the same implementation on VCF text/path/records; no Hail MatrixTable object is read. Variants are classified per ALT allele by the allele lengths after trimming shared bases: `SNV`, `MNV`, `indel` (subtype insertion, deletion or complex) and `SV` (symbolic alleles or a length change of at least 50 bases). Ti/Tv counts SNVs only (transitions A<->G and C<->T). The per-sample QC reports call rate, hom-ref/het/hom-alt counts, het/hom-alt ratio, SNV and indel counts, per-sample Ti/Tv and mean `DP` and `GQ`. These are plain counts of what is in the file; no filtering or reference is used.

- **API:** `dnakit.variants.summarize_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `n_records`, `n_alt_alleles`, `n_multiallelic_sites`, `n_samples`, `samples`, `type_counts`, `subtype_counts`, `by_chromosome`, `filter_counts`, `transitions`, `transversions`, `ti_tv`, `mean_qual`, `indel_length_counts`, `mean_sample_call_rate` and `sample_qc` (per sample: `call_rate`, `n_hom_ref`, `n_het`, `n_hom_alt`, `het_hom_alt_ratio`, `n_snv`, `n_indel`, `ti_tv`, `mean_depth`, `mean_genotype_quality`, ...). Ratios are `null` when the denominator is 0.

- **Example code:**

```python
import json

from dnakit.variants import summarize_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = summarize_mt(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "variants",
  "operation": "summarize_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 4,
    "n_alt_alleles": 5,
    "n_multiallelic_sites": 1,
    "n_samples": 3,
    "samples": ["S1", "S2", "S3"],
    "type_counts": {"SNV": 4, "indel": 1},
    "subtype_counts": {"indel:deletion": 1},
    "by_chromosome": {"chr1": 4},
    "filter_counts": {"PASS": 3, "q10": 1},
    "transitions": 3,
    "transversions": 1,
    "ti_tv": 3.0,
    "mean_qual": 40.0,
    "indel_length_counts": {"-1": 1},
    "mean_sample_call_rate": 0.916667,
    "sample_qc": [
      {
        "sample": "S1",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 2,
        "n_het": 2,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 1,
        "ti_tv": 2.0,
        "mean_depth": 10.0,
        "mean_genotype_quality": 50.0
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 0,
        "n_het": 4,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 1,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 40.0
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 2,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 30.0
      }
    ]
  }
}
```
