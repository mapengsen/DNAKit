# Clinical genomics

This category contains **18** documented functions. 12 of them compute their result from the evidence passed in (trio genotypes, plasma read counts, mutation matrices, ACMG/AMP criteria, VCF or ClinVar-style records); 6 always return `status: "not_implemented"` because they need reference databases, a trained model, an external caller or program, or have no defined meaning (`annotate`, `calculate`, `call`, `compare`, `predict`, `type_hla`). Every computing function works only on the evidence the caller supplies and none is a validated diagnostic. In particular, `classify` applies the ACMG/AMP 2015 combining rules to criteria that the caller has already decided are met (the criteria are not evaluated from data, and contradictory pathogenic and benign evidence gives Uncertain significance); `detect_circulating_tumor_dna` and `detect_minimal_residual_disease` use a default background error rate of 0.001 that is only a placeholder and has to be measured for the assay; `generate_clinical_variant_report` and `report` summarise supplied classifications and do not replace review by a clinical geneticist. Malformed input or a missing required option raises `InputFormatError` or `ConfigurationError` (with a `code` and a `hint`) instead of returning an envelope. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) analyze_trio_variants

- **Purpose:** Classify the variants of a child-father-mother trio by inheritance pattern: de novo, Mendelian error, inherited from one or both parents, homozygous with carrier parents, and compound-heterozygous candidates.

- **Method:** Multi-allelic sites are split first; for every ALT the child's diploid genotype is compared with the parents'. If the child genotype cannot be formed from one allele of each parent the site is `de_novo` (child carries the ALT and both parents are homozygous without it) or otherwise `mendelian_error` (including a child homozygous reference with a homozygous-ALT parent). Consistent sites with the ALT are `homozygous_both_parents_carriers`, `inherited_biparental`, `inherited_paternal` or `inherited_maternal`; sites where a parent genotype is missing and the child carries the ALT are `undetermined_missing_parent`. Compound-heterozygous candidates are pairs, within one gene, of a child-heterozygous variant carried only by the father and one carried only by the mother (gene symbols come from a `gene` key or INFO `GENE`/`SYMBOL`/`GENEINFO`; variants without a gene are not considered). Sites failing `min_depth`/`min_genotype_quality` in any of the three samples are excluded when `require_quality_pass` is true. Autosomal diploid logic only: X-linked and mitochondrial inheritance are not modelled, and phase is not used, so compound heterozygosity is a candidate list, not a confirmed in-trans call.

- **API:** `dnakit.clinical.analyze_trio_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`); the three trio samples are required. `options`: `child`, `father`, `mother` (sample names; when absent and the input has exactly three samples the column order child, father, mother is assumed, with a warning), `min_depth` (0), `min_genotype_quality` (0), `require_quality_pass` (True).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `category_counts`, `variants` (each with `site`, `chrom`, `pos`, `ref`, `alt`, `gene`, `category`, `quality_pass`, `child_gt`, `father_gt`, `mother_gt`), `compound_heterozygous_candidates` (`gene`, `paternal_variant`, `maternal_variant`) and `n_filtered_low_quality`. Reference-only child sites consistent with the parents are not listed. Wrong or missing sample names raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.clinical import analyze_trio_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=GENE1\tGT:DP:GQ\t0/1:30:60\t0/0:30:60\t0/0:30:60",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=GENE2\tGT:DP:GQ\t0/1:30:60\t0/1:30:60\t0/0:30:60",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=GENE2\tGT:DP:GQ\t0/1:30:60\t0/0:30:60\t0/1:30:60",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=GENE3\tGT:DP:GQ\t1/1:30:60\t0/1:30:60\t0/1:30:60",
    "chr1\t500\t.\tA\tC\t50\tPASS\tGENE=GENE4\tGT:DP:GQ\t0/0:30:60\t1/1:30:60\t0/0:30:60",
]) + "\n"

result = analyze_trio_variants(vcf, options={"child": "child", "father": "father", "mother": "mother"})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "analyze_trio_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 7,
    "text_length": 452,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eb0bb0a2db5f091abacfdcce919bbe3dd13373aaf43148678619612ad63bf6df"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-variant inheritance category from child/parent genotypes (after splitting multi-alleli...",
    "category_counts": {
      "de_novo": 1,
      "homozygous_both_parents_carriers": 1,
      "inherited_maternal": 1,
      "inherited_paternal": 1,
      "mendelian_error": 1
    },
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": "GENE1",
        "category": "de_novo",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/0",
        "mother_gt": "0/0"
      },
      {
        "site": "chr1:200:C:T",
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": "T",
        "gene": "GENE2",
        "category": "inherited_paternal",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/1",
        "mother_gt": "0/0"
      },
      "..."
    ],
    "compound_heterozygous_candidates": [{"gene": "GENE2", "paternal_variant": "chr1:200:C:T", "maternal_variant": "chr1:300:G:A"}],
    "n_filtered_low_quality": 0
  }
}
```

## 2) annotate

- **Purpose:** Annotate clinical variants with ClinVar, gnomAD and gene-disease information. Not implemented: it only reports what is missing (ClinVar, gnomAD and gene-disease reference databases).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Use the look-up functions of `dnakit.databases` / `dnakit.annotation` for individual variants; `read` extracts ClinVar fields that are already in a VCF.

- **API:** `dnakit.clinical.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("ClinVar / gnomAD / gene-disease databases (local files or network access)") and `message` ("clinical annotation needs reference databases that are not bundled."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import annotate

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
  "module": "clinical",
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
  "required_dependency": "ClinVar / gnomAD / gene-disease databases (local files or network access)",
  "message": "clinical annotation needs reference databases that are not bundled."
}
```

## 3) build_cancer_evolution_tree

- **Purpose:** Build a rooted perfect-phylogeny evolution tree of somatic mutations from a binary mutation-by-sample (or cell) matrix.

- **Method:** Each mutation is the set of samples/cells that carry it. Mutations with identical sets share a node; a node's parent is the node with the smallest strictly larger sample set, and the root is the germline. The tree exists only if the sample sets are pairwise nested or disjoint; if two mutations' sets partially overlap (a four-gamete violation) the result reports `perfect_phylogeny: false` with the violating pairs and builds no tree. It therefore assumes the infinite-sites model (no back mutation or recurrent mutation) and error-free calls; there is no likelihood model, error tolerance or clone-frequency estimation. Mutations absent from all samples are ignored with a warning.

- **API:** `dnakit.clinical.build_cancer_evolution_tree(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'mutations': {'TP53': [1, 1, 0], ...}, 'samples': [...]}` (0/1 presence per sample; or `mutation -> {sample: 0/1}`; at most 5000 mutations, all vectors of equal length; `samples` defaults to `sample_1`, ...). `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `perfect_phylogeny`, `violations` (`mutation_1`, `mutation_2`, `reason`), `n_mutations` and, when a tree exists, `nodes` (`id`, `label`, `mutations`, `samples`, `parent`, `children`, `depth`; the first is `root` = `germline`), `newick` and `sample_to_node` (the deepest node each sample belongs to). When there are violations, `tree` is `null` and `nodes` is absent.

- **Example code:**

```python
import json

from dnakit.clinical import build_cancer_evolution_tree

result = build_cancer_evolution_tree({
    "samples": ["S1", "S2", "S3", "S4"],
    "mutations": {
        "TP53": [1, 1, 1, 1],
        "KRAS": [1, 1, 0, 0],
        "PIK3CA": [1, 0, 0, 0],
        "APC": [0, 0, 1, 1],
    },
})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "build_cancer_evolution_tree",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 141,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d564184bd9487dc45696fc6e0f66a5f2e7510cc777268e6351a60433e7f8d170"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "perfect-phylogeny (infinite-sites) tree from a binary mutation x sample matrix: mutations...",
    "perfect_phylogeny": true,
    "violations": [],
    "n_mutations": 4,
    "nodes": [
      {
        "id": "root",
        "label": "germline",
        "mutations": [],
        "samples": ["S1", "S2", "..."],
        "parent": null,
        "children": ["n1"],
        "depth": 0
      },
      {
        "id": "n1",
        "label": "TP53",
        "mutations": ["TP53"],
        "samples": ["S1", "S2", "..."],
        "parent": "root",
        "children": ["n2", "n3"],
        "depth": 1
      },
      "..."
    ],
    "newick": "(((PIK3CA)KRAS,APC)TP53)germline;",
    "sample_to_node": {"S1": "n4", "S2": "n2", "S3": "n3", "S4": "n3"}
  }
}
```

## 4) calculate

- **Purpose:** Calculate an unspecified clinical quantity. Not implemented: it only reports what is missing (a specific clinical calculation; the generic name `calculate` has no defined meaning).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Call a specific function instead (`classify`, `qc`, `detect_de_novo_variants`, ...).

- **API:** `dnakit.clinical.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("a specific clinical calculation (this generic name has no defined meaning)") and `message` ("'calculate' is too unspecific to have a defined meaning; call a specific clinical function."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import calculate

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = calculate(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "clinical",
  "operation": "calculate",
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
  "required_dependency": "a specific clinical calculation (this generic name has no defined meaning)",
  "message": "'calculate' is too unspecific to have a defined meaning; call a specific clinical function."
}
```

## 5) call

- **Purpose:** Call clinical variants from sequencing reads. Not implemented: it only reports what is missing (an external variant caller and aligned reads).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. Use `dnakit.variants.call_freebayes`, `call_gatk_haplotypecaller` or `call_deepvariant`, which run the external caller when it is installed.

- **API:** `dnakit.clinical.call(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("variant caller binary and aligned reads") and `message` ("variant calling needs an external caller and aligned reads."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import call

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
  "module": "clinical",
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
  "required_dependency": "variant caller binary and aligned reads",
  "message": "variant calling needs an external caller and aligned reads."
}
```

## 6) classify

- **Purpose:** Combine ACMG/AMP 2015 evidence criteria into a five-tier variant classification (Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign).

- **Method:** The criteria that apply to the variant are supplied by the caller and are not evaluated from data. Each code (PVS1; PS1-4; PM1-6; PP1-5; BA1; BS1-4; BP1-7) has its default strength (very strong, strong, moderate, supporting, stand-alone), which a suffix changes (for example `PM2_Supporting`, `PVS1_Strong`); a code given twice with different strengths, an unknown code or a strength not defined for that direction is an error. The combination applies the rules of Table 5 in Richards et al. 2015: Pathogenic (PVS1 + >=1 PS; PVS1 + >=2 PM; PVS1 + PM + PP; PVS1 + >=2 PP; >=2 PS; PS + >=3 PM; PS + 2 PM + >=2 PP; PS + PM + >=4 PP), Likely pathogenic (PVS1 + PM; PS + 1-2 PM; PS + >=2 PP; >=3 PM; 2 PM + >=2 PP; PM + >=4 PP), Benign (BA1; >=2 BS) and Likely benign (BS + BP; >=2 BP). The pathogenic side and the benign side are evaluated separately. If only one side reaches a classification that is the result; if neither side does, the result is Uncertain significance; if both sides reach a classification (contradictory evidence, for example PVS1 + PS1 together with two BS), the result is Uncertain significance and `rule` names both conflicting rules (the 2015 guideline's default for conflicting evidence). Only these discrete count rules are implemented: the later Bayesian points framework (Tavtigian et al. 2018/2020), gene- or disease-specific ClinGen specifications and automatic criterion evaluation are not, and the result is not a clinical interpretation.

- **API:** `dnakit.clinical.classify(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the criteria met, as a list of codes (`['PVS1', 'PM2_Supporting']`, entries may also be `{'code', 'strength'}`), `{'criteria': [...], 'id': ...}`, or `{'variants': [{'id': ..., 'criteria': [...]}, ...]}` for several variants. `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `framework` and, for one variant, `classification`, `rule` (the rule that applied, or the description of the conflict, or "criteria do not meet any pathogenic or benign combination"), `evidence_counts` (counts per direction and strength), `criteria` (normalised `code`, `strength`, `direction`) and `id` when given; for several variants `results` holds one such entry per variant. An empty criteria list gives Uncertain significance with a warning. Invalid codes raise `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.clinical import classify

result = classify({
    "variants": [
        {"id": "v1", "criteria": ["PVS1", "PS1"]},
        {"id": "v2", "criteria": ["PM2_Supporting", "PP3"]},
        {"id": "v3", "criteria": ["PVS1", "PS1", "BS1", "BS2"]},  # pathogenic and benign rules both met
    ]
})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "classify",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 167,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2e728bced9ea44abed7aabb390830829cd366055831cd1a3492f622e97162582"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "framework": "ACMG/AMP 2015 (Richards et al., Genet Med 17:405-424) Table 5 combining rules applied to c...",
    "results": [
      {
        "id": "v1",
        "classification": "Pathogenic",
        "rule": "Pathogenic (i)(a): 1 Very strong + >=1 Strong",
        "evidence_counts": {
          "pathogenic_very_strong": 1,
          "pathogenic_strong": 1,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 0,
          "benign_stand_alone": 0,
          "benign_strong": 0,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PVS1", "strength": "very_strong", "direction": "pathogenic"},
          {"code": "PS1", "strength": "strong", "direction": "pathogenic"}
        ]
      },
      {
        "id": "v2",
        "classification": "Uncertain significance",
        "rule": "criteria do not meet any pathogenic or benign combination",
        "evidence_counts": {
          "pathogenic_very_strong": 0,
          "pathogenic_strong": 0,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 2,
          "benign_stand_alone": 0,
          "benign_strong": 0,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PM2", "strength": "supporting", "direction": "pathogenic"},
          {"code": "PP3", "strength": "supporting", "direction": "pathogenic"}
        ]
      },
      {
        "id": "v3",
        "classification": "Uncertain significance",
        "rule": "contradictory evidence: Pathogenic (i)(a): 1 Very strong + >=1 Strong vs Benign (ii): >=2 Strong",
        "evidence_counts": {
          "pathogenic_very_strong": 1,
          "pathogenic_strong": 1,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 0,
          "benign_stand_alone": 0,
          "benign_strong": 2,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PVS1", "strength": "very_strong", "direction": "pathogenic"},
          {"code": "PS1", "strength": "strong", "direction": "pathogenic"},
          "..."
        ]
      }
    ]
  }
}
```

## 7) compare

- **Purpose:** Compare clinical data sets. Not implemented: it only reports what is missing (a specific comparison definition; the generic name `compare` has no defined meaning here).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. For variant call sets use `dnakit.variants.compare_vcfs`.

- **API:** `dnakit.clinical.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("a specific comparison definition (this generic name has no defined meaning)") and `message` ("'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import compare

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = compare(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "clinical",
  "operation": "compare",
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
  "required_dependency": "a specific comparison definition (this generic name has no defined meaning)",
  "message": "'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets."
}
```

## 8) detect_circulating_tumor_dna

- **Purpose:** Test for circulating tumour DNA in plasma at tracked tumour variants (tumour-informed detection by a pooled binomial test).

- **Method:** The plasma ALT reads and depths at the tracked tumour variants are pooled; the pooled ALT count is tested against the background error rate with a one-sided binomial test, P(X >= alt | depth, error). Detected = pooled p <= `alpha` AND at least `min_positive_sites` sites individually positive (ALT >= `min_alt_reads` and p <= `site_alpha`). The estimated tumour fraction is `vaf_to_fraction_factor` * max(0, pooled VAF - background error), capped at 1 (a factor of 2 corresponds to clonal heterozygous, copy-neutral SNVs). The default `background_error_rate` of 0.001 is only a placeholder: it must be measured for the assay (for example from matched healthy plasma or from the non-reference bases at the tracked positions) and passed in, because the detection call and the p-values depend on it directly. Reads are assumed independent with one error rate (a per-site `error_rate` can be supplied), and there is no UMI/duplex error suppression, strand or context model, so this is a statistical sketch and not a validated clinical assay.

- **API:** `dnakit.clinical.detect_circulating_tumor_dna(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: tracked sites `[{'alt_count', 'depth'[, 'id', 'error_rate']}, ...]` or `{'sites': [...]}` (plasma read counts, integers with 0 <= alt_count <= depth), or a VCF/records with `FORMAT/AD` for the plasma sample (`options.sample`, default the first sample; multi-allelic sites are split). `options`: `background_error_rate` (0.001; a placeholder to be measured), `alpha` (0.01), `site_alpha` (0.01), `min_alt_reads` (2), `min_positive_sites` (1), `vaf_to_fraction_factor` (2.0).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `detected`, `n_sites`, `n_positive_sites`, `pooled_alt_reads`, `pooled_depth`, `pooled_vaf`, `pooled_p_value`, `estimated_tumor_fraction` and `sites` (per site: `id`, `alt_count`, `depth`, `error_rate`, `vaf`, `p_value`, `positive`).

- **Example code:**

```python
import json

from dnakit.clinical import detect_circulating_tumor_dna

# plasma reads at 4 tumour-informed sites: alt_count out of depth
sites = [
    {"id": "TP53:R175H", "alt_count": 12, "depth": 5000},
    {"id": "KRAS:G12D", "alt_count": 9, "depth": 4800},
    {"id": "APC:R1450*", "alt_count": 1, "depth": 5200},
    {"id": "PIK3CA:E545K", "alt_count": 10, "depth": 4900},
]

result = detect_circulating_tumor_dna({"sites": sites}, options={"background_error_rate": 0.0005})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "detect_circulating_tumor_dna",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 226,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d40e8a937f2b6f0485bbf2cd42cad1a6a7f68bf85efc96365e6faad9baa1ebba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "ctDNA detection: tumour-informed detection. The plasma ALT reads and depths at the tracked...",
    "parameters": {"background_error_rate": 0.0005, "alpha": 0.01, "site_alpha": 0.01, "min_alt_reads": 2, "min_positive_sites": 1},
    "detected": true,
    "n_sites": 4,
    "n_positive_sites": 3,
    "pooled_alt_reads": 32,
    "pooled_depth": 19900,
    "pooled_vaf": 0.00160804,
    "pooled_p_value": 2.17455e-08,
    "estimated_tumor_fraction": 0.00221608,
    "sites": [
      {
        "id": "TP53:R175H",
        "alt_count": 12,
        "depth": 5000,
        "error_rate": 0.0005,
        "vaf": 0.0024,
        "p_value": 1.24946e-05,
        "positive": true
      },
      {
        "id": "KRAS:G12D",
        "alt_count": 9,
        "depth": 4800,
        "error_rate": 0.0005,
        "vaf": 0.001875,
        "p_value": 0.000858515,
        "positive": true
      },
      "..."
    ]
  }
}
```

## 9) detect_de_novo_variants

- **Purpose:** Detect de novo variant candidates in a child-father-mother trio from genotypes, depths and allelic depths.

- **Method:** A candidate is a site (after splitting multi-allelic sites) where the child carries an ALT allele and both parents are called homozygous without it. Candidates are then filtered: child depth >= `min_child_depth`, parent depth >= `min_parent_depth`, GQ >= `min_genotype_quality` in all three samples, and, when `FORMAT/AD` exists, child ALT fraction >= `min_child_vaf` and parental ALT reads <= `max_parent_alt_reads`. Each rejected candidate is counted once in `rejected_by`, under the FIRST check it fails in a fixed order (the genotype checks that define a candidate come first, then: child, father and mother depth; child, father and mother GQ; child ALT fraction; father and mother ALT reads), so the reason is deterministic (for example `child_low_depth`, `child_low_vaf`, `father_alt_reads`). The fixed thresholds are generic defaults and the function has no knowledge of regions of low mappability, parental mosaicism or sequencing-error profiles, so the output is a candidate list that needs orthogonal validation, not a validated de novo call.

- **API:** `dnakit.clinical.detect_de_novo_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`); the three trio samples are required. `options`: `child`, `father`, `mother` (sample names; with exactly three samples the order child, father, mother is assumed, with a warning), `min_child_depth` (10), `min_parent_depth` (10), `min_genotype_quality` (20), `min_child_vaf` (0.25, applied when AD exists), `max_parent_alt_reads` (0, applied when AD exists).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `n_candidates_before_filters`, `rejected_by`, `n_de_novo` and `variants` (each with `site`, `chrom`, `pos`, `ref`, `alt`, `gene`, `category` = `de_novo`, `quality_pass`, `child_gt`, `father_gt`, `mother_gt`, `child_ad`).

- **Example code:**

```python
import json

from dnakit.clinical import detect_de_novo_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:35:80:20,15\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:8:80:5,3\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:40:80:36,4\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t400\t.\tT\tC\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:35:80:20,15\t0/0:30:60:28,2\t0/0:32:60:32,0",
]) + "\n"

result = detect_de_novo_variants(vcf, options={"child": "child", "father": "father", "mother": "mother"})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "detect_de_novo_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 415,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c5fffaae813fd43f7d1b4c5850ada93d33e140429b1303180591f6048a810982"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "de novo candidates = child carries the ALT allele while both parents are called homozygous...",
    "parameters": {
      "min_child_depth": 10.0,
      "min_parent_depth": 10.0,
      "min_genotype_quality": 20.0,
      "min_child_vaf": 0.25,
      "max_parent_alt_reads": 0
    },
    "n_candidates_before_filters": 4,
    "rejected_by": {"child_low_depth": 1, "child_low_vaf": 1, "father_alt_reads": 1},
    "n_de_novo": 1,
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": null,
        "category": "de_novo",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/0",
        "mother_gt": "0/0",
        "child_ad": [20, 15]
      }
    ]
  }
}
```

## 10) detect_minimal_residual_disease

- **Purpose:** Test for minimal residual disease (MRD) from plasma read counts at tracked tumour variants with a stricter tumour-informed pooled binomial test.

- **Method:** Same statistics as `detect_circulating_tumor_dna` with stricter defaults (`alpha` 0.001 and at least two individually positive sites). The plasma ALT reads and depths at the tracked tumour variants are pooled; the pooled ALT count is tested against the background error rate with a one-sided binomial test, P(X >= alt | depth, error). Detected = pooled p <= `alpha` AND at least `min_positive_sites` sites individually positive (ALT >= `min_alt_reads` and p <= `site_alpha`). The estimated tumour fraction is `vaf_to_fraction_factor` * max(0, pooled VAF - background error), capped at 1 (a factor of 2 corresponds to clonal heterozygous, copy-neutral SNVs). The default `background_error_rate` of 0.001 is only a placeholder: it must be measured for the assay (for example from matched healthy plasma or from the non-reference bases at the tracked positions) and passed in, because the detection call and the p-values depend on it directly. Reads are assumed independent with one error rate (a per-site `error_rate` can be supplied), and there is no UMI/duplex error suppression, strand or context model, so this is a statistical sketch and not a validated clinical assay.

- **API:** `dnakit.clinical.detect_minimal_residual_disease(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: tracked sites `[{'alt_count', 'depth'[, 'id', 'error_rate']}, ...]` or `{'sites': [...]}` (plasma read counts, integers with 0 <= alt_count <= depth), or a VCF/records with `FORMAT/AD` for the plasma sample (`options.sample`, default the first sample; multi-allelic sites are split). `options`: `background_error_rate` (0.001; a placeholder to be measured), `alpha` (0.001), `site_alpha` (0.01), `min_alt_reads` (2), `min_positive_sites` (2), `vaf_to_fraction_factor` (2.0).

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `detected`, `n_sites`, `n_positive_sites`, `pooled_alt_reads`, `pooled_depth`, `pooled_vaf`, `pooled_p_value`, `estimated_tumor_fraction` and `sites` (per site: `id`, `alt_count`, `depth`, `error_rate`, `vaf`, `p_value`, `positive`).

- **Example code:**

```python
import json

from dnakit.clinical import detect_minimal_residual_disease

# plasma reads at 4 tumour-informed sites: alt_count out of depth
sites = [
    {"id": "TP53:R175H", "alt_count": 12, "depth": 5000},
    {"id": "KRAS:G12D", "alt_count": 9, "depth": 4800},
    {"id": "APC:R1450*", "alt_count": 1, "depth": 5200},
    {"id": "PIK3CA:E545K", "alt_count": 10, "depth": 4900},
]

result = detect_minimal_residual_disease({"sites": sites}, options={"background_error_rate": 0.0005})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "detect_minimal_residual_disease",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 226,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d40e8a937f2b6f0485bbf2cd42cad1a6a7f68bf85efc96365e6faad9baa1ebba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "MRD detection: tumour-informed detection. The plasma ALT reads and depths at the tracked t...",
    "parameters": {"background_error_rate": 0.0005, "alpha": 0.001, "site_alpha": 0.01, "min_alt_reads": 2, "min_positive_sites": 2},
    "detected": true,
    "n_sites": 4,
    "n_positive_sites": 3,
    "pooled_alt_reads": 32,
    "pooled_depth": 19900,
    "pooled_vaf": 0.00160804,
    "pooled_p_value": 2.17455e-08,
    "estimated_tumor_fraction": 0.00221608,
    "sites": [
      {
        "id": "TP53:R175H",
        "alt_count": 12,
        "depth": 5000,
        "error_rate": 0.0005,
        "vaf": 0.0024,
        "p_value": 1.24946e-05,
        "positive": true
      },
      {
        "id": "KRAS:G12D",
        "alt_count": 9,
        "depth": 4800,
        "error_rate": 0.0005,
        "vaf": 0.001875,
        "p_value": 0.000858515,
        "positive": true
      },
      "..."
    ]
  }
}
```

## 11) generate_clinical_variant_report

- **Purpose:** Generate a clinical variant report (Markdown or plain text) from classified variants.

- **Method:** Each variant is given a tier: a supplied `classification` (`Pathogenic`, `Likely pathogenic`, `Uncertain significance`, `Likely benign`, `Benign`, or `P`/`LP`/`VUS`/`LB`/`B`) is used as is; when only `criteria` are given they are combined with the ACMG/AMP 2015 rules of `classify` (so contradictory criteria give Uncertain significance); otherwise the variant is `Not classified`. For VCF input the gene, ClinVar terms and consequence are read from the records and shown as evidence, but the ClinVar significance is not converted into a classification, so such variants are `Not classified` unless a classification is supplied; `hgvs_accession` (and `reference`) only adds an HGVS `g.` string. Variants are sorted by tier and gene, and the report lists the reportable findings (Pathogenic and Likely pathogenic), the variants of uncertain significance or not classified (unless `include_vus` is False), the number of benign/likely benign variants that are not listed, and a fixed statement of method and limitations. It summarises the evidence that was supplied; no database is queried and the criteria are not re-evaluated, so the text is not a validated diagnostic report and has to be reviewed by a qualified clinical geneticist.

- **API:** `dnakit.clinical.generate_clinical_variant_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'variants': [{'gene', 'hgvs', 'zygosity', 'classification' or 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}` (`sample`, `patient`, `patient_id`, `indication` and `date` are printed in the header), a list of such rows, or a VCF/records (gene, ClinVar and consequence fields are read from the records). `options`: `title`, `format` (`'markdown'` default or `'text'`), `include_vus` (True), `hgvs_accession` (string or `{contig: accession}`) and `reference` to build HGVS `g.` strings for VCF input, `output` (path to write the report), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `format`, `report` (the full report text), `n_variants`, `classification_counts`, `n_reportable`, `variants` (the rows with the resolved `classification`, and `classification_basis` when it was derived from criteria) and `output_path`. An unknown classification label or an empty variant list raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.clinical import generate_clinical_variant_report

result = generate_clinical_variant_report(
    {
        "sample": "P001",
        "indication": "hereditary cancer panel",
        "variants": [
            {"gene": "BRCA1", "hgvs": "NM_007294.4:c.68_69del", "zygosity": "heterozygous", "criteria": ["PVS1", "PS3"]},
            {"gene": "MSH6", "hgvs": "NM_000179.3:c.1A>G", "zygosity": "heterozygous", "classification": "VUS"},
            {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "LB"},
        ],
    },
    options={"title": "Hereditary cancer panel"},
)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "generate_clinical_variant_report",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 388,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3a1fb883b0a897f14eb0628ee788f198be4b6e220d8bea57bdbd15627135a355"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "report": "# Hereditary cancer panel\n\nSample: P001\nIndication: hereditary cancer panel\nVariants evalu...",
    "n_variants": 3,
    "classification_counts": {"Likely benign": 1, "Pathogenic": 1, "Uncertain significance": 1},
    "n_reportable": 1,
    "variants": [
      {
        "gene": "BRCA1",
        "hgvs": "NM_007294.4:c.68_69del",
        "zygosity": "heterozygous",
        "criteria": ["PVS1", "PS3"],
        "classification_basis": "ACMG/AMP 2015 combination of supplied criteria",
        "classification": "Pathogenic"
      },
      {
        "gene": "MSH6",
        "hgvs": "NM_000179.3:c.1A>G",
        "zygosity": "heterozygous",
        "classification": "Uncertain significance"
      },
      {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "Likely benign"}
    ],
    "output_path": null
  }
}
```

The `result["report"]` string printed in full:

    # Hereditary cancer panel

    Sample: P001
    Indication: hereditary cancer panel
    Variants evaluated: 3 (1 Likely benign, 1 Pathogenic, 1 Uncertain significance)

    ## Reportable findings (Pathogenic / Likely pathogenic)

    | Gene | Variant | Zygosity | Classification | Evidence |
    |---|---|---|---|---|
    | BRCA1 | NM_007294.4:c.68_69del | heterozygous | Pathogenic | PVS1, PS3 |

    ## Variants of uncertain significance / not classified

    | Gene | Variant | Zygosity | Classification | Evidence |
    |---|---|---|---|---|
    | MSH6 | NM_000179.3:c.1A>G | heterozygous | Uncertain significance | - |

    ## Not reported

    1 benign / likely benign variants and 0 uncertain variants are not listed.

    ## Method and limitations

    Classifications are taken from the input or computed from the supplied ACMG/AMP criteria with the 2015 combining rules (Richards et al., Genet Med 2015). The criteria were not re-evaluated against databases. This automatically generated summary is not a validated diagnostic report and must be reviewed by a qualified clinical geneticist.

## 12) infer_tumor_clones

- **Purpose:** Infer tumour clones by clustering mutation ALT fractions with a binomial mixture model (EM, number of clones chosen by BIC).

- **Method:** The (alt_count, depth) pairs are modelled as a mixture of k binomials fitted by EM from a deterministic quantile start; k = 1..`max_clones` (at most the number of mutations) is compared by BIC = -2 log L + (2k - 1) log n and the smallest BIC wins. Mutations are assigned to the component with the highest responsibility. Clusters are ranked by ALT fraction and the highest is labelled clonal, the others subclonal. The cancer cell fraction per mutation is CCF = VAF * (purity * CN + 2 * (1 - purity)) / (purity * multiplicity), with CN defaulting to 2 and multiplicity to 1 (both can be given per mutation), and each cluster reports the mean CCF. Copy-number changes are not inferred (they must be supplied), clone phylogeny/ordering is not modelled, and with few mutations or low depth the BIC choice is unreliable; it is not a replacement for PyClone, SciClone or similar tools.

- **API:** `dnakit.clinical.infer_tumor_clones(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: mutations `[{'alt_count', 'depth'[, 'id', 'copy_number', 'multiplicity']}, ...]` or `{'sites': [...]}`, or a VCF/records with `FORMAT/AD` (`options.sample`); at most 20,000 mutations. `options`: `purity` (1.0, 0.01-1), `max_clones` (4, at most 8), `copy_number` (default CN 2), `sample`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `purity`, `n_mutations`, `n_clones`, `model_selection` (`k`, `log_likelihood`, `bic` for each tried k), `clones` (per cluster: `cluster`, `role`, `cluster_vaf`, `pooled_vaf`, `weight`, `n_mutations`, `mean_ccf`, `mutations`) and `assignments` (`id`, `cluster`).

- **Example code:**

```python
import json

from dnakit.clinical import infer_tumor_clones

# (alt_count, depth) per mutation: a clonal group near VAF 0.40 and a subclonal group near 0.15
data = [(40, 100), (38, 100), (42, 100), (41, 100), (15, 100), (14, 100), (16, 100)]
sites = [{"id": f"m{i + 1}", "alt_count": a, "depth": d} for i, (a, d) in enumerate(data)]

result = infer_tumor_clones({"sites": sites}, options={"purity": 0.8})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "infer_tumor_clones",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 326,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b9fd33439bd546877ffc7c62a5baa0e90ee2b2d603c91c538205fc2ebcd69939"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "binomial mixture model of ALT/depth counts fitted by EM (deterministic quantile start), nu...",
    "purity": 0.8,
    "n_mutations": 7,
    "n_clones": 2,
    "model_selection": [{"k": 1, "log_likelihood": -424.161, "bic": 850.268}, {"k": 2, "log_likelihood": -401.198, "bic": 808.234}, "..."],
    "clones": [
      {
        "cluster_vaf": 0.4025,
        "pooled_vaf": 0.4025,
        "weight": 0.571429,
        "n_mutations": 4,
        "mean_ccf": 1.00625,
        "mutations": ["m1", "m2", "..."],
        "cluster": 1,
        "role": "clonal (highest ALT fraction)"
      },
      {
        "cluster_vaf": 0.15,
        "pooled_vaf": 0.15,
        "weight": 0.428571,
        "n_mutations": 3,
        "mean_ccf": 0.375,
        "mutations": ["m5", "m6", "m7"],
        "cluster": 2,
        "role": "subclonal"
      }
    ],
    "assignments": [{"id": "m1", "cluster": 1}, {"id": "m2", "cluster": 1}, "..."]
  }
}
```

## 13) predict

- **Purpose:** Predict clinical outcomes or variant effects with a trained model. Not implemented: it only reports what is missing (a trained clinical or pathogenicity predictor and its weights).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No model is loaded or run.

- **API:** `dnakit.clinical.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("trained clinical/pathogenicity prediction model") and `message` ("prediction needs trained model weights that are not available."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import predict

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = predict(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "clinical",
  "operation": "predict",
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
  "required_dependency": "trained clinical/pathogenicity prediction model",
  "message": "prediction needs trained model weights that are not available."
}
```

## 14) prioritize

- **Purpose:** Rank clinical candidate variants with a transparent additive score (consequence, ClinVar, rarity, CADD/REVEL, gene panel).

- **Method:** The same scoring as `dnakit.variants.prioritize_candidate_variants`. Annotations already on the records are read (consequence terms, `clinvar`/`CLNSIG`, `gnomad_af`, `cadd`, `revel`, `gene`); nothing is looked up. Score = consequence severity from a fixed Sequence Ontology table (0-1) + ClinVar (pathogenic 1.0, likely pathogenic 0.8, uncertain 0.1, benign -1.0) + rarity (0.5 when the population AF is absent or <= 1e-4, 0.3 up to 1e-3, 0.1 up to 1e-2) + 0.5 * max(CADD / 40, REVEL) + 0.5 for genes in the panel. Variants with population AF above `max_af`, and ClinVar benign/likely benign variants when `exclude_benign` is true, are listed as excluded. The weights are arbitrary, so this is a heuristic ranking and not a validated classifier; a variant with no population-frequency annotation keeps the full rarity score (0.5) but is flagged `frequency_unknown: true` and counted in one summary warning, so it can be told apart from a variant that is known to be rare. Every score component is returned.

- **API:** `dnakit.clinical.prioritize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). Annotations are read from record keys or INFO (`consequence`, `clinvar`/`CLNSIG`, `gnomad_af`, `cadd`, `revel`, `gene`). `options`: `max_af` (0.05; variants above it are excluded), `exclude_benign` (True), `genes` (panel list, +0.5), `top_n`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `parameters`, `ranked` (each with `chrom`, `pos`, `ref`, `alt`, `id`, `score`, `components`, `gene`, `rank`), `excluded` (with `excluded_reason`), `n_ranked`, `n_excluded` and `n_frequency_unknown` (number of variants, ranked or excluded, without a population-frequency annotation; the envelope `warnings` then holds one message with that count). Every ranked and excluded row has the boolean `frequency_unknown`.

- **Example code:**

```python
import json

from dnakit.clinical import prioritize

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=BRCA1;consequence=stop_gained;clinvar=Pathogenic;gnomad_af=0.00001;cadd=35",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=TTN;consequence=missense_variant;gnomad_af=0.002;revel=0.6",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=OBSCN;consequence=synonymous_variant;gnomad_af=0.2",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=APOE;consequence=missense_variant;clinvar=Benign",
]) + "\n"

result = prioritize(vcf, options={"genes": ["BRCA1"], "top_n": 3})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "prioritize",
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

## 15) qc

- **Purpose:** Quality-control a clinical call set: per-sample call rate, het/hom ratio, Ti/Tv and depth, per-site call rate and exact HWE test, with threshold flags.

- **Method:** Per-sample metrics are those of `dnakit.variants.summarize` (call rate, hom-ref/het/hom-alt counts, het/hom-alt ratio, Ti/Tv, mean depth, mean GQ); per-site metrics are call rate, allele counts/frequencies and the exact HWE p-value of `calculate_hwe_pvalue`. A sample is flagged `low_call_rate` below `min_sample_call_rate`, and, only when the option is given, `low_mean_depth`, `high_het_hom_ratio` or `ti_tv_out_of_range`; a site is flagged `low_call_rate` below `min_site_call_rate` or `hwe_violation` below `min_hwe_pvalue`. `passed` is true when no sample is flagged (site flags do not affect it). This covers called genotypes only, not raw sequencing-run metrics (base quality, duplication, contamination), and the thresholds are generic defaults, not laboratory-validated criteria.

- **API:** `dnakit.clinical.qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: `min_sample_call_rate` (0.95), `min_site_call_rate` (0.9), `min_hwe_pvalue` (1e-6), `min_mean_depth`, `max_het_hom_alt_ratio`, `ti_tv_range` (`[low, high]`); the last three are off unless given.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `method`, `thresholds`, `n_sites`, `n_samples`, `samples` (per-sample metrics plus `flags`), `samples_failing`, `n_sites_flagged`, `flagged_sites` (at most 1000, with per-site metrics and `flags`) and `passed`.

- **Example code:**

```python
import json

from dnakit.clinical import qc

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/0:10\t0/1:12\t1/1:11",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP\t0/1:9\t0/1:10\t./.:0",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:DP\t0/0:10\t./.:0\t./.:0",
    "chr1\t400\t.\tT\tC\t50\tPASS\t.\tGT:DP\t0/1:11\t1/1:12\t0/0:10",
]) + "\n"

result = qc(vcf, options={"min_sample_call_rate": 0.8})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "qc",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 280,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9a7e1e3b47679b6be206402306d35eae65b6c2089bdcde71517c09d310d2e168"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "call-set QC: per-sample call rate, het/hom-alt ratio, Ti/Tv and depth; per-site call rate...",
    "thresholds": {
      "min_sample_call_rate": 0.8,
      "min_site_call_rate": 0.9,
      "min_hwe_pvalue": 1e-06,
      "min_mean_depth": null,
      "max_het_hom_alt_ratio": null,
      "ti_tv_range": null
    },
    "n_sites": 4,
    "n_samples": 3,
    "samples": [
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
        "n_snv": 2,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": null,
        "flags": []
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 0,
        "n_het": 2,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 2.0,
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 8.5,
        "mean_genotype_quality": null,
        "flags": ["low_call_rate"]
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 2,
        "n_missing": 2,
        "call_rate": 0.5,
        "n_hom_ref": 1,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 5.25,
        "mean_genotype_quality": null,
        "flags": ["low_call_rate"]
      }
    ],
    "samples_failing": ["S2", "S3"],
    "n_sites_flagged": 2,
    "flagged_sites": [
      {
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": ["T"],
        "n_samples": 3,
        "call_rate": 0.666667,
        "an": 4,
        "ac": [2],
        "af": [0.5],
        "n_het": 2,
        "n_hom_alt": 0,
        "hwe_pvalue": 1.0,
        "mean_depth": 6.33333,
        "qual": 50.0,
        "filter": "PASS",
        "flags": ["low_call_rate"]
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": ["A"],
        "n_samples": 3,
        "call_rate": 0.333333,
        "an": 2,
        "ac": [0],
        "af": [0.0],
        "n_het": 0,
        "n_hom_alt": 0,
        "hwe_pvalue": 1.0,
        "mean_depth": 3.33333,
        "qual": 50.0,
        "filter": "PASS",
        "flags": ["low_call_rate"]
      }
    ],
    "passed": false
  }
}
```

## 16) read

- **Purpose:** Read a (ClinVar-style) VCF or variant records into clinical fields: gene, ClinVar significance and review status, disease, consequence and genotypes.

- **Method:** A parser, not an interpreter. For each record it reads INFO `GENEINFO`/`GENE` (gene symbol), `CLNSIG` (kept raw in `clinvar_significance_raw` and normalised to lower-case underscore terms in `clinvar_terms`, for example `likely_pathogenic`), `CLNREVSTAT`, `CLNDN` (split at `|` into diseases) and `MC` (Sequence Ontology term, or a `consequence` key), plus the genotype of every sample. Record keys of the same name (`gene`, `clinvar`, `consequence`, ...) are accepted for list-of-mappings input. No external look-up is made and no classification is derived.

- **API:** `dnakit.clinical.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: VCF text (tab-separated, optionally gzip bytes), a path to a VCF or VCF.gz, or a list of variant mappings (`chrom`, `pos`, `ref`, `alt`, optional `qual`, `filter`, `info`, and `genotypes`/`samples`). `options`: none.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `samples`, `n_variants`, `method` and `variants` (per record: `chrom`, `pos`, `id`, `ref`, `alt`, `gene`, `clinvar_significance_raw`, `clinvar_terms`, `review_status`, `disease`, `consequence`, `filter`, and `genotypes` when the VCF has samples).

- **Example code:**

```python
import json

from dnakit.clinical import read

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr17\t43094464\trs80357906\tC\tCT\t.\t.\tGENEINFO=BRCA1:672;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Hereditary_breast_ovarian_cancer_syndrome;MC=SO:0001589|frameshift_variant",
    "chr2\t47806520\trs1\tA\tG\t.\t.\tGENEINFO=MSH6:2956;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Lynch_syndrome|not_provided;MC=SO:0001583|missense_variant",
]) + "\n"

result = read(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "read",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 436,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b9a5b4d38401d514fb565af6b5e4ea890fab99311f3ea90124d1512a3de8c6d0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [],
    "n_variants": 2,
    "variants": [
      {
        "chrom": "chr17",
        "pos": 43094464,
        "id": "rs80357906",
        "ref": "C",
        "alt": ["CT"],
        "gene": "BRCA1",
        "clinvar_significance_raw": "Pathogenic",
        "clinvar_terms": ["pathogenic"],
        "review_status": "reviewed_by_expert_panel",
        "disease": ["Hereditary_breast_ovarian_cancer_syndrome"],
        "consequence": ["frameshift_variant"],
        "filter": "."
      },
      {
        "chrom": "chr2",
        "pos": 47806520,
        "id": "rs1",
        "ref": "A",
        "alt": ["G"],
        "gene": "MSH6",
        "clinvar_significance_raw": "Uncertain_significance",
        "clinvar_terms": ["uncertain_significance"],
        "review_status": "criteria_provided,_single_submitter",
        "disease": ["Lynch_syndrome", "not_provided"],
        "consequence": ["missense_variant"],
        "filter": "."
      }
    ],
    "method": "VCF/records parsed into clinical fields: ClinVar CLNSIG (normalised to lower-case terms),..."
  }
}
```

## 17) report

- **Purpose:** Generate a clinical variant report (Markdown or plain text) from classified variants.

- **Method:** Each variant is given a tier: a supplied `classification` (`Pathogenic`, `Likely pathogenic`, `Uncertain significance`, `Likely benign`, `Benign`, or `P`/`LP`/`VUS`/`LB`/`B`) is used as is; when only `criteria` are given they are combined with the ACMG/AMP 2015 rules of `classify` (so contradictory criteria give Uncertain significance); otherwise the variant is `Not classified`. For VCF input the gene, ClinVar terms and consequence are read from the records and shown as evidence, but the ClinVar significance is not converted into a classification, so such variants are `Not classified` unless a classification is supplied; `hgvs_accession` (and `reference`) only adds an HGVS `g.` string. Variants are sorted by tier and gene, and the report lists the reportable findings (Pathogenic and Likely pathogenic), the variants of uncertain significance or not classified (unless `include_vus` is False), the number of benign/likely benign variants that are not listed, and a fixed statement of method and limitations. It summarises the evidence that was supplied; no database is queried and the criteria are not re-evaluated, so the text is not a validated diagnostic report and has to be reviewed by a qualified clinical geneticist.

- **API:** `dnakit.clinical.report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'variants': [{'gene', 'hgvs', 'zygosity', 'classification' or 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}` (`sample`, `patient`, `patient_id`, `indication` and `date` are printed in the header), a list of such rows, or a VCF/records (gene, ClinVar and consequence fields are read from the records). `options`: `title`, `format` (`'markdown'` default or `'text'`), `include_vus` (True), `hgvs_accession` (string or `{contig: accession}`) and `reference` to build HGVS `g.` strings for VCF input, `output` (path to write the report), `overwrite`.

- **Output:** Returns the standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` has `format`, `report` (the full report text), `n_variants`, `classification_counts`, `n_reportable`, `variants` (the rows with the resolved `classification`, and `classification_basis` when it was derived from criteria) and `output_path`. An unknown classification label or an empty variant list raises `InputFormatError`.

- **Example code:**

```python
import json

from dnakit.clinical import report

rows = [
    {"gene": "CFTR", "hgvs": "NM_000492.4:c.1521_1523del", "zygosity": "homozygous", "classification": "P"},
    {"gene": "SCN5A", "hgvs": "NM_198056.3:c.5A>G", "zygosity": "heterozygous", "criteria": ["PM2", "PP3", "PP2"]},
    {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "B"},
]

result = report(rows, options={"format": "text", "include_vus": False, "title": "Cardiac panel"})

print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings are shortened with `...` and numbers are rounded to 6 significant digits; the output was produced by running the code above.

```json
{
  "module": "clinical",
  "operation": "report",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 320,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "083c03cb00fc96f648d730193daf78d3db4b45b82aff7d4b708e50d754ecbe62"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "text",
    "report": "Cardiac panel\n\nVariants evaluated: 3 (1 Benign, 1 Pathogenic, 1 Uncertain significance)\n\nR...",
    "n_variants": 3,
    "classification_counts": {"Benign": 1, "Pathogenic": 1, "Uncertain significance": 1},
    "n_reportable": 1,
    "variants": [
      {"gene": "CFTR", "hgvs": "NM_000492.4:c.1521_1523del", "zygosity": "homozygous", "classification": "Pathogenic"},
      {
        "gene": "SCN5A",
        "hgvs": "NM_198056.3:c.5A>G",
        "zygosity": "heterozygous",
        "criteria": ["PM2", "PP3", "PP2"],
        "classification_basis": "ACMG/AMP 2015 combination of supplied criteria",
        "classification": "Uncertain significance"
      },
      {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "Benign"}
    ],
    "output_path": null
  }
}
```

The `result["report"]` string printed in full:

```text
Cardiac panel

Variants evaluated: 3 (1 Benign, 1 Pathogenic, 1 Uncertain significance)

Reportable findings (Pathogenic / Likely pathogenic)

Gene=CFTR; Variant=NM_000492.4:c.1521_1523del; Zygosity=homozygous; Classification=Pathogenic; Evidence=-

Not reported

1 benign / likely benign variants and 1 uncertain variants are not listed.

Method and limitations

Classifications are taken from the input or computed from the supplied ACMG/AMP criteria with the 2015 combining rules (Richards et al., Genet Med 2015). The criteria were not re-evaluated against databases. This automatically generated summary is not a validated diagnostic report and must be reviewed by a qualified clinical geneticist.
```

## 18) type_hla

- **Purpose:** Type HLA alleles from sequencing reads. Not implemented: it only reports what is missing (an HLA typing program (OptiType, HLA-LA, arcasHLA, T1K), the IPD-IMGT/HLA reference and reads).

- **Method:** No computation is performed and no result is invented. The call validates nothing about the data and returns an honest `not_implemented` status. No external program is run.

- **API:** `dnakit.clinical.type_hla(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any value; it is only summarized in `input_summary`. `options` is accepted and ignored.

- **Output:** Always returns `status: "not_implemented"` with `required_dependency` ("HLA typing program (OptiType/HLA-LA/arcasHLA) with IPD-IMGT/HLA reference and reads") and `message` ("HLA typing needs an external program and the IMGT/HLA reference."); there is no `result`. The envelope also carries `module`, `operation`, `provenance` (`implementation: "not_implemented"`), `input_summary` and `warnings`.

- **Example code:**

```python
import json

from dnakit.clinical import type_hla

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = type_hla(vcf)

print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "clinical",
  "operation": "type_hla",
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
  "required_dependency": "HLA typing program (OptiType/HLA-LA/arcasHLA) with IPD-IMGT/HLA reference and reads",
  "message": "HLA typing needs an external program and the IMGT/HLA reference."
}
```
