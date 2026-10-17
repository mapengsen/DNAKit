# Forensic and reproductive genomics

This category contains **9** documented capabilities. **7** compute results natively from allele frequencies and single-source STR/SNP profiles using the standard library: `calculate` (random match probability, likelihood ratio, paternity index, kinship index), `compare`, `infer` (relationship posterior), `plot`, `read`, `report` and `validate`. **2** return `status: "not_implemented"` in every call because they need validated external models or software: `deconvolute` (probabilistic genotyping of mixtures) and `predict` (phenotype and ancestry prediction). There is no mutation, drop-in/dropout or mixture model. The functions are exported through the DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) calculate

- **Purpose:** Calculate a forensic statistic from allele frequencies: random match probability, likelihood ratio, paternity index or kinship index.

- **Method:** `random_match_probability`: the product rule over loci, p^2 for a homozygote and 2pq for a heterozygote; with `theta` > 0 the NRC II recommendation 4.10 formulas [2θ+(1-θ)p][3θ+(1-θ)p]/[(1+θ)(1+2θ)] and 2[θ+(1-θ)pi][θ+(1-θ)pj]/[(1+θ)(1+2θ)]. `likelihood_ratio`: Hp the suspect is the source, Hd an unrelated random person is; LR = 1/P(genotype) per matching locus, 0 on a mismatch, multiplied over loci (`theta` supported). `paternity_index`: PI = P(child | mother, alleged father) / P(child | mother, random man), each parent transmitting either allele with probability 1/2, combined index CPI = product of PI and probability of paternity W = CPI*prior/(CPI*prior + 1 - prior); duo (no mother) or trio; no mutation model, so one exclusion gives CPI = 0; theta = 0. `kinship_index`: LR = P(g2 | g1, relationship)/P(g2) with P(g2 | g1) = k0 P(g2) + k1 P(g2 | one IBD allele) + k2 [g2 = g1], where (k0, k1, k2) are the Cotterman IBD probabilities of the relationship (theta = 0, no mutation). Verbal equivalents use the ENFSI 2015 scale. Pure standard library.

- **API:** `dnakit.forensic_reproductive.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `options['quantity']`: `random_match_probability` (default; input `{'profile', 'frequencies'}`; `theta` 0 to 0.5), `likelihood_ratio` (`{'evidence', 'suspect', 'frequencies'}`; `theta`), `paternity_index` (`{'child', 'alleged_father', 'mother' (optional), 'frequencies'}`; `prior` 0.5) or `kinship_index` (`{'profile_1', 'profile_2', 'frequencies'}`; `relationship` = `parent_child`, `full_sibling` (default), `half_sibling`, `grandparent_grandchild`, `uncle_nephew`, `first_cousin`, `identical_twin`, or custom `ibd=[k0, k1, k2]` summing to 1). `frequencies` may also be passed in `options`. An allele missing from the frequency table raises an error unless `min_frequency` is set. Profiles are `{locus: [allele1, allele2]}` (a single allele is treated as homozygous); allele frequencies are `{locus: {allele: frequency}}` or records with columns `locus`, `allele`, `frequency`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` always holds `quantity`, `per_locus` (one record per locus), `combined`, `combined_log10` and `formula`. Random match probability: `theta`, per-locus `genotype_probability` and `log10`, `one_in`. Likelihood ratio: `theta`, per-locus `match`, `lr`, `log10_lr`, `excluded`, `verbal_equivalent`. Paternity index: `design` (`duo` or `trio`), per-locus `numerator`, `denominator`, `pi`, `log10_pi`, `exclusions`, `prior`, `probability_of_paternity`, `verbal_equivalent`. Kinship index: `relationship`, `ibd`, per-locus `lr`, `log10_lr`, `verbal_equivalent`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import calculate

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = calculate(
    {"profile": profile, "frequencies": frequencies},
    options={"quantity": "random_match_probability"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 264,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "604f5e370a70a7f4e8a5dad7a9412e6c919b171c7c709b6d9e46d955bf60fed1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "random_match_probability",
    "theta": 0.0,
    "per_locus": [
      {
        "locus": "D3S1358",
        "genotype": [
          "15",
          "16"
        ],
        "genotype_probability": 0.15,
        "log10": -0.8239087409443188
      },
      {
        "locus": "vWA",
        "genotype": [
          "17",
          "17"
        ],
        "genotype_probability": 0.12249999999999998,
        "log10": -0.9118639112994488
      },
      {
        "locus": "FGA",
        "genotype": [
          "22",
          "24"
        ],
        "genotype_probability": 0.12,
        "log10": -0.9208187539523752
      }
    ],
    "combined": 0.002205,
    "combined_log10": -2.6565914061961426,
    "one_in": 453.51473922902494,
    "formula": "product rule: homozygote p^2, heterozygote 2pq; with theta>0 the NRC II 4.10 formulas [2θ+(1-θ)p][3θ+(1-θ)p]/[(1+θ)(1+2θ)] and 2[θ+(1-θ)pi][θ+(1-θ)pj]/[(1+θ)(1+2θ)]"
  }
}
```

## 2) compare

- **Purpose:** Compare two single-source genotype profiles locus by locus (identical, one allele shared, no allele shared).

- **Method:** For every locus present in both profiles, the genotypes are compared as unordered pairs: identical (2 alleles shared), one allele shared (counting alleles with multiplicity) or none. The profiles are `consistent_with_same_source` when every compared locus is identical; no statistical weight is computed unless `frequencies` is given, in which case the random match probability of profile 1 over the compared loci is added.

- **API:** `dnakit.forensic_reproductive.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'profile_1': {...}, 'profile_2': {...}}` and optionally `'frequencies'` (or `options['frequencies']`). Profiles are `{locus: [allele1, allele2]}` (a single allele is treated as homozygous); allele frequencies are `{locus: {allele: frequency}}` or records with columns `locus`, `allele`, `frequency`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `per_locus` (`locus`, both genotypes, `status`, `alleles_shared`), `n_loci_compared`, the counts `identical`, `one_allele_shared`, `no_allele_shared`, `loci_only_in_profile_1`, `loci_only_in_profile_2`, `consistent_with_same_source`, `mean_alleles_shared` and, with frequencies, `random_match_probability_profile_1`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import compare

profile_1 = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
profile_2 = {"D3S1358": ["15", "16"], "vWA": ["17", "18"], "FGA": ["21", "23"]}
result = compare({"profile_1": profile_1, "profile_2": profile_2})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 164,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9a2cd44364446204b1f0719418ea04084042426d9b10b0dcd968f5421393ffe4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "per_locus": [
      {
        "locus": "D3S1358",
        "profile_1": [
          "15",
          "16"
        ],
        "profile_2": [
          "15",
          "16"
        ],
        "status": "identical",
        "alleles_shared": 2
      },
      {
        "locus": "vWA",
        "profile_1": [
          "17",
          "17"
        ],
        "profile_2": [
          "17",
          "18"
        ],
        "status": "one_allele_shared",
        "alleles_shared": 1
      },
      {
        "locus": "FGA",
        "profile_1": [
          "22",
          "24"
        ],
        "profile_2": [
          "21",
          "23"
        ],
        "status": "no_allele_shared",
        "alleles_shared": 0
      }
    ],
    "n_loci_compared": 3,
    "identical": 1,
    "one_allele_shared": 1,
    "no_allele_shared": 1,
    "loci_only_in_profile_1": [],
    "loci_only_in_profile_2": [],
    "consistent_with_same_source": false,
    "mean_alleles_shared": 1.0
  }
}
```

## 3) deconvolute

- **Purpose:** Deconvolute a DNA mixture into contributor genotypes. Not implemented in DNAKit.

- **Method:** No computation is done. Mixture deconvolution needs a validated probabilistic genotyping model with peak heights, drop-in/dropout and stutter parameters (for example EuroForMix or STRmix); the function only reports this requirement.

- **API:** `dnakit.forensic_reproductive.deconvolute(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any `input` and optional `options`; they are only summarised in `input_summary`.

- **Output:** Does not compute a result. A normal call returns the envelope with `status: "not_implemented"` (`provenance.implementation` is `not_implemented`), plus `required_dependency` and `message` stating what is missing: validated probabilistic genotyping software (e.g. EuroForMix, STRmix) and peak-height data.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import deconvolute

mixture = {"D3S1358": ["15", "16", "17"], "vWA": ["16", "17", "18"]}
result = deconvolute({"mixture": mixture})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "deconvolute",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 71,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "249f54c134c4b9460798aba692cf46ed0ed2bfc4fff3969693bb051e5987f4ef"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated probabilistic genotyping software (e.g. EuroForMix, STRmix) and peak-height data",
  "message": "Mixture deconvolution needs a validated probabilistic genotyping model with peak heights, drop-in/dropout and stutter parameters; DNAKit does not implement it."
}
```

## 4) infer

- **Purpose:** Infer the most probable biological relationship between two profiles from IBD-model likelihood ratios.

- **Method:** For every candidate relationship the kinship likelihood ratio against `unrelated` is computed from the Cotterman IBD probabilities (k0, k1, k2) (theta = 0, no mutation), and posterior(r) is proportional to prior(r) * LR(r vs unrelated), normalised over the candidates; priors default to uniform. Relationships with identical (k0, k1, k2) (half_sibling, grandparent_grandchild, uncle_nephew) cannot be distinguished by this model: their likelihood ratios are always equal (their posteriors differ only through the priors), and `most_probable` among them reflects list order, not evidence. Whenever two or more of the evaluated relationships share the same (k0, k1, k2), one warning in `warnings` lists them. An error is raised when the profiles are inconsistent with every candidate.

- **API:** `dnakit.forensic_reproductive.infer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'profile_1', 'profile_2', 'frequencies'}`. `options`: `relationships` (default `unrelated`, `parent_child`, `full_sibling`, `half_sibling`, `first_cousin`; any names accepted by `calculate`), `priors` ({relationship: prior}, one per listed relationship), `min_frequency`. Profiles are `{locus: [allele1, allele2]}` (a single allele is treated as homozygous); allele frequencies are `{locus: {allele: frequency}}` or records with columns `locus`, `allele`, `frequency`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `relationships` (sorted by posterior; each with `relationship`, `likelihood_ratio_vs_unrelated`, `log10_lr`, `prior`, `posterior`), `most_probable` and `method`; `warnings` holds one message per group of mathematically indistinguishable relationships among those evaluated. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import infer

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile_1 = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
profile_2 = {"D3S1358": ["15", "17"], "vWA": ["17", "18"], "FGA": ["24", "23"]}
result = infer(
    {"profile_1": profile_1, "profile_2": profile_2, "frequencies": frequencies},
    options={"relationships": ["unrelated", "parent_child", "full_sibling"]},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "infer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 348,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3adbb8c661a53b053f3e78105cc3481951e6d5846c902b6a1ba2d867d3e90554"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "relationships": [
      {
        "relationship": "parent_child",
        "likelihood_ratio_vs_unrelated": 1.7857142857142858,
        "log10_lr": 0.2518119729937996,
        "prior": 0.3333333333333333,
        "posterior": 0.5223636957231472
      },
      {
        "relationship": "unrelated",
        "likelihood_ratio_vs_unrelated": 1.0,
        "log10_lr": 0.0,
        "prior": 0.3333333333333333,
        "posterior": 0.29252366960496246
      },
      {
        "relationship": "full_sibling",
        "likelihood_ratio_vs_unrelated": 0.6328125000000001,
        "log10_lr": -0.19872495076921853,
        "prior": 0.3333333333333333,
        "posterior": 0.18511263467189032
      }
    ],
    "most_probable": "parent_child",
    "method": "posterior(r) ∝ prior(r) * LR(r vs unrelated), LR from the Cotterman IBD model (theta = 0, no mutation); priors default to uniform"
  }
}
```

## 5) plot

- **Purpose:** Draw the per-locus log10 genotype probability or likelihood ratio of a forensic statistic as an SVG bar chart.

- **Method:** Computes the statistic exactly as `calculate` does and draws one bar per locus (log10 genotype probability for `random_match_probability`, log10 LR / PI otherwise) on a fixed 720 x 400 SVG with a zero line, tick labels and rotated locus names. Loci with a zero LR (exclusions) have no bar. Nothing is written unless `path` is given.

- **API:** `dnakit.forensic_reproductive.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` and `options` as for `calculate` (`quantity` selects the statistic); `path` writes the SVG (an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `quantity`, `svg` (the SVG text), `n_loci` and `path` (null unless written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import plot

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
profile = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = plot(
    {"profile": profile, "frequencies": frequencies},
    options={"quantity": "random_match_probability"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "forensic_reproductive",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 264,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "604f5e370a70a7f4e8a5dad7a9412e6c919b171c7c709b6d9e46d955bf60fed1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "random_match_probability",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"400\" viewBox=\"0 0 720 400\"><title>random_match_probability p...",
    "n_loci": 3,
    "path": null
  }
}
```

## 6) predict

- **Purpose:** Predict phenotype or ancestry from a DNA profile. Not implemented in DNAKit.

- **Method:** No computation is done. Forensic DNA phenotyping and ancestry prediction require validated trained models and reference populations (for example HIrisPlex-S or a reference panel), which are not bundled; the function only reports this requirement.

- **API:** `dnakit.forensic_reproductive.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any `input` and optional `options`; they are only summarised in `input_summary`.

- **Output:** Does not compute a result. A normal call returns the envelope with `status: "not_implemented"` (`provenance.implementation` is `not_implemented`), plus `required_dependency` and `message` stating what is missing: validated forensic phenotype/ancestry prediction models and reference data (e.g. HIrisPlex-S, a reference panel).

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import predict

result = predict({"profile": {"D3S1358": ["15", "16"], "vWA": ["17", "17"]}})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "36756aa017f6c6033830b17c63bfe422f26f28d09aaa87be8526a8ab01acb794"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated forensic phenotype/ancestry prediction models and reference data (e.g. HIrisPlex-S, a reference panel)",
  "message": "Forensic DNA phenotyping and ancestry prediction require validated trained models and reference populations that are not bundled."
}
```

## 7) read

- **Purpose:** Read a genotype profile or an allele-frequency table (path, delimited text, records or mappings) into the DNAKit structures.

- **Method:** Delimited text, files, records or column mappings are loaded as a table. With `kind='auto'` the table is a frequency table when it has a `frequency`/`freq` column, otherwise a profile. Profile columns are `locus`, `allele1`, `allele2`; frequency columns are `locus`, `allele`, `frequency`. A mapping with `profile` and/or `frequencies` keys is normalised directly (allele labels become strings). Pure standard library.

- **API:** `dnakit.forensic_reproductive.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a path, delimited text, records or a mapping with `profile` / `frequencies`. `options`: `kind` = `auto` (default), `profile` or `frequencies`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `kind`, `profile` ({locus: [allele1, allele2]}) or `frequencies` ({locus: {allele: frequency}}) and `n_loci`; a mapping input returns just the normalised `profile` and/or `frequencies` entries. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import read

text = "locus,allele1,allele2\nD3S1358,15,16\nvWA,17,17\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "forensic_reproductive",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 46,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8491b2344053456425e5c6a64b5e7257c361b1691765ccfd9662f0422ca24a98"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "profile",
    "profile": {
      "D3S1358": [
        "15",
        "16"
      ],
      "vWA": [
        "17",
        "17"
      ]
    },
    "n_loci": 2
  }
}
```

## 8) report

- **Purpose:** Generate a Markdown statistical report (per-locus table, combined value, ENFSI verbal scale, limitations) for a forensic statistic.

- **Method:** Computes the statistic as `calculate` does and formats it as Markdown: title, statistic and design, the formula used, a table of per-locus values with log10, the combined value, the random match probability as 1 in N or the ENFSI 2015 verbal equivalent (plus probability of paternity and exclusions for a paternity index) and a fixed limitations paragraph (single-source profiles, no mutation, drop-in/dropout or mixture model). It is a formatting aid, not a validated casework report.

- **API:** `dnakit.forensic_reproductive.report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` and `options` as for `calculate` (`quantity` selects the statistic); `title` (default `Forensic statistical report`); `path` writes the Markdown (an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `format` (`markdown`), `text` (the report), `statistic` (the full result of `calculate`) and `path` (null unless written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import report

frequencies = {
    "D3S1358": {"15": 0.25, "16": 0.30, "17": 0.25, "18": 0.20},
    "vWA": {"14": 0.10, "16": 0.30, "17": 0.35, "18": 0.25},
    "FGA": {"21": 0.20, "22": 0.30, "23": 0.30, "24": 0.20},
}
evidence = {"D3S1358": ["15", "16"], "vWA": ["17", "17"], "FGA": ["22", "24"]}
result = report(
    {"evidence": evidence, "suspect": evidence, "frequencies": frequencies},
    options={"quantity": "likelihood_ratio", "title": "Case 1"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "forensic_reproductive",
  "operation": "report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 345,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "95080364deae6a65d99264fa39751ecd91ce28574936d2073d6f7aca4effb96c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "text": "# Case 1\n\nStatistic: **likelihood ratio**\n\nMethod: Hp: the suspect is the source of the evidence; Hd: an unrelated rando...",
    "statistic": {
      "quantity": "likelihood_ratio",
      "theta": 0.0,
      "per_locus": [
        {
          "locus": "D3S1358",
          "match": true,
          "lr": 6.666666666666667,
          "log10_lr": 0.8239087409443188
        },
        {
          "locus": "vWA",
          "match": true,
          "lr": 8.163265306122451,
          "log10_lr": 0.9118639112994489
        },
        {
          "locus": "FGA",
          "match": true,
          "lr": 8.333333333333334,
          "log10_lr": 0.9208187539523752
        }
      ],
      "combined": 453.5147392290254,
      "combined_log10": 2.656591406196143,
      "excluded": false,
      "verbal_equivalent": "moderately strong support for the first (prosecution/relationship) proposition",
      "formula": "Hp: the suspect is the source of the evidence; Hd: an unrelated random person is the source. Per locus LR = 1/P(genotype) when the genotypes match, 0 otherwise; combined = product over loci"
    },
    "path": null
  }
}
```

## 9) validate

- **Purpose:** Validate forensic profiles and allele-frequency tables: frequency sums, allele counts per locus, and loci or alleles missing from the frequency table.

- **Method:** Errors: the frequencies of a locus do not sum to 1 within `tolerance`; a locus has more than two alleles (a mixture, not a single-source profile); a locus is absent from the frequency table. Warnings: a single allele is listed (treated as homozygous, possible dropout); an allele is absent from the frequency table (calculations then need `min_frequency`). Content problems give `valid: false` rather than raising.

- **API:** `dnakit.forensic_reproductive.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping of any of `frequencies`, `profile`, `evidence`, `suspect`, `child`, ... (every profile-like entry is checked). `options`: `tolerance` (0.02).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `valid`, `n_errors`, `n_warnings` and `issues` (records with `severity`, `where` and `message`). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.forensic_reproductive import validate

frequencies = {"D3S1358": {"15": 0.5, "16": 0.3}, "vWA": {"17": 0.5, "18": 0.5}}
profile = {"D3S1358": ["15", "16", "17"], "vWA": ["17", "19"], "FGA": ["22", "24"]}
result = validate({"frequencies": frequencies, "profile": profile})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "forensic_reproductive",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 169,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1d22dcd2d8b8e87b5a0373b7cc09edab13b7cbd26fa968b9172c8e9b52a50855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_errors": 3,
    "n_warnings": 1,
    "issues": [
      {
        "severity": "error",
        "where": "D3S1358",
        "message": "allele frequencies sum to 0.8000, not 1 (tolerance 0.02)"
      },
      {
        "severity": "error",
        "where": "profile:D3S1358",
        "message": "3 alleles: a single-source profile has at most two"
      },
      "..."
    ]
  }
}
```
