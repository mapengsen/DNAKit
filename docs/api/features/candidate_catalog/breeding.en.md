# Breeding

This category contains **7** documented capabilities, all computed natively in DNAKit (none returns `not_implemented` in a normal call): breeder's-equation and heritability calculations, pedigree and genomic relationships, ridge/GBLUP fitting and prediction, truncation and Smith-Hazel selection, optimum-contribution optimisation, a forward selection simulation and an SVG progress plot. `fit`, `predict` and `calculate(quantity='genomic_relationship')` need numpy and return `status: "not_implemented"` without it. The functions are exported through the DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) calculate

- **Purpose:** Calculate a quantitative-genetics quantity: the breeder's-equation response to selection, heritability, the pedigree (additive) relationship matrix or the VanRaden genomic relationship matrix.

- **Method:** `options['quantity']` selects the quantity. `breeders_equation` (default): R = h2 * S with S = i * sigma_p; the selection intensity of truncation selection is i = phi(z)/p with z = Phi^-1(1-p) (infinite population, normal trait), and `response_per_year` = R / generation_interval when an interval is given. `heritability`: the slope of offspring on mid-parent (h2 = b) or on a single parent (h2 = 2b) from parent/offspring lists, or the realised heritability R / S, or Vg / Vp. `pedigree_relationship`: Henderson's tabular method, A_ij = (A_sj + A_dj)/2, A_ii = 1 + A_sd/2, inbreeding F = A_ii - 1. `genomic_relationship`: the VanRaden GRM (numpy needed; `not_implemented` without it).

- **API:** `dnakit.breeding.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` plus `options`. `breeders_equation`: a mapping with `heritability` (0 to 1) and either `selection_differential`, or `selection_intensity` / `proportion_selected` with `phenotypic_sd`; optional `generation_interval` (values may also be given in `options`). `heritability`: a mapping with `parent` and `offspring` lists (at least 3 pairs; `options['relationship']` = `midparent` (default) or `single_parent`), or `response` with `selection_differential`, or `genetic_variance` with `phenotypic_variance`. `pedigree_relationship`: records with `id`, `sire`, `dam` (or a mapping with `pedigree`). `genomic_relationship`: a genotype matrix (samples x markers dosages). `options['quantity']` defaults to `breeders_equation`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. For `breeders_equation`, `result` holds `quantity`, `heritability`, `selection_differential`, `selection_intensity` (null when not derivable), `response`, `response_per_year` (only with a generation interval) and `definition`. The other quantities return `h2`, `slope`, `slope_se`, `n_pairs` and `method` (heritability), `ids`, `matrix` and `inbreeding` ({id: F}) (pedigree), or `sample_ids`, `relationship_matrix`, `method`, `n_markers_used`, `allele_frequencies` and marker-filtering details (genomic). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import calculate

result = calculate(
    {"heritability": 0.5, "selection_differential": 2.0},
    options={"quantity": "breeders_equation"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "breeding",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 52,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d78fc459837b52d2ed27a19a773b13dd1a63f39d2106ec21454b9fa6df84268"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "breeders_equation",
    "heritability": 0.5,
    "selection_differential": 2.0,
    "selection_intensity": null,
    "response": 1.0,
    "definition": "R = h2 * S, S = i * sigma_p; truncation intensity i = phi(z)/p with z = Phi^-1(1-p) (infinite population, normal trait)"
  }
}
```

## 2) fit

- **Purpose:** Fit a genomic ridge-regression (SNP-BLUP, equivalent to GBLUP) model of a phenotype on marker dosages.

- **Method:** Ridge regression of the phenotype on centred dosages, minimising sum((y - mu - Zb)^2) + lambda * sum(b^2); this equals GBLUP with G = ZZ'/(2*sum(pq)) when lambda = 2*sum(pq)*(1-h2)/h2, with p estimated from the training sample. Breeding values are the fitted values minus the intercept. The training correlation is computed on the training data and is therefore optimistic (a warning says so). Needs numpy; returns `not_implemented` without it.

- **API:** `dnakit.breeding.fit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'genotypes': samples x markers dosages, 'phenotypes': [...]}` (samples with a null phenotype are dropped). `options`: `ridge_lambda` (penalty on the centred dosages, > 0), or `h2` (0 to 1, then lambda = 2*sum(p*q)*(1-h2)/h2); one of the two is required.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `model` (`rrblup_ridge`), `ridge_lambda` and `ridge_lambda_source`, `intercept`, `raw_intercept`, `marker_effects`, `allele_frequencies`, `n_markers_used`, `n_training_samples`, `training_predictions`, `training_breeding_values`, `training_correlation`, `sample_ids` and `definition`. The warning notes that the training correlation overstates accuracy. The result can be passed as `model` to `predict`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import fit

genotypes = [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 1, 0], [1, 0, 2, 1], [2, 1, 0, 2]]
phenotypes = [1.0, 2.1, 2.9, 0.8, 3.2, 4.0]
result = fit({"genotypes": genotypes, "phenotypes": phenotypes}, options={"ridge_lambda": 1.0})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "breeding",
  "operation": "fit",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 145,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eef533a99c44baec24b7857c098f6dfa38f9dd4eb4bcdf24fec13da360ff651f"
  },
  "warnings": [
    "training_correlation is computed on the training data and overstates prediction accuracy"
  ],
  "status": "ok",
  "result": {
    "model": "rrblup_ridge",
    "ridge_lambda": 1.0,
    "ridge_lambda_source": "ridge_lambda",
    "intercept": 2.3333333333333335,
    "raw_intercept": 1.2766513056835636,
    "marker_effects": [
      0.6798771121351764,
      -0.26359447004608283,
      0.12319508448540714,
      0.47327188940092196
    ],
    "allele_frequencies": [
      0.5,
      0.4166666666666667,
      0.5,
      0.5
    ],
    "n_markers_used": 4,
    "n_training_samples": 6,
    "training_predictions": [
      1.259447004608295,
      2.6394777265745013,
      3.2328725038402455,
      0.8726574500768052,
      2.6761904761904765,
      3.3193548387096774
    ],
    "training_breeding_values": [
      -1.0738863287250384,
      0.3061443932411677,
      0.8995391705069121,
      -1.4606758832565283,
      0.34285714285714286,
      0.9860215053763441
    ],
    "training_correlation": 0.9290475470728303,
    "sample_ids": [
      "sample1",
      "sample2",
      "sample3",
      "sample4",
      "sample5",
      "sample6"
    ],
    "definition": "ridge regression of the phenotype on centred dosages (SNP-BLUP), objective sum((y-mu-Zb)^2) + lambda*sum(b^2); equivalen..."
  }
}
```

## 3) optimize

- **Purpose:** Optimise parental contributions (optimum-contribution selection): maximise the mean breeding value minus a coancestry penalty.

- **Method:** Solves the convex quadratic programme max c'u - (penalty/2) c'Ac subject to sum(c) = 1 and c >= 0, where u holds the candidates' breeding values and A is the relationship matrix, by an active-set method; the KKT conditions are checked at the solution. Group coancestry is reported as c'Ac/2. The relationship matrix is either given or computed from a pedigree with Henderson's tabular method.

- **API:** `dnakit.breeding.optimize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'candidates': records with id and ebv, 'relationship': {'ids': [...], 'matrix': [[...]]}}` or `{'candidates': ..., 'pedigree': records id/sire/dam}`. `options`: `penalty` (default 1.0), `trait` (candidate column holding the breeding value, default `ebv`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `contributions` ({id: contribution}, non-zero ones), `mean_ebv`, `mean_relationship` (c'Ac), `group_coancestry` (c'Ac/2), `penalty`, `objective` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import optimize

candidates = [{"id": "A", "ebv": 3.0}, {"id": "B", "ebv": 2.5}, {"id": "C", "ebv": 1.0}]
relationship = {"ids": ["A", "B", "C"], "matrix": [[1.0, 0.8, 0.0], [0.8, 1.0, 0.0], [0.0, 0.0, 1.0]]}
result = optimize({"candidates": candidates, "relationship": relationship}, options={"penalty": 10.0})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "breeding",
  "operation": "optimize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 196,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ecac1283d9753f119385ad279d1a07411bc14b7fc9f605e69b801c825d0af4bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "contributions": {
      "A": 0.43421052631578966,
      "B": 0.18421052631578932,
      "C": 0.3815789473684211
    },
    "mean_ebv": 2.1447368421052633,
    "mean_relationship": 0.49605263157894747,
    "group_coancestry": 0.24802631578947373,
    "penalty": 10.0,
    "objective": -0.3355263157894739,
    "definition": "maximise c'u - (penalty/2) c'Ac subject to sum(c) = 1, c >= 0 (convex QP solved by an active-set method; KKT conditions verified at the solution)"
  }
}
```

## 4) plot

- **Purpose:** Draw breeding progress (a value per generation) as an SVG line chart.

- **Method:** A fixed 640 x 400 SVG with a bordered plot area, rounded axis ticks, one polyline per series (colours cycle through a five-colour palette), a series label and the x axis labelled by generation index. The y range is the min and max over all series; a constant series gets a unit range. Pure standard library; nothing is written unless `path` is given.

- **API:** `dnakit.breeding.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'series': {name: [values per generation]}}`, or a `breeding.simulate` result (a mapping with `generations`), where `options['metric']` (default `mean_genetic_value`) selects the plotted column. `options`: `path` (write the SVG there; an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `svg` (the SVG text), `series` ({name: number of points}) and `path` (null unless written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import plot

result = plot({"series": {"mean_genetic_value": [0.0, 0.8, 1.5, 2.1]}})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "breeding",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a3ac612626547dbf535b0080025c610c92cd77554db08194a08ac88dfaef3fd7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"400\" viewBox=\"0 0 640 400\"><title>breeding progress</title><...",
    "series": {
      "mean_genetic_value": 4
    },
    "path": null
  }
}
```

## 5) predict

- **Purpose:** Predict genomic estimated breeding values of new individuals with a ridge/GBLUP model.

- **Method:** With a fitted `model` (the `result` of `fit`), prediction = raw_intercept + sum_j b_j * dosage_j and the breeding value is the prediction minus the model intercept; a missing dosage is replaced by 2 x the training allele frequency. Without a model, the same ridge model as `fit` is fitted on `train_genotypes` and `phenotypes` (numpy needed; `not_implemented` without it) and a missing dosage is replaced by the training mean.

- **API:** `dnakit.breeding.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'model': <result of breeding.fit>, 'genotypes': new dosages}` (same markers in the same order) or `{'train_genotypes': ..., 'phenotypes': [...], 'genotypes': new dosages}`. `options` (refit path only): `ridge_lambda` or `h2`, as in `fit`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `sample_ids`, `predictions` and `breeding_values`; the refit path also reports `ridge_lambda`, `ridge_lambda_source` and `n_markers_used`, the model path a `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import predict

result = predict({
    "train_genotypes": [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 1, 0], [1, 0, 2, 1], [2, 1, 0, 2]],
    "phenotypes": [1.0, 2.1, 2.9, 0.8, 3.2, 4.0],
    "genotypes": [[1, 1, 1, 1], [2, 0, 2, 0]],
}, options={"ridge_lambda": 1.0})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "breeding",
  "operation": "predict",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 194,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b87bce67e2f2c7d47c89d18b80153af40dae06f8663c05579f7b6b60bb332f59"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sample_ids": [
      "sample1",
      "sample2"
    ],
    "predictions": [
      2.2894009216589866,
      2.8827956989247308
    ],
    "breeding_values": [
      -0.04393241167434713,
      0.5494623655913973
    ],
    "ridge_lambda": 1.0,
    "ridge_lambda_source": "ridge_lambda",
    "n_markers_used": 4
  }
}
```

## 6) select

- **Purpose:** Select breeding candidates by truncation on a trait or EBV, or by a Smith-Hazel selection index.

- **Method:** `truncation` keeps the k candidates with the best trait value. `index` computes the Smith-Hazel index I = b'x with b = P^-1 G a (P phenotypic and G genetic covariance matrices, a economic weights), with accuracy r_IH = sqrt(b'Pb / a'Ga); a singular P is rejected. The selection differential is the mean of the selected values minus the mean of all candidates. Standard library only.

- **API:** `dnakit.breeding.select(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: records with `id` (unique) and trait columns. `options`: `method` = `truncation` (default) or `index`; `trait` (truncation column, default `ebv`); for `index`: `traits`, `economic_weights`, `genetic_covariance`, `phenotypic_covariance`; exactly one of `n_select` or `proportion`; `higher_is_better` (default True).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `method`, `n_select`, `selected` (ids), `ranking` (id, score, rank for every candidate), `selection_differential`, `trait` (truncation) or `index_weights`, `traits`, `index_variance` and `accuracy` (index), and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import select

records = [{"id": "a", "ebv": 1.2}, {"id": "b", "ebv": 3.4}, {"id": "c", "ebv": 2.2}, {"id": "d", "ebv": 0.4}]
result = select(records, options={"method": "truncation", "trait": "ebv", "n_select": 2})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "breeding",
  "operation": "select",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "92fceebe3688fb1f2ee7418794195107fa73df082d1b93eace96de08af4e055f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "truncation",
    "n_select": 2,
    "selected": [
      "b",
      "c"
    ],
    "ranking": [
      {
        "id": "b",
        "score": 3.4,
        "rank": 1
      },
      {
        "id": "c",
        "score": 2.2,
        "rank": 2
      },
      "..."
    ],
    "selection_differential": 0.9999999999999998,
    "trait": "ebv",
    "definition": "truncation: keep the k best by the trait value; index: Smith-Hazel index I = b'x with b = P^-1 G a, accuracy r_IH = sqrt(b'Pb / a'Ga)"
  }
}
```

## 7) simulate

- **Purpose:** Simulate recurrent truncation selection on a polygenic trait (unlinked additive loci) and report the response per generation.

- **Method:** Biallelic additive loci with base frequency 0.5 and effects ~ N(0,1); phenotype = genetic value + N(0, Ve) with Ve fixed from the base-population genetic variance so that the base heritability equals `heritability`. Each generation the top `selection_proportion` of phenotypes are kept, parents are drawn at random (no selfing) and gametes recombine freely. For every generation it reports the mean genetic value, variances, selection differential, the breeder's-equation prediction (Vg/Vp) * S for the parent generation, and the observed response (change in mean genetic value). Deterministic given `seed` (Python `random`); this is a small model, not a linked-genome forward simulator such as AlphaSimR. Limits: 10 to 50,000 individuals, up to 2,000 loci, up to 100 generations.

- **API:** `dnakit.breeding.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: None or a mapping of the options below (it overrides `options`). `options`: `n_individuals` (1000), `n_loci` (100), `heritability` (0.5, base population), `selection_proportion` (0.2), `generations` (5), `seed` (0).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `generations` (generations + 1 rows: `generation`, `mean_genetic_value`, `genetic_variance`, `phenotypic_variance`, `realised_heritability_in_generation`, `selection_differential`, `predicted_response`, and from generation 1 `observed_response`), `parameters` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.breeding import simulate

result = simulate(None, options={"n_individuals": 200, "n_loci": 50, "generations": 3, "seed": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "breeding",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "generations": [
      {
        "generation": 0,
        "mean_genetic_value": 5.112519418461102,
        "genetic_variance": 18.871869146667038,
        "phenotypic_variance": 37.95812172495276,
        "realised_heritability_in_generation": 0.49717605321501257,
        "selection_differential": 8.570761111862272,
        "predicted_response": 4.261177182644397
      },
      {
        "generation": 1,
        "mean_genetic_value": 8.822494073920149,
        "genetic_variance": 18.800333738727726,
        "phenotypic_variance": 43.15957790675255,
        "realised_heritability_in_generation": 0.43560050052728416,
        "selection_differential": 9.147702786858458,
        "predicted_response": 3.9847439126303765,
        "observed_response": 3.709974655459047
      },
      "..."
    ],
    "parameters": {
      "n_individuals": 200,
      "n_loci": 50,
      "heritability": 0.5,
      "selection_proportion": 0.2,
      "environmental_variance": 18.871869146667038,
      "seed": 1
    },
    "definition": "unlinked additive biallelic loci (base frequency 0.5, effects ~ N(0,1)); phenotype = genetic value + N(0, Ve) with Ve fi..."
  }
}
```
