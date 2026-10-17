# Experimental design

This category contains **9** documented capabilities, all computed natively in DNAKit with the standard library (none returns `not_implemented` in a normal call): sample-size, power, coverage and replicate calculations, randomised layouts (CRD, RCBD, Latin square, factorial), standard-order factorial matrices, run-order randomisation, treatment assignment, design checks and validation, CSV/TSV/JSON export and a Monte Carlo power simulation. Randomisation is reproducible through `seed`. The functions are exported through the DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) assign

- **Purpose:** Randomly assign experimental units to treatment groups by (stratified) permuted-block randomisation.

- **Method:** Within each stratum (all units form one stratum without `stratify_by`), the units are shuffled and cut into consecutive blocks of sum(ratio) units; each block receives a random permutation of the group pattern (each group name repeated `ratio` times). Group sizes therefore differ from the target ratio by less than one block per stratum. Randomisation uses Python `random.Random(seed)` and is reproducible.

- **API:** `dnakit.experimental_design.assign(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of unit ids, an integer number of units (ids `u1`..`un`), or records (dicts, optionally with a stratification column). `options`: `groups` (list of at least two unique names; default `G1..Gk` with `n_groups`=2), `ratio` (positive integers, one per group; default equal), `stratify_by` (column name), `seed` (0).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `seed`, `n_units`, `assignments` (each unit record with an added `group`), `group_sizes` and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import assign

units = [
    {"unit": "m1", "sex": "F"}, {"unit": "m2", "sex": "F"}, {"unit": "m3", "sex": "F"},
    {"unit": "m4", "sex": "F"}, {"unit": "m5", "sex": "M"}, {"unit": "m6", "sex": "M"},
    {"unit": "m7", "sex": "M"}, {"unit": "m8", "sex": "M"},
]
result = assign(units, options={"groups": ["control", "treated"], "stratify_by": "sex", "seed": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "experimental_design",
  "operation": "assign",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 8,
    "text_length": 224,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "004f3d304d047d6db323ca275724203cedd8d41f637060e2a047d422cc41db24"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 3,
    "n_units": 8,
    "assignments": [
      {
        "unit": "m1",
        "sex": "F",
        "group": "treated"
      },
      {
        "unit": "m2",
        "sex": "F",
        "group": "treated"
      },
      "..."
    ],
    "group_sizes": {
      "control": 4,
      "treated": 4
    },
    "method": "permuted-block randomisation: within each stratum, consecutive blocks of sum(ratio) units receive a random permutation o..."
  }
}
```

## 2) calculate

- **Purpose:** Calculate sample size, power, sequencing coverage (Lander-Waterman) or the number of replicates for common designs.

- **Method:** `quantity` selects the calculation. `sample_size` / `power` for `two_sample_t`, `one_sample_t` and `paired_t` use the exact non-central t distribution (delta = d*sqrt(n1*n2/(n1+n2)) or d*sqrt(n); effect size is Cohen's d); `two_proportions` uses the normal approximation with pooled null variance as in R `power.prop.test`; `chi_square` uses a non-central chi-square with lambda = n*w^2 (Cohen's w) and the central critical value at 1-alpha. The sample size is the smallest n reaching the target power (`larger` and `smaller` are both treated as one-sided). `coverage` follows Lander & Waterman: C = L*N/G, covered fraction 1 - exp(-C), expected contigs N*exp(-C*(1-theta)), expected contig length L*((exp(C(1-theta)) - 1)/C + 1 - theta), and the fraction covered at depth >= k from Poisson(C). `replicates` is the smallest n whose t confidence half-width t * sd / sqrt(n) is <= `margin` (the z-based n is reported too).

- **API:** `dnakit.experimental_design.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Parameters come from the `input` mapping and `options` (`input` wins on a name clash). `quantity` = `sample_size` (default), `power`, `coverage` or `replicates`. Tests: `test` = `two_sample_t` (default), `one_sample_t`, `paired_t` (`effect_size` = Cohen's d), `two_proportions` (`p1`, `p2`) or `chi_square` (`effect_size` = Cohen's w, `df`); `alpha` (0.05), `power` (0.8), `alternative` = `two-sided` (default) / `larger` / `smaller`, `ratio` (n2/n1, 1), `n` (required to compute power). Coverage: `genome_size`, `read_length`, `n_reads` or `coverage`, optional `target_fraction`, `min_depth` (1), `overlap_fraction` (0). Replicates: `sd`, `margin`, `confidence` (0.95).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` depends on the quantity. Sample size: `quantity`, `test`, `alpha`, `alternative`, `target_power`, `n_per_group`, `n`, `n_group2`, `total` (t tests) or `n_unrounded` (proportions) or `n_total` (chi-square), `achieved_power` and `method`. Power: `n` and `power`. Coverage: `coverage`, `fraction_covered`, `fraction_covered_min_depth`, `expected_uncovered_bases`, `expected_contigs`, `expected_contig_length` and, with `target_fraction`, `coverage_needed_for_target_fraction` and `reads_needed_for_target_fraction`. Replicates: `replicates`, `replicates_normal_approximation`, `sd`, `margin`, `confidence`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import calculate

result = calculate(
    {"effect_size": 0.5, "alpha": 0.05, "power": 0.8},
    options={"quantity": "sample_size", "test": "two_sample_t"},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "experimental_design",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 49,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3228c308633c7dde297f035b87166f4969c8f503a7f49fa9c7d59a78fbd49aca"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "sample_size",
    "test": "two_sample_t",
    "alpha": 0.05,
    "alternative": "two-sided",
    "target_power": 0.8,
    "n_per_group": 64,
    "n": 64,
    "n_group2": 64,
    "total": 128,
    "achieved_power": 0.8014595579222591,
    "method": "exact non-central t power; delta = d*sqrt(n1*n2/(n1+n2)) (two-sample) or d*sqrt(n) (one-sample/paired), df = n1+n2-2 or n-1"
  }
}
```

## 3) check

- **Purpose:** Check the statistical properties of a design table: replication balance, complete blocks, the Latin-square property and orthogonality of coded factor columns.

- **Method:** Counts the replications per treatment (balanced when all are equal; a treatment replicated fewer than twice is flagged). With a `block` column it checks that every block holds every treatment equally often; with `row` and `column` columns it checks that each treatment appears once per row and per column; for `factor_columns` it computes all pairwise dot products of the numeric columns and calls the design orthogonal when all are zero (within 1e-9). Every failed check is added to `issues`; `ok` is true when there are none.

- **API:** `dnakit.experimental_design.check(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: design records, the whole envelope returned by `design()` (a mapping with `status` and `result`; its `result` is used and a `status` other than `ok` raises `InputFormatError` with code `INVALID_DESIGN_ENVELOPE`), or the inner `result` mapping of `design()` (the one holding `layout`), with a `treatment` column (`options['treatment_column']`), optionally `block`, or `row` and `column`. `options`: `treatment_column` (`treatment`), `factor_columns` (list of numeric, e.g. -1/+1 coded, columns).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `n_plots`, `n_treatments`, `replications`, `balanced`, then `n_blocks` and `complete_blocks` (with a `block` column), `latin_square` (with `row` and `column`), `column_dot_products` and `orthogonal` (with `factor_columns`), `issues` (list of messages) and `ok`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import check, design

layout = design({"design_type": "rcbd", "treatments": ["A", "B", "C"], "replicates": 3, "seed": 1})
result = check(layout)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "experimental_design",
  "operation": "check",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 82,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50d77680785df2e143b723bc829268992fcc49786029d1d325a320e054846b10"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_plots": 9,
    "n_treatments": 3,
    "replications": {
      "B": 3,
      "C": 3,
      "A": 3
    },
    "balanced": true,
    "n_blocks": 3,
    "complete_blocks": true,
    "issues": [],
    "ok": true
  }
}
```

## 4) design

- **Purpose:** Build a randomised experimental layout: completely randomised, randomised complete block, Latin square or factorial.

- **Method:** `crd`: each treatment repeated `replicates` times and the units shuffled onto plots. `rcbd`: every block holds every treatment once in a freshly shuffled order (number of blocks = `blocks`, default `replicates`). `latin_square`: a cyclic Latin square whose rows, columns and treatment labels are randomly permuted. `factorial`: the treatments are all combinations of the factor levels (labels such as `A=1|B=x`), laid out as a CRD or, when `blocks` is given, as an RCBD, with one column per factor. Shuffling uses Python `random.Random(seed)`; layouts of more than 100,000 plots are refused.

- **API:** `dnakit.experimental_design.design(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping (merged with `options`). `design_type` = `crd` (default), `rcbd`, `latin_square` or `factorial`; `treatments` (list) or `n_treatments`; `replicates` (2; the number of replicates for CRD, of blocks for RCBD unless `blocks` is given); `factors` ({name: levels}) for factorial; `seed` (0).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `design_type`, `seed`, `n_plots`, `layout` (records with `plot`, `block` or `row`/`column`, `treatment`, and factor columns for factorial designs), `replications` ({treatment: count}) and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import design

result = design({"design_type": "rcbd", "treatments": ["A", "B", "C"], "replicates": 2, "seed": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "experimental_design",
  "operation": "design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 82,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "37fea3b5d78613b96f698d7b9dbf24c5dd634e12891148cf0671e29b399d67a1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "design_type": "rcbd",
    "seed": 1,
    "n_plots": 6,
    "layout": [
      {
        "plot": 1,
        "block": 1,
        "treatment": "B"
      },
      {
        "plot": 2,
        "block": 1,
        "treatment": "C"
      },
      "..."
    ],
    "replications": {
      "B": 2,
      "C": 2,
      "A": 2
    },
    "definition": "each block holds every treatment once in random order (randomised complete block)"
  }
}
```

## 5) export

- **Purpose:** Export a design layout as CSV, TSV or JSON text, and optionally write it to a file.

- **Method:** Serialises the design records with the chosen columns (all columns of the first record by default); CSV/TSV cells that are missing are written empty, JSON keeps `null`. A file is written only when `path` is given, and an existing file is not overwritten unless `overwrite=True`. The text is always returned.

- **API:** `dnakit.experimental_design.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: design records, the whole envelope returned by `design()` (a mapping with `status` and `result`; its `result` is used and a `status` other than `ok` raises `InputFormatError` with code `INVALID_DESIGN_ENVELOPE`), or the inner `result` mapping of `design()` (the one holding `layout`). `options`: `format` = `csv` (default), `tsv` or `json`; `columns` (subset and order); `path` (or `output`); `overwrite` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `format`, `n_rows`, `columns`, `text` (the serialised table) and `path` (null unless a file was written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import design, export

layout = design({"design_type": "crd", "treatments": ["A", "B"], "replicates": 2, "seed": 1})
result = export(layout, options={"format": "csv", "columns": ["plot", "treatment"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "experimental_design",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d8dbc7378da676ba95efe489f195176d930d39dd9ba2373fd6e2742835265bd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "csv",
    "n_rows": 4,
    "columns": [
      "plot",
      "treatment"
    ],
    "text": "plot,treatment\n1,B\n2,A\n3,B\n4,A\n",
    "path": null
  }
}
```

## 6) generate

- **Purpose:** Generate a standard-order design matrix: full factorial, two-level fractional factorial or cyclic Latin square.

- **Method:** `full_factorial`: all level combinations in standard (Yates) order, the first factor varying fastest (at most 100,000 runs). `fractional_factorial`: a 2^(k-p) design for `n_factors` (2 to 12, named A, B, ...) where each generator such as `D=ABC` defines a generated column as the product of base columns, coded -1/+1; the base factors are run in standard order. `latin_square`: the cyclic square t[(r+c) mod k] without randomisation (use `design` for a randomised one). The matrices are not randomised.

- **API:** `dnakit.experimental_design.generate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping (merged with `options`). `kind` = `full_factorial` (default; `factors` = {name: levels}), `fractional_factorial` (`n_factors` and `generators`, e.g. `['D=ABC']`, expressed in base factors only) or `latin_square` (`treatments`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `kind` and, for the factorial kinds, `n_runs`, `runs` (records with `run` and one column per factor), `order` and (fractional) `generators`; for `latin_square` it holds `square` (list of rows) and `order`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import generate

result = generate(
    {"kind": "fractional_factorial", "n_factors": 4, "generators": ["D=ABC"]}
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "experimental_design",
  "operation": "generate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 73,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "feb488336a79b02b2e5ed7e00157348fe2840f9536dcde128d34618844093248"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "fractional_factorial",
    "n_runs": 8,
    "runs": [
      {
        "run": 1,
        "A": -1,
        "B": -1,
        "C": -1,
        "D": -1
      },
      {
        "run": 2,
        "A": 1,
        "B": -1,
        "C": -1,
        "D": 1
      },
      "..."
    ],
    "generators": [
      "D=ABC"
    ],
    "order": "standard order of the base factors (first varies fastest); generated columns are products of base columns"
  }
}
```

## 7) randomize

- **Purpose:** Randomise the run order of a list of runs, optionally within groups, and add a 1-based `run_order`.

- **Method:** A Fisher-Yates shuffle with Python `random.Random(seed)`. With `within`, the runs are shuffled inside each group of that column and the groups keep their first-appearance order. Between 1 and 1,000,000 runs are accepted.

- **API:** `dnakit.experimental_design.randomize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of items (each becomes `{'item': x}`) or records. `options`: `seed` (0), `within` (column name).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `seed`, `n_runs`, `within`, `runs` (the records in the new order, each with `run_order`) and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import randomize

runs = [{"sample": "s1", "batch": "b1"}, {"sample": "s2", "batch": "b1"},
        {"sample": "s3", "batch": "b2"}, {"sample": "s4", "batch": "b2"}]
result = randomize(runs, options={"seed": 7, "within": "batch"})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "experimental_design",
  "operation": "randomize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 132,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "455cf2af29029a45016dc21756ea96ad764597784664f0e65a9305c082fe8872"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "seed": 7,
    "n_runs": 4,
    "within": "batch",
    "runs": [
      {
        "sample": "s1",
        "batch": "b1",
        "run_order": 1
      },
      {
        "sample": "s2",
        "batch": "b1",
        "run_order": 2
      },
      "..."
    ],
    "method": "Fisher-Yates shuffle (Python random.Random(seed)); with `within`, runs are shuffled inside each group and groups keep first-appearance order"
  }
}
```

## 8) simulate

- **Purpose:** Estimate statistical power by Monte Carlo simulation of a t-test experiment and compare it with the exact power.

- **Method:** Draws `n_simulations` data sets of normal data with unit SD (group 1 mean 0, group 2 mean = `effect_size`; or one sample with mean `effect_size`), applies a two-sided pooled-variance t-test (one-sample t-test for `one_sample_t`) and reports the fraction with p < alpha as the empirical power, with its binomial standard error, next to the exact non-central-t power. Deterministic given `seed` (Python `random`).

- **API:** `dnakit.experimental_design.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping (merged with `options`) with `n` (per group), `effect_size` (Cohen's d), `test` = `two_sample_t` (default) or `one_sample_t`, `alpha` (0.05), `n_simulations` (2000), `seed` (0).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `test`, `n`, `effect_size`, `alpha`, `n_simulations`, `seed`, `empirical_power`, `standard_error`, `analytic_power` and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import simulate

result = simulate({"n": 20, "effect_size": 0.8, "n_simulations": 500, "seed": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "experimental_design",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 62,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dd66b43353d0d0b3ff488fa6d81c75710a1925dd884f65c2e3d613f7f6e6beca"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "two_sample_t",
    "n": 20,
    "effect_size": 0.8,
    "alpha": 0.05,
    "n_simulations": 500,
    "seed": 1,
    "empirical_power": 0.704,
    "standard_error": 0.020414896521902825,
    "analytic_power": 0.6934041966153844,
    "method": "Monte Carlo: normal data with unit SD, two-sided pooled t-test; compare with the exact non-central-t power"
  }
}
```

## 9) validate

- **Purpose:** Validate the structure of a design table: required columns present, no empty cells, unique plot ids and consistent rows.

- **Method:** Checks that the table has rows, that each required column exists and has no empty values, that the unit column (if present) has no duplicates and that all rows have the same number of columns. Content problems are collected in `errors` and give `valid: false`; they do not raise.

- **API:** `dnakit.experimental_design.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: design records, the whole envelope returned by `design()` (a mapping with `status` and `result`; its `result` is used and a `status` other than `ok` raises `InputFormatError` with code `INVALID_DESIGN_ENVELOPE`), or the inner `result` mapping of `design()` (the one holding `layout`). `options`: `required_columns` (default `['treatment']`), `unit_column` (default `plot`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `valid`, `n_rows` and `errors` (list of messages). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.experimental_design import validate

rows = [{"plot": 1, "treatment": "A"}, {"plot": 2, "treatment": ""}, {"plot": 2, "treatment": "B"}]
result = validate(rows)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "experimental_design",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 92,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ad0afe16f1236e53b11be4318dc7efee802a144ac7090c559777a3a886c87510"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_rows": 3,
    "errors": [
      "column 'treatment' has 1 empty value(s)",
      "unit column 'plot' has duplicated values"
    ]
  }
}
```
