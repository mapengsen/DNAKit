# Population genetics

This category contains **8** documented capabilities (the optional dadi/tskit backend wrappers that `dnakit.population` also exports are not described on this page). **7** compute results natively: `calculate` (allele frequency, heterozygosity, pairwise LD), `compare` (Fst), `estimate` (GRM, kinship, inbreeding, PCA), `export`, `plot`, `simulate` (Balding-Nichols and Wright-Fisher) and `test` (Hardy-Weinberg, differentiation); only `estimate` (grm, kinship, pca) and `plot(kind='pca')` need numpy and return `status: "not_implemented"` without it. **1**, `infer` (ancestry and structure inference), returns `status: "not_implemented"` in every call because it needs a reference panel or a model-fitting program. The functions are exported through the DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) calculate

- **Purpose:** Calculate per-variant allele frequencies, heterozygosity or pairwise linkage disequilibrium from a genotype matrix.

- **Method:** `allele_frequency`: alternate-allele frequency = sum(dosage)/(2 x called samples), MAF = min(p, 1-p), number called and missing rate, overall and (with `populations`) per population; fractional imputed dosages are allowed. `heterozygosity` (hard calls only): Ho = fraction of called samples with dosage 1, He = 2pq from the sample frequency, the unbiased He = 2pq x 2n/(2n-1) and F = 1 - Ho/He. `ld`: composite (genotype-dosage) LD for all variant pairs, r^2 = squared Pearson correlation of dosages, D = cov(dosage_a, dosage_b)/2 and D' = |D|/Dmax with Dmax = min(pA qB, qA pB) if D > 0 else min(pA pB, qA qB); pairs with fewer than 3 complete samples or a monomorphic variant are skipped. Pure standard library.

- **API:** `dnakit.population.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). `options`: `statistic` = `allele_frequency` (default), `heterozygosity` or `ld`; `max_variants` (100, at most 300; for `ld`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `statistic` and, for `allele_frequency`, `overall` (per variant: `variant_id`, `n_called`, `alt_allele_frequency`, `maf`, `missing_rate`) and `by_population`; for `heterozygosity`, `variants` (`n_called`, `alt_allele_frequency`, `observed_heterozygosity`, `expected_heterozygosity`, `expected_heterozygosity_unbiased`, `f_statistic`) and `definition`; for `ld`, `pairs` (`variant_1`, `variant_2`, `n`, `r2`, `d_composite`, `d_prime`) and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import calculate

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = calculate(genotypes, options={"statistic": "allele_frequency"})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "population",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "allele_frequency",
    "overall": [
      {
        "variant_id": "v1",
        "n_called": 6,
        "alt_allele_frequency": 0.5,
        "maf": 0.5,
        "missing_rate": 0.0
      },
      {
        "variant_id": "v2",
        "n_called": 6,
        "alt_allele_frequency": 0.5833333333333334,
        "maf": 0.41666666666666663,
        "missing_rate": 0.0
      },
      "..."
    ],
    "by_population": {
      "A": [
        {
          "variant_id": "v1",
          "n_called": 3,
          "alt_allele_frequency": 0.16666666666666666,
          "maf": 0.16666666666666666,
          "missing_rate": 0.0
        },
        {
          "variant_id": "v2",
          "n_called": 3,
          "alt_allele_frequency": 0.3333333333333333,
          "maf": 0.3333333333333333,
          "missing_rate": 0.0
        },
        "..."
      ],
      "B": [
        {
          "variant_id": "v1",
          "n_called": 3,
          "alt_allele_frequency": 0.8333333333333334,
          "maf": 0.16666666666666663,
          "missing_rate": 0.0
        },
        {
          "variant_id": "v2",
          "n_called": 3,
          "alt_allele_frequency": 0.8333333333333334,
          "maf": 0.16666666666666663,
          "missing_rate": 0.0
        },
        "..."
      ]
    }
  }
}
```

## 2) compare

- **Purpose:** Compare populations by Fst (Weir & Cockerham 1984 or Hudson), globally, per variant and pairwise.

- **Method:** `weir_cockerham` (default): the diploid variance components a, b, c of Weir & Cockerham (1984) per biallelic variant, theta = sum(a)/sum(a+b+c) over variants (ratio of sums), computed for all populations jointly and for every population pair. `hudson`: the pairwise estimator of Bhatia et al. (2013), the ratio of sums over variants of [(p1-p2)^2 - p1q1/(n1-1) - p2q2/(n2-1)] / [p1 q2 + p2 q1] with haploid sample sizes n; the per-variant table is given only when there are exactly two populations. Variants with undefined estimates are null and excluded from the ratio. Needs hard calls (0/1/2) and at least two populations.

- **API:** `dnakit.population.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). A `populations` list (one label per sample) is required. `options`: `method` = `weir_cockerham` (default) or `hudson`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `method`, `populations` ({name: sample count}), `fst` and `per_variant` ({variant_id: Fst}; for Hudson only with two populations), `pairwise` ({`popA|popB`: Fst}) and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import compare

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = compare(genotypes, options={"method": "weir_cockerham"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "population",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "weir_cockerham",
    "populations": {
      "A": 3,
      "B": 3
    },
    "fst": 0.459016393442623,
    "per_variant": {
      "v1": 0.5384615384615384,
      "v2": 0.31818181818181834,
      "v3": 0.5384615384615384,
      "v4": 0.5384615384615384,
      "v5": 0.3181818181818181
    },
    "definition": "Weir & Cockerham (1984) theta = sum(a) / sum(a + b + c) over variants (ratio of sums), diploid genotype-based components",
    "pairwise": {
      "A|B": 0.459016393442623
    }
  }
}
```

## 3) estimate

- **Purpose:** Estimate genomic relationships (VanRaden GRM or kinship), per-sample inbreeding F, or genotype principal components.

- **Method:** `grm` (default): VanRaden method 1 genomic relationship matrix. `kinship`: the same matrix divided by 2 (coancestry). `pca`: PCA of frequency-standardised dosages X = (dosage - 2p)/sqrt(2p(1-p)) by thin SVD, scores U*S, covariance eigenvalues S^2/(n-1); entirely missing or invariant markers are dropped and missing dosages are filled with the marker mean. These three use the numpy-based helpers of `dnakit.genomics` and return `not_implemented` without numpy. `inbreeding`: the PLINK `--het` method-of-moments F = (O_hom - E_hom)/(N - E_hom) per sample with E_hom = sum(1 - 2pq x 2n/(2n-1)) over called markers (hard calls, standard library).

- **API:** `dnakit.population.estimate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). `options`: `quantity` = `grm` (default), `kinship`, `inbreeding` or `pca`; `components` (2; for `pca`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `quantity` and, for `grm` and `kinship`, `sample_ids`, `relationship_matrix`, `method`, `n_markers_used`, `allele_frequencies` and marker-filtering details; for `pca`, `scores`, `loadings`, `eigenvalues`, `explained_variance_ratio`, `sample_ids`, `parameters` and marker details; for `inbreeding`, `samples` (`sample_id`, `n_markers`, `observed_homozygous`, `expected_homozygous`, `f`) and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import estimate

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = estimate(genotypes, options={"quantity": "grm"})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "population",
  "operation": "estimate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "quantity": "grm",
    "sample_ids": [
      "s1",
      "s2",
      "s3",
      "s4",
      "s5",
      "s6"
    ],
    "method": "vanraden_method_1",
    "ploidy": 2,
    "relationship_matrix": [
      [
        1.2359550561797752,
        0.8988764044943819,
        0.49438202247191015,
        -1.2584269662921348,
        -0.9213483146067415,
        -0.449438202247191
      ],
      [
        0.8988764044943819,
        1.3707865168539324,
        0.15730337078651685,
        -1.1910112359550562,
        -0.8539325842696629,
        -0.38202247191011235
      ],
      "..."
    ],
    "n_samples": 6,
    "denominator": 2.4722222222222223,
    "imputed_value_count": 0,
    "allele_frequencies": [
      0.5,
      0.5833333333333334,
      0.5,
      0.5,
      0.4166666666666667
    ],
    "marker_indices": [
      0,
      1,
      2,
      3,
      4
    ],
    "dropped_all_missing_indices": [],
    "dropped_monomorphic_indices": [],
    "n_markers_used": 5,
    "citations": [
      "https://doi.org/10.3168/jds.2007-0980",
      "https://doi.org/10.1534/g3.112.004259"
    ]
  }
}
```

## 4) export

- **Purpose:** Export a genotype matrix as TSV/CSV, PLINK text PED+MAP, or an allele-frequency table.

- **Method:** `tsv`/`csv`: one row per sample with `sample_id` and one column per variant. `ped_map` (hard calls): PLINK text files, genotype cells `ref ref`, `ref alt`, `alt alt` or `0 0` for missing, with the population label as family id, and a MAP with chromosome, variant id, 0 and position; `ref`/`alt`/`chrom`/`pos` lists in the input name the alleles (otherwise `A`/`G` and chromosome 0). `frq`: columns `SNP`, `A1`, `A2`, `MAF`, `NCHROBS` (2 x called samples), with `ALT`/`REF` as allele names when the input has none. The text is always returned; files are written only when `path` is given and existing files are not overwritten unless `overwrite=True`.

- **API:** `dnakit.population.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). `options`: `format` = `tsv` (default), `csv`, `ped_map` or `frq`; `path` (for `ped_map` the file stem: `<path>.ped` and `<path>.map` are written); `overwrite` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `format`, `n_samples`, `n_variants`, `text` (tsv, csv, frq) or `ped` and `map` (ped_map), and `path` (only when a file was written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import export

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = export(genotypes, options={"format": "frq"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "population",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "frq",
    "n_samples": 6,
    "n_variants": 5,
    "text": "SNP\tA1\tA2\tMAF\tNCHROBS\nv1\tALT\tREF\t0.5\t12\nv2\tALT\tREF\t0.41666666666666663\t12\nv3\tALT\tREF\t0.5\t12\nv4\tALT\tREF\t0.5\t12\nv5\tALT\tREF\t0.4166666666666667\t12\n"
  }
}
```

## 5) infer

- **Purpose:** Infer genetic ancestry or population structure. Not implemented in DNAKit.

- **Method:** No computation is done. Ancestry or structure inference needs labelled reference genotypes (for example 1000 Genomes) or a model-fitting program of the ADMIXTURE/STRUCTURE type, neither of which is bundled. The message points to `estimate(quantity='pca')` for a PCA of your own samples.

- **API:** `dnakit.population.infer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any `input` and optional `options`; they are only summarised in `input_summary`.

- **Output:** Does not compute a result. A normal call returns the envelope with `status: "not_implemented"` (`provenance.implementation` is `not_implemented`), plus `required_dependency` and `message` stating what is missing: an ancestry reference panel (e.g. 1000 Genomes) or an ADMIXTURE/STRUCTURE-type program.

- **Example code:**

```python
import json
from dnakit.population import infer

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = infer(genotypes)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "population",
  "operation": "infer",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "ancestry reference panel (e.g. 1000 Genomes) or an ADMIXTURE/STRUCTURE-type program",
  "message": "Ancestry/structure inference needs labelled reference genotypes or a model-fitting backend, neither of which is bundled; use estimate(quantity='pca') for a PCA of your own samples."
}
```

## 6) plot

- **Purpose:** Draw a population-genetics plot (MAF histogram, PCA scatter or per-variant Fst) as an SVG string.

- **Method:** `maf_histogram` (default): MAFs binned over [0, 0.5] into `bins` equal bins. `pca`: the first two principal components of `estimate(quantity='pca')`, points coloured by `populations` (numpy needed; `not_implemented` without it). `fst`: per-variant Fst from `compare` (`method` as there; needs `populations`) as a scatter over the variant index. Axes with rounded ticks are drawn with pure standard-library SVG; nothing is written unless `path` is given.

- **API:** `dnakit.population.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). `options`: `kind` = `maf_histogram` (default), `pca` or `fst`; `bins` (10); `method` (for `fst`, default `weir_cockerham`); `width` (640), `height` (420); `path` (an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `kind`, `svg`, `path` (null unless written) and, per kind, `bin_counts`, `bin_width`, `n_variants` (histogram), `n_points` (fst) or `n_points` and `populations` (pca). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import plot

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = plot(genotypes, options={"kind": "maf_histogram", "bins": 5})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "population",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "maf_histogram",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"420\" viewBox=\"0 0 640 420\"><title>population maf_histogram p...",
    "path": null,
    "bin_counts": [
      0,
      0,
      0,
      0,
      5
    ],
    "bin_width": 0.1,
    "n_variants": 5
  }
}
```

## 7) simulate

- **Purpose:** Simulate population genotypes (Balding-Nichols structured populations) or neutral Wright-Fisher allele-frequency drift.

- **Method:** `balding_nichols` (default): an ancestral frequency p per variant is drawn uniformly from [maf_min, 1-maf_min]; each population draws p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F) and each sample's genotype is Binomial(2, p_k). `wright_fisher`: for each replicate a biallelic locus with `initial_frequency` drifts for `generations` generations as p' = Binomial(2N, p)/2N (stopping at fixation or loss); the expected variance p0 q0 (1 - (1 - 1/2N)^t) is reported next to the observed one. Deterministic given `seed` (Python `random`). It is not a coalescent or forward simulator such as msprime or SLiM.

- **API:** `dnakit.population.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: None or a mapping of the options (it overrides `options`). `model` = `balding_nichols` (default): `n_populations` (2), `n_samples_per_population` (50), `n_variants` (100), `fst` (0.1), `maf_min` (0.1). `model` = `wright_fisher`: `population_size` (100 diploid), `generations` (50), `initial_frequency` (0.5), `replicates` (100). `seed` (0) for both.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. For `balding_nichols`, `result` holds `model`, `genotypes` (a genotype mapping usable by the other functions: `genotypes`, `populations`, `variant_ids`, `sample_ids`), `ancestral_frequencies`, `population_frequencies`, `fst`, `seed` and `definition`. For `wright_fisher`: `model`, `population_size`, `generations`, `initial_frequency`, `final_frequencies`, `mean_final_frequency`, `variance_final_frequency`, `expected_variance`, `trajectories` (the first 20 replicates), `seed` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import simulate

result = simulate(
    None,
    options={"model": "balding_nichols", "n_populations": 2, "n_samples_per_population": 3,
             "n_variants": 4, "fst": 0.1, "seed": 1},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "population",
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
    "model": "balding_nichols",
    "genotypes": {
      "genotypes": [
        [
          0,
          2,
          2,
          0
        ],
        [
          0,
          1,
          2,
          0
        ],
        "..."
      ],
      "populations": [
        "pop1",
        "pop1",
        "pop1",
        "pop2",
        "pop2",
        "pop2"
      ],
      "variant_ids": [
        "var1",
        "var2",
        "var3",
        "var4"
      ],
      "sample_ids": [
        "s1",
        "s2",
        "s3",
        "s4",
        "s5",
        "s6"
      ]
    },
    "ancestral_frequencies": [
      0.207491395289921,
      0.7779469895497861,
      0.7110196951812913,
      0.30405522059153733
    ],
    "population_frequencies": [
      [
        0.17901345124622997,
        0.8619574841214532,
        0.910923497312446,
        0.20385620033598056
      ],
      [
        0.39824481339193973,
        0.926723126702485,
        0.8023538517295546,
        0.3022080001595772
      ]
    ],
    "fst": 0.1,
    "seed": 1,
    "definition": "Balding-Nichols model: p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F) per population, genotypes ~ Binomial(2, p_k)"
  }
}
```

## 8) test

- **Purpose:** Test Hardy-Weinberg equilibrium per variant (exact or chi-square) or allele-frequency differentiation between populations.

- **Method:** `hwe` with `method='exact'` (default): the exact test of Wigginton, Cutler & Abecasis (2005), p = sum of the probabilities of heterozygote counts no more likely than the observed one. `method='chisq'`: Pearson chi-square (1 df) of the genotype counts against the expectations from the sample frequency (null when an expected count is 0). `differentiation`: Pearson chi-square test of homogeneity of the 2 x k allele-count table across populations (df = k - 1; null when a variant is monomorphic). Every tested variant also gets a Bonferroni-adjusted p (p x m) and a Benjamini-Hochberg q-value. Needs hard calls.

- **API:** `dnakit.population.test(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Genotypes are `{'genotypes': samples x variants alternate-allele dosages (0/1/2; null = missing), 'sample_ids': [...], 'variant_ids': [...], 'populations': [...]}` or just the matrix (samples are then named `sample1`.., variants `variant1`..). `populations` is required for `differentiation`. `options`: `test` = `hwe` (default) or `differentiation`; `method` = `exact` (default) or `chisq` (for `hwe`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `test`, `records` (per variant: `variant_id`, genotype counts `n_hom_ref`, `n_het`, `n_hom_alt` for HWE or `df` for differentiation, `chi2`, `p`, `p_bonferroni`, `q_bh`), `n_tested` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.population import test

genotypes = {
    "genotypes": [
        [0, 1, 2, 0, 1],
        [0, 1, 1, 0, 2],
        [1, 0, 2, 1, 1],
        [2, 1, 0, 2, 0],
        [2, 2, 0, 1, 0],
        [1, 2, 1, 2, 1],
    ],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "variant_ids": ["v1", "v2", "v3", "v4", "v5"],
    "populations": ["A", "A", "A", "B", "B", "B"],
}
result = test(genotypes, options={"test": "hwe", "method": "exact"})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "population",
  "operation": "test",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 263,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a04e58eeb71de4f8571c17160d4ba013e2cf61b5ada4c433ddea22fce45bdc9d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "hwe",
    "records": [
      {
        "variant_id": "v1",
        "n_hom_ref": 2,
        "n_het": 2,
        "n_hom_alt": 2,
        "chi2": null,
        "p": 0.4805194805194805,
        "p_bonferroni": 1.0,
        "q_bh": 0.8008658008658008
      },
      {
        "variant_id": "v2",
        "n_hom_ref": 1,
        "n_het": 3,
        "n_hom_alt": 2,
        "chi2": null,
        "p": 1.0,
        "p_bonferroni": 1.0,
        "q_bh": 1.0
      },
      "..."
    ],
    "n_tested": 5,
    "definition": "exact test of Wigginton et al. (2005): p = sum of probabilities of heterozygote counts no more likely than observed"
  }
}
```
