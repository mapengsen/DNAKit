# Regulation

This category contains **11** documented candidate capabilities; 10 compute real results natively (PWM motif scanning and scoring, motif comparison and validation, enrichment tests, variant effect proxies, in-silico mutagenesis, position-based region classification and distance-based enhancer-gene linking) and 1 (`predict`) always returns `status: "not_implemented"` because it needs trained sequence-model weights. All scoring is explicitly heuristic: motif scores are PWM log-odds (sum of log2(p / background)) with no hit p-values and no trained model, `link`/`annotate(entity="target_gene")` assign regions to genes by distance to the TSS only, and `classify` uses a position rule only. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Annotate sequences with PWM motif hits, or regions with their nearest target gene (`options["entity"]`).

- **Method:** `motif` (default; alias `tf_binding_site`): every window of every sequence is scored with each PWM. PWM scoring is the standard log-odds model: a matrix of counts or probabilities (columns summing to 1 within 1e-3 are taken as probabilities, otherwise counts; `matrix_type` can force either) is turned into probabilities p (counts: (count + `pseudocount`) / (column total + 4 `pseudocount`), `pseudocount` default 0; probabilities are floored at 1e-6, so a zero entry costs about -17.9 bits against a 0.25 background), and the score of a window is the sum over positions of log2(p(base) / background(base)) with `background` default 0.25 for each base. The relative score is (score - min) / (max - min) of the matrix. Both strands are scanned (the reverse strand scores the reverse-complemented window); windows containing N are skipped. This is a simple affinity proxy: there is no hit p-value, no genome-wide background model and no trained model. A window is a hit when its score is >= `threshold` (absolute bits) or, if no `threshold` is given, when its relative score is >= `min_relative_score` (default 0.8). Hits are sorted by position and capped at `max_hits` per sequence (`truncated` tells whether the cap was reached). `target_gene`: identical to `link` with `method="nearest"`. The distance is |TSS - midpoint| with midpoint = (start + end) // 2 of the region, and the TSS is the explicit `tss` if given, otherwise `start` for + or unstranded genes and `end - 1` for - genes. `method="nearest"` links each region to the single nearest TSS on the same chromosome (ties go to the alphabetically first gene name) if it lies within `max_distance`; `method="window"` links every gene whose TSS is within `max_distance`. `direction` is the position of the region relative to the gene orientation (`upstream`, `downstream` or `at_tss`). This is a proximity heuristic, not evidence of regulation: no contact, expression or TAD information is used, and the call adds a warning to this effect.

- **API:** `dnakit.regulation.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `motif`: `input` is DNA (string, FASTA text/path, list or `{id: sequence}`; ACGT and N, U is read as T); `options`: `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`); `threshold` (absolute bits) or `min_relative_score` (default 0.8); `strand` (`"both"` default | `"forward"` | `"reverse"`); `max_hits` (default 10000). `target_gene`: `input` is regions as BED text/path, rows or mappings (0-based half-open); `options`: `genes` (required: rows or mappings with `chrom`, `start`, `end`, `name` and optional `strand` or explicit `tss`), `max_distance` (default 1,000,000 bp). `options["entity"]` is `"motif"` (default), `"tf_binding_site"` or `"target_gene"`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. For `motif`, `result` has `entity` and `sequences` (per sequence `id`, `length`, `hit_count`, `hits`, `truncated`; each hit has `motif`, 0-based half-open forward `start`/`end`, `strand`, `score` (bits), `relative_score` and `matched`, the motif-oriented matched sequence). For `target_gene` the result is that of `link`: `result` has `method`, `max_distance` and `links` (one entry per region: `region`, `midpoint` and `genes`, each with `gene`, `tss`, `distance`, `direction`, `strand`; the list is empty when no gene is in range).

- **Example code:**

```python
import json
from dnakit.regulation import annotate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
result = annotate("TTGATAAGGCTTATCAC",
                  options={"motifs": {"GATA": gata}, "pseudocount": 0.5, "min_relative_score": 0.9})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 17,
    "sequence_bases": 17,
    "lengths": [
      17
    ],
    "sha256": "ef32207549d8a6455d7b256ce50fb87992ed624845b1ae6279f906be3d504032"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "motif",
    "sequences": [
      {
        "id": "seq1",
        "length": 17,
        "hit_count": 2,
        "hits": [
          {
            "motif": "GATA",
            "start": 2,
            "end": 7,
            "strand": "+",
            "score": 8.897784945974175,
            "relative_score": 1.0,
            "matched": "GATAA"
          },
          {
            "motif": "GATA",
            "start": 10,
            "end": 15,
            "strand": "-",
            "score": 8.897784945974175,
            "relative_score": 1.0,
            "matched": "GATAA"
          }
        ],
        "truncated": false
      }
    ]
  }
}
```

## 2) calculate

- **Purpose:** Calculate PWM-based regulatory effect deltas of variants, or the information content of motifs (`options["entity"]`).

- **Method:** `regulatory_effect_delta` (default): for each motif the best log2-odds score over all windows and both strands is computed for the reference and the alternate sequence, and `delta` = alt_score - ref_score in bits (negative means weaker predicted binding); `max_abs_delta` is the largest |delta| over the motifs. PWM scoring is the standard log-odds model: a matrix of counts or probabilities (columns summing to 1 within 1e-3 are taken as probabilities, otherwise counts; `matrix_type` can force either) is turned into probabilities p (counts: (count + `pseudocount`) / (column total + 4 `pseudocount`), `pseudocount` default 0; probabilities are floored at 1e-6, so a zero entry costs about -17.9 bits against a 0.25 background), and the score of a window is the sum over positions of log2(p(base) / background(base)) with `background` default 0.25 for each base. The relative score is (score - min) / (max - min) of the matrix. Both strands are scanned (the reverse strand scores the reverse-complemented window); windows containing N are skipped. This is a simple affinity proxy: there is no hit p-value, no genome-wide background model and no trained model. It is a position-weight-matrix proxy and not a trained model. `pwm_information_content`: per-position information content sum over bases of p log2(p / background) and its total in bits (2 bits for a fully conserved column with a uniform background).

- **API:** `dnakit.regulation.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `regulatory_effect_delta`: `input` is variants as a record or list, either `{"id", "ref_sequence", "alt_sequence"}` or `{"id", "sequence", "position", "ref", "alt"}` (0-based `position`; `ref` must match the sequence, otherwise an error is raised), at most 10000; `options`: `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`). `pwm_information_content`: `input` is one PWM or several (the `motifs` forms above); `options`: `background`, `pseudocount`, `matrix_type`. `options["entity"]` is `"regulatory_effect_delta"` (default) or `"pwm_information_content"`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `regulatory_effect_delta`: `result` has `entity` and `variants` (per variant `id`, `motifs` with `motif`, `ref_score`, `alt_score`, `delta`, and `max_abs_delta`; motifs without any scorable window are omitted), and `warnings` states that PWM log-odds is a simple proxy. `pwm_information_content`: `result` has `entity` and `motifs` (`motif`, `length`, `per_position`, `total_bits`).

- **Example code:**

```python
import json
from dnakit.regulation import calculate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
variant = {"id": "v1", "sequence": "TTGATAAGG", "position": 4, "ref": "T", "alt": "C"}
result = calculate(variant, options={"entity": "regulatory_effect_delta", "motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 9,
    "lengths": [
      9
    ],
    "sha256": "d79c2457c84e2ff160d36cf926e7da44d70d976b7e8eee511dbfc530fac2bd4f"
  },
  "warnings": [
    "PWM log-odds is a simple binding-affinity proxy."
  ],
  "status": "ok",
  "result": {
    "entity": "regulatory_effect_delta",
    "variants": [
      {
        "id": "v1",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 3.6123827271119264,
            "delta": -5.285402218862249
          }
        ],
        "max_abs_delta": 5.285402218862249
      }
    ]
  }
}
```

## 3) classify

- **Purpose:** Classify regions as promoter, proximal or distal by their position relative to gene TSSs (position rule only).

- **Method:** The TSS is the explicit `tss` if given, otherwise `start` for + or unstranded genes and `end - 1` for - genes. The promoter window is [TSS - `promoter_upstream`, TSS + `promoter_downstream`) for + genes and [TSS - `promoter_downstream` + 1, TSS + `promoter_upstream` + 1) for - genes. A region overlapping any promoter window is `promoter`; otherwise it is `proximal` when its midpoint (start + end) // 2 is within `proximal_distance` of a TSS on the same chromosome, else `distal` (also when the chromosome has no genes). The nearest gene is reported by midpoint-to-TSS distance. The thresholds are conventional defaults, and no chromatin, conservation or expression evidence is used (the call adds a warning).

- **API:** `dnakit.regulation.classify(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: regions as BED text/path, rows or mappings (0-based half-open). `options`: `genes` (required: rows or mappings with `chrom`, `start`, `end`, `name` and optional `strand` or explicit `tss`); `promoter_upstream` (default 1000), `promoter_downstream` (default 100), `proximal_distance` (default 10000); `entity` (only `"promoter"`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.regions` has one entry per region: `chrom`, `start`, `end`, `class` (`promoter`, `proximal` or `distal`), `nearest_gene` and `distance_to_tss` (both `null` when the chromosome has no genes).

- **Example code:**

```python
import json
from dnakit.regulation import classify

genes = [{"name": "GENE_A", "chrom": "chr1", "start": 5000, "end": 9000, "strand": "+"}]
regions = "chr1\t4600\t4800\tpromoter_like\nchr1\t9500\t9700\tnear\nchr1\t80000\t80400\tfar\n"
result = classify(regions, options={"genes": genes})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "classify",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ae3771f589ac85420ee7b88632e87b7bebaffb31e60bde9d3553caf7ca010ee8"
  },
  "warnings": [
    "Classification uses genomic position only."
  ],
  "status": "ok",
  "result": {
    "regions": [
      {
        "chrom": "chr1",
        "start": 4600,
        "end": 4800,
        "class": "promoter",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 300
      },
      {
        "chrom": "chr1",
        "start": 9500,
        "end": 9700,
        "class": "proximal",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 4600
      },
      {
        "chrom": "chr1",
        "start": 80000,
        "end": 80400,
        "class": "distal",
        "nearest_gene": "GENE_A",
        "distance_to_tss": 75200
      }
    ]
  }
}
```

## 4) compare

- **Purpose:** Compare PWM motifs pairwise by their best ungapped alignment (mean per-column Pearson correlation).

- **Method:** Motifs are converted to probability matrices. For every pair a, b all offsets and both orientations of b (forward and reverse complement) with at least `min_overlap` aligned columns are tried; the similarity of an alignment is the mean over the aligned columns of the Pearson correlation between the A/C/G/T probability vectors (a column with zero variance counts as 0). The best alignment gives `similarity` (1.0 = identical), `offset` (position of b's first column relative to a's), `orientation` and `overlap`. This is a simple column-correlation score; it is not Tomtom's p-value/E-value and gives no significance.

- **API:** `dnakit.regulation.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: two or more motifs (a list of matrices or `{name: matrix}`; matrices as counts or probabilities in the forms of `annotate`). `options`: `entity` (only `"motif"`); `min_overlap` (default 3, reduced to the shortest motif length); `background`, `pseudocount`, `matrix_type` as in `annotate`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.pairs` has one entry per unordered pair: `a`, `b`, `similarity`, `offset`, `orientation` (`forward` or `reverse_complement`) and `overlap`.

- **Example code:**

```python
import json
from dnakit.regulation import compare

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
ttatc = [[0, 0, 0, 10], [0, 0, 0, 10], [10, 0, 0, 0], [0, 0, 0, 10], [0, 10, 0, 0]]   # rows A, C, G, T per position
result = compare({"GATA": gata, "TTATC": ttatc}, options={"min_overlap": 4})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 189,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8983b08a91514321197f33dac0b3688b1e0ef8afdc615b50ecc871d94d24680e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": [
      {
        "a": "GATA",
        "b": "TTATC",
        "similarity": 0.9989409928538651,
        "offset": 0,
        "orientation": "reverse_complement",
        "overlap": 5
      }
    ]
  }
}
```

## 5) enrich

- **Purpose:** Test motif enrichment in foreground versus background sequences (hypergeometric or binomial upper tail).

- **Method:** The unit is the number of sequences with at least one motif hit. `hypergeometric` (default): with N = fg_total + bg_total, K = fg_hits + bg_hits and n = fg_total, p = P(X >= fg_hits) for X ~ Hypergeometric(N, K, n). `binomial`: p = P(X >= fg_hits) for X ~ Binomial(fg_total, bg_hits / bg_total). `fold_enrichment` = (fg_hits / fg_total) / (bg_hits / bg_total) (`null` when the background rate is 0), and `q_value` is the Benjamini-Hochberg adjustment over the tested motifs. With sequences, a hit is decided by the PWM scan of `annotate` (`threshold` or `min_relative_score`), so the result depends on that cut-off, and the background is whatever sequences are supplied (no matching for GC or length is done).

- **API:** `dnakit.regulation.enrich(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping, either counts `{"foreground_hits", "foreground_total", "background_hits", "background_total"}` (non-negative integers; hits cannot exceed totals, totals must be positive) or sequences `{"foreground": [...], "background": [...]}`. `options`: `entity` (`"motif"` default | `"tf_binding_site"`); `test` (`"hypergeometric"` default | `"binomial"`); `name` (label for the counts form, default `motif`); with sequences also `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`), `threshold`/`min_relative_score`, `strand`; each motif is tested.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `test` and `motifs` (per motif `motif`, `foreground_hits`, `foreground_total`, `background_hits`, `background_total`, `fold_enrichment`, `p_value`, `q_value`).

- **Example code:**

```python
import json
from dnakit.regulation import enrich

counts = {"foreground_hits": 30, "foreground_total": 100, "background_hits": 10, "background_total": 100}
result = enrich(counts, options={"test": "hypergeometric", "name": "GATA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "enrich",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 96,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "74bb3c6e8fb50effdc77642b2707471638887481358940928a82775ea5670cf4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "test": "hypergeometric",
    "motifs": [
      {
        "motif": "GATA",
        "foreground_hits": 30,
        "foreground_total": 100,
        "background_hits": 10,
        "background_total": 100,
        "fold_enrichment": 2.9999999999999996,
        "p_value": 0.000325205353801686,
        "q_value": 0.000325205353801686
      }
    ]
  }
}
```

## 6) link

- **Purpose:** Assign regulatory regions (enhancers and so on) to genes by the distance from the region midpoint to the gene TSS.

- **Method:** Distance-based assignment. The distance is |TSS - midpoint| with midpoint = (start + end) // 2 of the region, and the TSS is the explicit `tss` if given, otherwise `start` for + or unstranded genes and `end - 1` for - genes. `method="nearest"` links each region to the single nearest TSS on the same chromosome (ties go to the alphabetically first gene name) if it lies within `max_distance`; `method="window"` links every gene whose TSS is within `max_distance`. `direction` is the position of the region relative to the gene orientation (`upstream`, `downstream` or `at_tss`). This is a proximity heuristic, not evidence of regulation: no contact, expression or TAD information is used, and the call adds a warning to this effect. The `entity` value (`enhancer`, `promoter`, `silencer`, `insulator`, `target_gene`) is only a label: the rule is the same for all of them.

- **API:** `dnakit.regulation.link(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: regions as BED text/path, rows or mappings (0-based half-open). `options`: `entity` (`"enhancer"` default | `"promoter"` | `"silencer"` | `"insulator"` | `"target_gene"`); `genes` (required: rows or mappings with `chrom`, `start`, `end`, `name` and optional `strand` or explicit `tss`); `method` (`"nearest"` default | `"window"`); `max_distance` (default 1,000,000 bp).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `method`, `max_distance` and `links` (one entry per region: `region`, `midpoint` and `genes`, each with `gene`, `tss`, `distance`, `direction`, `strand`; the list is empty when no gene is in range).

- **Example code:**

```python
import json
from dnakit.regulation import link

enhancers = "chr1\t10000\t10400\tenh1\nchr1\t90000\t90200\tenh2\n"
genes = [
    {"name": "GENE_A", "chrom": "chr1", "start": 12000, "end": 20000, "strand": "+"},
    {"name": "GENE_B", "chrom": "chr1", "start": 60000, "end": 91000, "strand": "-"},
]
result = link(enhancers, options={"entity": "enhancer", "genes": genes, "method": "nearest", "max_distance": 50000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "link",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b35ce0cf833e4ae26ca3cce75f823f7c9a74d1755641cf0d0f0eb456361e590c"
  },
  "warnings": [
    "Distance-based assignment is a proximity heuristic, not evidence of regulation."
  ],
  "status": "ok",
  "result": {
    "method": "nearest",
    "max_distance": 50000,
    "links": [
      {
        "region": {
          "chrom": "chr1",
          "start": 10000,
          "end": 10400,
          "name": "enh1"
        },
        "midpoint": 10200,
        "genes": [
          {
            "gene": "GENE_A",
            "tss": 12000,
            "distance": 1800,
            "direction": "upstream",
            "strand": "+"
          }
        ]
      },
      {
        "region": {
          "chrom": "chr1",
          "start": 90000,
          "end": 90200,
          "name": "enh2"
        },
        "midpoint": 90100,
        "genes": [
          {
            "gene": "GENE_B",
            "tss": 90999,
            "distance": 899,
            "direction": "downstream",
            "strand": "-"
          }
        ]
      }
    ]
  }
}
```

## 7) perform_in_silico_mutagenesis

- **Purpose:** In-silico saturation mutagenesis of a sequence scored by a PWM (every single-base substitution).

- **Method:** For every position i and every alternative base b the effect is best(mutated) - best(reference), where best is the maximum PWM log2-odds score over all windows and both strands (the log-odds scoring of `score`: sum over positions of log2(p(base) / background(base)), background 0.25 by default). The model is the PWM itself, not a trained neural network, and a substitution has effect 0 when it neither lowers the best-scoring window nor creates a better one. A sequence shorter than a motif is an error.

- **API:** `dnakit.regulation.perform_in_silico_mutagenesis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one DNA sequence (string, FASTA text/path or `{"sequence": ...}`), at most 2000 bases, ACGT only. `options`: `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `sequence` (id) and `motifs`; per motif `motif`, `reference_score` and `positions`, a list with one entry per position: `position`, `ref` and `effects` (a mapping from each alternative base to the change in bits). `warnings` states that the model is a PWM.

- **Example code:**

```python
import json
from dnakit.regulation import perform_in_silico_mutagenesis

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
result = perform_in_silico_mutagenesis("AGATAAC", options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "regulation",
  "operation": "perform_in_silico_mutagenesis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 7,
    "sequence_bases": 7,
    "lengths": [
      7
    ],
    "sha256": "06ce638d4a9efc08740e01618439ea14f498bb9df12324428d7885607a97afcd"
  },
  "warnings": [
    "The scoring model is a PWM, not a trained sequence model."
  ],
  "status": "ok",
  "result": {
    "sequence": "seq1",
    "motifs": [
      {
        "motif": "GATA",
        "reference_score": 8.897784945974175,
        "positions": [
          {
            "position": 0,
            "ref": "A",
            "effects": {
              "C": 0.0,
              "G": 0.0,
              "T": 0.0
            }
          },
          {
            "position": 1,
            "ref": "G",
            "effects": {
              "A": -3.5443205162238103,
              "C": -3.5443205162238103,
              "T": -3.5443205162238103
            }
          },
          {
            "position": 2,
            "ref": "A",
            "effects": {
              "C": -5.285402218862249,
              "G": -5.285402218862249,
              "T": -3.7004397181410926
            }
          },
          {
            "position": 3,
            "ref": "T",
            "effects": {
              "A": -3.7004397181410926,
              "C": -5.285402218862249,
              "G": -5.285402218862249
            }
          },
          "..."
        ]
      }
    ]
  }
}
```

## 8) predict

- **Purpose:** Predict regulatory activity from sequence. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Sequence-to-function prediction needs a trained model such as Enformer, Borzoi, DeepSEA or Sei together with its weights, which DNAKit does not bundle. For a PWM-based proxy use `score`, `calculate` or `perform_in_silico_mutagenesis`. Nothing is fabricated.

- **API:** `dnakit.regulation.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a trained regulatory sequence model and weights (Enformer, Borzoi, DeepSEA, Sei)` and `message` is `Regulatory-activity prediction needs trained model weights that are not bundled.`

- **Example code:**

```python
import json
from dnakit.regulation import predict

result = predict("ACGTACGTACGTACGT")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 16,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "cf573e65038d08ff910a3345642ffd1e8329844633c2dcb15964b324ebdba4d0"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained regulatory sequence model and weights (Enformer, Borzoi, DeepSEA, Sei)",
  "message": "Regulatory-activity prediction needs trained model weights that are not bundled."
}
```

## 9) prioritize

- **Purpose:** Rank variants by the largest absolute PWM score change across motifs (a proxy for regulatory impact).

- **Method:** Each variant is scored as in `calculate` (`regulatory_effect_delta`): the best log2-odds score over all windows and both strands for the reference and alternate sequence, with `delta` = alt_score - ref_score in bits. Variants are sorted by `max_abs_delta` descending, ties by id; variants without any scorable motif window go last. The ranking is by this PWM proxy only; conservation, chromatin state and trained models are not used (the call adds a warning).

- **API:** `dnakit.regulation.prioritize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: variants as a record or list, either `{"id", "ref_sequence", "alt_sequence"}` or `{"id", "sequence", "position", "ref", "alt"}` (0-based `position`; `ref` must match the sequence, otherwise an error is raised), at most 10000. `options`: `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`); `entity` (only `"regulatory_variant"`); `top_n` (optional).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.ranked` has the variant entries of `calculate` (`id`, `motifs`, `max_abs_delta`) plus a 1-based `rank`.

- **Example code:**

```python
import json
from dnakit.regulation import prioritize

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
variants = [
    {"id": "in_motif", "sequence": "TTGATAAGG", "position": 4, "ref": "T", "alt": "C"},
    {"id": "outside", "sequence": "TTGATAAGG", "position": 8, "ref": "G", "alt": "A"},
]
result = prioritize(variants, options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "prioritize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 167,
    "sequence_bases": 18,
    "lengths": [
      9,
      9
    ],
    "sha256": "e31572e30c31d4dff3508cf3b45e845f40532b2bb237769b2d2ac7f3b605f6f8"
  },
  "warnings": [
    "Ranking uses a PWM score change proxy only."
  ],
  "status": "ok",
  "result": {
    "ranked": [
      {
        "id": "in_motif",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 3.6123827271119264,
            "delta": -5.285402218862249
          }
        ],
        "max_abs_delta": 5.285402218862249,
        "rank": 1
      },
      {
        "id": "outside",
        "motifs": [
          {
            "motif": "GATA",
            "ref_score": 8.897784945974175,
            "alt_score": 8.897784945974175,
            "delta": 0.0
          }
        ],
        "max_abs_delta": 0.0,
        "rank": 2
      }
    ]
  }
}
```

## 10) score

- **Purpose:** Score sequences against PWM motifs: the best log2-odds window per sequence and motif.

- **Method:** Every window on the allowed strands is scored and the highest-scoring one is reported with its position. PWM scoring is the standard log-odds model: a matrix of counts or probabilities (columns summing to 1 within 1e-3 are taken as probabilities, otherwise counts; `matrix_type` can force either) is turned into probabilities p (counts: (count + `pseudocount`) / (column total + 4 `pseudocount`), `pseudocount` default 0; probabilities are floored at 1e-6, so a zero entry costs about -17.9 bits against a 0.25 background), and the score of a window is the sum over positions of log2(p(base) / background(base)) with `background` default 0.25 for each base. The relative score is (score - min) / (max - min) of the matrix. Both strands are scanned (the reverse strand scores the reverse-complemented window); windows containing N are skipped. This is a simple affinity proxy: there is no hit p-value, no genome-wide background model and no trained model. The best score is the maximum over windows, so it grows with sequence length and is not a probability or a significance.

- **API:** `dnakit.regulation.score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: DNA (string, FASTA text/path, list or `{id: sequence}`; ACGT and N, U is read as T). `options`: `options["motifs"]` (required; `motif` is accepted as an alias): one PWM or several (a list or `{name: pwm}`), each given as counts or probabilities in the forms `{"A": [...], "C": [...], "G": [...], "T": [...]}`, rows `[[a, c, g, t], ...]` or JASPAR-style text (at most 1000 positions and 500 motifs); `background` (mapping A/C/G/T, default 0.25 each), `pseudocount` (default 0), `matrix_type` (`"auto"` default | `"counts"` | `"probabilities"`); `entity` (`"motif"` default | `"tf_binding_site"`); `strand` (`"both"` default | `"forward"` | `"reverse"`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.sequences` has one entry per sequence: `id` and `motifs`, each with `motif` and `best` (`start`, `end` as 0-based half-open forward coordinates, `strand`, `score` in bits, `relative_score`), or `best: null` when the sequence is shorter than the motif or has only windows containing N.

- **Example code:**

```python
import json
from dnakit.regulation import score

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
sequences = {"s1": "CCCGATAAGGCC", "s2": "CCCCCCCCCCCC"}
result = score(sequences, options={"motifs": {"GATA": gata}, "pseudocount": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "db0e3c9cf6a6c083e8cbd0b701de8166870199076e7dbe5f73b61e223ee0f361"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequences": [
      {
        "id": "s1",
        "motifs": [
          {
            "motif": "GATA",
            "best": {
              "start": 3,
              "end": 8,
              "strand": "+",
              "score": 8.897784945974175,
              "relative_score": 1.0
            }
          }
        ]
      },
      {
        "id": "s2",
        "motifs": [
          {
            "motif": "GATA",
            "best": {
              "start": 0,
              "end": 5,
              "strand": "-",
              "score": -12.16787507624152,
              "relative_score": 0.1440196391333088
            }
          }
        ]
      }
    ]
  }
}
```

## 11) validate

- **Purpose:** Validate PWM motifs: shape, finite non-negative entries, non-zero columns and whether the columns sum to 1.

- **Method:** Each matrix is parsed (A/C/G/T lists, an Lx4 list of rows or JASPAR-style text). It is invalid when it cannot be parsed, has negative or non-finite entries or a column that sums to zero (the error text is returned in `errors`). `matrix_type` is `probabilities` when every column sums to 1 within `tolerance`, otherwise `counts`. For valid matrices `total_information_bits` is the PWM information content sum over positions of sum_b p log2(p / 0.25) with the probability floor of 1e-6. A formatting and sanity check only: it does not judge whether a motif is biologically meaningful.

- **API:** `dnakit.regulation.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one PWM or several (a list or `{name: matrix}`) as counts or probabilities in the forms accepted by `annotate`. `options`: `entity` (only `"motif"`); `tolerance` (default 1e-3) for the column-sum test; `background`, `pseudocount` as in `annotate`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"regulation"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `motifs` (per motif `motif`, `valid`, `errors`, `length`, and for valid matrices `matrix_type` and `total_information_bits`) and `all_valid`. Malformed matrices are reported with `valid: false` rather than raising an error.

- **Example code:**

```python
import json
from dnakit.regulation import validate

gata = {"A": [1, 19, 1, 19, 18], "C": [1, 0, 0, 0, 1], "G": [17, 0, 0, 0, 0], "T": [1, 1, 19, 1, 1]}   # counts, GATAA-like
bad = [[1, 0, 0, 0], [0, -1, 0, 0]]
result = validate({"GATA": gata, "bad": bad})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "regulation",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 141,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a8ead67ae3934d9b11ccd1049efcc1a88aac7adb7a3ccf0433059478c62c867e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "motifs": [
      {
        "motif": "GATA",
        "valid": true,
        "errors": [],
        "length": 5,
        "matrix_type": "counts",
        "total_information_bits": 7.7240815230348225
      },
      {
        "motif": "bad",
        "valid": false,
        "errors": [
          "bad: negative or non-finite entries at positions [1]."
        ],
        "length": 2
      }
    ],
    "all_valid": false
  }
}
```
