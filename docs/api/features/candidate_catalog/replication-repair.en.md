# Replication and DNA repair

This category contains **10** documented capabilities. All 10 compute a result natively in pure Python (no external program or optional library) and return `status: "ok"`; none returns `not_implemented`. The scientific models are simple and some are explicitly heuristic: GC-skew origin/terminus estimates, DnaA-box density, an indel-signature rule for the repair pathway and a deterministic fork model are not substitutes for experimental replication or repair data. The functions are exported through the formal DNAKit namespace and also exposed through the Agent/MCP registry where the JSON call plan is compatible. The example outputs below are real outputs; long arrays and strings are abbreviated with `...`.

## 1) annotate

- **Purpose:** Annotates genomic features as leading-strand or lagging-strand relative to the replication origin and terminus of a circular chromosome.

- **Method:** The chromosome is treated as circular with 0-based, half-open coordinates. The right replichore runs from the origin to the terminus in increasing coordinates and the left replichore from the terminus back to the origin. A feature (taken at its midpoint) is co-directional with the fork (leading strand) when it lies on the + strand in the right replichore or on the - strand in the left replichore, and head-on (lagging strand) otherwise. Origin and terminus come from `origin`/`terminus` in the input or options; if no origin is given it is estimated from the cumulative GC skew of `sequence` (a heuristic, see `calculate`), and the terminus defaults to the antipode of the origin. Pure Python; no external program.

- **API:** `dnakit.replication_repair.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping with `features` (list of `{id, start, end, strand}`; `strand` is `+` or `-`) and either `length` (chromosome length) or `sequence` (DNA string), plus optional `origin` and `terminus`. Options: `origin`, `terminus` (override the input; defaults: estimated from GC skew of `sequence`, and the antipode of the origin). `options` is optional.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `length`, `origin`, `terminus`, `origin_source`, the two replichore lengths, `features` (each with `replichore`, `orientation`, `leading_strand` and `distance_from_origin_along_replichore`), `leading_strand_features`, `lagging_strand_features`, `leading_strand_fraction` and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import annotate

chromosome = {
    "length": 1000,
    "origin": 0,
    "terminus": 500,
    "features": [
        {"id": "g1", "start": 100, "end": 200, "strand": "+"},
        {"id": "g2", "start": 300, "end": 400, "strand": "-"},
        {"id": "g3", "start": 700, "end": 800, "strand": "-"},
    ],
}
result = annotate(chromosome)
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 225,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "453b7bb62ffd20dd42f2b582bfd5ca17c525b4b2bd305999ddc7c9b05f6b1f77"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 1000,
    "origin": 0,
    "terminus": 500,
    "origin_source": "supplied",
    "right_replichore_length": 500,
    "left_replichore_length": 500,
    "features": [
      {
        "id": "g1",
        "start": 100,
        "end": 200,
        "strand": "+",
        "replichore": "right",
        "orientation": "co-directional (leading strand)",
        "leading_strand": true,
        "distance_from_origin_along_replichore": 150
      },
      {
        "id": "g2",
        "start": 300,
        "end": 400,
        "strand": "-",
        "replichore": "right",
        "orientation": "head-on (lagging strand)",
        "leading_strand": false,
        "distance_from_origin_along_replichore": 350
      },
      {
        "id": "g3",
        "start": 700,
        "end": 800,
        "strand": "-",
        "replichore": "left",
        "orientation": "co-directional (leading strand)",
        "leading_strand": true,
        "distance_from_origin_along_replichore": 250
      }
    ],
    "leading_strand_features": 2,
    "lagging_strand_features": 1,
    "leading_strand_fraction": 0.666667,
    "definition": "Circular chromosome; the right replichore runs origin -> terminus in increasing coordinates, the left replichore terminus -> origin. A feature (position = midpo..."
  }
}
```

## 2) build_recombination_map

- **Purpose:** Builds a genetic (recombination) map from ordered markers: centimorgan distances between markers and cM/Mb recombination rates per interval.

- **Method:** For each marker after the first, the distance to the previous marker is taken from a given `cm` (cumulative cM), from a recombination fraction `r`, or from `r = recombinants / meioses`. A recombination fraction is converted with the Haldane map function d = -1/2 ln(1 - 2r) (default), the Kosambi map function d = 1/4 ln((1 + 2r) / (1 - 2r)), or left as d = r (`none`); r must be below 0.5 for Haldane and Kosambi. The distance is multiplied by 100 to give cM. The rate of an interval is delta cM divided by the interval length in Mb; the result also gives cumulative cM, the total map length and the mean rate (total cM over span in Mb). Pure Python; no external program.

- **API:** `dnakit.replication_repair.build_recombination_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a list (or `{"markers": [...]}`) of at least two markers sorted by position, each `{id, position}` (bp, non-negative integer). Every marker after the first needs one of `cm` (cumulative cM), `recombination_fraction` (r to the previous marker) or `recombinants` together with `meioses`. The first marker may carry a starting `cm`. Options: `map_function` (`haldane` default, `kosambi`, `none`), applied only to r-based markers.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `map_function`, `markers`, `total_cm`, `span_bp`, `mean_rate_cm_per_mb`, `map` (one row per marker with `cumulative_cm` and, from the second marker on, `interval_start`, `interval_end`, `recombination_fraction`, `delta_cm` and `rate_cm_per_mb`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import build_recombination_map

markers = [
    {"id": "m1", "position": 0},
    {"id": "m2", "position": 1_000_000, "recombination_fraction": 0.1},
    {"id": "m3", "position": 3_000_000, "recombinants": 5, "meioses": 50},
]
result = build_recombination_map(markers, options={"map_function": "kosambi"})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "build_recombination_map",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 164,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d84d990f59e801eaea2a20335fc868493341e19a259092cf7d02ba067d4f388d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "map_function": "kosambi",
    "markers": 3,
    "total_cm": 20.273256,
    "span_bp": 3000000,
    "mean_rate_cm_per_mb": 6.757752,
    "map": [
      {
        "id": "m1",
        "position": 0,
        "cumulative_cm": 0.0
      },
      {
        "id": "m2",
        "position": 1000000,
        "interval_start": 0,
        "interval_end": 1000000,
        "recombination_fraction": 0.1,
        "delta_cm": 10.136628,
        "cumulative_cm": 10.136628,
        "rate_cm_per_mb": 10.136628
      },
      {
        "id": "m3",
        "position": 3000000,
        "interval_start": 1000000,
        "interval_end": 3000000,
        "recombination_fraction": 0.1,
        "delta_cm": 10.136628,
        "cumulative_cm": 20.273255,
        "rate_cm_per_mb": 5.068314
      }
    ],
    "definition": "Interval rate = delta cM / (interval length in Mb). Distances from recombination fractions r use the Haldane (d = -1/2 ln(1-2r)) or Kosambi (d = 1/4 ln((1+2r)/(..."
  }
}
```

## 3) calculate

- **Purpose:** Calculates windowed GC (or AT) skew and cumulative skew of a DNA sequence and a heuristic estimate of the replication origin and terminus.

- **Method:** Skew per window is (G - C) / (G + C) (`skew_type` `GC`) or (A - T) / (A + T) (`AT`); the cumulative skew is the running sum of the window skews. Heuristic (Lobry 1996, Grigoriev 1998): the origin is the centre of the window with the minimum cumulative skew and the terminus the centre of the window with the maximum. This holds only for bacterial-like chromosomes with a strand composition bias and has window-size resolution. Pure Python; no external program.

- **API:** `dnakit.replication_repair.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a DNA string, single-record FASTA text, a file path or `{"sequence": ...}` (up to 20 Mb). Options: `window` (default min(1000, length)), `step` (default = `window`), `skew_type` (`GC` default or `AT`), `include_windows` (default `True`), `max_windows_returned` (default 20000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `skew_type`, `length`, `window`, `step`, `windows`, `overall_skew`, `origin_estimate` with `origin_window`, `terminus_estimate` with `terminus_window`, `cumulative_skew_range`, `heuristic`, `definition` and `window_table` (per window `start`, `end`, `skew`, `cumulative_skew` and the base counts). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import calculate

sequence = "GGGC" * 125 + "CCCG" * 125
result = calculate(sequence, options={"window": 100})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 1000,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "3490f072d0775230a9eb8cb9e2bc13b5ea9c4ae2878fa5265873210e8ad3c480"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "skew_type": "GC",
    "length": 1000,
    "window": 100,
    "step": 100,
    "windows": 10,
    "overall_skew": 0.0,
    "origin_estimate": 950,
    "origin_window": [
      900,
      1000
    ],
    "terminus_estimate": 450,
    "terminus_window": [
      400,
      500
    ],
    "cumulative_skew_range": 2.5,
    "heuristic": "Origin = centre of the window with the minimum cumulative skew, terminus = maximum (Lobry 1996; Grigoriev 1998). Valid only for circular bacterial-like chromoso...",
    "definition": "skew = (G-C)/(G+C) per window; cumulative skew = running sum of the window skews.",
    "window_table": [
      {
        "start": 0,
        "end": 100,
        "skew": 0.5,
        "cumulative_skew": 0.5,
        "count_G": 75,
        "count_C": 25
      },
      {
        "start": 100,
        "end": 200,
        "skew": 0.5,
        "cumulative_skew": 1.0,
        "count_G": 75,
        "count_C": 25
      },
      "..."
    ]
  }
}
```

## 4) classify_repair_pathway

- **Purpose:** Classifies the likely DNA double-strand-break repair pathway (NHEJ, MMEJ, SSA or HDR) from the reference and the repaired sequence with an indel-signature rule.

- **Method:** This is a heuristic rule, not a measurement of the pathway actually used. If the edited sequence equals the supplied `donor` the call is HDR. A pure deletion whose two junctions share at least `ssa_min_homology` identical bp is SSA, at least `min_microhomology` bp is MMEJ, otherwise NHEJ. A pure insertion of up to `max_nhej_insertion` bp is NHEJ. Substitutions, complex indels and longer insertions are `undetermined`, and identical sequences give `no_edit`. Junction homology is how far the deletion can slide: prefix match + suffix match - edited length. Pure Python; no external program.

- **API:** `dnakit.replication_repair.classify_repair_pathway(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping `{"reference": unedited sequence, "edited": repaired sequence}` with an optional `donor` (the expected HDR product). Options: `min_microhomology` (default 2), `ssa_min_homology` (default 25), `max_nhej_insertion` (default 20).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `pathway` (`NHEJ`, `MMEJ`, `SSA`, `HDR`, `no_edit` or `undetermined`), `reason`, `rule` (the rule text) and the sequence-change fields `reference_length`, `edited_length`, `deleted_bases`, `inserted_bases`, `first_difference` and, for deletions, `microhomology_length`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import classify_repair_pathway

pair = {"reference": "GGATGCATGCCGG", "edited": "GGATGCCGG"}
result = classify_repair_pathway(pair)
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "classify_repair_pathway",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4ed8687da434749ff2110833c98a49bd59d81d81f13b2294402bfa1059c84aa9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reference_length": 13,
    "edited_length": 9,
    "deleted_bases": 4,
    "inserted_bases": 0,
    "first_difference": 6,
    "microhomology_length": 4,
    "pathway": "MMEJ",
    "reason": "deletion of 4 bp with 4 bp microhomology at the junction (>= 2): microhomology-mediated end joining signature.",
    "rule": "Heuristic indel-signature rule, not a measurement of the pathway used: edited == donor -> HDR; pure deletion whose two junctions share >= ssa_min_homology (25) ..."
  }
}
```

## 5) compare

- **Purpose:** Compares the cumulative GC (or AT) skew profiles and the origin estimates of two DNA sequences.

- **Method:** The cumulative skew of each sequence (as in `calculate`) is resampled to a common number of points (`points`) over relative position and the two profiles are compared by Pearson correlation: +1 means the same shape and -1 an inverted skew, for example for a complementary sequence. The relative origin positions (cumulative-skew minima, a heuristic) are compared as well. Pure Python; no external program.

- **API:** `dnakit.replication_repair.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is `{"a": sequence, "b": sequence}` or a list of exactly two DNA strings/FASTA texts. Options: `window` (default length/20, at most 1000), `skew_type` (`GC` default or `AT`), `points` (default 100).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `skew_type`, `sequences` (per sequence `length`, `windows`, origin and terminus estimates with their relative positions), `cumulative_skew_correlation`, `origin_relative_position_difference` and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import compare

seq_a = "GGGC" * 125 + "CCCG" * 125
seq_b = seq_a.translate(str.maketrans("ACGT", "TGCA"))
result = compare({"a": seq_a, "b": seq_b}, options={"window": 50, "points": 20})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 2018,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "426e3a41b1350557a532a1a49d3da3c8d80283e11497373df65cc2cdb90834ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "skew_type": "GC",
    "sequences": {
      "a": {
        "length": 1000,
        "windows": 20,
        "origin_estimate": 975,
        "origin_relative_position": 0.975,
        "terminus_estimate": 475,
        "terminus_relative_position": 0.475
      },
      "b": {
        "length": 1000,
        "windows": 20,
        "origin_estimate": 475,
        "origin_relative_position": 0.475,
        "terminus_estimate": 975,
        "terminus_relative_position": 0.975
      }
    },
    "cumulative_skew_correlation": -1.0,
    "origin_relative_position_difference": 0.5,
    "definition": "Cumulative GC skew of each sequence is resampled to 20 points over relative position and compared with Pearson correlation (+1 = identical shape, -1 = inverted ..."
  }
}
```

## 6) detect

- **Purpose:** Detects mismatch-repair-relevant motifs in a DNA sequence: GATC (Dam) sites and perfect mono- to tetranucleotide microsatellite tracts.

- **Method:** An exact string scan. GATC is the methylation site that E. coli MutHLS uses to tell the parental strand from the new one; it is palindromic, so each site is one position. Microsatellite tracts are perfect tandem repeats of 1 to 4 bp units at or above the minimum copy numbers; these are the loci that mismatch-repair-deficient cells expand or contract. The scan does not model instability, and the density of GATC sites is reported next to the expected random density. Pure Python; no external program.

- **API:** `dnakit.replication_repair.detect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a DNA string, single-record FASTA text, a file path or `{"sequence": ...}`. Options: `min_copies` (mapping from unit length 1-6 to the minimum number of copies; default `{1: 8, 2: 5, 3: 4, 4: 3}`), `max_results` (default 1000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `length`, `gatc_sites` (`count`, `density_per_kb`, `expected_random_per_kb`, `positions`, `note`), `microsatellite_tracts` (`count`, `min_copies`, `tracts` with `start`, `end`, `unit`, `copies`, `length`, and a `note`) and a `heuristic` description. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import detect

sequence = "AAGATCTT" + "A" * 10 + "G" + "CA" * 6 + "GATC" + "TTT"
result = detect(sequence)
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 38,
    "sequence_bases": 38,
    "lengths": [
      38
    ],
    "sha256": "b3adc1d06006d25e6bd4d7958addc957746ee5779f68906991301c079daef47d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 38,
    "gatc_sites": {
      "count": 2,
      "density_per_kb": 52.631579,
      "expected_random_per_kb": 3.9062,
      "positions": [
        2,
        31
      ],
      "note": "GATC is the Dam-methylation site used by E. coli MutHLS to tell the parental from the new strand; it is palindromic, so each site is one position on both strand..."
    },
    "microsatellite_tracts": {
      "count": 2,
      "min_copies": {
        "1": 8,
        "2": 5,
        "3": 4,
        "4": 3
      },
      "tracts": [
        {
          "start": 8,
          "end": 18,
          "unit": "A",
          "copies": 10,
          "length": 10
        },
        {
          "start": 19,
          "end": 31,
          "unit": "CA",
          "copies": 6,
          "length": 12
        }
      ],
      "note": "Perfect repeats of 1-4 bp units are the loci that mismatch-repair-deficient cells expand or contract (microsatellite instability); presence does not imply insta..."
    },
    "heuristic": "Exact string scan for the GATC motif and perfect mono- to tetranucleotide tandem repeats at or above the minimum copy numbers; no statistical model of instabili..."
  }
}
```

## 7) map

- **Purpose:** Maps replication timing along the genome as the median-normalised log2 ratio of S-phase to G1 read depth.

- **Method:** Per window, score = log2((S / median S) / (G1 / median G1)) (in the style of Koren et al. 2012); a higher score means earlier replication. `pseudocount` is added to each median-normalised depth before the ratio, and is needed if windows can have zero depth. An optional running median over `smoothing_windows` windows is applied to the scores. Windows above `center` are labelled `early`, those below `late` (the rest `intermediate`), and interior windows that are local maxima above `center` are reported as early-replication peaks. Pure Python; no external program.

- **API:** `dnakit.replication_repair.map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping `{"s_phase": [depth per window], "g1": [depth per window]}` over the same fixed-size windows. Options: `pseudocount` (default 0; must be > 0 if windows can be empty), `smoothing_windows` (odd integer, default 1 = none), `window_size` (bp, default 1, used for coordinates), `center` (default 0), `max_windows_returned` (default 20000).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `windows`, `pseudocount`, `smoothing_windows`, `early_fraction`, `late_fraction`, `early_peaks` (`window`, `start`, `log2_s_over_g1`), `window_table` (per window `start`, `end`, `log2_s_over_g1`, `timing`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import map

depths = {"s_phase": [20, 10, 5, 10, 40, 10], "g1": [10, 10, 10, 10, 10, 10]}
result = map(depths, options={"window_size": 1000})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "map",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 68,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "22da7963a7d6b1d21d72ddc8d686c2e9f526341676ebc8e99db9f63470a79320"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "windows": 6,
    "pseudocount": 0.0,
    "smoothing_windows": 1,
    "early_fraction": 0.333333,
    "late_fraction": 0.166667,
    "early_peaks": [
      {
        "window": 4,
        "start": 4000,
        "log2_s_over_g1": 2.0
      }
    ],
    "window_table": [
      {
        "start": 0,
        "end": 1000,
        "log2_s_over_g1": 1.0,
        "timing": "early"
      },
      {
        "start": 1000,
        "end": 2000,
        "log2_s_over_g1": 0.0,
        "timing": "intermediate"
      },
      {
        "start": 2000,
        "end": 3000,
        "log2_s_over_g1": -1.0,
        "timing": "late"
      },
      "..."
    ],
    "definition": "Replication timing score = log2( (S-phase depth / median S depth) / (G1 depth / median G1 depth) ) per window (Koren et al. 2012 style): higher = replicated ear..."
  }
}
```

## 8) plot

- **Purpose:** Plots the cumulative GC (or AT) skew of a DNA sequence with the heuristic origin and terminus marked, as an SVG string.

- **Method:** The curve is the cumulative skew computed as in `calculate`; the origin and terminus marks are the cumulative-skew minimum and maximum (a heuristic that only holds for bacterial-like chromosomes). The SVG is returned as text and written to a file only when `output_path` is given. Pure Python; no plotting library.

- **API:** `dnakit.replication_repair.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` and the skew options as in `calculate` (`window`, `step`, `skew_type`). Extra option: `output_path` (optional `.svg` file; an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `svg` (the SVG text), `output_path` (`null` unless a file was written), `skew_type`, `windows`, `origin_estimate`, `terminus_estimate` and a `note`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import plot

sequence = "GGGC" * 125 + "CCCG" * 125
result = plot(sequence, options={"window": 100})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 1000,
    "sequence_bases": 1000,
    "lengths": [
      1000
    ],
    "sha256": "3490f072d0775230a9eb8cb9e2bc13b5ea9c4ae2878fa5265873210e8ad3c480"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"760\" height=\"380\" viewBox=\"0 0 760 380\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Cumulative GC skew</title...",
    "output_path": null,
    "skew_type": "GC",
    "windows": 10,
    "origin_estimate": 950,
    "terminus_estimate": 450,
    "note": "Origin/terminus markers are cumulative-skew extrema (heuristic)."
  }
}
```

## 9) predict

- **Purpose:** Predicts candidate bacterial replication-origin (oriC) regions by the density of DnaA boxes (TTATNCACA).

- **Method:** Heuristic: DnaA boxes (consensus TTATNCACA, either strand, Hamming distance at most `max_mismatches`) are located and the non-overlapping windows with the most boxes are returned. The GC-skew origin estimate (see `calculate`) and each window's circular distance to it are reported for cross-checking. Box density alone is a weak, bacteria-specific predictor. Pure Python; no external program.

- **API:** `dnakit.replication_repair.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a DNA string, single-record FASTA text, a file path or `{"sequence": ...}`. Options: `max_mismatches` (default 0), `window` (default 500 bp), `top_n` (default 3).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `length`, `dnaa_box_consensus`, `max_mismatches`, `boxes_found`, `candidate_regions` (per window `start`, `end`, `dnaa_boxes`, `box_positions` as `[position, strand]` pairs and `circular_distance_to_gc_skew_origin`), `gc_skew_origin_estimate` and a `heuristic` description. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import predict

filler = "CGGC" * 150
sequence = filler[:600] + "TTATCCACAGG" * 5 + filler[600:]
result = predict(sequence, options={"window": 100, "top_n": 2})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "predict",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 655,
    "sequence_bases": 655,
    "lengths": [
      655
    ],
    "sha256": "d0db177af5706fbac1f08bbbfff162d69116eb8f2169590c43d9b7a5eb7a7e84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 655,
    "dnaa_box_consensus": "TTATNCACA",
    "max_mismatches": 0,
    "boxes_found": 5,
    "candidate_regions": [
      {
        "start": 600,
        "end": 655,
        "dnaa_boxes": 5,
        "box_positions": [
          [
            600,
            "+"
          ],
          [
            611,
            "+"
          ],
          [
            622,
            "+"
          ],
          [
            633,
            "+"
          ],
          "..."
        ],
        "circular_distance_to_gc_skew_origin": 57
      }
    ],
    "heuristic": "oriC is predicted as the non-overlapping windows (default 500 bp) with the most DnaA boxes (TTATNCACA on either strand, Hamming <= max_mismatches). DnaA-box den...",
    "gc_skew_origin_estimate": 570
  }
}
```

## 10) simulate

- **Purpose:** Simulates replication-fork progression from one or more origins and returns the replication-timing profile.

- **Method:** A deterministic bidirectional-fork model: every origin fires at its `firing_time` and forks move in both directions at `fork_speed`. The replication time of position x is the minimum over origins of firing_time + |x - origin| / fork_speed (passive replication by a neighbour's fork included; circular distance when `circular`). Termination sites are the points where converging forks meet. Times are evaluated at window centres, so `s_phase_duration` is the largest window-centre time and can be smaller than a termination time when windows are wide. There is no stochastic firing, fork stalling or origin licensing; random origins (`n_origins`, `seed`) are uniform positions fired at time 0.

- **API:** `dnakit.replication_repair.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is the genome length in bp (an integer, or `{"length": n}`) or a DNA sequence (its length is used). Options: `origins` (list of positions or `{position, firing_time}`; default `n_origins` random positions drawn with `seed`), `n_origins` (1), `seed` (0), `fork_speed` (bp per time unit, default 1000), `circular` (default `False`), `window` (default length/100), `time_points` (5).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `length`, `circular`, `fork_speed_bp_per_time`, `origins`, `s_phase_duration`, `mean_replication_time`, `termination_sites` (`position`, `time`), `window`, `timing_profile` (per window `start`, `end`, `time`), `replication_curve` (`time`, `fraction_replicated`) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.replication_repair import simulate

result = simulate(1000, options={"origins": [250, 750], "fork_speed": 10, "window": 100, "time_points": 3})
print(result)
```

- **Example output:**

```json
{
  "module": "replication_repair",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 4,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "40510175845988f13f6162ed8526f0b09f73384467fa855e1e79b44a56562a58"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "length": 1000,
    "circular": false,
    "fork_speed_bp_per_time": 10.0,
    "origins": [
      {
        "position": 250.0,
        "firing_time": 0.0
      },
      {
        "position": 750.0,
        "firing_time": 0.0
      }
    ],
    "s_phase_duration": 20.0,
    "mean_replication_time": 12.0,
    "termination_sites": [
      {
        "position": 500.0,
        "time": 25.0
      }
    ],
    "window": 100,
    "timing_profile": [
      {
        "start": 0,
        "end": 100,
        "time": 20.0
      },
      {
        "start": 100,
        "end": 200,
        "time": 10.0
      },
      {
        "start": 200,
        "end": 300,
        "time": 0.0
      },
      "..."
    ],
    "replication_curve": [
      {
        "time": 0.0,
        "fraction_replicated": 0.2
      },
      {
        "time": 6.666667,
        "fraction_replicated": 0.2
      },
      {
        "time": 13.333333,
        "fraction_replicated": 0.6
      },
      "..."
    ],
    "definition": "Deterministic bidirectional-fork model: all origins fire at their firing_time and forks move at fork_speed in both directions; the replication time of position ..."
  }
}
```
