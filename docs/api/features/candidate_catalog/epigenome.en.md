# Epigenomics

This category contains **13** documented candidate capabilities; 11 compute real results natively (BED/bedGraph/Bismark parsing, interval merging and consensus peaks, annotation, methylation levels and QC, Fisher-test differential methylation, CpG statistics, normalisation, plotting and enhancer-gene linking) and 2 (`align`, `predict`) always return `status: "not_implemented"` because they need an aligner or trained model weights. Some results are explicit heuristics: `link_enhancers_to_genes` is a distance-to-TSS proximity rule, not evidence of regulation, and `call` and `build_chromatin_state_segments` only combine peak or state calls that are already given (they do not call peaks or learn states). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) align

- **Purpose:** Align epigenomic reads to a reference. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed and no aligner is wired into this function. Bisulfite data would need Bismark or bwa-meth, ChIP/ATAC data Bowtie2 or BWA, each with an indexed reference. Nothing is fabricated.

- **API:** `dnakit.epigenome.align(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a read aligner (Bismark/bwa-meth for bisulfite, Bowtie2/BWA for ChIP/ATAC) and an indexed reference` and `message` is `Read alignment is not implemented natively; no aligner is wired into this function.`

- **Example code:**

```python
import json
from dnakit.epigenome import align

result = align(["ACGTACGTACGTACGT"], options={"reference": "genome.fa"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "align",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "87cd4b17356c3dee4b13d54fba3beefea5678a7b651616a23b27c782da631af2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a read aligner (Bismark/bwa-meth for bisulfite, Bowtie2/BWA for ChIP/ATAC) and an indexed reference",
  "message": "Read alignment is not implemented natively; no aligner is wired into this function."
}
```

## 2) annotate

- **Purpose:** Annotate genomic regions (for example peaks) with the overlapping and the nearest features from a supplied annotation.

- **Method:** For every region the features on the same chromosome are scanned (quadratic search, at most 50 million region x feature pairs). A feature overlaps when `feature.start < region.end` and `feature.end > region.start`; `nearest` is the feature with the smallest gap (0 when overlapping or touching, otherwise the number of bases between the two intervals; the first feature in sorted order wins ties), with `side` `overlap`, `upstream` (feature lies before the region) or `downstream` (after it) in forward-strand coordinates. Features on other chromosomes are ignored and gene strand is not used. The annotation is whatever is supplied; nothing is looked up.

- **API:** `dnakit.epigenome.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: regions, intervals as BED text/path, rows `[chrom, start, end, ...]` or mappings with `chrom`/`start`/`end` (0-based half-open, `end` > `start` >= 0; `track`, `browser` and `#` lines and blanks are skipped; at most 2 million records). `options`: `features` (required: the annotation as BED with a name column, rows or mappings with `name`); `max_distance` (optional: a nearest feature farther than this is not reported); `entity` (only `"peaks"`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result.regions` has one entry per region: `chrom`, `start`, `end`, `overlapping` (feature names) and `nearest` (`name`, `distance`, `side`, or `null` if none exists or it is beyond `max_distance`).

- **Example code:**

```python
import json
from dnakit.epigenome import annotate

peaks = "chr1\t100\t200\tpeak1\nchr1\t1000\t1100\tpeak2\n"
features = "chr1\t150\t400\tgeneA\nchr1\t5000\t6000\tgeneB\n"
result = annotate(peaks, options={"features": features, "max_distance": 10000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1e0b18cc52249734849452cd1247dd17a4569b367d0b45b2f578d88f9e3f1d25"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "regions": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "overlapping": [
          "geneA"
        ],
        "nearest": {
          "name": "geneA",
          "distance": 0,
          "side": "overlap"
        }
      },
      {
        "chrom": "chr1",
        "start": 1000,
        "end": 1100,
        "overlapping": [],
        "nearest": {
          "name": "geneA",
          "distance": 600,
          "side": "upstream"
        }
      }
    ]
  }
}
```

## 3) build_chromatin_state_segments

- **Purpose:** Collapse per-bin chromatin-state calls into contiguous state segments.

- **Method:** Bins are sorted by chromosome, start and end; consecutive bins of the same state on one chromosome are joined when the gap to the previous segment is <= `max_gap` bases (overlapping bins of the same state are also joined). This only segments states that are already given; it does not learn or call states (it is not ChromHMM).

- **API:** `dnakit.epigenome.build_chromatin_state_segments(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: bins as BED text/path, rows or mappings `chrom, start, end, name` where `name` is the state label (`state` is accepted as an alias key in mappings), or `{chrom: [state, state, ...]}` with `options["bin_size"]` (bins start at 0 and are consecutive). `options`: `max_gap` (default 0), `bin_size` (required for the `{chrom: [states]}` form).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `bin_count`, `segment_count` and `segments` (`chrom`, `start`, `end`, `state`, `n_bins`) sorted by position.

- **Example code:**

```python
import json
from dnakit.epigenome import build_chromatin_state_segments

states = {"chr1": ["Quiescent", "Enhancer", "Enhancer", "Promoter", "Promoter", "Promoter"]}
result = build_chromatin_state_segments(states, options={"bin_size": 200})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "build_chromatin_state_segments",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 83,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b24877c056115a3c3ed7a4314d1dc9a56d57e546ddbf588405cab52d2608c84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_count": 6,
    "segment_count": 3,
    "segments": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 200,
        "state": "Quiescent",
        "n_bins": 1
      },
      {
        "chrom": "chr1",
        "start": 200,
        "end": 600,
        "state": "Enhancer",
        "n_bins": 2
      },
      {
        "chrom": "chr1",
        "start": 600,
        "end": 1200,
        "state": "Promoter",
        "n_bins": 3
      }
    ]
  }
}
```

## 4) calculate

- **Purpose:** Calculate methylation levels from counts, CpG observed/expected ratios from DNA, or the overlap of two interval sets (`options["entity"]`).

- **Method:** `methylation_level` (default): level = methylated / (methylated + unmethylated) per site; sites with coverage below `min_coverage` are dropped; `one_based` shifts Bismark's 1-based positions to 0-based. With `regions`, sites inside each 0-based half-open region are aggregated: `weighted_level` = sum(methylated) / sum(coverage) and `mean_site_level` = unweighted mean of the site levels. `cpg_observed_expected`: counts C, G and CpG, GC fraction and CpG o/e = CpG * L / (C * G) (`null` when C*G = 0), also per sliding window when `window` is set. `interval_overlap`: each set is merged first, then `overlap_bp`, `union_bp` and Jaccard = overlap / union are computed in bases. `a_intervals_overlapping_b` counts the original (unmerged) input intervals that share at least one base with any interval of `b` (half-open, so book-ended intervals do not touch); each input interval is counted once, however many `b` intervals it spans.

- **API:** `dnakit.epigenome.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `methylation_level`: `input` is per-site counts: Bismark `.cov` text/path or rows `chrom, start, end, percent, methylated, unmethylated`, or mappings with `chrom`, `start`, `methylated` and `unmethylated` (or `coverage`); `options`: `min_coverage` (default 1), `one_based` (default `false`; use `true` for `.cov`), `regions` (optional BED). `cpg_observed_expected`: `input` is DNA (string, FASTA text/path or `{id: sequence}`, ACGT; at most 10000 sequences and 50 Mb); `options`: `window` (optional), `step` (default = window). `interval_overlap`: `input` and `options["b"]` (required) are interval sets (intervals as BED text/path, rows `[chrom, start, end, ...]` or mappings with `chrom`/`start`/`end`). `options["entity"]` is `"methylation_level"` (default), `"cpg_observed_expected"` or `"interval_overlap"`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result.entity` echoes the choice. `methylation_level` without `regions`: `sites` (`chrom`, `pos`, `coverage`, `methylated`, `level`), `site_count`, `dropped_low_coverage`; with `regions`: `regions` (`chrom`, `start`, `end`, `n_sites`, `coverage`, `weighted_level`, `mean_site_level`). `cpg_observed_expected`: `sequences` (`id`, `length`, `c`, `g`, `cpg`, `gc_fraction`, `observed_expected`, and `windows` when requested). `interval_overlap`: `overlap_bp`, `a_bp`, `b_bp`, `union_bp`, `jaccard`, `a_intervals_overlapping_b`.

- **Example code:**

```python
import json
from dnakit.epigenome import calculate

cov = "chr1\t101\t101\t66.7\t2\t1\nchr1\t120\t120\t100\t5\t0\nchr1\t5000\t5000\t0\t0\t4\n"   # Bismark .cov rows
result = calculate(cov, options={"entity": "methylation_level", "one_based": True,
                                 "regions": "chr1\t0\t1000\tregionA"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "22911cbbf5fca89bad46f6d2e400272575fc4bd02d8613cfdf73a18d901f0f50"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "methylation_level",
    "regions": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 1000,
        "n_sites": 2,
        "coverage": 8,
        "weighted_level": 0.875,
        "mean_site_level": 0.8333333333333333
      }
    ]
  }
}
```

## 5) call

- **Purpose:** Call consensus peaks across samples: regions covered by the peak sets of at least `min_samples` samples.

- **Method:** Each sample's own overlapping peaks are merged first so that a sample votes at most once per base; the result is every maximal segment whose vote count is >= `min_samples` (0-based half-open), computed by a sweep over interval start and end events. This only combines existing peak calls; it does not call peaks from signal (it is not MACS), and it uses no statistics.

- **API:** `dnakit.epigenome.call(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: peak sets as a mapping `{sample: intervals}` (each intervals as BED text/path, rows `[chrom, start, end, ...]` or mappings with `chrom`/`start`/`end`) or a list of interval sets. `options`: `entity` (only `"consensus_peaks"`); `min_samples` (default 2, at most the number of samples).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `samples`, `min_samples`, `peak_count` and `peaks` (`chrom`, `start`, `end`).

- **Example code:**

```python
import json
from dnakit.epigenome import call

peak_sets = {
    "rep1": "chr1\t100\t300\nchr1\t1000\t1200\n",
    "rep2": "chr1\t200\t400\n",
    "rep3": "chr1\t250\t350\nchr1\t1100\t1300\n",
}
result = call(peak_sets, options={"entity": "consensus_peaks", "min_samples": 2})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "call",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 120,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6e5f5041ddc207500e99845fb2c45b47675903c48b38032e44d579799a510e75"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "rep1",
      "rep2",
      "rep3"
    ],
    "min_samples": 2,
    "peak_count": 2,
    "peaks": [
      {
        "chrom": "chr1",
        "start": 200,
        "end": 350
      },
      {
        "chrom": "chr1",
        "start": 1100,
        "end": 1200
      }
    ]
  }
}
```

## 6) find

- **Purpose:** Find differentially methylated regions by Fisher's exact test, or CpG islands in DNA (`options["entity"]`).

- **Method:** `differential_regions` (default): per region the pooled counts of group A and B give the group levels, `delta` = level_b - level_a and a two-sided Fisher exact test on the 2x2 table [[methylated_a, unmethylated_a], [methylated_b, unmethylated_b]] (the scipy definition: sum of table probabilities not larger than the observed one), the sample odds ratio ad/bc (`null` if bc = 0) and Benjamini-Hochberg `q_value` over the tested regions. Regions with coverage below `min_coverage` in either group are excluded; `significant` means q <= `alpha` and |delta| >= `min_delta`. Counts are pooled over replicates and cytosines, so replicate variability is not modelled (it is not a beta-binomial/DSS-style model). `cpg_islands`: delegates to `dnakit.patterns.find_cpg_islands`, a sliding-window rule (defaults: window 200, step 1, GC >= 0.5, CpG observed/expected >= 0.6, region length >= 200).

- **API:** `dnakit.epigenome.find(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `differential_regions`: `input` is rows or mappings with `chrom`, `start`, `end` (optional) and the pooled counts `methylated_a`, `unmethylated_a`, `methylated_b`, `unmethylated_b`; with `options["regions"]` (BED) the rows are single sites that are summed inside each region. `options`: `min_coverage` (default 10, per group), `alpha` (default 0.05), `min_delta` (default 0.0). `cpg_islands`: `input` is one DNA sequence (string, FASTA text/path or `{id: sequence}`); `options`: `window_size` (200), `step` (1), `min_gc` (0.5), `min_observed_expected` (0.6), `min_region_length` (200). `options["entity"]` is `"differential_regions"` (default) or `"cpg_islands"`.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `differential_regions`: `result` has `entity`, `tested`, `excluded_low_coverage` and `regions` sorted by p-value (`chrom`, `start`, `end`, `level_a`, `level_b`, `delta`, `coverage_a`, `coverage_b`, `p_value`, `odds_ratio`, `q_value`, `significant`). `cpg_islands`: `entity`, `sequence` and `islands` (`start`, `end`, `length`, `observed_expected`, `gc_fraction`).

- **Example code:**

```python
import json
from dnakit.epigenome import find

regions = [
    {"chrom": "chr1", "start": 100, "end": 200,
     "methylated_a": 40, "unmethylated_a": 60, "methylated_b": 80, "unmethylated_b": 20},
    {"chrom": "chr1", "start": 500, "end": 600,
     "methylated_a": 50, "unmethylated_a": 50, "methylated_b": 52, "unmethylated_b": 48},
]
result = find(regions, options={"entity": "differential_regions"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "find",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 258,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "759876ad6dfe28bf2acab3945ddd4e23a745b2f8098b355e76502f0a3913a0a7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "differential_regions",
    "tested": 2,
    "excluded_low_coverage": 0,
    "regions": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "level_a": 0.4,
        "level_b": 0.8,
        "delta": 0.4,
        "coverage_a": 100,
        "coverage_b": 100,
        "p_value": 1.0636025781775309e-08,
        "odds_ratio": 0.16666666666666666,
        "q_value": 2.1272051563550617e-08,
        "significant": true
      },
      {
        "chrom": "chr1",
        "start": 500,
        "end": 600,
        "level_a": 0.5,
        "level_b": 0.52,
        "delta": 0.020000000000000018,
        "coverage_a": 100,
        "coverage_b": 100,
        "p_value": 0.8875621593770635,
        "odds_ratio": 0.9230769230769231,
        "q_value": 0.8875621593770635,
        "significant": false
      }
    ]
  }
}
```

## 7) link_enhancers_to_genes

- **Purpose:** Assign enhancers (or any regions) to genes by genomic distance to the gene transcription start site (TSS).

- **Method:** Distance-based assignment: the distance is |TSS - midpoint| with midpoint = (start + end) // 2 of the enhancer. The TSS is the explicit `tss` if given, otherwise `start` for + or unstranded genes and `end - 1` for - genes. `method="nearest"` links each enhancer to the single nearest TSS on the same chromosome (ties go to the alphabetically first gene name) if it lies within `max_distance`; `method="window"` links every gene whose TSS lies within `max_distance`. `direction` tells whether the enhancer is upstream, downstream or at the TSS with respect to the gene orientation. This is a proximity heuristic, not evidence of regulation: no Hi-C, expression correlation or TAD boundary is used, and the call adds a warning to this effect.

- **API:** `dnakit.epigenome.link_enhancers_to_genes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: enhancers, intervals as BED text/path, rows `[chrom, start, end, ...]` or mappings with `chrom`/`start`/`end` (0-based half-open, `end` > `start` >= 0; `track`, `browser` and `#` lines and blanks are skipped; at most 2 million records). `options`: `genes` (required: rows or mappings with `chrom`, `start`, `end`, `name` and optional `strand` or explicit `tss`); `method` (`"nearest"` default | `"window"`); `max_distance` (default 1,000,000 bp; for `nearest`, farther genes are not linked).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `method`, `max_distance` and `links` (one entry per enhancer: `region`, `midpoint` and `genes`, each with `gene`, `tss`, `distance`, `direction`, `strand`; the list is empty when no gene is within range).

- **Example code:**

```python
import json
from dnakit.epigenome import link_enhancers_to_genes

enhancers = "chr1\t10000\t10400\tenh1\nchr1\t90000\t90200\tenh2\n"
genes = [
    {"name": "GENE_A", "chrom": "chr1", "start": 12000, "end": 20000, "strand": "+"},
    {"name": "GENE_B", "chrom": "chr1", "start": 60000, "end": 91000, "strand": "-"},
]
result = link_enhancers_to_genes(enhancers, options={"genes": genes, "method": "nearest", "max_distance": 50000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "link_enhancers_to_genes",
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

## 8) merge

- **Purpose:** Merge overlapping or nearby intervals (peaks, regions) per chromosome, like `bedtools merge`.

- **Method:** Intervals are sorted per chromosome and joined when the gap to the current merged interval is <= `distance` (so with 0 overlapping and book-ended intervals merge). Strand and names are ignored. Output is sorted by chromosome name and start, with `n_merged` counting the input intervals combined.

- **API:** `dnakit.epigenome.merge(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: intervals as BED text/path, rows `[chrom, start, end, ...]` or mappings with `chrom`/`start`/`end` (0-based half-open, `end` > `start` >= 0; `track`, `browser` and `#` lines and blanks are skipped; at most 2 million records). `options`: `entity` (only `"peaks"`); `distance` (default 0).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `distance`, `intervals` (`chrom`, `start`, `end`, `n_merged`) and `count`.

- **Example code:**

```python
import json
from dnakit.epigenome import merge

peaks = "chr1\t100\t200\nchr1\t150\t300\nchr1\t310\t400\nchr2\t5\t50\n"
result = merge(peaks, options={"distance": 20})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "merge",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 49,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c1640184d9b52456608b48132dcfc54c62503331851a4b0255849c89789202c3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "distance": 20,
    "intervals": [
      {
        "chrom": "chr1",
        "start": 100,
        "end": 400,
        "n_merged": 3
      },
      {
        "chrom": "chr2",
        "start": 5,
        "end": 50,
        "n_merged": 1
      }
    ],
    "count": 2
  }
}
```

## 9) normalize

- **Purpose:** Normalise signal values (bedGraph records or a plain list of numbers).

- **Method:** `cpm` (default): value / sum(values) * 1e6; `rpkm`: cpm / (interval length / 1000), which needs intervals; `log2cpm`: log2(cpm + `pseudocount`); `zscore`: (v - mean) / sd with the population sd unless `ddof` is 1; `minmax`: (v - min) / (max - min). The cpm family needs non-negative values and a positive sum, z-score a non-zero variance and min-max different values, otherwise an input error is raised. The normalisation is applied to the supplied values only, with no library-size or GC correction.

- **API:** `dnakit.epigenome.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: bedGraph text/path, rows `[chrom, start, end, value]` or mappings with a `value` (or `signal`) field, or a list of numbers (not for `rpkm`). `options`: `entity` (only `"signal"`); `method` (`"cpm"` default | `"rpkm"` | `"log2cpm"` | `"zscore"` | `"minmax"`); `pseudocount` (default 1, `log2cpm`); `ddof` (0 or 1, default 0, `zscore`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `method` and `records` (the input records with `value` replaced by the normalised value and `raw_value` keeping the original).

- **Example code:**

```python
import json
from dnakit.epigenome import normalize

signal = "chr1\t0\t100\t5\nchr1\t100\t200\t15\nchr1\t200\t400\t30\n"
result = normalize(signal, options={"method": "rpkm"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 45,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fa1a617de98e7c794ff977eaf8ac2571d21211bfe540f311c20841acfddc1a81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "rpkm",
    "records": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 100,
        "value": 1000000.0,
        "raw_value": 5.0
      },
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "value": 3000000.0,
        "raw_value": 15.0
      },
      {
        "chrom": "chr1",
        "start": 200,
        "end": 400,
        "value": 3000000.0,
        "raw_value": 30.0
      }
    ]
  }
}
```

## 10) plot

- **Purpose:** Draw a signal track (bedGraph bins as bars) for one genomic window as an SVG string.

- **Method:** A hand-written SVG (no plotting library): each bedGraph bin that overlaps the window becomes a `<rect class="bin">` with `data-start`, `data-end` and `data-value` attributes, positioned linearly in the window and with a height scaled to the largest value in the window. The window label shows the maximum. Only the first chromosome (or `chrom`) is drawn.

- **API:** `dnakit.epigenome.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: bedGraph text/path, rows `[chrom, start, end, value]` or mappings with a `value` (or `signal`) field. `options`: `entity` (only `"tracks"`); `chrom`, `start`, `end` (default: the first chromosome and the span of its bins); `output` (optional `.svg` path).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `svg` (the SVG text) and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.epigenome import plot

signal = "chr1\t0\t100\t5\nchr1\t100\t200\t15\nchr1\t200\t400\t30\n"
result = plot(signal, options={"chrom": "chr1", "start": 0, "end": 400})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "epigenome",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 45,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fa1a617de98e7c794ff977eaf8ac2571d21211bfe540f311c20841acfddc1a81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"150\"><text x=\"20\" y=\"12\" font-size=\"10\">chr1:0-400 (max 30)</text><rec...",
    "output_path": null
  }
}
```

## 11) predict

- **Purpose:** Predict epigenomic signal from sequence. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed. Sequence-to-epigenome prediction needs a trained model such as Enformer, Basenji or DeepSEA together with its weights, which DNAKit does not bundle. Nothing is fabricated.

- **API:** `dnakit.epigenome.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarised). `options`: ignored.

- **Output:** The call returns `status: "not_implemented"` (and `provenance.implementation` is `"not_implemented"`) with the usual envelope fields but no `result`. `required_dependency` is `a trained sequence-to-epigenome model and its weights (Enformer, Basenji, DeepSEA)` and `message` is `Epigenome prediction needs trained model weights that are not bundled.`

- **Example code:**

```python
import json
from dnakit.epigenome import predict

result = predict("ACGTACGTACGTACGT")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
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
  "required_dependency": "a trained sequence-to-epigenome model and its weights (Enformer, Basenji, DeepSEA)",
  "message": "Epigenome prediction needs trained model weights that are not bundled."
}
```

## 12) qc

- **Purpose:** Quality-check bisulfite methylation counts: coverage distribution, global methylation and conversion rate.

- **Method:** Control chromosomes (`control_chroms`, for example an unmethylated lambda phage or chrL) are separated from the sample sites. For the sample sites: number of sites, total, mean and median coverage, the fraction of sites with coverage >= `min_coverage` and `global_methylation` = sum(methylated) / sum(coverage). With controls, `conversion_rate` = 1 - methylated / coverage on the control sites, i.e. the share of unmethylated cytosines read as T. Plain summary statistics; no thresholds are applied beyond reporting.

- **API:** `dnakit.epigenome.qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: per-site counts: Bismark `.cov` text/path or rows `chrom, start, end, percent, methylated, unmethylated`, or mappings with `chrom`, `start`, `methylated` and `unmethylated` (or `coverage`). `options`: `entity` (only `"methylation"`); `min_coverage` (default 10); `control_chroms` (optional list of control chromosome names).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `sites`, `total_coverage`, `mean_coverage`, `median_coverage`, `min_coverage`, `fraction_sites_min_coverage`, `global_methylation` and, with controls, `control_coverage` and `conversion_rate`. If no sample sites remain an input error is raised.

- **Example code:**

```python
import json
from dnakit.epigenome import qc

cov = ("chr1\t101\t101\t66.7\t20\t10\nchr1\t120\t120\t100\t30\t0\nchr1\t5000\t5000\t0\t0\t4\n"
       "chrL\t10\t10\t1.0\t1\t99\n")
result = qc(cov, options={"min_coverage": 10, "control_chroms": ["chrL"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 87,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d643d3abac8dc0545d26854a692a340ad0c42126f53306d4fc3057691324ee32"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sites": 3,
    "total_coverage": 64,
    "mean_coverage": 21.333333333333332,
    "median_coverage": 30.0,
    "min_coverage": 10,
    "fraction_sites_min_coverage": 0.6666666666666666,
    "global_methylation": 0.78125,
    "control_coverage": 100,
    "conversion_rate": 0.99
  }
}
```

## 13) read

- **Purpose:** Read BED, bedGraph or Bismark coverage (per-site methylation counts) text, a file path or in-memory rows.

- **Method:** `track`, `browser` and `#` lines and blank lines are skipped; cells are split at tabs (or whitespace). Records are validated (integer coordinates, `end` > `start` >= 0, finite values, strand `+`/`-`/`.`) and malformed lines raise an error with the line number. bed/bedgraph coordinates are 0-based half-open; `bismark_cov` positions are returned as given (1-based) and only chrom, start, methylated and unmethylated are kept.

- **API:** `dnakit.epigenome.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text, a path, rows of cells or mappings. `options`: `format` (`"bed"` default: chrom, start, end, optional name/score/strand | `"bedgraph"`: chrom, start, end, value | `"bismark_cov"`: chrom, start, end, percent, methylated, unmethylated).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"epigenome"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `format`, `count` and `records`.

- **Example code:**

```python
import json
from dnakit.epigenome import read

bedgraph = "track type=bedGraph\nchr1\t0\t100\t1.5\nchr1\t100\t200\t3\n"
result = read(bedgraph, options={"format": "bedgraph"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "epigenome",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50e7b9157f4acdeffe27e5b4d0ac53b9255932da665686b2a6158b7f2dad752d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "bedgraph",
    "count": 2,
    "records": [
      {
        "chrom": "chr1",
        "start": 0,
        "end": 100,
        "value": 1.5
      },
      {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "value": 3.0
      }
    ]
  }
}
```
