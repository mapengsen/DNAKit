# Visualization

This category documents **7** functions of the `visualization` module. 5 compute a real result: `plot` draws seven kinds of sequence figures as standalone SVG, `generate_analysis_summary` computes sequence statistics, `generate_qc_html_report` and `generate_html_genome_report` build self-contained HTML reports, and `export_figure` saves an SVG (conversion to PNG/JPG/TIFF/PDF additionally needs the optional libraries `cairosvg` and `Pillow`, otherwise it returns `not_implemented`). `generate_clinical_pdf_report` always returns `status: "not_implemented"` because DNAKit cannot produce validated clinical interpretations. `save_figure` is the existing Matplotlib figure saver. Example outputs are real; `input_summary` is abbreviated to `{...}`.

## 1) export_figure

- **Purpose:** Write an SVG figure to a file, or convert it to PNG/JPG/TIFF/PDF.

- **Method:** The SVG must be well-formed XML with an `<svg>` root and numeric width/height (or a viewBox), and must not contain scripts, event-handler attributes, `foreignObject` or `javascript:` (`ConfigurationError` `UNSAFE_SVG`); DOCTYPE/ENTITY declarations are refused by the SVG artifact writer. A `.svg` target writes the exact bytes. Other suffixes (`.png`, `.jpg`/`.jpeg`, `.tif`/`.tiff`, `.pdf`) are converted with the optional libraries `cairosvg` and `Pillow`; if they are missing the call returns `status: "not_implemented"` with `required_dependency: "pip install cairosvg pillow"`. Any other suffix raises `ConfigurationError`.

- **API:** `dnakit.visualization.export_figure(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: SVG text, a `plot` result mapping (with `svg`), a whole `plot` envelope, or an `SVGArtifact`. `options`: `output` (required target path), `overwrite` (default false; an existing file raises `ConfigurationError`), `dpi` (raster export, default 96).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `path`, `format`, `bytes` and `sha256` of the written file.

- **Example code:**

```python
import json
from dnakit.visualization import export_figure, plot

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
figure = plot(fasta, options={"kind": "gc_distribution", "bins": 4})
result = export_figure(figure, options={"output": "gc.svg"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "export_figure",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "gc.svg",
    "format": "svg",
    "bytes": 3778,
    "sha256": "d9141462188a87dd8ebbaf0526446a1e705268ff8ed3d744b18c190e019186cb"
  }
}
```

## 2) generate_analysis_summary

- **Purpose:** Summarize a set of DNA sequences (counts, length statistics, GC content, base counts) as data and Markdown text.

- **Method:** `gc_percent` is `100 * (G+C) / (A+C+G+T)` over all sequences (`null` if there is no unambiguous base). `length` holds min, max, mean, median and N50 (the smallest length L such that sequences of length >= L hold at least half of all bases). `non_acgt_bases` counts every other symbol and `n_fraction` the share of `N`. `per_sequence` lists id, length and GC percent of at most the first 1000 sequences.

- **API:** `dnakit.visualization.generate_analysis_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text or a FASTA file path, a single sequence string, a list of sequences, `(id, sequence)` pairs or mappings with `id`/`sequence`, or a mapping `{id: sequence}` (a mapping with `records` is unwrapped). Sequences are upper-cased, whitespace is removed, `U` is read as `T`, and only letters and `-*.` are accepted (`InputFormatError` otherwise); at most 100,000 sequences and 500,000,000 bases. `options`: `title` (default `"Sequence analysis summary"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `sequence_count`, `total_bases`, `length`, `gc_percent`, `base_counts`, `non_acgt_bases`, `n_fraction`, `per_sequence` and `text` (the Markdown summary).

- **Example code:**

```python
import json
from dnakit.visualization import generate_analysis_summary

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = generate_analysis_summary(fasta, options={"title": "Demo"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "generate_analysis_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_count": 3,
    "total_bases": 30,
    "length": {
      "min": 8,
      "max": 12,
      "mean": 10.0,
      "median": 10.0,
      "n50": 10
    },
    "gc_percent": 46.42857142857143,
    "base_counts": {
      "A": 9,
      "C": 7,
      "G": 6,
      "N": 2,
      "T": 6
    },
    "non_acgt_bases": 2,
    "n_fraction": 0.06666666666666667,
    "per_sequence": [
      {
        "id": "a",
        "length": 10,
        "gc_percent": 50.0
      },
      {
        "id": "b",
        "length": 12,
        "gc_percent": 80.0
      },
      {
        "id": "c",
        "length": 8,
        "gc_percent": 0.0
      }
    ],
    "text": "# Demo\n\n- sequences: 3\n- total bases: 30\n- length min / median / mean / max: 8 / 10 / 10 / 12\n- N50: 10\n- GC content: 46.43 %\n- non-ACGT bases: 2\n"
  }
}
```

## 3) generate_clinical_pdf_report

- **Purpose:** Not implemented: a clinical report needs validated clinical interpretation data and a certified PDF report template.

- **Method:** No computation is performed. DNAKit cannot produce validated clinical interpretations here, so it does not generate a clinical report. For non-clinical summaries use `generate_analysis_summary`, `generate_qc_html_report` or `generate_html_genome_report`.

- **API:** `dnakit.visualization.generate_clinical_pdf_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("validated clinical variant interpretation (ClinVar/ACMG rules) and a certified report template") and a `message`; no `result`.

- **Example code:**

```python
import json
from dnakit.visualization import generate_clinical_pdf_report

result = generate_clinical_pdf_report({"variants": [{"gene": "BRCA1", "hgvs": "c.68_69del"}]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "generate_clinical_pdf_report",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated clinical variant interpretation (ClinVar/ACMG rules) and a certified report template",
  "message": "Clinical reports must contain validated interpretations, which DNAKit cannot produce here."
}
```

## 4) generate_html_genome_report

- **Purpose:** Build a self-contained HTML report of genome or assembly sequences with assembly statistics and an optional annotation summary.

- **Method:** The report is built with the package's HTML report builder from the sequences and the statistics of `generate_analysis_summary` (count, total bases, N50, GC). If `options["annotation"]` is given, `dnakit.annotation.summarize` is applied to it and its result is added under `annotation`. The HTML is returned in `html` and written only when `output` is given.

- **API:** `dnakit.visualization.generate_html_genome_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text or a FASTA file path, a single sequence string, a list of sequences, `(id, sequence)` pairs or mappings with `id`/`sequence`, or a mapping `{id: sequence}` (a mapping with `records` is unwrapped). Sequences are upper-cased, whitespace is removed, `U` is read as `T`, and only letters and `-*.` are accepted (`InputFormatError` otherwise); at most 100,000 sequences and 500,000,000 bases (contigs or chromosomes). `options`: `annotation` (GFF3/GTF text, path or features), `title`, `max_records` (default 1000), `max_total_sequence_symbols` (default 10,000,000), `output` (`.html` path; nothing is written without it), `overwrite` (default false; an existing file raises `ConfigurationError`), `include_html` (default true; false leaves out the `html` field).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `record_count`, `result_names`, `sha256` (of the HTML), `path`, `summary`, `annotation` (only with an annotation) and `html` (unless `include_html` is false; omitted in the example).

- **Example code:**

```python
import json
from dnakit.visualization import generate_html_genome_report

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
gff3 = "##gff-version 3\nchr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
result = generate_html_genome_report(fasta, options={"annotation": gff3, "include_html": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "generate_html_genome_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "record_count": 3,
    "result_names": [
      "annotation",
      "summary"
    ],
    "sha256": "c4d6a5fbf9302ce2e5d1ec158cfa3072ca43160d5192cd45a483bcd83d12f1a6",
    "path": null,
    "summary": {
      "sequence_count": 3,
      "total_bases": 30,
      "length": {
        "min": 8,
        "max": 12,
        "mean": 10.0,
        "median": 10.0,
        "n50": 10
      },
      "gc_percent": 46.42857142857143,
      "base_counts": {
        "A": 9,
        "C": 7,
        "G": 6,
        "N": 2,
        "T": 6
      },
      "non_acgt_bases": 2,
      "n_fraction": 0.06666666666666667
    },
    "annotation": {
      "format": "gff3",
      "feature_count": 1,
      "types": {
        "gene": {
          "count": 1,
          "total_bp": 4001,
          "mean_length": 4001.0
        }
      },
      "seqids": {
        "chr1": 1
      },
      "strands": {
        "+": 1
      },
      "gene_count": 1,
      "transcript_count": 0,
      "exon_count": 0,
      "intron_count": 0,
      "coding_transcripts": 0,
      "multi_exon_transcripts": 0,
      "transcripts_per_gene": {
        "mean": 0.0,
        "max": 0
      },
      "exons_per_transcript": {
        "mean": 0.0,
        "max": 0
      }
    }
  }
}
```

## 5) generate_qc_html_report

- **Purpose:** Build a self-contained HTML quality-control report of DNA sequences with per-sequence pass/fail flags.

- **Method:** A sequence fails when it is shorter than `min_length`, has an N fraction (`N` count / length) above `max_n_fraction`, or a GC percent (`100 * (G+C) / (A+C+G+T)`) outside [`min_gc_percent`, `max_gc_percent`]; a sequence without A/C/G/T has no GC value and fails any non-default range. Each failure lists its reasons. The report holds the summary statistics and the `qc` verdict. The thresholds are plain user-set cut-offs, not validated standards.

- **API:** `dnakit.visualization.generate_qc_html_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text or a FASTA file path, a single sequence string, a list of sequences, `(id, sequence)` pairs or mappings with `id`/`sequence`, or a mapping `{id: sequence}` (a mapping with `records` is unwrapped). Sequences are upper-cased, whitespace is removed, `U` is read as `T`, and only letters and `-*.` are accepted (`InputFormatError` otherwise); at most 100,000 sequences and 500,000,000 bases. `options`: `min_length` (default 1), `max_n_fraction` (default 0.1), `min_gc_percent` (default 0), `max_gc_percent` (default 100), `title`, `max_records` (default 1000), `max_total_sequence_symbols` (default 10,000,000), `output` (`.html` path; nothing is written without it), `overwrite` (default false; an existing file raises `ConfigurationError`), `include_html` (default true; false leaves out the `html` field).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `record_count`, `result_names`, `sha256`, `path`, `summary`, `qc` (`passed`, `failed_count`, `failed` with `id` and `reasons`, `thresholds`) and `html` (unless `include_html` is false; omitted in the example).

- **Example code:**

```python
import json
from dnakit.visualization import generate_qc_html_report

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = generate_qc_html_report(fasta, options={"min_length": 9, "include_html": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "generate_qc_html_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "record_count": 3,
    "result_names": [
      "qc",
      "summary"
    ],
    "sha256": "38211562920fb60a5494c1b34c47b0000a3d83a29af1d6bce34dea2d8c9681f6",
    "path": null,
    "summary": {
      "sequence_count": 3,
      "total_bases": 30,
      "length": {
        "min": 8,
        "max": 12,
        "mean": 10.0,
        "median": 10.0,
        "n50": 10
      },
      "gc_percent": 46.42857142857143,
      "base_counts": {
        "A": 9,
        "C": 7,
        "G": 6,
        "N": 2,
        "T": 6
      },
      "non_acgt_bases": 2,
      "n_fraction": 0.06666666666666667
    },
    "qc": {
      "passed": 1,
      "failed_count": 2,
      "failed": [
        {
          "id": "b",
          "reasons": [
            "N fraction 0.167 > 0.1"
          ]
        },
        {
          "id": "c",
          "reasons": [
            "length 8 < 9"
          ]
        }
      ],
      "thresholds": {
        "min_length": 9,
        "max_n_fraction": 0.1,
        "min_gc_percent": 0.0,
        "max_gc_percent": 100.0
      }
    }
  }
}
```

## 6) save_figure

- **Purpose:** Document a reusable save figure capability for Visualization workflows. Candidate CAP-0828 comes from 180 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the save figure operation, and return structured results with provenance and diagnostics. The optional `plot_kind` parameter selects one supported variant (the supported values recorded in the source catalog). The formal package reuses the existing Matplotlib figure persistence API; plot-kind and output-format variants remain options of the existing callable. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.visualization.save_figure|save_figure`.

- **API:** `dnakit.visualization.save_figure(figure[required], target[required], overwrite[optional], tight_layout[optional], close[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `figure`, `target`; optional arguments: `overwrite`, `tight_layout`, `close`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The existing formal `save_figure` function returns the staged output `Path`. It accepts a Matplotlib figure and chooses the output format from the target suffix.

- **Example code:**

```python
from pathlib import Path

import matplotlib.pyplot as plt

from dnakit.visualization import save_figure

figure, axis = plt.subplots()
axis.plot([0, 1], [0, 1])
result = save_figure(figure, Path("results/save_figure.svg"), overwrite=True)
print(result)
```

- **Example output:**

```json
{
  "status": "formal_existing",
  "capability_id": "CAP-0828",
  "operation": "save_figure",
  "result_type": "PosixPath",
  "result_preview": "results/save_figure.svg",
  "provenance": {
    "schema": "formal_dnakit",
    "implementation": "existing_visualization_api"
  }
}
```

## 7) plot

- **Purpose:** Draw a sequence-based figure as a standalone SVG (valid XML, no scripts) and return it together with the plotted numbers.

- **Method:** `kind` selects the figure, drawn by the dependency-free SVG renderers of the package: `length_histogram` (default; sequence lengths in nt), `gc_distribution` (per-sequence `100 * (G+C) / (A+C+G+T)` in `bins` equal-width bins over 0-100; sequences without an unambiguous base are skipped and counted in `skipped`), `base_composition` (proportion of each symbol over all sequences, A, C, G, T first), `sequence` (the symbols of the first sequence or the one named by `id`), `linear_map` and `circular_map` (first sequence, or the one named by `id`, with `features` given as `{type, start, end}` in 0-based half-open coordinates; the circular map always treats the sequence as circular: it forces topology `circular` whatever `topology` says, so no circular sequence needs to be requested) and `label_comparison` (paired proportion bars for `true` and `predicted` class labels). `data` holds the exact counts, edges or values drawn. Nothing is written unless `output` is given.

- **API:** `dnakit.visualization.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text or a FASTA file path, a single sequence string, a list of sequences, `(id, sequence)` pairs or mappings with `id`/`sequence`, or a mapping `{id: sequence}` (a mapping with `records` is unwrapped). Sequences are upper-cased, whitespace is removed, `U` is read as `T`, and only letters and `-*.` are accepted (`InputFormatError` otherwise); at most 100,000 sequences and 500,000,000 bases; for `label_comparison` a mapping `{"true": [...], "predicted": [...]}`. `options`: `kind` (default `"length_histogram"`; unknown kinds raise `ConfigurationError`), `bins` (default 30 for lengths, 20 for GC, 1 to 1000), `id`, `features`, `topology` (`"linear"` default, for `sequence`/`linear_map`), `title`, `output` (`.svg` path to write), `overwrite` (default false).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` is a dict with `kind`, `svg` (the SVG text), `width`, `height`, `sha256` (of the SVG), `data` (the numbers drawn; its keys depend on `kind`) and `path` (the written file or `null`). In the example the long `svg` text is abbreviated to `<svg ... </svg>`.

- **Example code:**

```python
import json
from dnakit.visualization import plot

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = plot(fasta, options={"kind": "gc_distribution", "bins": 4, "title": "GC content"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "visualization",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "gc_distribution",
    "svg": "<svg ... </svg>",
    "width": 720,
    "height": 720,
    "sha256": "78bc7424b731541585975c3e833d831387d12f7cabb3de97d995009360bc6861",
    "data": {
      "count": 3,
      "skipped": 0,
      "gc_percent": [
        50.0,
        80.0,
        0.0
      ],
      "bin_edges": [
        0.0,
        25.0,
        50.0,
        75.0,
        100.0
      ],
      "counts": [
        1,
        0,
        1,
        1
      ]
    },
    "path": null
  }
}
```
