# Annotation

This category documents **18** functions of the `annotation` module; 14 of them are covered by the text below (the sections of `convert_annotation_coordinates`, `make_intergenic_regions`, `make_promoter_regions` and `merge_gene_models` are separate). 12 of these 14 compute a real result from GFF3/GTF text, files or feature mappings, without any database: `parse`, `validate`, `check_annotation_consistency`, `build_gene_model`, `index`, `summarize`, `compare`, `export`, `extract`, `annotate` and the two lookups `assign_gene_ontology` / `assign_pathways`. Among them `annotate` needs `options["annotation"]`, and `assign_gene_ontology` and `assign_pathways` need `options["mapping"]` or `options["gene_sets"]`; without them they return `status: "not_implemented"`, because DNAKit bundles no annotation, Gene Ontology or pathway database. `assign_protein_domains` and `find_homologous_genes` always return `status: "not_implemented"` (they would need Pfam/InterPro with an HMM search, or an orthology database / BLAST). All coordinates are 1-based and closed. Example outputs are real; `input_summary` is abbreviated to `{...}`.

## 1) annotate

- **Purpose:** Annotate genomic regions with the overlapping features of a supplied annotation and, when none overlaps, the nearest feature.

- **Method:** Features of the selected types are indexed per sequence by start coordinate; a region overlaps a feature when `feature.start <= region.end` and `feature.end >= region.start`, and `overlap_bp` is the length of the shared interval. If nothing overlaps, the nearest feature to the left or right is reported with its `distance` in bp (region start minus feature end, or feature start minus region end); `relation` is `upstream`/`downstream` relative to the feature strand, or `before`/`after` for unstranded features; features farther than `max_distance` are not reported. All coordinates are 1-based and closed. Without `options["annotation"]` the function returns `status: "not_implemented"`, because DNAKit bundles no annotation database.

- **API:** `dnakit.annotation.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: regions as strings `"chr1:100-200"` / `"chr1:150"` (commas allowed in numbers), mappings `{"seqid"|"chrom", "start", "end", "id"}`, a single one of these, or a mapping with `regions` (at most 100,000). `options`: `annotation` (required: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`)), `feature_types` (default `["gene"]`, case-insensitive), `max_distance` (default unlimited), `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `feature_types` and `annotations`, one row per region with `region`, `overlapping` (feature `id`, `name`, `type`, `seqid`, `start`, `end`, `strand`, `overlap_bp`) and `nearest` (the same fields plus `distance` and `relation`, or `null`). Without `annotation`: `status: "not_implemented"` with `required_dependency` and `message`.

- **Example code:**

```python
import json
from dnakit.annotation import annotate

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = annotate(["chr1:1100-1150", "chr1:6000-6100"], options={"annotation": gff3, "feature_types": ["gene"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_types": [
      "gene"
    ],
    "annotations": [
      {
        "region": {
          "id": "chr1:1100-1150",
          "seqid": "chr1",
          "start": 1100,
          "end": 1150
        },
        "overlapping": [
          {
            "id": "g1",
            "name": "GENE1",
            "type": "gene",
            "seqid": "chr1",
            "start": 1000,
            "end": 5000,
            "strand": "+",
            "overlap_bp": 51
          }
        ],
        "nearest": null
      },
      {
        "region": {
          "id": "chr1:6000-6100",
          "seqid": "chr1",
          "start": 6000,
          "end": 6100
        },
        "overlapping": [],
        "nearest": {
          "id": "g1",
          "name": "GENE1",
          "type": "gene",
          "seqid": "chr1",
          "start": 1000,
          "end": 5000,
          "strand": "+",
          "distance": 1000,
          "relation": "downstream"
        }
      }
    ]
  }
}
```

## 2) assign_gene_ontology

- **Purpose:** Look up the Gene Ontology terms of genes in a gene-to-GO mapping that the caller supplies; no GO database is bundled.

- **Method:** This is a table lookup, not an enrichment or prediction: the terms come entirely from `mapping` and/or `gene_sets`, which are merged into one gene -> terms table. Every term must look like `GO:` followed by seven digits (`InputFormatError` `INVALID_GO_TERM` otherwise); the GO hierarchy is not consulted. Genes without an entry are listed under `unannotated`. Without `mapping`/`gene_sets` the function returns `status: "not_implemented"`.

- **API:** `dnakit.annotation.assign_gene_ontology(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a gene id, a comma/whitespace separated string, a list of gene ids, or a mapping with `genes`. `options`: `mapping` (`{gene: [terms]}`) and/or `gene_sets` (`{term: [genes]}`), at least one of them is required for a result; `term_names` (`{term: name}`, optional labels), `case_insensitive` (default false: gene names match case-sensitively).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `assignments` (gene -> sorted list of terms), `terms` (term -> `name` and the `genes` carrying it), `unannotated` (genes without any term) and `annotated_fraction`. Without `mapping` and `gene_sets` the envelope has `status: "not_implemented"`, `required_dependency` and `message` instead.

- **Example code:**

```python
import json
from dnakit.annotation import assign_gene_ontology

result = assign_gene_ontology(
    "TP53, BRCA1, UNKNOWN",
    options={
        "mapping": {"TP53": ["GO:0006915", "GO:0003700"], "BRCA1": ["GO:0006281"]},
        "term_names": {"GO:0006281": "DNA repair"},
    },
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "assign_gene_ontology",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "assignments": {
      "TP53": [
        "GO:0003700",
        "GO:0006915"
      ],
      "BRCA1": [
        "GO:0006281"
      ],
      "UNKNOWN": []
    },
    "terms": {
      "GO:0003700": {
        "name": null,
        "genes": [
          "TP53"
        ]
      },
      "GO:0006281": {
        "name": "DNA repair",
        "genes": [
          "BRCA1"
        ]
      },
      "GO:0006915": {
        "name": null,
        "genes": [
          "TP53"
        ]
      }
    },
    "unannotated": [
      "UNKNOWN"
    ],
    "annotated_fraction": 0.6666666666666666
  }
}
```

## 3) assign_pathways

- **Purpose:** Look up the pathways of genes in a gene-to-pathway mapping that the caller supplies; no pathway database is bundled.

- **Method:** The same table lookup as `assign_gene_ontology`, except that pathway identifiers are free-form strings (for example KEGG or Reactome ids or names) and are not format-checked. It is not an enrichment test. Without `mapping`/`gene_sets` the function returns `status: "not_implemented"`.

- **API:** `dnakit.annotation.assign_pathways(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a gene id, a comma/whitespace separated string, a list of gene ids, or a mapping with `genes`. `options`: `mapping` (`{gene: [terms]}`) and/or `gene_sets` (`{term: [genes]}`), at least one of them is required for a result; `term_names` (`{term: name}`, optional labels), `case_insensitive` (default false: gene names match case-sensitively).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `assignments` (gene -> sorted list of terms), `terms` (term -> `name` and the `genes` carrying it), `unannotated` (genes without any term) and `annotated_fraction`. Without `mapping` and `gene_sets` the envelope has `status: "not_implemented"`, `required_dependency` and `message` instead.

- **Example code:**

```python
import json
from dnakit.annotation import assign_pathways

result = assign_pathways(
    ["TP53", "BRCA1", "UNKNOWN"],
    options={"gene_sets": {"p53 signaling": ["TP53", "MDM2"], "DNA repair": ["BRCA1", "TP53"]}},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "assign_pathways",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "assignments": {
      "TP53": [
        "DNA repair",
        "p53 signaling"
      ],
      "BRCA1": [
        "DNA repair"
      ],
      "UNKNOWN": []
    },
    "terms": {
      "DNA repair": {
        "name": null,
        "genes": [
          "TP53",
          "BRCA1"
        ]
      },
      "p53 signaling": {
        "name": null,
        "genes": [
          "TP53"
        ]
      }
    },
    "unannotated": [
      "UNKNOWN"
    ],
    "annotated_fraction": 0.6666666666666666
  }
}
```

## 4) assign_protein_domains

- **Purpose:** Not implemented: protein domains need profile-HMM or domain databases (Pfam, InterPro) and a search program.

- **Method:** No computation is performed; neither a domain database nor an HMM search is bundled.

- **API:** `dnakit.annotation.assign_protein_domains(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("Pfam/InterPro profile databases and hmmscan or InterProScan") and a `message` (Domain assignment needs a domain database and an HMM search; neither is bundled); no `result`.

- **Example code:**

```python
import json
from dnakit.annotation import assign_protein_domains

result = assign_protein_domains("MKTAYIAKQRQISFVKSHFSRQ")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "assign_protein_domains",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Pfam/InterPro profile databases and hmmscan or InterProScan",
  "message": "Domain assignment needs a domain database and an HMM search; neither is bundled."
}
```

## 5) build_gene_model

- **Purpose:** Assemble gene models (gene, transcripts, exons, introns, CDS, UTRs) from GFF3/GTF text or feature mappings.

- **Method:** Transcripts are mRNA/transcript/`*RNA` features, linked to their exons and CDS through GFF3 `Parent` or GTF `transcript_id`. Introns are the gaps between consecutive exons. If a transcript has CDS (and, for GTF, start/stop codon) features but no annotated UTR, the 5' and 3' UTRs are derived as the exon parts outside the CDS span and `utr_derived` is true (in GTF a stop codon outside the CDS is taken into account). `protein_length_aa` is `cds_length / 3` when the CDS length is divisible by 3, otherwise `null`. Intervals are `[start, end]`, 1-based closed, ascending genomic order. No sequence is read, so the protein length is a length only.

- **API:** `dnakit.annotation.build_gene_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`). `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `gene_count`, `transcript_count` and `genes`, each with `gene_id`, `name`, `type`, `seqid`, `start`, `end`, `strand`, `biotype` and `transcripts` (`transcript_id`, `type`, `seqid`, `start`, `end`, `strand`, `name`, `exons`, `introns`, `cds`, `five_prime_utr`, `three_prime_utr`, `utr_derived`, `exon_count`, `spliced_length`, `cds_length`, `protein_length_aa`).

- **Example code:**

```python
import json
from dnakit.annotation import build_gene_model

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = build_gene_model(gff3)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "build_gene_model",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "gene_count": 1,
    "transcript_count": 1,
    "genes": [
      {
        "gene_id": "g1",
        "name": "GENE1",
        "type": "gene",
        "seqid": "chr1",
        "start": 1000,
        "end": 5000,
        "strand": "+",
        "biotype": null,
        "transcripts": [
          {
            "transcript_id": "t1",
            "type": "mRNA",
            "seqid": "chr1",
            "start": 1000,
            "end": 5000,
            "strand": "+",
            "name": null,
            "exons": [
              [
                1000,
                1200
              ],
              [
                2000,
                2300
              ]
            ],
            "introns": [
              [
                1201,
                1999
              ]
            ],
            "cds": [
              [
                1101,
                1200
              ],
              [
                2000,
                2199
              ]
            ],
            "five_prime_utr": [
              [
                1000,
                1100
              ]
            ],
            "three_prime_utr": [
              [
                2200,
                2300
              ]
            ],
            "utr_derived": true,
            "exon_count": 2,
            "spliced_length": 502,
            "cds_length": 300,
            "protein_length_aa": 100
          }
        ]
      }
    ]
  }
}
```

## 6) check_annotation_consistency

- **Purpose:** Check the logical consistency of a gene annotation: children inside parents, exon and CDS geometry, and CDS phases.

- **Method:** Rules (code, severity): `DUPLICATE_ID` (GFF3, error), `MISSING_PARENT` (error), `CHILD_OUTSIDE_PARENT` (transcript outside gene, exon/CDS/UTR outside transcript; error), `SEQID_MISMATCH` and `STRAND_MISMATCH` between child and parent (error), `EXON_OVERLAP` (exons of one transcript overlap; error), `CDS_OUTSIDE_EXON` (a CDS interval not contained in one exon; error), `CDS_PHASE_MISMATCH` (the phase differs from the one implied by the lengths of the preceding CDS parts in transcript order, `(3 - preceding_length mod 3) mod 3` relative to the first phase; error), `CDS_LENGTH_NOT_MULTIPLE_OF_3` (CDS plus stop codon length; warning), `GENE_WITHOUT_TRANSCRIPT` and `TRANSCRIPT_WITHOUT_EXON` (warning), and `BEYOND_SEQUENCE_END` (needs `contig_lengths`; error). It checks internal consistency only, not whether the annotation is biologically correct. Syntax is checked by `validate`.

- **API:** `dnakit.annotation.check_annotation_consistency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`). `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`), `contig_lengths` (`{seqid: length}`, enables `BEYOND_SEQUENCE_END`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `consistent` (no error), `issues` (each `{severity, code, line, message}`), `error_count`, `warning_count` and `checked` (`format`, `features`, `genes`, `transcripts`). In the example the second CDS carries phase 0 instead of the implied 2.

- **Example code:**

```python
import json
from dnakit.annotation import check_annotation_consistency

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t0\tID=c2;Parent=t1\n"
)
result = check_annotation_consistency(gff3)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "check_annotation_consistency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "consistent": false,
    "issues": [
      {
        "severity": "error",
        "code": "CDS_PHASE_MISMATCH",
        "line": 7,
        "message": "transcript 't1': CDS 2000-2199 has phase 0, expected 2"
      }
    ],
    "error_count": 1,
    "warning_count": 0,
    "checked": {
      "format": "gff3",
      "features": 6,
      "genes": 1,
      "transcripts": 1
    }
  }
}
```

## 7) compare

- **Purpose:** Compare two annotations by feature coordinates `(seqid, type, start, end, strand)`.

- **Method:** Each annotation becomes a multiset of `(seqid, type, start, end, strand)` keys; `shared`, `only_left` and `only_right` count the distinct keys in the intersection and in each difference, and `jaccard` is `shared / union` (1.0 for two empty sets). Attributes, ids and scores are ignored, so renamed but identically placed features count as shared. `by_type` splits the counts per feature type.

- **API:** `dnakit.annotation.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the first annotation (GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`)); the second is `options["other"]`. A pair `[left, right]` or a mapping `{"left", "right"}` is also accepted; a missing second annotation raises `ConfigurationError` (`OTHER_REQUIRED`). `options`: `other`, `feature_types` (restrict to these types, case-insensitive), `max_examples` (default 20), `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `left_count`, `right_count`, `shared`, `only_left`, `only_right`, `jaccard`, `by_type` and the sample lists `only_left_examples` / `only_right_examples` (`seqid`, `type`, `start`, `end`, `strand`).

- **Example code:**

```python
import json
from dnakit.annotation import compare

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
other = gff3.replace("1000\t1200", "1000\t1250", 1)
result = compare(gff3, options={"other": other, "feature_types": ["exon", "CDS"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "left_count": 4,
    "right_count": 4,
    "shared": 3,
    "only_left": 1,
    "only_right": 1,
    "jaccard": 0.6,
    "by_type": {
      "CDS": {
        "shared": 2,
        "only_left": 0,
        "only_right": 0
      },
      "exon": {
        "shared": 1,
        "only_left": 1,
        "only_right": 1
      }
    },
    "only_left_examples": [
      {
        "seqid": "chr1",
        "type": "exon",
        "start": 1000,
        "end": 1200,
        "strand": "+"
      }
    ],
    "only_right_examples": [
      {
        "seqid": "chr1",
        "type": "exon",
        "start": 1000,
        "end": 1250,
        "strand": "+"
      }
    ]
  }
}
```

## 8) convert_annotation_coordinates

- **Purpose:** Document a reusable convert annotation coordinates capability for Annotation workflows. Candidate CAP-0042 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert annotation coordinates operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.convert_annotation_coordinates(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `convert_annotation_coordinates` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import convert_annotation_coordinates

result = convert_annotation_coordinates(
    input={"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
    options={"coordinate_system": "1-based-closed", "target_coordinate_system": "0-based-half-open"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0042",
  "operation": "convert_annotation_coordinates",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'convert_annotation_coordinates', 'result': {'records': [{'chrom': 'chr1', 'start': 1, 'end': 5, 'id': 'g1', 'coordinate_system': '0-based-half-open'}], 'record_count': 1, 'source_coordinate_system': '1-based-closed', 'target_coordinate_system': '0-based-half-...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) export

- **Purpose:** Convert an annotation to GFF3, GTF, BED6, BED12 or JSON text, optionally writing a file.

- **Method:** GFF3 to GTF resolves `gene_id`/`transcript_id` through `Parent` links; GTF to GFF3 creates the missing gene and transcript lines from their children (`synthesize_parents`). BED uses 0-based half-open starts (`start - 1`); BED6 exports the features of one type; BED12 has one line per transcript (blocks are the exons, the thick region is the CDS plus start/stop codons, empty when non-coding). The file is written atomically and only when `output` is given.

- **API:** `dnakit.annotation.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`). `options`: `format` (target: `"gff3"` default, `"gtf"`, `"bed"`, `"bed12"`, `"json"`), `source_format` (format of the input: `"auto"`, `"gff3"`, `"gtf"`), `feature_type` (BED6 only, default `"gene"`), `synthesize_parents` (GTF to GFF3, default true), `output` (file path; nothing is written without it), `overwrite` (default false; an existing file raises `ConfigurationError`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `text`, `line_count` (data lines, without comments) and `path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.annotation import export

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = export(gff3, options={"format": "bed12"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "bed12",
    "text": "chr1\t999\t2300\tt1\t0\t+\t1100\t2199\t0\t2\t201,301,\t0,1000,\n",
    "line_count": 1,
    "path": null
  }
}
```

## 10) extract

- **Purpose:** Extract the features of one biological entity from GFF3/GTF annotation, or scan a DNA sequence for open reading frames.

- **Method:** `entity` (default `gene`) is one of `gene`, `transcript`, `mrna`, `exon`, `cds`, `start_codon`, `stop_codon`, `mirna` (miRNA/pre_miRNA types, or genes/transcripts with a miRNA biotype), `intron` (derived as the gaps between consecutive exons and numbered in transcript order), `utr`/`five_prime_utr`/`three_prime_utr` (annotated, or derived from the exons and the CDS span as in `build_gene_model`), any other value (matched case-insensitively against the feature type) or `orf`. With `genome`, the genomic sequence of each feature is added (reverse-complemented on the minus strand). For `orf` the sequence is scanned in all frames (both strands unless `both_strands` is false): an ORF runs from the first start codon after a stop codon to the next in-frame stop codon, using the standard genetic code; coordinates are 1-based on the forward strand. ORFs are reported by this definition only, with no gene-likelihood filtering beyond `min_length`.

- **API:** `dnakit.annotation.extract(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: for annotation entities GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`); for `orf` a DNA sequence, FASTA text, a mapping `{id: sequence}` or a list. `options`: `entity`, `seqid`, `strand`, `start`/`end` (keep features overlapping that range), `limit` (default 100000), `genome` (FASTA text/path or `{seqid: sequence}`), `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`); for `orf`: `min_length` (nucleotides including the stop codon, default 90), `start_codons` (default `["ATG"]`), `stop_codons` (default `["TAA", "TAG", "TGA"]`), `require_stop` (default true), `both_strands` (default true).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `entity`, `format`, `count`, `truncated` and `features` (feature records with `seqid`, `type`, `start`, `end`, `strand`, `phase`, `attributes`, `derived`, and `sequence` when `genome` is given); for `orf` `result` has `entity`, `count` and `features` (the ORFs: `seqid`, `strand`, `frame` (+1..+3 / -1..-3), `start`, `end`, `length` in nucleotides including the stop, `has_stop`, `protein` and `protein_length`).

- **Example code:**

```python
import json
from dnakit.annotation import extract

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = extract(gff3, options={"entity": "intron"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "extract",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "intron",
    "format": "gff3",
    "count": 1,
    "truncated": false,
    "features": [
      {
        "line": null,
        "seqid": "chr1",
        "source": "dnakit",
        "type": "intron",
        "start": 1201,
        "end": 1999,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "Parent": "t1",
          "intron_number": 1
        },
        "derived": true
      }
    ]
  }
}
```

## 11) find_homologous_genes

- **Purpose:** Not implemented: homology search needs an orthology database or a sequence-similarity search program.

- **Method:** No computation is performed; neither reference data (Ensembl Compara, OrthoDB) nor BLAST/DIAMOND with reference proteomes is bundled.

- **API:** `dnakit.annotation.find_homologous_genes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarized in `input_summary`); `options` is unused.

- **Output:** Always returns the envelope with `status: "not_implemented"`, `required_dependency` ("an orthology database (Ensembl Compara, OrthoDB) or BLAST/DIAMOND with reference proteomes") and a `message` (Finding homologous genes needs reference data or a similarity search; neither is bundled); no `result`.

- **Example code:**

```python
import json
from dnakit.annotation import find_homologous_genes

result = find_homologous_genes("BRCA1", options={"species": "mouse"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "find_homologous_genes",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an orthology database (Ensembl Compara, OrthoDB) or BLAST/DIAMOND with reference proteomes",
  "message": "Finding homologous genes needs reference data or a similarity search; neither is bundled."
}
```

## 12) index

- **Purpose:** Build lookup tables for a GFF3/GTF annotation: by id, by sequence (sorted intervals) and parent-to-children.

- **Method:** One pass over the features. `by_id` keeps the first feature of each id; `by_seqid` lists `[start, end, id_or_null, type]` sorted by start then end; `children` maps each parent id to the ids of the features that list it as parent, in file order (features without an id, such as exon/CDS lines, appear as `"<type>@<line>"`); `gene_transcripts` maps gene ids to transcript ids. These are plain dictionaries, not an interval tree.

- **API:** `dnakit.annotation.index(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`). `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `feature_count`, `by_id`, `by_seqid`, `children`, `gene_transcripts` and `type_counts`.

- **Example code:**

```python
import json
from dnakit.annotation import index

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
)
result = index(gff3)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "index",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 2,
    "by_id": {
      "g1": {
        "seqid": "chr1",
        "type": "gene",
        "start": 1000,
        "end": 5000,
        "strand": "+"
      },
      "t1": {
        "seqid": "chr1",
        "type": "mRNA",
        "start": 1000,
        "end": 5000,
        "strand": "+"
      }
    },
    "by_seqid": {
      "chr1": [
        [
          1000,
          5000,
          "g1",
          "gene"
        ],
        [
          1000,
          5000,
          "t1",
          "mRNA"
        ]
      ]
    },
    "children": {
      "g1": [
        "t1"
      ]
    },
    "gene_transcripts": {
      "g1": [
        "t1"
      ]
    },
    "type_counts": {
      "gene": 1,
      "mRNA": 1
    }
  }
}
```

## 13) make_intergenic_regions

- **Purpose:** Document a reusable make intergenic regions capability for Annotation workflows. Candidate CAP-0047 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the make intergenic regions operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.make_intergenic_regions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `make_intergenic_regions` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import make_intergenic_regions

result = make_intergenic_regions(
    input={"genes": [{"chrom": "chr1", "start": 2, "end": 4}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0047",
  "operation": "make_intergenic_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_intergenic_regions', 'result': {'regions': [{'id': 'intergenic_chr1_1', 'chrom': 'chr1', 'start': 1, 'end': 1, 'feature_type': 'intergenic', 'coordinate_system': '1-based-closed'}, {'id': 'intergenic_chr1_2', 'chrom': 'chr1', 'start': 5, 'end': 10, 'feat...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) make_promoter_regions

- **Purpose:** Document a reusable make promoter regions capability for Annotation workflows. Candidate CAP-0048 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the make promoter regions operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.make_promoter_regions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `make_promoter_regions` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import make_promoter_regions

result = make_promoter_regions(
    input={"genes": [{"id": "g1", "chrom": "chr1", "start": 3, "end": 6, "strand": "+"}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0048",
  "operation": "make_promoter_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_promoter_regions', 'result': {'promoters': [{'id': 'promoter_g1', 'gene_id': 'g1', 'chrom': 'chr1', 'start': 1, 'end': 3, 'strand': '+', 'feature_type': 'promoter', 'coordinate_system': '1-based-closed', 'upstream': 2, 'downstream': 1}], 'promoter_count'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) merge_gene_models

- **Purpose:** Document a reusable merge gene models capability for Annotation workflows. Candidate CAP-0049 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge gene models operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.merge_gene_models(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `merge_gene_models` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import merge_gene_models

result = merge_gene_models(
    input={"segments": [{"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3}, {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5}]} ,
    options={"coordinate_system": "1-based-closed", "merge_adjacent": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0049",
  "operation": "merge_gene_models",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'merge_gene_models', 'result': {'models': [{'gene_id': 'g1', 'gene_name': 'g1', 'chrom': 'chr1', 'strand': '+', 'start': 1, 'end': 5, 'segments': [{'start': 1, 'end': 5}], 'segment_count': 1, 'coordinate_system': '1-based-closed'}], 'model_count': 1, 'input_se...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) parse

- **Purpose:** Parse GFF3 or GTF text (or a `.gff`/`.gtf`/`.gz` file) into feature records.

- **Method:** The format is detected from a `##gff-version` line or the attribute style unless `format` is given. GFF3 attributes are percent-decoded and comma lists become lists; GTF attributes are `key "value";` pairs with repeated keys becoming lists. Coordinates stay 1-based and closed. With `strict` (default) a malformed line raises `InputFormatError` (`INVALID_ANNOTATION`); with `strict: false` it is skipped and listed under `skipped`. A list of feature mappings is normalized instead of parsed.

- **API:** `dnakit.annotation.parse(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the annotation text, a file path (`.gz` allowed), or a list of feature mappings. `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`), `strict` (default true), `limit` (maximum features returned, default 100000).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `feature_count`, `truncated`, `features` (each `line`, `seqid`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`), `directives` (the `##` lines) and `skipped`.

- **Example code:**

```python
import json
from dnakit.annotation import parse

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
)
result = parse(gff3)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "parse",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 2,
    "truncated": false,
    "features": [
      {
        "line": 2,
        "seqid": "chr1",
        "source": "src",
        "type": "gene",
        "start": 1000,
        "end": 5000,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "ID": "g1",
          "Name": "GENE1"
        }
      },
      {
        "line": 3,
        "seqid": "chr1",
        "source": "src",
        "type": "mRNA",
        "start": 1000,
        "end": 5000,
        "score": null,
        "strand": "+",
        "phase": null,
        "attributes": {
          "ID": "t1",
          "Parent": "g1"
        }
      }
    ],
    "directives": [
      "##gff-version 3"
    ],
    "skipped": []
  }
}
```

## 17) summarize

- **Purpose:** Summarize a GFF3/GTF annotation: feature counts and lengths per type, per sequence and strand, and gene-structure statistics.

- **Method:** Lengths are `end - start + 1`. Gene structure comes from the same model as `build_gene_model`: `intron_count` is derived from the gaps between exons, `coding_transcripts` counts transcripts with CDS, `multi_exon_transcripts` those with more than one exon, and `transcripts_per_gene` / `exons_per_transcript` give mean and maximum (0 when there are no values).

- **API:** `dnakit.annotation.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: GFF3 or GTF text, the path of a GFF3/GTF file (`.gz` allowed), or a list of feature mappings (as returned by `parse`; each needs `seqid`, integer `start` <= `end` >= 1, and optionally `type`, `strand`, `phase`, `attributes`). `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `feature_count`, `types` (type -> `count`, `total_bp`, `mean_length`), `seqids`, `strands`, `gene_count`, `transcript_count`, `exon_count`, `intron_count`, `coding_transcripts`, `multi_exon_transcripts`, `transcripts_per_gene` and `exons_per_transcript` (`mean`, `max`).

- **Example code:**

```python
import json
from dnakit.annotation import summarize

gff3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tsrc\tCDS\t2000\t2199\t.\t+\t2\tID=c2;Parent=t1\n"
)
result = summarize(gff3)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "summarize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "gff3",
    "feature_count": 6,
    "types": {
      "CDS": {
        "count": 2,
        "total_bp": 300,
        "mean_length": 150.0
      },
      "exon": {
        "count": 2,
        "total_bp": 502,
        "mean_length": 251.0
      },
      "gene": {
        "count": 1,
        "total_bp": 4001,
        "mean_length": 4001.0
      },
      "mRNA": {
        "count": 1,
        "total_bp": 4001,
        "mean_length": 4001.0
      }
    },
    "seqids": {
      "chr1": 6
    },
    "strands": {
      "+": 6
    },
    "gene_count": 1,
    "transcript_count": 1,
    "exon_count": 2,
    "intron_count": 1,
    "coding_transcripts": 1,
    "multi_exon_transcripts": 1,
    "transcripts_per_gene": {
      "mean": 1.0,
      "max": 1
    },
    "exons_per_transcript": {
      "mean": 2.0,
      "max": 2
    }
  }
}
```

## 18) validate

- **Purpose:** Validate the syntax of a GFF3 or GTF file: column count, coordinates, strand, phase, attributes, ids and parents.

- **Method:** Errors: wrong column count, non-integer or non-positive coordinates, end < start, strand not in `+ - . ?`, phase not in 0/1/2/`.`, CDS without phase, unparsable attributes; GFF3 adds duplicate `ID` of non-CDS features and `Parent` pointing to an undefined `ID`; GTF adds a missing `gene_id` (all lines) or `transcript_id` (all but gene lines). Warning: GFF3 without a `##gff-version 3` directive. Malformed lines are reported, not raised. Logical checks of gene structures are in `check_annotation_consistency`.

- **API:** `dnakit.annotation.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: annotation text or a file path (`.gz` allowed). `options`: `format` (`"auto"` default: a `##gff-version` line or the attribute style decides; `"gff3"`, `"gtf"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `valid` (no error), `format`, `feature_count` (lines that parsed), `errors` and `warnings` (each `{line, severity, code, message}`, sorted by line).

- **Example code:**

```python
import json
from dnakit.annotation import validate

bad = "chr1\tsrc\tgene\t50\t10\t.\t+\t.\tID=g1\nchr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=nope\n"
result = validate(bad)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "annotation",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "format": "gff3",
    "feature_count": 1,
    "errors": [
      {
        "line": 1,
        "severity": "error",
        "code": "COORDINATE_ORDER",
        "message": "end (10) is smaller than start (50)"
      },
      {
        "line": 2,
        "severity": "error",
        "code": "UNDEFINED_PARENT",
        "message": "Parent 'nope' is not defined"
      }
    ],
    "warnings": [
      {
        "line": 1,
        "severity": "warning",
        "code": "MISSING_VERSION",
        "message": "no '##gff-version 3' directive"
      }
    ]
  }
}
```
