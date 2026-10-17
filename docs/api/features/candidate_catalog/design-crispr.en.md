# CRISPR design

This category contains **19** documented candidate capabilities; 18 compute real results natively (using the standard library and, for some, `dnakit.patterns` and `dnakit.molbio`) and 1 (`predict_editing_outcomes`) always returns `status: "not_implemented"` because it needs trained repair-outcome model weights. Several steps are explicit heuristics rather than validated predictors: the guide `proxy_score` is a rule-based proxy that was not fitted to editing data, off-target search is an ungapped Hamming search without CFD/MIT scoring, and the primer checks count contiguous Watson-Crick pairs without free-energy models. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) assemble_dna_fragments

- **Purpose:** Assemble DNA fragments in the given order by exact terminal overlaps (in-silico Gibson/overlap joining).

- **Method:** Delegates to `dnakit.molbio.simulate_assembly(method="gibson")`: each fragment is joined to the next at the longest exact suffix/prefix overlap whose length lies in `[min_overlap, max_overlap]`; with `circularize` the last fragment is also joined back to the first. Only exact sequence overlaps are used; reaction kinetics, secondary structure and mis-annealing are not modelled. If two consecutive fragments share no admissible overlap the call raises a DNAKit error instead of returning a product.

- **API:** `dnakit.design_crispr.assemble_dna_fragments(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: ordered fragments as a list of DNA strings, a mapping `{id: sequence}` (insertion order is the assembly order) or FASTA text/path; at least two fragments (at most 1000), characters limited to ACGT (U is read as T). `options`: `min_overlap` (default 15), `max_overlap` (default 200), `circularize` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `product` (assembled sequence), `length`, `circular` and `junctions` (one entry per junction with `left`, `right`, `overlap` and `overlap_length`).

- **Example code:**

```python
import json
from dnakit.design_crispr import assemble_dna_fragments

overlap = "ACGTACGTACGTACGTAAAA"
fragments = {"left": "AAAAAC" + overlap, "right": overlap + "CCCCCT"}
result = assemble_dna_fragments(fragments, options={"min_overlap": 15})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "assemble_dna_fragments",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 77,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c9c510fcdf0c197ea95010181d5c1150e5924cb4cb3ea553a5a3e7054f95e81f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "product": "AAAAACACGTACGTACGTACGTAAAACCCCCT",
    "length": 32,
    "circular": false,
    "junctions": [
      {
        "left": "left",
        "right": "right",
        "overlap": "ACGTACGTACGTACGTAAAA",
        "overlap_length": 20
      }
    ]
  }
}
```

## 2) calculate_primer_tm

- **Purpose:** Melting temperature (Tm) of one or more primers, by the Wallace rule for short oligos or the SantaLucia nearest-neighbour model, together with each primer's length and GC fraction.

- **Method:** Each primer is upper-cased and stripped of whitespace; GC fraction = (G + C) / length. `method` `"wallace"`: Tm = 2 x (A + T) + 4 x (G + C) in degrees Celsius, a rule of thumb for oligos of about 14 nt or fewer that ignores salt, concentration and sequence context and overestimates Tm for longer oligos (a warning is added for a Wallace Tm of a primer longer than 14 nt). `"nearest_neighbor"`: `dnakit.thermodynamics.melting_temperature` (SantaLucia 1998 unified DNA/DNA parameters, perfectly matched primer/complement duplex, 2-60 nt) at the given `sodium_molar`, `potassium_molar` and `strand_concentration_molar`, rounded to 0.01 degrees; Mg2+, dNTP and DMSO corrections are not modelled, so the value in a Mg2+ PCR buffer is typically higher. `"auto"` (default): Wallace for length <= 13, nearest-neighbour above. Ambiguous bases (N, IUPAC codes) are rejected because Tm is undefined for them; `U` is rejected unless `rna` is true, in which case it is read as `T`.

- **API:** `dnakit.design_crispr.calculate_primer_tm(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a primer string (5'->3'), `{"id": ..., "sequence": ...}`, a mapping `{name: primer}`, or a list of strings / `{"id", "sequence"}` mappings; 1 to 1000 primers of at most 200 bases (a string with several lines is one primer: whitespace is removed). `options`: `method` (`"auto"` default | `"wallace"` | `"nearest_neighbor"`), `sodium_molar` (default 0.05), `potassium_molar` (default 0.0), `strand_concentration_molar` (default 250e-9), all in mol/L and only used by the nearest-neighbour Tm, and `rna` (default `false`).

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"design_crispr"`), `operation` (`"calculate_primer_tm"`), `provenance` (`implementation`, `backend`), `input_summary` (generic text summary: `record_count`, `text_length`, `sequence_bases`, `lengths`, `sha256`; computed from the input as text), `warnings` and `result`. `result.primers` has one entry per primer in input order: `id` (`primer1`, `primer2`, ... for unnamed ones), `sequence` (cleaned), `length`, `gc_fraction`, `tm_celsius` (a float) and `method` (`"wallace"` or `"nearest_neighbor"`, the method actually used). `result.conditions` (`sodium_molar`, `potassium_molar`, `strand_concentration_molar`) is present when any primer used the nearest-neighbour model. An empty or ambiguous sequence raises `InputFormatError` (`CRISPR_SEQUENCE_INVALID`); a nearest-neighbour request for a primer outside 2-60 nt raises `ConfigurationError` (`CRISPR_PRIMER_LENGTH`).

- **Example code:**

```python
import json
from dnakit.design_crispr import calculate_primer_tm

result = calculate_primer_tm({"short": "ATGCGTAC", "long": "GACGTTAGCATCAGGCTAAC"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "calculate_primer_tm",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d1cef6d866bd036122f644676770ae56c33d025666a8977ff9212d2c0653bfd5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "id": "short",
        "sequence": "ATGCGTAC",
        "length": 8,
        "gc_fraction": 0.5,
        "tm_celsius": 24.0,
        "method": "wallace"
      },
      {
        "id": "long",
        "sequence": "GACGTTAGCATCAGGCTAAC",
        "length": 20,
        "gc_fraction": 0.5,
        "tm_celsius": 52.86,
        "method": "nearest_neighbor"
      }
    ],
    "conditions": {
      "sodium_molar": 0.05,
      "potassium_molar": 0.0,
      "strand_concentration_molar": 2.5e-07
    }
  }
}
```

## 3) check_primer_dimer

- **Purpose:** Screen primers for self- and cross-dimers by the longest perfectly complementary ungapped stretch.

- **Method:** Every primer is compared with itself and with every other primer. For primers A and B, `longest_complementary_run` is the length of the longest common substring of A and the reverse complement of B (dynamic programming, O(len(A)*len(B))); `three_prime_run_a`/`three_prime_run_b` is the longest suffix (3' end) of one primer that is reverse-complementary to some part of the other. A pair is flagged `three_prime_dimer` when either 3' run exceeds `max_3prime_run` and `long_complementary_run` when the longest stretch exceeds `max_run`. This counts contiguous Watson-Crick pairs only; it is not a free-energy (nearest-neighbour) calculation and the thresholds are rules of thumb.

- **API:** `dnakit.design_crispr.check_primer_dimer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one primer string (self-dimer only), a list of primers (at most 100) or a mapping such as `{"forward": ..., "reverse": ...}`; primers are ACGT only, up to 200 nt. `options`: `max_3prime_run` (default 3), `max_run` (default 6).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.pairs` has one entry per unordered primer pair including self pairs: `a`, `b`, `self`, `longest_complementary_run`, `run_start_in_a` (0-based, `null` if none), `three_prime_run_a`, `three_prime_run_b`, `three_prime_dimer`, `long_complementary_run`.

- **Example code:**

```python
import json
from dnakit.design_crispr import check_primer_dimer

result = check_primer_dimer({"forward": "AAAAGGGG", "reverse": "CCCCTTTT"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_dimer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 46,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5d740a084f018177281992c5b4d9cb8ee12cee423c0c05d77fe37c9bf5920439"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": [
      {
        "a": "forward",
        "b": "forward",
        "self": true,
        "longest_complementary_run": 0,
        "run_start_in_a": null,
        "three_prime_run_a": 0,
        "three_prime_run_b": 0,
        "three_prime_dimer": false,
        "long_complementary_run": false
      },
      {
        "a": "forward",
        "b": "reverse",
        "self": false,
        "longest_complementary_run": 8,
        "run_start_in_a": 0,
        "three_prime_run_a": 8,
        "three_prime_run_b": 8,
        "three_prime_dimer": true,
        "long_complementary_run": true
      },
      {
        "a": "reverse",
        "b": "reverse",
        "self": true,
        "longest_complementary_run": 0,
        "run_start_in_a": null,
        "three_prime_run_a": 0,
        "three_prime_run_b": 0,
        "three_prime_dimer": false,
        "long_complementary_run": false
      }
    ]
  }
}
```

## 4) check_primer_secondary_structure

- **Purpose:** Screen primers for hairpins: maximal perfect Watson-Crick stems closed by a loop of at least `min_loop` bases.

- **Method:** All position pairs (i, j) are tested: a stem is a run of nested A-T/C-G pairs (seq[i+k], seq[j-k]); only maximal stems are reported (the outer pair cannot be extended outwards, and the stem is extended inwards while the loop stays >= `min_loop`). Stems shorter than `min_stem` are dropped and hairpins are sorted by stem length. There are no mismatches, bulges or G-U pairs and no energy model, so this is a stem-count screen and not a thermodynamic fold (it is not NUPACK or mfold).

- **API:** `dnakit.design_crispr.check_primer_secondary_structure(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one primer string, a list of primers (at most 1000) or a mapping `{name: primer}`; ACGT only, up to 200 nt. `options`: `min_stem` (default 4, range 2-50) minimum stem length in base pairs; `min_loop` (default 3, range 0-50) minimum unpaired loop length.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.primers` has one entry per primer: `primer`, `sequence`, `hairpins` (each with `stem_length`, `left_start`, `right_start`, `loop_length`, `stem_5prime`, `stem_3prime`, `involves_3prime_end`), `longest_stem` and `three_prime_hairpin`.

- **Example code:**

```python
import json
from dnakit.design_crispr import check_primer_secondary_structure

result = check_primer_secondary_structure({"p1": "GCGCAAAAGCGC", "p2": "ACACACACACAC"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_secondary_structure",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 44,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "82ba75d2df38faa9ec3b57f89e3989391a1df3dfe43ef1202550a994c32f5885"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "primer": "p1",
        "sequence": "GCGCAAAAGCGC",
        "hairpins": [
          {
            "stem_length": 4,
            "left_start": 0,
            "right_start": 8,
            "loop_length": 4,
            "stem_5prime": "GCGC",
            "stem_3prime": "GCGC",
            "involves_3prime_end": true
          }
        ],
        "longest_stem": 4,
        "three_prime_hairpin": true
      },
      {
        "primer": "p2",
        "sequence": "ACACACACACAC",
        "hairpins": [],
        "longest_stem": 0,
        "three_prime_hairpin": false
      }
    ]
  }
}
```

## 5) check_primer_specificity

- **Purpose:** Count where a primer binds (ungapped, mismatches allowed) on each template sequence, on both strands.

- **Method:** Delegates to `dnakit.molbio.match_primer` for every template: an ungapped search of the primer on both strands that allows up to `max_mismatches` mismatches and none within the last `strict_three_prime_bases` bases of the 3' end. All hits of all templates are pooled; `specific` is true when exactly one site is found in total. Binding energy and extension competence are not evaluated, so a hit is a sequence match and not a prediction of amplification.

- **API:** `dnakit.design_crispr.check_primer_specificity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: one primer string, a list of primers or a mapping `{name: primer}` (ACGT, up to 200 nt). `options`: `templates` (required: mapping `{id: sequence}`, list of sequences, or FASTA text/path; total at most 20 Mb), `max_mismatches` (default 0, range 0-10), `strict_three_prime_bases` (default 0, range 0-50). Without `templates` the call raises a configuration error.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result.primers` has one entry per primer: `primer`, `sequence`, `total_hits`, `specific` and `hits` (each with `template`, `strand`, 0-based half-open forward `start`/`end`, `mismatch_positions` from the 5' end).

- **Example code:**

```python
import json
from dnakit.design_crispr import check_primer_specificity

primer = "GAATTCGGGACG"
template = "TTT" + primer + "AAAA" + "CGTCCCGAATTC" + "CC"   # second site is on the reverse strand
result = check_primer_specificity(primer, options={"templates": {"t1": template}})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "check_primer_specificity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 12,
    "sequence_bases": 12,
    "lengths": [
      12
    ],
    "sha256": "7994e92a211084996bdf46d96ef74833b03dd460002c7633153f9953445df383"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "primers": [
      {
        "primer": "primer1",
        "sequence": "GAATTCGGGACG",
        "total_hits": 2,
        "specific": false,
        "hits": [
          {
            "template": "t1",
            "strand": "forward",
            "start": 3,
            "end": 15,
            "mismatch_positions": []
          },
          {
            "template": "t1",
            "strand": "reverse",
            "start": 19,
            "end": 31,
            "mismatch_positions": []
          }
        ]
      }
    ]
  }
}
```

## 6) design

- **Purpose:** Design CRISPR guide RNAs for a target sequence, or an HDR donor template around a cut site (`options["entity"]` = `"grna"` or `"donor_template"`).

- **Method:** `grna`: guides are enumerated next to a PAM on both strands with `dnakit.patterns.scan_pam_candidates` (SpCas9 `NGG` 3' PAM with 20 nt guide, SaCas9 `NNGRRT` with 21 nt, AsCas12a `TTTV` 5' PAM with 20 nt, or a custom PAM), then filtered by GC range, homopolymer length, `TTTT` (RNA Pol III terminator) and excluded motifs/restriction sites (both strands), and sorted by `proxy_score`, then position. `proxy_score` is a rule-based heuristic starting at 100 and subtracting 15 for GC outside 40-70 % (40 instead outside 20-80 %), 40 for a `TTTT` run and 15 for an A/C/G homopolymer of length >= 4; it was not fitted to editing data and is not a predicted efficiency. With `references`, each guide gets a Hamming off-target count (ungapped, both strands, PAM required): `exact_sites` (perfect matches), `near_sites` (1..`max_mismatches` mismatches) and `off_target_sites` = near_sites + max(exact_sites - 1, 0). `cut_site` is the boundary 3 bp from the PAM for SpCas9/SaCas9 and `null` otherwise. `donor_template`: the donor is left homology arm + `insert` + right homology arm, with `delete_length` bases after the cut removed; no HDR efficiency or PAM-blocking mutation is considered.

- **API:** `dnakit.design_crispr.design(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the target DNA (string, FASTA text/path or `{"sequence": ...}`; one sequence, ACGT). `options`: `entity` (`"grna"` default | `"donor_template"`). For `grna`: `nuclease` (`"SpCas9"` default, `"SaCas9"`, `"AsCas12a"`) or custom `pam`/`pam_side` (`"3prime"` | `"5prime"`)/`guide_length`; `strand` (`"both"` default | `"forward"` | `"reverse"`); `min_gc` 0.4, `max_gc` 0.7, `max_homopolymer` 4 (longer runs rejected), `avoid_poly_t` `true`, `exclude_motifs` (IUPAC list), `avoid_enzymes` (built-in restriction enzyme names), `max_guides` 100, optional `references` (mapping/FASTA) and `max_mismatches` (default 3). For `donor_template`: `cut_site` (required, 0-based boundary), `insert` (default empty), `delete_length` (default 0), `left_arm`/`right_arm` (default 60).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. For `grna`, `result` has `target`, `nuclease`, `candidates_scanned` and `guides` (each with `id`, `sequence`, `pam`, `strand`, 0-based half-open forward `start`/`end`, `pam_start`, `cut_site`, `gc`, `longest_homopolymer`, `proxy_score`, and with `references` also `exact_sites`, `near_sites`, `off_target_sites`); `warnings` states that `proxy_score` is a heuristic. For `donor_template`, `result` has `target`, `donor`, `donor_length`, `left_arm`/`right_arm` (`start`, `end`, `sequence`), `insert`, `deleted` and `edited_sequence`. Homology arms that leave the target raise a configuration error.

- **Example code:**

```python
import json
from dnakit.design_crispr import design

target = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
result = design(target, options={"nuclease": "SpCas9", "max_guides": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 46,
    "sequence_bases": 46,
    "lengths": [
      46
    ],
    "sha256": "899ca71f07ea07f2bf5addd556b05562c9a7750d20636de9e9f49275ea9b61f9"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "target": "seq1",
    "nuclease": "SpCas9",
    "candidates_scanned": 4,
    "guides": [
      {
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "pam": "AGG",
        "strand": "forward",
        "start": 0,
        "end": 20,
        "pam_start": 20,
        "cut_site": 17,
        "gc": 0.5,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g1"
      },
      {
        "sequence": "AGCTTCAAGGCCAGCTAGCT",
        "pam": "AGG",
        "strand": "forward",
        "start": 13,
        "end": 33,
        "pam_start": 33,
        "cut_site": 30,
        "gc": 0.55,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g2"
      },
      {
        "sequence": "GGCCAGCTAGCTAGGCTAAT",
        "pam": "CGG",
        "strand": "forward",
        "start": 21,
        "end": 41,
        "pam_start": 41,
        "cut_site": 38,
        "gc": 0.55,
        "longest_homopolymer": 2,
        "proxy_score": 100.0,
        "id": "seq1_g3"
      }
    ]
  }
}
```

## 7) encode_data_into_dna

- **Purpose:** Encode bytes as a DNA string with a fixed 2-bits-per-base code (4 bases per byte). Raw bytes and files are encoded as they are. This is a minimal toy mapping, not a synthesis-ready DNA-storage scheme.

- **Method:** The input is turned into bytes (see Input) and every byte is written as four bases, most significant bit pair first, with A = 00, C = 01, G = 10, T = 11 (`"Hi!"` = bytes 0x48 0x69 0x21 gives `CAGA CGGC AGAC`; bytes 0x00 0x48 0xFF 0xB2 give `AAAA CAGA TTTT GTAG`). The code is lossless for every byte value 0-255: `dnakit.design_crispr._impl_local.dna_to_bytes` is the exact inverse (it needs ACGT and a length that is a multiple of 4), and the same mapping is what the shared runtime's `decode_dna_storage` assumes (that runtime branch is not exposed by any public module and decodes to UTF-8 text only; use `dna_to_bytes` for binary data). The code has no index, header, checksum, error correction, primers or constraints on GC content or homopolymer runs, so the output is not suitable for synthesis as is. `dnakit.genomics.encode_dna_storage` and `decode_dna_storage` provide indexed segments with an integrity digest (they need the optional `dnakit[dna-storage]` extra).

- **API:** `dnakit.design_crispr.encode_data_into_dna(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `bytes`/`bytearray`/`memoryview` (encoded as they are, any byte values), a `Path` of an existing file (its bytes), a `str` (its UTF-8 bytes, so `"é"` gives 8 bases; with `options["source"] = "path"` the string is read as a file path instead), or a mapping/list (sorted-key JSON text, `ensure_ascii=False`, UTF-8). At most 1,000,000 bytes (`ConfigurationError`, `CRISPR_INPUT_TOO_LARGE`); a missing file, an unsupported type (such as an `int`) or non-JSON mapping content raises `InputFormatError`. `options`: `source` (`"auto"` default | `"text"` | `"path"`); with `"auto"` and `"text"` a string is text, never a path.

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"design_crispr"`), `operation` (`"encode_data_into_dna"`), `provenance` (`implementation`, `backend`), `input_summary` (generic text summary: `record_count`, `text_length`, `sequence_bases`, `lengths`, `sha256`; computed from the input as text), `warnings` and `result`. `result` is a string of `A`, `C`, `G`, `T` of length 4 x (number of bytes). Note that `input_summary` describes the input as text (raw bytes appear as a UTF-8 decoding or, when not valid UTF-8, as base64 JSON).

- **Example code:**

```python
import json
from dnakit.design_crispr import encode_data_into_dna

result = encode_data_into_dna(bytes([0x00, 0x48, 0xFF, 0xB2]))   # raw bytes, encoded as they are
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "encode_data_into_dna",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 4,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d909adf4e65a93c9f454770461199a5eff63c53a49a20aeff3418ae597a85212"
  },
  "warnings": [],
  "status": "ok",
  "result": "AAAACAGATTTTGTAG"
}
```

## 8) evaluate_primer_pair

- **Purpose:** Evaluate a forward/reverse primer pair: Tm of each primer and their difference, GC fraction, 3' GC clamp, self- and cross-complementarity, and a list of flags against documented thresholds.

- **Method:** Each primer gets the Tm of `calculate_primer_tm` (same `method` for both; with `auto` a pair that straddles 13 nt gets two different methods and a warning, because `delta_tm_celsius` is then not meaningful). `delta_tm_celsius` = |Tm(forward) - Tm(reverse)|. `gc_clamp` = number of G/C among the last 5 bases of the primer (3' end). `longest_self_complementary` = longest contiguous stretch of the primer whose reverse complement also occurs in the same primer (the self-dimer case of two copies annealing antiparallel; a palindromic primer such as `GAATTC` gives its full length; hairpins are covered by `check_primer_secondary_structure`). `longest_cross_complementary` = longest contiguous stretch of the forward primer that is the reverse complement of a stretch of the reverse primer. Complementarity is exact contiguous Watson-Crick matching; no free energy, mismatches or bulges are modelled (use `check_primer_dimer` for 3'-end runs). Flags (default thresholds in brackets, all overridable): `tm_mismatch` when delta Tm > `max_delta_tm` [5.0 C]; per primer `length_out_of_range` outside `min_length`-`max_length` [18-30], `gc_out_of_range` outside `min_gc`-`max_gc` [0.40-0.60], `gc_clamp_out_of_range` outside `min_gc_clamp`-`max_gc_clamp` [1-3], `self_complementarity` when the self run > `max_self_complementary` [4]; `cross_complementarity` when the cross run > `max_cross_complementary` [4]. These are common design rules of thumb, not guarantees of PCR success.

- **API:** `dnakit.design_crispr.evaluate_primer_pair(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{"forward": ..., "reverse": ...}` (exactly these two keys) or a list of two primers `[forward, reverse]`, both as strings written 5'->3' (the reverse primer is given as ordered, not reverse-complemented); ACGT only (`U` with `rna`). Anything else raises `InputFormatError`. `options`: `method`, `rna`, `sodium_molar`, `potassium_molar`, `strand_concentration_molar` as in `calculate_primer_tm`, and the thresholds `max_delta_tm`, `min_length`, `max_length`, `min_gc`, `max_gc`, `min_gc_clamp`, `max_gc_clamp`, `max_self_complementary`, `max_cross_complementary`.

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"design_crispr"`), `operation` (`"evaluate_primer_pair"`), `provenance` (`implementation`, `backend`), `input_summary` (generic text summary: `record_count`, `text_length`, `sequence_bases`, `lengths`, `sha256`; computed from the input as text), `warnings` and `result`. `result` has `forward` and `reverse` (each `sequence`, `length`, `gc_fraction`, `tm_celsius`, `method`, `gc_clamp`, `longest_self_complementary`), `delta_tm_celsius`, `longest_cross_complementary`, `thresholds` (the values used), `flags` (a list of `{"flag", "primer", "value", "threshold"}`; `primer` is `"forward"`, `"reverse"` or `null` for pair-level flags), `pass` (true when there are no flags) and, when a nearest-neighbour Tm was computed, `conditions`.

- **Example code:**

```python
import json
from dnakit.design_crispr import evaluate_primer_pair

result = evaluate_primer_pair({"forward": "GACGTTAGCATCAGGCTAAC", "reverse": "CCATGTTCAGCTAGGATCAT"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "evaluate_primer_pair",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3fdd37d2399c743d84981ed823a66973cf44c6c25c083abdbdc8bdd2613dda0b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "forward": {
      "sequence": "GACGTTAGCATCAGGCTAAC",
      "length": 20,
      "gc_fraction": 0.5,
      "tm_celsius": 52.86,
      "method": "nearest_neighbor",
      "gc_clamp": 2,
      "longest_self_complementary": 6
    },
    "reverse": {
      "sequence": "CCATGTTCAGCTAGGATCAT",
      "length": 20,
      "gc_fraction": 0.45,
      "tm_celsius": 50.99,
      "method": "nearest_neighbor",
      "gc_clamp": 1,
      "longest_self_complementary": 4
    },
    "delta_tm_celsius": 1.87,
    "longest_cross_complementary": 4,
    "thresholds": {
      "max_delta_tm": 5.0,
      "min_length": 18,
      "max_length": 30,
      "min_gc": 0.4,
      "max_gc": 0.6,
      "min_gc_clamp": 1,
      "max_gc_clamp": 3,
      "max_self_complementary": 4,
      "max_cross_complementary": 4
    },
    "flags": [
      {
        "flag": "self_complementarity",
        "primer": "forward",
        "value": 6,
        "threshold": 4
      }
    ],
    "pass": false,
    "conditions": {
      "sodium_molar": 0.05,
      "potassium_molar": 0.0,
      "strand_concentration_molar": 2.5e-07
    }
  }
}
```

## 9) export

- **Purpose:** Serialise guide RNAs as TSV, CSV, FASTA or JSON text, optionally writing a file.

- **Method:** Guides are parsed and validated (ACGT), then written in the chosen format. TSV/CSV have columns `id`, `sequence` and then the sorted union of any extra record fields (for example the `start`, `strand`, `proxy_score` fields returned by `design`); records lacking a field get an empty cell. FASTA writes `>id` plus the sequence; JSON dumps the records with sorted keys. An existing output file is replaced only when `overwrite` is true.

- **API:** `dnakit.design_crispr.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field (as returned by `design`; a `{"guides": [...]}` mapping is accepted). `options`: `format` (`"tsv"` default | `"csv"` | `"fasta"` | `"json"`); `output` (optional file path); `overwrite` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `format`, `count`, `text` (the serialised guides) and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.design_crispr import export

guides = [{"id": "g1", "sequence": "GACGTACGTCGTAAGCTTCA", "start": 0, "strand": "forward"},
          {"id": "g2", "sequence": "ACGTACGTTTTTACGTACGT"}]
result = export(guides, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 133,
    "sequence_bases": 40,
    "lengths": [
      20,
      20
    ],
    "sha256": "b99bcb771c0c55eaedd9ee50be13be0a112cdaf706aaa2e103450482ec3eca20"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "count": 2,
    "text": "id\tsequence\tstart\tstrand\ng1\tGACGTACGTCGTAAGCTTCA\t0\tforward\ng2\tACGTACGTTTTTACGTACGT\t\t\n",
    "output_path": null
  }
}
```

## 10) filter

- **Purpose:** Filter guide RNAs by GC content, homopolymer runs, poly-T, length and excluded motifs, reporting why each rejected guide was dropped.

- **Method:** Each guide is checked against fixed design rules of thumb: GC fraction below `min_gc` (`gc_below_min`) or above `max_gc` (`gc_above_max`); longest homopolymer longer than `max_homopolymer` (`homopolymer_run`); a `TTTT` run, the RNA Pol III terminator (`poly_t`, when `avoid_poly_t`); an excluded motif or restriction site on either strand (`excluded_motif:<name>`); and a length different from `length` (`length`). These are common design heuristics, not a model of editing efficiency.

- **API:** `dnakit.design_crispr.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field. `options`: `min_gc` (default 0.4), `max_gc` (0.7), `max_homopolymer` (4: longer runs fail), `avoid_poly_t` (`true`), `exclude_motifs` (IUPAC list, both strands), `avoid_enzymes` (built-in restriction enzyme names), `length` (optional exact guide length).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `kept` and `rejected` (each record with `id`, `sequence`, `gc`, `longest_homopolymer` and the list `failed` of violated rules) plus `kept_count` and `rejected_count`.

- **Example code:**

```python
import json
from dnakit.design_crispr import filter as filter_guides

guides = {"ok": "GACGTACGTCGTAAGCTTCA", "at_rich": "ATATATATATATATATATAT",
          "poly_t": "GACTTTTCGACTGACTGACT", "ecori": "GACGAATTCGTAAGCTTCAG"}
result = filter_guides(guides, options={"avoid_enzymes": ["EcoRI"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "filter",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 132,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "02704497d75ca487b918a3590bdae63dbfa789f44f7116c23d2bba76448b6562"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kept": [
      {
        "id": "ok",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "gc": 0.5,
        "longest_homopolymer": 2,
        "failed": []
      }
    ],
    "rejected": [
      {
        "id": "at_rich",
        "sequence": "ATATATATATATATATATAT",
        "gc": 0.0,
        "longest_homopolymer": 1,
        "failed": [
          "gc_below_min"
        ]
      },
      {
        "id": "poly_t",
        "sequence": "GACTTTTCGACTGACTGACT",
        "gc": 0.45,
        "longest_homopolymer": 4,
        "failed": [
          "poly_t"
        ]
      },
      {
        "id": "ecori",
        "sequence": "GACGAATTCGTAAGCTTCAG",
        "gc": 0.45,
        "longest_homopolymer": 2,
        "failed": [
          "excluded_motif:EcoRI"
        ]
      }
    ],
    "kept_count": 1,
    "rejected_count": 3
  }
}
```

## 11) plot

- **Purpose:** Draw guide RNA positions along a target as an SVG map (forward-strand guides above the line, reverse-strand guides below).

- **Method:** A plain hand-written SVG (no plotting library): a baseline of the target length with one `<rect class="guide">` per guide, positioned linearly from `start` to `end`, blue above the line for forward guides and red below for the others. Each rectangle carries `data-start`, `data-end` and `data-strand` attributes. The picture shows positions only, no scores.

- **API:** `dnakit.design_crispr.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping `{"length": N, "guides": [...]}` or `{"sequence": ..., "guides": [...]}`; each guide needs 0 <= `start` < `end` (0-based half-open) and `strand` (`forward`/`reverse`), for example the guides in the `result` of `design`. Without a length, the largest `end` is used. `options`: `output` (optional `.svg` path).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `svg` (the SVG text) and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.design_crispr import plot

guides = [{"start": 0, "end": 20, "strand": "forward"}, {"start": 30, "end": 50, "strand": "reverse"}]
result = plot({"length": 60, "guides": guides})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "design_crispr",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 119,
    "sequence_bases": 60,
    "lengths": [
      60
    ],
    "sha256": "8f8a17b9efbb550b409b4f24f119aadff3263fc6863ad7764bffe66e5c0b8b80"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"90\"><line x1=\"20\" y1=\"45\" x2=\"780\" y2=\"45\" stroke=\"#333\" stroke-width=\"2\"/><text x=\"20\" y=\"85\" font-...",
    "output_path": null
  }
}
```

## 12) predict_editing_outcomes

- **Purpose:** Predict CRISPR repair outcomes. Not implemented: the function always returns `status: "not_implemented"`.

- **Method:** No computation is performed and no result is fabricated. Indel and repair-outcome prediction needs a trained model such as inDelphi, FORECasT or Lindel together with its weights, which DNAKit does not bundle. For HDR geometry without any outcome model use `design` (`entity="donor_template"`) and `simulate`.

- **API:** `dnakit.design_crispr.predict_editing_outcomes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: any value (it is only summarised). `options`: ignored.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. The status is `"not_implemented"` (and `provenance.implementation` is `"not_implemented"`); there is no `result`. `required_dependency` is `a trained repair-outcome model (inDelphi, FORECasT, Lindel) and its weights` and `message` is `Indel/repair outcome prediction needs trained model weights that are not bundled.`

- **Example code:**

```python
import json
from dnakit.design_crispr import predict_editing_outcomes

result = predict_editing_outcomes("GACGTACGTCGTAAGCTTCAAGG")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "predict_editing_outcomes",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 23,
    "sequence_bases": 23,
    "lengths": [
      23
    ],
    "sha256": "e17f729834a93bd40c59a9c547a2621d3e85a9c76868193434768ff96cb3b1d2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained repair-outcome model (inDelphi, FORECasT, Lindel) and its weights",
  "message": "Indel/repair outcome prediction needs trained model weights that are not bundled."
}
```

## 13) rank

- **Purpose:** Rank guide RNAs by a numeric field; the default `proxy_score` is a rule-based heuristic (see `score`).

- **Method:** Records are sorted by the chosen field, descending unless `ascending`; ties keep the input order (stable sort) and `top_n` truncates the list. When `by` is `proxy_score` and the field is absent it is computed with the same rule-based proxy as `score` (start at 100; -15 for GC outside 40-70 %, -40 instead outside 20-80 %; -40 for `TTTT`; -15 for an A/C/G homopolymer >= 4). `proxy_score` is a rule-based heuristic and is not fitted to or validated against editing data; it is not a predicted editing efficiency (use Rule Set 2, DeepSpCas9 or similar for that). Any other field (for example `off_target_sites` from `design`) must be numeric in every record.

- **API:** `dnakit.design_crispr.rank(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field. `options`: `by` (default `"proxy_score"`), `ascending` (default `false`), `top_n` (optional). A record without a numeric `by` field raises an input error.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `by`, `ascending` and `ranked` (the records, with `proxy_score` added when used, each with a 1-based `rank`). When `by` is `proxy_score`, `warnings` states that it is a heuristic.

- **Example code:**

```python
import json
from dnakit.design_crispr import rank

guides = ["GGGGGGGGGGGGGGGGGGGG", "GACGTACGTCGTAAGCTTCA", "ACGTACGTTTTTACGTACGT"]
result = rank(guides)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "rank",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 72,
    "sequence_bases": 60,
    "lengths": [
      20,
      20,
      20
    ],
    "sha256": "9bb60ae6e4e1ebb70a61019350afe0fb603f219157932ffb6cd98efba658b434"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "by": "proxy_score",
    "ascending": false,
    "ranked": [
      {
        "id": "guide2",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "proxy_score": 100.0,
        "rank": 1
      },
      {
        "id": "guide3",
        "sequence": "ACGTACGTTTTTACGTACGT",
        "proxy_score": 60.0,
        "rank": 2
      },
      {
        "id": "guide1",
        "sequence": "GGGGGGGGGGGGGGGGGGGG",
        "proxy_score": 45.0,
        "rank": 3
      }
    ]
  }
}
```

## 14) remove_restriction_sites

- **Purpose:** Remove restriction-enzyme recognition sites from a DNA sequence by minimal base substitutions.

- **Method:** Delegates to `dnakit.molbio.optimize_sequence_rules`, a greedy, seed-deterministic optimizer, with the recognition sequences as forbidden motifs on both strands; sites are counted before and after (the count accepts IUPAC patterns). With `preserve_translation` the input must be an in-frame coding sequence and only synonymous changes are made. Sites that cannot be removed under the constraints stay and `removed_all` is false. The edit is not checked for creating new, unlisted sites or for effects on expression.

- **API:** `dnakit.design_crispr.remove_restriction_sites(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: DNA string, FASTA text/path or `{"sequence": ...}` (one sequence). `options`: `enzymes` (names from the built-in catalog such as `"EcoRI"`) and/or `sites` (IUPAC recognition sequences), at least one is required; `preserve_translation` (default `false`); `seed` (default 0).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `sequence_id`, `sequence` (edited), `sites_before`, `sites_after`, `removed_all` and `changes` (each with 0-based `position`, `before`, `after`).

- **Example code:**

```python
import json
from dnakit.design_crispr import remove_restriction_sites

result = remove_restriction_sites("TTTTGAATTCAAAACCCGGATCCTTTT", options={"enzymes": ["EcoRI", "BamHI"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "remove_restriction_sites",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 27,
    "sequence_bases": 27,
    "lengths": [
      27
    ],
    "sha256": "f370022c244baaa8f4faedb685da22ddd8021fff5f05b2aaa52b195cdc90afe7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_id": "seq1",
    "sequence": "TTTTGCATTCAAAACCCAGATCCTTTT",
    "sites_before": 2,
    "sites_after": 0,
    "removed_all": true,
    "changes": [
      {
        "position": 5,
        "before": "A",
        "after": "C"
      },
      {
        "position": 17,
        "before": "G",
        "after": "A"
      }
    ]
  }
}
```

## 15) run_crispr_screen_design

- **Purpose:** Design a pooled CRISPR guide library: the best filtered guides per target sequence.

- **Method:** For each target, guides are enumerated and filtered exactly as in `design`, sorted by the heuristic `proxy_score`, and picked greedily until `guides_per_target` are chosen, skipping a guide whose cut site is closer than `min_cut_spacing` bp to an already chosen one (the guide start is used when there is no Cas9 cut site) and, with `references` and `max_off_target_sites`, a guide with more non-self Hamming sites than allowed. `proxy_score` is a rule-based heuristic and is not fitted to or validated against editing data; it is not a predicted editing efficiency (use Rule Set 2, DeepSpCas9 or similar for that). No non-targeting controls are generated and no gene- or exon-level annotation is used.

- **API:** `dnakit.design_crispr.run_crispr_screen_design(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: targets as a mapping `{target_id: DNA}`, a list of sequences or FASTA text/path (for example exons), up to 10000. `options`: `guides_per_target` (default 4); the enumeration and filter options of `design` (`nuclease`/`pam`, `min_gc` 0.4, `max_gc` 0.7, `max_homopolymer` 4, `avoid_poly_t`, `exclude_motifs`, `avoid_enzymes`); `min_cut_spacing` (default 0); optional `references`, `max_mismatches` (default 3) and `max_off_target_sites` (default: no limit).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `guides_per_target`, `target_count`, `guide_count`, `library` (rows with `target`, `guide_id`, the guide fields of `design` and, with `references`, the off-target counts) and `targets_short` (targets that got fewer guides than requested, with `requested` and `designed`).

- **Example code:**

```python
import json
from dnakit.design_crispr import run_crispr_screen_design

exon = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
result = run_crispr_screen_design({"exon1": exon}, options={"guides_per_target": 2, "min_cut_spacing": 5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "run_crispr_screen_design",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7e584ed8cb48bf6be7db6302c06095c9fde603d3140b4f062b450b482e52a1fb"
  },
  "warnings": [
    "proxy_score is a rule-based heuristic, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "guides_per_target": 2,
    "target_count": 1,
    "guide_count": 2,
    "library": [
      {
        "target": "exon1",
        "guide_id": "exon1_g1",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "pam": "AGG",
        "strand": "forward",
        "start": 0,
        "end": 20,
        "pam_start": 20,
        "cut_site": 17,
        "gc": 0.5,
        "proxy_score": 100.0
      },
      {
        "target": "exon1",
        "guide_id": "exon1_g2",
        "sequence": "AGCTTCAAGGCCAGCTAGCT",
        "pam": "AGG",
        "strand": "forward",
        "start": 13,
        "end": 33,
        "pam_start": 33,
        "cut_site": 30,
        "gc": 0.55,
        "proxy_score": 100.0
      }
    ],
    "targets_short": []
  }
}
```

## 16) scan_grna_offtargets

- **Purpose:** Find off-target sites of guide RNAs in supplied sequences by exact ungapped Hamming search on both strands.

- **Method:** Every window of the guide length on either strand of each reference is compared with the guide; a site is reported when it differs at <= `max_mismatches` positions and, if `pam` is set, the adjacent PAM (3' or 5' side) matches the IUPAC pattern. The guide's own target is reported too (0 mismatches). There are no bulges or indels and no CFD/MIT-style position-weighted scoring, so the mismatch count is the only ranking; the search is brute force and refuses jobs with guides x reference bases x 2 above 50 million. A reference with ambiguity codes is accepted, an ambiguity code in the reference counts as a mismatch.

- **API:** `dnakit.design_crispr.scan_grna_offtargets(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field (at most 1000). `options`: `references` (required: mapping `{id: sequence}`, list of sequences or FASTA text/path), `max_mismatches` (default 3, range 0-6), `pam` (IUPAC, default `"NGG"`; `null` ignores the PAM), `pam_side` (`"3prime"` default | `"5prime"`). Without `references` the call raises a configuration error.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `max_mismatches`, `pam` and `guides` (per guide `id`, `sequence`, `hit_count`, `counts_by_mismatch` and `hits` sorted by mismatches, each with `reference`, `strand` (`+`/`-`), 0-based half-open forward `start`/`end`, `mismatches`, 1-based `mismatch_positions` from the guide 5' end and the matched `site`).

- **Example code:**

```python
import json
from dnakit.design_crispr import scan_grna_offtargets

guide = "GACGTACGTCGTAAGCTTCA"
one_mismatch = guide[:9] + "A" + guide[10:]
reference = guide + "TGG" + "TT" + one_mismatch + "CGG" + "TT" + one_mismatch + "AAA"   # last site has no NGG PAM
result = scan_grna_offtargets(guide, options={"references": {"chr1": reference}, "max_mismatches": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "scan_grna_offtargets",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 20,
    "lengths": [
      20
    ],
    "sha256": "b015b64a47743134d8b718d27eacebc4a00519b9304d5f0f27607012ff7efa0c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "max_mismatches": 3,
    "pam": "NGG",
    "guides": [
      {
        "id": "seq1",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "hit_count": 2,
        "counts_by_mismatch": {
          "0": 1,
          "1": 1
        },
        "hits": [
          {
            "reference": "chr1",
            "strand": "+",
            "start": 0,
            "end": 20,
            "mismatches": 0,
            "mismatch_positions": [],
            "site": "GACGTACGTCGTAAGCTTCA"
          },
          {
            "reference": "chr1",
            "strand": "+",
            "start": 25,
            "end": 45,
            "mismatches": 1,
            "mismatch_positions": [
              10
            ],
            "site": "GACGTACGTAGTAAGCTTCA"
          }
        ]
      }
    ]
  }
}
```

## 17) score

- **Purpose:** Score guide RNAs with a rule-based heuristic on-target proxy in [0, 100]; it is not a trained efficiency model.

- **Method:** The proxy starts at 100 and subtracts 15 when the GC fraction lies outside 40-70 % (40 instead when outside 20-80 %), 40 for a `TTTT` run (RNA Pol III terminator) and 15 for an A/C/G homopolymer of length >= 4, with a floor of 0. It encodes only common design rules of thumb. `proxy_score` is a rule-based heuristic and is not fitted to or validated against editing data; it is not a predicted editing efficiency (use Rule Set 2, DeepSpCas9 or similar for that). Guides must be 10-40 nt of ACGT (U is read as T).

- **API:** `dnakit.design_crispr.score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field. `options`: none used.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `method` (`"heuristic_rule_proxy"`) and `scores` (per guide `id`, `sequence`, `proxy_score`, `gc`, `longest_homopolymer`, `poly_t`); `warnings` states that the score is a heuristic. A guide outside 10-40 nt or with other characters raises an input error.

- **Example code:**

```python
import json
from dnakit.design_crispr import score

result = score({"balanced": "GACGTACGTCGTAAGCTTCA", "poly_t": "ACGTACGTTTTTACGTACGT",
                "all_g": "GGGGGGGGGGGGGGGGGGGG"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 103,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6c8cb8de2dc6b2a2809628762c7b49bc5a2f66b8c7654fb54181ebb6d99e9f63"
  },
  "warnings": [
    "Heuristic rule-based proxy, not a trained on-target efficiency model."
  ],
  "status": "ok",
  "result": {
    "method": "heuristic_rule_proxy",
    "scores": [
      {
        "id": "balanced",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "proxy_score": 100.0,
        "gc": 0.5,
        "longest_homopolymer": 2,
        "poly_t": false
      },
      {
        "id": "poly_t",
        "sequence": "ACGTACGTTTTTACGTACGT",
        "proxy_score": 60.0,
        "gc": 0.4,
        "longest_homopolymer": 5,
        "poly_t": true
      },
      {
        "id": "all_g",
        "sequence": "GGGGGGGGGGGGGGGGGGGG",
        "proxy_score": 45.0,
        "gc": 1.0,
        "longest_homopolymer": 20,
        "poly_t": false
      }
    ]
  }
}
```

## 18) simulate

- **Purpose:** Simulate homology-directed repair in silico: swap the region between a donor's homology arms into the target.

- **Method:** The donor is split into left arm, payload and right arm by the given arm lengths. Both arms must occur exactly once on the forward strand of the target, the left one before the right one; the target sequence between them is replaced by the payload. Repair-pathway competition, cutting efficiency, NHEJ indels and mutations in the arms are not modelled, so this is string substitution and not an efficiency prediction. A missing, repeated or mis-ordered arm raises an input error.

- **API:** `dnakit.design_crispr.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: target DNA (string, FASTA text/path or `{"sequence": ...}`; one sequence). `options`: `donor` (required, left arm + payload + right arm), `left_arm` and `right_arm` (required arm lengths that fit in the donor).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `target`, `edited_sequence`, `region` (replaced 0-based half-open `start`/`end`), `deleted`, `inserted` and `length_change`.

- **Example code:**

```python
import json
from dnakit.design_crispr import simulate

target = "AAAACCCCTTGGGGTTTT"
donor = "CCCC" + "AA" + "GGGG"          # left arm + payload + right arm
result = simulate(target, options={"donor": donor, "left_arm": 4, "right_arm": 4})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 18,
    "lengths": [
      18
    ],
    "sha256": "82435c45a7de0432d66f77bce4656ba8124463ba82f4c4688556c6b902bf5651"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "target": "seq1",
    "edited_sequence": "AAAACCCCAAGGGGTTTT",
    "region": {
      "start": 8,
      "end": 10
    },
    "deleted": "TT",
    "inserted": "AA",
    "length_change": 0
  }
}
```

## 19) validate

- **Purpose:** Validate guide RNAs: alphabet, length, PAM adjacency in a target (optional) and design-rule warnings.

- **Method:** Errors make a guide invalid: characters outside ACGT (`invalid_characters`), a length different from the expected one (`length_<n>_expected_<m>`) and, when a `target` is given, no occurrence of the guide next to a matching PAM on either strand (`no_pam_adjacent_site_in_target`). Warnings do not invalidate: GC outside 40-70 % (`gc_outside_40_70`), a homopolymer run longer than 4 (`homopolymer_run_gt_4`) and `TTTT` (`poly_t`). Target sites are found by exact string matching on both strands; there is no mismatch tolerance and no efficiency or off-target assessment.

- **API:** `dnakit.design_crispr.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: guides as a list of strings, a FASTA text/path, a mapping `{id: guide}` or a list of records with a `guide`/`sequence` field. `options`: `length` (expected guide length, default the nuclease default, 20 for SpCas9); `nuclease` (`"SpCas9"` default, `"SaCas9"`, `"AsCas12a"`) or custom `pam`/`pam_side`/`guide_length` as in `design`; `target` (optional DNA; the guide must then occur in it on either strand next to a matching PAM).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"design_crispr"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256; a mapping input counts as one record), `warnings` and `result`. `result` has `nuclease`, `pam`, `expected_length`, `all_valid` and `guides` (per guide `id`, `sequence`, `valid`, `errors`, `warnings` and `target_sites` with `strand`, 0-based half-open forward `start`/`end`, `pam`).

- **Example code:**

```python
import json
from dnakit.design_crispr import validate

target = "GACGTACGTCGTAAGCTTCA" + "AGG" + "CCA" + "GCTAGCTAGGCTAATCGGAT"
guides = {"good": "GACGTACGTCGTAAGCTTCA", "short": "GACGTACGTCGTAAGCTTC", "no_pam": "CATCATCATCATCATCATCA"}
result = validate(guides, options={"target": target})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "design_crispr",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 98,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d888e1a0acbc03ba69ad873c3fd9569ce541e899e96b5a9b738889ecb86b1d1c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "nuclease": "SpCas9",
    "pam": "NGG",
    "expected_length": 20,
    "guides": [
      {
        "id": "good",
        "sequence": "GACGTACGTCGTAAGCTTCA",
        "valid": true,
        "errors": [],
        "warnings": [],
        "target_sites": [
          {
            "strand": "forward",
            "start": 0,
            "end": 20,
            "pam": "AGG"
          }
        ]
      },
      {
        "id": "short",
        "sequence": "GACGTACGTCGTAAGCTTC",
        "valid": false,
        "errors": [
          "length_19_expected_20"
        ],
        "warnings": [],
        "target_sites": []
      },
      {
        "id": "no_pam",
        "sequence": "CATCATCATCATCATCATCA",
        "valid": false,
        "errors": [
          "no_pam_adjacent_site_in_target"
        ],
        "warnings": [
          "gc_outside_40_70"
        ],
        "target_sites": []
      }
    ],
    "all_valid": false
  }
}
```
