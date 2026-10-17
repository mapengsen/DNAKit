# Molecular laboratory

This category contains **10** documented capabilities. All 10 compute a result natively (pure Python; Tm and extinction coefficients come from `dnakit.thermodynamics`) and return `status: "ok"`; none returns `not_implemented`, and none needs an external program or optional library. They are bench calculators, protocol planners and small sequence tools (Sanger trace base calling and read assembly, primer analysis, optimization and validation), not replacements for dedicated assemblers or primer-design suites. The functions are exported through the formal DNAKit namespace and also exposed through the Agent/MCP registry where the JSON call plan is compatible. The example outputs below are real outputs; long arrays and strings are abbreviated with `...`.

## 1) analyze

- **Purpose:** Analyzes oligonucleotide sequences: composition, GC fraction, reverse complement, molecular weight, A260 extinction coefficient and melting temperature.

- **Method:** Per sequence it counts bases, computes the GC fraction, the reverse complement, the longest homopolymer and whether the 3' end has a G/C clamp. The molecular weight is that of the anhydrous single strand with a 5'-OH, and the 260 nm extinction coefficient is the theoretical nearest-neighbour value (with hypochromicity) from `dnakit.thermodynamics`, which also gives nmol and ug per OD260. Tm uses the SantaLucia (1998) nearest-neighbour model by default (`sodium_mm`, `primer_nm`); `wallace` (2 per A/T plus 4 per G/C) applies to 2-13 nt and `basic` (64.9 + 41 (G+C-16.4)/N) to sequences of at least 14 nt, otherwise Tm is `null` with a `tm_note`. Pure Python; no external program.

- **API:** `dnakit.molecular_lab.analyze(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a DNA string (A, C, G, T only), FASTA text, a list of sequences or `{id, sequence}` mappings, or a mapping name -> sequence (up to 1000 sequences of 100 kb). Options: `tm_method` (`nearest_neighbor` default, `wallace`, `basic`), `sodium_mm` (default 50), `primer_nm` (default 250).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `sequences` (per sequence: `length`, `composition`, `gc_fraction`, `reverse_complement`, `is_self_complementary`, `longest_homopolymer`, `ends_with_gc_clamp`, `molecular_weight_da`, `extinction_coefficient_260nm`, `nmol_per_od260`, `ug_per_od260`, `tm_celsius`, `tm_method`), `count` and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import analyze

result = analyze({"primer": "CGTTGACCGTAGCTAGGCTA"}, options={"sodium_mm": 50, "primer_nm": 250})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "analyze",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "874a8555eeae8766308c0b460b1ba35a68a63691bbf0f8bb6ce8790a59570a39"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequences": {
      "primer": {
        "length": 20,
        "composition": {
          "A": 4,
          "C": 5,
          "G": 6,
          "T": 5
        },
        "gc_fraction": 0.55,
        "reverse_complement": "TAGCCTAGCTACGGTCAACG",
        "is_self_complementary": false,
        "longest_homopolymer": 2,
        "ends_with_gc_clamp": false,
        "molecular_weight_da": 6133.04,
        "extinction_coefficient_260nm": 191900.0,
        "nmol_per_od260": 5.21105,
        "ug_per_od260": 31.9596,
        "tm_celsius": 55.432,
        "tm_method": "nearest_neighbor"
      }
    },
    "count": 1,
    "definition": "Composition, GC fraction, reverse complement, anhydrous single-strand molecular weight, theoretical 260 nm extinction coefficient (nearest-neighbour hypochromic..."
  }
}
```

## 2) assemble_sanger_reads

- **Purpose:** Assembles overlapping Sanger reads into contigs by greedy suffix/prefix overlap merging, considering both strands.

- **Method:** Greedy overlap-layout: the pair of contigs (or reads) with the longest suffix/prefix overlap of at least `min_overlap` bases and ungapped identity of at least `min_identity` is merged, repeatedly. With `both_strands` a read may be reverse complemented to fit. In an overlap the consensus keeps the earlier contig's base: there is no quality weighting and no gap handling, so reads should be quality-trimmed first. Contained reads are absorbed; reads that join nothing are returned as singlets. Pure Python; no external assembler.

- **API:** `dnakit.molecular_lab.assemble_sanger_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is FASTA text, a list of sequences or `{id, sequence}` mappings, or a mapping name -> sequence (up to 100 reads of 20 kb; IUPAC letters allowed). Options: `min_overlap` (default 20 bases), `min_identity` (default 0.95, ungapped), `both_strands` (default `True`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `reads`, `contigs` (per contig `id`, `length`, `sequence` and `reads` with `id`, `strand`, `offset`, `length`), `singlets`, the options used and a `definition` string. The contig orientation is arbitrary: it can be the reverse complement of the original template. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import assemble_sanger_reads

template = "ATGGCTAGCTAGGATCCGATTACAGGCTTAACGGATCGCTAGGCATCGATTCGGCAAGTCTTGACCATGCAGG"
reverse_complement = template[40:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
reads = {"f1": template[:30], "f2": template[20:50], "r3": reverse_complement}
result = assemble_sanger_reads(reads, options={"min_overlap": 10})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "assemble_sanger_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 123,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b7ac51ed6b5dfc15af47ad4da1d341477c59246037d2cf25852f57ad54a4d0b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "contigs": [
      {
        "id": "contig1",
        "length": 73,
        "sequence": "CCTGCATGGTCAAGACTTGCCGAATCGATGCCTAGCGATCCGTTAAGCCTGTAATCGGATCCTAGCTAGCCAT",
        "reads": [
          {
            "id": "r3",
            "strand": "+",
            "offset": 0,
            "length": 33
          },
          {
            "id": "f2",
            "strand": "-",
            "offset": 23,
            "length": 30
          },
          {
            "id": "f1",
            "strand": "-",
            "offset": 43,
            "length": 30
          }
        ]
      }
    ],
    "singlets": [],
    "min_overlap": 10,
    "min_identity": 0.95,
    "both_strands": true,
    "definition": "Greedy overlap-layout: repeatedly merge the pair of contigs with the longest suffix/prefix overlap (>= min_overlap bases, identity >= min_identity, reads or the..."
  }
}
```

## 3) calculate

- **Purpose:** Performs a bench calculation selected by `options["calculation"]`: dilution, molarity, copy number, ligation, A260, Tm, PCR yield or gel size.

- **Method:** Parameters are read from the input mapping and from `options` (options win). `dilution`: three of c1, v1, c2, v2 (C1 V1 = C2 V2) give the fourth and the diluent volume. `molarity`: converts mass and molar concentration using 660 g/mol per bp (`ds`) or 330 g/mol per nt (`ss`), or the exact weight when `sequence` is given. `copy_number`: copies = mass / MW x 6.02214076e23. `ligation`: insert mass = ratio x vector mass x insert bp / vector bp, with volumes when concentrations are given. `a260`: concentration from A260 with 50 (dsDNA), 33 (ssDNA) or 40 (RNA) ng/uL per A260 unit, or from the oligo extinction coefficient when `sequence` is given. `tm`: nearest-neighbour (default), `wallace` or `basic`. `pcr_yield`: N0 x (1 + E)^cycles. `gel_size`: least-squares fit of log10(size) against migration distance. Pure Python; Tm and extinction coefficients use `dnakit.thermodynamics`.

- **API:** `dnakit.molecular_lab.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping of parameters; `options["calculation"]` is required and selects the calculation. `dilution`: three of `c1`, `v1`, `c2`, `v2`. `molarity`: `length` (bp) or `sequence`, `strand` (`ds` default or `ss`), and one of `mass_concentration_ng_per_ul` or `molar_concentration_nm`. `copy_number`: `mass_ng`, `length` or `sequence`, `strand`. `ligation`: `vector_bp`, `insert_bp`, `vector_ng`, `insert_to_vector_ratio` (default 3), optional `vector_conc_ng_per_ul` and `insert_conc_ng_per_ul`. `a260`: `a260`, `nucleic_acid` (`dsDNA`, `ssDNA`, `RNA`) or `sequence`, `dilution_factor`, `path_length_cm`. `tm`: `sequence`, `method` (`nearest_neighbor` default, `wallace`, `basic`), `sodium_mm` (50), `primer_nm` (250). `pcr_yield`: `template_copies`, `cycles`, `efficiency` (default 1.0). `gel_size`: `ladder` (list of `{bp, distance}`) and `distances`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `calculation`, a `formula` string and fields that depend on the calculation: `dilution` the four quantities plus `solved`, `diluent_volume` and `dilution_factor`; `molarity` `molecular_weight`, `mw_basis`, `mass_concentration_ng_per_ul`, `molar_concentration_nm`, `pmol_per_ul`; `copy_number` `molecular_weight`, `mass_ng`, `copies`; `ligation` `insert_ng`, `vector_ng`, `insert_to_vector_ratio`, `vector_pmol`, `insert_pmol` (and volumes when concentrations are given); `a260` `method`, `mass_concentration_ng_per_ul`; `tm` `tm_celsius`, `gc_fraction`, `conditions`, `model`; `pcr_yield` `product_copies`, `fold_amplification`; `gel_size` `slope_log10bp_per_distance`, `intercept_log10bp`, `r_squared`, `estimated_bp`, `extrapolated_points`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import calculate

parameters = {"vector_bp": 3000, "insert_bp": 1000, "vector_ng": 50, "vector_conc_ng_per_ul": 25, "insert_conc_ng_per_ul": 10}
result = calculate(parameters, options={"calculation": "ligation"})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 113,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50280e11e9bf70b18a72becc28c313b3cb8c918d937adb6ffd3f004a5dd1f027"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "calculation": "ligation",
    "insert_ng": 50.0,
    "vector_ng": 50.0,
    "insert_to_vector_ratio": 3.0,
    "vector_pmol": 0.0252525,
    "insert_pmol": 0.0757576,
    "formula": "insert_ng = ratio * vector_ng * insert_bp / vector_bp (equal ratio of ends, 660 g/mol per bp)",
    "vector_volume_ul": 2.0,
    "insert_volume_ul": 5.0
  }
}
```

## 4) call_bases_from_sanger

- **Purpose:** Calls bases from Sanger chromatogram traces by picking the strongest channel at each peak, with IUPAC codes for mixed peaks.

- **Method:** At each peak location the base is the channel with the highest signal. If the second-highest channel is at least `mixed_peak_ratio` of the highest, the two-base IUPAC code is reported and the position is listed under `mixed_positions`. For an ABIF `.ab1` file the four channels are DATA9-12 ordered by FWO_1 and the peaks come from PLOC2; the result also gives the agreement with the instrument's own base calls. When no peak locations are supplied they are detected as local maxima of the per-sample maximum across the four channels, above `min_peak_height` and at least `min_peak_spacing` samples apart. No Phred-style quality values are estimated. Pure Python; no external base caller.

- **API:** `dnakit.molecular_lab.call_bases_from_sanger(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is an ABIF `.ab1` file path, or a mapping `{"traces": {"A": [...], "C": [...], "G": [...], "T": [...]}, "peak_locations": [...]}` where `peak_locations` is optional. Options: `mixed_peak_ratio` (default 0.5), `min_peak_height` (default 0.1 of the global maximum, used when peaks are detected), `min_peak_spacing` (default 3 samples).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `source`, `samples`, `peaks`, `peaks_detected_automatically`, `sequence`, `secondary_to_primary_peak_ratio` (per peak), `mixed_positions` (`index`, `sample`, `bases`, `code`, `secondary_to_primary`), `peak_locations` and a `definition` string; for `.ab1` files also the channel order and `agreement_with_instrument_basecalls`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import call_bases_from_sanger

import math

def peak(n, centres, sigma=2.0):
    return [sum(h * math.exp(-((i - c) ** 2) / (2 * sigma**2)) for c, h in centres.items()) for i in range(n)]

traces = {
    "A": peak(100, {10: 100, 90: 100, 50: 60}),
    "C": peak(100, {30: 100}),
    "G": peak(100, {50: 100}),
    "T": peak(100, {70: 100}),
}
result = call_bases_from_sanger({"traces": traces, "peak_locations": [10, 30, 50, 70, 90]})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "call_bases_from_sanger",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 9056,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8570c3e0ee1c1dfb6b9f2ea92f12e1b21c7946c975fc500d682a7d0a69c963e9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "mapping",
    "samples": 100,
    "peaks": 5,
    "peaks_detected_automatically": false,
    "sequence": "ACRTA",
    "secondary_to_primary_peak_ratio": [
      0.0,
      0.0,
      0.6,
      0.0,
      "..."
    ],
    "mixed_positions": [
      {
        "index": 2,
        "sample": 50,
        "bases": "GA",
        "code": "R",
        "secondary_to_primary": 0.6
      }
    ],
    "peak_locations": [
      10,
      30,
      50,
      70,
      "..."
    ],
    "definition": "At each peak location the base is the channel with the highest signal; if the second-highest channel is at least mixed_peak_ratio (default 0.5) of it, the IUPAC..."
  }
}
```

## 5) export

- **Purpose:** Serialises lab records (samples, primers, sequences) to CSV, TSV, JSON or FASTA text, optionally writing a file.

- **Method:** Each record becomes a row (CSV/TSV, standard quoting), a JSON object or a FASTA entry whose header is the record's `id` or `name` and whose sequence is wrapped at `line_width`. Columns default to the fields in order of appearance. A file is written only when `output_path` is given, and an existing file is never overwritten unless `overwrite=True`. Pure Python.

- **API:** `dnakit.molecular_lab.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a list of record mappings (or sequence strings), a single record, or a mapping id -> sequence. Options: `format` (`csv` default, `tsv`, `json`, `fasta`), `columns` (field order; default: fields in order of appearance), `line_width` (FASTA, default 80), `output_path` (optional file), `overwrite` (default `False`). FASTA needs a `sequence` (or `seq`) field.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `format`, `records`, `columns`, `text` (the serialised text) and `output_path` (`null` unless a file was written). Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import export

records = [{"id": "p1", "sequence": "ACGT", "conc": 5}, {"id": "p2", "sequence": "GG", "conc": None}]
result = export(records, options={"format": "fasta"})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 91,
    "sequence_bases": 6,
    "lengths": [
      4,
      2
    ],
    "sha256": "540186960d98208a3e6d57ad53aa6c394164e5559c0e9ac130c9bdb7d63aa38e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "records": 2,
    "columns": [
      "id",
      "sequence",
      "conc"
    ],
    "text": ">p1\nACGT\n>p2\nGG\n",
    "output_path": null
  }
}
```

## 6) optimize

- **Purpose:** Optimizes primer length to reach a target melting temperature.

- **Method:** The primer is the first n bases of the template (`forward`) or the reverse complement of its last n bases (`reverse`), for n from `min_length` to `max_length`. The Tm of every candidate is computed (SantaLucia 1998 nearest-neighbour via `dnakit.thermodynamics` with 50 mM Na+ and 250 nM primer, or the `basic` formula) and the length whose Tm is closest to `target_tm` is chosen; ties go to the shorter primer. With `require_gc_clamp` only candidates whose 3' base is G or C are eligible (an error is raised if none is). Secondary structure and specificity are not evaluated.

- **API:** `dnakit.molecular_lab.optimize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is the template DNA string. Options: `target_tm` (default 60), `min_length` (18), `max_length` (30), `direction` (`forward` default or `reverse`), `require_gc_clamp` (default `False`), `tm_method` (`nearest_neighbor` default or `basic`; `sodium_mm` 50 and `primer_nm` 250 for the former).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `direction`, `target_tm`, `tm_method`, `best` (the chosen candidate: `length`, `primer`, `tm_celsius`, `tm_error`, `gc_fraction`, `gc_clamp`), `candidates` (the same fields for every length tried) and a `definition` string. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import optimize

template = "CGTTGACCGTAGCTAGGCTAATCGGATTCGCAGGT"
result = optimize(template, options={"target_tm": 60})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "optimize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 35,
    "sequence_bases": 35,
    "lengths": [
      35
    ],
    "sha256": "3de911088164c6583899211652ff4c24cfac73b79ceab8df82c00b738ee6f019"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "direction": "forward",
    "target_tm": 60.0,
    "tm_method": "nearest_neighbor",
    "best": {
      "length": 25,
      "primer": "CGTTGACCGTAGCTAGGCTAATCGG",
      "tm_celsius": 60.835,
      "tm_error": 0.83537,
      "gc_fraction": 0.56,
      "gc_clamp": true
    },
    "candidates": [
      {
        "length": 18,
        "primer": "CGTTGACCGTAGCTAGGC",
        "tm_celsius": 54.425,
        "tm_error": 5.5755,
        "gc_fraction": 0.611111,
        "gc_clamp": true
      },
      {
        "length": 19,
        "primer": "CGTTGACCGTAGCTAGGCT",
        "tm_celsius": 55.786,
        "tm_error": 4.2139,
        "gc_fraction": 0.578947,
        "gc_clamp": false
      },
      {
        "length": 20,
        "primer": "CGTTGACCGTAGCTAGGCTA",
        "tm_celsius": 55.432,
        "tm_error": 4.5675,
        "gc_fraction": 0.55,
        "gc_clamp": false
      },
      "..."
    ],
    "definition": "Primers are the first (forward) or reverse complement of the last (reverse) min_length..max_length bases of the template; the length whose Tm is closest to targ..."
  }
}
```

## 7) parse

- **Purpose:** Parses laboratory quantity strings (`50 ng/uL`, `10 mM`, `200 uL`) into SI values, or expands plate well ranges (`A1:B3`).

- **Method:** A quantity is `<number> <unit>` with the units g, L, mol and M, the prefixes k, m, u (or the micro sign), n, p, f, mass/volume or amount/volume ratios, and `X` for a fold concentration; the value is converted to g, L, mol, g/L or mol/L and its dimension is reported. For `kind` `wells`, single wells and rectangular ranges on a 96-well (8x12) or 384-well (16x24) plate are expanded in row-major or column-major order. Pure Python.

- **API:** `dnakit.molecular_lab.parse(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a quantity string or a list of them (`kind` `quantity`), or well text such as `A1:B3, H12` (`kind` `wells`). Options: `kind` (`quantity` default or `wells`), `plate` (96 default or 384) and `order` (`row` default or `column`) for wells.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds for quantities `kind`, `quantities` (per item `text`, `value`, `unit`, `dimension`, `si_value`, `si_unit`) and a `definition` string; for wells `kind`, `plate`, `order`, `wells`, `count` and a `definition`. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import parse

result = parse(["50 ng/uL", "1.5 mM", "200 \u00b5L", "5X"])
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "parse",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 38,
    "sequence_bases": 8,
    "lengths": [4, 2, 1, 1],
    "sha256": "825a57d9c055f721e499e4afc93c58da3805cc480a7ed180994258ebc4f722f0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "quantity",
    "quantities": [
      {
        "text": "50 ng/uL",
        "value": 50.0,
        "unit": "ng/uL",
        "dimension": "mass_concentration",
        "si_value": 0.05,
        "si_unit": "g/L"
      },
      {
        "text": "1.5 mM",
        "value": 1.5,
        "unit": "mM",
        "dimension": "molar_concentration",
        "si_value": 0.0015,
        "si_unit": "mol/L"
      },
      {
        "text": "200 µL",
        "value": 200.0,
        "unit": "µL",
        "dimension": "volume",
        "si_value": 0.0002,
        "si_unit": "L"
      },
      "..."
    ],
    "definition": "'<number> <unit>' with units g, L, mol, M and prefixes k, m, u (or the micro sign), n, p, f; ratios mass/volume and amount/volume; 'X' for fold concentration. S..."
  }
}
```

## 8) plan

- **Purpose:** Plans a bench protocol with exact volumes: PCR master mix, serial dilution or ligation setup.

- **Method:** Parameters are read from the input mapping and from `options` (options win), and `protocol` selects the plan. `pcr`: each component volume = final / stock x reaction volume, water makes up the remainder, and the master mix (everything except the template) is scaled by n_reactions x (1 + overage). `serial_dilution`: each tube is the previous one diluted by `dilution_factor` into `final_volume_ul`, with transfer and diluent volumes per tube. `ligation`: insert mass from the molar ratio (as in `calculate`), then volumes of vector, insert, 10x ligase buffer, ligase and water. A plan whose components exceed the reaction volume is rejected. Pure Python.

- **API:** `dnakit.molecular_lab.plan(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping of parameters (may be empty); `options["protocol"]` is `pcr` (default), `serial_dilution` or `ligation`. `pcr`: `reaction_volume_ul` (25), `n_reactions` (1), `overage_fraction` (0.1), `template_volume_ul` (1), `components` (list of `{name, stock, final, unit}`; default 1x buffer from 5x, 0.2 mM dNTP, 0.5 uM primers, 0.025 U/uL polymerase). `serial_dilution`: `start_concentration`, `dilution_factor` (10), `steps`, `final_volume_ul` (100). `ligation`: `vector_bp`, `insert_bp`, `vector_ng`, `insert_to_vector_ratio` (3), optional concentrations, `reaction_volume_ul` (20).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `protocol` and a `definition` string plus protocol-specific fields: `pcr` `reaction_volume_ul`, `n_reactions`, `overage_fraction`, `master_mix_for_reactions`, `components` (per component `stock`, `final`, `unit`, `volume_per_reaction_ul`, `master_mix_volume_ul`, including water), `master_mix_total_ul`, `master_mix_aliquot_per_tube_ul`, `template_per_tube_ul`; `serial_dilution` `dilution_factor`, `steps`, `tubes`, `total_diluent_ul`; `ligation` the ligation amounts of `calculate` plus `reaction_volume_ul`, `ligase_buffer_10x_ul`, `ligase_ul` and, with concentrations, the vector, insert and water volumes. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import plan

result = plan({}, options={"protocol": "pcr", "n_reactions": 10})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "plan",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 2,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "protocol": "pcr",
    "reaction_volume_ul": 25.0,
    "n_reactions": 10,
    "overage_fraction": 0.1,
    "master_mix_for_reactions": 11.0,
    "components": [
      {
        "component": "reaction buffer",
        "stock": 5.0,
        "final": 1.0,
        "unit": "X",
        "volume_per_reaction_ul": 5.0,
        "master_mix_volume_ul": 55.0
      },
      {
        "component": "dNTP mix",
        "stock": 10.0,
        "final": 0.2,
        "unit": "mM each",
        "volume_per_reaction_ul": 0.5,
        "master_mix_volume_ul": 5.5
      },
      {
        "component": "forward primer",
        "stock": 10.0,
        "final": 0.5,
        "unit": "uM",
        "volume_per_reaction_ul": 1.25,
        "master_mix_volume_ul": 13.75
      },
      "..."
    ],
    "master_mix_total_ul": 264.0,
    "master_mix_aliquot_per_tube_ul": 24.0,
    "template_per_tube_ul": 1.0,
    "definition": "volume = final / stock * reaction volume per component; water makes up the rest; the master mix (everything but the template) is scaled by n_reactions * (1 + ov..."
  }
}
```

## 9) plot

- **Purpose:** Draws a schematic virtual agarose gel as an SVG string.

- **Method:** Each lane is a list of fragment sizes in bp. Migration distance is linear in log10(size) between `min_bp` and `max_bp` (the usual semilog behaviour of agarose) and band brightness grows with fragment mass. It is a drawing, not a simulation of a particular gel, buffer or voltage. The SVG is returned as text and written to a file only when `output_path` is given. Pure Python; no plotting library.

- **API:** `dnakit.molecular_lab.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a mapping lane name -> list of fragment sizes in bp (or a list of such lists), for example `{"ladder": [1000, 500, 100], "PCR": [650]}`. Options: `min_bp` (default 100), `max_bp` (default 10000), `output_path` (optional `.svg`; an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `svg` (the SVG text), `output_path` (`null` unless a file was written), `lanes`, `min_bp`, `max_bp` and a `model` description. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import plot

result = plot({"ladder": [1000, 500, 100], "PCR": [650]})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1363816d06183a2cf60426b59ffcdf3e55f048b4076c468c821cc930afe14db7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"188\" height=\"400\" viewBox=\"0 0 188 400\" font-family=\"sans-serif\" font-size=\"10\">\n<title>Virtual agarose gel</titl...",
    "output_path": null,
    "lanes": 2,
    "min_bp": 100.0,
    "max_bp": 10000.0,
    "model": "Schematic gel: migration distance is linear in log10(size) between min_bp and max_bp (the usual semilog behaviour of agarose); band brightness grows with fragme..."
  }
}
```

## 10) validate

- **Purpose:** Validates PCR primers against rule-of-thumb design criteria.

- **Method:** Each primer is checked for length, GC content, a 3' G/C clamp (3' base G or C and at most 3 G/C among the last 5 bases), homopolymer runs, a dinucleotide unit repeated 4 or more times, and Tm (SantaLucia 1998 nearest-neighbour by default); for a forward/reverse pair the Tm difference is checked too. Every check reports `passed` and a `detail`. Passing does not predict amplification: secondary structure, primer dimers and off-target binding are not assessed. Pure Python.

- **API:** `dnakit.molecular_lab.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` is a primer sequence, a list or mapping of primers, or a pair `{"forward": ..., "reverse": ...}`. Options: `min_length` (18), `max_length` (30), `min_gc` (0.40), `max_gc` (0.60), `min_tm` (52), `max_tm` (65), `max_homopolymer` (4), `max_tm_difference` (pair, 5), `tm_method` (`nearest_neighbor` default), `sodium_mm`, `primer_nm`.

- **Output:** Returns a JSON-serializable envelope with `status` (`"ok"`), `module`, `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result`. `result` holds `kind` (`primer_pair` for a forward/reverse pair, otherwise `primers`), the overall `passed`, the per-primer reports (`forward` and `reverse` for a pair, a `primers` mapping name -> report otherwise; each with `sequence`, `tm_celsius`, `passed` and `checks`, every check having `check`, `passed` and `detail`), for pairs `pair_checks` (the Tm difference) and a `rules` description of the thresholds used. Invalid input or options raise `InputFormatError` or `ConfigurationError`.

- **Example code:**

```python
from dnakit.molecular_lab import validate

result = validate({"forward": "CGTTGACCGTAGCTAGGCTA", "reverse": "GCGCGCGCGCGCGCGCGCGC"})
print(result)
```

- **Example output:**

```json
{
  "module": "molecular_lab",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6a15be3b2d065bb6549f70ec94b2bba749a5f7594bfe14e633230b4ed427261e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "primer_pair",
    "passed": false,
    "forward": {
      "sequence": "CGTTGACCGTAGCTAGGCTA",
      "tm_celsius": 55.432,
      "passed": false,
      "checks": [
        {
          "check": "length",
          "passed": true,
          "detail": "20 nt (allowed 18-30)"
        },
        {
          "check": "gc_content",
          "passed": true,
          "detail": "55.0% (allowed 40-60%)"
        },
        {
          "check": "three_prime_gc_clamp",
          "passed": false,
          "detail": "last 5 nt GGCTA has 3 G/C; 3' base A (want 3' G/C and <= 3 G/C in the last 5)"
        },
        "..."
      ]
    },
    "reverse": {
      "sequence": "GCGCGCGCGCGCGCGCGCGC",
      "tm_celsius": 80.078,
      "passed": false,
      "checks": [
        {
          "check": "length",
          "passed": true,
          "detail": "20 nt (allowed 18-30)"
        },
        {
          "check": "gc_content",
          "passed": false,
          "detail": "100.0% (allowed 40-60%)"
        },
        {
          "check": "three_prime_gc_clamp",
          "passed": false,
          "detail": "last 5 nt CGCGC has 5 G/C; 3' base C (want 3' G/C and <= 3 G/C in the last 5)"
        },
        "..."
      ]
    },
    "pair_checks": [
      {
        "check": "tm_difference",
        "passed": false,
        "detail": "|Tm1 - Tm2| = 24.6 C (allowed <= 5)"
      }
    ],
    "rules": "Rule-of-thumb PCR primer criteria (adjustable): length 18-30 nt, GC 40-60%, Tm 52-65 C (nearest-neighbour, 50 mM Na+, 250 nM), 3' terminal G/C with at most 3 G/..."
  }
}
```
