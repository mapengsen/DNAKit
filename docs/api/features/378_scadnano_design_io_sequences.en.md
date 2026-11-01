# scadnano designs: files, exports, sequences and modifications

Functions that read and write scadnano designs (`.sc`), convert them to cadnano v2, oxView and IDT order formats, build plate maps, and work with DNA sequences and modifications of the strands.

These functions are wrappers of [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package) (MIT license, tested with version 0.21.1), a Python package for designing DNA nanostructures. scadnano is optional and is imported only when a function runs; install it with `pip install "dnakit[dna-nanotechnology]"` (the Excel order also needs openpyxl and the plate tables need tabulate, which scadnano installs). A design is a scadnano JSON mapping (the content of a `.sc` file). Editing functions never change their argument: they return the new design under `design`. Strands are addressed by their 0-based position in `design["strands"]` and domains by their position inside the strand; helix indices and offsets are scadnano's own. The functions do not predict folding, yield or stability. A design with sequences that were assigned before an insertion, deletion or domain change keeps the old sequence, so assign sequences last.

## 1) `dnakit.genomics.load_scadnano_design()`

- **Purpose:** Load and validate a scadnano design from a `.sc` file, a JSON string or a JSON mapping.
- **Method:** Calls `Design.from_scadnano_file`, `from_scadnano_json_str` or `from_scadnano_json_map` depending on the argument given, then serialises the validated design again.
- **API:** `dnakit.genomics.load_scadnano_design(design: 'Mapping | None' = None, *, path: 'str | None' = None, json_text: 'str | None' = None) -> 'dict[str, Any]'`.
- **Input:** Exactly one of `design` (mapping), `path` (file) or `json_text`.
- **Output:** `design` (normalised JSON mapping), grid, helix, strand and scaffold counts, strand names and lengths, `source_kind`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    load_scadnano_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = load_scadnano_design(path="design.sc")
print(r["source_kind"], r["helices"], r["strands"], r["strand_names"])
r = load_scadnano_design(json_text=open("design.sc").read())
print(r["source_kind"], r["strand_lengths"])
r = load_scadnano_design(design)
print(r["source_kind"], r["grid"])
```

- **Example output:**

```text
file 2 2 [None, 'stapA']
json_string [48, 24]
json_map square
```

## 2) `dnakit.genomics.export_scadnano_design()`

- **Purpose:** Serialise a design to scadnano `.sc` JSON text and optionally write it to a new file.
- **Method:** Calls `Design.to_json` (`suppress_indent` is scadnano's compact layout) and writes the text like `write_scadnano_file`; an existing output file is refused.
- **API:** `dnakit.genomics.export_scadnano_design(design: 'Mapping', *, output_path: 'str | None' = None, suppress_indent: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** A design mapping, an optional new `output_path`, `suppress_indent`.
- **Output:** `json_text`, `output_path`, design counts.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_scadnano_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

import os, tempfile
out = os.path.join(tempfile.mkdtemp(), "out.sc")
r = export_scadnano_design(design, output_path=out)
print(r["json_text"].splitlines()[0:3])
print(os.path.getsize(out) == len(r["json_text"].encode()))
```

- **Example output:**

```text
['{', '  "version": "0.21.1",', '  "grid": "square",']
True
```

## 3) `dnakit.genomics.import_cadnano_v2_design()`

- **Purpose:** Convert a cadnano v2 design to a scadnano design.
- **Method:** Calls `Design.from_cadnano_v2` with either the parsed mapping (`json_dict`) or a `.json` file.
- **API:** `dnakit.genomics.import_cadnano_v2_design(source: 'str | Mapping', *) -> 'dict[str, Any]'`.
- **Input:** A cadnano v2 mapping or the path of a cadnano v2 JSON file.
- **Output:** The converted `design` and its counts.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    export_cadnano_v2_design,
    import_cadnano_v2_design,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

cad = export_cadnano_v2_design(design)["cadnano_v2"]
r = import_cadnano_v2_design(cad)
print(r["grid"], r["helices"], r["strands"])
```

- **Example output:**

```text
square 4 21
```

## 4) `dnakit.genomics.export_cadnano_v2_design()`

- **Purpose:** Export a design as cadnano v2 JSON.
- **Method:** Calls `Design.to_cadnano_v2_json` (also the content `write_cadnano_v2_file` writes) and returns the text and its parsed form; optionally writes a new file.
- **API:** `dnakit.genomics.export_cadnano_v2_design(design: 'Mapping', *, name: 'str' = '', whitespace: 'bool' = True, output_path: 'str | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `name`, `whitespace`, optional new `output_path`.
- **Output:** `json_text`, `cadnano_v2`, `output_path`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    export_cadnano_v2_design,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = export_cadnano_v2_design(design, name="rect")
print(sorted(r["cadnano_v2"])[:4], len(r["cadnano_v2"]["vstrands"]))
```

- **Example output:**

```text
['name', 'vstrands'] 4
```

## 5) `dnakit.genomics.export_oxview_design()`

- **Purpose:** Export a design as an oxView file.
- **Method:** Calls `Design.to_oxview_format` (the same content as `to_oxview_json` and `write_oxview_file`) and counts strands and nucleotides in the result.
- **API:** `dnakit.genomics.export_oxview_design(design: 'Mapping', *, use_strand_colors: 'bool' = True, output_path: 'str | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `use_strand_colors`, optional new `output_path`.
- **Output:** `json_text`, `strands`, `nucleotides`, `output_path`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_oxview_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = export_oxview_design(design)
print(r["strands"], r["nucleotides"])
```

- **Example output:**

```text
2 72
```

## 6) `dnakit.genomics.export_idt_order()`

- **Purpose:** Export the strands as an IDT bulk-input text or as a plate Excel workbook.
- **Method:** `order_format="bulk_text"` calls `Design.to_idt_bulk_input_format`; `"plate_excel"` calls `write_idt_plate_excel_file` and needs openpyxl and a new `.xlsx` path. Modification codes are written into the sequences.
- **API:** `dnakit.genomics.export_idt_order(design: 'Mapping', *, order_format: 'str' = 'bulk_text', output_path: 'str | None' = None, delimiter: 'str' = ',', domain_delimiter: 'str' = '', only_strands_with_vendor_fields: 'bool' = False, export_scaffold: 'bool' = False, export_non_modified_strand_version: 'bool' = False, plate_type: 'int' = 96) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `order_format`, delimiters, scaffold and vendor-field switches, `plate_type` 96 or 384.
- **Output:** Bulk: `text`, `lines`. Plate: `output_path`, `bytes`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_idt_order,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = export_idt_order(design)
print(r["text"])
print(r["lines"])
```

- **Example output:**

```text
stapA,ACGTACGTACGTACGTACGTACGT,25nm,STD
1
```

## 7) `dnakit.genomics.build_plate_maps()`

- **Purpose:** Build 96 or 384 well synthesis plate maps.
- **Method:** Calls `Design.plate_maps` and `PlateMap.to_table`. Only strands whose vendor fields name a plate and a well are placed.
- **API:** `dnakit.genomics.build_plate_maps(design: 'Mapping', *, plate_type: 'int' = 96, strand_indices: 'Sequence[int] | None' = None, table_format: 'str' = 'pipe') -> 'dict[str, Any]'`.
- **Input:** A design mapping, `plate_type`, optional `strand_indices`, `table_format`.
- **Output:** `plates`: list with `plate_name`, `wells` (well to strand name) and the printable `table`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_plate_maps,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = build_plate_maps(design)
print(r["plates"][0]["plate_name"], r["plates"][0]["wells"])
```

- **Example output:**

```text
P1 {'A1': 'stapA'}
```

## 8) `dnakit.genomics.get_strand_sequence()`

- **Purpose:** Return the DNA sequence of a strand or of one of its domains, with domain delimiters or vendor formatting.
- **Method:** `Strand.dna_sequence_delimited`, `Strand.vendor_dna_sequence` and the `dna_sequence` / `vendor_dna_sequence` of a domain, loopout or extension, chosen by the arguments. Unassigned bases are `?`.
- **API:** `dnakit.genomics.get_strand_sequence(design: 'Mapping', strand_index: 'int', *, domain_index: 'int | None' = None, delimiter: 'str' = '', vendor_format: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, optional `domain_index`, `delimiter`, `vendor_format`.
- **Output:** `sequence` (or `None` when no DNA is assigned).
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_strand_sequence,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(get_strand_sequence(design, 0, delimiter="|")["sequence"])
print(get_strand_sequence(design, 0, domain_index=1)["sequence"])
print(get_strand_sequence(design, 1, vendor_format=True)["sequence"])
```

- **Example output:**

```text
ACGTACGTACGTACGTACGTACGT|ACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGT
```

## 9) `dnakit.genomics.get_domain_subsequence()`

- **Purpose:** Return the DNA subsequence of a domain between two helix offsets.
- **Method:** Calls `Domain.dna_sequence_in`; both offsets are inclusive and deletions and insertions are respected.
- **API:** `dnakit.genomics.get_domain_subsequence(design: 'Mapping', strand_index: 'int', domain_index: 'int', offset_left: 'int', offset_right: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `domain_index`, `offset_left`, `offset_right`.
- **Output:** `sequence`, the two offsets.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_domain_subsequence,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(get_domain_subsequence(design, 0, 0, 2, 5)["sequence"])
```

- **Example output:**

```text
GTAC
```

## 10) `dnakit.genomics.assign_dna_sequence()`

- **Purpose:** Assign DNA to a strand or to one domain, optionally writing the complement to the bound strands.
- **Method:** Calls `Design.assign_dna` (this covers `Strand.set_dna_sequence` and the builder's `with_sequence`). `?` is the unknown-base wildcard.
- **API:** `dnakit.genomics.assign_dna_sequence(design: 'Mapping', sequence: 'str', *, strand_index: 'int', domain_index: 'int | None' = None, assign_complement: 'bool' = True, check_length: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `sequence`, `strand_index`, optional `domain_index`, `assign_complement`, `check_length`.
- **Output:** The new `design`, counts and the `assigned` strand summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import assign_dna_sequence, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)
print(r["assigned"]["dna_length"], r["strand_lengths"])
print(r["design"]["strands"][1]["sequence"])
r = assign_dna_sequence(design, "AAAA", strand_index=0, domain_index=0, assign_complement=False)
print(r["design"]["strands"][0]["sequence"])
```

- **Example output:**

```text
48 [48, 24]
ACGTACGTACGTACGTACGTACGT
AAAA????????????????????????????????????????????
```

## 11) `dnakit.genomics.assign_complement_sequence()`

- **Purpose:** Give a strand the Watson-Crick complement of the strand it is bound to.
- **Method:** Calls `Strand.assign_dna_complement_from`; only the bound regions get bases.
- **API:** `dnakit.genomics.assign_complement_sequence(design: 'Mapping', strand_index: 'int', other_strand_index: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index` (receiver), `other_strand_index` (source).
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_complement_sequence,
    assign_dna_sequence,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0, assign_complement=False)["design"]
r = assign_complement_sequence(design, 1, 0)
print(r["assigned"]["dna_sequence"])
```

- **Example output:**

```text
ACGTACGTACGTACGTACGTACGT
```

## 12) `dnakit.genomics.assign_helix_domain_sequences()`

- **Purpose:** Assign a strand the concatenated DNA of the domains on one helix.
- **Method:** Calls `Design.assign_dna_sequence`: the domains of the opposite direction on the helix are joined in helix order and written as they are.
- **API:** `dnakit.genomics.assign_helix_domain_sequences(design: 'Mapping', strand_index: 'int', *, helix: 'int', forward: 'bool') -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `helix`, `forward` (direction of the receiving strand).
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    assign_helix_domain_sequences,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0, assign_complement=False)["design"]
r = assign_helix_domain_sequences(design, 1, helix=0, forward=False)
print(r["assigned"]["dna_sequence"])
```

- **Example output:**

```text
ACGTACGTACGTACGTACGTACGT
```

## 13) `dnakit.genomics.check_complementarity()`

- **Purpose:** Test Watson-Crick complementarity of two bases or two sequences.
- **Method:** `mode="bases"` calls `bases_complementary`; `mode="sequences"` calls `reverse_complementary`, which compares the two sequences position by position (scadnano's definition).
- **API:** `dnakit.genomics.check_complementarity(first: 'str', second: 'str', *, mode: 'str' = 'bases', allow_wildcard: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** `first`, `second`, `mode`, `allow_wildcard`.
- **Output:** `complementary` (bool) and `mode`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import check_complementarity

print(check_complementarity("A", "T")["complementary"])
print(check_complementarity("A", "G")["complementary"])
print(check_complementarity("ACG", "TGC", mode="sequences")["complementary"])
print(check_complementarity("A?", "TG", mode="sequences", allow_wildcard=True)["complementary"])
```

- **Example output:**

```text
True
False
True
True
```

## 14) `dnakit.genomics.get_m13_sequence()`

- **Purpose:** Return the M13mp18 scaffold sequence.
- **Method:** Calls `scadnano.m13` with the rotation (default 5587) and variant (`p7249`, `p7560`, `p8064`, `p8634`).
- **API:** `dnakit.genomics.get_m13_sequence(*, rotation: 'int' = 5587, variant: 'str' = 'p7249') -> 'dict[str, Any]'`.
- **Input:** `rotation`, `variant`.
- **Output:** `sequence`, `length`, `variant`, `rotation`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import get_m13_sequence

r = get_m13_sequence()
print(r["length"], r["sequence"][:20])
print(get_m13_sequence(rotation=0, variant="p8064")["length"])
```

- **Example output:**

```text
7249 TTCCCTTCCTTTCTCGCCAC
8064
```

## 15) `dnakit.genomics.assign_m13_scaffold()`

- **Purpose:** Assign M13mp18 to the scaffold strand and complements to the staples.
- **Method:** Calls `Design.assign_m13_to_scaffold`.
- **API:** `dnakit.genomics.assign_m13_scaffold(design: 'Mapping', *, rotation: 'int' = 5587, variant: 'str' = 'p7249') -> 'dict[str, Any]'`.
- **Input:** A design mapping with one scaffold, `rotation`, `variant`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import assign_m13_scaffold, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = assign_m13_scaffold(design)
print(r["design"]["strands"][0]["sequence"][:20], r["design"]["strands"][1]["sequence"][:20])
```

- **Example output:**

```text
TTCCCTTCCTTTCTCGCCAC TAGAGCTTGACGGGGAAAGC
```

## 16) `dnakit.genomics.get_strand_without_modifications()`

- **Purpose:** Return a copy of a strand with its modifications removed; the design is untouched.
- **Method:** Calls `Strand.no_modifications_version`.
- **API:** `dnakit.genomics.get_strand_without_modifications(design: 'Mapping', strand_index: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and `strand_index`.
- **Output:** `strand` (JSON mapping) and `dna_sequence`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_strand_without_modifications,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")["design"]
print(design["strands"][1].get("5prime_modification"))
r = get_strand_without_modifications(design, 1)
print(r["strand"].get("5prime_modification"), r["dna_sequence"])
```

- **Example output:**

```text
/5Biosg/
None ACGTACGTACGTACGTACGTACGT
```

## 17) `dnakit.genomics.inline_deletions_insertions()`

- **Purpose:** Fold all deletions and insertions into the domain lengths.
- **Method:** Calls `Design.inline_deletions_insertions`; the number of bases stays the same.
- **API:** `dnakit.genomics.inline_deletions_insertions(design: 'Mapping', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping.
- **Output:** The new `design` and `strand_lengths_before`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_insertion,
    build_strand,
    inline_deletions_insertions,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = add_design_insertion(design, 0, 10, 2)["design"]
print(design["strands"][0]["domains"][0])
r = inline_deletions_insertions(design)
print(r["design"]["strands"][0]["domains"][0], r["strand_lengths"] == r["strand_lengths_before"])
```

- **Example output:**

```text
{'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'insertions': [[10, 2]]}
{'helix': 0, 'forward': True, 'start': 0, 'end': 26} True
```

## 18) `dnakit.genomics.set_strand_modification()`

- **Purpose:** Attach a 5', 3' or internal DNA modification to a strand.
- **Method:** Calls `Strand.set_modification_5p`, `set_modification_3p` or `set_modification_internal` according to `location`.
- **API:** `dnakit.genomics.set_strand_modification(design: 'Mapping', strand_index: 'int', *, location: 'str', display_text: 'str', vendor_code: 'str', connector_length: 'int' = 4, index: 'int | None' = None, allowed_bases: 'Sequence[str] | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `location` (`5p`, `3p`, `internal`), `display_text`, `vendor_code`, `connector_length`, `index` and `allowed_bases` for internal ones.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")
print(r["design"]["strands"][1]["5prime_modification"])
r = set_strand_modification(r["design"], 1, location="internal", display_text="T*", vendor_code="/iBiodT/", index=3, allowed_bases=["T"])
print(r["design"]["strands"][1]["internal_modifications"])
```

- **Example output:**

```text
/5Biosg/
{'3': '/iBiodT/'}
```

## 19) `dnakit.genomics.remove_strand_modification()`

- **Purpose:** Remove the 5', 3' or an internal modification from a strand.
- **Method:** Calls `Strand.remove_modification_5p`, `remove_modification_3p` or `remove_modification_internal`.
- **API:** `dnakit.genomics.remove_strand_modification(design: 'Mapping', strand_index: 'int', *, location: 'str', index: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `location`, `index` for internal ones.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    remove_strand_modification,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="3p", display_text="Q", vendor_code="/3IABkFQ/")["design"]
print("3prime_modification" in design["strands"][1])
r = remove_strand_modification(design, 1, location="3p")
print("3prime_modification" in r["design"]["strands"][1])
```

- **Example output:**

```text
True
False
```

## 20) `dnakit.genomics.list_design_modifications()`

- **Purpose:** List the distinct DNA modifications used in a design.
- **Method:** Calls `Design.modifications`, optionally for one modification type.
- **API:** `dnakit.genomics.list_design_modifications(design: 'Mapping', *, location: 'str | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping and optional `location`.
- **Output:** `modifications` (list of mappings) and `count`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    list_design_modifications,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")["design"]
r = list_design_modifications(design)
print(r["count"], r["modifications"])
print(list_design_modifications(design, location="3p")["count"])
```

- **Example output:**

```text
1 [{'display_text': 'B', 'vendor_code': '/5Biosg/', 'location': "5'"}]
0
```

## 21) `dnakit.genomics.compute_base_pairs()`

- **Purpose:** Compute the base pairs of a design.
- **Method:** Calls `Design.base_pairs` (Watson-Crick pairs, or any pairs with `allow_mismatches=True`).
- **API:** `dnakit.genomics.compute_base_pairs(design: 'Mapping', *, allow_mismatches: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** A design mapping with assigned sequences, `allow_mismatches`.
- **Output:** `base_pairs` (helix index to offsets) and `total`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    compute_base_pairs,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = compute_base_pairs(design)
print(r["total"], r["base_pairs"]["1"][:5])
```

- **Example output:**

```text
24 [0, 1, 2, 3, 4]
```
