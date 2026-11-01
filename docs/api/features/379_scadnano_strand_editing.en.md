# scadnano designs: strand and domain editing

Functions that add or remove helices and strands, build strands domain by domain, edit domain lists, add crossovers, deletions, insertions, nicks and ligations, and shift the whole design.

These functions are wrappers of [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package) (MIT license, tested with version 0.21.1), a Python package for designing DNA nanostructures. scadnano is optional and is imported only when a function runs; install it with `pip install "dnakit[dna-nanotechnology]"` (the Excel order also needs openpyxl and the plate tables need tabulate, which scadnano installs). A design is a scadnano JSON mapping (the content of a `.sc` file). Editing functions never change their argument: they return the new design under `design`. Strands are addressed by their 0-based position in `design["strands"]` and domains by their position inside the strand; helix indices and offsets are scadnano's own. The functions do not predict folding, yield or stability. A design with sequences that were assigned before an insertion, deletion or domain change keeps the old sequence, so assign sequences last.

## 1) `dnakit.genomics.add_design_helix()`

- **Purpose:** Add a helix to a design.
- **Method:** Calls `Design.add_helix` and sets the helix index so that the saved JSON reloads (scadnano leaves it at -1).
- **API:** `dnakit.genomics.add_design_helix(design: 'Mapping', *, idx: 'int', max_offset: 'int', min_offset: 'int' = 0, grid_position: 'Sequence[int] | None' = None, roll: 'float' = 0.0, group: 'str' = 'default_group') -> 'dict[str, Any]'`.
- **Input:** A design mapping, `idx`, `max_offset`, `min_offset`, `grid_position`, `roll`, `group`.
- **Output:** The new `design` and counts.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_helix,
    assign_dna_sequence,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = add_design_helix(design, idx=2, max_offset=64, grid_position=[0, 2])
print(r["helices"], r["design"]["helices"][2])
```

- **Example output:**

```text
3 {'max_offset': 64, 'grid_position': [0, 2]}
```

## 2) `dnakit.genomics.add_design_strand()`

- **Purpose:** Add a strand, given as scadnano JSON, to a design.
- **Method:** Appends the mapping to the strand list and reloads the design, which validates it (`Design.add_strand`).
- **API:** `dnakit.genomics.add_design_strand(design: 'Mapping', strand: 'Mapping', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and a strand mapping with a `domains` list.
- **Output:** The new `design`, `added_strand_index`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_strand, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_strand(design, {"domains": [{"helix": 0, "forward": True, "start": 30, "end": 40}], "name": "extra"})
print(r["added_strand_index"], r["strand_names"])
```

- **Example output:**

```text
2 [None, 'stapA', 'extra']
```

## 3) `dnakit.genomics.build_strand()`

- **Purpose:** Build a strand domain by domain (crossovers, loopouts, extensions, deletions, insertions, sequence, scaffold, circular, vendor fields).
- **Method:** Replays the steps on scadnano's `StrandBuilder` obtained from `Design.draw_strand`; only the listed operations are allowed and unknown keys are refused.
- **API:** `dnakit.genomics.build_strand(design: 'Mapping', *, helix: 'int', offset: 'int', steps: 'Sequence[Mapping]') -> 'dict[str, Any]'`.
- **Input:** A design mapping, start `helix` and `offset`, and `steps`, a list of mappings with `op` (`to`, `move`, `cross`, `loopout`, `extension_3p`, `extension_5p`, `with_deletions`, `with_insertions`, `with_sequence`, `with_domain_sequence`, `with_name`, `with_label`, `with_domain_name`, `as_scaffold`, `as_circular`, `with_vendor_fields`).
- **Output:** The new `design`, `added_strand_index` and the `strand` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

steps = [{"op": "to", "offset": 40}, {"op": "loopout", "helix": 1, "length": 3}, {"op": "move", "delta": -16}, {"op": "with_deletions", "deletions": [30]}, {"op": "with_name", "name": "hairpin"}]
r = build_strand(design, helix=0, offset=24, steps=steps)
print(r["added_strand_index"], r["strand"]["name"], r["strand"]["domains"], r["strand"]["dna_length"])
```

- **Example output:**

```text
2 hairpin 3 34
```

## 4) `dnakit.genomics.remove_design_strand()`

- **Purpose:** Remove a strand from a design.
- **Method:** Calls `Design.remove_strand`.
- **API:** `dnakit.genomics.remove_design_strand(design: 'Mapping', strand_index: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and `strand_index`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    remove_design_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = remove_design_strand(design, 1)
print(r["strands"], r["removed_strand_index"])
```

- **Example output:**

```text
1 1
```

## 5) `dnakit.genomics.find_strand_by_name()`

- **Purpose:** Find a strand by name.
- **Method:** Searches the strands like `Design.strand_with_name`; a missing name returns `found=False` instead of an error.
- **API:** `dnakit.genomics.find_strand_by_name(design: 'Mapping', name: 'str', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and `name`.
- **Output:** `found` and the strand summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_strand_by_name,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(find_strand_by_name(design, "stapA")["strand_index"])
print(find_strand_by_name(design, "nope")["found"])
```

- **Example output:**

```text
1
False
```

## 6) `dnakit.genomics.list_strands_at_helix_end()`

- **Purpose:** List the strands whose 5' or 3' end lies on a helix.
- **Method:** Calls `Design.strands_starting_on_helix` (`end="5p"`) or `strands_ending_on_helix` (`end="3p"`).
- **API:** `dnakit.genomics.list_strands_at_helix_end(design: 'Mapping', helix: 'int', *, end: 'str' = '5p') -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `end`.
- **Output:** `strands`: summaries with `strand_index`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    list_strands_at_helix_end,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print([s["strand_index"] for s in list_strands_at_helix_end(design, 1, end="5p")["strands"]])
print([s["strand_index"] for s in list_strands_at_helix_end(design, 1, end="3p")["strands"]])
```

- **Example output:**

```text
[1]
[0, 1]
```

## 7) `dnakit.genomics.find_domains_at_offset()`

- **Purpose:** Find the domains that cover a helix offset.
- **Method:** Calls `Design.domains_at`, or `Design.domain_at` when `forward` is given.
- **API:** `dnakit.genomics.find_domains_at_offset(design: 'Mapping', helix: 'int', offset: 'int', *, forward: 'bool | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`, optional `forward`.
- **Output:** `domains`: rows with strand and domain index, helix, direction, start and end.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_domains_at_offset,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

rows = find_domains_at_offset(design, 1, 10)["domains"]
print([(d["strand_index"], d["domain_index"], d["forward"]) for d in rows])
print(find_domains_at_offset(design, 1, 10, forward=False)["domains"][0]["strand_index"])
```

- **Example output:**

```text
[(0, 1, False), (1, 0, True)]
0
```

## 8) `dnakit.genomics.insert_strand_domain()`

- **Purpose:** Insert (or append) a domain, loopout or extension into a strand.
- **Method:** Calls `Design.insert_domain` (or `Design.append_domain` when `order` is omitted).
- **API:** `dnakit.genomics.insert_strand_domain(design: 'Mapping', strand_index: 'int', domain: 'Mapping', *, order: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, a domain mapping, optional `order`.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    insert_strand_domain,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = insert_strand_domain(design, 1, {"helix": 1, "forward": True, "start": 30, "end": 38})
print(r["assigned"]["domains"], r["design"]["strands"][1]["domains"][-1])
r = insert_strand_domain(design, 1, {"loopout": 4}, order=0)
print(r["design"]["strands"][1]["domains"][0])
```

- **Example output:**

```text
2 {'helix': 1, 'forward': True, 'start': 30, 'end': 38}
{'loopout': 4}
```

## 9) `dnakit.genomics.remove_strand_domain()`

- **Purpose:** Remove a domain from a strand.
- **Method:** Calls `Design.remove_domain`.
- **API:** `dnakit.genomics.remove_strand_domain(design: 'Mapping', strand_index: 'int', domain_index: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `domain_index`.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    remove_strand_domain,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = remove_strand_domain(design, 3, 0)
print(design["strands"][3]["domains"])
print(r["design"]["strands"][3]["domains"])
```

- **Example output:**

```text
[{'helix': 0, 'forward': False, 'start': 160, 'end': 176}, {'helix': 1, 'forward': True, 'start': 160, 'end': 176}]
[{'helix': 1, 'forward': True, 'start': 160, 'end': 176}]
```

## 10) `dnakit.genomics.set_strand_domains()`

- **Purpose:** Replace the whole domain list of a strand.
- **Method:** Calls `Strand.set_domains` with domains built from the mappings.
- **API:** `dnakit.genomics.set_strand_domains(design: 'Mapping', strand_index: 'int', domains: 'Sequence[Mapping]', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, a list of domain mappings.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    set_strand_domains,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = set_strand_domains(design, 3, [{"helix": 3, "forward": True, "start": 16, "end": 24}])
print(r["assigned"]["domains"], r["assigned"]["dna_length"])
```

- **Example output:**

```text
1 8
```

## 11) `dnakit.genomics.rotate_strand_domains()`

- **Purpose:** Rotate the order of the domains of a strand.
- **Method:** Calls `Strand.rotate_domains`.
- **API:** `dnakit.genomics.rotate_strand_domains(design: 'Mapping', strand_index: 'int', rotation: 'int', *, forward: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `rotation`, `forward`.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    rotate_strand_domains,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = rotate_strand_domains(design, 3, 1)
print([d["helix"] for d in design["strands"][3]["domains"]], [d["helix"] for d in r["design"]["strands"][3]["domains"]])
```

- **Example output:**

```text
[0, 1] [1, 0]
```

## 12) `dnakit.genomics.reverse_strands()`

- **Purpose:** Reverse the polarity of one strand or of all strands.
- **Method:** Calls `Strand.reverse` for one strand, or `Design.reverse_all` when `strand_index` is `None`.
- **API:** `dnakit.genomics.reverse_strands(design: 'Mapping', *, strand_index: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping and optional `strand_index`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    reverse_strands,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = reverse_strands(design, strand_index=3)
print(design["strands"][3]["domains"][0]["helix"], r["design"]["strands"][3]["domains"][0]["helix"])
print(reverse_strands(design)["strands"])
```

- **Example output:**

```text
0 1
21
```

## 13) `dnakit.genomics.set_strand_properties()`

- **Purpose:** Mark a strand as scaffold, make it circular or linear, or set its vendor fields.
- **Method:** Calls `Strand.set_scaffold`, `set_circular` / `set_linear` and sets `VendorFields`; arguments left as `None` are not changed.
- **API:** `dnakit.genomics.set_strand_properties(design: 'Mapping', strand_index: 'int', *, scaffold: 'bool | None' = None, circular: 'bool | None' = None, vendor_fields: 'Mapping[str, str | None] | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `scaffold`, `circular`, `vendor_fields`.
- **Output:** The new `design` and the `assigned` summary.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    set_strand_properties,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = set_strand_properties(design, 1, scaffold=True, vendor_fields={"scale": "100nm", "purification": "HPLC"})
print(r["assigned"]["scaffold"], r["design"]["strands"][1]["vendor"])
r = set_strand_properties(design, 0, circular=False)
print(r["assigned"]["circular"])
```

- **Example output:**

```text
True {'scale': '100nm', 'purification': 'HPLC'}
False
```

## 14) `dnakit.genomics.add_crossover()`

- **Purpose:** Add a half crossover or a full crossover between two helices.
- **Method:** Calls `Design.add_half_crossover`, or `Design.add_full_crossover` (two half crossovers) when `full=True`.
- **API:** `dnakit.genomics.add_crossover(design: 'Mapping', *, helix: 'int', helix2: 'int', offset: 'int', forward: 'bool', offset2: 'int | None' = None, forward2: 'bool | None' = None, full: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `helix2`, `offset`, `forward`, optional `offset2`, `forward2`, `full`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import add_crossover, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
design = build_strand(empty, helix=0, offset=0, steps=[{"op": "to", "offset": 24}])["design"]
design = build_strand(design, helix=1, offset=24, steps=[{"op": "move", "delta": -24}])["design"]

r = add_crossover(design, helix=0, helix2=1, offset=23, forward=True)
print(len(design["strands"]), "->", r["strands"], r["design"]["strands"][0]["domains"])
r = add_crossover(design, helix=0, helix2=1, offset=10, forward=True, full=True)
print(r["strands"])
```

- **Example output:**

```text
2 -> 1 [{'helix': 0, 'forward': True, 'start': 0, 'end': 24}, {'helix': 1, 'forward': False, 'start': 0, 'end': 24}]
2
```

## 15) `dnakit.genomics.add_design_deletion()`

- **Purpose:** Add a single-base deletion at a helix offset.
- **Method:** Calls `Design.add_deletion` (also what the builder's `with_deletions` does).
- **API:** `dnakit.genomics.add_design_deletion(design: 'Mapping', helix: 'int', offset: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_deletion, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_deletion(design, 0, 10)
print(r["strand_lengths"], r["design"]["strands"][0]["domains"][0])
```

- **Example output:**

```text
[47, 24] {'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'deletions': [10]}
```

## 16) `dnakit.genomics.add_design_insertion()`

- **Purpose:** Add an insertion of extra bases at a helix offset.
- **Method:** Calls `Design.add_insertion` (also the builder's `with_insertions`).
- **API:** `dnakit.genomics.add_design_insertion(design: 'Mapping', helix: 'int', offset: 'int', length: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`, `length`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_insertion, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_insertion(design, 0, 10, 2)
print(r["strand_lengths"], r["design"]["strands"][0]["domains"][0])
```

- **Example output:**

```text
[50, 24] {'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'insertions': [[10, 2]]}
```

## 17) `dnakit.genomics.add_design_nick()`

- **Purpose:** Nick a domain, splitting its strand.
- **Method:** Calls `Design.add_nick`; the cut lies between `offset - 1` and `offset`.
- **API:** `dnakit.genomics.add_design_nick(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *, new_color: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`, `forward`, `new_color`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_nick, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_nick(design, 0, 12, True)
print(r["strands"], r["strand_lengths"])
```

- **Example output:**

```text
3 [24, 12, 36]
```

## 18) `dnakit.genomics.ligate_design_nick()`

- **Purpose:** Join two strands that meet in a nick.
- **Method:** Calls `Design.ligate`.
- **API:** `dnakit.genomics.ligate_design_nick(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`, `forward`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_nick,
    build_strand,
    ligate_design_nick,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

nicked = add_design_nick(design, 0, 12, True)["design"]
r = ligate_design_nick(nicked, 0, 12, True)
print(r["strands"], r["strand_lengths"])
```

- **Example output:**

```text
2 [24, 48]
```

## 19) `dnakit.genomics.move_strands()`

- **Purpose:** Shift all strands along the helices or across helices.
- **Method:** `axis="offset"` calls `Design.move_strand_offsets`; `axis="helix"` calls `Design.move_strands_on_helices`.
- **API:** `dnakit.genomics.move_strands(design: 'Mapping', delta: 'int', *, axis: 'str' = 'offset') -> 'dict[str, Any]'`.
- **Input:** A design mapping, `delta`, `axis`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    move_strands,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = move_strands(design, 5)
print(design["strands"][0]["domains"][0], r["design"]["strands"][0]["domains"][0])
r = move_strands(design, 0, axis="helix")
print(r["strands"])
```

- **Example output:**

```text
{'helix': 0, 'forward': True, 'start': 0, 'end': 24} {'helix': 0, 'forward': True, 'start': 5, 'end': 29}
2
```
