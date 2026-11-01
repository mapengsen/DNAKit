# scadnano designs: overlaps, crossovers, helix geometry and DNA origami

Functions that test overlaps of domains, list crossovers and potential crossover sites, run autostaple and roll relaxation, report backbone angles, and generate a DNA origami rectangle with twist-correction deletions.

These functions are wrappers of [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package) (MIT license, tested with version 0.21.1), a Python package for designing DNA nanostructures. scadnano is optional and is imported only when a function runs; install it with `pip install "dnakit[dna-nanotechnology]"` (the Excel order also needs openpyxl and the plate tables need tabulate, which scadnano installs). A design is a scadnano JSON mapping (the content of a `.sc` file). Editing functions never change their argument: they return the new design under `design`. Strands are addressed by their 0-based position in `design["strands"]` and domains by their position inside the strand; helix indices and offsets are scadnano's own. The functions do not predict folding, yield or stability. A design with sequences that were assigned before an insertion, deletion or domain change keeps the old sequence, so assign sequences last.

## 1) `dnakit.genomics.check_domain_overlap()`

- **Purpose:** Test whether two domains or two strands overlap on a helix, and get the overlap interval.
- **Method:** Calls `Domain.overlaps` (opposite directions), `Domain.compute_overlap` and `Domain.overlaps_illegally` (same direction), or `Strand.overlaps` for whole strands.
- **API:** `dnakit.genomics.check_domain_overlap(design: 'Mapping', strand_index: 'int', other_strand_index: 'int', *, domain_index: 'int | None' = None, other_domain_index: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, two strand indices and optionally two domain indices.
- **Output:** `overlaps`, `overlap_interval` (start inclusive, end exclusive), `overlaps_illegally`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    check_domain_overlap,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = check_domain_overlap(design, 0, 0, domain_index=0, other_domain_index=0)
print(r["overlap_interval"], r["overlaps"], r["overlaps_illegally"])
print(check_domain_overlap(design, 0, 1)["overlaps"])
```

- **Example output:**

```text
[0, 24] False True
True
```

## 2) `dnakit.genomics.find_overlapping_domains()`

- **Purpose:** Find the domains whose offset range overlaps that of a given domain.
- **Method:** Calls `Domain.find_overlapping_ranges` on the other domains of the design. It compares ranges only, not helices or directions.
- **API:** `dnakit.genomics.find_overlapping_domains(design: 'Mapping', strand_index: 'int', domain_index: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `domain_index`.
- **Output:** `domains`: located rows.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_overlapping_domains,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = find_overlapping_domains(design, 0, 0)
print([(d["strand_index"], d["domain_index"]) for d in r["domains"]])
```

- **Example output:**

```text
[(0, 1), (1, 0)]
```

## 3) `dnakit.genomics.find_helix_domain_overlap_pairs()`

- **Purpose:** List forward/reverse domain pairs that overlap on a helix.
- **Method:** Calls `scadnano.find_overlapping_domains_on_helix` for one or all helices.
- **API:** `dnakit.genomics.find_helix_domain_overlap_pairs(design: 'Mapping', *, helix: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping and optional `helix`.
- **Output:** `pairs` and `count`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_helix_domain_overlap_pairs,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = find_helix_domain_overlap_pairs(design)
print(r["count"], r["pairs"])
```

- **Example output:**

```text
1 [{'helix': 1, 'first': {'strand_index': 1, 'domain_index': 0}, 'second': {'strand_index': 0, 'domain_index': 1}}]
```

## 4) `dnakit.genomics.domain_has_crossover()`

- **Purpose:** Test whether a domain has a crossover at an offset.
- **Method:** Calls `Domain.has_crossover_at`.
- **API:** `dnakit.genomics.domain_has_crossover(design: 'Mapping', strand_index: 'int', domain_index: 'int', offset: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `domain_index`, `offset`.
- **Output:** `has_crossover`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    domain_has_crossover,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(domain_has_crossover(design, 0, 0, 23)["has_crossover"])
print(domain_has_crossover(design, 0, 0, 10)["has_crossover"])
```

- **Example output:**

```text
True
False
```

## 5) `dnakit.genomics.map_domain_offset_to_dna_index()`

- **Purpose:** Map domain positions to indices in the strand DNA sequence.
- **Method:** Calls `Strand.dna_index_start_domain` (equal to `Domain.get_seq_start_idx` and `Loopout.get_seq_start_idx`) and, with `offset`, `Domain.domain_offset_to_strand_dna_idx`.
- **API:** `dnakit.genomics.map_domain_offset_to_dna_index(design: 'Mapping', strand_index: 'int', domain_index: 'int', *, offset: 'int | None' = None, offset_closer_to_5p: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `strand_index`, `domain_index`, optional `offset`, `offset_closer_to_5p`.
- **Output:** `start_index` and, with `offset`, `dna_index`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    map_domain_offset_to_dna_index,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = map_domain_offset_to_dna_index(design, 0, 1, offset=10)
print(r["start_index"], r["dna_index"])
```

- **Example output:**

```text
24 37
```

## 6) `dnakit.genomics.list_helix_crossovers()`

- **Purpose:** List the crossovers incident to a helix.
- **Method:** Calls `Helix.crossover_addresses`.
- **API:** `dnakit.genomics.list_helix_crossovers(design: 'Mapping', helix: 'int', *, allow_intrahelix: 'bool' = True, allow_intergroup: 'bool' = True) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `allow_intrahelix`, `allow_intergroup`.
- **Output:** `crossovers`: other helix, offset, forward; `count`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_helix_crossovers,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_helix_crossovers(design, 0)
print(r["count"], r["crossovers"][:3])
```

- **Example output:**

```text
10 [{'helix': 1, 'offset': 16, 'forward': True}, {'helix': 1, 'offset': 175, 'forward': True}, {'helix': 1, 'offset': 31, 'forward': False}]
```

## 7) `dnakit.genomics.split_helix_strands()`

- **Purpose:** Split the domains on a helix into staple and scaffold domains.
- **Method:** Calls `Helix.get_strand_sets`, cadnano's parity rule.
- **API:** `dnakit.genomics.split_helix_strands(design: 'Mapping', helix: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and `helix`.
- **Output:** `staple_domains` and `scaffold_domains`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    split_helix_strands,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = split_helix_strands(design, 0)
print(len(r["staple_domains"]), len(r["scaffold_domains"]))
```

- **Example output:**

```text
9 1
```

## 8) `dnakit.genomics.list_crossover_points()`

- **Purpose:** Return cadnano's crossover-point tables for the design grid.
- **Method:** Calls `Design.get_crossover_points` (square and honeycomb grids).
- **API:** `dnakit.genomics.list_crossover_points(design: 'Mapping', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping.
- **Output:** `scaffold_low`, `scaffold_high`, `staple_low`, `staple_high`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_crossover_points,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_crossover_points(design)
print(r["grid"], r["scaffold_low"][0], r["staple_low"])
```

- **Example output:**

```text
square [4, 26, 15] [[31], [23], [15], [7]]
```

## 9) `dnakit.genomics.list_potential_crossovers()`

- **Purpose:** List the offsets where a crossover to a neighbouring helix is possible.
- **Method:** Calls `Design.potential_crossover_list`.
- **API:** `dnakit.genomics.list_potential_crossovers(design: 'Mapping', helix: 'int', *, index: 'int | None' = None) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, optional `index`.
- **Output:** `sites`: neighbour helix, offset, `scaffold` or `staple`, `is_low`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_potential_crossovers,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_potential_crossovers(design, 0)
print(r["count"], r["sites"][:2])
```

- **Example output:**

```text
46 [{'neighbor_helix': 1, 'offset': 18, 'strand_type': 'scaffold', 'is_low': True}, {'neighbor_helix': 1, 'offset': 28, 'strand_type': 'scaffold', 'is_low': True}]
```

## 10) `dnakit.genomics.get_helix_neighbors()`

- **Purpose:** Return the neighbouring helices of a helix on the lattice.
- **Method:** Calls `Design.get_helix_neighbors`; `None` marks an empty lattice position.
- **API:** `dnakit.genomics.get_helix_neighbors(design: 'Mapping', helix: 'int', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping and `helix`.
- **Output:** `neighbors`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    get_helix_neighbors,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

print(get_helix_neighbors(design, 0)["neighbors"])
print(get_helix_neighbors(design, 1)["neighbors"])
```

- **Example output:**

```text
[None, 1, None, None]
[None, 0, None, 2]
```

## 11) `dnakit.genomics.autostaple_design()`

- **Purpose:** Generate staples over the scaffold with cadnano's autostaple algorithm.
- **Method:** Calls `Design.autostaple`, which removes the old staples first.
- **API:** `dnakit.genomics.autostaple_design(design: 'Mapping', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping without loopouts or extensions.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    autostaple_design,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = autostaple_design(design)
print(len(design["strands"]), "->", r["strands"], r["scaffolds"])
```

- **Example output:**

```text
21 -> 9 1
```

## 12) `dnakit.genomics.relax_helix_rolls()`

- **Purpose:** Relax the helix roll angles, or only compute the roll corrections.
- **Method:** `compute_only=False` calls `Design.relax_helix_rolls`; `True` only reports `Helix.compute_relaxed_roll_delta` per helix and leaves the design unchanged.
- **API:** `dnakit.genomics.relax_helix_rolls(design: 'Mapping', *, compute_only: 'bool' = False) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `compute_only`.
- **Output:** `applied`, and `roll_deltas` or the new `design` with `rolls`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    relax_helix_rolls,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = relax_helix_rolls(design, compute_only=True)
print(r["applied"], {k: round(v, 2) for k, v in r["roll_deltas"].items()})
r = relax_helix_rolls(design)
print(r["applied"], {k: round(v, 2) for k, v in r["rolls"].items()})
```

- **Example output:**

```text
False {'0': 349.71, '1': 331.71, '2': 349.71, '3': 331.71}
True {'0': 349.71, '1': 331.71, '2': 349.71, '3': 331.71}
```

## 13) `dnakit.genomics.get_backbone_angle()`

- **Purpose:** Return the backbone angle at a helix offset.
- **Method:** Calls `Helix.backbone_angle_at_offset` with the design geometry.
- **API:** `dnakit.genomics.get_backbone_angle(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `offset`, `forward`.
- **Output:** `angle_degrees`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    get_backbone_angle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

print(round(get_backbone_angle(design, 0, 10, True)["angle_degrees"], 2))
print(round(get_backbone_angle(design, 0, 10, False)["angle_degrees"], 2))
```

- **Example output:**

```text
342.86
132.86
```

## 14) `dnakit.genomics.create_origami_rectangle()`

- **Purpose:** Generate a DNA origami rectangle design.
- **Method:** Calls `origami_rectangle.create` with M13mp18 assigned by default.
- **API:** `dnakit.genomics.create_origami_rectangle(*, num_helices: 'int', num_cols: 'int', assign_seq: 'bool' = True, seam_left_column: 'int' = -1, nick_pattern: 'str' = 'staggered', twist_correction_deletion_spacing: 'int' = 0, twist_correction_start_col: 'int' = 1, twist_correction_deletion_offset: 'int' = -1, num_flanking_columns: 'int' = 1, num_flanking_helices: 'int' = 0, custom_scaffold: 'str | None' = None, edge_staples: 'bool' = True, scaffold_nick_offset: 'int' = -1) -> 'dict[str, Any]'`.
- **Input:** `num_helices`, `num_cols` and the upstream options (`nick_pattern`, twist-correction spacing, flanking columns, `custom_scaffold`, ...).
- **Output:** The `design` and its counts.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import create_origami_rectangle

r = create_origami_rectangle(num_helices=4, num_cols=10)
print(r["helices"], r["strands"], r["scaffolds"], r["strand_lengths"][:3])
r = create_origami_rectangle(num_helices=4, num_cols=10, twist_correction_deletion_spacing=3)
print(r["strand_lengths"][0])
```

- **Example output:**

```text
4 21 1 [640, 32, 32]
628
```

## 15) `dnakit.genomics.add_twist_correction_deletions()`

- **Purpose:** Add twist-correction deletions to an origami rectangle.
- **Method:** Calls `origami_rectangle.add_twist_correction_deletions`: one deletion per helix every `deletion_spacing` columns.
- **API:** `dnakit.genomics.add_twist_correction_deletions(design: 'Mapping', *, offset_start: 'int', deletion_spacing: 'int', deletion_start_col: 'int', deletion_offset: 'int', num_helices: 'int', num_cols: 'int', num_flanking_helices: 'int' = 0) -> 'dict[str, Any]'`.
- **Input:** A design mapping and the rectangle layout parameters.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_twist_correction_deletions,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["design"]

r = add_twist_correction_deletions(design, offset_start=0, deletion_spacing=2, deletion_start_col=1, deletion_offset=-1, num_helices=4, num_cols=10)
print(sum(r["strand_lengths"]) - sum(create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["strand_lengths"]))
```

- **Example output:**

```text
-40
```

## 16) `dnakit.genomics.add_deletion_in_range()`

- **Purpose:** Add one deletion at a free offset inside a range of a helix.
- **Method:** Calls `origami_rectangle.add_deletion_in_range`; free offsets are covered by two domains and are not domain ends.
- **API:** `dnakit.genomics.add_deletion_in_range(design: 'Mapping', *, helix: 'int', start: 'int', end: 'int', deletion_offset: 'int' = -1) -> 'dict[str, Any]'`.
- **Input:** A design mapping, `helix`, `start`, `end`, `deletion_offset`.
- **Output:** The new `design`.
- **Example code:**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_deletion_in_range,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["design"]

r = add_deletion_in_range(design, helix=1, start=17, end=30)
print(sum(r["strand_lengths"]) - sum(create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["strand_lengths"]))
```

- **Example output:**

```text
-2
```
