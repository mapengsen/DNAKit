# tskit and toytree: tree sequences and tree drawing

This group of functions is a set of adapters over two optional packages. [tskit](https://github.com/tskit-dev/tskit) is the tree-sequence toolkit: `dnakit.population` uses it to draw tree sequences (SVG, HTML and text), to generate random nucleotide strings, and to mirror the commands of `python -m tskit` on a `.trees` file (info, tree listing, one text export per table, FASTA and VCF). [toytree](https://github.com/eaton-lab/toytree) (with its renderer [toyplot](https://github.com/sandialabs/toyplot)) is a Python library for phylogenetic trees: `dnakit.phylogeny` uses it to draw single trees, tree grids and tree clouds, to colour or size nodes by a numeric feature, to export the node data table, and to prune and root trees. The tskit functions import from `dnakit.population` and the toytree functions from `dnakit.phylogeny`, for example `from dnakit.population import tree_file_info`. None of them predicts or simulates data; they read, convert and draw what you give them.

tskit is MIT-licensed. toytree is GPL-3.0 (its package metadata says "GNU General Public License v3"); DNAKit imports it only inside these functions and copies none of its code. toyplot is BSD-licensed.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises `BackendUnavailableError` with the install command (`dnakit[bio-tskit]` or `dnakit[bio-toytree]`). The drawing functions never overwrite a file: the output path must be new and its folder must exist.

## Example data

The functions that read a `.trees` file (sections 6 to 17) need a binary tree-sequence file, which cannot be shown as text. Their examples therefore load the following seven small text tables with `tskit.load_text()`, add one provenance record, and dump the result to a temporary `example.trees`. The tree sequence has 4 sample nodes (two diploid individuals), length 100, two trees (`[0, 50)` and `[50, 100)`), two sites and two mutations. The tskit drawing examples (sections 1 to 5 and 18) pass the same tree sequence inline as table rows, and the toytree examples pass Newick text, so they need no file.

**Example file:** `nodes.txt`

```text
is_sample time population individual
1 0 0 0
1 0 0 0
1 0 0 1
1 0 0 1
0 1 0 -1
0 1 0 -1
0 2 0 -1
0 1 0 -1
0 1 0 -1
0 3 0 -1
```

**Example file:** `edges.txt`

```text
left right parent child
0 50 4 0
0 50 4 1
0 50 5 2
0 50 5 3
0 50 6 4
0 50 6 5
50 100 7 0
50 100 7 2
50 100 8 1
50 100 8 3
50 100 9 7
50 100 9 8
```

**Example file:** `sites.txt`

```text
position ancestral_state
10 A
60 C
```

**Example file:** `mutations.txt`

```text
site node derived_state
0 4 G
1 7 T
```

**Example file:** `individuals.txt`

```text
flags location parents
0 0.5 -1,-1
0 1.5 -1,-1
```

**Example file:** `populations.txt`

```text
metadata
POP1
```

**Example file:** `migrations.txt`

```text
source dest node time left right
0 0 4 0.5 0 50
```

## 1) `dnakit.population.draw_tree_html()`

- **Purpose:** Write tskit's HTML view of one tree of a tree sequence to a new file.
- **Method:** Build the tree sequence from explicit table rows with tskit, pick the tree with `ts.at_index(tree_index)`, and write the text of `tskit.tree_html(tree)` as UTF-8. The destination must not exist yet and its folder must exist.
- **API:** `draw_tree_html(tables: 'Mapping[str, Any]', tree_index: 'int', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tables` is a mapping with `sequence_length`, `nodes` (`time`, optional `sample`), `edges` (`left`, `right`, `parent`, `child`), `sites` (`position`, `ancestral_state`) and `mutations` (`site`, `node`, `derived_state`); `max_rows` bounds each table (at most 1,000,000). `tree_index` is an integer from 0 to `num_trees - 1`. `destination` is a new file path.
- **Output:** A dictionary with `status`, `operation`, `image_path` (the HTML file), `format` (`html`), `bytes`, `tree_index`, `interval` (`[left, right)`), `num_roots`, `num_sites` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "tree1.html"
    result = draw_tree_html(tables, 1, out, verbose=False)
    html = out.read_text(encoding="utf-8")
    print(result["format"], result["tree_index"], result["interval"])
    print(result["num_roots"], result["num_sites"], "<div>" in html)
```

- **Example output:**

```text
html 1 [50.0, 100.0]
1 1 True
```

## 2) `dnakit.population.draw_tree_sequence_html()`

- **Purpose:** Write tskit's HTML summary of a whole tree sequence to a new file.
- **Method:** Build the tree sequence from explicit table rows and write the text of `tskit.tree_sequence_html(ts)` as UTF-8. The summary lists the table sizes, provenance records and metadata. The destination must not exist yet.
- **API:** `draw_tree_sequence_html(tables: 'Mapping[str, Any]', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tables` as in `draw_tree_html()`; `destination` is a new file path; `max_rows` bounds each table.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `format` (`html`), `bytes`, `sequence_length`, `num_trees`, `num_samples`, `num_sites`, `num_mutations`, `num_provenances` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_sequence_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "summary.html"
    result = draw_tree_sequence_html(tables, out, verbose=False)
    html = out.read_text(encoding="utf-8")
    print(result["format"], result["sequence_length"], result["num_trees"])
    print(result["num_samples"], result["num_sites"], result["num_mutations"], result["num_provenances"])
    print("Provenances" in html)
```

- **Example output:**

```text
html 100.0 2
4 2 2 0
True
```

## 3) `dnakit.population.draw_tree_sequence_svg()`

- **Purpose:** Draw a tree sequence, or one of its trees, as an SVG file.
- **Method:** Build the tree sequence from explicit table rows. With `tree_index=None` it calls `TreeSequence.draw_svg(path=...)` (all trees, optionally limited by `x_lim`); with an index it calls `Tree.draw_svg(path=...)` on `ts.at_index(index)`. The file must not exist yet.
- **API:** `draw_tree_sequence_svg(tables: 'Mapping[str, Any]', destination: 'str | Path', *, tree_index: 'int | None' = None, x_lim: 'Sequence[float] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tables` as in `draw_tree_html()`; `destination` is a new file path; `tree_index` is null or an integer from 0 to `num_trees - 1`; `x_lim` is null or two numbers `[left, right]` with `0 <= left < right <= sequence_length`, allowed only for the whole sequence (with `tree_index` it raises `ConfigurationError`).
- **Output:** A dictionary with `status`, `operation`, `image_path`, `format` (`svg`), `scope` (`tree_sequence` or `tree`), `tree_index`, `x_lim`, `bytes`, `sequence_length`, `num_trees`, `intervals` (one `[left, right]` per tree), `num_samples`, `num_sites`, `num_mutations` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_sequence_svg

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    whole = draw_tree_sequence_svg(tables, Path(tmp) / "all.svg", x_lim=[0, 60], verbose=False)
    single = draw_tree_sequence_svg(tables, Path(tmp) / "tree1.svg", tree_index=1, verbose=False)
    print(whole["scope"], whole["x_lim"], whole["intervals"], whole["num_sites"])
    print(single["scope"], single["tree_index"], single["num_trees"])
    print("<svg" in (Path(tmp) / "tree1.svg").read_text(encoding="utf-8"))
```

- **Example output:**

```text
tree_sequence [0.0, 60.0] [[0.0, 50.0], [50.0, 100.0]] 2
tree 1 2
True
```

## 4) `dnakit.population.draw_variant_html()`

- **Purpose:** Write tskit's HTML view of one variant (one site) to a new file.
- **Method:** Build the tree sequence from explicit table rows, take the `site_index`-th item of `ts.variants()` (sites in position order) and write the text of `tskit.variant_html(variant)` as UTF-8. The destination must not exist yet.
- **API:** `draw_variant_html(tables: 'Mapping[str, Any]', site_index: 'int', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tables` as in `draw_tree_html()`; `site_index` is an integer from 0 to `num_sites - 1`; `destination` is a new file path.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `format` (`html`), `bytes`, `site_index`, `position`, `alleles` (null for a missing allele), `num_samples` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_variant_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "site1.html"
    result = draw_variant_html(tables, 1, out, verbose=False)
    print(result["site_index"], result["position"], result["alleles"], result["num_samples"])
    print("<div>" in out.read_text(encoding="utf-8"))
```

- **Example output:**

```text
1 60.0 ['C', 'T'] 4
True
```

## 5) `dnakit.population.tskit_random_nucleotides()`

- **Purpose:** Draw a random nucleotide string of a given length, reproducible with a seed.
- **Method:** Call `tskit.random_nucleotides(length, seed=seed)`, which draws uniformly from `ACTG` with a numpy `RandomState`. The same seed gives tskit's string exactly; it is not the same stream as DNAKit's `generate_random_dna()`.
- **API:** `tskit_random_nucleotides(*, length: 'int', seed: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `length` is an integer from 0 to 10,000,000 (keyword only); `seed` is null (fresh random string) or an integer from 0 to 2^32 - 1.
- **Output:** A dictionary with `status`, `operation`, `length`, `seed`, `alphabet` (`ACTG`), `sequence` and `provenance`.
- **Example code:**

```python
from dnakit.population import tskit_random_nucleotides

first = tskit_random_nucleotides(length=16, seed=7, verbose=False)
second = tskit_random_nucleotides(length=16, seed=7, verbose=False)
print(first["sequence"])
print(first["sequence"] == second["sequence"], first["length"], first["alphabet"])
```

- **Example output:**

```text
GACTGGGGACTGACTT
True 16 ACTG
```

## 6) `dnakit.population.export_tree_file_edges_table()`

- **Purpose:** Write the edges table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(edges=handle, precision=precision)` into the new file. This is what tskit's `edges` command prints.
- **API:** `export_tree_file_edges_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_edges_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "edges.txt"
    result = export_tree_file_edges_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
edges 12
left | right | parent | child | metadata
0.000000 | 50.000000 | 4 | 0 |
0.000000 | 50.000000 | 4 | 1 |
0.000000 | 50.000000 | 5 | 2 |
0.000000 | 50.000000 | 5 | 3 |
50.000000 | 100.000000 | 7 | 0 |
50.000000 | 100.000000 | 7 | 2 |
50.000000 | 100.000000 | 8 | 1 |
50.000000 | 100.000000 | 8 | 3 |
0.000000 | 50.000000 | 6 | 4 |
0.000000 | 50.000000 | 6 | 5 |
50.000000 | 100.000000 | 9 | 7 |
50.000000 | 100.000000 | 9 | 8 |
```

## 7) `dnakit.population.export_tree_file_fasta()`

- **Purpose:** Write the sample haplotypes of a `.trees` file as FASTA.
- **Method:** Load the file with `tskit.load()` and call `TreeSequence.write_fasta(path, wrap_width=wrap_width)`, the operation of tskit's `run_fasta()` handler. Each sample node becomes a record named `n<node id>`. Positions without a site are written as `N` because the file has no reference sequence.
- **API:** `export_tree_file_fasta(tree_sequence_path: 'str | Path', destination: 'str | Path', *, wrap_width: 'int' = 60, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes`; `destination` is a new file path; `wrap_width` is an integer of at least 0 (default 60, 0 writes each sequence on one line).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `file_format` (`fasta`), `wrap_width`, `bytes`, `num_samples`, `sequence_length`, `num_sites` and `provenance`.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_fasta

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "samples.fa"
    result = export_tree_file_fasta(tree_file, out, wrap_width=0, verbose=False)
    print(result["file_format"], result["num_samples"], result["sequence_length"], result["num_sites"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line[:12], len(line))
```

- **Example output:**

```text
fasta 4 100.0 2
>n0 3
NNNNNNNNNNGN 100
>n1 3
NNNNNNNNNNGN 100
>n2 3
NNNNNNNNNNAN 100
>n3 3
NNNNNNNNNNAN 100
```

## 8) `dnakit.population.export_tree_file_individuals_table()`

- **Purpose:** Write the individuals table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(individuals=handle, precision=precision)` into the new file. This is what tskit's `individuals` command prints.
- **API:** `export_tree_file_individuals_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_individuals_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "individuals.txt"
    result = export_tree_file_individuals_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
individuals 2
id | flags | location | parents | metadata
0 | 0 | 0.5 | -1,-1 |
1 | 0 | 1.5 | -1,-1 |
```

## 9) `dnakit.population.export_tree_file_migrations_table()`

- **Purpose:** Write the migrations table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(migrations=handle, precision=precision)` into the new file. This is what tskit's `migrations` command prints.
- **API:** `export_tree_file_migrations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_migrations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "migrations.txt"
    result = export_tree_file_migrations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
migrations 1
left | right | node | source | dest | time | metadata
0.0 | 50.0 | 4 | 0 | 0 | 0.5 |  |
```

## 10) `dnakit.population.export_tree_file_mutations_table()`

- **Purpose:** Write the mutations table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(mutations=handle, precision=precision)` into the new file. This is what tskit's `mutations` command prints.
- **API:** `export_tree_file_mutations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_mutations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "mutations.txt"
    result = export_tree_file_mutations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
mutations 2
site | node | time | derived_state | parent | metadata
0 | 4 | unknown | G | -1 |
1 | 7 | unknown | T | -1 |
```

## 11) `dnakit.population.export_tree_file_nodes_table()`

- **Purpose:** Write the nodes table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(nodes=handle, precision=precision)` into the new file. This is what tskit's `nodes` command prints.
- **API:** `export_tree_file_nodes_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_nodes_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "nodes.txt"
    result = export_tree_file_nodes_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
nodes 10
id | is_sample | time | population | individual | metadata
0 | 1 | 0.000000 | 0 | 0 |
1 | 1 | 0.000000 | 0 | 0 |
2 | 1 | 0.000000 | 0 | 1 |
3 | 1 | 0.000000 | 0 | 1 |
4 | 0 | 1.000000 | 0 | -1 |
5 | 0 | 1.000000 | 0 | -1 |
6 | 0 | 2.000000 | 0 | -1 |
7 | 0 | 1.000000 | 0 | -1 |
8 | 0 | 1.000000 | 0 | -1 |
9 | 0 | 3.000000 | 0 | -1 |
```

## 12) `dnakit.population.export_tree_file_populations_table()`

- **Purpose:** Write the populations table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(populations=handle)` into the new file. This is what tskit's `populations` command prints.
- **API:** `export_tree_file_populations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; this function has no `precision` argument.
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_populations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "populations.txt"
    result = export_tree_file_populations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
populations 1
id | metadata
0 | POP1
```

## 13) `dnakit.population.export_tree_file_provenances_table()`

- **Purpose:** Write the provenance records of a `.trees` file as table text or as readable blocks.
- **Method:** Load the file with `tskit.load()`. By default it calls `TreeSequence.dump_text(provenances=handle)` for the tabular text. With `human_readable=True` it writes one block per record in the form `id=..., timestamp=..., record=<JSON with 4-space indent>`, as tskit's `-H` option does; a record that is not valid JSON raises `ConfigurationError`.
- **API:** `export_tree_file_provenances_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, human_readable: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes`; `destination` is a new file path; `human_readable` is a boolean (default false).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (`provenances`), `row_count` (records), `bytes` and `provenance`; the human-readable form also has `human_readable` set to true.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_provenances_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    table = export_tree_file_provenances_table(tree_file, Path(tmp) / "table.txt", verbose=False)
    blocks = export_tree_file_provenances_table(
        tree_file, Path(tmp) / "blocks.txt", human_readable=True, verbose=False
    )
    print(table["row_count"], blocks["row_count"], blocks["human_readable"])
    for line in (Path(tmp) / "table.txt").read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
    print((Path(tmp) / "blocks.txt").read_text(encoding="utf8"), end="")
```

- **Example output:**

```text
1 1 True
id | timestamp | record
0 | 2024-01-01T00:00:00 | {"software": "demo"} |
id=0, timestamp=2024-01-01T00:00:00, record={
    "software": "demo"
}
```

## 14) `dnakit.population.export_tree_file_sites_table()`

- **Purpose:** Write the sites table of a `.trees` file as text.
- **Method:** Load the file with `tskit.load()` after a size check, then call `TreeSequence.dump_text(sites=handle, precision=precision)` into the new file. This is what tskit's `sites` command prints.
- **API:** `export_tree_file_sites_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000); `destination` is a new file path whose folder exists; `precision` is an integer from 0 to 16 (default 6).
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `table` (the table name), `row_count` (rows without the header line), `bytes` and `provenance`. The text file has one header line and one line per row.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_sites_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "sites.txt"
    result = export_tree_file_sites_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
sites 2
position | ancestral_state | metadata
10.000000 | A |
60.000000 | C |
```

## 15) `dnakit.population.export_tree_file_vcf()`

- **Purpose:** Write the variants of a `.trees` file as a VCF file.
- **Method:** Load the file with `tskit.load()` and call `TreeSequence.write_vcf(handle, ploidy=..., contig_id=..., allow_position_zero=...)`, the operation of tskit's `run_vcf()` handler. tskit writes one VCF column per individual (or per `ploidy` samples when the file has no individuals). A tskit error is re-raised as `ConfigurationError` and the partly written file is removed.
- **API:** `export_tree_file_vcf(tree_sequence_path: 'str | Path', destination: 'str | Path', *, ploidy: 'int | None' = None, contig_id: 'str' = '1', allow_position_zero: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes`; `destination` is a new file path; `ploidy` is null (use tskit's default) or an integer of at least 1 that divides the number of samples (tskit refuses it when the file already has individuals); `contig_id` is the CHROM text (default `1`); `allow_position_zero` is a boolean.
- **Output:** A dictionary with `status`, `operation`, `source`, `destination`, `file_format` (`vcf`), `ploidy`, `contig_id`, `allow_position_zero`, `bytes`, `num_samples`, `num_records` (lines that are not header lines) and `provenance`.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_vcf

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "variants.vcf"
    result = export_tree_file_vcf(tree_file, out, contig_id="chr1", verbose=False)
    print(result["file_format"], result["contig_id"], result["num_samples"], result["num_records"])
    lines = out.read_text(encoding="utf8").splitlines()
    for line in lines:
        if not line.startswith("##"):
            print(line.replace("\t", " | ").rstrip())
```

- **Example output:**

```text
vcf chr1 4 2
#CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | tsk_0 | tsk_1
chr1 | 10 | 0 | A | G | . | PASS | . | GT | 1|1 | 0|0
chr1 | 60 | 1 | C | T | . | PASS | . | GT | 1|0 | 1|0
```

## 16) `dnakit.population.tree_file_info()`

- **Purpose:** Summarise a `.trees` file: the tskit summary table and the table counts.
- **Method:** Load the file with `tskit.load()` after a size check and return `str(TreeSequence)` together with the counts read from the tree sequence, as tskit's `info` command does.
- **API:** `tree_file_info(tree_sequence_path: 'str | Path', *, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes` (default 100,000,000; at most 1,000,000,000).
- **Output:** A dictionary with `status`, `operation`, `source`, `summary_text` (the multi-line tskit table, which includes the in-memory size of the file), `sequence_length`, `num_trees`, `num_samples`, `num_nodes`, `num_edges`, `num_sites`, `num_mutations`, `num_individuals`, `num_populations`, `num_migrations`, `num_provenances` and `provenance`.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import tree_file_info

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    result = tree_file_info(tree_file, verbose=False)
    print(result["sequence_length"], result["num_trees"], result["num_samples"], result["num_nodes"], result["num_edges"])
    print(result["num_sites"], result["num_mutations"], result["num_individuals"], result["num_populations"])
    print(result["num_migrations"], result["num_provenances"])
    print(result["summary_text"].splitlines()[3])
```

- **Example output:**

```text
100.0 2 4 10 12
2 2 2 1
1 1
║Trees          │      2║
```

## 17) `dnakit.population.tree_file_trees()`

- **Purpose:** List the trees of a `.trees` file with index, interval and number of sites.
- **Method:** Load the file with `tskit.load()` and iterate `TreeSequence.trees()`. The interval is rounded to `precision` decimals, as in tskit's `trees` command. With `draw=True` each record also holds `Tree.draw(format="unicode")`.
- **API:** `tree_file_trees(tree_sequence_path: 'str | Path', *, precision: 'int' = 6, draw: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tree_sequence_path` is an existing `.trees` file no larger than `max_bytes`; `precision` is an integer from 0 to 16 (default 6); `draw` is a boolean (default false).
- **Output:** A dictionary with `status`, `operation`, `source`, `precision`, `num_trees`, `trees` (a list of `index`, `num_sites`, `interval` and, with `draw`, `drawing`) and `provenance`.
- **Example code:**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import tree_file_trees

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    result = tree_file_trees(tree_file, precision=2, draw=True, verbose=False)
    print(result["num_trees"], result["precision"])
    for tree in result["trees"]:
        print(tree["index"], tree["interval"], tree["num_sites"])
        print(tree["drawing"])
```

- **Example output:**

```text
2 2
0 [0.0, 50.0] 1
   6   
 ┏━┻━┓ 
 4   5 
┏┻┓ ┏┻┓
0 1 2 3
1 [50.0, 100.0] 1
   9   
 ┏━┻━┓ 
 7   8 
┏┻┓ ┏┻┓
0 2 1 3
```

## 18) `dnakit.population.draw_tree_text()`

- **Purpose:** Write one tree of a tree sequence as tskit's text drawing to a new file.
- **Method:** Build the tree sequence from explicit table rows, pick the tree with `ts.at_index(tree_index)`, draw it with `Tree.draw(format=text_format)` and write the text as UTF-8. The result also carries the drawing, the tree's Newick string and its interval.
- **API:** `draw_tree_text(tables: 'Mapping[str, Any]', tree_index: 'int', destination: 'str | Path', *, text_format: 'str' = 'unicode', max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package tskit (extra `dnakit[bio-tskit]`); no external program, internet access or model.
- **Input:** `tables` as in `draw_tree_html()`; `tree_index` is an integer from 0 to `num_trees - 1`; `destination` is a new file path; `text_format` is `unicode` (default, box characters) or `ascii`, anything else raises `ConfigurationError`.
- **Output:** A dictionary with `status`, `operation`, `image_path`, `format`, `bytes`, `tree_index`, `interval`, `newick`, `num_roots`, `num_sites`, `text` (the drawing) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_text

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    unicode_result = draw_tree_text(tables, 0, Path(tmp) / "tree0.txt", verbose=False)
    ascii_result = draw_tree_text(
        tables, 1, Path(tmp) / "tree1.txt", text_format="ascii", verbose=False
    )
    print(unicode_result["interval"], unicode_result["newick"])
    print(unicode_result["text"])
    print(ascii_result["format"], ascii_result["interval"])
    print(ascii_result["text"])
```

- **Example output:**

```text
[0.0, 50.0] ((1:1.00000000000000,2:1.00000000000000):1.00000000000000,(3:1.00000000000000,4:1.00000000000000):1.00000000000000);
   6   
 ┏━┻━┓ 
 4   5 
┏┻┓ ┏┻┓
0 1 2 3
ascii [50.0, 100.0]
   9   
 +-+-+ 
 7   8 
+++ +++
0 2 1 3
```

## 19) `dnakit.phylogeny.draw_tree_cloud_figure()`

- **Purpose:** Overlay several trees that share tip labels as one tree cloud and save it as an image.
- **Method:** Parse every tree with `toytree.tree()`, check that all have the same tip labels, build `toytree.mtree(trees)`, draw it with `MultiTree.draw_cloud_tree()` and save the canvas with `toytree.save()`. Only the first tree carries tip labels and the tip order follows the consensus tree.
- **API:** `draw_tree_cloud_figure(trees: 'Sequence[str]', output_path: 'str | Path', *, image_format: 'str' = 'svg', width: 'int | None' = None, height: 'int | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `trees` is a list of 1 to 100 Newick strings or paths of Newick files with identical, unique tip labels; `output_path` must be a new file whose suffix matches `image_format` (`svg`, `png`, `html` or `pdf`); `width` and `height` are null or integers from 50 to 10,000; `max_tips` and `max_file_bytes` bound each tree.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `image_path`, `image_format`, `image_bytes`, `tree_count`, `tip_count`, `tip_labels` (sorted), `newick` (one string per input tree) and `consensus_newick`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_cloud_figure

trees = [
    "((A:1,B:1):1,(C:2,D:2):1);",
    "((A:1,C:1):1,(B:2,D:2):1);",
    "((A:1,B:1):1,(C:1,D:1):1);",
]

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "cloud.svg"
    result = draw_tree_cloud_figure(trees, out, verbose=False)["result"]
    print(result["image_format"], result["tree_count"], result["tip_count"], result["tip_labels"])
    print(result["newick"][1])
    print(result["consensus_newick"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **Example output:**

```text
svg 3 4 ['A', 'B', 'C', 'D']
((A:1,C:1):1,(B:2,D:2):1);
((A:2,B:2.66666666667)0.666666666667:1,C:2.66666666667,D:3.33333333333);
True
```

## 20) `dnakit.phylogeny.draw_tree_figure()`

- **Purpose:** Draw one tree with toytree to an SVG, PNG, HTML or PDF file and return the data it was drawn from.
- **Method:** Read the Newick text (or file) with `toytree.tree()`, draw it with `ToyTree.draw()` using the chosen layout, style and edge type, add the requested `toytree.annotate` or `toytree.utils` annotations to the axes, and save the canvas with `toytree.save()`. When `tip_labels` is among the annotations, the tip labels come from `add_tip_labels()` (labelled by `feature` or the tip name) instead of `tree.draw`. The node table holds the layout coordinates `x` and `y` of every node.
- **API:** `draw_tree_figure(tree: 'object', output_path: 'str | Path', *, image_format: 'str' = 'svg', layout: 'str' = 'r', tree_style: 'str' = 'n', edge_type: 'str | None' = None, width: 'int | None' = None, height: 'int | None' = None, tip_labels_align: 'bool' = False, annotations: 'Sequence[str]' = (), feature: 'str | None' = None, pie_proportions: 'Sequence[Sequence[float]] | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `tree` is Newick text (it must contain `;`) or a path to a one-tree Newick file; `output_path` is a new file whose suffix matches `image_format` (`svg`, `png`, `html` or `pdf`); `layout` is one of `r`, `l`, `u`, `d`, `c`, `c0-180`, `un`; `tree_style` is one of `n`, `r`, `c`, `s`, `o`, `b`; `edge_type` is null, `p`, `c` or `b`; `annotations` lists names from `tip_labels`, `tip_markers`, `node_markers`, `edge_markers`, `node_labels`, `edge_labels`, `node_bars`, `scale_bar`, `box_outline`, `box_outline_spines`, `external_ticks`, `node_pie_charts`, `edge_pie_charts` without repeats (`node_bars` needs `feature`; the two pie charts need `pie_proportions`, one row per node index whose values sum to 1); `width` and `height` are null or integers from 50 to 10,000.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `image_path`, `image_format`, `image_bytes`, `layout`, `tree_style`, `edge_type`, `annotations`, `tip_count`, `tip_labels` (sorted), `newick`, `nodes` (a list of `idx`, `name`, `is_tip`, `height`, `dist`, `support`, `x`, `y`) and `tip_label_angles`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_figure

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "tree.svg"
    result = draw_tree_figure(
        tree, out, annotations=["tip_labels", "node_markers", "scale_bar"], verbose=False
    )["result"]
    print(result["image_format"], result["layout"], result["annotations"])
    print(result["tip_count"], result["tip_labels"], result["newick"])
    for node in result["nodes"][:3]:
        print(node["idx"], node["name"], node["is_tip"], node["x"], node["y"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **Example output:**

```text
svg r ['tip_labels', 'node_markers', 'scale_bar']
4 ['A', 'B', 'C', 'D'] ((A:1,B:1):1,(C:2,D:2):1);
0 A True -1.0 0.0
1 B True -1.0 1.0
2 C True 0.0 2.0
True
```

## 21) `dnakit.phylogeny.draw_tree_grid_figure()`

- **Purpose:** Draw several trees side by side in a rows by columns grid and save the image.
- **Method:** Parse every tree with `toytree.tree()`, check that all have the same tip labels, build `toytree.mtree(trees)`, draw it with `MultiTree.draw(shape=(rows, columns))` so each tree has its own cell, and save the canvas with `toytree.save()`.
- **API:** `draw_tree_grid_figure(trees: 'Sequence[str]', output_path: 'str | Path', *, image_format: 'str' = 'svg', rows: 'int' = 1, columns: 'int' = 4, width: 'int | None' = None, height: 'int | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `trees` is a list of 1 to 100 Newick strings or file paths with identical, unique tip labels, and at most `rows * columns` of them; `rows` and `columns` are integers from 1 to 20 (defaults 1 and 4); `output_path`, `image_format`, `width`, `height`, `max_tips` and `max_file_bytes` are as in `draw_tree_cloud_figure()`.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `image_path`, `image_format`, `image_bytes`, `shape` (`[rows, columns]`), `tree_count` and `trees` (the Newick string and sorted tip labels of each tree).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_grid_figure

trees = [
    "((A:1,B:1):1,(C:2,D:2):1);",
    "((A:1,C:1):1,(B:2,D:2):1);",
]

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "grid.svg"
    result = draw_tree_grid_figure(trees, out, rows=1, columns=2, verbose=False)["result"]
    print(result["shape"], result["tree_count"])
    for tree in result["trees"]:
        print(tree["newick"], tree["tip_labels"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **Example output:**

```text
[1, 2] 2
((A:1,B:1):1,(C:2,D:2):1); ['A', 'B', 'C', 'D']
((A:1,C:1):1,(B:2,D:2):1); ['A', 'B', 'C', 'D']
True
```

## 22) `dnakit.phylogeny.draw_tree_node_feature_figure()`

- **Purpose:** Draw a tree with its nodes coloured or sized by a numeric node feature and return the mapped values.
- **Method:** Read the tree with `toytree.tree()` and take the feature column from `tree.get_node_data()`. With `encoding="color"` it maps the values with `toytree.style.get_color_mapped_feature()`; with `encoding="size"` it uses `get_range_mapped_feature()` between `min_size` and `max_size` (missing values become 0 unless `nan_value` is set). The tree is drawn with `ToyTree.draw()` using those node colours or sizes and saved with `toytree.save()`.
- **API:** `draw_tree_node_feature_figure(tree: 'object', feature: 'str', output_path: 'str | Path', *, encoding: 'str' = 'color', image_format: 'str' = 'svg', layout: 'str' = 'r', tree_style: 'str' = 'n', tips_only: 'bool' = False, cmap: 'str | None' = None, domain_min: 'float | None' = None, domain_max: 'float | None' = None, nan_value: 'float | None' = None, min_size: 'float' = 5.0, max_size: 'float' = 15.0, reverse: 'bool' = False, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `tree` is Newick text or a one-tree file path; `feature` must be a column of the node data such as `dist`, `height` or `support`, otherwise `ConfigurationError`; `output_path` is a new file whose suffix matches `image_format`; `encoding` is `color` or `size`; `layout` and `tree_style` as in `draw_tree_figure()`; for colours `cmap` (a toyplot colour map name; an unknown name raises `ConfigurationError`), `domain_min`, `domain_max`, `nan_value`, `tips_only` and `reverse`; for sizes `min_size` and `max_size` (0 to 1000, with `max_size >= min_size`).
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `image_path`, `image_format`, `image_bytes`, `feature`, `encoding`, `tips_only`, `tip_labels`, `newick`, `nodes` (as in `draw_tree_figure()`), `feature_values` (one per node index) and `mapped` (hex colours or sizes, one per node index).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_node_feature_figure

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    colored = draw_tree_node_feature_figure(tree, "height", Path(tmp) / "color.svg", verbose=False)["result"]
    sized = draw_tree_node_feature_figure(
        tree, "height", Path(tmp) / "size.svg", encoding="size", verbose=False
    )["result"]
    print(colored["feature_values"])
    print(colored["mapped"])
    print([round(value, 2) for value in sized["mapped"]])
```

- **Example output:**

```text
[1.0, 1.0, 0.0, 0.0, 2.0, 2.0, 3.0]
['#bfe5a0', '#bfe5a0', '#5e4fa2', '#5e4fa2', '#fdbf6f', '#fdbf6f', '#9e0142']
[8.33, 8.33, 5.0, 5.0, 11.67, 11.67, 15.0]
```

## 23) `dnakit.phylogeny.export_tree_node_data()`

- **Purpose:** Return the node data table of a tree and optionally write it to a delimited file.
- **Method:** Read the tree with `toytree.tree()` and take `tree.get_node_data(features)`, or `tree.get_tip_data(features)` when `tips_only` is true. Selected `nodes` are returned in the given order and a name that is not in the tree raises `ConfigurationError`. With `output_path` the table is written with `to_csv(sep=separator, index=False)`, `idx` column first, missing values empty. This is the operation of toytree's `get-node-data` command. Internal node names are available here even though `tree.write()` omits them.
- **API:** `export_tree_node_data(tree: 'object', *, features: 'Sequence[str] | None' = None, nodes: 'Sequence[str] | None' = None, tips_only: 'bool' = False, separator: 'str' = '\t', output_path: 'str | Path | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `tree` is Newick text or a one-tree file path; `features` is null (all columns: `idx`, `name`, `height`, `dist`, `support` and any other feature) or a list of column names; `nodes` is null or a list of distinct node names; `tips_only` is a boolean; `separator` is one character (default tab); `output_path` is null or a new file.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `tips_only`, `columns`, `row_count`, `rows` (dictionaries; missing numbers are null), `output_path`, `output_bytes` and `separator`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import export_tree_node_data

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "nodes.csv"
    result = export_tree_node_data(
        tree, features=["name", "height"], nodes=["n2", "A"], separator=",", output_path=out, verbose=False
    )["result"]
    print(result["columns"], result["row_count"])
    for row in result["rows"]:
        print(row)
    print(out.read_text(encoding="utf-8"), end="")
    tips = export_tree_node_data(tree, tips_only=True, verbose=False)["result"]
    print(tips["row_count"], [row["name"] for row in tips["rows"]])
```

- **Example output:**

```text
['idx', 'name', 'height'] 2
{'idx': 5, 'name': 'n2', 'height': 2.0}
{'idx': 0, 'name': 'A', 'height': 1.0}
idx,name,height
5,n2,2.0
0,A,1.0
4 ['A', 'B', 'C', 'D']
```

## 24) `dnakit.phylogeny.map_node_feature_values()`

- **Purpose:** Map a list of numbers to colours or sizes with toytree's value mappers, without drawing.
- **Method:** Turn `None` entries into NaN, then call `toytree.style.get_color_mapped_values()` (hex colours, toytree's default colour map when `cmap` is null) or `get_range_mapped_values()` (sizes from `min_size` to `max_size`). NaN is filled with `nan_value`; for sizes the default fill is 0.0, and for colours a missing value with no `nan_value` comes out as `#000000`.
- **API:** `map_node_feature_values(values: 'Sequence[float | None]', *, encoding: 'str' = 'color', cmap: 'str | None' = None, domain_min: 'float | None' = None, domain_max: 'float | None' = None, nan_value: 'float | None' = None, min_size: 'float' = 5.0, max_size: 'float' = 15.0, reverse: 'bool' = False, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `values` is a non-empty list of up to 1,000,000 numbers or null; `encoding` is `color` or `size`; `cmap` is null or a toyplot colour map name (unknown names raise `ConfigurationError`); `domain_min`, `domain_max` and `nan_value` are null or finite numbers; `min_size` and `max_size` are numbers from 0 to 1000 with `max_size >= min_size`; `reverse` is a boolean.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `encoding`, `count` and `mapped` (hex colour strings or floats, in input order).
- **Example code:**

```python
from dnakit.phylogeny import map_node_feature_values

values = [0, 1, None, 2]
colors = map_node_feature_values(values, verbose=False)["result"]
filled = map_node_feature_values(values, nan_value=1, domain_min=0, domain_max=2, verbose=False)["result"]
sizes = map_node_feature_values(values, encoding="size", verbose=False)["result"]
print(colors["count"], colors["mapped"])
print(filled["mapped"])
print(sizes["encoding"], sizes["mapped"])
```

- **Example output:**

```text
4 ['#5e4fa2', '#ffffbf', '#000000', '#9e0142']
['#5e4fa2', '#ffffbf', '#ffffbf', '#9e0142']
size [5.0, 10.0, 0.0, 15.0]
```

## 25) `dnakit.phylogeny.prune_tree_tips()`

- **Purpose:** Keep the smallest subtree that connects the named nodes and return it as Newick.
- **Method:** Read the tree with `toytree.tree()`, check that every name is a node name of the tree, and call `toytree.mod.prune(tree, *names, preserve_dists=..., require_root=...)`. Other tips are removed and nodes that become unary are dropped. This is the operation of toytree's `prune` command.
- **API:** `prune_tree_tips(tree: 'object', names: 'Sequence[str]', *, preserve_dists: 'bool' = True, require_root: 'bool' = False, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `tree` is Newick text or a one-tree file path; `names` is a list of at least two distinct names of tips or named internal nodes; `preserve_dists` (default true) keeps the path lengths between retained nodes, and false keeps each retained edge's own length; `require_root` keeps the original root as a node of the result.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `newick`, `tip_count`, `tips`, `retained_names`, `preserve_dists` and `require_root`.
- **Example code:**

```python
from dnakit.phylogeny import prune_tree_tips

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"
merged = prune_tree_tips(tree, ["A", "C"], verbose=False)["result"]
own = prune_tree_tips(tree, ["A", "C"], preserve_dists=False, verbose=False)["result"]
print(merged["newick"], merged["tip_count"], merged["retained_names"])
print(own["newick"], own["preserve_dists"])
```

- **Example output:**

```text
(A:2,C:3); 2 ['A', 'C']
(A:1,C:2); False
```

## 26) `dnakit.phylogeny.root_toytree_on_outgroup()`

- **Purpose:** Root a tree on the edge above an outgroup given by node names and return it as Newick.
- **Method:** Read the tree with `toytree.tree()`, check that every outgroup name is a node name, and call `toytree.mod.root(tree, *outgroup, edge_features=...)`. One name puts the root on that node's edge; several names put it on the edge above their most recent common ancestor. In the source module the function is called `root_tree_on_outgroup`, but it is exported as `root_toytree_on_outgroup` because `dnakit.phylogeny` already has a different `root_tree_on_outgroup` that writes a file.
- **API:** `root_toytree_on_outgroup(tree: 'object', outgroup: 'Sequence[str]', *, edge_features: 'Sequence[str] | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`.
- **Requirements:** needs the optional package toytree (extra `dnakit[bio-toytree]`, which brings toyplot, numpy and pandas); no external program, internet access or model.
- **Input:** `tree` is Newick text or a one-tree file path; `outgroup` is a list of one or more distinct names of tips or named internal nodes; `edge_features` is null or a list of node feature names that toytree splits between the two edges of the new root.
- **Output:** A dictionary with `status`, `operation`, `result` and `provenance`. `result` holds `newick` (the rooted tree), `tip_count`, `tips`, `outgroup` and `edge_features`.
- **Example code:**

```python
from dnakit.phylogeny import root_toytree_on_outgroup

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"
single = root_toytree_on_outgroup(tree, ["D"], verbose=False)["result"]
pair = root_toytree_on_outgroup(tree, ["A", "B"], verbose=False)["result"]
print(single["newick"], single["outgroup"])
print(pair["newick"], pair["tips"])
```

- **Example output:**

```text
(D:1,(C:2,(A:1,B:1):2):1); ['D']
((A:1,B:1):1,(C:2,D:2):1); ['A', 'B', 'C', 'D']
```
