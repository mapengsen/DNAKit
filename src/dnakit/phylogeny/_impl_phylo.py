"""Native Newick, distance-matrix and tree-building helpers behind the phylogeny candidate functions.

Pure standard library and deterministic.  Definitions:

* p-distance = mismatches / compared sites over sites where both sequences have A, C, G or T
  (``U`` counts as ``T``, case-insensitive; gaps and ambiguity codes are skipped pairwise).
* Jukes-Cantor (JC69) distance = -3/4 * ln(1 - 4/3 * p); undefined for p >= 0.75 (saturation).
* UPGMA = average-linkage clustering, branch length = (merge distance / 2) - child height.
* Neighbour-joining = Saitou and Nei (1987) with the Q criterion of Studier and Keppler (1988).
* Robinson-Foulds = size of the symmetric difference of the sets of non-trivial splits.
"""

from __future__ import annotations

import math
import random
import re
from collections import Counter
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_TREE_NODES = 200_000
MAX_TEXT_BYTES = 50_000_000
MAX_TAXA_MATRIX = 400
MAX_TAXA_ALIGNMENT = 400
MAX_ALIGNMENT_COLUMNS = 5_000_000
MAX_BOOTSTRAP_WORK = 30_000_000


class Node:
    """A tree node (root has ``parent`` None)."""

    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name: str = "", length: float | None = None) -> None:
        self.name = name
        self.length = length
        self.children: list[Node] = []
        self.parent: Node | None = None

    def add(self, child: Node) -> Node:
        child.parent = self
        self.children.append(child)
        return child

    @property
    def is_leaf(self) -> bool:
        return not self.children


# --------------------------------------------------------------------------- options


def options_dict(options: Mapping[str, object] | None) -> dict[str, object]:
    return dict(options or {})


def get_int(settings: Mapping[str, object], key: str, default: int, *, low: int, high: int) -> int:
    value = settings.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"options.{key} must be an integer in [{low}, {high}].", code="PHYLOGENY_OPTION", hint=f"Set options['{key}'] to an integer in [{low}, {high}].")
    return value


def get_choice(settings: Mapping[str, object], key: str, default: str, allowed: Sequence[str]) -> str:
    value = settings.get(key, default)
    if not isinstance(value, str) or value.lower() not in allowed:
        raise ConfigurationError(f"options.{key} must be one of {', '.join(allowed)}; got {value!r}.", code="PHYLOGENY_OPTION")
    return value.lower()


def write_output(settings: Mapping[str, object], text: str) -> str | None:
    target = settings.get("output")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or "\x00" in str(target):
        raise ConfigurationError("options.output must be a file path.", code="PHYLOGENY_OPTION")
    path = Path(target)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- text input


def read_text(value: object, what: str = "input") -> str | None:
    """Text of a string/bytes/Path input (existing single-line path strings are read); None for other types."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    path: Path | None = None
    if isinstance(value, Path):
        path = value
    elif isinstance(value, str) and "\n" not in value and 0 < len(value) < 4096 and not value.lstrip().startswith(("(", ">", "#")):
        try:
            if Path(value).is_file():
                path = Path(value)
        except OSError:
            path = None
    if path is not None:
        if not path.is_file():
            raise InputFormatError(f"{what} file does not exist: {path}", code="PHYLOGENY_INPUT_FILE")
        if path.stat().st_size > MAX_TEXT_BYTES:
            raise InputFormatError(f"{what} file is larger than {MAX_TEXT_BYTES} bytes.", code="PHYLOGENY_INPUT_TOO_LARGE")
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, str):
        if len(value) > MAX_TEXT_BYTES:
            raise InputFormatError(f"{what} is larger than {MAX_TEXT_BYTES} characters.", code="PHYLOGENY_INPUT_TOO_LARGE")
        return value
    return None


# --------------------------------------------------------------------------- Newick


_UNQUOTED_STOP = set("(),:;[]'")


def parse_newick(text: str) -> Node:
    """Parse a single Newick tree (comments in [] are dropped; quoted labels supported)."""

    if len(text) > MAX_TEXT_BYTES:
        raise InputFormatError("Newick text is too large.", code="PHYLOGENY_INPUT_TOO_LARGE")
    text = re.sub(r"\[[^\]]*\]", "", text).strip()
    if not text:
        raise InputFormatError("empty Newick string.", code="NEWICK_FORMAT")
    pos = 0
    n = len(text)
    root = Node()
    stack: list[Node] = []
    current = root
    count = 1
    expect_label = False  # True after ')' or at a leaf start

    def skip_space() -> None:
        nonlocal pos
        while pos < n and text[pos].isspace():
            pos += 1

    def read_label() -> str:
        nonlocal pos
        skip_space()
        if pos < n and text[pos] == "'":
            pos += 1
            chars = []
            while pos < n:
                if text[pos] == "'":
                    if pos + 1 < n and text[pos + 1] == "'":
                        chars.append("'")
                        pos += 2
                        continue
                    pos += 1
                    return "".join(chars)
                chars.append(text[pos])
                pos += 1
            raise InputFormatError("unterminated quoted label in Newick.", code="NEWICK_FORMAT")
        start = pos
        while pos < n and text[pos] not in _UNQUOTED_STOP:
            pos += 1
        return text[start:pos].strip()

    def read_length() -> float | None:
        nonlocal pos
        skip_space()
        if pos < n and text[pos] == ":":
            pos += 1
            skip_space()
            start = pos
            while pos < n and text[pos] not in _UNQUOTED_STOP and not text[pos].isspace():
                pos += 1
            try:
                value = float(text[start:pos])
            except ValueError:
                raise InputFormatError(f"invalid branch length {text[start:pos]!r} in Newick.", code="NEWICK_FORMAT") from None
            if not math.isfinite(value):
                raise InputFormatError("non-finite branch length in Newick.", code="NEWICK_FORMAT")
            return value
        return None

    skip_space()
    if pos < n and text[pos] == "(":
        pos += 1
        stack.append(root)
        current = Node()
        root.add(current)
        count += 1
        expect_label = False
    else:
        # a single-leaf tree such as "A;"
        root.name = read_label()
        root.length = read_length()
        skip_space()
        if pos < n and text[pos] == ";":
            pos += 1
        if text[pos:].strip():
            raise InputFormatError("unexpected text after the Newick tree.", code="NEWICK_FORMAT")
        if not root.name:
            raise InputFormatError("empty Newick string.", code="NEWICK_FORMAT")
        return root

    done = False
    while pos < n and not done:
        skip_space()
        if pos >= n:
            break
        ch = text[pos]
        if ch == "(":
            pos += 1
            stack.append(current)
            child = Node()
            current.add(child)
            count += 1
            current = child
        elif ch == ",":
            pos += 1
            if not stack:
                raise InputFormatError("',' outside parentheses in Newick.", code="NEWICK_FORMAT")
            current = stack[-1].add(Node())
            count += 1
        elif ch == ")":
            pos += 1
            if not stack:
                raise InputFormatError("unbalanced ')' in Newick.", code="NEWICK_FORMAT")
            current = stack.pop()
            current.name = read_label()
            current.length = read_length()
            if not stack:
                skip_space()
                if pos < n and text[pos] == ";":
                    pos += 1
                done = True
        elif ch == ";":
            pos += 1
            done = True
        else:
            current.name = read_label()
            current.length = read_length()
            if pos < n and text[pos] not in ",()" and not text[pos].isspace() and text[pos] != ";":
                raise InputFormatError(f"unexpected character {text[pos]!r} in Newick.", code="NEWICK_FORMAT")
        if count > MAX_TREE_NODES:
            raise InputFormatError("Newick tree has too many nodes.", code="PHYLOGENY_INPUT_TOO_LARGE")
    if stack:
        raise InputFormatError("unbalanced '(' in Newick.", code="NEWICK_FORMAT", hint="Check that every '(' has a matching ')'.")
    if text[pos:].strip().strip(";"):
        raise InputFormatError("unexpected text after the Newick tree.", code="NEWICK_FORMAT", hint="Pass one tree per call.")
    names = [leaf.name for leaf in leaves(root)]
    if any(not name for name in names):
        raise InputFormatError("every leaf needs a name.", code="NEWICK_FORMAT")
    duplicates = sorted(name for name, c in Counter(names).items() if c > 1)
    if duplicates:
        raise InputFormatError(f"duplicate leaf names: {', '.join(duplicates[:10])}.", code="NEWICK_FORMAT")
    return root


def load_tree(value: object, what: str = "tree") -> Node:
    text = read_text(value, what)
    if text is None:
        raise InputFormatError(f"{what} must be Newick text or a path to a Newick file.", code="NEWICK_FORMAT", hint="Pass a string such as '((A:1,B:1):1,C:2);'.")
    return parse_newick(text)


def preorder(root: Node) -> Iterator[Node]:
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def postorder(root: Node) -> Iterator[Node]:
    return iter(reversed(list(_reverse_preorder(root))))


def _reverse_preorder(root: Node) -> Iterator[Node]:
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node.children)


def leaves(root: Node) -> list[Node]:
    return [node for node in preorder(root) if node.is_leaf]


def _quote(name: str) -> str:
    if name and not any(c in "(),:;[]' \t\n" for c in name):
        return name
    return "'" + name.replace("'", "''") + "'"


def _fmt(length: float) -> str:
    return format(length, ".12g")


def to_newick(root: Node, *, internal_labels: bool = True) -> str:
    # iterative rendering to avoid recursion limits on caterpillar trees
    order = list(_reverse_preorder(root))
    rendered: dict[int, str] = {}
    for node in reversed(order):
        out = ""
        if node.children:
            out += "(" + ",".join(rendered.pop(id(child)) for child in node.children) + ")"
        label = node.name if node.is_leaf or internal_labels else ""
        out += _quote(label) if label else ""
        if node.length is not None and node.parent is not None:
            out += ":" + _fmt(node.length)
        rendered[id(node)] = out
    return rendered[id(root)] + ";"


# --------------------------------------------------------------------------- tree statistics


def clade_leaf_sets(root: Node) -> dict[int, frozenset[str]]:
    sets: dict[int, frozenset[str]] = {}
    for node in postorder(root):
        if node.is_leaf:
            sets[id(node)] = frozenset([node.name])
        else:
            merged: set[str] = set()
            for child in node.children:
                merged |= sets[id(child)]
            sets[id(node)] = frozenset(merged)
    return sets


def depths(root: Node) -> tuple[dict[int, int], dict[int, float | None]]:
    """Depth in edges and root distance (None if any length on the path is missing)."""

    edge_depth: dict[int, int] = {id(root): 0}
    distance: dict[int, float | None] = {id(root): 0.0}
    for node in preorder(root):
        for child in node.children:
            edge_depth[id(child)] = edge_depth[id(node)] + 1
            base = distance[id(node)]
            distance[id(child)] = None if base is None or child.length is None else base + child.length
    return edge_depth, distance


def tree_statistics(root: Node) -> dict[str, Any]:
    nodes = list(preorder(root))
    leaf_nodes = [n for n in nodes if n.is_leaf]
    internal = [n for n in nodes if not n.is_leaf]
    edge_depth, distance = depths(root)
    non_root = [n for n in nodes if n is not root]
    lengths = [n.length for n in non_root if n.length is not None]
    all_lengths = len(lengths) == len(non_root)
    if any(length < 0 for length in lengths):
        negative = sum(1 for length in lengths if length < 0)
    else:
        negative = 0
    leaf_depths = [edge_depth[id(n)] for n in leaf_nodes]
    stats: dict[str, Any] = {
        "leaf_count": len(leaf_nodes),
        "internal_node_count": len(internal),
        "node_count": len(nodes),
        "leaves": [n.name for n in leaf_nodes],
        "root_degree": len(root.children),
        "rooted_binary": all(len(n.children) == 2 for n in internal),
        "polytomy_count": sum(1 for n in internal if len(n.children) > 2 and not (n is root and len(n.children) == 3)),
        "max_depth_edges": max(leaf_depths),
        "mean_leaf_depth_edges": sum(leaf_depths) / len(leaf_depths),
        "sackin_index": sum(leaf_depths),
        "has_branch_lengths": bool(lengths),
        "all_branches_have_lengths": all_lengths,
        "negative_branch_count": negative,
        "total_branch_length": sum(lengths) if lengths else None,
    }
    if all_lengths and non_root:
        tips = [distance[id(n)] for n in leaf_nodes]
        tip_values = [float(t) for t in tips if t is not None]
        stats["max_root_to_tip"] = max(tip_values)
        stats["min_root_to_tip"] = min(tip_values)
        stats["mean_root_to_tip"] = sum(tip_values) / len(tip_values)
        scale = max(abs(max(tip_values)), 1e-300)
        stats["ultrametric"] = (max(tip_values) - min(tip_values)) <= 1e-9 * scale
        internal_length = sum(n.length or 0.0 for n in non_root if not n.is_leaf)
        total = sum(lengths)
        stats["treeness"] = internal_length / total if total else None
    if stats["rooted_binary"] and len(leaf_nodes) > 1:
        sizes = {k: len(v) for k, v in clade_leaf_sets(root).items()}
        stats["colless_index"] = sum(abs(sizes[id(n.children[0])] - sizes[id(n.children[1])]) for n in internal)
    return stats


def patristic_matrix(root: Node) -> tuple[list[str], list[list[float]]]:
    leaf_nodes = leaves(root)
    if len(leaf_nodes) > 1_000:
        raise InputFormatError("patristic matrices are limited to 1,000 leaves.", code="PHYLOGENY_INPUT_TOO_LARGE")
    nonroot = [n for n in preorder(root) if n is not root]
    if any(n.length is None for n in nonroot):
        raise InputFormatError("every branch needs a length for patristic distances.", code="NEWICK_FORMAT", hint="Provide a Newick tree with branch lengths.")
    _, dist = depths(root)
    # distance via lowest common ancestor on parent chains
    chains: list[list[Node]] = []
    for leaf in leaf_nodes:
        chain = []
        node: Node | None = leaf
        while node is not None:
            chain.append(node)
            node = node.parent
        chains.append(chain)
    names = [leaf.name for leaf in leaf_nodes]
    matrix = [[0.0] * len(names) for _ in names]
    for i in range(len(names)):
        ancestors = {id(n): k for k, n in enumerate(chains[i])}
        for j in range(i + 1, len(names)):
            for node in chains[j]:
                if id(node) in ancestors:
                    d = (dist[id(leaf_nodes[i])] or 0.0) + (dist[id(leaf_nodes[j])] or 0.0) - 2 * (dist[id(node)] or 0.0)
                    matrix[i][j] = matrix[j][i] = d
                    break
    return names, matrix


def annotate_tree(root: Node, metadata: Mapping[str, object] | None, max_clade: int) -> dict[str, Any]:
    nodes = list(preorder(root))
    index = {id(n): i for i, n in enumerate(nodes)}
    edge_depth, distance = depths(root)
    sets = clade_leaf_sets(root)
    if metadata is not None:
        unknown = sorted(set(metadata) - {n.name for n in nodes if n.is_leaf})
        if unknown:
            raise ConfigurationError(f"options.metadata names unknown leaves: {', '.join(map(str, unknown[:10]))}.", code="PHYLOGENY_OPTION")
    labelled = 0
    for node in nodes:
        if not node.is_leaf and not node.name:
            node.name = f"N{index[id(node)]}"
            labelled += 1
    rows = []
    for node in nodes:
        row: dict[str, Any] = {
            "id": index[id(node)],
            "name": node.name,
            "parent": None if node.parent is None else index[id(node.parent)],
            "is_leaf": node.is_leaf,
            "child_count": len(node.children),
            "leaf_count": len(sets[id(node)]),
            "depth_edges": edge_depth[id(node)],
            "branch_length": node.length if node.parent is not None else None,
            "root_distance": distance[id(node)],
        }
        if not node.is_leaf and len(sets[id(node)]) <= max_clade:
            row["clade"] = sorted(sets[id(node)])
        if metadata is not None:
            values = sorted({str(metadata[leaf]) for leaf in sets[id(node)] if leaf in metadata})
            row["metadata_values"] = values
            row["metadata_homogeneous"] = len(values) == 1 if values else None
        rows.append(row)
    return {"nodes": rows, "annotated_newick": to_newick(root), "internal_nodes_labelled": labelled}


# --------------------------------------------------------------------------- splits, Robinson-Foulds


def split_set(root: Node, rooted: bool) -> tuple[set[frozenset[str]], frozenset[str]]:
    sets = clade_leaf_sets(root)
    universe = sets[id(root)]
    result: set[frozenset[str]] = set()
    reference = min(universe)
    for node in preorder(root):
        if node is root or node.is_leaf:
            continue
        clade = sets[id(node)]
        if rooted:
            if 1 < len(clade) < len(universe):
                result.add(clade)
            continue
        side = universe - clade if reference in clade else clade
        if len(side) >= 2 and len(universe) - len(side) >= 2:
            result.add(frozenset(side))
    return result, universe


def robinson_foulds(a: Node, b: Node, rooted: bool) -> dict[str, Any]:
    sa, ua = split_set(a, rooted)
    sb, ub = split_set(b, rooted)
    if ua != ub:
        only_a, only_b = sorted(ua - ub), sorted(ub - ua)
        raise InputFormatError(
            f"the trees have different leaf sets (only in first: {only_a[:5]}, only in second: {only_b[:5]}).",
            code="RF_LEAF_MISMATCH", hint="Prune both trees to the shared leaves first.",
        )
    rf = len(sa ^ sb)
    denominator = len(sa) + len(sb)
    n = len(ua)
    return {
        "robinson_foulds": rf,
        "normalized_rf": rf / denominator if denominator else 0.0,
        "shared_splits": len(sa & sb),
        "splits_only_first": len(sa - sb),
        "splits_only_second": len(sb - sa),
        "splits_first": len(sa),
        "splits_second": len(sb),
        "max_rf_binary": (2 * (n - 2) if rooted else 2 * (n - 3)) if n >= 3 else 0,
        "leaf_count": n,
        "rooted": rooted,
        "definition": "RF = |S1 symmetric-difference S2| over non-trivial splits (bipartitions); normalized = RF / (|S1| + |S2|)",
    }


# --------------------------------------------------------------------------- alignments and distances


def parse_alignment(value: object) -> tuple[list[str], list[str]]:
    text = read_text(value, "alignment")
    pairs: list[tuple[str, str]] = []
    if text is not None:
        if not text.lstrip().startswith(">"):
            raise InputFormatError("alignment text must be FASTA.", code="ALIGNMENT_FORMAT", hint="Pass FASTA text, {name: sequence} or a list of {'name','sequence'}.")
        name: str | None = None
        parts: list[str] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    pairs.append((name, "".join(parts)))
                header = line[1:].split()
                if not header:
                    raise InputFormatError("FASTA header without a name.", code="ALIGNMENT_FORMAT")
                name, parts = header[0], []
            else:
                parts.append(line)
        if name is not None:
            pairs.append((name, "".join(parts)))
    elif isinstance(value, Mapping):
        for key in ("alignment", "sequences", "records"):
            if key in value:
                return parse_alignment(value[key])
        pairs = [(str(k), v) for k, v in value.items()]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for index, item in enumerate(value):
            if isinstance(item, Mapping) and isinstance(item.get("sequence", item.get("seq")), str):
                pairs.append((str(item.get("name", item.get("id", f"seq_{index + 1}"))), item.get("sequence", item.get("seq"))))  # type: ignore[arg-type]
            else:
                raise InputFormatError(f"alignment record {index} must be a mapping with 'name' and 'sequence'.", code="ALIGNMENT_FORMAT")
    else:
        raise InputFormatError("alignment must be FASTA text, a mapping or a list of records.", code="ALIGNMENT_FORMAT")
    if any(not isinstance(seq, str) for _, seq in pairs):
        raise InputFormatError("alignment sequences must be strings.", code="ALIGNMENT_FORMAT")
    if len(pairs) < 2:
        raise InputFormatError("an alignment needs at least two sequences.", code="ALIGNMENT_FORMAT")
    if len(pairs) > MAX_TAXA_ALIGNMENT:
        raise InputFormatError(f"alignments are limited to {MAX_TAXA_ALIGNMENT} sequences.", code="PHYLOGENY_INPUT_TOO_LARGE")
    names = [name for name, _ in pairs]
    if len(set(names)) != len(names):
        raise InputFormatError("sequence names must be unique.", code="ALIGNMENT_FORMAT")
    seqs = [seq.upper().replace("U", "T") for _, seq in pairs]
    if len({len(s) for s in seqs}) != 1 or not seqs[0]:
        raise InputFormatError("all sequences must have the same non-zero length (aligned).", code="ALIGNMENT_FORMAT")
    if len(seqs[0]) * len(seqs) > MAX_ALIGNMENT_COLUMNS:
        raise InputFormatError("alignment is too large.", code="PHYLOGENY_INPUT_TOO_LARGE")
    return names, seqs


_BASES = frozenset("ACGT")


def jc_distance(p: float) -> float | None:
    """Jukes-Cantor distance, None when p >= 0.75."""

    argument = 1.0 - 4.0 * p / 3.0
    if argument <= 0.0:
        return None
    return -0.75 * math.log(argument)


def pairwise_counts(a: str, b: str) -> tuple[int, int]:
    mismatches = compared = 0
    for x, y in zip(a, b):
        if x in _BASES and y in _BASES:
            compared += 1
            mismatches += x != y
    return mismatches, compared


def distance_from_counts(mismatches: int, compared: int, kind: str, saturation: float | None, pair: str) -> float:
    if compared == 0:
        if saturation is None:
            raise InputFormatError(f"sequences {pair} share no comparable (A/C/G/T) sites.", code="DISTANCE_UNDEFINED", hint="Remove one of the sequences or set options.saturation_distance.")
        return saturation
    p = mismatches / compared
    if kind == "p":
        return p
    d = jc_distance(p)
    if d is None:
        if saturation is None:
            raise InputFormatError(f"Jukes-Cantor distance is undefined for sequences {pair} (p = {p:.3f} >= 0.75).", code="DISTANCE_SATURATED", hint="Use options.distance='p' or set options.saturation_distance.")
        return saturation
    return d


def distance_matrix_from_alignment(names: Sequence[str], seqs: Sequence[str], kind: str, saturation: float | None) -> tuple[list[list[float]], list[list[int]]]:
    n = len(names)
    matrix = [[0.0] * n for _ in range(n)]
    compared_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        compared_matrix[i][i] = sum(1 for c in seqs[i] if c in _BASES)
        for j in range(i + 1, n):
            mismatches, compared = pairwise_counts(seqs[i], seqs[j])
            matrix[i][j] = matrix[j][i] = distance_from_counts(mismatches, compared, kind, saturation, f"{names[i]}/{names[j]}")
            compared_matrix[i][j] = compared_matrix[j][i] = compared
    return matrix, compared_matrix


def parse_distance_matrix(value: object) -> tuple[list[str], list[list[float]]]:
    names: list[str]
    matrix: list[list[float]]
    if isinstance(value, Mapping) and "matrix" in value:
        raw = value["matrix"]
        names_raw = value.get("names", value.get("labels"))
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
            raise InputFormatError("'matrix' must be a list of rows.", code="DISTANCE_MATRIX")
        names = [str(x) for x in names_raw] if isinstance(names_raw, Sequence) and not isinstance(names_raw, (str, bytes)) else [f"t{i + 1}" for i in range(len(raw))]
        matrix = [[_number(x) for x in row] for row in raw]  # type: ignore[union-attr]
    elif isinstance(value, Mapping) and value and all(isinstance(v, Mapping) for v in value.values()):
        names = [str(k) for k in value]
        matrix = []
        for a in value:
            row = value[a]
            matrix.append([0.0 if a == b else _number(row[b] if b in row else value[b][a] if a in value[b] else None) for b in value])  # type: ignore[index]
    else:
        raise InputFormatError("distance matrix must be {'names': [...], 'matrix': [[...]]} or a nested {a: {b: d}} mapping.", code="DISTANCE_MATRIX")
    n = len(names)
    if n < 2:
        raise InputFormatError("a distance matrix needs at least two taxa.", code="DISTANCE_MATRIX")
    if n > MAX_TAXA_MATRIX:
        raise InputFormatError(f"distance matrices are limited to {MAX_TAXA_MATRIX} taxa.", code="PHYLOGENY_INPUT_TOO_LARGE")
    if len(set(names)) != n:
        raise InputFormatError("taxon names must be unique.", code="DISTANCE_MATRIX")
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise InputFormatError("the matrix must be square with one row per name.", code="DISTANCE_MATRIX")
    for i in range(n):
        if abs(matrix[i][i]) > 1e-12:
            raise InputFormatError("the diagonal of a distance matrix must be zero.", code="DISTANCE_MATRIX")
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-9 * max(1.0, abs(matrix[i][j])):
                raise InputFormatError(f"the matrix is not symmetric at ({names[i]}, {names[j]}).", code="DISTANCE_MATRIX")
            if matrix[i][j] < 0:
                raise InputFormatError("distances must be non-negative.", code="DISTANCE_MATRIX")
    return names, matrix


def _number(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise InputFormatError("distance matrix entries must be finite numbers.", code="DISTANCE_MATRIX")
    return float(value)


# --------------------------------------------------------------------------- tree building


def neighbor_joining(names: Sequence[str], matrix: Sequence[Sequence[float]]) -> Node:
    n = len(names)
    if n == 2:
        root = Node()
        half = matrix[0][1] / 2.0
        root.add(Node(names[0], half))
        root.add(Node(names[1], half))
        return root
    nodes: list[Node] = [Node(name) for name in names]
    dist = {(i, j): float(matrix[i][j]) for i in range(n) for j in range(n) if i != j}
    active = list(range(n))
    next_id = n

    def d(i: int, j: int) -> float:
        return dist[(i, j)]

    while len(active) > 3:
        m = len(active)
        totals = {i: sum(d(i, j) for j in active if j != i) for i in active}
        best: tuple[float, int, int] | None = None
        for a in range(m):
            for b in range(a + 1, m):
                i, j = active[a], active[b]
                q = (m - 2) * d(i, j) - totals[i] - totals[j]
                if best is None or q < best[0] - 1e-12:
                    best = (q, i, j)
        assert best is not None
        _, i, j = best
        dij = d(i, j)
        li = dij / 2.0 + (totals[i] - totals[j]) / (2.0 * (m - 2))
        lj = dij - li
        parent = Node()
        nodes.append(parent)
        nodes[i].length, nodes[j].length = li, lj
        parent.add(nodes[i])
        parent.add(nodes[j])
        u = next_id
        next_id += 1
        for k in active:
            if k not in (i, j):
                value = (d(i, k) + d(j, k) - dij) / 2.0
                dist[(u, k)] = dist[(k, u)] = value
        active = [k for k in active if k not in (i, j)] + [u]
    i, j, k = active
    root = Node()
    lengths = (
        (d(i, j) + d(i, k) - d(j, k)) / 2.0,
        (d(i, j) + d(j, k) - d(i, k)) / 2.0,
        (d(i, k) + d(j, k) - d(i, j)) / 2.0,
    )
    for index, length in zip((i, j, k), lengths):
        nodes[index].length = length
        root.add(nodes[index])
    return root


def upgma(names: Sequence[str], matrix: Sequence[Sequence[float]]) -> Node:
    n = len(names)
    nodes: dict[int, Node] = {i: Node(names[i]) for i in range(n)}
    height = {i: 0.0 for i in range(n)}
    size = {i: 1 for i in range(n)}
    dist = {(i, j): float(matrix[i][j]) for i in range(n) for j in range(n) if i != j}
    active = list(range(n))
    next_id = n
    while len(active) > 1:
        best: tuple[float, int, int] | None = None
        for a in range(len(active)):
            for b in range(a + 1, len(active)):
                i, j = active[a], active[b]
                if best is None or dist[(i, j)] < best[0] - 1e-12:
                    best = (dist[(i, j)], i, j)
        assert best is not None
        dij, i, j = best
        new_height = dij / 2.0
        parent = Node()
        for child in (i, j):
            nodes[child].length = new_height - height[child]
            parent.add(nodes[child])
        u = next_id
        next_id += 1
        nodes[u], height[u], size[u] = parent, new_height, size[i] + size[j]
        for k in active:
            if k not in (i, j):
                value = (dist[(i, k)] * size[i] + dist[(j, k)] * size[j]) / (size[i] + size[j])
                dist[(u, k)] = dist[(k, u)] = value
        active = [k for k in active if k not in (i, j)] + [u]
    return nodes[active[0]]


def build_tree(names: Sequence[str], matrix: Sequence[Sequence[float]], method: str) -> Node:
    return neighbor_joining(names, matrix) if method == "nj" else upgma(names, matrix)


# --------------------------------------------------------------------------- bootstrap


def bootstrap_support(names: Sequence[str], seqs: Sequence[str], method: str, kind: str, replicates: int, seed: int, saturation: float) -> dict[str, Any]:
    n, length = len(names), len(seqs[0])
    matrix, _ = distance_matrix_from_alignment(names, seqs, kind, saturation)
    reference = build_tree(names, matrix, method)
    ref_splits, _ = split_set(reference, rooted=(method == "upgma"))
    # group identical columns into patterns
    columns: Counter[str] = Counter("".join(seq[c] for seq in seqs) for c in range(length))
    patterns = sorted(columns)
    pattern_weights = [columns[p] for p in patterns]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    work = len(pairs) * len(patterns) * replicates
    if work > MAX_BOOTSTRAP_WORK:
        raise InputFormatError(f"bootstrap work ({work}) exceeds the limit {MAX_BOOTSTRAP_WORK}; reduce replicates or taxa.", code="PHYLOGENY_INPUT_TOO_LARGE")
    mis: dict[tuple[int, int], list[int]] = {}
    cmp_: dict[tuple[int, int], list[int]] = {}
    for i, j in pairs:
        mis[(i, j)] = [1 if (p[i] in _BASES and p[j] in _BASES and p[i] != p[j]) else 0 for p in patterns]
        cmp_[(i, j)] = [1 if (p[i] in _BASES and p[j] in _BASES) else 0 for p in patterns]
    rng = random.Random(seed)
    hits: Counter[frozenset[str]] = Counter()
    cumulative = []
    running = 0
    for w in pattern_weights:
        running += w
        cumulative.append(running)
    import bisect

    for _ in range(replicates):
        draws = Counter(bisect.bisect_right(cumulative, rng.randrange(length)) for _ in range(length))
        weights = [draws.get(k, 0) for k in range(len(patterns))]
        rep = [[0.0] * n for _ in range(n)]
        for i, j in pairs:
            m = sum(w * x for w, x in zip(weights, mis[(i, j)]))
            c = sum(w * x for w, x in zip(weights, cmp_[(i, j)]))
            rep[i][j] = rep[j][i] = distance_from_counts(m, c, kind, saturation, f"{names[i]}/{names[j]}")
        tree = build_tree(names, rep, method)
        splits, _ = split_set(tree, rooted=(method == "upgma"))
        for split in ref_splits & splits:
            hits[split] += 1
    clades = [
        {"clade": sorted(split), "support": hits[split] / replicates}
        for split in sorted(ref_splits, key=lambda s: (len(s), sorted(s)))
    ]
    return {"reference": reference, "clades": clades, "split_hits": hits, "ref_splits": ref_splits}


def label_support(root: Node, support: Mapping[frozenset[str], float], rooted: bool) -> None:
    sets = clade_leaf_sets(root)
    universe = sets[id(root)]
    reference = min(universe)
    for node in preorder(root):
        if node is root or node.is_leaf:
            continue
        clade = sets[id(node)]
        key = clade if rooted else frozenset(universe - clade if reference in clade else clade)
        if key in support:
            node.name = str(round(100 * support[key]))


# --------------------------------------------------------------------------- export and plot


def tree_to_json(root: Node) -> dict[str, Any]:
    out: dict[int, dict[str, Any]] = {}
    for node in postorder(root):
        entry: dict[str, Any] = {"name": node.name or None, "branch_length": node.length if node.parent is not None else None}
        if node.children:
            entry["children"] = [out[id(c)] for c in node.children]
        out[id(node)] = entry
    return out[id(root)]


def edge_list(root: Node) -> list[dict[str, Any]]:
    nodes = list(preorder(root))
    index = {id(n): i for i, n in enumerate(nodes)}
    return [
        {"parent": index[id(n.parent)], "child": index[id(n)], "child_name": n.name or None, "length": n.length}
        for n in nodes if n.parent is not None
    ]


def to_nexus(root: Node) -> str:
    names = [leaf.name for leaf in leaves(root)]
    rooted = "[&R] "
    return (
        "#NEXUS\n\nBEGIN TAXA;\n"
        f"  DIMENSIONS NTAX={len(names)};\n  TAXLABELS {' '.join(_quote(n) for n in names)};\nEND;\n\n"
        f"BEGIN TREES;\n  TREE tree1 = {rooted}{to_newick(root)}\nEND;\n"
    )


def tree_svg(root: Node, use_lengths: bool, width: int = 640) -> tuple[str, int]:
    leaf_nodes = leaves(root)
    if len(leaf_nodes) > 2_000:
        raise InputFormatError("plot is limited to 2,000 leaves.", code="PHYLOGENY_INPUT_TOO_LARGE")
    edge_depth, distance = depths(root)
    if use_lengths:
        x_of = {k: (v if v is not None else 0.0) for k, v in distance.items()}
        if any(v is None for v in distance.values()):
            raise InputFormatError("branch lengths are missing; use options.branch_lengths=False for a cladogram.", code="NEWICK_FORMAT")
    else:
        x_of = {k: float(v) for k, v in edge_depth.items()}
    max_x = max(x_of.values()) or 1.0
    label_space = 8 * max(len(leaf.name) for leaf in leaf_nodes) + 12
    left, top, row = 10.0, 10.0, 18.0
    plot_width = max(width - label_space - left, 50.0)
    height = int(top * 2 + row * len(leaf_nodes))
    y: dict[int, float] = {}
    for i, leaf in enumerate(leaf_nodes):
        y[id(leaf)] = top + row * (i + 0.5)
    for node in postorder(root):
        if node.children:
            ys = [y[id(c)] for c in node.children]
            y[id(node)] = (min(ys) + max(ys)) / 2.0
    px = lambda node: left + plot_width * x_of[id(node)] / max_x  # noqa: E731
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    parts.append('<g stroke="#333" stroke-width="1" fill="none">')
    for node in preorder(root):
        if node.children:
            ys = [y[id(c)] for c in node.children]
            parts.append(f'<line x1="{px(node):.1f}" y1="{min(ys):.1f}" x2="{px(node):.1f}" y2="{max(ys):.1f}"/>')
        if node.parent is not None:
            parts.append(f'<line x1="{px(node.parent):.1f}" y1="{y[id(node)]:.1f}" x2="{px(node):.1f}" y2="{y[id(node)]:.1f}"/>')
    parts.append("</g>")
    for leaf in leaf_nodes:
        parts.append(f'<text x="{px(leaf) + 4:.1f}" y="{y[id(leaf)] + 4:.1f}" font-size="12" font-family="sans-serif">{escape(leaf.name)}</text>')
    parts.append("</svg>")
    return "\n".join(parts), len(leaf_nodes)
