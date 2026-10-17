"""Newick tree engine (also used for the TREES block of NEXUS files).

Trees are parsed into :class:`Node` objects by an iterative scanner (no recursion, so very deep
"caterpillar" trees are fine) and serialised back by an iterative writer.  Comments in square
brackets are discarded; underscores in unquoted labels are kept as they are.
"""

# ruff: noqa: E501, RUF012, B905

from __future__ import annotations

import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from dnakit.exceptions import ConfigurationError
from dnakit.io_formats._impl_common import (
    Format,
    Opts,
    Payload,
    bad_input,
    decode_text,
    fmt_number,
    json_payload,
    to_float,
)

MAX_NODES = 2_000_000


class Node:
    __slots__ = ("children", "length", "name")

    def __init__(
        self,
        name: str | None = None,
        length: float | None = None,
        children: list[Node] | None = None,
    ) -> None:
        self.name = name
        self.length = length
        self.children: list[Node] = children if children is not None else []


@dataclass
class TreeEntry:
    root: Node
    name: str | None = None
    rooted: bool | None = None
    offset: int = 0
    length: int = 0


class NewickSyntaxError(ValueError):
    pass


_STOP = set("(),:;[]")


def _skip_comment(text: str, i: int) -> int:
    depth = 0
    n = len(text)
    while i < n:
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise NewickSyntaxError("unterminated [comment]")


def _read_label(text: str, i: int) -> tuple[str, int]:
    n = len(text)
    if text[i] == "'":
        out = []
        i += 1
        while i < n:
            if text[i] == "'":
                if i + 1 < n and text[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                return "".join(out), i + 1
            out.append(text[i])
            i += 1
        raise NewickSyntaxError("unterminated quoted label")
    j = i
    while j < n and not text[j].isspace() and text[j] not in _STOP:
        j += 1
    return text[i:j], j


def parse_one(text: str, i: int, budget: int = MAX_NODES) -> tuple[Node, int]:
    """Parse one tree starting at ``text[i]``; return the root and the index after its ';'."""

    n = len(text)
    root = cur = Node()
    stack: list[Node] = []
    count = 1
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c == "[":
            i = _skip_comment(text, i)
        elif c == "(":
            if cur.children or cur.name is not None or cur.length is not None:
                raise NewickSyntaxError(f"unexpected '(' at position {i}")
            child = Node()
            cur.children.append(child)
            stack.append(cur)
            cur = child
            count += 1
            i += 1
        elif c == ",":
            if not stack:
                raise NewickSyntaxError(f"',' outside parentheses at position {i}")
            cur = Node()
            stack[-1].children.append(cur)
            count += 1
            i += 1
        elif c == ")":
            if not stack:
                raise NewickSyntaxError(f"unbalanced ')' at position {i}")
            cur = stack.pop()
            i += 1
        elif c == ":":
            if cur.length is not None:
                raise NewickSyntaxError(f"two branch lengths on one node at position {i}")
            i += 1
            while i < n and (text[i].isspace() or text[i] == "["):
                i = _skip_comment(text, i) if text[i] == "[" else i + 1
            j = i
            while j < n and not text[j].isspace() and text[j] not in _STOP:
                j += 1
            try:
                cur.length = to_float(text[i:j])
            except ValueError as exc:
                raise NewickSyntaxError(f"bad branch length at position {i}: {exc}") from exc
            i = j
        elif c == ";":
            if stack:
                raise NewickSyntaxError("unbalanced parentheses: missing ')'")
            return root, i + 1
        else:
            if cur.name is not None or cur.length is not None:
                raise NewickSyntaxError(f"unexpected label at position {i}")
            label, i = _read_label(text, i)
            cur.name = label or None
        if count > budget:
            raise ConfigurationError(
                f"Tree has more than {budget} nodes.", code="IO_FORMATS_TOO_LARGE"
            )
    raise NewickSyntaxError("tree is not terminated by ';'")


def _quote_label(name: str | None) -> str:
    if not name:
        return ""
    if any(ch.isspace() or ch in "()[]':;," for ch in name):
        return "'" + name.replace("'", "''") + "'"
    return name


def to_newick(root: Node, quote: Any = _quote_label) -> str:
    out: list[str] = []
    stack: list[tuple[str, Any]] = [("n", root)]
    while stack:
        kind, value = stack.pop()
        if kind == "t":
            out.append(value)
            continue
        suffix = quote(value.name) + (
            ":" + fmt_number(value.length) if value.length is not None else ""
        )
        kids = value.children
        if kids:
            out.append("(")
            stack.append(("t", ")" + suffix))
            for idx in range(len(kids) - 1, -1, -1):
                stack.append(("n", kids[idx]))
                if idx:
                    stack.append(("t", ","))
        else:
            out.append(suffix)
    return "".join(out) + ";"


def preorder(root: Node) -> Iterator[Node]:
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def postorder(root: Node) -> list[Node]:
    out: list[Node] = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(node.children)
    out.reverse()
    return out


def tip_names(root: Node) -> list[str | None]:
    return [n.name for n in preorder(root) if not n.children]


def clone(root: Node) -> Node:
    copies: dict[int, Node] = {}
    for node in postorder(root):
        copies[id(node)] = Node(node.name, node.length, [copies[id(c)] for c in node.children])
    return copies[id(root)]


def _join_length(a: float | None, b: float | None) -> float | None:
    return None if a is None and b is None else (a or 0.0) + (b or 0.0)


def prune(root: Node, keep: set[str]) -> Node | None:
    """Drop tips not in ``keep``; remove emptied clades and merge the branches of unary nodes."""

    removed: set[int] = set()
    replaced: dict[int, Node] = {}
    for node in postorder(root):
        if not node.children:
            if node.name not in keep:
                removed.add(id(node))
            continue
        kids = [replaced.get(id(c), c) for c in node.children if id(c) not in removed]
        node.children = kids
        if not kids:
            removed.add(id(node))
        elif len(kids) == 1 and node is not root:
            kids[0].length = _join_length(node.length, kids[0].length)
            replaced[id(node)] = kids[0]
    if id(root) in removed:
        return None
    while len(root.children) == 1:
        child = root.children[0]
        child.length = _join_length(root.length, child.length)
        root = child
    return root


def collapse_short(root: Node, threshold: float) -> Node:
    """Contract internal branches shorter than ``threshold``; descendants keep their own lengths."""

    for node in postorder(root):
        if not node.children:
            continue
        kids: list[Node] = []
        for child in node.children:
            if child.children and child.length is not None and child.length < threshold:
                kids.extend(child.children)
            else:
                kids.append(child)
        node.children = kids
    return root


def ladderize(root: Node, order: str) -> Node:
    tips: dict[int, int] = {}
    first: dict[int, str] = {}
    for node in postorder(root):
        if not node.children:
            tips[id(node)] = 1
            first[id(node)] = node.name or ""
        else:
            tips[id(node)] = sum(tips[id(c)] for c in node.children)
            first[id(node)] = min(first[id(c)] for c in node.children)
            if order == "alphabetical":
                node.children.sort(key=lambda c: first[id(c)])
            elif order == "ascending":
                node.children.sort(key=lambda c: (tips[id(c)], first[id(c)]))
            else:
                node.children.sort(key=lambda c: (-tips[id(c)], first[id(c)]))
    return root


def flat_nodes(root: Node) -> list[dict[str, object]]:
    ids: dict[int, int] = {}
    out: list[dict[str, object]] = []
    parents: dict[int, int | None] = {id(root): None}
    for node in preorder(root):
        ids[id(node)] = len(out)
        out.append(
            {"id": len(out), "parent": parents[id(node)], "name": node.name, "length": node.length}
        )
        for child in node.children:
            parents[id(child)] = ids[id(node)]
    return out


def from_flat(nodes: Sequence[Mapping[str, Any]]) -> Node:
    built: dict[Any, Node] = {}
    root: Node | None = None
    for item in nodes:
        length = item.get("length")
        node = Node(item.get("name") or None, None if length is None else float(length))
        built[item["id"]] = node
        parent = item.get("parent")
        if parent is None:
            if root is not None:
                raise ValueError("more than one root (parent=null)")
            root = node
        else:
            if parent not in built:
                raise ValueError(f"node {item['id']} appears before its parent {parent}")
            built[parent].children.append(node)
    if root is None:
        raise ValueError("no root node")
    return root


def from_nested(spec: Mapping[str, Any]) -> Node:
    root = Node()
    work: list[tuple[Mapping[str, Any], Node]] = [(spec, root)]
    while work:
        item, node = work.pop()
        node.name = str(item["name"]) if item.get("name") not in (None, "") else None
        length = item.get("length")
        node.length = None if length is None else float(length)
        for child in item.get("children", []) or []:
            sub = Node()
            node.children.append(sub)
            work.append((child, sub))
    return root


def tree_stats(root: Node) -> dict[str, object]:
    depth = {id(root): 0}
    dist = {id(root): 0.0}
    tips = internal = 0
    total = 0.0
    max_depth = 0
    height = 0.0
    names: list[str | None] = []
    negative = 0
    has_lengths = False
    binary = True
    for node in preorder(root):
        for child in node.children:
            depth[id(child)] = depth[id(node)] + 1
            dist[id(child)] = dist[id(node)] + (child.length or 0.0)
        if node.length is not None:
            has_lengths = True
            total += node.length
            negative += node.length < 0
        if node.children:
            internal += 1
            allowed = (2, 3) if node is root and len(node.children) == 3 else (2,)
            binary = binary and len(node.children) in allowed
        else:
            tips += 1
            names.append(node.name)
            max_depth = max(max_depth, depth[id(node)])
            height = max(height, dist[id(node)])
    named = [x for x in names if x]
    return {
        "tip_count": tips,
        "internal_node_count": internal,
        "node_count": tips + internal,
        "root_degree": len(root.children),
        "rooted": len(root.children) == 2,
        "fully_resolved": binary,
        "max_depth_edges": max_depth,
        "height": height,
        "total_branch_length": total if has_lengths else None,
        "has_branch_lengths": has_lengths,
        "negative_branch_count": negative,
        "unnamed_tip_count": tips - len(named),
        "duplicate_tip_count": len(named) - len(set(named)),
    }


# NEXUS helpers ------------------------------------------------------------------------------

_NEXUS_SAFE = re.compile(r"[A-Za-z0-9_.|]+\Z")


def nexus_token(name: str) -> str:
    if _NEXUS_SAFE.match(name):
        return name
    return "'" + name.replace("'", "''") + "'"


def nexus_taxa_block(labels: Sequence[str]) -> str:
    return (
        "BEGIN TAXA;\n"
        + f"    DIMENSIONS NTAX={len(labels)};\n"
        + "    TAXLABELS "
        + " ".join(nexus_token(t) for t in labels)
        + ";\nEND;\n"
    )


def nexus_trees_block(entries: Sequence[TreeEntry]) -> str:
    lines = ["BEGIN TREES;"]
    for number, entry in enumerate(entries, 1):
        flag = {True: "[&R] ", False: "[&U] ", None: ""}[entry.rooted]
        lines.append(
            f"    TREE {nexus_token(entry.name or f'tree{number}')} = {flag}{to_newick(entry.root, lambda n: nexus_token(n) if n else '')}"
        )
    lines.append("END;")
    return "\n".join(lines) + "\n"


# =============================================================================================
# the Newick format engine
# =============================================================================================


class NewickFormat(Format):
    name = "newick"
    title = "Newick tree file"
    noun = "trees"
    suffix = ".nwk"
    targets = ("newick", "nexus", "tsv", "json")
    options = {
        "read": {"max_nodes": "int, default 2000000: maximum nodes per tree"},
        "filter": {
            "keep_tips": "list of tip names to keep (other tips are pruned; unary nodes are merged, summing branch lengths)",
            "drop_tips": "list of tip names to prune",
            "collapse_below": "float: contract internal branches shorter than this into polytomies",
            "min_tips": "int: drop trees with fewer tips (evaluated after pruning)",
            "max_tips": "int: drop trees with more tips (evaluated after pruning)",
        },
        "sort": {
            "order": "'ascending' (default: clades with fewer tips first), 'descending' or 'alphabetical' "
            "(by the smallest tip name of each clade); ties are broken by that name"
        },
        "merge": {
            "deduplicate": "bool, default False: drop trees whose normalised Newick equals an earlier tree"
        },
        "split": {"trees_per_part": "int, default 1: trees in each output part"},
        "validate": {},
        "summarize": {},
        "index": {},
        "convert": {},
        "stream": {},
        "write": {},
    }
    semantics = {
        "validate": "Each ';'-terminated tree must have balanced parentheses, well-formed labels and finite numeric branch "
        "lengths. Duplicate or missing tip names are errors; negative branch lengths and single-child nodes are warnings.",
        "summarize": "Per tree: tip/internal/node counts, root degree, rooted (= root has exactly 2 children), fully-resolved flag, "
        "depth in edges, height (longest root-to-tip path length; missing lengths count as 0), total branch length.",
        "index": "Gives each tree its byte offset and length in the text and lists, for every tip name, the trees containing it.",
        "convert": "Targets: 'newick' (normalised), 'nexus' (TREES block, tree names tree1..), 'tsv' (node table: tree, node_id, "
        "parent_id, name, length, is_tip) and 'json' (records with flat node tables).",
        "read": "Records carry the normalised Newick string, tip names and a flat node table (id, parent, name, length in preorder).",
    }

    # -- parse ----------------------------------------------------------------------------
    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> list[TreeEntry]:
        text = decode_text(data)
        budget = opts.integer("max_nodes", MAX_NODES, 1) or MAX_NODES
        entries: list[TreeEntry] = []
        i, n = 0, len(text)
        byte_pos = char_pos = 0
        while True:
            while i < n and text[i].isspace():
                i += 1
            while i < n and text[i] == "[":  # leading comment between trees
                try:
                    i = _skip_comment(text, i)
                except NewickSyntaxError as exc:
                    self.bad(errors, "", str(exc))
                    return entries
                while i < n and text[i].isspace():
                    i += 1
            if i >= n:
                break
            try:
                root, end = parse_one(text, i, budget)
            except NewickSyntaxError as exc:
                self.bad(errors, f"tree {len(entries) + 1}", str(exc))
                break
            byte_pos += len(text[char_pos:i].encode("utf-8"))
            span = len(text[i:end].encode("utf-8"))
            entries.append(TreeEntry(root, None, None, byte_pos, span))
            byte_pos += span
            char_pos = end
            i = end
        if not entries and errors is None:
            raise bad_input(
                "No tree found in the Newick input.",
                code="IO_FORMATS_BAD_NEWICK",
                hint="A Newick tree looks like '(A:1,B:1);'.",
            )
        if not entries and errors is not None and not errors:
            errors.append("no tree found in the input")
        return entries

    def from_records(self, records: object, opts: Opts) -> list[TreeEntry]:
        if isinstance(records, Mapping) and not any(
            k in records for k in ("newick", "nodes", "children", "name")
        ):
            for key in ("trees", "records", "items"):
                if key in records:
                    records = records[key]
                    break
        items = [records] if isinstance(records, Mapping) else records
        if isinstance(items, str) or not isinstance(items, Sequence):
            raise bad_input("write_newick needs Newick text, a tree dict or a list of them.")
        entries: list[TreeEntry] = []
        for number, item in enumerate(items, 1):
            try:
                if isinstance(item, str):
                    entries.extend(self.parse(item.encode("utf-8"), None, opts))
                elif isinstance(item, Mapping) and "newick" in item:
                    entries.extend(self.parse(str(item["newick"]).encode("utf-8"), None, opts))
                elif isinstance(item, Mapping) and "nodes" in item:
                    entries.append(TreeEntry(from_flat(item["nodes"])))
                elif isinstance(item, Mapping):
                    entries.append(TreeEntry(from_nested(item)))
                else:
                    raise TypeError(f"unsupported tree description of type {type(item).__name__}")
            except (KeyError, TypeError, ValueError) as exc:
                if isinstance(exc, ConfigurationError):
                    raise
                raise bad_input(
                    f"tree {number} is malformed: {exc}",
                    hint="Use Newick strings, {'newick': ...}, {'nodes': [...]} or nested {'name','length','children'} dicts.",
                ) from exc
        return entries

    # -- records / dump -------------------------------------------------------------------
    def records(self, doc: list[TreeEntry]) -> list[dict[str, object]]:
        out = []
        for i, entry in enumerate(doc):
            out.append(
                {
                    "index": i,
                    "newick": to_newick(entry.root),
                    "tips": tip_names(entry.root),
                    "nodes": flat_nodes(entry.root),
                }
            )
        return out

    def dump(self, doc: list[TreeEntry], opts: Opts) -> Payload:
        return Payload(
            "newick", text="".join(to_newick(e.root) + "\n" for e in doc), count=len(doc)
        )

    # -- operations -----------------------------------------------------------------------
    def validate_doc(
        self, doc: list[TreeEntry], opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        tips_total = 0
        for i, entry in enumerate(doc, 1):
            stats = tree_stats(entry.root)
            tips_total += stats["tip_count"]  # type: ignore[operator]
            if stats["unnamed_tip_count"]:
                errors.append(f"tree {i}: {stats['unnamed_tip_count']} tip(s) without a name")
            if stats["duplicate_tip_count"]:
                errors.append(f"tree {i}: {stats['duplicate_tip_count']} duplicated tip name(s)")
            if stats["negative_branch_count"]:
                warnings.append(
                    f"tree {i}: {stats['negative_branch_count']} negative branch length(s)"
                )
            if any(len(n.children) == 1 for n in preorder(entry.root)):
                warnings.append(f"tree {i}: contains single-child (unary) internal nodes")
        return {"tree_count": len(doc), "tip_count": tips_total}

    def filter(self, doc: list[TreeEntry], opts: Opts) -> list[TreeEntry]:
        keep, drop = opts.names("keep_tips"), opts.names("drop_tips")
        threshold = opts.number("collapse_below", None, 0.0)
        min_tips, max_tips = opts.integer("min_tips", None, 0), opts.integer("max_tips", None, 0)
        out = []
        for entry in doc:
            root: Node | None = clone(entry.root)
            if keep is not None or drop is not None:
                names = {n for n in tip_names(root) if n}
                wanted = (set(keep) if keep is not None else set(names)) - set(drop or [])
                root = prune(root, wanted)
                if root is None:
                    continue
            if threshold is not None:
                root = collapse_short(root, threshold)
            count = len(tip_names(root))
            if (min_tips is not None and count < min_tips) or (
                max_tips is not None and count > max_tips
            ):
                continue
            out.append(TreeEntry(root, entry.name, entry.rooted))
        return out

    def sort(self, doc: list[TreeEntry], opts: Opts) -> list[TreeEntry]:
        order = str(opts.text("order", "ascending", ["ascending", "descending", "alphabetical"]))
        return [TreeEntry(ladderize(clone(e.root), order), e.name, e.rooted) for e in doc]

    def merge(self, docs: list[list[TreeEntry]], opts: Opts) -> list[TreeEntry]:
        dedupe = opts.boolean("deduplicate", False)
        out, seen = [], set()
        for doc in docs:
            for entry in doc:
                if dedupe:
                    text = to_newick(entry.root)
                    if text in seen:
                        continue
                    seen.add(text)
                out.append(entry)
        return out

    def split(self, doc: list[TreeEntry], opts: Opts) -> list[tuple[str, list[TreeEntry]]]:
        size = opts.integer("trees_per_part", 1, 1) or 1
        return [(f"{i // size + 1:04d}", doc[i : i + size]) for i in range(0, len(doc), size)]

    def summarize(self, doc: list[TreeEntry], opts: Opts) -> dict[str, object]:
        per_tree = [tree_stats(e.root) for e in doc]
        return {
            "tree_count": len(doc),
            "trees": per_tree,
            "tip_count_min": min((s["tip_count"] for s in per_tree), default=None),  # type: ignore[type-var]
            "tip_count_max": max((s["tip_count"] for s in per_tree), default=None),
        }  # type: ignore[type-var]

    def index(self, doc: list[TreeEntry], opts: Opts) -> dict[str, object]:
        trees = []
        tips: dict[str, list[int]] = {}
        for i, entry in enumerate(doc):
            names = [t for t in tip_names(entry.root) if t]
            trees.append(
                {
                    "index": i,
                    "offset": entry.offset,
                    "length": entry.length,
                    "tip_count": len(names),
                }
            )
            for name in dict.fromkeys(names):
                tips.setdefault(name, []).append(i)
        return {"key": "tree index", "record_count": len(doc), "trees": trees, "tips": tips}

    def convert(self, doc: list[TreeEntry], target: str, opts: Opts) -> Payload:
        if target == "newick":
            return self.dump(doc, opts)
        if target == "nexus":
            labels = list(dict.fromkeys(t for e in doc for t in tip_names(e.root) if t))
            return Payload(
                "nexus",
                text="#NEXUS\n\n" + nexus_taxa_block(labels) + "\n" + nexus_trees_block(doc),
                count=len(doc),
            )
        if target == "tsv":
            lines = ["tree\tnode_id\tparent_id\tname\tlength\tis_tip"]
            for i, entry in enumerate(doc):
                kids = {id(n) for n in preorder(entry.root) if not n.children}
                for node, row in zip(preorder(entry.root), flat_nodes(entry.root)):
                    lines.append(
                        "\t".join(
                            [
                                str(i),
                                str(row["id"]),
                                "" if row["parent"] is None else str(row["parent"]),
                                node.name or "",
                                "" if node.length is None else fmt_number(node.length),
                                "1" if id(node) in kids else "0",
                            ]
                        )
                    )
            return Payload("tsv", text="\n".join(lines) + "\n", count=len(doc))
        return json_payload("newick", self.records(doc), len(doc))
