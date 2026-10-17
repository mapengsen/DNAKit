"""Pure-Python property-graph operations behind the ``knowledge_graph`` candidate functions.

A graph is ``{"directed": bool, "nodes": [{"id": str, ...properties}], "edges": [{"source", "relation",
"target", ...properties}]}``. Edge identity is ``(source, relation, target)``; for undirected graphs the two
endpoints are unordered. Everything is deterministic and bounded; nothing touches the network or a database.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from collections import deque
from collections.abc import Mapping, Sequence
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_FILE_BYTES = 100_000_000
DEFAULT_MAX_EDGES = 500_000
DEFAULT_RELATION = "related_to"
_SOURCE_KEYS = ("subject", "source", "head", "from", "s", "src", "id1")
_RELATION_KEYS = ("predicate", "relation", "relationship", "type", "p", "rel", "label")
_TARGET_KEYS = ("object", "target", "tail", "to", "o", "dst", "id2")
_HEADER_WORDS = {
    "subject",
    "source",
    "head",
    "from",
    "src",
    "predicate",
    "relation",
    "object",
    "target",
    "tail",
    "to",
}


def bad_input(
    message: str, *, code: str = "INVALID_GRAPH_INPUT", hint: str | None = None
) -> InputFormatError:
    return InputFormatError(message, code=code, hint=hint)


def bad_option(
    message: str, *, code: str = "INVALID_GRAPH_OPTION", hint: str | None = None
) -> ConfigurationError:
    return ConfigurationError(message, code=code, hint=hint)


class Graph:
    """Mutable in-memory property graph with unique node ids and unique edge keys."""

    def __init__(self, directed: bool = True) -> None:
        self.directed = directed
        self.nodes: dict[str, dict[str, object]] = {}
        self.edges: dict[tuple[str, str, str], dict[str, object]] = {}
        self.metadata: dict[str, object] = {}

    def key(self, source: str, relation: str, target: str) -> tuple[str, str, str]:
        if not self.directed and target < source:
            source, target = target, source
        return (source, relation, target)

    def add_node(self, node_id: str, properties: Mapping[str, object] | None = None) -> bool:
        """Add a node or merge properties into it (existing keys are kept); True when new."""
        new = node_id not in self.nodes
        current = self.nodes.setdefault(node_id, {})
        for k, v in (properties or {}).items():
            current.setdefault(k, v)
        return new

    def add_edge(
        self,
        source: str,
        relation: str,
        target: str,
        properties: Mapping[str, object] | None = None,
    ) -> str:
        """Add an edge; returns 'added' or 'duplicate' (a duplicate gets a ``count`` and new properties)."""
        self.add_node(source)
        self.add_node(target)
        key = self.key(source, relation, target)
        if key in self.edges:
            edge = self.edges[key]
            edge["count"] = int(edge.get("count", 1)) + 1  # type: ignore[call-overload]
            for k, v in (properties or {}).items():
                edge.setdefault(k, v)
            return "duplicate"
        edge = {"source": source, "relation": relation, "target": target}
        edge.update({k: v for k, v in (properties or {}).items() if k not in edge})
        self.edges[key] = edge
        return "added"

    def to_dict(self) -> dict[str, object]:
        out: dict[str, object] = {}
        out.update(self.metadata)
        out["directed"] = self.directed
        out["nodes"] = [{"id": nid, **props} for nid, props in self.nodes.items()]
        out["edges"] = [dict(e) for e in self.edges.values()]
        return out

    def copy(self) -> Graph:
        clone = Graph(self.directed)
        clone.metadata = dict(self.metadata)
        clone.nodes = {k: dict(v) for k, v in self.nodes.items()}
        clone.edges = {k: dict(v) for k, v in self.edges.items()}
        return clone

    # adjacency ------------------------------------------------------------------
    def adjacency(
        self, direction: str, relations: set[str] | None
    ) -> dict[str, list[tuple[str, dict[str, object]]]]:
        """Neighbour lists in edge insertion order for 'out', 'in' or 'both'."""
        adjacency: dict[str, list[tuple[str, dict[str, object]]]] = {n: [] for n in self.nodes}
        forward = not self.directed or direction in {"out", "both"}
        backward = not self.directed or direction in {"in", "both"}
        for edge in self.edges.values():
            if relations is not None and edge["relation"] not in relations:
                continue
            s, t = str(edge["source"]), str(edge["target"])
            if forward:
                adjacency[s].append((t, edge))
            if backward and not (s == t and forward):
                adjacency[t].append((s, edge))
        return adjacency


# ------------------------------------------------------------------ parsing


def _text_id(value: object, *, what: str) -> str:
    if value is None or isinstance(value, (Mapping, list, tuple, set)):
        raise bad_input(
            f"{what} must be a string or number, got {value!r}.", code="INVALID_NODE_ID"
        )
    text = str(value).strip()
    if not text:
        raise bad_input(f"{what} is empty.", code="INVALID_NODE_ID")
    return text


def _first(mapping: Mapping[str, object], keys: Sequence[str]) -> str | None:
    for key in keys:
        if key in mapping:
            return key
    return None


def _read_text_file(path: Path) -> str:
    try:
        if path.stat().st_size > MAX_FILE_BYTES:
            raise bad_input(
                f"file is larger than {MAX_FILE_BYTES} bytes: {path.name}",
                code="GRAPH_INPUT_TOO_LARGE",
            )
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise bad_input(
            f"could not read UTF-8 file {path}: {exc}", code="GRAPH_INPUT_UNREADABLE"
        ) from exc


def _is_path(value: object) -> bool:
    if isinstance(value, Path):
        return True
    return (
        isinstance(value, str)
        and "\n" not in value
        and "\t" not in value
        and 0 < len(value) < 4096
        and os.path.isfile(value)
    )


class RawGraph:
    """Parsed but not yet de-duplicated graph content."""

    def __init__(self) -> None:
        self.directed: bool | None = None
        self.nodes: list[dict[str, object]] = []
        self.edges: list[dict[str, object]] = []
        self.explicit_nodes = False
        self.metadata: dict[str, object] = {}


def _edge_from(item: object, default_relation: str, where: str) -> dict[str, object]:
    if isinstance(item, Mapping):
        sk, tk = _first(item, _SOURCE_KEYS), _first(item, _TARGET_KEYS)
        if sk is None or tk is None:
            raise bad_input(
                f"{where}: edge needs a source/subject and a target/object field.",
                code="INVALID_EDGE",
            )
        rk = _first(item, [k for k in _RELATION_KEYS if k not in (sk, tk)])
        used = {sk, tk, rk}
        relation = default_relation if rk is None else str(item[rk]).strip()
        edge: dict[str, object] = {
            "source": _text_id(item[sk], what=f"{where} source"),
            "relation": relation,
            "target": _text_id(item[tk], what=f"{where} target"),
        }
        for k, v in item.items():
            if k not in used and k not in edge:
                edge[str(k)] = v
        return edge
    if isinstance(item, Sequence) and not isinstance(item, (str, bytes)):
        if len(item) == 2:
            return {
                "source": _text_id(item[0], what=f"{where} source"),
                "relation": default_relation,
                "target": _text_id(item[1], what=f"{where} target"),
            }
        if len(item) in (3, 4):
            edge = {
                "source": _text_id(item[0], what=f"{where} source"),
                "relation": str(item[1]).strip(),
                "target": _text_id(item[2], what=f"{where} target"),
            }
            if len(item) == 4:
                edge["weight"] = item[3]
            return edge
    raise bad_input(
        f"{where}: expected [subject, relation, object], [source, target] or an edge mapping.",
        code="INVALID_EDGE",
    )


def _parse_delimited(
    text: str, default_relation: str, delimiter: str | None
) -> list[dict[str, object]]:
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []
    if delimiter is None:
        delimiter = "\t" if "\t" in lines[0] else ","
    rows = list(csv.reader(lines, delimiter=delimiter))
    if {c.strip().lower() for c in rows[0]} & _HEADER_WORDS:
        header = [c.strip() for c in rows[0]]
        return [
            _edge_from(dict(zip(header, r, strict=False)), default_relation, f"row {i + 2}")
            for i, r in enumerate(rows[1:])
        ]
    return [
        _edge_from([c.strip() for c in r], default_relation, f"row {i + 1}")
        for i, r in enumerate(rows)
    ]


def _parse_graphml(text: str) -> RawGraph:
    if re.search(r"<!(DOCTYPE|ENTITY)", text, re.I):
        raise bad_input(
            "GraphML with DOCTYPE/ENTITY declarations is rejected for safety.", code="UNSAFE_XML"
        )
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise bad_input(f"invalid GraphML: {exc}", code="INVALID_GRAPHML") from exc

    def tag(element: ET.Element) -> str:
        return element.tag.rsplit("}", 1)[-1]

    keys: dict[str, tuple[str, str]] = {}
    for el in root.iter():
        if tag(el) == "key":
            keys[el.get("id", "")] = (el.get("attr.name", el.get("id", "")), el.get("for", "all"))
    graph_el = next((el for el in root if tag(el) == "graph"), None)
    if graph_el is None:
        raise bad_input("GraphML has no <graph> element.", code="INVALID_GRAPHML")
    raw = RawGraph()
    raw.directed = graph_el.get("edgedefault", "directed") != "undirected"
    raw.explicit_nodes = True

    def data(el: ET.Element) -> dict[str, object]:
        return {
            keys.get(d.get("key", ""), (d.get("key", ""), ""))[0]: (d.text or "")
            for d in el
            if tag(d) == "data"
        }

    for el in graph_el:
        if tag(el) == "node":
            raw.nodes.append({"id": el.get("id", ""), **data(el)})
        elif tag(el) == "edge":
            props = data(el)
            relation = str(props.pop("relation", DEFAULT_RELATION))
            raw.edges.append(
                {
                    "source": el.get("source", ""),
                    "relation": relation,
                    "target": el.get("target", ""),
                    **props,
                }
            )
    return raw


def raw_graph(value: object, options: Mapping[str, object], *, what: str = "input") -> RawGraph:
    """Parse any supported graph input without normalising it.

    Accepted: a graph mapping (``nodes``/``edges``/``triples``), a list of triples ``[s, p, o]``, 2-item
    edges or edge mappings (aliases subject/predicate/object, source/relation/target, head/relation/tail,
    from/type/to), delimited text or a file (.json, .jsonl, .tsv, .csv, .txt, .graphml).
    """

    default_relation = str(options.get("default_relation", DEFAULT_RELATION))
    limit = options.get("max_edges", DEFAULT_MAX_EDGES)
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise bad_option("option 'max_edges' must be a positive integer.")
    raw = RawGraph()
    if _is_path(value):
        path = Path(str(value)).expanduser()
        text = _read_text_file(path)
        suffix = path.suffix.lower()
        if suffix == ".graphml":
            raw = _parse_graphml(text)
        elif suffix == ".json":
            try:
                raw = raw_graph(json.loads(text), options, what=what)
            except json.JSONDecodeError as exc:
                raise bad_input(f"invalid JSON in {path.name}: {exc}") from exc
        elif suffix == ".jsonl":
            items = []
            for number, line in enumerate(text.splitlines(), 1):
                if line.strip():
                    try:
                        items.append(json.loads(line))
                    except json.JSONDecodeError as exc:
                        raise bad_input(f"invalid JSON on line {number} of {path.name}.") from exc
            if any(isinstance(i, Mapping) and "_kind" in i for i in items):
                # files written by knowledge_graph.export(format="jsonl"): one object per node or edge
                nodes = [{k: v for k, v in i.items() if k != "_kind"} for i in items if i.get("_kind") == "node"]
                edges = [{k: v for k, v in i.items() if k != "_kind"} for i in items if i.get("_kind") == "edge"]
                raw = raw_graph({"nodes": nodes, "edges": edges}, options, what=what)
            else:
                raw = raw_graph(items, options, what=what)
        elif suffix in {".tsv", ".csv", ".txt", ".tab", ""}:
            raw.edges = _parse_delimited(text, default_relation, "," if suffix == ".csv" else None)
        else:
            raise bad_input(
                f"unsupported graph file type '{suffix}'.",
                code="UNSUPPORTED_GRAPH_FORMAT",
                hint="Use .json, .jsonl, .tsv, .csv or .graphml.",
            )
    elif isinstance(value, str):
        stripped = value.lstrip()
        if stripped.startswith("<"):
            raw = _parse_graphml(value)
        elif stripped.startswith(("{", "[")):
            try:
                raw = raw_graph(json.loads(value), options, what=what)
            except json.JSONDecodeError as exc:
                raise bad_input(f"invalid JSON text: {exc}") from exc
        else:
            raw.edges = _parse_delimited(value, default_relation, None)
    elif isinstance(value, Mapping):
        if not any(k in value for k in ("nodes", "edges", "triples", "links", "relationships")):
            raise bad_input(
                f"{what} mapping needs 'nodes' and/or 'edges' (or 'triples').",
                code="INVALID_GRAPH_INPUT",
            )
        directed = value.get("directed")
        if directed is not None and not isinstance(directed, bool):
            raise bad_input("'directed' must be a boolean.")
        raw.directed = directed
        for k, v in value.items():
            if k in {"name", "description", "version", "metadata"}:
                raw.metadata[str(k)] = v
        nodes = value.get("nodes")
        if nodes is not None:
            raw.explicit_nodes = True
            if isinstance(nodes, Mapping):
                nodes = [
                    {"id": k, **(v if isinstance(v, Mapping) else {})} for k, v in nodes.items()
                ]
            if not isinstance(nodes, Sequence) or isinstance(nodes, (str, bytes)):
                raise bad_input("'nodes' must be a list.")
            for i, n in enumerate(nodes):
                if isinstance(n, Mapping):
                    if "id" not in n:
                        raise bad_input(f"nodes[{i}] has no 'id'.", code="INVALID_NODE_ID")
                    raw.nodes.append({str(k): v for k, v in n.items()})
                else:
                    raw.nodes.append({"id": n})
        for key in ("edges", "triples", "links", "relationships"):
            items = value.get(key)
            if items is not None:
                if not isinstance(items, Sequence) or isinstance(items, (str, bytes)):
                    raise bad_input(f"'{key}' must be a list.")
                raw.edges.extend(
                    _edge_from(it, default_relation, f"{key}[{i}]") for i, it in enumerate(items)
                )
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        for i, item in enumerate(value):
            raw.edges.append(_edge_from(item, default_relation, f"item {i}"))
            if len(raw.edges) > limit:
                break
    else:
        raise bad_input(
            f"{what} must be triples, an edge list, a graph mapping or a file path.",
            code="INVALID_GRAPH_INPUT",
            hint="Example: [['BRCA1', 'associated_with', 'breast cancer']].",
        )
    if len(raw.edges) > limit:
        raise bad_input(
            f"more than max_edges={limit} edges.",
            code="GRAPH_INPUT_TOO_LARGE",
            hint="Raise options['max_edges'] deliberately or split the input.",
        )
    for edge in raw.edges:
        edge["source"] = _text_id(edge["source"], what="edge source")
        edge["target"] = _text_id(edge["target"], what="edge target")
    return raw


def directed_option(options: Mapping[str, object], raw: RawGraph | None = None) -> bool:
    value = options.get("directed")
    if value is None:
        return True if raw is None or raw.directed is None else raw.directed
    if not isinstance(value, bool):
        raise bad_option("option 'directed' must be a boolean.")
    return value


def build(raw: RawGraph, options: Mapping[str, object]) -> tuple[Graph, dict[str, int]]:
    """Normalise raw content: merge duplicate nodes/edges, create nodes mentioned by edges."""

    graph = Graph(directed_option(options, raw))
    graph.metadata = dict(raw.metadata)
    stats = {"duplicate_nodes": 0, "duplicate_edges": 0, "self_loops": 0, "dropped_self_loops": 0}
    for node in raw.nodes:
        node_id = _text_id(node["id"], what="node id")
        if not graph.add_node(node_id, {k: v for k, v in node.items() if k != "id"}):
            stats["duplicate_nodes"] += 1
    drop_loops = bool(options.get("drop_self_loops", False))
    for edge in raw.edges:
        s, t = str(edge["source"]), str(edge["target"])
        relation = str(edge["relation"]).strip() or DEFAULT_RELATION
        if s == t:
            if drop_loops:
                stats["dropped_self_loops"] += 1
                continue
            stats["self_loops"] += 1
        props = {k: v for k, v in edge.items() if k not in {"source", "relation", "target"}}
        if graph.add_edge(s, relation, t, props) == "duplicate":
            stats["duplicate_edges"] += 1
    return graph, stats


def load_graph(value: object, options: Mapping[str, object], *, what: str = "input") -> Graph:
    """Parse and normalise a graph from any accepted input form."""

    return build(raw_graph(value, options, what=what), options)[0]


def summary(graph: Graph) -> dict[str, object]:
    relations: dict[str, int] = {}
    for e in graph.edges.values():
        relations[str(e["relation"])] = relations.get(str(e["relation"]), 0) + 1
    return {
        "directed": graph.directed,
        "n_nodes": len(graph.nodes),
        "n_edges": len(graph.edges),
        "relation_counts": dict(sorted(relations.items())),
    }


# ------------------------------------------------------------------ helpers


def str_list(options: Mapping[str, object], key: str) -> list[str] | None:
    value = options.get(key)
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if not isinstance(value, Sequence):
        raise bad_option(f"option '{key}' must be a list of strings.")
    return [str(v) for v in value]


def get_int(options: Mapping[str, object], key: str, default: int, low: int, high: int) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise bad_option(f"option '{key}' must be an integer in [{low}, {high}].")
    return value


def get_choice(
    options: Mapping[str, object], key: str, default: str, choices: Sequence[str]
) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or value.lower() not in choices:
        raise bad_option(f"option '{key}' must be one of {list(choices)}.")
    return value.lower()


def write_output(text: str, options: Mapping[str, object]) -> dict[str, object]:
    """Write ``text`` atomically to ``options['output']`` (no silent overwrite) or return it inline."""

    data = text.encode("utf-8")
    info: dict[str, object] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    output = options.get("output")
    if output is None:
        if len(text) > 2_000_000:
            raise bad_option(
                "content is too large to return inline; set options['output'].",
                code="GRAPH_OUTPUT_TOO_LARGE",
            )
        info.update(content=text, output_path=None)
        return info
    if not isinstance(output, (str, os.PathLike)) or not str(output).strip() or "\0" in str(output):
        raise bad_option("option 'output' must be a file path.")
    path = Path(output).expanduser()
    if path.exists() and not options.get("overwrite", False):
        raise bad_option(
            f"output file already exists: {path}",
            code="GRAPH_OUTPUT_EXISTS",
            hint="Set options['overwrite']=True or choose another path.",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(dir=str(path.parent), prefix=".dnakit-", suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
        os.replace(temp_name, path)
    except BaseException:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
        raise
    info["output_path"] = str(path)
    return info


# ------------------------------------------------------------------ operations


def op_build(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph, stats = build(raw_graph(value, options), options)
    return {**summary(graph), **stats, "graph": graph.to_dict()}


def op_create(value: object, options: Mapping[str, object]) -> dict[str, object]:
    raw = raw_graph(value, options) if value is not None else RawGraph()
    extra_nodes = options.get("nodes")
    if extra_nodes is not None:
        if not isinstance(extra_nodes, Sequence) or isinstance(extra_nodes, (str, bytes)):
            raise bad_option("option 'nodes' must be a list of ids or node mappings.")
        raw.nodes.extend(dict(n) if isinstance(n, Mapping) else {"id": n} for n in extra_nodes)
        raw.explicit_nodes = True
    extra_edges = options.get("edges")
    if extra_edges is not None:
        if not isinstance(extra_edges, Sequence) or isinstance(extra_edges, (str, bytes)):
            raise bad_option("option 'edges' must be a list of triples or edge mappings.")
        raw.edges.extend(
            _edge_from(e, str(options.get("default_relation", DEFAULT_RELATION)), f"edges[{i}]")
            for i, e in enumerate(extra_edges)
        )
    graph, stats = build(raw, options)
    for key in ("name", "description"):
        if options.get(key) is not None:
            graph.metadata[key] = str(options[key])
    return {**summary(graph), **stats, "graph": graph.to_dict()}


def _relationships(options: Mapping[str, object]) -> list[dict[str, object]]:
    items = options.get("relationships")
    default_relation = str(options.get("default_relation", DEFAULT_RELATION))
    if items is not None:
        if not isinstance(items, Sequence) or isinstance(items, (str, bytes)):
            raise bad_option("option 'relationships' must be a list.")
        return [
            _edge_from(it, default_relation, f"relationships[{i}]") for i, it in enumerate(items)
        ]
    if "subject" in options and "object" in options:
        props = options.get("properties", {})
        if not isinstance(props, Mapping):
            raise bad_option("option 'properties' must be a mapping.")
        edge = _edge_from(
            {
                "subject": options["subject"],
                "predicate": options.get("predicate", default_relation),
                "object": options["object"],
            },
            default_relation,
            "relationship",
        )
        edge.update(
            {str(k): v for k, v in props.items() if k not in {"source", "relation", "target"}}
        )
        return [edge]
    raise bad_option(
        "give options['subject'], ['predicate'], ['object'] or options['relationships'].",
        hint="Example: options={'subject': 'TP53', 'predicate': 'regulates', 'object': 'MDM2'}.",
    )


def op_add_relationship(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options) if value is not None else Graph(directed_option(options))
    create_missing = options.get("create_missing_nodes", True) is not False
    on_duplicate = get_choice(options, "on_duplicate", "ignore", ("ignore", "error", "update"))
    added, skipped = [], []
    for edge in _relationships(options):
        s, t, r = (
            str(edge["source"]),
            str(edge["target"]),
            str(edge["relation"]).strip() or DEFAULT_RELATION,
        )
        if not create_missing and (s not in graph.nodes or t not in graph.nodes):
            raise bad_input(
                f"relationship {s!r} -{r}-> {t!r} references an unknown node.",
                code="UNKNOWN_NODE",
                hint="Set options['create_missing_nodes']=True or add the nodes first.",
            )
        props = {k: v for k, v in edge.items() if k not in {"source", "relation", "target"}}
        key = graph.key(s, r, t)
        if key in graph.edges:
            if on_duplicate == "error":
                raise bad_input(
                    f"relationship {s!r} -{r}-> {t!r} already exists.", code="DUPLICATE_EDGE"
                )
            if on_duplicate == "update":
                graph.edges[key].update(props)
            skipped.append({"source": s, "relation": r, "target": t})
            continue
        graph.add_edge(s, r, t, props)
        added.append({"source": s, "relation": r, "target": t})
    return {**summary(graph), "added": added, "already_present": skipped, "graph": graph.to_dict()}


def _edge_ref(item: object, what: str) -> tuple[str, str, str]:
    if isinstance(item, Mapping):
        sk, tk, rk = (
            _first(item, _SOURCE_KEYS),
            _first(item, _TARGET_KEYS),
            _first(item, _RELATION_KEYS),
        )
        if sk is None or tk is None or rk is None:
            raise bad_option(f"{what} needs source, relation and target.")
        return str(item[sk]).strip(), str(item[rk]).strip(), str(item[tk]).strip()
    if isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 3:
        return str(item[0]).strip(), str(item[1]).strip(), str(item[2]).strip()
    raise bad_option(f"{what} must be [source, relation, target] or a mapping.")


def op_update(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    ignore_missing = bool(options.get("ignore_missing", False))
    changes = {
        "edges_removed": 0,
        "nodes_removed": 0,
        "node_properties_set": 0,
        "edge_properties_set": 0,
    }
    known = {
        "remove_edges",
        "remove_nodes",
        "set_node_properties",
        "set_edge_properties",
        "ignore_missing",
        "directed",
        "default_relation",
        "max_edges",
        "backend",
        "drop_self_loops",
    }
    unknown = sorted(set(options) - known)
    if unknown:
        raise bad_option(
            f"unknown option(s) {unknown}.",
            code="UNKNOWN_GRAPH_OPTION",
            hint=f"Supported: {sorted(known)}.",
        )
    if not any(
        k in options
        for k in ("remove_edges", "remove_nodes", "set_node_properties", "set_edge_properties")
    ):
        raise bad_option(
            "no change requested.",
            hint="Use remove_edges, remove_nodes, set_node_properties or set_edge_properties.",
        )
    for item in options.get("remove_edges", []) or []:  # type: ignore[union-attr]
        s, r, t = _edge_ref(item, "remove_edges item")
        key = graph.key(s, r, t)
        if key in graph.edges:
            del graph.edges[key]
            changes["edges_removed"] += 1
        elif not ignore_missing:
            raise bad_input(f"edge {s!r} -{r}-> {t!r} does not exist.", code="UNKNOWN_EDGE")
    for node_id in str_list(options, "remove_nodes") or []:
        if node_id not in graph.nodes:
            if ignore_missing:
                continue
            raise bad_input(f"node {node_id!r} does not exist.", code="UNKNOWN_NODE")
        del graph.nodes[node_id]
        changes["nodes_removed"] += 1
        for key in [
            k for k, e in graph.edges.items() if e["source"] == node_id or e["target"] == node_id
        ]:
            del graph.edges[key]
            changes["edges_removed"] += 1
    node_props = options.get("set_node_properties", {}) or {}
    if not isinstance(node_props, Mapping):
        raise bad_option("option 'set_node_properties' must map node ids to property mappings.")
    for node_id, props in node_props.items():
        if str(node_id) not in graph.nodes:
            if ignore_missing:
                continue
            raise bad_input(f"node {node_id!r} does not exist.", code="UNKNOWN_NODE")
        if not isinstance(props, Mapping) or "id" in props:
            raise bad_option(
                "node properties must be a mapping without 'id' (use knowledge_graph.map to rename)."
            )
        graph.nodes[str(node_id)].update({str(k): v for k, v in props.items()})
        changes["node_properties_set"] += len(props)
    for item in options.get("set_edge_properties", []) or []:  # type: ignore[union-attr]
        if not isinstance(item, Mapping) or not isinstance(item.get("properties"), Mapping):
            raise bad_option(
                "set_edge_properties items need source, relation, target and a 'properties' mapping."
            )
        s, r, t = _edge_ref(item, "set_edge_properties item")
        key = graph.key(s, r, t)
        if key not in graph.edges:
            if ignore_missing:
                continue
            raise bad_input(f"edge {s!r} -{r}-> {t!r} does not exist.", code="UNKNOWN_EDGE")
        props = {
            str(k): v
            for k, v in item["properties"].items()
            if k not in {"source", "relation", "target"}
        }
        graph.edges[key].update(props)
        changes["edge_properties_set"] += len(props)
    return {**summary(graph), **changes, "graph": graph.to_dict()}


def op_link(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    relation = str(options.get("relation", "linked_to")).strip() or "linked_to"
    created: list[dict[str, object]] = []
    pairs = options.get("pairs")
    by = options.get("by")
    if (pairs is None) == (by is None):
        raise bad_option(
            "give exactly one of options['pairs'] (list of [a, b]) or options['by'] (node property name)."
        )
    if pairs is not None:
        if not isinstance(pairs, Sequence) or isinstance(pairs, (str, bytes)):
            raise bad_option("option 'pairs' must be a list of [source, target] pairs.")
        for i, pair in enumerate(pairs):
            if not isinstance(pair, Sequence) or isinstance(pair, str) or len(pair) != 2:
                raise bad_option(f"pairs[{i}] must be [source, target].")
            s, t = _text_id(pair[0], what="pair source"), _text_id(pair[1], what="pair target")
            if graph.key(s, relation, t) not in graph.edges:
                graph.add_edge(s, relation, t)
                created.append({"source": s, "relation": relation, "target": t})
    else:
        if not isinstance(by, str) or not by:
            raise bad_option("option 'by' must be a property name.")
        max_group = get_int(options, "max_group_size", 100, 2, 10_000)
        groups: dict[str, list[str]] = {}
        for node_id, props in graph.nodes.items():
            if props.get(by) not in (None, ""):
                groups.setdefault(str(props[by]).strip().casefold(), []).append(node_id)
        for members in groups.values():
            if len(members) < 2:
                continue
            if len(members) > max_group:
                raise bad_input(
                    f"{len(members)} nodes share one '{by}' value (limit {max_group}).",
                    code="GRAPH_GROUP_TOO_LARGE",
                    hint="Raise options['max_group_size'] to allow all pairs inside a group.",
                )
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    if graph.key(members[i], relation, members[j]) not in graph.edges:
                        graph.add_edge(members[i], relation, members[j])
                        created.append(
                            {"source": members[i], "relation": relation, "target": members[j]}
                        )
    return {
        **summary(graph),
        "n_created": len(created),
        "created": created[:1000],
        "graph": graph.to_dict(),
    }


def _mapping_table(value: object) -> dict[str, str]:
    if isinstance(value, Mapping):
        return {str(k).strip(): str(v).strip() for k, v in value.items()}
    if _is_path(value):
        text = _read_text_file(Path(str(value)).expanduser())
        value = text
    if isinstance(value, str):
        rows = [
            r
            for r in csv.reader(io.StringIO(value), delimiter="\t" if "\t" in value else ",")
            if r and not r[0].startswith("#")
        ]
        if (
            rows
            and rows[0][0].strip().lower() in {"from", "source", "old", "id", "key"}
            and len(rows[0]) >= 2
        ):
            rows = rows[1:]
        return {r[0].strip(): r[1].strip() for r in rows if len(r) >= 2}
    if isinstance(value, Sequence):
        table = {}
        for item in value:
            if isinstance(item, Mapping) and "from" in item and "to" in item:
                table[str(item["from"]).strip()] = str(item["to"]).strip()
            elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                table[str(item[0]).strip()] = str(item[1]).strip()
            else:
                raise bad_option("mapping rows must be {'from', 'to'} or [from, to].")
        return table
    raise bad_option("option 'mapping' must be a dict, rows or a two-column table path.")


def op_map(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    if "mapping" not in options and "relation_mapping" not in options:
        raise bad_option(
            "give options['mapping'] (old id -> new id) and/or options['relation_mapping']."
        )
    table = _mapping_table(options["mapping"]) if "mapping" in options else {}
    relation_table = (
        _mapping_table(options["relation_mapping"]) if "relation_mapping" in options else {}
    )
    unmapped = get_choice(options, "unmapped", "keep", ("keep", "drop", "error"))
    new = Graph(graph.directed)
    new.metadata = dict(graph.metadata)
    dropped: set[str] = set()
    missing: list[str] = []
    rename: dict[str, str | None] = {}
    for node_id in graph.nodes:
        if node_id in table:
            rename[node_id] = table[node_id]
        elif unmapped == "keep" or not table:
            rename[node_id] = node_id
        elif unmapped == "drop":
            rename[node_id] = None
            dropped.add(node_id)
        else:
            missing.append(node_id)
    if missing:
        raise bad_input(
            f"{len(missing)} node(s) have no mapping, e.g. {missing[:5]}.", code="UNMAPPED_NODE"
        )
    merged = 0
    for node_id, props in graph.nodes.items():
        target = rename[node_id]
        if target is not None and not new.add_node(target, props):
            merged += 1
    for edge in graph.edges.values():
        s, t = rename[str(edge["source"])], rename[str(edge["target"])]
        if s is None or t is None:
            continue
        relation = relation_table.get(str(edge["relation"]), str(edge["relation"]))
        if s == t and options.get("drop_self_loops", False):
            continue
        new.add_edge(
            s,
            relation,
            t,
            {k: v for k, v in edge.items() if k not in {"source", "relation", "target"}},
        )
    return {
        **summary(new),
        "n_ids_renamed": sum(1 for k, v in rename.items() if v is not None and v != k),
        "n_nodes_merged": merged,
        "n_nodes_dropped": len(dropped),
        "graph": new.to_dict(),
    }


def degree_table(graph: Graph) -> dict[str, tuple[int, int]]:
    deg = {n: [0, 0] for n in graph.nodes}  # in, out
    for e in graph.edges.values():
        deg[str(e["source"])][1] += 1
        deg[str(e["target"])][0] += 1
    return {n: (v[0], v[1]) for n, v in deg.items()}


def components(graph: Graph) -> list[list[str]]:
    """Weakly connected components (insertion order of first node, nodes in insertion order)."""

    parent = {n: n for n in graph.nodes}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in graph.edges.values():
        a, b = find(str(e["source"])), find(str(e["target"]))
        if a != b:
            parent[b] = a
    groups: dict[str, list[str]] = {}
    for n in graph.nodes:
        groups.setdefault(find(n), []).append(n)
    return list(groups.values())


def op_aggregate(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    n, m = len(graph.nodes), len(graph.edges)
    top_k = get_int(options, "top_k", 10, 1, 1000)
    deg = degree_table(graph)
    totals = {k: (v[0] + v[1]) for k, v in deg.items()}  # a self-loop counts twice (in + out)
    values = sorted(totals.values())

    def stat(seq: Sequence[int]) -> dict[str, object]:
        if not seq:
            return {"min": None, "max": None, "mean": None, "median": None}
        s = sorted(seq)
        mid = len(s) // 2
        median = s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2
        return {"min": s[0], "max": s[-1], "mean": sum(s) / len(s), "median": median}

    comps = components(graph)
    type_key = str(options.get("type_property", "type"))
    types: dict[str, int] = {}
    for props in graph.nodes.values():
        label = str(props.get(type_key, "unknown"))
        types[label] = types.get(label, 0) + 1
    possible = n * (n - 1) if graph.directed else n * (n - 1) / 2
    loops = sum(1 for e in graph.edges.values() if e["source"] == e["target"])
    hubs = sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))[:top_k]
    return {
        **summary(graph),
        "degree": {
            "total": stat(values),
            "in": stat([v[0] for v in deg.values()]),
            "out": stat([v[1] for v in deg.values()]),
        },
        "density": (m - loops) / possible if possible else None,
        "n_isolated_nodes": sum(1 for v in totals.values() if v == 0),
        "n_self_loops": loops,
        "n_weakly_connected_components": len(comps),
        "largest_component_size": max((len(c) for c in comps), default=0),
        "node_type_counts": dict(sorted(types.items())),
        "top_hubs": [
            {"id": k, "degree": v, "in_degree": deg[k][0], "out_degree": deg[k][1]} for k, v in hubs
        ],
        "degree_note": "degree = in + out edges of the graph's edge set; a self-loop counts twice",
    }


def op_neighbors(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    node = options.get("node")
    if node is None:
        raise bad_option("option 'node' (the query node id) is required.")
    node = str(node).strip()
    if node not in graph.nodes:
        raise bad_input(f"node {node!r} is not in the graph.", code="UNKNOWN_NODE")
    depth = get_int(options, "depth", 1, 1, 20)
    direction = get_choice(options, "direction", "both", ("out", "in", "both"))
    relations = set(str_list(options, "relations") or []) or None
    limit = get_int(options, "limit", 1000, 1, 1_000_000)
    adjacency = graph.adjacency(direction, relations)
    parent: dict[str, tuple[str, dict[str, object]] | None] = {node: None}
    distance = {node: 0}
    queue: deque[str] = deque([node])
    order: list[str] = []
    while queue:
        current = queue.popleft()
        if distance[current] == depth:
            continue
        for other, edge in adjacency[current]:
            if other not in distance:
                distance[other] = distance[current] + 1
                parent[other] = (current, edge)
                order.append(other)
                queue.append(other)

    def path_to(target: str) -> list[str]:
        out = [target]
        while parent[out[-1]] is not None:
            out.append(parent[out[-1]][0])  # type: ignore[index]
        return out[::-1]

    neighbors = []
    for other in order[:limit]:
        via = parent[other][1]  # type: ignore[index]
        neighbors.append(
            {
                "id": other,
                "distance": distance[other],
                "path": path_to(other),
                "via_edge": {
                    "source": via["source"],
                    "relation": via["relation"],
                    "target": via["target"],
                },
                **{k: v for k, v in graph.nodes[other].items()},
            }
        )
    result: dict[str, object] = {
        "node": node,
        "depth": depth,
        "direction": direction if graph.directed else "both",
        "n_neighbors": len(order),
        "neighbors": neighbors,
        "truncated": len(order) > limit,
    }
    target = options.get("target")
    if target is not None:
        target = str(target).strip()
        if target not in graph.nodes:
            raise bad_input(f"target node {target!r} is not in the graph.", code="UNKNOWN_NODE")
        # unbounded BFS for the shortest path (the depth limit applies to the neighbourhood only)
        max_path = get_int(options, "max_path_length", 50, 1, 1000)
        prev: dict[str, tuple[str, dict[str, object]] | None] = {node: None}
        dist = {node: 0}
        q: deque[str] = deque([node])
        while q and target not in dist:
            cur = q.popleft()
            if dist[cur] >= max_path:
                continue
            for other, edge in adjacency[cur]:
                if other not in dist:
                    dist[other] = dist[cur] + 1
                    prev[other] = (cur, edge)
                    q.append(other)
        if target in dist:
            nodes_path = [target]
            edges_path = []
            while prev[nodes_path[-1]] is not None:
                p, e = prev[nodes_path[-1]]  # type: ignore[misc]
                edges_path.append(
                    {"source": e["source"], "relation": e["relation"], "target": e["target"]}
                )
                nodes_path.append(p)
            result["shortest_path"] = {
                "found": True,
                "length": dist[target],
                "nodes": nodes_path[::-1],
                "edges": edges_path[::-1],
            }
        else:
            result["shortest_path"] = {"found": False, "length": None, "nodes": [], "edges": []}
    if options.get("include_subgraph", False):
        keep = set(distance)
        sub = Graph(graph.directed)
        for nid in graph.nodes:
            if nid in keep:
                sub.add_node(nid, graph.nodes[nid])
        for e in graph.edges.values():
            if (
                e["source"] in keep
                and e["target"] in keep
                and (relations is None or e["relation"] in relations)
            ):
                sub.add_edge(
                    str(e["source"]),
                    str(e["relation"]),
                    str(e["target"]),
                    {k: v for k, v in e.items() if k not in {"source", "relation", "target"}},
                )
        result["subgraph"] = sub.to_dict()
    return result


def op_search(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    query = options.get("query")
    if not isinstance(query, str) or not query.strip():
        raise bad_option("option 'query' (non-empty text) is required.")
    match = get_choice(options, "match", "contains", ("contains", "exact", "prefix"))
    case = bool(options.get("case_sensitive", False))
    fields = str_list(options, "fields")
    node_type = options.get("node_type")
    type_key = str(options.get("type_property", "type"))
    limit = get_int(options, "limit", 100, 1, 100_000)
    q = query.strip() if case else query.strip().casefold()

    def hit(text: str) -> bool:
        t = text if case else text.casefold()
        return t == q if match == "exact" else t.startswith(q) if match == "prefix" else q in t

    matches = []
    for nid, props in graph.nodes.items():
        if node_type is not None and str(props.get(type_key)) != str(node_type):
            continue
        candidates = {
            "id": nid,
            **{k: v for k, v in props.items() if isinstance(v, (str, int, float))},
        }
        matched = [
            k for k, v in candidates.items() if (fields is None or k in fields) and hit(str(v))
        ]
        if matched:
            matches.append({"id": nid, "matched_fields": matched, **props})
    return {
        "query": query,
        "match": match,
        "n_matches": len(matches),
        "matches": matches[:limit],
        "truncated": len(matches) > limit,
    }


def op_read(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph, stats = build(raw_graph(value, options), options)
    return {**summary(graph), **stats, "graph": graph.to_dict()}


def to_graphml(graph: Graph) -> str:
    node_keys = sorted({k for p in graph.nodes.values() for k in p})
    edge_keys = sorted(
        {k for e in graph.edges.values() for k in e if k not in {"source", "target"}}
    )
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
    ]
    for i, k in enumerate(node_keys):
        lines.append(f'  <key id="n{i}" for="node" attr.name={quoteattr(k)} attr.type="string"/>')
    for i, k in enumerate(edge_keys):
        lines.append(f'  <key id="e{i}" for="edge" attr.name={quoteattr(k)} attr.type="string"/>')

    def val(v: object) -> str:
        return escape(v if isinstance(v, str) else json.dumps(v, default=str))

    lines.append(f'  <graph id="G" edgedefault="{"directed" if graph.directed else "undirected"}">')
    for nid, props in graph.nodes.items():
        lines.append(f"    <node id={quoteattr(nid)}>")
        lines.extend(
            f'      <data key="n{node_keys.index(k)}">{val(v)}</data>' for k, v in props.items()
        )
        lines.append("    </node>")
    for i, e in enumerate(graph.edges.values()):
        lines.append(
            f'    <edge id="edge{i}" source={quoteattr(str(e["source"]))} target={quoteattr(str(e["target"]))}>'
        )
        lines.extend(
            f'      <data key="e{edge_keys.index(k)}">{val(v)}</data>'
            for k, v in e.items()
            if k not in {"source", "target"}
        )
        lines.append("    </edge>")
    lines += ["  </graph>", "</graphml>", ""]
    return "\n".join(lines)


def op_export(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    fmt = options.get("format")
    if fmt is None:
        suffix = (
            Path(str(options["output"])).suffix.lower().lstrip(".") if options.get("output") else ""
        )
        fmt = suffix if suffix in {"json", "tsv", "graphml", "jsonl"} else "json"
    if not isinstance(fmt, str) or fmt.lower() not in {"json", "tsv", "graphml", "jsonl"}:
        raise bad_option("option 'format' must be one of ['json', 'tsv', 'graphml', 'jsonl'].")
    fmt = fmt.lower()
    note = None
    if fmt == "json":
        text = json.dumps(graph.to_dict(), indent=2, ensure_ascii=False, default=str) + "\n"
    elif fmt == "jsonl":
        data = graph.to_dict()
        text = "".join(
            json.dumps({"_kind": "node", **n}, ensure_ascii=False, default=str) + "\n"
            for n in data["nodes"]
        )  # type: ignore[union-attr]
        text += "".join(
            json.dumps({"_kind": "edge", **e}, ensure_ascii=False, default=str) + "\n"
            for e in data["edges"]
        )  # type: ignore[union-attr]
    elif fmt == "tsv":
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "relation", "target"])
        for e in graph.edges.values():
            writer.writerow([e["source"], e["relation"], e["target"]])
        text = buffer.getvalue()
        isolated = sum(
            1
            for n in graph.nodes
            if n not in {str(e[k]) for e in graph.edges.values() for k in ("source", "target")}
        )
        note = (
            "TSV holds edge triples only: node/edge properties"
            + (f" and {isolated} isolated node(s)" if isolated else "")
            + " are not written."
        )
    else:
        text = to_graphml(graph)
        note = "GraphML stores all property values as strings."
    info = write_output(text, options)
    result: dict[str, object] = {**summary(graph), "format": fmt, **info}
    if note:
        result["note"] = note
    return result


def op_validate(value: object, options: Mapping[str, object]) -> dict[str, object]:
    raw = raw_graph(value, options)
    errors: list[dict[str, object]] = []
    warnings: list[dict[str, object]] = []

    def add(bucket: list[dict[str, object]], code: str, message: str) -> None:
        if len(bucket) < 100:
            bucket.append({"code": code, "message": message})

    allowed = set(str_list(options, "allowed_relations") or []) or None
    allow_loops = options.get("allow_self_loops", True) is not False
    directed = directed_option(options, raw)
    ids: dict[str, int] = {}
    for node in raw.nodes:
        nid = _text_id(node["id"], what="node id")
        ids[nid] = ids.get(nid, 0) + 1
    for nid, count in ids.items():
        if count > 1:
            add(errors, "DUPLICATE_NODE", f"node id {nid!r} is defined {count} times")
    seen: dict[tuple[str, str, str], int] = {}
    touched: set[str] = set()
    for edge in raw.edges:
        s, r, t = str(edge["source"]), str(edge["relation"]).strip(), str(edge["target"])
        touched.update((s, t))
        if not r:
            add(errors, "MISSING_RELATION", f"edge {s!r} -> {t!r} has an empty relation")
        if allowed is not None and r and r not in allowed:
            add(
                errors,
                "UNKNOWN_RELATION",
                f"relation {r!r} is not in allowed_relations ({s!r} -> {t!r})",
            )
        if raw.explicit_nodes and not options.get("allow_implicit_nodes", False):
            for end in (s, t):
                if end not in ids:
                    add(
                        errors,
                        "DANGLING_EDGE",
                        f"edge {s!r} -{r}-> {t!r} references undefined node {end!r}",
                    )
        if s == t:
            add(errors if not allow_loops else warnings, "SELF_LOOP", f"self-loop on {s!r} ({r})")
        key = (s, r, t) if directed or s <= t else (t, r, s)
        seen[key] = seen.get(key, 0) + 1
    for (s, r, t), count in seen.items():
        if count > 1:
            add(warnings, "DUPLICATE_EDGE", f"edge {s!r} -{r}-> {t!r} appears {count} times")
    isolated = [n for n in ids if n not in touched]
    if isolated:
        add(
            warnings, "ISOLATED_NODE", f"{len(isolated)} node(s) have no edges, e.g. {isolated[:5]}"
        )
    return {
        "valid": not errors,
        "n_nodes": len(set(ids) | touched),
        "n_edges": len(raw.edges),
        "errors": errors,
        "warnings": warnings,
    }


def canonical_form(graph: Graph) -> dict[str, object]:
    """Order-independent representation used for hashing and diffs."""

    nodes = sorted(
        (
            [nid, json.loads(json.dumps(props, sort_keys=True, default=str))]
            for nid, props in graph.nodes.items()
        ),
        key=lambda x: x[0],
    )
    edges = sorted(
        (
            [
                e["source"],
                e["relation"],
                e["target"],
                json.loads(
                    json.dumps(
                        {k: v for k, v in e.items() if k not in {"source", "relation", "target"}},
                        sort_keys=True,
                        default=str,
                    )
                ),
            ]
            for e in graph.edges.values()
        ),
        key=lambda x: (str(x[0]), str(x[1]), str(x[2])),
    )
    return {"directed": graph.directed, "nodes": nodes, "edges": edges}


def version_id(graph: Graph) -> str:
    payload = json.dumps(
        canonical_form(graph), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def op_version(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    vid = version_id(graph)
    result: dict[str, object] = {
        "version_id": vid,
        "short_id": vid[7:19],
        **summary(graph),
        "definition": "SHA-256 of the canonical (sorted) JSON of nodes, edges and their properties",
    }
    if options.get("label") is not None:
        result["label"] = str(options["label"])
    previous = options.get("previous_version")
    if previous is not None:
        result["previous_version"] = str(previous)
        result["changed_since_previous"] = str(previous) != vid
    return result


def op_snapshot(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    vid = version_id(graph)
    snapshot: dict[str, object] = {
        "format": "dnakit.knowledge_graph.snapshot",
        "format_version": 1,
        "version_id": vid,
        "label": str(options["label"]) if options.get("label") is not None else None,
    }
    if options.get("timestamp") is not None:
        snapshot["timestamp"] = str(options["timestamp"])
    ordered = graph.copy()
    ordered.nodes = dict(sorted(ordered.nodes.items()))
    ordered.edges = dict(sorted(ordered.edges.items(), key=lambda kv: kv[0]))
    snapshot["graph"] = ordered.to_dict()
    info = write_output(json.dumps(snapshot, indent=2, ensure_ascii=False) + "\n", options)
    return {
        "version_id": vid,
        **summary(graph),
        "snapshot": None if info.get("output_path") else snapshot,
        **info,
    }


def _pair(value: object, options: Mapping[str, object]) -> tuple[Graph, Graph]:
    if isinstance(value, Mapping) and "old" in value and "new" in value:
        return load_graph(value["old"], options, what="old"), load_graph(
            value["new"], options, what="new"
        )
    if options.get("other") is not None:
        return load_graph(value, options, what="old"), load_graph(
            options["other"], options, what="new"
        )
    if (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == 2
        and all(isinstance(v, (Mapping, Path)) or _is_path(v) for v in value)
    ):
        return load_graph(value[0], options, what="old"), load_graph(value[1], options, what="new")
    raise bad_input(
        "provide {'old': graph, 'new': graph}, a pair [old, new] of graph mappings/files, or input plus options['other'].",
        code="INVALID_GRAPH_INPUT",
    )


def op_compare(value: object, options: Mapping[str, object]) -> dict[str, object]:
    old, new = _pair(value, options)
    if old.directed != new.directed:
        raise bad_input(
            "cannot compare a directed graph with an undirected one.", code="GRAPH_KIND_MISMATCH"
        )
    limit = get_int(options, "max_items", 1000, 1, 1_000_000)
    added_nodes = [n for n in new.nodes if n not in old.nodes]
    removed_nodes = [n for n in old.nodes if n not in new.nodes]
    changed_nodes = [
        {"id": n, "old": old.nodes[n], "new": new.nodes[n]}
        for n in old.nodes
        if n in new.nodes and old.nodes[n] != new.nodes[n]
    ]
    added_edges = [dict(e) for k, e in new.edges.items() if k not in old.edges]
    removed_edges = [dict(e) for k, e in old.edges.items() if k not in new.edges]
    changed_edges = [
        {
            "source": k[0],
            "relation": k[1],
            "target": k[2],
            "old": {
                a: b for a, b in old.edges[k].items() if a not in {"source", "relation", "target"}
            },
            "new": {
                a: b for a, b in new.edges[k].items() if a not in {"source", "relation", "target"}
            },
        }
        for k in old.edges
        if k in new.edges and old.edges[k] != new.edges[k]
    ]
    old_id, new_id = version_id(old), version_id(new)
    counts = {
        "nodes_added": len(added_nodes),
        "nodes_removed": len(removed_nodes),
        "nodes_changed": len(changed_nodes),
        "edges_added": len(added_edges),
        "edges_removed": len(removed_edges),
        "edges_changed": len(changed_edges),
    }
    return {
        "identical": old_id == new_id,
        "old_version_id": old_id,
        "new_version_id": new_id,
        "counts": counts,
        "added_nodes": added_nodes[:limit],
        "removed_nodes": removed_nodes[:limit],
        "changed_nodes": changed_nodes[:limit],
        "added_edges": added_edges[:limit],
        "removed_edges": removed_edges[:limit],
        "changed_edges": changed_edges[:limit],
        "truncated": any(v > limit for v in counts.values()),
    }


def op_lineage(value: object, options: Mapping[str, object]) -> dict[str, object]:
    out_keys = ("output", "outputs", "target", "dataset", "id", "artifact", "name")
    in_keys = ("inputs", "input", "sources", "source", "derived_from", "parents", "from")
    if isinstance(value, Mapping) and "records" in value:
        value = value["records"]
    if isinstance(value, Mapping):
        value = [value]
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise bad_input(
            "input must be a non-empty list of lineage records {'output', 'inputs', 'step'}.",
            code="INVALID_LINEAGE",
        )
    graph = Graph(True)
    steps: dict[str, list[str]] = {}
    for i, rec in enumerate(value):
        if not isinstance(rec, Mapping):
            raise bad_input(f"record {i} is not a mapping.", code="INVALID_LINEAGE")
        ok = _first(rec, out_keys)
        if ok is None:
            raise bad_input(f"record {i} has no output/target/dataset id.", code="INVALID_LINEAGE")
        outputs = (
            rec[ok] if isinstance(rec[ok], Sequence) and not isinstance(rec[ok], str) else [rec[ok]]
        )
        ik = _first(rec, [k for k in in_keys if k != ok])
        raw_in = [] if ik is None else rec[ik]
        inputs = (
            [raw_in]
            if isinstance(raw_in, str)
            else list(raw_in)
            if isinstance(raw_in, Sequence)
            else [raw_in]
        )
        step = next(
            (
                str(rec[k])
                for k in ("step", "tool", "process", "activity")
                if rec.get(k) is not None
            ),
            None,
        )
        extra = {
            k: v
            for k, v in rec.items()
            if k not in {ok, ik, "step", "tool", "process", "activity"}
            and isinstance(v, (str, int, float, bool))
        }
        for out in outputs:
            out_id = _text_id(out, what="output id")
            graph.add_node(out_id, extra)
            for inp in inputs:
                in_id = _text_id(inp, what="input id")
                if in_id == out_id:
                    raise bad_input(f"{out_id!r} is listed as its own input.", code="LINEAGE_CYCLE")
                graph.add_edge(out_id, "derived_from", in_id, {"step": step} if step else None)
                if step:
                    steps.setdefault(out_id, []).append(step)
    indeg = {
        n: 0 for n in graph.nodes
    }  # consumers -> "uses"; compute topological order inputs first
    uses: dict[str, list[str]] = {n: [] for n in graph.nodes}
    for e in graph.edges.values():
        indeg[str(e["source"])] += 1  # outputs depend on inputs
        uses[str(e["target"])].append(str(e["source"]))
    queue = deque(n for n in graph.nodes if indeg[n] == 0)
    order: list[str] = []
    remaining = dict(indeg)
    while queue:
        n = queue.popleft()
        order.append(n)
        for consumer in uses[n]:
            remaining[consumer] -= 1
            if remaining[consumer] == 0:
                queue.append(consumer)
    if len(order) != len(graph.nodes):
        cyclic = sorted(n for n in graph.nodes if n not in set(order))
        raise bad_input(f"lineage contains a cycle involving {cyclic[:10]}.", code="LINEAGE_CYCLE")
    roots = [n for n in graph.nodes if indeg[n] == 0]
    leaves = [n for n in graph.nodes if not uses[n]]
    result: dict[str, object] = {
        **summary(graph),
        "roots": roots,
        "leaves": leaves,
        "topological_order": order,
        "edge_meaning": "source was derived_from target (output -> input)",
        "graph": graph.to_dict(),
    }
    trace = options.get("trace")
    if trace is not None:
        start = str(trace).strip()
        if start not in graph.nodes:
            raise bad_input(f"trace node {start!r} is not in the lineage.", code="UNKNOWN_NODE")
        parents: dict[str, list[str]] = {n: [] for n in graph.nodes}
        for e in graph.edges.values():
            parents[str(e["source"])].append(str(e["target"]))
        seen: dict[str, int] = {}
        queue2: deque[str] = deque([start])
        depth = {start: 0}
        while queue2:
            n = queue2.popleft()
            for p in parents[n]:
                if p not in depth:
                    depth[p] = depth[n] + 1
                    queue2.append(p)
        seen = {k: v for k, v in depth.items() if k != start}
        result["upstream_of"] = {
            "node": start,
            "ancestors": [
                {"id": k, "distance": v}
                for k, v in sorted(seen.items(), key=lambda kv: (kv[1], kv[0]))
            ],
        }
    return result


def op_access_policy(value: object, options: Mapping[str, object]) -> dict[str, object]:
    graph = load_graph(value, options)
    levels = str_list(options, "levels") or ["public", "internal", "restricted", "secret"]
    if len(set(levels)) != len(levels):
        raise bad_option("option 'levels' must list distinct levels from lowest to highest.")
    clearance = str(options.get("clearance", levels[0]))
    if clearance not in levels:
        raise bad_option(
            f"clearance {clearance!r} is not one of the levels {levels}.", code="UNKNOWN_CLEARANCE"
        )
    prop = str(options.get("classification_property", "classification"))
    unclassified = get_choice(options, "unclassified", "deny", ("deny", "allow"))
    deny_nodes = set(str_list(options, "deny_nodes") or [])
    redact = set(str_list(options, "redact_properties") or [])
    rank = {name: i for i, name in enumerate(levels)}

    def allowed(props: Mapping[str, object], what: str) -> bool:
        label = props.get(prop)
        if label in (None, ""):
            return unclassified == "allow"
        if str(label) not in rank:
            raise bad_input(
                f"{what} has unknown classification {label!r} (levels: {levels}).",
                code="UNKNOWN_CLASSIFICATION",
            )
        return rank[str(label)] <= rank[clearance]

    kept = Graph(graph.directed)
    kept.metadata = dict(graph.metadata)
    removed_nodes: list[str] = []
    for nid, props in graph.nodes.items():
        if nid in deny_nodes or not allowed(props, f"node {nid!r}"):
            removed_nodes.append(nid)
        else:
            kept.add_node(nid, {k: v for k, v in props.items() if k not in redact})
    removed_edges = 0
    for e in graph.edges.values():
        if (
            e["source"] not in kept.nodes
            or e["target"] not in kept.nodes
            or not allowed(e, f"edge {e['source']!r}->{e['target']!r}")
        ):
            removed_edges += 1
            continue
        kept.add_edge(
            str(e["source"]),
            str(e["relation"]),
            str(e["target"]),
            {k: v for k, v in e.items() if k not in {"source", "relation", "target"} | redact},
        )
    return {
        **summary(kept),
        "clearance": clearance,
        "levels": levels,
        "n_nodes_removed": len(removed_nodes),
        "removed_nodes": removed_nodes[:1000],
        "n_edges_removed": removed_edges,
        "policy": f"keep items with {prop} <= {clearance}; unclassified items are {'kept' if unclassified == 'allow' else 'removed'}",
        "graph": kept.to_dict(),
    }
