"""Real local implementations behind the knowledge graph capability family.

Graphs are plain triples/edge lists processed in pure Python (no graph database or web service). Building,
merging, querying, statistics, validation, versioning, access filtering and export are native;
``pseudonymize_genomic_metadata_fields`` uses the shared metadata pseudonymisation (key-name and pattern based).
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_ok
from dnakit._metadata_pseudonymize import pseudonymize_envelope
from dnakit.exceptions import ConfigurationError
from dnakit.knowledge_graph import _impl_graph as _graph

_MODULE = 'knowledge_graph'


def _run(operation: str, input: object, options: Mapping[str, object] | None, implementation: object) -> object:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping or None.", code="INVALID_GRAPH_OPTIONS")
    opts = dict(options or {})
    result = implementation(input, opts)  # type: ignore[operator]
    return candidate_ok(_MODULE, operation, input, result, options=opts)


def add_graph_relationship(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Add one or more relationships (edges) to a knowledge graph and return the new graph (the input is not modified).

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``subject``, ``predicate``, ``object`` and optional ``properties`` for one edge, or ``relationships``
    (list of triples / edge mappings); ``create_missing_nodes`` (True; if False unknown nodes raise),
    ``on_duplicate`` ('ignore' | 'error' | 'update'). Input may be None to start from an empty graph.
    """
    return _run('add_graph_relationship', input, options, _graph.op_add_relationship)

def aggregate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Aggregate a knowledge graph into summary statistics: relation counts, degree statistics, density, components, node types and top hubs.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``top_k`` (10 hubs), ``type_property`` ('type'). Degree = in + out edges (a self-loop counts twice);
    density = edges / possible edges excluding self-loops; components are weakly connected.
    """
    return _run('aggregate', input, options, _graph.op_aggregate)

def apply_access_policy(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Apply a clearance-level access policy to a knowledge graph, removing nodes and edges classified above the clearance and redacting properties.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``clearance`` (default the lowest level), ``levels`` (low to high; default ['public', 'internal',
    'restricted', 'secret']), ``classification_property`` ('classification', read from node and edge properties),
    ``unclassified`` ('deny' default: items without a classification are removed | 'allow'), ``deny_nodes`` (ids
    removed regardless), ``redact_properties`` (property names dropped from kept items). Edges touching removed
    nodes are removed. This filters data; it is not an authentication system.
    """
    return _run('apply_access_policy', input, options, _graph.op_access_policy)

def build_data_lineage(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a data-lineage graph (output derived_from input) from provenance records and report roots, leaves and a topological order.

    Input: a list of records ``{'output': id (or 'target' / 'dataset' / 'id'), 'inputs': [ids] (or 'sources' /
    'derived_from'), 'step': name (or 'tool' / 'process')}``. Edges point from each output to its inputs with the
    step stored on the edge. A cycle raises an error. Option ``trace`` (node id) adds all upstream ancestors with
    their distance.
    """
    return _run('build_data_lineage', input, options, _graph.op_lineage)

def build_knowledge_graph(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a knowledge graph (unique nodes, merged duplicate edges) from triples or an edge list.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Nodes mentioned only by edges are created; repeated (source, relation, target) edges are merged and
    counted in ``count``. Returns node/edge counts, relation counts, duplicate and self-loop counts and the
    ``graph`` ({'directed', 'nodes', 'edges'}).
    """
    return _run('build_knowledge_graph', input, options, _graph.op_build)

def compare_database_versions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two versions of a knowledge graph: added, removed and changed nodes and edges.

    Input: ``{'old': graph, 'new': graph}``, a pair ``[old, new]`` of graph mappings or files, or one graph
    plus ``options['other']`` as the new version. Graph forms as in build_knowledge_graph. Edges are identified by
    (source, relation, target); a node or edge is "changed" when its properties differ. Option ``max_items``
    (1000) bounds each list.
    """
    return _run('compare_database_versions', input, options, _graph.op_compare)

def create(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Create a new knowledge graph, optionally seeded with nodes and edges.

    Input: None for an empty graph, or any graph input accepted by build_knowledge_graph. Options: ``name``,
    ``description``, ``directed`` (True), ``nodes`` (list of ids or {'id': ..., properties}), ``edges`` (list of
    triples or edge mappings), ``default_relation``. Returns the graph and its counts.
    """
    return _run('create', input, options, _graph.op_create)

def pseudonymize_genomic_metadata_fields(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Pseudonymise sensitive metadata fields (redact, drop or keyed-hash); key-name and pattern based, NOT de-identification.

    Input: a mapping, a list of mappings, or JSON text of those; anything else raises InputFormatError. Nested mappings
    and lists are searched. A key is sensitive when its name, lower-cased with everything except letters and digits
    removed (``Patient_ID`` -> ``patientid``), is in a default list (name, first/last/full name, patient, patient id/name,
    mrn, medical record number, ssn, social security number, email, phone, telephone, mobile, address, street, city, zip,
    postal code, dob, date of birth, birth date, national id, passport, insurance, account number, ip address, device id,
    photo, biometric and a few variants).
    Options: ``mode`` ('redact' default: value -> "[REDACTED]"; 'drop': key removed; 'hash': HMAC-SHA256 token),
    ``salt`` (required for 'hash', string of >= 8 characters, else ConfigurationError), ``digest_chars`` (hash token
    length, 8-64, default 16), ``keys`` (extra sensitive key names), ``keep`` (key names never treated as sensitive),
    ``scan_text`` (True: e-mail, phone and SSN-like substrings of every string value become "[REDACTED:email]" etc.).
    The input is not modified. The result has ``metadata``, ``mode``, ``n_redacted``, ``redacted_paths`` (key paths only,
    never values) and ``warnings``; ``input_summary`` has only ``record_count`` and ``n_keys``. This is NOT a guarantee
    of de-identification: dates, ages over 89, small geographic units and names in free text are not handled.
    """
    return pseudonymize_envelope('knowledge_graph', 'pseudonymize_genomic_metadata_fields', input, options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise a knowledge graph as JSON, JSONL, TSV (edge triples only) or GraphML text, optionally writing it to a file.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``format`` ('json' | 'jsonl' | 'tsv' | 'graphml'; default from the ``output`` suffix, else 'json'),
    ``output`` (path; without it the text is returned under ``content``), ``overwrite`` (False). GraphML stores
    property values as strings; TSV drops properties and isolated nodes.
    """
    return _run('export', input, options, _graph.op_export)

def link(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Link entities: create edges between explicit node pairs, or between all nodes that share the same value of a property.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``relation`` ('linked_to'); exactly one of ``pairs`` (list of [a, b]; missing nodes are created) or
    ``by`` (node property name; nodes whose values are equal ignoring case are linked pairwise, at most
    ``max_group_size`` (100) nodes per group). Existing edges are not duplicated.
    """
    return _run('link', input, options, _graph.op_link)

def map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Map (rename) node identifiers and relation names through a lookup table, merging nodes that map to the same identifier.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``mapping`` (dict old id -> new id, rows {'from', 'to'} / [from, to], or a two-column TSV/CSV path),
    ``relation_mapping`` (same forms, for relation names), ``unmapped`` ('keep' | 'drop' | 'error' for nodes
    absent from ``mapping``). Merged nodes keep the properties of the first node; merged edges are counted.
    """
    return _run('map', input, options, _graph.op_map)

def query_graph_neighbors(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Query the neighbourhood of a node by breadth-first search, optionally with the shortest path to a target and the induced subgraph.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``node`` (required), ``depth`` (1, max 20), ``direction`` ('both' | 'out' | 'in'; ignored for
    undirected graphs), ``relations`` (restrict to these relation names), ``limit`` (1000), ``target`` (also return
    the shortest path from ``node`` under the same direction/relation rules; ``max_path_length`` 50),
    ``include_subgraph`` (False).
    """
    return _run('query_graph_neighbors', input, options, _graph.op_neighbors)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a knowledge graph from a file or text (JSON, JSONL, TSV/CSV triples, GraphML) and return it with its summary.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    GraphML with DOCTYPE/ENTITY declarations is rejected. Files are limited to 100 MB.
    """
    return _run('read', input, options, _graph.op_read)

def search(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Search the nodes of a knowledge graph by text in their id and string/number properties.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``query`` (required), ``match`` ('contains' | 'exact' | 'prefix'), ``case_sensitive`` (False),
    ``fields`` (restrict to these property names, 'id' included), ``node_type`` with ``type_property`` ('type'),
    ``limit`` (100). Matches keep the graph's node order.
    """
    return _run('search', input, options, _graph.op_search)

def snapshot_database_version(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Create a snapshot of a knowledge graph (canonical content plus version id), optionally saved as a JSON file.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``label``, ``timestamp`` (stored only if given, so snapshots stay reproducible), ``output`` (path;
    otherwise the snapshot is returned), ``overwrite`` (False).
    """
    return _run('snapshot_database_version', input, options, _graph.op_snapshot)

def update(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Apply a change set to a knowledge graph: remove edges or nodes (with their edges) and set node or edge properties.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options (at least one required): ``remove_edges`` (list of [source, relation, target]), ``remove_nodes``
    (ids), ``set_node_properties`` ({id: {property: value}}), ``set_edge_properties`` (list of {source, relation,
    target, properties}), ``ignore_missing`` (False: unknown ids raise). Returns the new graph and change counts.
    """
    return _run('update', input, options, _graph.op_update)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate the integrity of a knowledge graph: duplicate nodes, dangling edges, empty or disallowed relations, self-loops, duplicate edges and isolated nodes.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    Options: ``allowed_relations`` (list; others are errors), ``allow_self_loops`` (True; False makes them errors),
    ``allow_implicit_nodes`` (False: when a node list is given, edges to undefined nodes are errors). Returns
    ``valid`` (no errors), ``errors`` and ``warnings`` (at most 100 each).
    """
    return _run('validate', input, options, _graph.op_validate)

def version(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute a content-based version identifier (SHA-256 of the canonical sorted graph) with node and edge counts.

    Input graph: a mapping ``{'nodes': [...], 'edges': [...], 'directed': bool}``, a list of triples
    ``[subject, relation, object]`` / 2-item edges / edge mappings (keys subject|source|head|from, predicate|relation|
    type, object|target|tail|to), delimited triple text, or a .json/.jsonl/.tsv/.csv/.graphml path. Common options:
    ``directed`` (default True; from the input mapping if present), ``default_relation`` ('related_to'),
    ``max_edges`` (500000), ``drop_self_loops`` (False).
    The id does not depend on node or edge order. Options: ``label`` (echoed), ``previous_version`` (a version id;
    adds ``changed_since_previous``).
    """
    return _run('version', input, options, _graph.op_version)

__all__ = ['add_graph_relationship', 'aggregate', 'apply_access_policy', 'build_data_lineage', 'build_knowledge_graph', 'compare_database_versions', 'create', 'export', 'link', 'map', 'pseudonymize_genomic_metadata_fields', 'query_graph_neighbors', 'read', 'search', 'snapshot_database_version', 'update', 'validate', 'version']

__candidate_adapter_operations__ = set(['add_graph_relationship', 'aggregate', 'apply_access_policy', 'build_data_lineage', 'build_knowledge_graph', 'compare_database_versions', 'create', 'export', 'link', 'map', 'pseudonymize_genomic_metadata_fields', 'query_graph_neighbors', 'read', 'search', 'snapshot_database_version', 'update', 'validate', 'version'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
