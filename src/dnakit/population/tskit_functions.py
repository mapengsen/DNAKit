"""tskit wrappers for DNA-level sample export, text tables, genomic intervals and topologies.

The functions here wrap upstream tskit functions that work on sequences, genomic intervals,
text-format tables, PLINK ``.fam`` files and phylogenetic topologies. Where a tree sequence is
needed, it is built from the explicit table mapping used by ``build_tree_sequence()``. Inputs are
local; no coalescent simulation is run and no network is used.

Output files are never overwritten, and tskit is imported only when a function runs.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import (
    _DEFAULT_MAX_ROWS,
    _MAX_ROWS_LIMIT,
    _TSKIT_EXTRA,
    _integer,
    _limit,
    _number,
    _provenance,
    _sample_sets,
    _summary,
    _text,
    _tree_sequence,
)

_MAX_PRECISION = 17
_MAX_LEAVES = 1000
_SEQUENCE_FORMATS = ("fasta", "nexus")
_TEXT_TABLES = ("nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations")
_TREE_SHAPES = ("random_binary", "balanced", "star", "comb")
_FORBIDDEN_LABEL_CHARACTERS = frozenset("(),:;[]'\t\n ")
_TEXT_CONVENTION = (
    "tskit text tables: whitespace or tab delimited with a header row; intervals are half-open "
    "[left, right); metadata columns are omitted from the returned rows"
)


def _new_file(destination: object, name: str) -> Path:
    """Resolve an output file path that must not exist yet; its directory must exist."""

    if not isinstance(destination, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(_text(str(destination), name)).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def _precision(value: object, name: str = "precision") -> int:
    return _integer(value, name, minimum=0, maximum=_MAX_PRECISION)


def _tree_index(value: object, num_trees: int) -> int:
    return _integer(value, "tree_index", minimum=0, maximum=num_trees - 1)


def _jsonable(value: Any) -> Any:
    """Convert numpy values and tuples into plain JSON-serializable Python values."""

    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    return value


def _newick(
    tskit: Any, tree: Any, precision: int, labels: Mapping[int, str], branches: bool
) -> str:
    try:
        root = tree.root
    except ValueError as exc:
        raise ConfigurationError(f"tskit tree has no single root: {exc}") from exc
    return tskit.build_newick(
        tree,
        root=root,
        precision=precision,
        node_labels=dict(labels),
        include_branch_lengths=branches,
    )


def _sample_labels(ts: Any, taxon_labels: object) -> dict[int, str]:
    """Label each sample node; the default label is ``n<node id>``, as in tskit's NEXUS output."""

    samples = [int(node) for node in ts.samples()]
    if taxon_labels is None:
        return {node: f"n{node}" for node in samples}
    if isinstance(taxon_labels, (str, bytes)) or not isinstance(taxon_labels, Sequence):
        raise ConfigurationError("taxon_labels must be a list of strings, one per sample.")
    if len(taxon_labels) != len(samples):
        raise ConfigurationError(
            f"taxon_labels has {len(taxon_labels)} entries for {len(samples)} samples."
        )
    labels = [_text(label, "taxon label") for label in taxon_labels]
    for label in labels:
        if set(label) & _FORBIDDEN_LABEL_CHARACTERS:
            raise ConfigurationError(f"taxon label {label!r} has Newick-reserved characters.")
    if len(set(labels)) != len(labels):
        raise ConfigurationError("taxon labels must be unique.")
    return dict(zip(samples, labels, strict=True))


def _plain_row(row: Any) -> dict[str, Any]:
    return {
        field.name: _jsonable(getattr(row, field.name))
        for field in dataclasses.fields(row)
        if field.name != "metadata" and not field.name.endswith("_offset")
    }


def export_sample_sequences(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    file_format: str = "fasta",
    wrap_width: int = 60,
    missing_data_character: str = "N",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write each sample's haplotype sequence to a FASTA or NEXUS file, using tskit.

    ``tables`` takes the same explicit rows as ``build_tree_sequence()``. Sequences are the
    alignments tskit derives from the sites and mutations; positions without a site are filled
    with ``missing_data_character``. ``wrap_width=0`` disables line wrapping (FASTA only).
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if file_format not in _SEQUENCE_FORMATS:
        raise ConfigurationError(f"file_format must be one of {_SEQUENCE_FORMATS}.")
    width = _integer(wrap_width, "wrap_width", minimum=0)
    if len(_text(missing_data_character, "missing_data_character")) != 1:
        raise ConfigurationError("missing_data_character must be a single character.")
    ts = _tree_sequence(tskit, tables, max_rows)
    target = _new_file(destination, "destination")
    if file_format == "fasta":
        tskit.write_fasta(
            ts,
            str(target),
            wrap_width=width,
            reference_sequence=None,
            missing_data_character=missing_data_character,
        )
    else:
        tskit.write_nexus(
            ts,
            str(target),
            precision=None,
            include_trees=False,
            include_alignments=None,
            reference_sequence=None,
            missing_data_character=missing_data_character,
            node_labels=None,
        )
    return {
        "status": "completed",
        "operation": "export_sample_sequences",
        "file_format": file_format,
        "destination": str(target),
        "bytes": target.stat().st_size,
        "num_samples": int(ts.num_samples),
        "sequence_length": float(ts.sequence_length),
        "num_sites": int(ts.num_sites),
        "provenance": _provenance(
            "tskit",
            "sample haplotypes from tskit alignments; positions without a site are "
            "missing_data_character",
        ),
    }


def export_ms_format(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    print_trees: bool = False,
    precision: int = 4,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the genotypes of a tree sequence in Hudson's ``ms`` output format, using tskit.

    Each row is one sample and each column one segregating site, with ``1`` for the derived
    allele. ``print_trees=True`` also writes the Newick tree; ``precision`` sets the digits of
    the site positions.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if not isinstance(print_trees, bool):
        raise ConfigurationError("print_trees must be true or false.")
    digits = _precision(precision)
    ts = _tree_sequence(tskit, tables, max_rows)
    target = _new_file(destination, "destination")
    with target.open("w", encoding="utf8") as handle:
        tskit.write_ms(
            ts,
            handle,
            print_trees=print_trees,
            precision=digits,
            num_replicates=1,
            write_header=True,
        )
    return {
        "status": "completed",
        "operation": "export_ms_format",
        "destination": str(target),
        "bytes": target.stat().st_size,
        "num_samples": int(ts.num_samples),
        "num_segregating_sites": int(ts.num_sites),
        "provenance": _provenance("tskit", "ms format; one row per sample, one column per site"),
    }


def complement_intervals(
    intervals: Sequence[Sequence[float]],
    *,
    start: float = 0.0,
    end: float,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the half-open ``[left, right)`` intervals in ``[start, end)`` not covered by others.

    ``intervals`` must be sorted, non-overlapping ``[left, right]`` pairs inside the range,
    as tskit requires. An empty list returns the whole range.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    lower = _number(start, "start", minimum=0.0)
    upper = _number(end, "end", positive=True)
    if upper <= lower:
        raise ConfigurationError("end must be greater than start.")
    pairs = _interval_pairs(intervals, lower, upper)
    complement = tskit.negate_intervals(pairs, lower, upper)
    rows = [[float(left), float(right)] for left, right in complement]
    complement_length = sum(right - left for left, right in rows)
    return {
        "status": "completed",
        "operation": "complement_intervals",
        "start": lower,
        "end": upper,
        "intervals": rows,
        "num_intervals": len(rows),
        "complement_length": complement_length,
        "covered_length": upper - lower - complement_length,
        "provenance": _provenance("tskit", "half-open [left, right) intervals within [start, end)"),
    }


def _interval_pairs(intervals: object, lower: float, upper: float) -> list[tuple[float, float]]:
    if isinstance(intervals, (str, bytes)) or not isinstance(intervals, Sequence):
        raise ConfigurationError("intervals must be a list of [left, right] pairs.")
    pairs: list[tuple[float, float]] = []
    previous = lower
    for item in intervals:
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 2:
            raise ConfigurationError("each interval must be a [left, right] pair.")
        left = _number(item[0], "interval left", minimum=lower)
        right = _number(item[1], "interval right", minimum=lower)
        if not left < right:
            raise ConfigurationError(f"interval [{left}, {right}) needs left < right.")
        if right > upper:
            raise ConfigurationError(f"interval [{left}, {right}) extends beyond end={upper}.")
        if left < previous:
            raise ConfigurationError("intervals must be sorted and must not overlap.")
        previous = right
        pairs.append((left, right))
    return pairs


def read_text_table(
    source: str | Path,
    *,
    table: str = "nodes",
    strict: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse one tskit text table file and return its rows, using tskit's text parsers.

    ``table`` selects the parser: ``nodes``, ``edges``, ``sites``, ``mutations``,
    ``individuals``, ``populations`` or ``migrations``. Rows keep the numeric and string
    columns; metadata columns are omitted. ``strict=True`` requires tab-delimited input.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if table not in _TEXT_TABLES:
        raise ConfigurationError(f"table must be one of {_TEXT_TABLES}.")
    if not isinstance(strict, bool):
        raise ConfigurationError("strict must be true or false.")
    limit = _limit(max_rows, "max_rows", _MAX_ROWS_LIMIT)
    path = input_path(source)
    parser = getattr(tskit, f"parse_{table}")
    try:
        with path.open(encoding="utf8") as handle:
            parsed = parser(handle, strict=strict)
    except (tskit.TskitException, ValueError, OSError) as exc:
        raise ConfigurationError(f"tskit could not parse {path.name} as {table}: {exc}") from exc
    count = int(parsed.num_rows)
    if count > limit:
        raise ConfigurationError(f"{path.name} has {count} rows, above max_rows={limit}.")
    rows = [_plain_row(parsed[index]) for index in range(count)]
    columns = [name for name in parsed.column_names if name != "metadata"]
    columns = [name for name in columns if not name.endswith("_offset")]
    return {
        "status": "completed",
        "operation": "read_text_table",
        "table": table,
        "source": str(path),
        "num_rows": count,
        "columns": columns,
        "rows": rows,
        "provenance": _provenance("tskit", _TEXT_CONVENTION),
    }


def read_text_tree_sequence(
    nodes: str | Path,
    edges: str | Path,
    *,
    sites: str | Path | None = None,
    mutations: str | Path | None = None,
    sequence_length: float = 0.0,
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a tree sequence from tskit text table files, using ``tskit.load_text()``.

    ``nodes`` and ``edges`` are required; ``sites`` and ``mutations`` are optional. With
    ``sequence_length=0`` tskit takes the length from the edges. The result is summarized
    like ``read_tree_sequence()``.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    length = _number(sequence_length, "sequence_length", minimum=0.0)
    if not isinstance(strict, bool):
        raise ConfigurationError("strict must be true or false.")
    paths = {
        "nodes": input_path(nodes),
        "edges": input_path(edges),
        "sites": None if sites is None else input_path(sites),
        "mutations": None if mutations is None else input_path(mutations),
    }
    try:
        with ExitStack() as stack:
            handles = {
                name: None if path is None else stack.enter_context(path.open(encoding="utf8"))
                for name, path in paths.items()
            }
            ts = tskit.load_text(
                nodes=handles["nodes"],
                edges=handles["edges"],
                sites=handles["sites"],
                mutations=handles["mutations"],
                sequence_length=length,
                strict=strict,
            )
    except (tskit.TskitException, ValueError, OSError) as exc:
        raise ConfigurationError(f"tskit rejected the text tables: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_text_tree_sequence",
        "sources": {name: None if path is None else str(path) for name, path in paths.items()},
        **_summary(ts),
        "provenance": _provenance("tskit", _TEXT_CONVENTION),
    }


def write_text_tree_sequence(
    tables: Mapping[str, Any],
    destination_dir: str | Path,
    *,
    precision: int = 6,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write nodes, edges, sites and mutations as tskit text tables in a new directory.

    The files are ``nodes.txt``, ``edges.txt``, ``sites.txt`` and ``mutations.txt``, which
    ``read_text_tree_sequence()`` reads back. The directory must not exist yet.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    digits = _precision(precision)
    if not isinstance(destination_dir, (str, Path)):
        raise ConfigurationError("destination_dir must be a path.")
    directory = Path(_text(str(destination_dir), "destination_dir")).expanduser().resolve()
    if directory.exists():
        raise ConfigurationError(f"Output directory already exists: {directory}")
    if not directory.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {directory.parent}")
    ts = _tree_sequence(tskit, tables, max_rows)
    directory.mkdir()
    names = ("nodes", "edges", "sites", "mutations")
    with ExitStack() as stack:
        handles = {
            name: stack.enter_context((directory / f"{name}.txt").open("w", encoding="utf8"))
            for name in names
        }
        tskit.dump_text(
            ts,
            nodes=handles["nodes"],
            edges=handles["edges"],
            sites=handles["sites"],
            mutations=handles["mutations"],
            individuals=None,
            populations=None,
            migrations=None,
            provenances=None,
            precision=digits,
            encoding="utf8",
            base64_metadata=True,
        )
    return {
        "status": "completed",
        "operation": "write_text_tree_sequence",
        "destination_dir": str(directory),
        "files": [f"{name}.txt" for name in names],
        "num_nodes": int(ts.num_nodes),
        "num_edges": int(ts.num_edges),
        "num_sites": int(ts.num_sites),
        "num_mutations": int(ts.num_mutations),
        "provenance": _provenance("tskit", _TEXT_CONVENTION),
    }


def count_sample_topologies(
    tables: Mapping[str, Any],
    sample_sets: Mapping[str, Sequence[int]],
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the embedded topologies of named sample sets in every tree, using tskit.

    ``sample_sets`` maps names to lists of sample node IDs; every set needs at least two
    samples. For each tree, the counts are tskit's topology counters for every combination of
    sets, labelled by the set names joined with ``+``. Each count is keyed by tskit's topology
    rank (``shape`` and ``label``), and entries are sorted by rank.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if sample_sets is None:
        raise ConfigurationError("sample_sets is required, for example {'A': [0, 1], 'B': [2, 3]}.")
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    names = list(sets)
    groups = [sets[name] for name in names]
    trees = []
    for index, counter in enumerate(tskit.combinatorics.treeseq_count_topologies(ts, groups)):
        entries: dict[str, list[dict[str, Any]]] = {}
        for key in sorted(counter.topologies):
            label = "+".join(names[position] for position in key)
            ranks = counter.topologies[key]
            ordered = sorted(ranks.items(), key=lambda item: (item[0].shape, item[0].label))
            entries[label] = [
                {"shape": _jsonable(rank.shape), "label": _jsonable(rank.label), "count": int(n)}
                for rank, n in ordered
            ]
        trees.append({"tree_index": index, "counts": entries})
    return {
        "status": "completed",
        "operation": "count_sample_topologies",
        "sample_sets": {name: len(sets[name]) for name in names},
        "num_trees": int(ts.num_trees),
        "trees": trees,
        "provenance": _provenance("tskit", "tskit treeseq_count_topologies, per tree"),
    }


def resolve_polytomies(
    tables: Mapping[str, Any],
    *,
    tree_index: int = 0,
    random_seed: int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resolve multifurcations of one tree into bifurcations, using tskit, and return Newick.

    ``random_seed`` makes the random resolution reproducible. A tree that is already binary
    is returned unchanged, with ``added_nodes`` equal to 0.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    index = _tree_index(tree_index, ts.num_trees)
    seed = None if random_seed is None else _integer(random_seed, "random_seed", minimum=0)
    tree = ts.at_index(index)
    resolved = tskit.combinatorics.split_polytomies(tree, random_seed=seed, record_provenance=False)
    resolved_ts = resolved.tree_sequence
    return {
        "status": "completed",
        "operation": "resolve_polytomies",
        "tree_index": index,
        "newick": _newick(tskit, resolved, 6, _sample_labels(resolved_ts, None), True),
        "added_nodes": int(resolved_ts.num_nodes - ts.num_nodes),
        "added_edges": int(resolved_ts.num_edges - ts.num_edges),
        "random_seed": seed,
        "provenance": _provenance("tskit", "split_polytomies; default labels n<node id>"),
    }


def tree_to_newick(
    tables: Mapping[str, Any],
    *,
    tree_index: int = 0,
    precision: int = 6,
    include_branch_lengths: bool = True,
    taxon_labels: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return one tree of a tree sequence as a Newick string, using tskit's ``build_newick()``.

    ``taxon_labels`` gives one label per sample in node-ID order (default ``n<node id>``).
    Labels must be unique and must not contain Newick-reserved characters.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    digits = _precision(precision)
    if not isinstance(include_branch_lengths, bool):
        raise ConfigurationError("include_branch_lengths must be true or false.")
    ts = _tree_sequence(tskit, tables, max_rows)
    index = _tree_index(tree_index, ts.num_trees)
    labels = _sample_labels(ts, taxon_labels)
    tree = ts.at_index(index)
    return {
        "status": "completed",
        "operation": "tree_to_newick",
        "tree_index": index,
        "num_samples": int(ts.num_samples),
        "newick": _newick(tskit, tree, digits, labels, include_branch_lengths),
        "provenance": _provenance("tskit", "tskit build_newick with taxon labels"),
    }


def random_tree_newick(
    num_leaves: int,
    *,
    shape: str = "random_binary",
    seed: int | None = None,
    arity: int = 2,
    branch_length: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a tree on ``num_leaves`` taxa of a named shape, using tskit, and return Newick.

    ``shape`` is ``random_binary`` (uniform over leaf-labelled binary trees, reproducible with
    ``seed``), ``balanced`` (``arity`` children per node), ``star`` or ``comb``. Taxa are named
    ``n0`` to ``n<num_leaves - 1>``. Each branch has length ``branch_length``.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if shape not in _TREE_SHAPES:
        raise ConfigurationError(f"shape must be one of {_TREE_SHAPES}.")
    leaves = _integer(num_leaves, "num_leaves", minimum=2, maximum=_MAX_LEAVES)
    length = _number(branch_length, "branch_length", positive=True)
    common: dict[str, Any] = {"span": 1, "branch_length": length, "record_provenance": False}
    combinatorics = tskit.combinatorics
    if shape == "random_binary":
        seed_value = None if seed is None else _integer(seed, "seed", minimum=0)
        tree = combinatorics.generate_random_binary(leaves, random_seed=seed_value, **common)
    else:
        if seed is not None:
            raise ConfigurationError("seed applies only to shape='random_binary'.")
        if shape == "balanced":
            width = _integer(arity, "arity", minimum=2)
            tree = combinatorics.generate_balanced(leaves, arity=width, **common)
        elif shape == "star":
            tree = combinatorics.generate_star(leaves, **common)
        else:
            tree = combinatorics.generate_comb(leaves, **common)
    labels = {int(node): f"n{int(node)}" for node in tree.tree_sequence.samples()}
    return {
        "status": "completed",
        "operation": "random_tree_newick",
        "shape": shape,
        "num_leaves": leaves,
        "seed": seed if shape == "random_binary" else None,
        "newick": _newick(tskit, tree, 6, labels, True),
        "provenance": _provenance("tskit", f"tskit combinatorics generator, shape={shape}"),
    }


def map_allele_indices(
    alleles_from: Sequence[str],
    alleles_to: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map each allele in ``alleles_from`` to its index in ``alleles_to``, using tskit.

    An allele that is missing from ``alleles_to`` gets a new index after the end of
    ``alleles_to``, as tskit's ``allele_remap()`` assigns; those entries are listed in
    ``unmapped_alleles``. ``alleles_to`` must not repeat an allele.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    source = _allele_list(alleles_from, "alleles_from")
    target = _allele_list(alleles_to, "alleles_to")
    if len(set(target)) != len(target):
        raise ConfigurationError("alleles_to must not repeat an allele.")
    indices = [int(value) for value in tskit.genotypes.allele_remap(source, target)]
    unmapped = [
        allele for allele, index in zip(source, indices, strict=True) if index >= len(target)
    ]
    return {
        "status": "completed",
        "operation": "map_allele_indices",
        "indices": indices,
        "num_alleles_to": len(target),
        "unmapped_alleles": unmapped,
        "provenance": _provenance("tskit", "first match in alleles_to; unmapped get new indices"),
    }


def _allele_list(value: object, name: str) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of allele strings.")
    return [_text(allele, f"{name} entry") for allele in value]


def read_plink_fam(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a PLINK ``.fam`` file into individuals, using tskit's ``parse_fam()``.

    The file has five columns: family ID, individual ID, father, mother and sex, with sex
    coded 1 for male and 2 for female. Other sex codes count as unknown. Individuals are
    returned in file order; ``parents`` holds tskit individual indices (-1 when unknown).
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    source = input_path(path)
    try:
        with source.open(encoding="utf8") as handle:
            table = tskit.parse_fam(handle)
    except (tskit.TskitException, ValueError, OSError) as exc:
        raise ConfigurationError(
            f"tskit could not parse {source.name} as PLINK .fam: {exc}"
        ) from exc
    individuals = []
    for index in range(int(table.num_rows)):
        row = table[index]
        meta = row.metadata
        individuals.append(
            {
                "fid": str(meta["plink_fid"]),
                "iid": str(meta["plink_iid"]),
                "sex": int(meta["sex"]),
                "parents": [int(parent) for parent in row.parents],
            }
        )
    males = sum(1 for item in individuals if item["sex"] == 1)
    females = sum(1 for item in individuals if item["sex"] == 2)
    return {
        "status": "completed",
        "operation": "read_plink_fam",
        "source": str(source),
        "num_individuals": len(individuals),
        "num_males": males,
        "num_females": females,
        "num_unknown_sex": len(individuals) - males - females,
        "individuals": individuals,
        "provenance": _provenance("tskit", "PLINK .fam: FID IID PAT MAT SEX; sex 1 male, 2 female"),
    }


__all__ = [
    "complement_intervals",
    "count_sample_topologies",
    "export_ms_format",
    "export_sample_sequences",
    "map_allele_indices",
    "random_tree_newick",
    "read_plink_fam",
    "read_text_table",
    "read_text_tree_sequence",
    "resolve_polytomies",
    "tree_to_newick",
    "write_text_tree_sequence",
]

_register_public_verbose(globals(), tuple(__all__))
