"""dendropy command-line applications as Python functions: format conversion, sumtrees and sumlabels.

``dendropy_convert_format()`` is the ``dendropy-format`` application: its ``convert`` command and the
``to_nexus``, ``to_newick``, ``to_nexml``, ``to_fasta`` and ``to_phylip`` shortcuts, which only fix the
output format. ``dendropy_sumtrees_summary()`` is the ``sumtrees`` consensus or clade-credibility summary,
and ``dendropy_sumlabels_summary()`` is the ``sumlabels`` label transfer. Each function calls the dendropy
library routines that the application calls, with explicit input and output paths and no shell, and returns
a JSON-serialisable dictionary. The splash, citation, usage, description and overwrite-prompt helpers are
command-line infrastructure and are not wrapped.

Outputs are new files: an existing output path is refused rather than overwritten.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path as _input_file
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_DENDROPY_EXTRA = "bio-dendropy"
_GREATER_THAN_HALF_NOTE = (
    "dendropy.utility.constants.GREATER_THAN_HALF equals 0.5 as a float, so the sumtrees default "
    "keeps clades present in exactly half of the trees; pass min_clade_freq=math.nextafter(0.5, 1.0) "
    "for a strict majority."
)

# dendropy-format --from choices -> (dendropy schema, DataSet.get keyword arguments).
_CONVERT_INPUT: Mapping[str, tuple[str, Mapping[str, bool]]] = {
    "fasta": ("fasta", {}),
    "newick": ("newick", {}),
    "nexml": ("nexml", {}),
    "nexus": ("nexus", {}),
    "phylip-strict": ("phylip", {"strict": True, "multispace_delimiter": False, "interleaved": False}),
    "phylip-relaxed-singlespace": (
        "phylip",
        {"strict": False, "multispace_delimiter": False, "interleaved": False},
    ),
    "phylip-relaxed-multispace": (
        "phylip",
        {"strict": False, "multispace_delimiter": True, "interleaved": False},
    ),
    "phylip-strict-interleaved": (
        "phylip",
        {"strict": True, "multispace_delimiter": False, "interleaved": True},
    ),
    "phylip-relaxed-singlespace-interleaved": (
        "phylip",
        {"strict": False, "multispace_delimiter": False, "interleaved": True},
    ),
    "phylip-relaxed-multispace-interleaved": (
        "phylip",
        {"strict": False, "multispace_delimiter": True, "interleaved": True},
    ),
}
# dendropy-format --to choices -> (dendropy schema, DataSet.write keyword arguments).
_CONVERT_OUTPUT: Mapping[str, tuple[str, Mapping[str, bool]]] = {
    "fasta": ("fasta", {}),
    "newick": ("newick", {}),
    "nexml": ("nexml", {}),
    "nexus": ("nexus", {}),
    "phylip": ("phylip", {}),
    "phylip-strict": ("phylip", {"strict": True}),
}
_DATA_TYPES = ("dna", "rna", "standard")
_RECODE = ("gap", "missing")
_SUMTREES_INPUT = ("nexus/newick", "nexus", "newick")
_SUMTREES_TARGETS = ("consensus", "mcc", "msct")
_TREE_OUTPUT = ("newick", "nexus")


def _dendropy() -> Any:
    return require("dendropy", _DENDROPY_EXTRA)


def _version() -> str | None:
    try:
        return metadata.version("dendropy")
    except metadata.PackageNotFoundError:  # pragma: no cover - dendropy is installed when used
        return None


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def _new_output(value: object, name: str = "output_path") -> Path:
    """A file that does not exist yet, in a directory that does; nothing is overwritten."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _fraction(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    number = float(value)
    if not 0.0 < number <= 1.0:
        raise ConfigurationError(f"{name} must be in (0, 1].")
    return number


def _sorted_clade(node: Any) -> list[str]:
    return sorted(leaf.taxon.label for leaf in node.leaf_iter() if leaf.taxon is not None)


def _tree_text(tree: Any) -> str:
    return tree.as_string(schema="newick", suppress_rooting=False)


# ---------------------------------------------------------------------------
# dendropy-format: convert and the to_* shortcuts
# ---------------------------------------------------------------------------


def dendropy_convert_format(
    source_path: str,
    output_path: str,
    *,
    input_format: str,
    output_format: str,
    data_type: str | None = None,
    unquoted_underscores: bool = False,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence or tree file between formats, as ``dendropy-format convert`` does.

    The ``to_nexus``, ``to_newick``, ``to_nexml``, ``to_fasta`` and ``to_phylip`` commands of the
    application are this function with ``output_format`` fixed. ``data_type`` applies to FASTA and
    PHYLIP input. ``recode_uncertain`` ("gap" or "missing") replaces polymorphic and ambiguous DNA
    states with a gap or a missing-data symbol, as the application's option does.
    """
    dendropy = _dendropy()
    source = _input_file(source_path)
    target = _new_output(output_path)
    in_schema, read_options = _CONVERT_INPUT[
        _choice(input_format, tuple(_CONVERT_INPUT), "input_format")
    ]
    out_schema, write_options = _CONVERT_OUTPUT[
        _choice(output_format, tuple(_CONVERT_OUTPUT), "output_format")
    ]
    read_kwargs: dict[str, Any] = dict(read_options)
    if in_schema in ("fasta", "phylip"):
        if data_type is not None:
            _choice(data_type, _DATA_TYPES, "data_type")
        read_kwargs["data_type"] = data_type
    write_kwargs: dict[str, Any] = dict(write_options)
    underscores = _flag(unquoted_underscores, "unquoted_underscores")
    if underscores and out_schema in ("nexus", "newick"):
        write_kwargs["unquoted_underscores"] = True
    with source.open(encoding="utf-8") as src:
        dataset = dendropy.DataSet.get(file=src, schema=in_schema, **read_kwargs)
    if recode_uncertain is not None:
        alphabet = dendropy.DNA_STATE_ALPHABET
        _choice(recode_uncertain, _RECODE, "recode_uncertain")
        convert_to = alphabet.gap if recode_uncertain == "gap" else alphabet.missing
        convert_from = set(alphabet._polymorphic_states) | set(alphabet._ambiguous_states)
        for char_matrix in dataset.char_matrices:
            for taxon in char_matrix:
                seq = char_matrix[taxon]
                for index, symbol in enumerate(seq):
                    if symbol in convert_from:
                        seq[index] = convert_to
    with target.open("w", encoding="utf-8") as dest:
        dataset.write(file=dest, schema=out_schema, **write_kwargs)
    labels = {taxon.label for namespace in dataset.taxon_namespaces for taxon in namespace}
    first = dataset.char_matrices[0] if dataset.char_matrices else None
    return {
        "status": "completed",
        "operation": "convert_format",
        "input_format": input_format,
        "output_format": output_format,
        "outputs": {"output": str(target)},
        "tree_list_count": len(dataset.tree_lists),
        "char_matrix_count": len(dataset.char_matrices),
        "taxon_count": len(labels),
        "sequence_length": getattr(first, "sequence_size", None) if first is not None else None,
        "provenance": {
            "package": "dendropy",
            "version": _version(),
            "application": "dendropy-format",
            "library_calls": ["dendropy.DataSet.get", "dendropy.DataSet.write"],
        },
    }


# ---------------------------------------------------------------------------
# sumtrees: consensus or clade-credibility summary of a tree file
# ---------------------------------------------------------------------------


def dendropy_sumtrees_summary(
    input_path: str,
    output_path: str,
    *,
    input_format: str = "nexus/newick",
    burnin: int = 0,
    summary_target: str = "consensus",
    min_clade_freq: float | None = None,
    is_source_trees_rooted: bool = False,
    preserve_underscores: bool = False,
    output_format: str = "newick",
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a file of sampled trees as ``sumtrees`` does: consensus, MCC or MSCT tree.

    ``summary_target`` is "consensus" (majority-rule; ``min_clade_freq`` defaults to the application's
    value), "mcc" (maximum product of clade credibilities) or "msct" (maximum sum). The first
    ``burnin`` trees are skipped. Support values are proportions; edge lengths are means over the
    trees. The result lists each internal clade with its support. ``output_format`` is "newick" or
    "nexus". The application's default ``min_clade_freq`` keeps clades present in exactly half of the
    trees (dendropy's ``GREATER_THAN_HALF`` constant equals 0.5 as a float); pass
    ``math.nextafter(0.5, 1.0)`` for a strict majority.
    """
    dendropy = _dendropy()
    constants = require("dendropy.utility.constants", _DENDROPY_EXTRA)
    source = _input_file(input_path)
    target = _new_output(output_path)
    schema = _choice(input_format, _SUMTREES_INPUT, "input_format")
    summary = _choice(summary_target, _SUMTREES_TARGETS, "summary_target")
    out_format = _choice(output_format, _TREE_OUTPUT, "output_format")
    if isinstance(burnin, bool) or not isinstance(burnin, int) or burnin < 0:
        raise ConfigurationError("burnin must be a non-negative integer.")
    _flag(is_source_trees_rooted, "is_source_trees_rooted")
    _flag(preserve_underscores, "preserve_underscores")
    freq = None if min_clade_freq is None else _fraction(min_clade_freq, "min_clade_freq")
    trees = dendropy.TreeArray(taxon_namespace=dendropy.TaxonNamespace())
    trees.read_from_files(
        files=[str(source)],
        schema=schema,
        rooting=dendropy.get_rooting_argument(is_rooted=is_source_trees_rooted),
        tree_offset=burnin,
        store_tree_weights=False,
        preserve_underscores=preserve_underscores,
        ignore_unrecognized_keyword_arguments=True,
    )
    tree_count = int(trees.total_trees_counted)
    if tree_count == 0:
        raise ConfigurationError("No trees remain after burn-in; the input has too few trees.")
    if summary == "consensus":
        used = constants.GREATER_THAN_HALF if freq is None else freq
        tree = trees.consensus_tree(min_freq=used, summarize_splits=False)
    elif summary == "mcc":
        used = None
        tree = trees.maximum_product_of_split_support_tree(
            include_external_splits=False, summarize_splits=False
        )
    else:
        used = None
        tree = trees.maximum_sum_of_split_support_tree(
            include_external_splits=False, summarize_splits=False
        )
    tree.encode_bipartitions()
    if freq is not None and summary != "consensus":
        trees.collapse_edges_with_less_than_minimum_support(tree=tree, min_freq=freq)
    trees.summarize_splits_on_tree(
        tree=tree,
        is_bipartitions_updated=True,
        support_as_percentages=False,
        support_label_decimals=4,
        set_support_as_node_label=True,
        add_support_as_node_attribute=True,
        add_support_as_node_annotation=True,
        set_edge_lengths="mean",
        error_on_negative_edge_lengths=False,
        add_edge_length_summaries_as_edge_attributes=True,
        add_edge_length_summaries_as_edge_annotations=True,
        add_node_age_summaries_as_node_attributes=False,
        add_node_age_summaries_as_node_annotations=False,
    )
    splits = []
    for node in tree.internal_nodes():
        if node is tree.seed_node or node.is_leaf():
            continue
        support = getattr(node, "support", None)
        splits.append(
            {
                "clade": _sorted_clade(node),
                "support": None if support is None else float(support),
            }
        )
    if out_format == "newick":
        target.write_text(_tree_text(tree) + "\n", encoding="utf-8")
    else:
        dataset = dendropy.DataSet()
        dataset.new_tree_list(taxon_namespace=trees.taxon_namespace).append(tree)
        with target.open("w", encoding="utf-8") as dest:
            dataset.write(file=dest, schema="nexus")
    return {
        "status": "completed",
        "operation": "sumtrees_summary",
        "summary_target": summary,
        "tree_count": tree_count,
        "min_clade_freq": used,
        "outputs": {"output": str(target)},
        "newick": _tree_text(tree),
        "splits": splits,
        "provenance": {
            "package": "dendropy",
            "version": _version(),
            "application": "sumtrees",
            "library_calls": [
                "dendropy.TreeArray.read_from_files",
                "dendropy.TreeArray.consensus_tree"
                if summary == "consensus"
                else "dendropy.TreeArray.maximum_*_split_support_tree",
                "dendropy.TreeArray.summarize_splits_on_tree",
            ],
            "note": _GREATER_THAN_HALF_NOTE,
        },
    }


# ---------------------------------------------------------------------------
# sumlabels: transfer labels from source trees onto a target tree
# ---------------------------------------------------------------------------


def dendropy_sumlabels_summary(
    target_path: str,
    source_paths: Sequence[str],
    output_path: str,
    *,
    rooted: bool = False,
    separator: str = "/",
    preserve_target_labels: bool = False,
    output_format: str = "nexus",
    verbose: bool = True,
) -> dict[str, Any]:
    """Label each split of a target tree with the labels the same splits carry in source trees.

    This is ``sumlabels``: the first tree of ``target_path`` is the target; every labelled edge of
    every source tree adds its label to the matching split (labels are joined with ``separator``).
    Splits that match no source are left unlabelled. ``output_format`` is "nexus" or "newick".
    """
    dendropy = _dendropy()
    target_file = _input_file(target_path)
    sources = [_input_file(path) for path in _sources(source_paths)]
    target = _new_output(output_path)
    out_format = _choice(output_format, ("nexus", "newick"), "output_format")
    _flag(rooted, "rooted")
    _flag(preserve_target_labels, "preserve_target_labels")
    separator = _text(separator, "separator")
    rooting = "force-rooted" if rooted else None
    namespace = dendropy.TaxonNamespace()
    target_tree = None
    for tree in dendropy.Tree.yield_from_files(
        files=[str(target_file)], schema="nexus/newick", taxon_namespace=namespace, rooting=rooting
    ):
        target_tree = tree
        break
    if target_tree is None:
        raise ConfigurationError("The target file contains no tree.")
    source_labels: dict[Any, list[str]] = {}
    source_tree_count = 0
    for source in sources:
        for tree in dendropy.Tree.yield_from_files(
            files=[str(source)], schema="nexus/newick", taxon_namespace=namespace, rooting=rooting
        ):
            source_tree_count += 1
            tree.encode_bipartitions()
            for bipartition, edge in tree.bipartition_edge_map.items():
                label = edge.head_node.label
                if label:
                    source_labels.setdefault(bipartition, []).append(label)
    target_tree.encode_bipartitions()
    for bipartition, edge in target_tree.bipartition_edge_map.items():
        collected: list[str] = []
        if preserve_target_labels and edge.head_node.label:
            collected.append(edge.head_node.label)
        elif not preserve_target_labels:
            edge.head_node.label = None
        collected.extend(source_labels.get(bipartition, []))
        if collected:
            edge.head_node.label = separator.join(collected)
    labelled = []
    for edge in target_tree.bipartition_edge_map.values():
        if edge.head_node.label:
            labelled.append({"clade": _sorted_clade(edge.head_node), "label": edge.head_node.label})
    dataset = dendropy.DataSet()
    dataset.new_tree_list(taxon_namespace=namespace).append(target_tree)
    with target.open("w", encoding="utf-8") as dest:
        if out_format == "newick":
            dataset.write(
                file=dest,
                schema="newick",
                suppress_rooting=False,
                suppress_edge_lengths=False,
                unquoted_underscores=False,
                preserve_spaces=False,
                store_tree_weights=False,
            )
        else:
            dataset.write(file=dest, schema="nexus")
    return {
        "status": "completed",
        "operation": "sumlabels_summary",
        "source_tree_count": source_tree_count,
        "labelled_split_count": len(labelled),
        "labels": labelled,
        "outputs": {"output": str(target)},
        "provenance": {
            "package": "dendropy",
            "version": _version(),
            "application": "sumlabels",
            "library_calls": ["dendropy.Tree.yield_from_files", "dendropy.DataSet.write"],
        },
    }


def _sources(values: object) -> list[str]:
    if isinstance(values, str) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("source_paths must be a non-empty list of file paths.")
    return [_text(value, "source_paths item") for value in values]


__all__ = [
    "dendropy_convert_format",
    "dendropy_sumlabels_summary",
    "dendropy_sumtrees_summary",
]

_register_public_verbose(globals(), tuple(__all__))
