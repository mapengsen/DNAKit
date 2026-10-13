"""Optional DendroPy wrappers for DNA alignment statistics and tree measures on taxa.

- ``count_segregating_sites()``, ``calculate_mean_pairwise_differences()``,
  ``calculate_nucleotide_diversity()``, ``calculate_tajimas_d()``,
  ``calculate_wattersons_theta()``, ``build_derived_state_matrix()`` and
  ``calculate_unfolded_site_frequency_spectrum()`` use ``dendropy.calculate.popgenstat``
  on a DNA alignment given as an ID-to-sequence mapping or a ``DNA`` batch.
- ``calculate_weighted_robinson_foulds_distance()``,
  ``count_false_positive_and_negative_splits()``, ``calculate_colless_tree_imbalance()``,
  ``calculate_sackin_index()``, ``calculate_b1_index()``,
  ``calculate_pybus_harvey_gamma()``, ``calculate_treeness()``,
  ``calculate_parsimony_score()``, ``calculate_reconciliation_discordance()`` and
  ``calculate_monophyletic_partition_discordance()`` use ``dendropy`` on Newick trees
  given as text or as a path to a one-tree Newick file.

Nothing is drawn and no external program or network service is called. Packages are
imported only when a function runs. A missing package raises ``BackendUnavailableError``
with the ``pip install "dnakit[bio-dendropy]"`` hint.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, named_sequence_items, require, sequence_text
from dnakit.core import DNA
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dendropy"
_MAX_FILE_BYTES_LIMIT = 200_000_000
_NORMALIZE_COLLESS = (None, "max", "yule", "pda")
_NORMALIZE_SACKIN = (None, "tips", "yule", "pda")
_ALIGNMENT_CONVENTION = (
    "Sequences are upper-cased, sorted by ID and read as DNA with gap and IUPAC symbols. "
    "Missing-data handling (ignore_uncertain) and the statistic itself follow dendropy."
)
_TREE_CONVENTION = (
    "Trees are read from Newick text or a one-tree Newick file. Tip labels must be unique and "
    "non-empty, underscores are kept as written, and the value is dendropy's."
)


def _int(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 200 or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string of at most 200 characters.")
    return value


def _package() -> str:
    return f"dendropy=={metadata.version('dendropy')}"


def _number(value: object) -> float | None:
    number = float(value)
    return round(number, 12) if math.isfinite(number) else None


def _alignment(sequences: Mapping[str, str] | DNA, *, max_sequences: int, max_length: int):
    """Validated (ID, sequence) pairs sorted by ID; all sequences share one length."""

    _int(max_sequences, "max_sequences", 2, 100_000)
    _int(max_length, "max_length", 1, 10_000_000)
    items = named_sequence_items(sequences)
    if len(items) < 2:
        raise ConfigurationError("At least two sequences are required.")
    if len(items) > max_sequences:
        raise ConfigurationError(f"{len(items)} sequences exceed max_sequences={max_sequences}.")
    pairs = sorted(
        (_name(name, "sequence ID"), sequence_text(value, gaps=True)) for name, value in items
    )
    names = [name for name, _ in pairs]
    if len(set(names)) != len(names):
        raise ConfigurationError("Sequence IDs must be unique.")
    lengths = {len(text) for _, text in pairs}
    if len(lengths) != 1:
        raise ConfigurationError("All sequences must have the same length.")
    if lengths.pop() > max_length:
        raise ConfigurationError(f"Sequences are longer than max_length={max_length}.")
    return pairs


def _matrix(dendropy: Any, pairs: Sequence[tuple[str, str]], namespace: Any = None) -> Any:
    options: dict[str, Any] = {} if namespace is None else {"taxon_namespace": namespace}
    return dendropy.DnaCharacterMatrix.from_dict({name: text for name, text in pairs}, **options)


def _alignment_result(
    operation: str,
    pairs: Sequence[tuple[str, str]],
    fields: Mapping[str, object],
    *,
    ignore_uncertain: bool,
) -> dict[str, object]:
    result: dict[str, object] = {
        "sequence_count": len(pairs),
        "site_count": len(pairs[0][1]),
        "ignore_uncertain": ignore_uncertain,
    }
    result.update(fields)
    return {
        "status": "completed",
        "operation": operation,
        "result": result,
        "provenance": {"package": _package(), "convention": _ALIGNMENT_CONVENTION},
    }


def _tree_source(value: object, max_file_bytes: int) -> tuple[str, Any]:
    """Return ("text", Newick string) or ("path", Path) after size checks."""

    max_bytes = _int(max_file_bytes, "max_file_bytes", 1, _MAX_FILE_BYTES_LIMIT)
    if isinstance(value, str) and ";" in value and "\x00" not in value:
        size = len(value.encode("utf-8"))
        if size > max_bytes:
            raise ConfigurationError(
                f"Newick text has {size} bytes; max_file_bytes is {max_bytes}."
            )
        return "text", value.strip()
    if not isinstance(value, (str, Path)) or "\x00" in str(value):
        raise ConfigurationError("tree must be Newick text containing ';' or a path to a file.")
    path = input_path(value)
    size = path.stat().st_size
    if size > max_bytes:
        raise ConfigurationError(f"Tree file has {size} bytes; max_file_bytes is {max_bytes}.")
    return "path", path


def _tip_labels(tree: Any) -> list[str]:
    labels = []
    for leaf in tree.leaf_node_iter():
        taxon = leaf.taxon
        if taxon is None or not taxon.label:
            raise ConfigurationError("Every tip needs a non-empty label.")
        labels.append(taxon.label)
    if len(set(labels)) != len(labels):
        raise ConfigurationError("Tip labels must be unique.")
    return sorted(labels)


def _read_tree(
    dendropy: Any,
    source: object,
    *,
    namespace: Any = None,
    rooting: str | None = None,
    max_tips: int,
    max_file_bytes: int,
) -> Any:
    _int(max_tips, "max_tips", 2, 100_000)
    kind, value = _tree_source(source, max_file_bytes)
    options: dict[str, Any] = {"schema": "newick", "preserve_underscores": True}
    if namespace is not None:
        options["taxon_namespace"] = namespace
    if rooting is not None:
        options["rooting"] = rooting
    try:
        if kind == "path":
            trees = dendropy.TreeList.get(path=str(value), **options)
        else:
            trees = dendropy.TreeList.get(data=value, **options)
    except Exception as exc:
        raise ConfigurationError(f"dendropy could not parse the Newick tree: {exc}") from exc
    if len(trees) != 1:
        raise ConfigurationError(f"Expected exactly one tree; found {len(trees)}.")
    tree = trees[0]
    count = len(_tip_labels(tree))
    if count < 2:
        raise ConfigurationError("A tree needs at least two tips.")
    if count > max_tips:
        raise ConfigurationError(f"Tree has {count} tips; max_tips is {max_tips}.")
    return tree


def _tree_pair(
    dendropy: Any, source_a: object, source_b: object, *, max_tips: int, max_file_bytes: int
) -> tuple[Any, Any]:
    namespace = dendropy.TaxonNamespace()
    tree_a = _read_tree(
        dendropy,
        source_a,
        namespace=namespace,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
    )
    tree_b = _read_tree(
        dendropy,
        source_b,
        namespace=namespace,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
    )
    if _tip_labels(tree_a) != _tip_labels(tree_b):
        raise ConfigurationError("Both trees must have the same tip labels.")
    return tree_a, tree_b


def _tree_result(
    operation: str, tree: Any, fields: Mapping[str, object], convention: str = _TREE_CONVENTION
) -> dict[str, object]:
    return {
        "status": "completed",
        "operation": operation,
        "result": {"tip_count": len(_tip_labels(tree)), **fields},
        "provenance": {"package": _package(), "convention": convention},
    }


def count_segregating_sites(
    sequences: Mapping[str, str] | DNA,
    *,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Count polymorphic sites in a DNA alignment with dendropy."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    value = popgen.num_segregating_sites(_matrix(dendropy, pairs), ignore_uncertain=flag)
    return _alignment_result(
        "count_segregating_sites",
        pairs,
        {"segregating_sites": int(value)},
        ignore_uncertain=flag,
    )


def calculate_mean_pairwise_differences(
    sequences: Mapping[str, str] | DNA,
    *,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Average number of pairwise differences (Tajima's k) in a DNA alignment."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    value = popgen.average_number_of_pairwise_differences(
        _matrix(dendropy, pairs), ignore_uncertain=flag
    )
    return _alignment_result(
        "calculate_mean_pairwise_differences",
        pairs,
        {"mean_pairwise_differences": _number(value)},
        ignore_uncertain=flag,
    )


def calculate_nucleotide_diversity(
    sequences: Mapping[str, str] | DNA,
    *,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Nucleotide diversity (pi) of a DNA alignment, as computed by dendropy."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    value = popgen.nucleotide_diversity(_matrix(dendropy, pairs), ignore_uncertain=flag)
    return _alignment_result(
        "calculate_nucleotide_diversity",
        pairs,
        {"nucleotide_diversity": _number(value)},
        ignore_uncertain=flag,
    )


def calculate_tajimas_d(
    sequences: Mapping[str, str] | DNA,
    *,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Tajima's D for a DNA alignment; None when there are no segregating sites."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    matrix = _matrix(dendropy, pairs)
    segregating = int(popgen.num_segregating_sites(matrix, ignore_uncertain=flag))
    value: float | None = None
    if segregating > 0:
        try:
            value = _number(popgen.tajimas_d(matrix, ignore_uncertain=flag))
        except ZeroDivisionError:
            value = None
    return _alignment_result(
        "calculate_tajimas_d",
        pairs,
        {"tajimas_d": value, "segregating_sites": segregating},
        ignore_uncertain=flag,
    )


def calculate_wattersons_theta(
    sequences: Mapping[str, str] | DNA,
    *,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Watterson's theta for a DNA alignment (segregating sites divided by a_n, not by length)."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    value = popgen.wattersons_theta(_matrix(dendropy, pairs), ignore_uncertain=flag)
    return _alignment_result(
        "calculate_wattersons_theta",
        pairs,
        {
            "wattersons_theta": _number(value),
            "convention_note": "segregating sites divided by a_n = sum(1/i) for i in 1..n-1",
        },
        ignore_uncertain=flag,
    )


def build_derived_state_matrix(
    sequences: Mapping[str, str] | DNA,
    *,
    ancestral_id: str,
    ignore_uncertain: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Mark each ingroup site 0 (ancestral state) or 1 (derived state) against an ancestral ID."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    ancestor = _name(ancestral_id, "ancestral_id")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    by_id = dict(pairs)
    if ancestor not in by_id:
        raise ConfigurationError(f"ancestral_id {ancestor!r} is not one of the sequence IDs.")
    ingroup = [(name, text) for name, text in pairs if name != ancestor]
    ancestral_matrix = _matrix(dendropy, [(ancestor, by_id[ancestor])])
    ancestral_sequence = next(iter(ancestral_matrix.values()))
    derived = popgen.derived_state_matrix(
        _matrix(dendropy, ingroup),
        ancestral_sequence=ancestral_sequence,
        ignore_uncertain=flag,
    )
    rows = sorted(
        (
            {"id": str(getattr(taxon, "label", taxon)), "states": sequence.symbols_as_string()}
            for taxon, sequence in derived.items()
        ),
        key=lambda row: str(row["id"]),
    )
    return _alignment_result(
        "build_derived_state_matrix",
        ingroup,
        {"ancestral_id": ancestor, "rows": rows},
        ignore_uncertain=flag,
    )


def calculate_unfolded_site_frequency_spectrum(
    sequences: Mapping[str, str] | DNA,
    *,
    ancestral_id: str,
    ignore_uncertain: bool = False,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Unfolded site frequency spectrum of the ingroup, relative to an ancestral ID."""

    popgen = require("dendropy.calculate.popgenstat", _EXTRA)
    dendropy = require("dendropy", _EXTRA)
    flag = _flag(ignore_uncertain, "ignore_uncertain")
    ancestor = _name(ancestral_id, "ancestral_id")
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    by_id = dict(pairs)
    if ancestor not in by_id:
        raise ConfigurationError(f"ancestral_id {ancestor!r} is not one of the sequence IDs.")
    ingroup = [(name, text) for name, text in pairs if name != ancestor]
    ancestral_matrix = _matrix(dendropy, [(ancestor, by_id[ancestor])])
    ancestral_sequence = next(iter(ancestral_matrix.values()))
    spectrum = popgen.unfolded_site_frequency_spectrum(
        _matrix(dendropy, ingroup),
        ancestral_sequence=ancestral_sequence,
        ignore_uncertain=flag,
        pad=True,
    )
    counts_by_derived = {int(key): int(value) for key, value in spectrum.items()}
    total = len(ingroup)
    return _alignment_result(
        "calculate_unfolded_site_frequency_spectrum",
        ingroup,
        {
            "ancestral_id": ancestor,
            "ingroup_count": total,
            "site_frequency_spectrum": [counts_by_derived.get(k, 0) for k in range(total + 1)],
        },
        ignore_uncertain=flag,
    )


def calculate_weighted_robinson_foulds_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Weighted Robinson-Foulds distance between two Newick trees with the same tips."""

    dendropy = require("dendropy", _EXTRA)
    treecompare = require("dendropy.calculate.treecompare", _EXTRA)
    first, second = _tree_pair(
        dendropy, tree_a, tree_b, max_tips=max_tips, max_file_bytes=max_file_bytes
    )
    for tree in (first, second):
        for node in tree.preorder_node_iter():
            if node is not tree.seed_node and node.edge.length is None:
                raise ConfigurationError(
                    "Weighted Robinson-Foulds needs a branch length on every edge."
                )
    distance = treecompare.weighted_robinson_foulds_distance(first, second)
    return _tree_result(
        "calculate_weighted_robinson_foulds_distance",
        first,
        {"distance": _number(distance)},
        convention=(
            "dendropy's weighted Robinson-Foulds distance with edge_weight_attr='length'; "
            "both trees are read as written, so their rooting must match the intended comparison."
        ),
    )


def count_false_positive_and_negative_splits(
    reference_tree: str | Path,
    comparison_tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Bipartitions found only in the comparison tree (false positives) or only in the reference."""

    dendropy = require("dendropy", _EXTRA)
    treecompare = require("dendropy.calculate.treecompare", _EXTRA)
    reference, comparison = _tree_pair(
        dendropy, reference_tree, comparison_tree, max_tips=max_tips, max_file_bytes=max_file_bytes
    )
    false_positives, false_negatives = treecompare.false_positives_and_negatives(
        reference, comparison
    )
    return _tree_result(
        "count_false_positive_and_negative_splits",
        reference,
        {"false_positives": int(false_positives), "false_negatives": int(false_negatives)},
        convention=(
            "false_positives are bipartitions found in the comparison tree but not the reference; "
            "false_negatives are the reverse."
        ),
    )


def calculate_colless_tree_imbalance(
    tree: str | Path,
    *,
    normalize: str | None = "max",
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Colless tree imbalance of a strictly bifurcating tree."""

    dendropy = require("dendropy", _EXTRA)
    measure = require("dendropy.calculate.treemeasure", _EXTRA)
    if normalize not in _NORMALIZE_COLLESS:
        raise ConfigurationError("normalize must be one of None, 'max', 'yule' or 'pda'.")
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    try:
        value = measure.colless_tree_imbalance(source, normalize=normalize)
    except TypeError as exc:
        raise ConfigurationError(f"Colless imbalance needs a bifurcating tree: {exc}") from exc
    return _tree_result(
        "calculate_colless_tree_imbalance",
        source,
        {"colless_imbalance": _number(value), "normalize": normalize or "none"},
    )


def calculate_sackin_index(
    tree: str | Path,
    *,
    normalize: str | None = "tips",
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Sackin index: the sum over tips of their number of ancestors, optionally normalized."""

    dendropy = require("dendropy", _EXTRA)
    measure = require("dendropy.calculate.treemeasure", _EXTRA)
    if normalize not in _NORMALIZE_SACKIN:
        raise ConfigurationError("normalize must be one of None, 'tips', 'yule' or 'pda'.")
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    upstream = True if normalize == "tips" else normalize
    value = measure.sackin_index(source, normalize=upstream)
    return _tree_result(
        "calculate_sackin_index",
        source,
        {"sackin_index": _number(value), "normalize": normalize or "none"},
    )


def calculate_b1_index(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """B1 tree-shape statistic: sum over non-root internal nodes of 1 / (1 + child depth)."""

    dendropy = require("dendropy", _EXTRA)
    measure = require("dendropy.calculate.treemeasure", _EXTRA)
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    return _tree_result(
        "calculate_b1_index",
        source,
        {"b1": _number(measure.B1(source))},
        convention=(
            "Value is the sum of 1/m over non-root internal nodes, where m is one more than the "
            "largest m of the node's children and tips have m = 0, as dendropy computes it."
        ),
    )


def calculate_pybus_harvey_gamma(
    tree: str | Path,
    *,
    precision: float = 1e-5,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Pybus and Harvey gamma statistic; the tree must be ultrametric within precision."""

    dendropy = require("dendropy", _EXTRA)
    measure = require("dendropy.calculate.treemeasure", _EXTRA)
    if isinstance(precision, bool) or not isinstance(precision, (int, float)):
        raise ConfigurationError("precision must be a positive number.")
    if not math.isfinite(float(precision)) or not 0 < float(precision) <= 1:
        raise ConfigurationError("precision must be in (0, 1].")
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    try:
        value = measure.pybus_harvey_gamma(source, prec=float(precision))
    except ValueError as exc:
        message = str(exc).splitlines()[0]
        raise ConfigurationError(
            f"Pybus-Harvey gamma needs an ultrametric tree: {message}"
        ) from exc
    return _tree_result("calculate_pybus_harvey_gamma", source, {"gamma": _number(value)})


def calculate_treeness(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Treeness: internal branch length as a fraction of total branch length."""

    dendropy = require("dendropy", _EXTRA)
    measure = require("dendropy.calculate.treemeasure", _EXTRA)
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    for node in source.preorder_node_iter():
        if node is not source.seed_node and node.edge.length is None:
            raise ConfigurationError("Treeness needs a branch length on every edge.")
    try:
        value = measure.treeness(source)
    except ZeroDivisionError as exc:
        raise ConfigurationError("Treeness needs a tree with positive total length.") from exc
    return _tree_result("calculate_treeness", source, {"treeness": _number(value)})


def calculate_parsimony_score(
    tree: str | Path,
    sequences: Mapping[str, str] | DNA,
    *,
    gaps_as_missing: bool = True,
    max_sequences: int = 1000,
    max_length: int = 100_000,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Fitch parsimony score of DNA characters on a tree whose tips match the sequence IDs."""

    dendropy = require("dendropy", _EXTRA)
    treescore = require("dendropy.calculate.treescore", _EXTRA)
    flag = _flag(gaps_as_missing, "gaps_as_missing")
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    pairs = _alignment(sequences, max_sequences=max_sequences, max_length=max_length)
    if [name for name, _ in pairs] != _tip_labels(source):
        raise ConfigurationError("Sequence IDs must match the tree tips exactly.")
    matrix = _matrix(dendropy, pairs, namespace=source.taxon_namespace)
    score = treescore.parsimony_score(source, matrix, gaps_as_missing=flag)
    return {
        "status": "completed",
        "operation": "calculate_parsimony_score",
        "result": {
            "tip_count": len(pairs),
            "site_count": len(pairs[0][1]),
            "parsimony_score": int(score),
            "gaps_as_missing": flag,
        },
        "provenance": {
            "package": _package(),
            "convention": (
                "Fitch parsimony: the minimum number of state changes on the tree, summed over "
                "sites. Gaps count as missing data when gaps_as_missing is true."
            ),
        },
    }


def calculate_reconciliation_discordance(
    gene_tree: str | Path,
    species_tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Gene duplications implied by reconciling a gene tree on a species tree (same tips)."""

    dendropy = require("dendropy", _EXTRA)
    reconcile = require("dendropy.model.reconcile", _EXTRA)
    namespace = dendropy.TaxonNamespace()
    gene = _read_tree(
        dendropy,
        gene_tree,
        namespace=namespace,
        rooting="force-rooted",
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
    )
    species = _read_tree(
        dendropy,
        species_tree,
        namespace=namespace,
        rooting="force-rooted",
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
    )
    if _tip_labels(gene) != _tip_labels(species):
        raise ConfigurationError("Gene and species trees must have the same tip labels.")
    gene.encode_bipartitions()
    species.encode_bipartitions()
    value = reconcile.reconciliation_discordance(gene, species)
    return _tree_result(
        "calculate_reconciliation_discordance",
        gene,
        {"duplications": int(value)},
        convention=(
            "Both trees are read as rooted. The count is the number of gene duplications implied "
            "by the gene tree reconciled on the species tree (Goodman et al. 1979; Maddison 1997)."
        ),
    )


def calculate_monophyletic_partition_discordance(
    tree: str | Path,
    groups: Sequence[Sequence[str]],
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Deep coalescences if the given tip groups were each forced to be monophyletic."""

    dendropy = require("dendropy", _EXTRA)
    reconcile = require("dendropy.model.reconcile", _EXTRA)
    source = _read_tree(dendropy, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    if not isinstance(groups, Sequence) or isinstance(groups, str) or not groups:
        raise ConfigurationError("groups must be a non-empty sequence of tip-label groups.")
    cleaned = []
    for group in groups:
        if not isinstance(group, Sequence) or isinstance(group, str) or not group:
            raise ConfigurationError("Each group must be a non-empty sequence of tip labels.")
        cleaned.append([_name(label, "tip label") for label in group])
    flat = [label for group in cleaned for label in group]
    if len(set(flat)) != len(flat) or sorted(set(flat)) != _tip_labels(source):
        raise ConfigurationError("groups must partition the tree tips: each tip exactly once.")
    namespace = source.taxon_namespace
    partition = dendropy.TaxonNamespacePartition(
        taxon_namespace=namespace,
        membership_lists=[[namespace.get_taxon(label) for label in group] for group in cleaned],
    )
    source.encode_bipartitions()
    value = reconcile.monophyletic_partition_discordance(source, partition)
    return _tree_result(
        "calculate_monophyletic_partition_discordance",
        source,
        {"deep_coalescences": int(value), "group_count": len(cleaned)},
        convention=(
            "Number of deep coalescences needed if each tip group formed one monophyletic clade "
            "(dendropy monophyletic_partition_discordance)."
        ),
    )


__all__ = [
    "build_derived_state_matrix",
    "calculate_b1_index",
    "calculate_colless_tree_imbalance",
    "calculate_mean_pairwise_differences",
    "calculate_monophyletic_partition_discordance",
    "calculate_nucleotide_diversity",
    "calculate_parsimony_score",
    "calculate_pybus_harvey_gamma",
    "calculate_reconciliation_discordance",
    "calculate_sackin_index",
    "calculate_tajimas_d",
    "calculate_treeness",
    "calculate_unfolded_site_frequency_spectrum",
    "calculate_wattersons_theta",
    "calculate_weighted_robinson_foulds_distance",
    "count_false_positive_and_negative_splits",
    "count_segregating_sites",
]

_register_public_verbose(globals(), tuple(__all__))
