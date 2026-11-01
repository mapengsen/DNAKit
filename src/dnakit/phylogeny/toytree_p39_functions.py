"""Tree comparison, consensus, inference and local genealogy (P39)."""
from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.bioinformatics._common import require
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.exceptions import ConfigurationError


def _toy() -> Any:
    return require("toytree", "bio-toytree")


def _result(**values: Any) -> dict[str, Any]:
    return {"status": "completed", "backend": "toytree", "version": _toy().__version__, **values}


def toytree_rf_distance(first: str, second: str, *, normalize: bool = False, verbose: bool = True) -> dict[str, Any]:
    """Robinson-Foulds bipartition distance between two Newick trees with the same tips."""
    toy = _toy()
    return _result(distance=float(toy.distance.get_treedist_rf(toy.tree(first), toy.tree(second), normalize=normalize)))


def toytree_consensus(trees: Sequence[str], *, min_frequency: float = 0.5, verbose: bool = True) -> dict[str, Any]:
    """Majority-rule consensus; min_frequency=0.5 keeps splits present in at least half the trees."""
    toy = _toy()
    if isinstance(trees, str) or not trees or not 0 <= min_frequency <= 1:
        raise ConfigurationError("Provide trees and a min_frequency between 0 and 1.")
    result = toy.infer.get_consensus_tree([toy.tree(t) for t in trees], min_freq=min_frequency)
    return _result(newick=result.write())


def _distances(distances: list[list[float]], names: Sequence[str], algorithm: str) -> dict[str, Any]:
    np, pd = require("numpy", "bio-analysis"), require("pandas", "bio-analysis")
    matrix = np.asarray(distances, dtype=float)
    if matrix.shape != (len(names), len(names)) or len(names) < 2 or len(set(names)) != len(names):
        raise ConfigurationError("Provide a square distance matrix and unique labels for at least two tips.")
    if algorithm == "neighbor_joining_tree" and len(names) < 3:
        raise ConfigurationError("toytree neighbor joining requires at least three tips.")
    if not np.isfinite(matrix).all() or (matrix < 0).any() or not np.allclose(matrix, matrix.T) or not np.allclose(np.diag(matrix), 0):
        raise ConfigurationError("Distances must be finite, nonnegative, symmetric and have a zero diagonal.")
    result = getattr(_toy().infer, algorithm)(pd.DataFrame(matrix, index=names, columns=names))
    return _result(newick=result.write())


def toytree_neighbor_joining(distances: list[list[float]], names: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Infer a neighbor-joining tree from labelled pairwise distances for at least three tips."""
    return _distances(distances, names, "neighbor_joining_tree")


def toytree_upgma(distances: list[list[float]], names: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Infer a UPGMA tree from labelled pairwise distances."""
    return _distances(distances, names, "upgma_tree")


def toytree_mrca(tree: str, tips: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Return the index and descendant tip labels of the most recent common ancestor."""
    obj = _toy().tree(tree)
    if isinstance(tips, str) or not tips:
        raise ConfigurationError("tips must be a nonempty list.")
    node = obj.get_mrca_node(*tips)
    return _result(node_index=int(node.idx), tips=sorted(node.get_leaf_names()))


def toytree_monophyly(tree: str, tips: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Test whether selected tips form a monophyletic group in a rooted Newick tree."""
    if isinstance(tips, str) or not tips:
        raise ConfigurationError("tips must be a nonempty list.")
    return _result(monophyletic=bool(_toy().tree(tree).is_monophyletic(*tips)))


def toytree_unique_topologies(trees: Sequence[str], *, include_root: bool = False, verbose: bool = True) -> dict[str, Any]:
    """Count distinct rooted or unrooted topologies in a Newick tree sample."""
    if isinstance(trees, str) or not trees:
        raise ConfigurationError("trees must be a nonempty list.")
    unique = _toy().mtree(list(trees)).get_unique_topologies(include_root=include_root)
    return _result(topologies=[{"newick": t.write(), "count": int(n)} for t, n in unique])


def toytree_local_genealogy(tree_sequence: str, *, index: int | None = None, position: int | None = None,
                          verbose: bool = True) -> dict[str, Any]:
    """Extract one local ToyTree genealogy from a local tskit file, by tree index or genomic site."""
    if (index is None) == (position is None):
        raise ConfigurationError("Supply exactly one of index or position.")
    ts = require("tskit", "bio-popgen").load(str(_h.existing(tree_sequence, "tree_sequence")))
    module = require("toytree.utils.src.toytree_sequence", "bio-toytree")
    obj = module.ToyTreeSequence(ts)
    tree = obj.at_index(_h.integer(index, "index", minimum=0)) if index is not None else obj.at_site(_h.integer(position, "position", minimum=0))
    return _result(newick=tree.write())


def toytree_k80_matrix(distance: float, *, kappa: float = 2.0, verbose: bool = True) -> dict[str, Any]:
    """K80 nucleotide transition matrix in toytree's upstream base order."""
    module = require("toytree.infer.src.maximum_likelihood", "bio-toytree")
    if distance < 0 or kappa <= 0:
        raise ConfigurationError("distance must be nonnegative and kappa positive.")
    return _result(matrix=module.K80().get_p_matrix_function(kappa=kappa)(distance).tolist(), base_order=list(module.BASE_ORDER))


def toytree_tn93_matrix(distance: float, *, kappas: list[float] = (2.0, 2.0),
                       frequencies: list[float] = (0.25, 0.25, 0.25, 0.25), verbose: bool = True) -> dict[str, Any]:
    """TN93 nucleotide transition matrix using upstream parameter and base ordering."""
    module = require("toytree.infer.src.maximum_likelihood", "bio-toytree")
    if distance < 0 or len(kappas) != 2 or min(kappas) <= 0 or len(frequencies) != 4 or min(frequencies) <= 0 or abs(sum(frequencies) - 1) > 1e-7:
        raise ConfigurationError("Require nonnegative distance, two positive kappas and four positive frequencies summing to one.")
    matrix = module.TN93().get_p_matrix_function(kappa_ls=tuple(kappas), pi_ls=tuple(frequencies))(distance)
    return _result(matrix=matrix.tolist(), base_order=list(module.BASE_ORDER))


__all__ = ["toytree_rf_distance", "toytree_consensus", "toytree_neighbor_joining", "toytree_upgma", "toytree_mrca",
           "toytree_monophyly", "toytree_unique_topologies", "toytree_local_genealogy", "toytree_k80_matrix", "toytree_tn93_matrix"]

_register_public_verbose(globals(), tuple(__all__))
