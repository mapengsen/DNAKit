"""FastTree facades (pack P16): distance matrix, and tree inference with starting trees, constraints, bootstrap and ME-only mode.

Runs FastTree (tested 2.2.0, double precision) without a shell through ``dnakit.bioinformatics._p16_cli``. The tree or matrix
goes to a new file; results are JSON-serialisable dictionaries. ``dnakit.phylogeny.external.infer_fast_tree`` covers the
plain run; these functions add the other modes.
"""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.exceptions import ConfigurationError

_TOOL = cli.Tool("fasttree", "FastTree", "Install FastTree (conda install -c bioconda fasttree) and put it on PATH or pass fasttree_path.",
                 "2.2.0", "https://doi.org/10.1371/journal.pone.0009490", ("-expert",), r"FastTree (?:Version )?(\d+(?:\.\d+)+)")
_ENV_NOTE = "FastTree prints its progress on standard error; the tail is returned in log_tail."


def compute_fasttree_distance_matrix(
    alignment: str,
    output_path: str,
    *,
    nucleotide: bool = True,
    raw_distance: bool = False,
    pseudocounts: float | None = None,
    fasttree_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise distance matrix of an alignment in PHYLIP format (``FastTree -makematrix``).

    ``nucleotide=True`` uses ``-nt`` (Jukes-Cantor distances; protein alignments use log-corrected BLOSUM45 distances).
    ``raw_distance`` is ``-rawdist`` (fraction different) and ``pseudocounts`` is ``-pseudo <weight>`` (for sequences
    with little overlap). The matrix is written to ``output_path`` and returned as ``labels`` and ``matrix``.
    """
    source = cli.existing_file(alignment, "alignment")
    target = cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    if cli.flag(nucleotide, "nucleotide"):
        arguments.append("-nt")
    if cli.flag(raw_distance, "raw_distance"):
        arguments.append("-rawdist")
    if pseudocounts is not None:
        arguments += ["-pseudo", str(cli.number(pseudocounts, "pseudocounts", minimum=0.0))]
    arguments += ["-makematrix", str(source)]
    record = cli.invoke(_TOOL, arguments, cwd=target.parent, stdout_to=target, explicit=fasttree_path)
    lines = [ln.split() for ln in target.read_text(encoding="utf-8").splitlines() if ln.strip()]
    count = int(lines[0][0])
    labels = [row[0] for row in lines[1 : count + 1]]
    matrix = [[float(value) for value in row[1:]] for row in lines[1 : count + 1]]
    return cli.finish(record, "compute_fasttree_distance_matrix", labels=labels, matrix=matrix)


def infer_fasttree_tree(
    alignment: str,
    output_path: str,
    *,
    nucleotide: bool = True,
    gtr: bool = False,
    gamma: bool = False,
    starting_tree: str | None = None,
    refine: str | None = None,
    constraints: str | None = None,
    constraint_weight: float | None = None,
    bootstrap_replicates: int | None = None,
    no_support: bool = False,
    seed: int | None = None,
    fasttree_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer or refine a tree with FastTree modes beyond a plain run.

    ``starting_tree`` is a Newick file used as the starting topology (``-intree``; branch lengths are ignored).
    ``refine`` selects what happens next: ``None`` (default FastTree search), ``"mllen"`` (``-mllen``: optimise
    maximum-likelihood branch lengths without topology moves; with a starting tree this is the ML refinement of that tree),
    ``"mllen_fixed_topology"`` (``-mllen -nome``: branch lengths on a fixed topology) or ``"noml"`` (``-noml``: a
    minimum-evolution tree without maximum-likelihood steps). ``constraints`` is an alignment of 0/1 split constraints
    (``-constraints``) with ``constraint_weight`` (``-constraintWeight``, default 100). ``bootstrap_replicates``
    (``-boot N``) adds resampled support values; ``no_support`` (``-nosupport``) omits support values. The Newick tree is
    written to ``output_path`` (``-out``).
    """
    source = cli.existing_file(alignment, "alignment")
    target = cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    if cli.flag(nucleotide, "nucleotide"):
        arguments.append("-nt")
    if cli.flag(gtr, "gtr"):
        arguments.append("-gtr")
    if cli.flag(gamma, "gamma"):
        arguments.append("-gamma")
    if starting_tree is not None:
        arguments += ["-intree", str(cli.existing_file(starting_tree, "starting_tree"))]
    modes = {None: [], "mllen": ["-mllen"], "mllen_fixed_topology": ["-mllen", "-nome"], "noml": ["-noml"]}
    if refine not in modes:
        raise ConfigurationError(f"refine must be one of {[m for m in modes if m]} or None.")
    if refine == "mllen_fixed_topology" and starting_tree is None:
        raise ConfigurationError("refine='mllen_fixed_topology' needs starting_tree.")
    arguments += modes[refine]
    if constraints is not None:
        arguments += ["-constraints", str(cli.existing_file(constraints, "constraints"))]
        if constraint_weight is not None:
            arguments += ["-constraintWeight", str(cli.number(constraint_weight, "constraint_weight", minimum=0.0))]
    elif constraint_weight is not None:
        raise ConfigurationError("constraint_weight needs constraints.")
    if bootstrap_replicates is not None and no_support:
        raise ConfigurationError("bootstrap_replicates and no_support are mutually exclusive.")
    if bootstrap_replicates is not None:
        arguments += ["-boot", str(cli.integer(bootstrap_replicates, "bootstrap_replicates", minimum=1))]
    if cli.flag(no_support, "no_support"):
        arguments.append("-nosupport")
    if seed is not None:
        arguments += ["-seed", str(cli.integer(seed, "seed", minimum=0))]
    arguments += ["-quiet", "-out", str(target), str(source)]
    record = cli.invoke(_TOOL, arguments, cwd=target.parent, expected=[target], explicit=fasttree_path)
    newick = target.read_text(encoding="utf-8").strip()
    return cli.finish(record, "infer_fasttree_tree", refine=refine, newick=newick[:100_000],
                      tips=newick.count(",") + 1 if newick else 0)


__all__ = ["compute_fasttree_distance_matrix", "infer_fasttree_tree"]

_register_public_verbose(globals(), tuple(__all__))
