"""Biopython phylogeny, NEXUS, tree-query and PAML operations added by integration pack P19.

The functions wrap ``Bio.Nexus`` (character-matrix editing and the NEXUS tree class), ``Bio.Phylo``
(tree queries, random trees, parsimony and consensus), ``Bio.codonalign`` / ``Bio.Align.analysis``
(dN/dS) and the PAML wrappers ``Bio.Phylo.PAML``. Alignments are passed as a mapping from taxon name
to aligned sequence; trees are passed and returned as Newick text. Biopython is imported only when a
function runs. Functions that use randomness take a ``seed``; the seed is applied to Python's global
random state for the duration of the call and the previous state is restored afterwards.

``run_paml`` is a facade over the PAML programs ``codeml``, ``baseml`` and ``yn00`` (installed by the
user; ``conda install -c bioconda paml``). It writes the control file with Biopython's option table,
runs the program with a fixed argument list and no shell, and returns the output files.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import io
import random
import re
import shutil
import warnings
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_MAX_TAXA = 5_000
_MAX_ALIGNMENT_LENGTH = 5_000_000
_MAX_TREE_TEXT = 5_000_000
_MAX_TREES = 10_000
_NAME = re.compile(r"[A-Za-z0-9_.\-]{1,100}")
_PARTITION_NAME = re.compile(r"[A-Za-z0-9_\-]{1,40}")
_SEQUENCE = re.compile(r"[A-Za-z?\-*.]+")
_DATATYPES = ("dna", "rna", "nucleotide")
_DISTANCE_MODELS = ("identity", "blastn", "megablast", "trans", "hoxd70")
_DN_DS_METHODS = ("NG86", "LWL85", "YN00", "ML")
_TREE_METHODS = ("UPGMA", "NJ")
_PAML_PROGRAMS = ("codeml", "baseml", "yn00")
_SUPPORT_MODES = ("branchlength_to_support", "absolute_to_relative")
_ORDERS = ("preorder", "postorder", "level")
_PAML_HINT = (
    "Install PAML (conda install -c bioconda paml, or http://abacus.gene.ucl.ac.uk/software/paml.html) "
    "and put codeml, baseml and yn00 on PATH, or pass executable."
)


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _bio(module: str) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Bio.codonalign and Bio.Nexus emit import-time warnings
        return require(module, _EXTRA)


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value != value:
        raise ConfigurationError(f"{name} must be a number.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _name_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or any(ch in value for ch in "\x00\n\r"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


@contextlib.contextmanager
def _seeded(seed: int | None) -> Iterator[None]:
    """Seed Python's global random state for one call and restore it afterwards."""

    if seed is None:
        yield
        return
    value = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    state = random.getstate()
    random.seed(value)
    try:
        yield
    finally:
        random.setstate(state)


def _alignment(value: object, name: str = "alignment", *, minimum_rows: int = 1) -> dict[str, str]:
    """Validate a mapping of taxon name to aligned sequence (equal lengths)."""

    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of taxon name to aligned sequence.")
    if not minimum_rows <= len(value) <= _MAX_TAXA:
        raise ConfigurationError(f"{name} must hold {minimum_rows} to {_MAX_TAXA} sequences.")
    cleaned: dict[str, str] = {}
    for taxon, sequence in value.items():
        if not isinstance(taxon, str) or not _NAME.fullmatch(taxon):
            raise ConfigurationError(
                f"{name}: taxon names must use letters, digits, '_', '.' or '-' (got {taxon!r})."
            )
        if not isinstance(sequence, str) or not _SEQUENCE.fullmatch(sequence):
            raise ConfigurationError(f"{name}: sequence for {taxon!r} has invalid characters.")
        if len(sequence) > _MAX_ALIGNMENT_LENGTH:
            raise ConfigurationError(f"{name}: sequences may be at most {_MAX_ALIGNMENT_LENGTH} long.")
        cleaned[taxon] = sequence
    if len({len(text) for text in cleaned.values()}) != 1:
        raise ConfigurationError(f"{name}: all sequences must have the same length.")
    return cleaned


def _nexus_text(alignment: Mapping[str, str], datatype: str) -> str:
    names = list(alignment)
    length = len(alignment[names[0]])
    lines = [
        "#NEXUS",
        "begin data;",
        f"dimensions ntax={len(names)} nchar={length};",
        f"format datatype={datatype} missing=? gap=-;",
        "matrix",
        *(f"{taxon} {alignment[taxon]}" for taxon in names),
        ";",
        "end;",
        "",
    ]
    return "\n".join(lines)


def _nexus(alignment: Mapping[str, str], datatype: str) -> Any:
    nexus = _bio("Bio.Nexus.Nexus")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return nexus.Nexus(_nexus_text(alignment, datatype))
    except nexus.NexusError as exc:
        raise ConfigurationError(f"Biopython could not read the alignment as NEXUS: {exc}") from exc


def _matrix_dict(matrix: Mapping[str, Any], order: Sequence[str] | None = None) -> dict[str, str]:
    keys = list(order) if order is not None else list(matrix)
    return {taxon: str(matrix[taxon]) for taxon in keys if taxon in matrix}


def _indices(values: object, name: str, limit: int) -> list[int]:
    if values is None:
        return []
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of integers.")
    result = []
    for item in values:
        result.append(_integer(item, name, minimum=0, maximum=max(limit - 1, 0)))
    return result


def _taxa_list(values: object, name: str, known: Sequence[str]) -> list[str]:
    if values is None:
        return []
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of taxon names.")
    unknown = [item for item in values if item not in known]
    if unknown:
        raise ConfigurationError(f"{name} names taxa that are not in the alignment: {unknown}.")
    return list(values)


def _tree_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be Newick text.")
    if len(value) > _MAX_TREE_TEXT:
        raise ConfigurationError(f"{name} is longer than {_MAX_TREE_TEXT} characters.")
    return value.strip()


def _phylo_tree(value: object, name: str = "newick") -> Any:
    phylo = _bio("Bio.Phylo")
    text = _tree_text(value, name)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return phylo.read(io.StringIO(text), "newick")
    except Exception as exc:  # Biopython raises several parser exception types
        raise ConfigurationError(f"{name} is not readable Newick text: {exc}") from exc


def _newick(tree: Any) -> str:
    return str(tree.format("newick")).strip()


def _clade_record(clade: Any) -> dict[str, Any]:
    return {
        "name": clade.name,
        "branch_length": clade.branch_length,
        "confidence": clade.confidence,
        "terminal": bool(clade.is_terminal()),
        "terminal_count": len(clade.get_terminals()),
    }


def _exact(name: str) -> str:
    return re.escape(name)


def _unique_clade(tree: Any, name: str, label: str) -> Any:
    matches = list(tree.find_clades(name=_exact(name)))
    if not matches:
        raise ConfigurationError(f"{label} {name!r} is not a clade name in the tree.")
    if len(matches) > 1:
        raise ConfigurationError(f"{label} {name!r} names {len(matches)} clades; use a unique name.")
    return matches[0]


def _jsonable(value: Any) -> Any:
    """Convert Biopython result objects (numpy numbers, tuples, sets) to plain JSON types."""

    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    if isinstance(value, bool) or value is None or isinstance(value, (int, float, str)):
        return value
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except (ValueError, TypeError):
            pass
    return str(value)


def _msa(alignment: Mapping[str, str]) -> Any:
    seq = _bio("Bio.Seq")
    record = _bio("Bio.SeqRecord")
    align = _bio("Bio.Align")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return align.MultipleSeqAlignment(
            [record.SeqRecord(seq.Seq(text), id=taxon) for taxon, text in alignment.items()]
        )


def _full_matrix(lower: Sequence[Sequence[float]], size: int) -> list[list[float]]:
    full = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row):
            full[row][column] = full[column][row] = float(lower[row][column])
    return full


# --------------------------------------------------------------------------------------------
# NEXUS character-matrix operations
# --------------------------------------------------------------------------------------------


def nexus_terminal_gaps_to_missing(
    alignment: Mapping[str, str],
    *,
    missing: str = "?",
    skip_n: bool = True,
    datatype: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the leading and trailing gaps of every sequence by the missing-data symbol.

    This is ``Bio.Nexus.Nexus.Nexus.terminal_gap_to_missing``: runs such as ``??---ACG--?`` become
    ``?????ACG???``; internal gaps stay gaps. With ``skip_n=False`` terminal ``N`` runs are also
    treated as terminal gaps. ``alignment`` maps taxon name to aligned sequence (gap ``-``,
    missing ``?``).
    """

    rows = _alignment(alignment)
    symbol = _name_text(missing, "missing")
    if len(symbol) != 1:
        raise ConfigurationError("missing must be a single character.")
    nexus = _nexus(rows, _choice(datatype, _DATATYPES, "datatype"))
    nexus.terminal_gap_to_missing(missing=symbol, skip_n=_flag(skip_n, "skip_n"))
    return {
        "status": "completed",
        "operation": "nexus_terminal_gaps_to_missing",
        "taxa": list(nexus.taxlabels),
        "alignment": _matrix_dict(nexus.matrix, nexus.taxlabels),
        "provenance": _provenance(parser="Bio.Nexus.Nexus.Nexus.terminal_gap_to_missing"),
    }


def nexus_combine_matrices(
    matrices: Mapping[str, Mapping[str, str]],
    *,
    datatype: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate several character matrices into one supermatrix (``Bio.Nexus.Nexus.combine``).

    ``matrices`` maps a data-set name to an alignment (taxon name to aligned sequence). Taxa that are
    absent from a data set receive the missing symbol ``?`` for its columns. The result lists the
    combined alignment, one character set per data set (0-based column indices) and the partition.
    """

    if not isinstance(matrices, Mapping) or not matrices:
        raise ConfigurationError("matrices must be a non-empty mapping of data-set name to alignment.")
    kind = _choice(datatype, _DATATYPES, "datatype")
    nexus_module = _bio("Bio.Nexus.Nexus")
    pairs = []
    for label, alignment in matrices.items():
        if not isinstance(label, str) or not _PARTITION_NAME.fullmatch(label):
            raise ConfigurationError(f"data-set name {label!r} must use letters, digits, '_' or '-'.")
        pairs.append((label, _nexus(_alignment(alignment, f"matrices[{label!r}]"), kind)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        combined = nexus_module.combine(pairs)
    return {
        "status": "completed",
        "operation": "nexus_combine_matrices",
        "taxa": list(combined.taxlabels),
        "character_count": int(combined.nchar),
        "alignment": _matrix_dict(combined.matrix, combined.taxlabels),
        "charsets": {key: [int(i) for i in value] for key, value in combined.charsets.items()},
        "partition": {
            key: [int(i) for i in value]
            for key, value in combined.charpartitions.get("combined", {}).items()
        },
        "provenance": _provenance(parser="Bio.Nexus.Nexus.combine"),
    }


def nexus_write_partition_files(
    alignment: Mapping[str, str],
    charpartition: Mapping[str, Sequence[int]],
    output_dir: str | Path,
    *,
    file_stem: str = "alignment",
    exclude: Sequence[int] | None = None,
    delete: Sequence[str] | None = None,
    interleave: bool = False,
    datatype: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one NEXUS data file per character partition (``write_nexus_data_partitions``).

    ``charpartition`` maps a partition name to its 0-based column indices (for example codon positions).
    Each file ``<file_stem>_<partition>.data`` holds only that partition's columns; ``exclude`` removes
    columns from all files and ``delete`` removes taxa. The output directory must not exist yet.
    """

    rows = _alignment(alignment)
    length = len(next(iter(rows.values())))
    if not isinstance(charpartition, Mapping) or not charpartition:
        raise ConfigurationError("charpartition must be a non-empty mapping of name to column indices.")
    partition: dict[str, list[int]] = {}
    for label, columns in charpartition.items():
        if not isinstance(label, str) or not _PARTITION_NAME.fullmatch(label):
            raise ConfigurationError(f"partition name {label!r} must use letters, digits, '_' or '-'.")
        partition[label] = _indices(columns, f"charpartition[{label!r}]", length)
    stem = _name_text(file_stem, "file_stem")
    if not _PARTITION_NAME.fullmatch(stem):
        raise ConfigurationError("file_stem must use letters, digits, '_' or '-'.")
    nexus = _nexus(rows, _choice(datatype, _DATATYPES, "datatype"))
    excluded = _indices(exclude, "exclude", length)
    deleted = _taxa_list(delete, "delete", nexus.taxlabels)
    directory = cli.new_directory(output_dir, "output_dir")
    target = directory / f"{stem}.nex"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        written = nexus.write_nexus_data_partitions(
            filename=str(target),
            charpartition=partition,
            exclude=excluded,
            delete=deleted,
            interleave=_flag(interleave, "interleave"),
        )
    files = {label: str(Path(path)) for label, path in written.items()}
    return {
        "status": "completed",
        "operation": "nexus_write_partition_files",
        "files": files,
        "partition_sizes": {label: len(set(columns) - set(excluded)) for label, columns in partition.items()},
        "provenance": _provenance(parser="Bio.Nexus.Nexus.Nexus.write_nexus_data_partitions"),
    }


def nexus_weighted_stepmatrix(
    alignment: Mapping[str, str],
    *,
    name: str = "your_name_here",
    exclude: Sequence[int] | None = None,
    delete: Sequence[str] | None = None,
    datatype: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a weighted-parsimony step matrix from an alignment (Wheeler 1990).

    This is ``Bio.Nexus.Nexus.Nexus.weighted_stepmatrix``: the observed co-occurrence of character
    states within sites is normalised and log-transformed. The result carries the NEXUS ``usertype``
    text and the matrix rows (``inf`` marks pairs that never co-occur).
    """

    rows = _alignment(alignment, minimum_rows=2)
    length = len(next(iter(rows.values())))
    nexus = _nexus(rows, _choice(datatype, _DATATYPES, "datatype"))
    label = _name_text(name, "name")
    if not re.fullmatch(r"[A-Za-z0-9_]+", label):
        raise ConfigurationError("name must use letters, digits and '_'.")
    text = nexus.weighted_stepmatrix(
        name=label,
        exclude=_indices(exclude, "exclude", length),
        delete=_taxa_list(delete, "delete", nexus.taxlabels),
    )
    symbols: list[str] = []
    matrix: dict[str, dict[str, float | None]] = {}
    for line in text.splitlines()[2:]:
        match = re.match(r"\[(.)\]\s+(.*)", line)
        if not match:
            continue
        symbols.append(match.group(1))
        matrix[match.group(1)] = {}
    header = text.splitlines()[1].split()
    for line in text.splitlines()[2:]:
        match = re.match(r"\[(.)\]\s+(.*)", line)
        if not match:
            continue
        values = match.group(2).split()
        for column, value in zip(header, values):
            if value == ".":
                matrix[match.group(1)][column] = 0.0
            elif value.startswith("inf"):
                matrix[match.group(1)][column] = None
            else:
                matrix[match.group(1)][column] = float(value)
    return {
        "status": "completed",
        "operation": "nexus_weighted_stepmatrix",
        "symbols": symbols,
        "stepmatrix": text,
        "weights": matrix,
        "note": "weights: null marks 'inf' (the two states never occur together in a site); 0.0 is the diagonal",
        "provenance": _provenance(parser="Bio.Nexus.Nexus.Nexus.weighted_stepmatrix"),
    }


def nexus_bootstrap_alignment(
    alignment: Mapping[str, str],
    *,
    seed: int | None = None,
    exclude: Sequence[int] | None = None,
    delete: Sequence[str] | None = None,
    datatype: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Resample the alignment columns with replacement (``Bio.Nexus.Nexus.Nexus.bootstrap``).

    Returns a bootstrapped alignment of the same length. ``exclude`` removes columns and ``delete``
    removes taxa before resampling. With ``seed`` the draw is reproducible.
    """

    rows = _alignment(alignment)
    length = len(next(iter(rows.values())))
    nexus = _nexus(rows, _choice(datatype, _DATATYPES, "datatype"))
    with _seeded(seed), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sample = nexus.bootstrap(
            delete=_taxa_list(delete, "delete", nexus.taxlabels),
            exclude=_indices(exclude, "exclude", length),
        )
    return {
        "status": "completed",
        "operation": "nexus_bootstrap_alignment",
        "seed": seed,
        "alignment": _matrix_dict(sample, nexus.taxlabels),
        "provenance": _provenance(
            parser="Bio.Nexus.Nexus.Nexus.bootstrap",
            randomness="the global random state is seeded for the call and restored afterwards",
        ),
    }


# --------------------------------------------------------------------------------------------
# Distances, parsimony and consensus
# --------------------------------------------------------------------------------------------


def biopython_alignment_distance_matrix(
    alignment: Mapping[str, str],
    *,
    model: str = "identity",
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise distance matrix of a DNA alignment with ``Bio.Phylo.TreeConstruction.DistanceCalculator``.

    ``model`` is ``identity`` (fraction of differing aligned positions) or one of the nucleotide
    substitution-matrix models ``blastn``, ``megablast``, ``trans`` and ``hoxd70`` (distance derived from
    the matrix scores). Gap columns are scored as Biopython does. Returns taxon ids and the full
    symmetric matrix.
    """

    rows = _alignment(alignment, minimum_rows=2)
    selected = _choice(model, _DISTANCE_MODELS, "model")
    construction = _bio("Bio.Phylo.TreeConstruction")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            matrix = construction.DistanceCalculator(selected).get_distance(_msa(rows))
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Biopython could not compute the distances: {exc}") from exc
    ids = list(matrix.names)
    return {
        "status": "completed",
        "operation": "biopython_alignment_distance_matrix",
        "model": selected,
        "ids": ids,
        "matrix": _full_matrix(matrix.matrix, len(ids)),
        "provenance": _provenance(parser=f"Bio.Phylo.TreeConstruction.DistanceCalculator({selected!r})"),
    }


def parsimony_tree_search(
    alignment: Mapping[str, str],
    *,
    starting_tree: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Maximum-parsimony tree search by nearest-neighbour interchange (NNI).

    This is ``ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()))`` of Biopython. Without
    ``starting_tree`` the search starts from a neighbour-joining tree built on identity distances.
    ``starting_tree`` must carry branch lengths. The result holds the Newick tree and its parsimony score
    (Fitch, unit costs).
    """

    rows = _alignment(alignment, minimum_rows=3)
    construction = _bio("Bio.Phylo.TreeConstruction")
    start = _phylo_tree(starting_tree, "starting_tree") if starting_tree is not None else None
    if start is not None:
        tips = {clade.name for clade in start.get_terminals()}
        if tips != set(rows):
            raise ConfigurationError("starting_tree tips must equal the alignment taxa.")
        if any(clade.branch_length is None for clade in start.find_clades() if clade is not start.root):
            raise ConfigurationError("starting_tree needs a branch length on every branch (NNI roots at the midpoint).")
    msa = _msa(rows)
    scorer = construction.ParsimonyScorer()
    searcher = construction.NNITreeSearcher(scorer)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tree = construction.ParsimonyTreeConstructor(searcher, start).build_tree(msa)
            score = scorer.get_score(tree, msa)
    except (ValueError, TypeError, KeyError) as exc:
        raise ConfigurationError(f"Biopython could not search for a parsimony tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parsimony_tree_search",
        "newick": _newick(tree),
        "parsimony_score": float(score),
        "starting_tree": "supplied" if start is not None else "neighbour joining (identity distance)",
        "provenance": _provenance(
            parser="Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor / NNITreeSearcher"
        ),
    }


def strict_consensus_tree(newick_trees: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Strict consensus of several trees: the clades present in every input tree.

    This is the result of ``Bio.Phylo.Consensus.strict_consensus``. Biopython compares clades by the
    position of the tips, which is only correct when every tree lists its tips in the same order; this
    function matches clades by tip name instead, so the order of the tips does not matter, and returns
    the same tree as Biopython when the order is identical. All trees must contain the same taxa.
    Branch lengths of the result are not meaningful.
    """

    if isinstance(newick_trees, (str, bytes)) or not isinstance(newick_trees, Sequence):
        raise ConfigurationError("newick_trees must be a list of Newick strings.")
    if not 2 <= len(newick_trees) <= _MAX_TREES:
        raise ConfigurationError(f"newick_trees must hold 2 to {_MAX_TREES} trees.")
    trees = [_phylo_tree(text, f"newick_trees[{i}]") for i, text in enumerate(newick_trees)]
    names = [[clade.name for clade in tree.get_terminals()] for tree in trees]
    if any(name is None or name == "" for row in names for name in row):
        raise ConfigurationError("every tip must have a name.")
    taxa = set(names[0])
    if any(set(row) != taxa or len(set(row)) != len(row) for row in names):
        raise ConfigurationError("all trees must contain the same, unique taxa.")
    shared: set[frozenset[str]] | None = None
    for tree in trees:
        clades = {
            frozenset(tip.name for tip in clade.get_terminals())
            for clade in tree.find_clades(terminal=False)
        }
        shared = clades if shared is None else shared & clades
    base = _bio("Bio.Phylo.BaseTree")
    tips = {name: base.Clade(name=name) for name in names[0]}
    built: list[tuple[frozenset[str], Any]] = []
    for group in sorted((g for g in (shared or set()) if 1 < len(g) < len(taxa)), key=len):
        node = base.Clade()
        covered: set[str] = set()
        for subset, child in built:
            if subset <= group and not subset & covered:
                node.clades.append(child)
                covered |= subset
        node.clades.extend(tips[name] for name in names[0] if name in group and name not in covered)
        built = [(s, c) for s, c in built if not s <= group]
        built.append((group, node))
    root = base.Clade()
    covered = set()
    for subset, child in built:
        root.clades.append(child)
        covered |= subset
    root.clades.extend(tips[name] for name in names[0] if name not in covered)
    result = base.Tree(root=root)
    return {
        "status": "completed",
        "operation": "strict_consensus_tree",
        "tree_count": len(trees),
        "shared_clade_count": len([g for g in (shared or set()) if 1 < len(g) < len(taxa)]),
        "newick": _newick(result),
        "provenance": _provenance(
            parser="Bio.Phylo.Consensus.strict_consensus semantics (clades matched by tip name)"
        ),
    }


# --------------------------------------------------------------------------------------------
# Tree editing and queries
# --------------------------------------------------------------------------------------------


def collapse_tree_genera(
    newick: str,
    *,
    space_equals_underscore: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collapse every subtree whose taxa share the same genus (``Bio.Nexus.Trees.Tree.collapse_genera``).

    The genus is the first word of the taxon name (``Genus_species``). A collapsed subtree is replaced
    by one tip called ``<Genus> <collapsed>``. With ``space_equals_underscore`` a space also separates
    genus and species.
    """

    trees = _bio("Bio.Nexus.Trees")
    text = _tree_text(newick, "newick")
    try:
        tree = trees.Tree(text if text.endswith(";") else text + ";")
    except Exception as exc:
        raise ConfigurationError(f"newick is not readable by Bio.Nexus.Trees: {exc}") from exc
    before = len(tree.get_terminals())
    tree.collapse_genera(space_equals_underscore=_flag(space_equals_underscore, "space_equals_underscore"))
    return {
        "status": "completed",
        "operation": "collapse_tree_genera",
        "newick": tree.to_string(plain=False, plain_newick=True) + ";",
        "terminal_count_before": before,
        "terminal_count_after": len(tree.get_terminals()),
        "provenance": _provenance(parser="Bio.Nexus.Trees.Tree.collapse_genera"),
    }


def split_tree_clade(
    newick: str,
    target: str,
    *,
    count: int = 2,
    branch_length: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a named clade into ``count`` new descendants (``Bio.Phylo.BaseTree.TreeMixin.split``).

    The new child clades are named ``<target>0``, ``<target>1``, ... and receive ``branch_length``.
    A tip becomes an internal node with the new tips below it.
    """

    tree = _phylo_tree(newick)
    clade = _unique_clade(tree, _name_text(target, "target"), "target")
    number = _integer(count, "count", minimum=1, maximum=1000)
    length = _number(branch_length, "branch_length", minimum=0.0)
    before = {id(item) for item in tree.find_clades()}
    clade.split(n=number, branch_length=length)
    created = [item.name for item in tree.find_clades() if id(item) not in before]
    return {
        "status": "completed",
        "operation": "split_tree_clade",
        "new_clades": created,
        "newick": _newick(tree),
        "provenance": _provenance(parser="Bio.Phylo.BaseTree.TreeMixin.split"),
    }


def random_tree(
    taxa: Sequence[str],
    *,
    branch_length: float = 1.0,
    branch_stdev: float | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a random bifurcating tree on the given taxa (``Bio.Phylo.BaseTree.Tree.randomized``).

    Terminal branches get ``branch_length``; with ``branch_stdev`` the lengths are drawn from a normal
    distribution (negative draws become 0). The algorithm repeatedly splits a random tip, so it is not a
    coalescent or birth-death model. Use ``seed`` for reproducible output.
    """

    if isinstance(taxa, (str, bytes)) or not isinstance(taxa, Sequence):
        raise ConfigurationError("taxa must be a list of names.")
    if not 2 <= len(taxa) <= 10_000 or len(set(taxa)) != len(taxa):
        raise ConfigurationError("taxa must be 2 to 10000 distinct names.")
    names = [_name_text(item, "taxa") for item in taxa]
    length = _number(branch_length, "branch_length", minimum=0.0)
    stdev = None if branch_stdev is None else _number(branch_stdev, "branch_stdev", minimum=0.0)
    base = _bio("Bio.Phylo.BaseTree")
    with _seeded(seed):
        tree = base.Tree.randomized(names, branch_length=length, branch_stdev=stdev)
    return {
        "status": "completed",
        "operation": "random_tree",
        "seed": seed,
        "newick": _newick(tree),
        "provenance": _provenance(
            parser="Bio.Phylo.BaseTree.Tree.randomized",
            randomness="the global random state is seeded for the call and restored afterwards",
        ),
    }


def tree_is_bifurcating(newick: str, *, verbose: bool = True) -> dict[str, Any]:
    """Test whether a tree is strictly bifurcating (``Bio.Phylo.BaseTree.TreeMixin.is_bifurcating``).

    A tree is bifurcating when every internal node has two children; a root with three children is
    accepted (unrooted form), as in Biopython.
    """

    tree = _phylo_tree(newick)
    internal = [clade for clade in tree.find_clades() if not clade.is_terminal()]
    return {
        "status": "completed",
        "operation": "tree_is_bifurcating",
        "is_bifurcating": bool(tree.is_bifurcating()),
        "terminal_count": len(tree.get_terminals()),
        "internal_child_counts": sorted({len(clade.clades) for clade in internal}),
        "provenance": _provenance(parser="Bio.Phylo.BaseTree.TreeMixin.is_bifurcating"),
    }


def tree_is_parent_of(
    newick: str, ancestor: str, descendant: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Test whether the clade ``ancestor`` contains the clade ``descendant`` (``is_parent_of``).

    Clades are selected by exact name; both names must be unique in the tree. As in Biopython, a clade
    counts as a descendant of itself.
    """

    tree = _phylo_tree(newick)
    parent = _unique_clade(tree, _name_text(ancestor, "ancestor"), "ancestor")
    child = _unique_clade(tree, _name_text(descendant, "descendant"), "descendant")
    return {
        "status": "completed",
        "operation": "tree_is_parent_of",
        "ancestor": ancestor,
        "descendant": descendant,
        "is_parent_of": bool(parent.is_parent_of(child)),
        "provenance": _provenance(parser="Bio.Phylo.BaseTree.TreeMixin.is_parent_of"),
    }


def find_tree_clades(
    newick: str,
    *,
    name: str | None = None,
    terminal: bool | None = None,
    branch_length: float | None = None,
    order: str = "preorder",
    first_only: bool = False,
    max_results: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find clades by attributes (``find_clades`` / ``find_any`` / ``find_elements`` of ``Bio.Phylo``).

    ``name`` is a regular expression that must match the whole clade name (as in Biopython);
    ``terminal`` selects tips (True) or internal nodes (False); ``branch_length`` matches by equality (within 1e-9).
    ``first_only=True`` returns only the first match (``find_any``). ``order`` is preorder, postorder
    or level.
    """

    tree = _phylo_tree(newick)
    criteria: dict[str, Any] = {}
    if name is not None:
        criteria["name"] = _name_text(name, "name")
    wanted_length = None if branch_length is None else _number(branch_length, "branch_length")
    if terminal is not None:
        criteria["terminal"] = _flag(terminal, "terminal")
    limit = _integer(max_results, "max_results", minimum=1, maximum=1_000_000)
    selected = _choice(order, _ORDERS, "order")
    single = _flag(first_only, "first_only")
    try:
        found = []
        for clade in tree.find_clades(order=selected, **criteria):
            if wanted_length is not None and (
                clade.branch_length is None or abs(clade.branch_length - wanted_length) > 1e-9
            ):
                continue  # Biopython cannot match floats; the branch length is compared here
            found.append(clade)
            if single or len(found) >= limit:
                break
    except re.error as exc:
        raise ConfigurationError(f"name is not a valid regular expression: {exc}") from exc
    return {
        "status": "completed",
        "operation": "find_tree_clades",
        "match_count": len(found),
        "clades": [_clade_record(clade) for clade in found],
        "provenance": _provenance(parser="Bio.Phylo.BaseTree.TreeMixin.find_clades / find_any"),
    }


def _nexus_tree(value: object, name: str) -> Any:
    trees = _bio("Bio.Nexus.Trees")
    text = _tree_text(value, name)
    try:
        return trees.Tree(text if text.endswith(";") else text + ";")
    except Exception as exc:
        raise ConfigurationError(f"{name} is not readable by Bio.Nexus.Trees: {exc}") from exc


def trees_are_identical(newick_a: str, newick_b: str, *, verbose: bool = True) -> dict[str, Any]:
    """Test whether two trees have identical topology and taxa (``Bio.Nexus.Trees.Tree.is_identical``).

    The comparison ignores branch lengths and the order of children.
    """

    first = _nexus_tree(newick_a, "newick_a")
    second = _nexus_tree(newick_b, "newick_b")
    return {
        "status": "completed",
        "operation": "trees_are_identical",
        "identical": bool(first.is_identical(second)),
        "provenance": _provenance(parser="Bio.Nexus.Trees.Tree.is_identical"),
    }


def trees_are_compatible(
    newick_a: str,
    newick_b: str,
    *,
    threshold: float,
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check two trees for conflicting supported clades (``Bio.Nexus.Trees.Tree.is_compatible``).

    Only internal branches with support at least ``threshold`` are compared (the trees must carry support
    values written as ``)support:length`` after the internal node, for example ``(A,B)0.9:1``;
    a single number after ``)`` is read as a branch length). ``strict=True`` requires the same taxa. The trees are compatible when no conflict is listed.
    """

    first = _nexus_tree(newick_a, "newick_a")
    second = _nexus_tree(newick_b, "newick_b")
    cutoff = _number(threshold, "threshold")
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            conflicts = first.is_compatible(second, cutoff, strict=_flag(strict, "strict"))
        except Exception as exc:
            raise ConfigurationError(f"Biopython could not compare the trees: {exc}") from exc
    return {
        "status": "completed",
        "operation": "trees_are_compatible",
        "compatible": not conflicts,
        "conflict_count": len(conflicts),
        "conflicts": [
            {
                "clade_a": sorted(item[0]),
                "support_a": item[1],
                "clade_b": sorted(item[2]),
                "support_b": item[3],
                "shared": sorted(item[4]),
            }
            for item in conflicts
        ],
        "provenance": _provenance(parser="Bio.Nexus.Trees.Tree.is_compatible"),
    }


def convert_tree_support(
    newick: str,
    *,
    mode: str = "absolute_to_relative",
    replicates: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert clade support values (``Bio.Nexus.Trees.Tree``).

    ``mode="absolute_to_relative"`` divides clade counts (as written by PHYLIP consense) by
    ``replicates`` (``convert_absolute_support``). ``mode="branchlength_to_support"`` moves values that
    were read as branch lengths (for example PAUP bootstrap output) to the support field and sets the
    branch lengths to 0 (``branchlength2support``). Returns Newick with ``support:length`` labels.
    """

    selected = _choice(mode, _SUPPORT_MODES, "mode")
    tree = _nexus_tree(newick, "newick")
    if selected == "absolute_to_relative":
        if replicates is None:
            raise ConfigurationError("replicates is required for mode 'absolute_to_relative'.")
        tree.convert_absolute_support(_integer(replicates, "replicates", minimum=1, maximum=10**9))
    else:
        tree.branchlength2support()
    nodes = [
        {
            "taxon": tree.node(node).data.taxon,
            "support": tree.node(node).data.support,
            "branch_length": tree.node(node).data.branchlength,
        }
        for node in tree.all_ids()
    ]
    return {
        "status": "completed",
        "operation": "convert_tree_support",
        "mode": selected,
        "newick": tree.to_string(plain=False, plain_newick=True) + ";",
        "nodes": nodes,
        "provenance": _provenance(
            parser="Bio.Nexus.Trees.Tree.convert_absolute_support"
            if selected == "absolute_to_relative"
            else "Bio.Nexus.Trees.Tree.branchlength2support"
        ),
    }


# --------------------------------------------------------------------------------------------
# dN/dS
# --------------------------------------------------------------------------------------------


def _codon_alignment_rows(value: object) -> dict[str, str]:
    rows = _alignment(value, "codon_alignment", minimum_rows=2)
    for taxon, text in rows.items():
        if len(text) % 3 or len(text) < 3:
            raise ConfigurationError(f"codon_alignment: {taxon!r} is not a whole number of codons.")
        if set(text.upper()) - set("ACGT"):
            raise ConfigurationError(
                f"codon_alignment: {taxon!r} must hold only A, C, G and T (no gaps or ambiguity)."
            )
    return {taxon: text.upper() for taxon, text in rows.items()}


def codon_dn_ds_estimates(
    codon_alignment: Mapping[str, str],
    *,
    method: str = "YN00",
    genetic_code: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise dN and dS of a codon alignment (``Bio.Align.analysis.calculate_dn_ds_matrix``).

    ``method`` is ``NG86`` (Nei-Gojobori), ``LWL85`` (Li-Wu-Luo), ``YN00`` (Yang-Nielsen) or ``ML``
    (maximum likelihood, Goldman-Yang). Rows of ``codon_alignment`` are in frame and ungapped. A value
    of null marks a saturated pair (Biopython reports -1). Returns full symmetric matrices.
    """

    rows = _codon_alignment_rows(codon_alignment)
    selected = _choice(method, _DN_DS_METHODS, "method")
    code = _integer(genetic_code, "genetic_code", minimum=1, maximum=33)
    tables = _bio("Bio.Data.CodonTable")
    analysis = _bio("Bio.Align.analysis")
    align = _bio("Bio.Align")
    seq = _bio("Bio.Seq")
    record = _bio("Bio.SeqRecord")
    np = require("numpy", "bio-analysis")
    try:
        table = tables.generic_by_id[code]
    except KeyError as exc:
        raise ConfigurationError(f"genetic_code {code} is not a NCBI translation table.") from exc
    length = len(next(iter(rows.values())))
    alignment = align.Alignment(
        [record.SeqRecord(seq.Seq(text), id=taxon) for taxon, text in rows.items()],
        np.array([[0, length] for _ in rows]),
    )
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dn, ds = analysis.calculate_dn_ds_matrix(alignment, method=selected, codon_table=table)
    except (ZeroDivisionError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"{selected} could not estimate dN/dS for this alignment: {exc}") from exc
    ids = list(rows)

    def full(lower: Sequence[Sequence[float]]) -> list[list[float | None]]:
        result: list[list[float | None]] = [[0.0] * len(ids) for _ in ids]
        for i in range(len(ids)):
            for j in range(i):
                value = lower[i][j]
                cell = None if value is None or value < 0 else float(value)
                result[i][j] = result[j][i] = cell
        return result

    return {
        "status": "completed",
        "operation": "codon_dn_ds_estimates",
        "method": selected,
        "genetic_code": code,
        "ids": ids,
        "dn": full(dn.matrix),
        "ds": full(ds.matrix),
        "provenance": _provenance(
            parser=f"Bio.Align.analysis.calculate_dn_ds_matrix ({selected})",
            saturation="null marks a negative value, which Biopython reports when changes saturate",
        ),
    }


def codon_dn_ds_trees(
    codon_alignment: Mapping[str, str],
    *,
    dn_ds_method: str = "NG86",
    tree_method: str = "UPGMA",
    genetic_code: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a dN tree and a dS tree from a codon alignment (``CodonAlignment.get_dn_ds_tree``).

    The pairwise dN and dS matrices are estimated with ``dn_ds_method`` (NG86, LWL85, YN00 or ML) and
    turned into trees with ``tree_method`` (UPGMA or NJ). Returns both trees as Newick.
    """

    rows = _codon_alignment_rows(codon_alignment)
    method = _choice(dn_ds_method, _DN_DS_METHODS, "dn_ds_method")
    builder = _choice(tree_method, _TREE_METHODS, "tree_method")
    code = _integer(genetic_code, "genetic_code", minimum=1, maximum=33)
    tables = _bio("Bio.Data.CodonTable")
    codonalign = _bio("Bio.codonalign.codonalignment")
    codonseq = _bio("Bio.codonalign.codonseq")
    record = _bio("Bio.SeqRecord")
    try:
        table = tables.generic_by_id[code]
    except KeyError as exc:
        raise ConfigurationError(f"genetic_code {code} is not a NCBI translation table.") from exc
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            alignment = codonalign.CodonAlignment(
                [record.SeqRecord(codonseq.CodonSeq(text), id=taxon) for taxon, text in rows.items()]
            )
            dn_tree, ds_tree = alignment.get_dn_ds_tree(
                dn_ds_method=method, tree_method=builder, codon_table=table
            )
    except (ZeroDivisionError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"{method} could not estimate dN/dS for this alignment: {exc}") from exc
    return {
        "status": "completed",
        "operation": "codon_dn_ds_trees",
        "dn_ds_method": method,
        "tree_method": builder,
        "genetic_code": code,
        "dn_tree_newick": _newick(dn_tree),
        "ds_tree_newick": _newick(ds_tree),
        "provenance": _provenance(parser="Bio.codonalign.codonalignment.CodonAlignment.get_dn_ds_tree"),
    }


# --------------------------------------------------------------------------------------------
# PAML
# --------------------------------------------------------------------------------------------


def _paml_class(program: str) -> Any:
    module = _bio(f"Bio.Phylo.PAML.{program}")
    return getattr(module, {"codeml": "Codeml", "baseml": "Baseml", "yn00": "Yn00"}[program])


def _paml_options(program: str, options: Mapping[str, Any] | None) -> list[tuple[str, Any]]:
    probe = _paml_class(program)()
    known = dict(probe.get_all_options())
    chosen: list[tuple[str, Any]] = []
    for key, value in (options or {}).items():
        if key not in known:
            raise ConfigurationError(f"{key!r} is not a {program} control-file option.")
        if key == "NSsites":
            if isinstance(value, int) and not isinstance(value, bool):
                value = [value]
            if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
                raise ConfigurationError("NSsites must be an integer or a list of integers.")
            value = [_integer(item, "NSsites", minimum=0, maximum=13) for item in value]
        elif isinstance(value, bool) or not isinstance(value, (int, float, str)):
            raise ConfigurationError(f"option {key} must be a number or text.")
        elif isinstance(value, str) and (not value.strip() or "\n" in value or "\r" in value):
            raise ConfigurationError(f"option {key} must be single-line text.")
        chosen.append((key, value))
    return chosen


def run_paml(
    program: str,
    alignment_path: str | Path,
    working_dir: str | Path,
    *,
    tree_path: str | Path | None = None,
    options: Mapping[str, Any] | None = None,
    executable: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run PAML ``codeml``, ``baseml`` or ``yn00`` on a PHYLIP/PAML sequence file.

    ``program`` selects the program (``Bio.Phylo.PAML.codeml.Codeml.run``, ``baseml.Baseml.run``,
    ``yn00.Yn00.run``). ``options`` are control-file options by their PAML names (for example
    ``{"runmode": -2, "seqtype": 1, "CodonFreq": 2}``); an unknown name is refused. The alignment
    (and the optional tree; a bare Newick tree gets the ``<taxa> 1`` header line that PAML expects) are
    copied into ``working_dir`` (a new directory), the control file is written there, and the program runs with the control file as its only argument. The results file
    ``results.out`` is parsed with ``read_paml_results`` when the run succeeds.
    """

    chosen = _choice(program, _PAML_PROGRAMS, "program")
    source = cli.existing_file(alignment_path, "alignment_path")
    tree_source = cli.existing_file(tree_path, "tree_path") if tree_path is not None else None
    if tree_source is not None and chosen == "yn00":
        raise ConfigurationError("yn00 does not use a tree file.")
    settings = _paml_options(chosen, options)
    limit = _number(timeout_seconds, "timeout_seconds", minimum=1.0)
    directory = cli.new_directory(working_dir, "working_dir")
    seq_name, tree_name, out_name, ctl_name = "seq.phy", "tree.nwk", "results.out", f"{chosen}.ctl"
    shutil.copyfile(source, directory / seq_name)
    if tree_source is not None:
        tree_text = tree_source.read_text(encoding="utf-8", errors="replace")
        header = source.read_text(encoding="utf-8", errors="replace").split(None, 1)
        if tree_text.lstrip().startswith("(") and header and header[0].isdigit():
            # PAML ends with an "end of tree file" error on a bare Newick tree: give it the "<taxa> <trees>" header
            tree_text = f"{header[0]} 1\n{tree_text.lstrip()}"
        (directory / tree_name).write_text(tree_text if tree_text.endswith("\n") else tree_text + "\n", encoding="utf-8")
    lines = [f"seqfile = {seq_name}", f"outfile = {out_name}"]
    if tree_source is not None:
        lines.append(f"treefile = {tree_name}")
    defaults = dict(_paml_class(chosen)().get_all_options())
    merged = {key: value for key, value in defaults.items() if value is not None}
    merged.update(dict(settings))
    for key, value in merged.items():
        text = " ".join(str(item) for item in value) if key == "NSsites" else str(value)
        lines.append(f"{key} = {text}")
    (directory / ctl_name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    outcome = cli.run(
        chosen,
        [ctl_name],
        cwd=directory,
        expected=[directory / out_name],
        explicit=executable,
        hint=_PAML_HINT,
        timeout_seconds=limit,
    )
    outcome.update(
        {
            "operation": "run_paml",
            "program": chosen,
            "working_dir": str(directory),
            "control_file": str(directory / ctl_name),
            "results_file": str(directory / out_name),
            "files": sorted(path.name for path in directory.iterdir()),
        }
    )
    if outcome["status"] == "completed":
        head = (directory / out_name).read_text(encoding="utf-8", errors="replace")[:2000]
        found = re.search(r"paml version ([0-9][0-9A-Za-z.\-]*)", head)
        outcome["version"] = found.group(1) if found else outcome.get("version")
        try:
            outcome["results"] = read_paml_results(chosen, directory / out_name, verbose=False)["results"]
        except ConfigurationError as exc:
            outcome["results"] = None
            outcome["parse_error"] = str(exc)
    outcome["provenance"] = _provenance(
        parser=f"Bio.Phylo.PAML.{chosen} option table; program run by DNAKit",
    )
    return outcome


def read_paml_results(
    program: str, results_file: str | Path, *, verbose: bool = True
) -> dict[str, Any]:
    """Parse a PAML results file of ``codeml``, ``baseml`` or ``yn00`` (``Bio.Phylo.PAML.<program>.read``).

    Returns the parsed content as nested dictionaries: ``codeml`` gives site-model and pairwise
    results, ``baseml`` the model parameters and tree lengths, ``yn00`` the Nei-Gojobori, Yang-Nielsen
    and LWL85 dN/dS tables.
    """

    chosen = _choice(program, _PAML_PROGRAMS, "program")
    path = input_path(results_file)
    module = _bio(f"Bio.Phylo.PAML.{chosen}")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parsed = module.read(str(path))
    except (ValueError, KeyError, IndexError, TypeError) as exc:
        raise ConfigurationError(f"{path.name} is not a valid {chosen} results file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_paml_results",
        "program": chosen,
        "results": _jsonable(parsed),
        "provenance": _provenance(parser=f"Bio.Phylo.PAML.{chosen}.read"),
    }


__all__ = [
    "biopython_alignment_distance_matrix",
    "codon_dn_ds_estimates",
    "codon_dn_ds_trees",
    "collapse_tree_genera",
    "convert_tree_support",
    "find_tree_clades",
    "nexus_bootstrap_alignment",
    "nexus_combine_matrices",
    "nexus_terminal_gaps_to_missing",
    "nexus_weighted_stepmatrix",
    "nexus_write_partition_files",
    "parsimony_tree_search",
    "random_tree",
    "read_paml_results",
    "run_paml",
    "split_tree_clade",
    "strict_consensus_tree",
    "tree_is_bifurcating",
    "tree_is_parent_of",
    "trees_are_compatible",
    "trees_are_identical",
]

_register_public_verbose(globals(), tuple(__all__))
