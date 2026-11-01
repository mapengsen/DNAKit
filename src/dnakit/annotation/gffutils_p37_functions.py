"""Feature look-ups, derived features, database edits, GFF writing and gene-model plots through gffutils (pack P37).

Every function reads a local GFF3/GTF file (loaded into an in-memory gffutils database, nothing is
written) or an existing gffutils ``.db`` file. Functions that change a database never touch the input:
they write a new ``.db`` file at ``output_path`` that must not exist yet. Coordinates are 1-based and
inclusive as stored in GFF. gffutils is imported only when a function runs; a missing package raises
``BackendUnavailableError`` with the ``dnakit[bio-gffutils]`` hint.
"""

from __future__ import annotations

import io
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import _DEFAULT_ROWS, _checked, _package, _rows_limit, _text
from dnakit.annotation.gffutils_review_functions import _feature_record, _open_database
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_DATABASE_SUFFIXES = frozenset({".db", ".sqlite", ".sqlite3"})
_IMAGE_SUFFIXES = frozenset({".png", ".svg", ".pdf"})
_MAX_LEVEL = 50


def _gffutils() -> Any:
    return require("gffutils", "bio-gffutils")


def _provenance(**notes: object) -> dict[str, Any]:
    return {"package": _package("gffutils"), **notes}


def _optional_text(value: object, name: str) -> str | None:
    return None if value is None else _text(value, name)


def _level(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_LEVEL:
        raise ConfigurationError(f"level must be None or an integer in [1, {_MAX_LEVEL}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _new_file(value: str | Path, suffixes: frozenset[str] | None, name: str) -> Path:
    path = Path(value).expanduser()
    if path.exists():
        raise ConfigurationError(f"{name} already exists: {path.name}; choose a new file.")
    if suffixes is not None and path.suffix.lower() not in suffixes:
        raise ConfigurationError(f"{name} must end with one of {sorted(suffixes)}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}.")
    return path


def _open(source: str | Path) -> tuple[Any, Any]:
    gffutils = _gffutils()
    return gffutils, _open_database(gffutils, input_path(source))


def _new_database(source: str | Path, output_path: str | Path) -> tuple[Any, Any, Path]:
    """Copy a database, or build one from a GFF/GTF file, into a new ``.db`` file and open it."""

    gffutils = _gffutils()
    path = input_path(source)
    target = _new_file(output_path, _DATABASE_SUFFIXES, "output_path")
    try:
        if path.suffix in _DATABASE_SUFFIXES:
            shutil.copyfile(path, target)
        else:
            gffutils.create_db(str(path), str(target), verbose=False)
        return gffutils, gffutils.FeatureDB(str(target)), target
    except Exception as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"gffutils could not build a database from {path.name}: {exc}") from exc


def _feature(db: Any, feature_id: str) -> Any:
    try:
        return db[_text(feature_id, "feature_id")]
    except Exception as exc:
        raise ConfigurationError(f"feature {feature_id!r} was not found in the annotation.") from exc


def _rows(features: Any, limit: int, what: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for feature in features:
        rows.append(_feature_record(feature) | {"id": feature.id})
        if len(rows) > limit:
            break
    return _checked(rows, limit, what)


def fetch_gffutils_feature(source: str | Path, feature_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Fetch one annotation feature by its ID from a GFF/GTF file or gffutils database (``gffutils-cli fetch``)."""

    _, db = _open(source)
    feature = _feature(db, feature_id)
    return {
        "status": "completed",
        "operation": "fetch_gffutils_feature",
        "feature": _feature_record(feature) | {"id": feature.id},
        "gff_line": str(feature),
        "provenance": _provenance(),
    }


def list_gffutils_features(
    source: str | Path,
    *,
    featuretype: str | None = None,
    seqid: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List annotation features, all of them or one feature type, optionally on one sequence.

    ``featuretype=None`` walks ``FeatureDB.all_features()``; a type name walks ``features_of_type()``.
    """

    limit = _rows_limit(max_rows)
    kind = _optional_text(featuretype, "featuretype")
    chrom = _optional_text(seqid, "seqid")
    _, db = _open(source)
    features = db.features_of_type(kind) if kind else db.all_features()
    if chrom:
        features = (feature for feature in features if feature.seqid == chrom)
    rows = _rows(features, limit, "feature list")
    return {
        "status": "completed",
        "operation": "list_gffutils_features",
        "featuretype": kind,
        "seqid": chrom,
        "count": len(rows),
        "features": rows,
        "provenance": _provenance(),
    }


def group_gffutils_parent_children(
    source: str | Path,
    *,
    parent_featuretype: str = "gene",
    level: int | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group every feature of ``parent_featuretype`` with its child features (``FeatureDB.iter_by_parent_childs``)."""

    limit = _rows_limit(max_rows)
    kind = _text(parent_featuretype, "parent_featuretype")
    _, db = _open(source)
    groups: list[dict[str, Any]] = []
    total = 0
    for group in db.iter_by_parent_childs(featuretype=kind, level=_level(level)):
        parent, children = group[0], group[1:]
        total += 1 + len(children)
        groups.append(
            {
                "parent": _feature_record(parent) | {"id": parent.id},
                "children": [_feature_record(child) | {"id": child.id} for child in children],
            }
        )
        if total > limit:
            break
    _checked([None] * total, limit, "parent-child groups")
    return {
        "status": "completed",
        "operation": "group_gffutils_parent_children",
        "parent_featuretype": kind,
        "group_count": len(groups),
        "groups": groups,
        "provenance": _provenance(),
    }


def fetch_gffutils_parents(
    source: str | Path,
    feature_id: str,
    *,
    level: int | None = None,
    featuretype: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch the parent features of a feature (``FeatureDB.parents``, ``gffutils-cli parents``)."""

    _, db = _open(source)
    feature = _feature(db, feature_id)
    parents = list(db.parents(feature, level=_level(level), featuretype=_optional_text(featuretype, "featuretype")))
    return {
        "status": "completed",
        "operation": "fetch_gffutils_parents",
        "feature_id": feature.id,
        "count": len(parents),
        "parents": [_feature_record(parent) | {"id": parent.id} for parent in parents],
        "provenance": _provenance(),
    }


def get_gffutils_interfeatures(
    source: str | Path,
    *,
    featuretype: str = "gene",
    new_featuretype: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the gaps between consecutive features of one type (``FeatureDB.interfeatures``).

    Features are sorted by sequence and start; gaps are made only between neighbours on the same sequence.
    """

    limit = _rows_limit(max_rows)
    kind = _text(featuretype, "featuretype")
    _, db = _open(source)
    ordered = db.features_of_type(kind, order_by=("seqid", "start"))
    gaps = []
    for gap in db.interfeatures(ordered, new_featuretype=_optional_text(new_featuretype, "new_featuretype")):
        gaps.append(_feature_record(gap))
        if len(gaps) > limit:
            break
    _checked(gaps, limit, "inter-feature gaps")
    return {
        "status": "completed",
        "operation": "get_gffutils_interfeatures",
        "featuretype": kind,
        "count": len(gaps),
        "interfeatures": gaps,
        "provenance": _provenance(),
    }


def delete_gffutils_features(
    source: str | Path,
    output_path: str | Path,
    *,
    feature_ids: Sequence[str] | None = None,
    featuretype: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete features by ID or all features of one type; the edited database is written to a new ``.db`` file."""

    ids = None if feature_ids is None else [_text(item, "feature_ids") for item in feature_ids]
    kind = _optional_text(featuretype, "featuretype")
    if (ids is None) == (kind is None):
        raise ConfigurationError("give exactly one of feature_ids or featuretype.")
    _, db, target = _new_database(source, output_path)
    before = db.count_features_of_type()
    doomed = [_feature(db, item) for item in ids] if ids is not None else list(db.features_of_type(kind))
    db.delete(doomed, make_backup=False)
    after = db.count_features_of_type()
    return {
        "status": "completed",
        "operation": "delete_gffutils_features",
        "output_path": str(target),
        "deleted": len(doomed),
        "features_before": int(before),
        "features_after": int(after),
        "provenance": _provenance(),
    }


def update_gffutils_features(
    source: str | Path,
    additions: Sequence[str],
    output_path: str | Path,
    *,
    merge_strategy: str = "error",
    verbose: bool = True,
) -> dict[str, Any]:
    """Add GFF lines to an annotation database (``FeatureDB.update``); the result is a new ``.db`` file.

    ``merge_strategy`` is a gffutils strategy: ``error``, ``warning``, ``merge``, ``create_unique`` or ``replace``.
    """

    if isinstance(additions, (str, bytes)) or not additions:
        raise ConfigurationError("additions must be a non-empty list of GFF3 lines.")
    strategies = ("error", "warning", "merge", "create_unique", "replace")
    if merge_strategy not in strategies:
        raise ConfigurationError(f"merge_strategy must be one of {list(strategies)}.")
    lines = [_text(item, "additions", maximum=100_000) for item in additions]
    gffutils, db, target = _new_database(source, output_path)
    before = int(db.count_features_of_type())
    try:
        features = [gffutils.feature.feature_from_line(line) for line in lines]
        db.update(features, make_backup=False, merge_strategy=merge_strategy)
    except Exception as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"gffutils could not add the features: {exc}") from exc
    return {
        "status": "completed",
        "operation": "update_gffutils_features",
        "output_path": str(target),
        "features_before": before,
        "features_after": int(db.count_features_of_type()),
        "merge_strategy": merge_strategy,
        "provenance": _provenance(),
    }


def add_gffutils_relation(
    source: str | Path,
    parent_id: str,
    child_id: str,
    level: int,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a parent-child relation of a given level between two features (``FeatureDB.add_relation``)."""

    if isinstance(level, bool) or not isinstance(level, int) or not 1 <= level <= _MAX_LEVEL:
        raise ConfigurationError(f"level must be an integer in [1, {_MAX_LEVEL}].")
    _, db, target = _new_database(source, output_path)
    parent, child = _feature(db, parent_id), _feature(db, child_id)
    before = len(list(db.children(parent)))
    db.add_relation(parent, child, level)
    children = [item.id for item in db.children(parent)]
    return {
        "status": "completed",
        "operation": "add_gffutils_relation",
        "output_path": str(target),
        "parent_id": parent.id,
        "child_id": child.id,
        "level": level,
        "children_before": before,
        "children_after": len(children),
        "children": children,
        "provenance": _provenance(),
    }


def _derived(
    method: str,
    source: str | Path,
    output_path: str | Path | None,
    limit: int,
    **kwargs: Any,
) -> dict[str, Any]:
    target = None if output_path is None else _new_file(output_path, None, "output_path")
    _, db = _open(source)
    features = []
    try:
        for feature in getattr(db, method)(**kwargs):
            features.append(feature)
            if len(features) > limit:
                break
    except KeyError as exc:
        raise ConfigurationError(
            f"gffutils needs an ID attribute on every exon to derive features; missing {exc}."
        ) from exc
    _checked(features, limit, "derived features")
    if target is not None:
        target.write_text("##gff-version 3\n" + "".join(f"{item}\n" for item in features), encoding="utf-8")
    return {
        "status": "completed",
        "count": len(features),
        "features": [_feature_record(item) for item in features],
        "output_path": None if target is None else str(target),
        "provenance": _provenance(),
    }


def create_gffutils_introns(
    source: str | Path,
    *,
    exon_featuretype: str = "exon",
    grandparent_featuretype: str | None = "gene",
    parent_featuretype: str | None = None,
    new_featuretype: str = "intron",
    merge_attributes: bool = True,
    numeric_sort: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Derive intron features from the gaps between the exons of each transcript (``FeatureDB.create_introns``).

    The new features are returned and, when ``output_path`` (a new file) is given, written as GFF3.
    """

    result = _derived(
        "create_introns",
        source,
        output_path,
        _rows_limit(max_rows),
        exon_featuretype=_text(exon_featuretype, "exon_featuretype"),
        grandparent_featuretype=_optional_text(grandparent_featuretype, "grandparent_featuretype"),
        parent_featuretype=_optional_text(parent_featuretype, "parent_featuretype"),
        new_featuretype=_text(new_featuretype, "new_featuretype"),
        merge_attributes=_flag(merge_attributes, "merge_attributes"),
        numeric_sort=_flag(numeric_sort, "numeric_sort"),
    )
    return {"operation": "create_gffutils_introns", **result}


def create_gffutils_splice_sites(
    source: str | Path,
    *,
    exon_featuretype: str = "exon",
    grandparent_featuretype: str | None = "gene",
    parent_featuretype: str | None = None,
    merge_attributes: bool = True,
    numeric_sort: bool = False,
    output_path: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Derive 5' and 3' splice-site features from the exon boundaries of each transcript (``FeatureDB.create_splice_sites``).

    Every exon must carry an ``ID`` attribute, as gffutils builds the new IDs from them.
    """

    result = _derived(
        "create_splice_sites",
        source,
        output_path,
        _rows_limit(max_rows),
        exon_featuretype=_text(exon_featuretype, "exon_featuretype"),
        grandparent_featuretype=_optional_text(grandparent_featuretype, "grandparent_featuretype"),
        parent_featuretype=_optional_text(parent_featuretype, "parent_featuretype"),
        merge_attributes=_flag(merge_attributes, "merge_attributes"),
        numeric_sort=_flag(numeric_sort, "numeric_sort"),
    )
    return {"operation": "create_gffutils_splice_sites", **result}


def gffutils_children_bp(
    source: str | Path,
    feature_id: str,
    *,
    child_featuretype: str = "exon",
    merge: bool = False,
    ignore_strand: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Total base pairs of the children of one type of a feature, for example exonic bp of a gene (``FeatureDB.children_bp``).

    With ``merge=True`` overlapping children are merged first (same sequence, strand and type).
    """

    _, db = _open(source)
    feature = _feature(db, feature_id)
    kind = _text(child_featuretype, "child_featuretype")
    merge = _flag(merge, "merge")
    ignore_strand = _flag(ignore_strand, "ignore_strand")
    kwargs: dict[str, Any] = {"child_featuretype": kind, "merge": merge}
    if ignore_strand:
        kwargs["ignore_strand"] = True
    base_pairs = int(db.children_bp(feature, **kwargs))
    return {
        "status": "completed",
        "operation": "gffutils_children_bp",
        "feature_id": feature.id,
        "child_featuretype": kind,
        "merge": merge,
        "base_pairs": base_pairs,
        "child_count": len(list(db.children(feature, featuretype=kind))),
        "provenance": _provenance(),
    }


def extract_gffutils_feature_sequence(
    source: str | Path,
    fasta: str | Path,
    feature_id: str,
    *,
    use_strand: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the genome sequence of a feature from a FASTA file (``Feature.sequence``).

    With ``use_strand=True`` a minus-strand feature is reverse-complemented, as gffutils does.
    """

    fasta_path = input_path(fasta)
    _, db = _open(source)
    feature = _feature(db, feature_id)
    try:
        sequence = feature.sequence(str(fasta_path), use_strand=_flag(use_strand, "use_strand"))
    except Exception as exc:
        raise ConfigurationError(f"gffutils could not read the sequence of {feature.id}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "extract_gffutils_feature_sequence",
        "feature_id": feature.id,
        "seqid": feature.seqid,
        "start": int(feature.start),
        "end": int(feature.end),
        "strand": feature.strand,
        "length": len(sequence),
        "sequence": str(sequence),
        "provenance": _provenance(),
    }


def write_gffutils_gff(
    source: str | Path,
    output_path: str | Path,
    *,
    featuretype: str | None = None,
    canonical_gene_order: bool = False,
    gene_featuretype: str = "gene",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write features as GFF records with ``GFFWriter`` (``write_rec`` / ``write_recs``).

    ``canonical_gene_order=True`` instead writes each gene followed by its mRNAs (longest first), their
    exons and the remaining children (``write_gene_recs``); ``featuretype`` is then ignored.
    """

    gffutils, db = _open(source)
    target = _new_file(output_path, None, "output_path")
    kind = _optional_text(featuretype, "featuretype")
    canonical = _flag(canonical_gene_order, "canonical_gene_order")
    genes = _text(gene_featuretype, "gene_featuretype")
    writer = gffutils.gffwriter.GFFWriter(str(target), with_header=True)
    written = 0
    try:
        if canonical:
            for gene in db.features_of_type(genes):
                writer.write_gene_recs(db, gene.id)
        else:
            writer.write_recs(db.features_of_type(kind) if kind else db.all_features())
    finally:
        writer.close()
    written = sum(1 for line in target.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))
    return {
        "status": "completed",
        "operation": "write_gffutils_gff",
        "output_path": str(target),
        "records_written": written,
        "canonical_gene_order": canonical,
        "provenance": _provenance(),
    }


def remove_gffutils_duplicates(
    source: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove duplicate GFF records by merging records that share an ID (``gffutils-cli rmdups``).

    The database is built with ``merge_strategy="merge"`` and every record is written once to ``output_path``.
    """

    gffutils = _gffutils()
    path = input_path(source)
    target = _new_file(output_path, None, "output_path")
    lines_before = sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))
    try:
        db = gffutils.create_db(str(path), ":memory:", verbose=False, merge_strategy="merge")
    except Exception as exc:
        raise ConfigurationError(f"gffutils could not read {path.name}: {exc}") from exc
    writer = gffutils.gffwriter.GFFWriter(str(target), with_header=True)
    try:
        for record in db.all_features():
            writer.write_rec(record)
    finally:
        writer.close()
    lines_after = sum(1 for line in target.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))
    return {
        "status": "completed",
        "operation": "remove_gffutils_duplicates",
        "output_path": str(target),
        "records_before": lines_before,
        "records_after": lines_after,
        "removed": lines_before - lines_after,
        "provenance": _provenance(merge_strategy="merge"),
    }


def plot_gffutils_gene_model(
    source: str | Path,
    gene_id: str,
    output_path: str | Path,
    *,
    transcripts: Sequence[str] = ("mRNA",),
    utrs: Sequence[str] = ("3'UTR", "5'UTR"),
    cds: Sequence[str] = ("CDS",),
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a gene model (transcript line, UTRs, CDS) with ``gffutils.contrib.plotting.Gene`` and save the image.

    ``utrs`` depends on the annotation: GTF uses ``3'UTR``/``5'UTR``, FlyBase GFF ``three_prime_UTR``/``five_prime_UTR``.
    Needs pybedtools and matplotlib. The result holds the image path and the drawn transcripts.
    """

    target = _new_file(output_path, _IMAGE_SUFFIXES, "output_path")
    kinds = {
        name: [_text(item, name) for item in values]
        for name, values in (("transcripts", transcripts), ("utrs", utrs), ("cds", cds))
    }
    require("pybedtools", "bio-pybedtools")
    matplotlib = require("matplotlib", "viz")
    plotting = require("gffutils.contrib.plotting", "bio-gffutils")
    _, db = _open(source)
    gene = _feature(db, gene_id)
    model = plotting.Gene(db, gene.id, transcripts=kinds["transcripts"], utrs=kinds["utrs"], cds=kinds["cds"])
    matplotlib.use("Agg", force=True)
    from matplotlib import pyplot

    figure, axes = pyplot.subplots(figsize=(8, max(1.5, 0.8 * max(1, len(model.tracks) // 3))))
    try:
        model.add_to_ax(axes)
        axes.set_xlim(int(gene.start) - 1, int(gene.end))
        axes.set_ylim(0, model.max_y)
        axes.set_yticks([])
        axes.set_title(gene.id)
        figure.savefig(target, dpi=150, bbox_inches="tight")
    finally:
        pyplot.close(figure)
    transcript_rows = []
    for child in db.children(gene.id, level=1):
        if child.featuretype in kinds["transcripts"]:
            transcript_rows.append(
                {
                    "id": child.id,
                    "start": int(child.start),
                    "end": int(child.end),
                    "utr_count": sum(1 for x in db.children(child, level=1) if x.featuretype in kinds["utrs"]),
                    "cds_count": sum(1 for x in db.children(child, level=1) if x.featuretype in kinds["cds"]),
                }
            )
    return {
        "status": "completed",
        "operation": "plot_gffutils_gene_model",
        "gene_id": gene.id,
        "image_path": str(target),
        "transcripts": transcript_rows,
        "provenance": _provenance(),
    }


__all__ = [
    "add_gffutils_relation",
    "create_gffutils_introns",
    "create_gffutils_splice_sites",
    "delete_gffutils_features",
    "extract_gffutils_feature_sequence",
    "fetch_gffutils_feature",
    "fetch_gffutils_parents",
    "get_gffutils_interfeatures",
    "gffutils_children_bp",
    "group_gffutils_parent_children",
    "list_gffutils_features",
    "plot_gffutils_gene_model",
    "remove_gffutils_duplicates",
    "update_gffutils_features",
    "write_gffutils_gff",
]

_register_public_verbose(globals(), tuple(__all__))
