"""scikit-bio file readers/writers, sequence queries and alignment helpers (pack P20, part 1).

Every function imports scikit-bio only when it runs. Sequences must be upper-case IUPAC text,
as scikit-bio requires; nothing is upper-cased silently. Positions are 0-based and intervals are
half-open ``[start, end)`` like Python slices.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _skbio_p20_common as common
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_SEQUENCE_TYPES = ("DNA", "RNA")
_READ_FORMATS = (
    "qseq",
    "clustal",
    "phylip",
    "newick",
    "lsmat",
    "phylip_dm",
    "binary_dm",
    "ordination",
    "sample_metadata",
    "biom",
    "taxdump",
)
_TAXDUMP_SCHEMES = (
    "nodes",
    "names",
    "citations",
    "delnodes",
    "gencode",
    "images",
    "merged",
    "division",
)
_MOTIF_TYPES = ("purine-run", "pyrimidine-run")
_POSITION_KINDS = ("gaps", "degenerates", "definites")
_MAX_SEQUENCE = 50_000_000
_MAX_ALIGNED = 500
_MAX_ALIGNMENT_LENGTH = 1_000_000


def _sequence_class(sequence_type: str) -> Any:
    kind = common.choice(sequence_type, _SEQUENCE_TYPES, "sequence_type")
    return getattr(common.skbio(), kind)


def _make(sequence: object, sequence_type: str, name: str = "sequence") -> Any:
    cls = _sequence_class(sequence_type)
    value = common.text(sequence, name, maximum=_MAX_SEQUENCE)
    try:
        return cls(value)
    except ValueError as exc:
        raise ConfigurationError(f"{name} is not valid {sequence_type}: {exc}") from exc


def _intervals(value: object, length: int, name: str) -> list[tuple[int, int]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of [start, end) pairs.")
    result: list[tuple[int, int]] = []
    for item in value:
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 2:
            raise ConfigurationError(f"{name} must be a list of [start, end) pairs.")
        start = common.integer(item[0], f"{name} start", minimum=0, maximum=length)
        end = common.integer(item[1], f"{name} end", minimum=start, maximum=length)
        result.append((start, end))
    return result


def _mask(intervals: list[tuple[int, int]], length: int) -> Any:
    import numpy as np

    mask = np.zeros(length, dtype=bool)
    for start, end in intervals:
        mask[start:end] = True
    return mask


def _runs(mask: Any) -> list[list[int]]:
    """Convert a boolean mask to half-open [start, end) runs."""

    runs: list[list[int]] = []
    start: int | None = None
    for index, value in enumerate(list(mask)):
        if value and start is None:
            start = index
        elif not value and start is not None:
            runs.append([start, index])
            start = None
    if start is not None:
        runs.append([start, len(mask)])
    return runs


def _aligned_rows(alignment: object, sequence_type: str) -> tuple[list[str], Any]:
    if not isinstance(alignment, Mapping) or not alignment:
        raise ConfigurationError("alignment must be a mapping of name to aligned sequence.")
    if len(alignment) > _MAX_ALIGNED:
        raise ConfigurationError(f"alignment may hold at most {_MAX_ALIGNED} sequences.")
    names = [str(key) for key in alignment]
    cls = _sequence_class(sequence_type)
    rows = []
    for name, value in alignment.items():
        text = common.text(value, f"alignment[{name!r}]", maximum=_MAX_ALIGNMENT_LENGTH)
        try:
            rows.append(cls(text))
        except ValueError as exc:
            raise ConfigurationError(
                f"alignment[{name!r}] is not valid {sequence_type}: {exc}"
            ) from exc
    if len({len(row) for row in rows}) != 1:
        raise ConfigurationError("aligned sequences must all have the same length.")
    msa = common.skbio().TabularMSA(rows, index=names)
    return names, msa


# --------------------------------------------------------------------------- readers / writers
def read_skbio_file(
    path: str,
    file_format: str,
    *,
    sequence_type: str = "DNA",
    qseq_variant: str = "illumina1.8",
    taxdump_scheme: str = "nodes",
    max_items: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a bioinformatics file with scikit-bio's reader and return its content as plain data.

    ``file_format`` selects the reader and therefore the rows of the scikit-bio I/O registry this
    covers: ``qseq`` (Illumina QSEQ reads, ``qseq_variant`` is the quality encoding),
    ``clustal`` and ``phylip`` (alignments, residues typed by ``sequence_type``), ``newick``
    (tree), ``lsmat``, ``phylip_dm`` and ``binary_dm`` (distance matrices), ``ordination``
    (ordination results), ``sample_metadata`` (QIIME 2 metadata TSV), ``biom`` (BIOM table) and
    ``taxdump`` (an NCBI ``*.dmp`` table; ``taxdump_scheme`` names its columns, for example
    ``nodes`` or ``names``). Lists are cut to ``max_items`` entries; ``total`` gives the full count.
    """

    file = common.input_path(path)
    fmt = common.choice(file_format, _READ_FORMATS, "file_format")
    limit = common.integer(max_items, "max_items", minimum=1, maximum=1_000_000)
    package = common.skbio()
    result: dict[str, Any] = {"status": "completed", "operation": "read_skbio_file", "format": fmt}
    try:
        data = _read_one(package, file, fmt, sequence_type, qseq_variant, taxdump_scheme, limit)
    except ConfigurationError:
        raise
    except Exception as exc:  # scikit-bio raises many format-specific error types
        raise ConfigurationError(f"scikit-bio could not read {file.name} as {fmt}: {exc}") from exc
    result.update(data)
    result["provenance"] = common.provenance(method=f"skbio.io.read(format={fmt!r})")
    return result


def _read_one(
    package: Any,
    file: Any,
    fmt: str,
    sequence_type: str,
    qseq_variant: str,
    taxdump_scheme: str,
    limit: int,
) -> dict[str, Any]:
    import itertools

    import pandas as pd

    if fmt == "qseq":
        constructor = _sequence_class(sequence_type)
        reads = list(
            itertools.islice(
                package.io.read(
                    str(file), format="qseq", variant=qseq_variant, constructor=constructor
                ),
                limit + 1,
            )
        )
        items = [
            {
                "sequence": str(read),
                "metadata": common.jsonable(read.metadata),
                "quality": common.jsonable(
                    read.positional_metadata["quality"].to_numpy()
                    if "quality" in read.positional_metadata
                    else None
                ),
            }
            for read in reads[:limit]
        ]
        return {"records": items, "truncated": len(reads) > limit, "variant": qseq_variant}
    if fmt in ("clustal", "phylip"):
        msa = package.io.read(
            str(file),
            format=fmt,
            into=package.TabularMSA,
            constructor=_sequence_class(sequence_type),
        )
        names = [str(i) for i in msa.index]
        return {
            "alignment": {
                n: str(s) for n, s in zip(names[:limit], list(msa)[:limit], strict=False)
            },
            "sequence_count": len(msa),
            "alignment_length": msa.shape.position,
        }
    if fmt == "newick":
        tree = package.io.read(str(file), format="newick", into=package.TreeNode)
        return {
            "newick": str(tree).strip(),
            "tip_count": tree.count(tips=True),
            "node_count": tree.count(),
        }
    if fmt in ("lsmat", "phylip_dm", "binary_dm"):
        matrix = package.io.read(str(file), format=fmt, into=package.DistanceMatrix)
        return {
            "ids": [str(i) for i in matrix.ids],
            "matrix": common.jsonable(matrix.data),
            "size": len(matrix.ids),
        }
    if fmt == "ordination":
        ordination = package.io.read(str(file), format="ordination", into=package.OrdinationResults)
        return {"ordination": common.ordination_to_dict(ordination)}
    if fmt == "sample_metadata":
        metadata = package.SampleMetadata.load(str(file))
        frame = metadata.to_dataframe()
        return {
            "id_column": metadata.id_header,
            "columns": [str(c) for c in frame.columns],
            "column_types": {str(k): v.type for k, v in metadata.columns.items()},
            "ids": [str(i) for i in frame.index[:limit]],
            "data": {
                str(c): common.jsonable(frame[c].iloc[:limit].tolist()) for c in frame.columns
            },
            "total": len(frame),
        }
    if fmt == "biom":
        biom = require("biom", common.EXTRA)
        table = package.io.read(str(file), format="biom", into=biom.Table)
        return {
            "sample_ids": [str(i) for i in table.ids()[:limit]],
            "feature_ids": [str(i) for i in table.ids(axis="observation")[:limit]],
            "shape": [int(table.shape[0]), int(table.shape[1])],
            "counts_features_by_samples": common.jsonable(
                table.matrix_data.toarray()[:limit, :limit]
            ),
            "total": [int(table.shape[0]), int(table.shape[1])],
        }
    scheme = common.choice(taxdump_scheme, _TAXDUMP_SCHEMES, "taxdump_scheme")
    frame = package.io.read(str(file), format="taxdump", into=pd.DataFrame, scheme=scheme)
    return {
        "scheme": scheme,
        "columns": [str(c) for c in frame.columns],
        "index": [str(i) for i in frame.index[:limit]],
        "data": common.jsonable(frame.iloc[:limit].to_numpy()),
        "total": len(frame),
    }


def write_skbio_distance_matrix(
    distance_matrix: Sequence[Sequence[float]],
    output_path: str,
    *,
    ids: Sequence[str] | None = None,
    file_format: str = "lsmat",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a symmetric distance matrix with scikit-bio.

    ``file_format`` is ``lsmat`` (tab-separated text) or ``binary_dm`` (HDF5); the result
    reports the written path, the ids and the size.
    """

    fmt = common.choice(file_format, ("lsmat", "binary_dm"), "file_format")
    matrix = common.distance_matrix(distance_matrix, ids)
    target = common.output_path(output_path)
    try:
        matrix.write(str(target), format=fmt)
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not write {fmt}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_skbio_distance_matrix",
        "format": fmt,
        "output": str(target),
        "ids": [str(i) for i in matrix.ids],
        "size": len(matrix.ids),
        "provenance": common.provenance(method=f"DistanceMatrix.write(format={fmt!r})"),
    }


def write_skbio_ordination(
    sample_coordinates: Sequence[Sequence[float]],
    eigenvalues: Sequence[float],
    output_path: str,
    *,
    sample_ids: Sequence[str] | None = None,
    proportion_explained: Sequence[float] | None = None,
    method: str = "PCoA",
    long_method_name: str = "Principal Coordinate Analysis",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write ordination results (sample coordinates and eigenvalues) in scikit-bio's text format.

    ``sample_coordinates`` has one row per sample and one column per axis; ``eigenvalues`` and
    the optional ``proportion_explained`` have one value per axis.
    """

    import pandas as pd

    package = common.skbio()
    coords = common.matrix(sample_coordinates, "sample_coordinates")
    ids = common.labels(sample_ids, coords.shape[0], "sample_ids", "S")
    axes = coords.shape[1]
    eig = (
        pd.Series(
            common.matrix([list(eigenvalues)], "eigenvalues")[0],
            index=[f"PC{i + 1}" for i in range(axes)],
        )
        if len(eigenvalues) == axes
        else None
    )
    if eig is None:
        raise ConfigurationError("eigenvalues must have one value per coordinate column.")
    explained = None
    if proportion_explained is not None:
        if len(proportion_explained) != axes:
            raise ConfigurationError("proportion_explained must have one value per axis.")
        explained = pd.Series(
            common.matrix([list(proportion_explained)], "proportion_explained")[0], index=eig.index
        )
    result = package.OrdinationResults(
        common.text(method, "method", maximum=64),
        common.text(long_method_name, "long_method_name", maximum=256),
        eig,
        pd.DataFrame(coords, index=ids, columns=eig.index),
        proportion_explained=explained,
    )
    target = common.output_path(output_path)
    result.write(str(target), format="ordination")
    return {
        "status": "completed",
        "operation": "write_skbio_ordination",
        "output": str(target),
        "sample_count": len(ids),
        "axes": axes,
        "provenance": common.provenance(method="OrdinationResults.write(format='ordination')"),
    }


def write_skbio_sample_metadata(
    columns: Mapping[str, Sequence[Any]],
    sample_ids: Sequence[str],
    output_path: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a QIIME 2 sample metadata TSV (with ``#sk:types`` row) via ``SampleMetadata.save``.

    ``columns`` maps a column name to one value per sample; numbers become numeric columns and
    strings categorical columns.
    """

    package = common.skbio()
    if isinstance(sample_ids, (str, bytes)) or not isinstance(sample_ids, Sequence):
        raise ConfigurationError("sample_ids must be a list of strings.")
    ids = common.labels(sample_ids, len(sample_ids), "sample_ids", "S")
    frame = common.dataframe(columns, ids, "columns")
    frame.index.name = "id"
    try:
        metadata = package.SampleMetadata(frame)
    except Exception as exc:
        raise ConfigurationError(f"columns are not valid sample metadata: {exc}") from exc
    target = common.output_path(output_path)
    metadata.save(str(target))
    return {
        "status": "completed",
        "operation": "write_skbio_sample_metadata",
        "output": str(target),
        "sample_count": len(ids),
        "columns": [str(c) for c in frame.columns],
        "provenance": common.provenance(method="SampleMetadata.save"),
    }


def write_skbio_biom_table(
    counts: Sequence[Sequence[float]],
    output_path: str,
    *,
    sample_ids: Sequence[str] | None = None,
    feature_ids: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a count table as a BIOM (HDF5) file with scikit-bio's ``biom`` writer.

    ``counts`` has one row per sample and one column per feature (the usual layout of DNAKit
    microbiome functions); the BIOM file stores features by samples.
    """

    biom = require("biom", common.EXTRA)
    package = common.skbio()
    array = common.matrix(counts, "counts", nonnegative=True)
    samples = common.labels(sample_ids, array.shape[0], "sample_ids", "S")
    features = common.labels(feature_ids, array.shape[1], "feature_ids", "F")
    table = biom.Table(array.T, features, samples)
    target = common.output_path(output_path)
    try:
        package.io.write(table, format="biom", into=str(target))
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not write the BIOM table: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_skbio_biom_table",
        "output": str(target),
        "sample_count": len(samples),
        "feature_count": len(features),
        "provenance": common.provenance(method="skbio.io.write(format='biom')"),
    }


# --------------------------------------------------------------------------- sequence queries
def lowercase_masked_sequence(
    sequence: str,
    mask_intervals: Sequence[Sequence[int]],
    *,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the sequence with the given ``[start, end)`` intervals written in lower case.

    Same as ``skbio.Sequence.lowercase`` with a boolean mask (soft-masking).
    """

    seq = _make(sequence, sequence_type)
    intervals = _intervals(mask_intervals, len(seq), "mask_intervals")
    rendered = seq.lowercase(_mask(intervals, len(seq)))
    return {
        "status": "completed",
        "operation": "lowercase_masked_sequence",
        "sequence": rendered,
        "masked_positions": int(_mask(intervals, len(seq)).sum()),
        "length": len(seq),
        "provenance": common.provenance(method="skbio.Sequence.lowercase"),
    }


def find_subsequence(
    sequence: str,
    subsequence: str,
    *,
    start: int | None = None,
    end: int | None = None,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the first occurrence of ``subsequence`` (``skbio.Sequence.index``).

    ``start`` and ``end`` bound the search like Python's ``str.index``. ``found`` is False and
    ``index`` None when the subsequence is absent.
    """

    seq = _make(sequence, sequence_type)
    sub = _make(subsequence, sequence_type, "subsequence")
    low = None if start is None else common.integer(start, "start", minimum=0, maximum=len(seq))
    high = None if end is None else common.integer(end, "end", minimum=0, maximum=len(seq))
    try:
        index: int | None = int(seq.index(sub, start=low, end=high))
    except ValueError:
        index = None
    return {
        "status": "completed",
        "operation": "find_subsequence",
        "found": index is not None,
        "index": index,
        "provenance": common.provenance(method="skbio.Sequence.index"),
    }


def split_contiguous_segments(
    sequence: str,
    *,
    included_intervals: Sequence[Sequence[int]] | None = None,
    included_characters: str | None = None,
    min_length: int = 1,
    invert: bool = False,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a sequence into contiguous segments selected by a mask (``Sequence.iter_contiguous``).

    The mask is either ``included_intervals`` (``[start, end)`` pairs) or every position whose
    character is in ``included_characters`` (for example ``"N"``). ``invert`` selects the
    complement and ``min_length`` drops shorter segments.
    """

    seq = _make(sequence, sequence_type)
    if (included_intervals is None) == (included_characters is None):
        raise ConfigurationError("give exactly one of included_intervals and included_characters.")
    if included_intervals is not None:
        mask = _mask(_intervals(included_intervals, len(seq), "included_intervals"), len(seq))
    else:
        chars = common.text(included_characters, "included_characters", maximum=64)
        mask = [c in chars for c in str(seq)]
    minimum = common.integer(min_length, "min_length", minimum=1, maximum=max(len(seq), 1))
    flag_invert = common.flag(invert, "invert")
    chosen = [bool(v) != flag_invert for v in list(mask)]
    runs = [r for r in _runs(chosen) if r[1] - r[0] >= minimum]
    pieces = [
        str(piece) for piece in seq.iter_contiguous(mask, min_length=minimum, invert=flag_invert)
    ]
    segments = [
        {"start": start, "end": end, "sequence": piece}
        for (start, end), piece in zip(runs, pieces, strict=True)
    ]
    return {
        "status": "completed",
        "operation": "split_contiguous_segments",
        "segments": segments,
        "segment_count": len(segments),
        "provenance": common.provenance(method="skbio.Sequence.iter_contiguous"),
    }


def locate_sequence_positions(
    sequence: str,
    kind: str,
    *,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate gap, degenerate or definite positions of an IUPAC sequence.

    ``kind`` is ``gaps`` (``-`` and ``.``), ``degenerates`` (IUPAC ambiguity codes such as N, R,
    Y) or ``definites`` (the unambiguous bases). The result lists 0-based ``positions``, their
    ``runs`` as ``[start, end)`` and ``has_definites`` (``GrammaredSequence.has_definites``).
    """

    seq = _make(sequence, sequence_type)
    chosen = common.choice(kind, _POSITION_KINDS, "kind")
    mask = getattr(seq, chosen)()
    positions = [int(i) for i, v in enumerate(list(mask)) if v]
    return {
        "status": "completed",
        "operation": "locate_sequence_positions",
        "kind": chosen,
        "positions": positions,
        "count": len(positions),
        "runs": _runs(mask),
        "length": len(seq),
        "has_definites": bool(seq.has_definites()),
        "provenance": common.provenance(method=f"skbio.GrammaredSequence.{chosen}"),
    }


def find_nucleotide_motifs(
    sequence: str,
    motif_type: str,
    *,
    min_length: int = 1,
    ignore_intervals: Sequence[Sequence[int]] | None = None,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Find runs of purines (``purine-run``) or pyrimidines (``pyrimidine-run``).

    Wraps ``GrammaredSequence.find_motifs`` (also available on DNA and RNA); ``ignore_intervals``
    marks ``[start, end)`` positions to skip and ``min_length`` filters short runs.
    """

    seq = _make(sequence, sequence_type)
    kind = common.choice(motif_type, _MOTIF_TYPES, "motif_type")
    minimum = common.integer(min_length, "min_length", minimum=1, maximum=max(len(seq), 1))
    ignore = None
    if ignore_intervals is not None:
        ignore = _mask(_intervals(ignore_intervals, len(seq), "ignore_intervals"), len(seq))
    motifs = [
        {"start": int(s.start), "end": int(s.stop), "sequence": str(seq[s])}
        for s in seq.find_motifs(kind, min_length=minimum, ignore=ignore)
    ]
    return {
        "status": "completed",
        "operation": "find_nucleotide_motifs",
        "motif_type": kind,
        "motifs": motifs,
        "count": len(motifs),
        "provenance": common.provenance(method="skbio.GrammaredSequence.find_motifs"),
    }


def is_reverse_complement_pair(
    sequence_a: str,
    sequence_b: str,
    *,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether two sequences are reverse complements of each other (IUPAC aware)."""

    first = _make(sequence_a, sequence_type, "sequence_a")
    second = _make(sequence_b, sequence_type, "sequence_b")
    return {
        "status": "completed",
        "operation": "is_reverse_complement_pair",
        "is_reverse_complement": bool(first.is_reverse_complement(second)),
        "provenance": common.provenance(method="skbio.NucleotideMixin.is_reverse_complement"),
    }


# --------------------------------------------------------------------------- alignments
def msa_consensus_sequence(
    alignment: Mapping[str, str],
    *,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Majority consensus of an aligned set (``TabularMSA.consensus``, scikit-bio tie rule)."""

    _, msa = _aligned_rows(alignment, sequence_type)
    return {
        "status": "completed",
        "operation": "msa_consensus_sequence",
        "consensus": str(msa.consensus()),
        "sequence_count": len(msa),
        "alignment_length": msa.shape.position,
        "provenance": common.provenance(method="skbio.TabularMSA.consensus"),
    }


def msa_conservation(
    alignment: Mapping[str, str],
    *,
    metric: str = "inverse_shannon_uncertainty",
    degenerate_mode: str = "error",
    gap_mode: str = "nan",
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-column conservation of an alignment (``TabularMSA.conservation``).

    Defaults reproduce scikit-bio. Columns that cannot be scored (gaps with ``gap_mode='nan'``)
    are reported as null.
    """

    _, msa = _aligned_rows(alignment, sequence_type)
    values = msa.conservation(
        metric=common.choice(metric, ("inverse_shannon_uncertainty",), "metric"),
        degenerate_mode=common.choice(degenerate_mode, ("nan", "error"), "degenerate_mode"),
        gap_mode=common.choice(gap_mode, ("nan", "ignore", "error", "include"), "gap_mode"),
    )
    return {
        "status": "completed",
        "operation": "msa_conservation",
        "metric": metric,
        "conservation": common.jsonable(values),
        "provenance": common.provenance(method="skbio.TabularMSA.conservation"),
    }


def msa_gap_frequencies(
    alignment: Mapping[str, str],
    *,
    axis: str = "sequence",
    relative: bool = False,
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Count (or give the fraction of) gap characters per column or per sequence.

    ``axis='sequence'`` returns one value per alignment position, ``axis='position'`` one value
    per sequence (``TabularMSA.gap_frequencies``).
    """

    names, msa = _aligned_rows(alignment, sequence_type)
    chosen = common.choice(axis, ("sequence", "position"), "axis")
    values = msa.gap_frequencies(axis=chosen, relative=common.flag(relative, "relative"))
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "msa_gap_frequencies",
        "axis": chosen,
        "relative": relative,
        "frequencies": common.jsonable(values),
        "provenance": common.provenance(method="skbio.TabularMSA.gap_frequencies"),
    }
    if chosen == "position":
        result["names"] = names
    return result


def parse_cigar_alignment_path(
    cigar: str,
    *,
    starts: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a CIGAR string into a pairwise alignment path (``PairAlignPath.from_cigar``).

    Returns run lengths, per-run states (0 match, 1 deletion from the first sequence, 2 insertion,
    in scikit-bio's encoding), the start offsets, the alignment length and the CIGAR round trip.
    """

    package = common.skbio()
    text = common.text(cigar, "cigar", maximum=1_000_000)
    begin = None
    if starts is not None:
        begin = [common.integer(v, "starts item", minimum=0, maximum=2**40) for v in starts]
        if len(begin) != 2:
            raise ConfigurationError("starts must hold two offsets.")
    try:
        path = package.alignment.PairAlignPath.from_cigar(text, starts=begin)
    except Exception as exc:
        raise ConfigurationError(f"cigar is not valid: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_cigar_alignment_path",
        "lengths": common.jsonable(path.lengths),
        "states": common.jsonable(path.states),
        "starts": common.jsonable(path.starts),
        "sequence_count": int(path.shape[0]),
        "alignment_length": int(path.shape[1]),
        "cigar": path.to_cigar(),
        "provenance": common.provenance(method="skbio.alignment.PairAlignPath.from_cigar"),
    }


def apply_alignment_path(
    sequences: Mapping[str, str],
    *,
    cigar: str | None = None,
    lengths: Sequence[int] | None = None,
    states: Sequence[Sequence[int]] | None = None,
    starts: Sequence[int] | None = None,
    gap_char: str = "-",
    sequence_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Rebuild the gapped alignment of the original sequences from an alignment path.

    The path is a pairwise ``cigar`` (two sequences) or, for any number of sequences, run
    ``lengths`` with ``states`` (one row of bit-coded run states per sequence; 1 marks a gap).
    ``starts`` gives each sequence's first aligned position. Equivalent to
    ``TabularMSA.from_path_seqs`` / ``AlignPath.to_aligned`` with ``gap_char``.
    """

    package = common.skbio()
    if not isinstance(sequences, Mapping) or not 1 <= len(sequences) <= _MAX_ALIGNED:
        raise ConfigurationError("sequences must be a mapping of name to unaligned sequence.")
    names = [str(k) for k in sequences]
    cls = _sequence_class(sequence_type)
    seqs = []
    for name, value in sequences.items():
        text = common.text(value, f"sequences[{name!r}]", maximum=_MAX_ALIGNMENT_LENGTH)
        try:
            seqs.append(cls(text))
        except ValueError as exc:
            raise ConfigurationError(
                f"sequences[{name!r}] is not valid {sequence_type}: {exc}"
            ) from exc
    if (cigar is None) == (lengths is None):
        raise ConfigurationError("give either cigar or lengths with states.")
    gap = common.text(gap_char, "gap_char", maximum=1)
    try:
        if cigar is not None:
            if len(seqs) != 2:
                raise ConfigurationError("a CIGAR path needs exactly two sequences.")
            begin = None if starts is None else list(starts)
            path = package.alignment.PairAlignPath.from_cigar(
                common.text(cigar, "cigar"), starts=begin
            )
        else:
            if lengths is None or states is None or starts is None:
                raise ConfigurationError("lengths requires states and starts.")
            runs = [list(run) for run in states]
            if len(runs) != len(lengths) or any(len(run) != len(seqs) for run in runs):
                raise ConfigurationError(
                    "states needs one 0/1 gap flag per sequence for every run."
                )
            rows = [
                [
                    sum((1 if run[8 * r + k] else 0) << k for k in range(min(8, len(seqs) - 8 * r)))
                    for run in runs
                ]
                for r in range((len(seqs) + 7) // 8)
            ]
            path = package.alignment.AlignPath(
                lengths=list(lengths), states=rows, starts=list(starts)
            )
        aligned = path.to_aligned(seqs, gap_char=gap)
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"the alignment path does not fit the sequences: {exc}") from exc
    return {
        "status": "completed",
        "operation": "apply_alignment_path",
        "alignment": dict(zip(names, aligned, strict=True)),
        "alignment_length": int(path.shape[1]),
        "provenance": common.provenance(method="skbio.alignment.AlignPath.to_aligned"),
    }


__all__ = [
    "apply_alignment_path",
    "find_nucleotide_motifs",
    "find_subsequence",
    "is_reverse_complement_pair",
    "locate_sequence_positions",
    "lowercase_masked_sequence",
    "msa_consensus_sequence",
    "msa_conservation",
    "msa_gap_frequencies",
    "parse_cigar_alignment_path",
    "read_skbio_file",
    "split_contiguous_segments",
    "write_skbio_biom_table",
    "write_skbio_distance_matrix",
    "write_skbio_ordination",
    "write_skbio_sample_metadata",
]

_register_public_verbose(globals(), tuple(__all__))
