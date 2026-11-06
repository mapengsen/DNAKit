"""DNA block discovery, records, annotation and figures using geneblocks (P39)."""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError


def _module() -> Any:
    return require("geneblocks", "bio-geneblocks")


def _result(**values: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        "backend": "geneblocks",
        "version": _module().__version__,
        **values,
    }


def _records(
    sequences: Mapping[str, str], annotations: Mapping[str, list[dict[str, Any]]] | None = None
) -> dict[str, Any]:
    biotools = require("geneblocks.biotools", "bio-geneblocks")
    result = {n: biotools.sequence_to_record(s, record_id=n, name=n) for n, s in sequences.items()}
    for name, features in (annotations or {}).items():
        for feature in features:
            biotools.annotate_record(
                result[name],
                location=(feature["start"], feature["end"], feature.get("strand", 1)),
                feature_type=feature.get("type", "misc_feature"),
                **feature.get("qualifiers", {}),
            )
    return result


def _feature(feature: Any) -> dict[str, Any]:
    return {
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "strand": feature.location.strand,
        "type": feature.type,
        "qualifiers": feature.qualifiers,
    }


def _record(record: Any) -> dict[str, Any]:
    return {
        "id": record.id,
        "sequence": str(record.seq),
        "features": [_feature(f) for f in record.features],
    }


def _common(
    sequences: Mapping[str, str],
    min_block_size: int,
    annotations: Mapping[str, list[dict[str, Any]]] | None = None,
) -> Any:
    module = _module()
    if shutil.which("blastn") is None:
        raise BackendUnavailableError(
            "geneblocks common-block discovery requires BLAST+ blastn.",
            hint="conda install -c bioconda blast",
        )
    return module.CommonBlocks.from_sequences(
        _records(sequences, annotations),
        min_block_size=_h.integer(min_block_size, "min_block_size", minimum=1),
    )


def geneblocks_find_common(
    sequences: Mapping[str, str], *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Discover common DNA blocks shared across sequences using upstream all-versus-all BLAST+."""
    blocks = _common(sequences, min_block_size)
    return _result(blocks=blocks.common_blocks)


def geneblocks_unique_segments(
    sequences: Mapping[str, str], *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Coordinates of non-shared segments using CommonBlocks.compute_unique_blocks."""
    return _result(
        segments={
            n: [list(x) for x in v]
            for n, v in _common(sequences, min_block_size).compute_unique_blocks().items()
        }
    )


def geneblocks_export_csv(
    sequences: Mapping[str, str], output_csv: str, *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Write the upstream semicolon-separated common-block table to a new file."""
    target = _h.new_file(output_csv, "output_csv")
    blocks = _common(sequences, min_block_size)
    blocks.common_blocks_to_csv(str(target))
    return _result(outputs=[str(target)], blocks=blocks.common_blocks)


def geneblocks_common_records(
    sequences: Mapping[str, str], *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Export common blocks as JSON representations of upstream sequence records."""
    return _result(
        records=[_record(r) for r in _common(sequences, min_block_size).common_blocks_records()]
    )


def geneblocks_unique_records(
    sequences: Mapping[str, str], *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Export non-shared segments as upstream sequence records including their DNA."""
    return _result(
        records=[_record(r) for r in _common(sequences, min_block_size).unique_blocks_records()]
    )


def geneblocks_annotate_common(
    sequences: Mapping[str, str], *, min_block_size: int = 80, verbose: bool = True
) -> dict[str, Any]:
    """Annotate original DNA sequences with discovered common-block feature coordinates."""
    records = _common(sequences, min_block_size).sequences_with_annotated_blocks()
    return _result(records={n: _record(r) for n, r in records.items()})


def geneblocks_transfer_annotations(
    sequences: Mapping[str, str],
    annotations: Mapping[str, list[dict[str, Any]]],
    *,
    min_block_size: int = 80,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transfer strand-aware features between homologous common blocks of annotated DNA records."""
    records = _common(sequences, min_block_size, annotations).copy_features_between_common_blocks(
        inplace=False
    )
    return _result(records={n: _record(r) for n, r in records.items()})


def geneblocks_plot_common(
    sequences: Mapping[str, str],
    output_image: str,
    *,
    min_block_size: int = 80,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot shared blocks along DNA sequences; return the image path and underlying block data."""
    target = _h.new_file(output_image, "output_image")
    blocks = _common(sequences, min_block_size)
    plt = require("matplotlib.pyplot", "bio-geneblocks")
    fig, axes = plt.subplots(len(sequences), 1, squeeze=False, figsize=(10, 2 * len(sequences)))
    try:
        blocks.plot_common_blocks(axes=list(axes[:, 0]))
        fig.savefig(target)
    finally:
        plt.close(fig)
    return _result(outputs=[str(target)], blocks=blocks.common_blocks)


def _diff(first: str, second: str) -> Any:
    return _module().DiffBlocks.from_sequences(first, second, blast_over=None)


def _blocks(blocks: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "operation": b.operation,
            "first": [int(b.s1_location.start), int(b.s1_location.end)],
            "second": [int(b.s2_location.start), int(b.s2_location.end)],
        }
        for b in blocks
    ]


def geneblocks_merge_differences(
    first: str, second: str, *, change_gap: int = 100, replace_gap: int = 10, verbose: bool = True
) -> dict[str, Any]:
    """Merge nearby edit blocks into larger changes using DiffBlocks.merged."""
    diff = _diff(first, second).merged(change_gap=change_gap, replace_gap=replace_gap)
    return _result(blocks=_blocks(diff.blocks))


def geneblocks_difference_features(
    first: str, second: str, *, sequence: str = "s2", verbose: bool = True
) -> dict[str, Any]:
    """Represent DNA differences as annotated features relative to s1 or s2."""
    _h.choice(sequence, ("s1", "s2"), "sequence")
    return _result(
        features=[_feature(f) for f in _diff(first, second).diffs_as_features(sequence=sequence)]
    )


def geneblocks_plot_differences(
    first: str, second: str, output_image: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Plot DNA edit blocks and return the image path and block coordinates."""
    target = _h.new_file(output_image, "output_image")
    diff = _diff(first, second)
    plt = require("matplotlib.pyplot", "bio-geneblocks")
    fig, axes = plt.subplots(2, 1, figsize=(8, 6))
    try:
        diff.plot(axes=tuple(axes))
        fig.savefig(target)
    finally:
        plt.close(fig)
    return _result(outputs=[str(target)], blocks=_blocks(diff.blocks))


def geneblocks_reconstruct_sequences(
    first: str, second: str, blocks: list[dict[str, Any]], *, verbose: bool = True
) -> dict[str, Any]:
    """Reconstruct both sequences from a complete user-supplied list of edit and equal blocks.

    Each block has operation, first=[start,end], second=[start,end], all 0-based half-open.
    Upstream short-sequence diff discovery omits equal blocks, so callers must supply those explicitly.
    """
    module = _module()
    Location = require("geneblocks.Location", "bio-geneblocks").Location
    DiffBlock = require("geneblocks.DiffBlocks.DiffBlock", "bio-geneblocks").DiffBlock
    native = [
        DiffBlock(
            _h.choice(
                b["operation"], ("equal", "replace", "change", "insert", "delete"), "operation"
            ),
            Location(*b["first"], sequence=first),
            Location(*b["second"], sequence=second),
        )
        for b in blocks
    ]
    a, b = module.DiffBlocks.reconstruct_sequences_from_blocks(native)
    return _result(first=str(a), second=str(b))


def geneblocks_sequence_record(
    sequence: str, *, name: str = "sequence", verbose: bool = True
) -> dict[str, Any]:
    """Wrap DNA as a sequence record using geneblocks.biotools.sequence_to_record."""
    return _result(record=_record(_records({name: sequence})[name]))


def geneblocks_circular_equality(
    sequences: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Test whether DNA strings are rotations of one circular sequence using upstream semantics."""
    module = _module()
    if shutil.which("blastn") is None:
        raise BackendUnavailableError(
            "geneblocks circular equality requires BLAST+ blastn.",
            hint="conda install -c bioconda blast",
        )
    return _result(equal=bool(module.sequences_are_circularly_equal(list(sequences))))


def geneblocks_extract_subsequence(
    sequence: str, start: int, end: int, *, strand: int = 1, verbose: bool = True
) -> dict[str, Any]:
    """Extract a strand-aware 0-based half-open DNA subsequence with geneblocks.Location."""
    _h.choice(strand, (1, -1), "strand")
    Location = require("geneblocks.Location", "bio-geneblocks").Location
    return _result(sequence=str(Location(start, end, strand, sequence=sequence).extract_sequence()))


__all__ = [
    "geneblocks_annotate_common",
    "geneblocks_circular_equality",
    "geneblocks_common_records",
    "geneblocks_difference_features",
    "geneblocks_export_csv",
    "geneblocks_extract_subsequence",
    "geneblocks_find_common",
    "geneblocks_merge_differences",
    "geneblocks_plot_common",
    "geneblocks_plot_differences",
    "geneblocks_reconstruct_sequences",
    "geneblocks_sequence_record",
    "geneblocks_transfer_annotations",
    "geneblocks_unique_records",
    "geneblocks_unique_segments",
]

_register_public_verbose(globals(), tuple(__all__))
