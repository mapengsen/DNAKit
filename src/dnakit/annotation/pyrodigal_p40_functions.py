"""pyrodigal gene-finder operations as DNAKit functions (pack P40).

The functions call the pyrodigal Python package (a Cython re-implementation of Prodigal) in process on
an explicit local FASTA file. Defaults follow Prodigal: single mode trains on the input (all records are
joined by Prodigal's linker), translation table 11, closed ends off, minimum gene length 90 nt. pyrodigal is
imported only when a function runs. Outputs written to disk go to new paths; nothing is overwritten.
"""

from __future__ import annotations

import io
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-pyrodigal"
_OUTPUTS = ("nucleotides", "translations", "genbank", "gff", "scores")
_OUTPUT_SUFFIX = {"nucleotides": "genes.fna", "translations": "proteins.faa", "genbank": "genes.gbk",
                  "gff": "genes.gff", "scores": "scores.tsv"}


def _lib() -> Any:
    return require("pyrodigal", _EXTRA)


def _records(fasta: str | Path) -> list[tuple[str, str]]:
    return [(rid, seq.upper()) for rid, _d, seq in _h.read_fasta(_h.existing(fasta, "fasta"))]


def _training(lib: Any, records: Sequence[tuple[str, str]], training_file: str | Path | None, *,
              translation_table: int, force_nonsd: bool, start_weight: float | None) -> Any:
    if training_file is not None:
        with open(_h.existing(training_file, "training_file"), "rb") as handle:
            return lib.TrainingInfo.load(handle)
    finder = lib.GeneFinder()
    kwargs: dict[str, Any] = {"force_nonsd": force_nonsd, "translation_table": translation_table}
    if start_weight is not None:
        kwargs["start_weight"] = float(start_weight)
    try:
        return finder.train(*[seq for _id, seq in records], **kwargs)
    except ValueError as exc:
        raise ConfigurationError(f"pyrodigal could not train on this input: {exc}") from exc


def _info_summary(info: Any) -> dict[str, Any]:
    return {
        "gc": float(info.gc),
        "translation_table": int(info.translation_table),
        "uses_shine_dalgarno": bool(info.uses_sd),
        "start_weight": float(info.start_weight),
        "bias": [float(x) for x in info.bias],
        "type_weights": [float(x) for x in info.type_weights],
    }


def pyrodigal_train_gene_finder(
    fasta: str | Path,
    *,
    translation_table: int = 11,
    force_nonsd: bool = False,
    start_weight: float | None = None,
    training_output: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train the Prodigal gene finder on a genome FASTA (pyrodigal GeneFinder.train, ``pyrodigal -t``).

    All records are treated as contigs of one genome. ``training_output`` optionally writes the training
    info in Prodigal's binary format (machine specific, as in Prodigal); the result lists the GC content,
    whether a Shine-Dalgarno motif is used, start weight, bias and start-type weights.
    """

    lib = _lib()
    table = _h.integer(translation_table, "translation_table", minimum=1, maximum=33)
    records = _records(fasta)
    info = _training(lib, records, None, translation_table=table, force_nonsd=bool(force_nonsd),
                     start_weight=start_weight)
    result: dict[str, Any] = {"status": "completed", "backend": "pyrodigal", "version": lib.__version__,
                              "record_count": len(records), "training": _info_summary(info), "outputs": []}
    if training_output is not None:
        target = _h.new_file(training_output, "training_output")
        with open(target, "wb") as handle:
            info.dump(handle)
        result["outputs"].append(str(target))
    return result


def pyrodigal_gc_frame_plot(
    fasta: str | Path,
    *,
    window_size: int = 120,
    record: str | None = None,
    max_segments: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Maximum-GC frame plot of one sequence (pyrodigal Sequence.max_gc_frame_plot).

    For every position the frame (0, 1 or 2) with the highest GC content inside a sliding window is
    recorded; the result gives the per-frame counts and the run-length encoded segments (bounded).
    Upstream -1 values (no frame, for example sequences shorter than the window)
    are preserved in segments and counted separately in ``unknown_frame_count``.
    """

    lib = _lib()
    window = _h.integer(window_size, "window_size", minimum=1)
    limit = _h.integer(max_segments, "max_segments", minimum=1, maximum=_h.MAX_ROWS)
    records = _records(fasta)
    chosen = _pick(records, record)
    plot = lib.Sequence(chosen[1]).max_gc_frame_plot(window)
    counts = [0, 0, 0]
    unknown_count = 0
    total_segments = 0
    segments: list[dict[str, int]] = []
    start = 0
    for index, frame in enumerate(plot):
        if frame == -1:
            unknown_count += 1
        else:
            counts[frame] += 1
        if index + 1 == len(plot) or plot[index + 1] != frame:
            total_segments += 1
            if len(segments) < limit:
                segments.append({"start": start, "end": index + 1, "frame": int(frame)})
            start = index + 1
    return {"status": "completed", "backend": "pyrodigal", "version": lib.__version__, "record": chosen[0],
            "window_size": window, "length": len(plot), "frame_counts": counts, "unknown_frame_count": unknown_count,
            "segments": segments, "segments_truncated": total_segments > len(segments)}


def _pick(records: Sequence[tuple[str, str]], record: str | None) -> tuple[str, str]:
    if record is None:
        return records[0]
    for item in records:
        if item[0] == record:
            return item
    raise ConfigurationError(f"record {record!r} is not in the FASTA file.")


def pyrodigal_shine_dalgarno(
    fasta: str | Path,
    *,
    start: int,
    position: int | None = None,
    strand: int = 1,
    exact: bool = False,
    record: str | None = None,
    training_file: str | Path | None = None,
    translation_table: int = 11,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score the Shine-Dalgarno motif upstream of a start codon (pyrodigal Sequence.shine_dalgarno).

    ``start`` is the 0-based position of the start codon on the sequence and ``position`` where to scan
    from (default: the base just upstream of ``start``). ``exact=True`` scores only an exact AGGAGG,
    otherwise one mismatch is allowed. The weights come from ``training_file`` (a Prodigal training file) or
    from training on the input. Returns the index of the best motif and its weight.
    """

    lib = _lib()
    records = _records(fasta)
    chosen = _pick(records, record)
    if strand not in (1, -1):
        raise ConfigurationError("strand must be 1 or -1.")
    begin = _h.integer(start, "start", minimum=0, maximum=len(chosen[1]) - 1)
    pos = begin - strand if position is None else _h.integer(position, "position", minimum=0)
    info = _training(lib, records, training_file, translation_table=translation_table, force_nonsd=False,
                     start_weight=None)
    try:
        index = int(lib.Sequence(chosen[1]).shine_dalgarno(pos, begin, info, strand=strand, exact=bool(exact)))
    except ValueError as exc:
        raise ConfigurationError(f"pyrodigal rejected the positions: {exc}") from exc
    weights = [float(x) for x in info.rbs_weights]
    return {"status": "completed", "backend": "pyrodigal", "version": lib.__version__, "record": chosen[0],
            "start": begin, "position": pos, "strand": strand, "exact": bool(exact),
            "motif_index": index, "motif_weight": weights[index] if index < len(weights) else None}


def pyrodigal_codon_probability(
    fasta: str | Path,
    *,
    kind: str = "start",
    record: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Start- or stop-codon probability expected from the GC content of a sequence.

    ``kind`` is ``"start"`` (pyrodigal Sequence.start_probability) or ``"stop"``
    (Sequence.stop_probability).
    """

    lib = _lib()
    which = _h.choice(kind, ("start", "stop"), "kind")
    chosen = _pick(_records(fasta), record)
    seq = lib.Sequence(chosen[1])
    value = seq.start_probability() if which == "start" else seq.stop_probability()
    return {"status": "completed", "backend": "pyrodigal", "version": lib.__version__, "record": chosen[0],
            "kind": which, "gc": float(seq.gc), "probability": float(value)}


def pyrodigal_candidate_nodes(
    fasta: str | Path,
    *,
    score: bool = True,
    record: str | None = None,
    translation_table: int = 11,
    closed: bool = False,
    training_file: str | Path | None = None,
    max_nodes: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract (and optionally score) the candidate start/stop nodes of a sequence.

    Covers pyrodigal Nodes.extract (``score=False``) and Nodes.extract followed by Nodes.score
    (``score=True``, needs the training info trained on the input or read from ``training_file``).
    """

    lib = _lib()
    table = _h.integer(translation_table, "translation_table", minimum=1, maximum=33)
    limit = _h.integer(max_nodes, "max_nodes", minimum=1, maximum=_h.MAX_ROWS)
    records = _records(fasta)
    chosen = _pick(records, record)
    sequence = lib.Sequence(chosen[1])
    nodes = lib.Nodes()
    count = nodes.extract(sequence, closed=bool(closed), translation_table=table)
    if score:
        info = _training(lib, records, training_file, translation_table=table, force_nonsd=False,
                         start_weight=None)
        nodes.score(sequence, info, closed=bool(closed))
    nodes.sort()
    rows = []
    for index in range(min(int(nodes.length), limit)):
        node = nodes[index]
        row: dict[str, Any] = {"index": int(node.index), "type": str(node.type), "strand": int(node.strand),
                               "edge": bool(node.edge)}
        if score:
            row.update({"coding_score": float(node.cscore), "rbs_score": float(node.rscore),
                        "start_score": float(node.sscore), "total_score": float(node.score)})
        rows.append(row)
    return {"status": "completed", "backend": "pyrodigal", "version": lib.__version__, "record": chosen[0],
            "scored": bool(score), "node_count": int(count), "nodes": rows,
            "truncated": int(nodes.length) > limit}


def pyrodigal_find_genes(
    fasta: str | Path,
    *,
    mode: str = "single",
    outputs: Sequence[str] = (),
    output_dir: str | Path | None = None,
    translation_table: int = 11,
    closed: bool = False,
    mask: bool = False,
    min_gene: int | None = None,
    training_file: str | Path | None = None,
    max_genes: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict genes with pyrodigal and write the chosen annotation files.

    ``mode`` is ``"single"`` (train on the input, or use ``training_file``) or ``"meta"`` (pre-trained
    metagenomic models). The result lists every gene with its coordinates (1-based, inclusive), strand,
    start codon, RBS motif, scores and the prediction confidence (pyrodigal Gene.confidence). ``outputs``
    selects files written into the new ``output_dir``: ``nucleotides`` (Gene.sequence / write_genes,
    ``pyrodigal -d``), ``translations`` (Gene.translate / write_translations, ``-a``), ``genbank``
    (write_genbank, ``-f gbk``), ``gff`` (write_gff, ``-f gff``) and ``scores`` (write_scores, ``-s``).
    """

    lib = _lib()
    which = _h.choice(mode, ("single", "meta"), "mode")
    wanted = [_h.choice(item, _OUTPUTS, "outputs") for item in outputs]
    limit = _h.integer(max_genes, "max_genes", minimum=1, maximum=_h.MAX_ROWS)
    table = _h.integer(translation_table, "translation_table", minimum=1, maximum=33)
    records = _records(fasta)
    if wanted and output_dir is None:
        raise ConfigurationError("output_dir is required when outputs are requested.")
    kwargs: dict[str, Any] = {"closed": bool(closed), "mask": bool(mask)}
    if min_gene is not None:
        kwargs["min_gene"] = _h.integer(min_gene, "min_gene", minimum=1)
    if which == "meta":
        finder = lib.GeneFinder(meta=True, **kwargs)
    else:
        info = _training(lib, records, training_file, translation_table=table, force_nonsd=False,
                         start_weight=None)
        finder = lib.GeneFinder(info, **kwargs)
    directory = _h.new_dir(output_dir, "output_dir") if wanted else None
    handles = {name: io.StringIO() for name in wanted}
    genes_out: list[dict[str, Any]] = []
    total = 0
    for number, (rid, seq) in enumerate(records):
        genes = finder.find_genes(seq)
        total += len(genes)
        for gene in genes:
            if len(genes_out) < limit:
                genes_out.append({
                    "record": rid, "begin": int(gene.begin), "end": int(gene.end), "strand": int(gene.strand),
                    "start_type": str(gene.start_type), "partial_begin": bool(gene.partial_begin),
                    "partial_end": bool(gene.partial_end), "rbs_motif": gene.rbs_motif,
                    "rbs_spacer": gene.rbs_spacer, "gc_content": float(gene.gc_cont),
                    "score": float(gene.score), "coding_score": float(gene.cscore),
                    "confidence": float(gene.confidence()), "length": int(gene.end - gene.begin + 1)})
        if "nucleotides" in handles:
            genes.write_genes(handles["nucleotides"], sequence_id=rid)
        if "translations" in handles:
            genes.write_translations(handles["translations"], sequence_id=rid)
        if "genbank" in handles:
            genes.write_genbank(handles["genbank"], sequence_id=rid)
        if "gff" in handles:
            genes.write_gff(handles["gff"], sequence_id=rid, header=number == 0)
        if "scores" in handles:
            genes.write_scores(handles["scores"], sequence_id=rid, header=number == 0)
    files = []
    for name, buffer in handles.items():
        target = directory / _OUTPUT_SUFFIX[name]  # type: ignore[operator]
        target.write_text(buffer.getvalue(), encoding="utf-8")
        files.append(str(target))
    return {"status": "completed", "backend": "pyrodigal", "version": lib.__version__, "mode": which,
            "record_count": len(records), "gene_count": total, "genes": genes_out,
            "truncated": total > len(genes_out), "outputs": files}


__all__ = [
    "pyrodigal_candidate_nodes",
    "pyrodigal_codon_probability",
    "pyrodigal_find_genes",
    "pyrodigal_gc_frame_plot",
    "pyrodigal_shine_dalgarno",
    "pyrodigal_train_gene_finder",
]

_register_public_verbose(globals(), tuple(__all__))
