"""Candidate visualization functions backed by the dependency-free SVG renderers of this package.

``plot`` draws sequence-based figures (length histogram, GC distribution, base composition, sequence view, linear/circular map,
label comparison) as valid XML SVG without scripts.  ``export_figure`` saves an SVG (or converts it, with the optional
``cairosvg``/``Pillow``).  ``generate_analysis_summary``, ``generate_qc_html_report`` and ``generate_html_genome_report`` compute
sequence statistics; the HTML reports reuse ``build_html_report``.  ``generate_clinical_pdf_report`` is not implemented.
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from collections.abc import Mapping
from math import ceil
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok, run_candidate_operation
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError
from dnakit.visualization import _impl_plot as _plots

_MODULE = "visualization"
_UNSAFE_SVG = re.compile(r"<\s*script|<\s*foreignObject|\son[a-z]+\s*=|javascript\s*:", re.IGNORECASE)


def _opts(options: Mapping[str, object] | None) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping.", code="INVALID_OPTION")
    return dict(options or {})


def _refuse_overwrite(error: FileExistsError) -> ConfigurationError:
    return ConfigurationError(str(error), code="OUTPUT_EXISTS", hint="Pass options={'overwrite': True} to replace the file.")


def _artifact_result(artifact: object, data: Mapping[str, object], settings: Mapping[str, object]) -> dict[str, object]:
    """Common ``plot`` result; writes the SVG only when ``options["output"]`` is given."""

    from dnakit.visualization.config import SaveConfig
    from dnakit.visualization.export import save_svg

    path = settings.get("output")
    saved = None
    if path is not None:
        overwrite = settings.get("overwrite", False)
        if not isinstance(overwrite, bool):
            raise ConfigurationError("options.overwrite must be a boolean.", code="INVALID_OPTION")
        try:
            save_svg(artifact, str(path), config=SaveConfig(overwrite=overwrite), verbose=False)  # type: ignore[arg-type]
        except FileExistsError as exc:
            raise _refuse_overwrite(exc) from exc
        saved = str(path)
    return {
        "kind": data["kind"], "svg": artifact.svg, "width": artifact.width, "height": artifact.height,  # type: ignore[attr-defined]
        "sha256": artifact.sha256, "data": {k: v for k, v in data.items() if k != "kind"}, "path": saved,  # type: ignore[attr-defined]
    }


def export_figure(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Write an SVG figure to a file, or convert it to PNG/JPG/TIFF/PDF (optional ``cairosvg`` and ``Pillow``).

    ``input``: SVG text, a ``plot`` result (mapping with ``svg``), a ``plot`` envelope, or an ``SVGArtifact``.  The SVG must be
    well-formed XML with an ``<svg>`` root, no DOCTYPE/ENTITY declarations and no scripts, event-handler attributes or
    ``foreignObject``.  Options: ``output`` (required target path; ``.svg`` writes the exact bytes, other suffixes convert),
    ``overwrite`` (default false), ``dpi`` (raster export, default 96).  If the conversion libraries are missing the result is
    ``not_implemented``.  Result: ``path``, ``format``, ``bytes`` and ``sha256`` of the written file.
    """

    from dnakit.visualization.config import ImageExportConfig, SaveConfig
    from dnakit.visualization.export import save_image, save_svg
    from dnakit.visualization.results import SVGArtifact

    settings = _opts(options)
    target = settings.get("output")
    if not isinstance(target, (str, os.PathLike)) or not str(target).strip():
        raise ConfigurationError("options.output (the target path) is required.", code="OUTPUT_REQUIRED")
    overwrite = settings.get("overwrite", False)
    if not isinstance(overwrite, bool):
        raise ConfigurationError("options.overwrite must be a boolean.", code="INVALID_OPTION")
    source = input
    if isinstance(source, Mapping) and isinstance(source.get("result"), Mapping):
        source = source["result"]
    if isinstance(source, Mapping):
        source = source.get("svg")
    if isinstance(source, SVGArtifact):
        artifact = source
    elif isinstance(source, str) and source.lstrip().startswith("<"):
        if _UNSAFE_SVG.search(source):
            raise ConfigurationError("The SVG contains scripts, event handlers or foreignObject.", code="UNSAFE_SVG")
        try:
            root = ET.fromstring(source)
        except ET.ParseError as exc:
            raise InputFormatError(f"The input is not well-formed SVG/XML: {exc}", code="INVALID_SVG") from exc
        if root.tag.rsplit("}", 1)[-1] != "svg":
            raise InputFormatError("The root element must be <svg>.", code="INVALID_SVG")
        width, height = _svg_size(root)
        artifact = SVGArtifact("exported-figure", source, width, height, sum(1 for _ in root.iter()))
    else:
        raise InputFormatError("input must be SVG text, a plot result or an SVGArtifact.", code="INVALID_INPUT")
    path = Path(target)
    suffix = path.suffix.lower().lstrip(".")
    try:
        if suffix == "svg":
            save_svg(artifact, path, config=SaveConfig(overwrite=overwrite), verbose=False)
        elif suffix in {"png", "jpg", "jpeg", "tif", "tiff", "pdf"}:
            save_image(artifact, path, config=ImageExportConfig(overwrite=overwrite, dpi=int(settings.get("dpi", 96))), verbose=False)
        else:
            raise ConfigurationError("The output suffix must be .svg, .png, .jpg, .tiff or .pdf.", code="INVALID_VISUALIZATION_FORMAT")
    except FileExistsError as exc:
        raise _refuse_overwrite(exc) from exc
    except BackendUnavailableError as exc:
        return candidate_not_implemented(_MODULE, "export_figure", input, options=options, requires="pip install cairosvg pillow", message=str(exc))
    import hashlib

    data = path.read_bytes()
    result = {"path": str(path), "format": "tiff" if suffix == "tif" else "jpg" if suffix == "jpeg" else suffix, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    return candidate_ok(_MODULE, "export_figure", input, result, options=options)


def _svg_size(root: ET.Element) -> tuple[int, int]:
    def number(text: str | None) -> float | None:
        match = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*(?:px)?\s*$", text or "")
        return float(match.group(1)) if match else None

    width, height = number(root.get("width")), number(root.get("height"))
    if (width is None or height is None) and root.get("viewBox"):
        parts = root.get("viewBox", "").replace(",", " ").split()
        if len(parts) == 4:
            width, height = width or float(parts[2]), height or float(parts[3])
    if not width or not height or width <= 0 or height <= 0:
        raise InputFormatError("The SVG needs numeric width/height attributes or a viewBox.", code="INVALID_SVG")
    return ceil(width), ceil(height)


def generate_analysis_summary(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Summarize a set of DNA sequences (counts, length statistics, GC content, base counts) as data and Markdown text.

    ``input``: FASTA text/path, a sequence, a list of sequences, ``(id, sequence)`` pairs or mappings, or ``{id: sequence}``.
    Result: ``sequence_count``, ``total_bases``, ``length`` (min, max, mean, median, n50), ``gc_percent`` (100 * (G+C)/(A+C+G+T) over
    all sequences), ``base_counts``, ``non_acgt_bases``, ``n_fraction``, ``per_sequence`` (id, length, gc_percent; at most 1000
    rows) and ``text`` (Markdown).  Option ``title`` (default "Sequence analysis summary").
    """

    settings = _opts(options)
    records = _plots.read_sequences(input)
    result = _plots.summarize(records)
    result["per_sequence"] = [{"id": name, "length": len(seq), "gc_percent": _plots.gc_percent(seq)} for name, seq in records[:1000]]
    result["text"] = _markdown(str(settings.get("title", "Sequence analysis summary")), result)
    return candidate_ok(_MODULE, "generate_analysis_summary", input, result, options=options)


def _markdown(title: str, summary: Mapping[str, object]) -> str:
    length = summary["length"]
    gc = summary["gc_percent"]
    return "\n".join([
        f"# {title}", "", f"- sequences: {summary['sequence_count']}", f"- total bases: {summary['total_bases']}",
        f"- length min / median / mean / max: {length['min']} / {length['median']:g} / {length['mean']:g} / {length['max']}",  # type: ignore[index]
        f"- N50: {length['n50']}", f"- GC content: {'n/a' if gc is None else format(gc, '.2f') + ' %'}",  # type: ignore[index]
        f"- non-ACGT bases: {summary['non_acgt_bases']}", "",
    ])


def _html_report(input: object, options: Mapping[str, object] | None, operation: str, title: str, extra: Mapping[str, object]) -> object:
    from dnakit.core import DNA
    from dnakit.visualization.config import SaveConfig
    from dnakit.visualization.export import save_html_report
    from dnakit.visualization.report import build_html_report

    settings = _opts(options)
    records = _plots.read_sequences(input)
    results: dict[str, object] = dict(extra)
    results["summary"] = _plots.summarize(records)
    try:
        dna_records = [DNA(sequence, id=name).record for name, sequence in records]
    except Exception as exc:  # DNAKit sequence errors: report them as an input problem
        raise InputFormatError(f"A sequence cannot be represented as DNA: {exc}", code="INVALID_SEQUENCE") from exc
    artifact = build_html_report(
        dna_records, results=results, title=str(settings.get("title", title)), max_records=int(settings.get("max_records", 1000)),
        max_total_sequence_symbols=int(settings.get("max_total_sequence_symbols", 10_000_000)), verbose=False,
    )
    path = settings.get("output")
    if path is not None:
        overwrite = settings.get("overwrite", False)
        if not isinstance(overwrite, bool):
            raise ConfigurationError("options.overwrite must be a boolean.", code="INVALID_OPTION")
        try:
            save_html_report(artifact, str(path), config=SaveConfig(overwrite=overwrite), verbose=False)
        except FileExistsError as exc:
            raise _refuse_overwrite(exc) from exc
    result = {"record_count": artifact.record_count, "result_names": list(artifact.result_names), "sha256": artifact.sha256,
              "path": None if path is None else str(path), "summary": results["summary"]}
    for key, value in extra.items():
        result[key] = value
    if settings.get("include_html", True) is not False:
        result["html"] = artifact.html
    return candidate_ok(_MODULE, operation, input, result, options=options)


def generate_clinical_pdf_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Not implemented: a clinical report needs validated clinical interpretation data and a PDF report template."""

    return candidate_not_implemented(
        _MODULE, "generate_clinical_pdf_report", input, options=options,
        requires="validated clinical variant interpretation (ClinVar/ACMG rules) and a certified report template",
        message="Clinical reports must contain validated interpretations, which DNAKit cannot produce here.",
    )


def generate_html_genome_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Build a self-contained HTML report of genome/assembly sequences with assembly statistics (and optional annotation summary).

    ``input``: FASTA text/path or sequences as for ``generate_analysis_summary`` (contigs/chromosomes).  The report lists the
    sequences and the statistics (count, total bases, N50, GC) and, if ``options["annotation"]`` (GFF3/GTF text, path or
    features) is given, the ``summarize`` result of that annotation.  Options: ``title``, ``max_records`` (default 1000),
    ``max_total_sequence_symbols`` (default 10,000,000), ``output`` (``.html`` path; nothing is written without it), ``overwrite``
    (default false), ``include_html`` (default true).  Result: ``record_count``, ``sha256``, ``summary``, ``annotation`` and ``html``.
    """

    settings = _opts(options)
    extra: dict[str, object] = {}
    if settings.get("annotation") is not None:
        from dnakit.annotation import summarize as summarize_annotation

        extra["annotation"] = summarize_annotation(settings["annotation"], verbose=False)["result"]  # type: ignore[index]
    return _html_report(input, options, "generate_html_genome_report", "DNAKit genome report", extra)


def generate_qc_html_report(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Build a self-contained HTML quality-control report of DNA sequences with per-sequence pass/fail flags.

    ``input``: FASTA text/path or sequences as for ``generate_analysis_summary``.  A sequence fails when it is shorter than
    ``min_length`` (default 1), has an N fraction above ``max_n_fraction`` (default 0.1), or a GC percent outside
    [``min_gc_percent``, ``max_gc_percent``] (default 0 to 100; a sequence without A/C/G/T has no GC value and fails a non-default
    range).  The report holds the summary statistics and the ``qc`` verdict (``passed``, ``failed`` with reasons).  Other options
    as in ``generate_html_genome_report`` (``title``, ``max_records``, ``output``, ``overwrite``, ``include_html``).
    """

    settings = _opts(options)
    records = _plots.read_sequences(input)
    minimum, max_n = int(settings.get("min_length", 1)), float(settings.get("max_n_fraction", 0.1))  # type: ignore[call-overload]
    low, high = float(settings.get("min_gc_percent", 0)), float(settings.get("max_gc_percent", 100))  # type: ignore[arg-type]
    failed = []
    for name, sequence in records:
        reasons = []
        if len(sequence) < minimum:
            reasons.append(f"length {len(sequence)} < {minimum}")
        n_fraction = sequence.count("N") / len(sequence) if sequence else 0.0
        if n_fraction > max_n:
            reasons.append(f"N fraction {n_fraction:.3f} > {max_n}")
        gc = _plots.gc_percent(sequence)
        if (gc is None and (low > 0 or high < 100)) or (gc is not None and not low <= gc <= high):
            reasons.append(f"GC percent {'n/a' if gc is None else format(gc, '.2f')} outside [{low}, {high}]")
        if reasons:
            failed.append({"id": name, "reasons": reasons})
    qc = {"passed": len(records) - len(failed), "failed_count": len(failed), "failed": failed,
          "thresholds": {"min_length": minimum, "max_n_fraction": max_n, "min_gc_percent": low, "max_gc_percent": high}}
    return _html_report(input, options, "generate_qc_html_report", "DNAKit sequence QC report", {"qc": qc})


def plot(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Draw a sequence-based figure as a standalone SVG (valid XML, no scripts) and return it with the plotted numbers.

    Option ``kind`` selects the figure; ``input`` is FASTA text/path, a sequence, a list of sequences or ``{id: sequence}``:
    "length_histogram" (default; sequence lengths in nt, ``bins`` default 30), "gc_distribution" (per-sequence
    100 * (G+C)/(A+C+G+T), ``bins`` equal-width bins over 0-100, default 20), "base_composition" (proportion of each base symbol
    over all sequences; A, C, G, T first), "sequence" (symbols of the first sequence, or the one named by ``id``), "linear_map" /
    "circular_map" (first sequence with ``features`` = list of ``{type, start, end}`` in 0-based half-open coordinates; "circular_map"
    always draws the sequence as circular (the sequence is forced to topology "circular" whatever ``topology`` says;
    ``topology`` is only used by the other maps)), and "label_comparison" (``input`` = ``{"true": [...],
    "predicted": [...]}`` class labels, paired proportion bars).  Common options: ``title``, ``output`` (``.svg`` path to write,
    nothing is written without it), ``overwrite`` (default false).  Result: ``kind``, ``svg``, ``width``, ``height``, ``sha256``,
    ``data`` (the exact counts/edges/values drawn) and ``path``.
    """

    from dnakit.visualization.config import DistributionPlotConfig
    from dnakit.visualization.distribution import plot_distribution, plot_label_comparison

    settings = _opts(options)
    kind = str(settings.get("kind", "length_histogram"))
    if kind not in {"length_histogram", "gc_distribution", "base_composition", "sequence", "linear_map", "circular_map", "label_comparison"}:
        raise ConfigurationError(
            f"Unknown plot kind {kind!r}.", code="INVALID_OPTION",
            hint="Use length_histogram, gc_distribution, base_composition, sequence, linear_map, circular_map or label_comparison.",
        )
    title = settings.get("title")

    def config(**values: object) -> DistributionPlotConfig:
        return DistributionPlotConfig(title=None if title is None else str(title), **values)  # type: ignore[arg-type]

    if kind == "label_comparison":
        if not isinstance(input, Mapping) or "true" not in input or "predicted" not in input:
            raise InputFormatError("label_comparison needs {'true': [...], 'predicted': [...]}.", code="INVALID_INPUT")
        artifact = plot_label_comparison(list(input["true"]), list(input["predicted"]), config=config(y_label="Sample proportion (%)"), verbose=False)  # type: ignore[call-overload]
        labels = list(dict.fromkeys([*map(str, input["true"]), *map(str, input["predicted"])]))  # type: ignore[call-overload]
        data = {"kind": kind, "labels": labels, "true_counts": {k: [str(v) for v in input["true"]].count(k) for k in labels},  # type: ignore[call-overload]
                "predicted_counts": {k: [str(v) for v in input["predicted"]].count(k) for k in labels}}  # type: ignore[call-overload]
        return candidate_ok(_MODULE, "plot", input, _artifact_result(artifact, data, settings), options=options)
    records = _plots.read_sequences(input)
    if kind == "length_histogram":
        lengths = [len(sequence) for _, sequence in records]
        artifact = plot_distribution(lengths, task="regression", config=config(bins=int(settings.get("bins", 30)), x_label="Sequence length (nt)", y_label="Sequence proportion (%)"), verbose=False)  # type: ignore[call-overload]
        meta = artifact.metadata
        data = {"kind": kind, "count": len(lengths), "min": min(lengths), "max": max(lengths), "mean": sum(lengths) / len(lengths),
                "bin_edges": list(meta["bin_edges"]), "counts": list(meta["counts"])}
    elif kind == "gc_distribution":
        values = [(name, _plots.gc_percent(sequence)) for name, sequence in records]
        gc = [value for _, value in values if value is not None]
        if not gc:
            raise InputFormatError("No sequence contains an unambiguous A/C/G/T base.", code="INVALID_SEQUENCE")
        bins = int(settings.get("bins", 20))
        if not 1 <= bins <= 1000:
            raise ConfigurationError("bins must be between 1 and 1000.", code="INVALID_OPTION")
        edges = [100.0 * i / bins for i in range(bins + 1)]
        artifact = plot_distribution(gc, task="regression", config=config(bin_edges=edges, x_label="GC content (%)", y_label="Sequence proportion (%)"), verbose=False)
        data = {"kind": kind, "count": len(gc), "skipped": len(values) - len(gc), "gc_percent": gc, "bin_edges": edges, "counts": list(artifact.metadata["counts"])}
    elif kind == "base_composition":
        bases = [symbol for _, sequence in records for symbol in sequence]
        observed = sorted(set(bases) - set("ACGT"))
        categories = [*"ACGT", *observed]
        artifact = plot_distribution(bases, task="classification", config=config(categories=categories, x_label="Base", y_label="Base proportion (%)"), verbose=False)
        counts = {label: count for label, count in zip(artifact.metadata["labels"], artifact.metadata["counts"])}
        data = {"kind": kind, "total_bases": len(bases), "counts": counts, "proportions": {k: v / len(bases) for k, v in counts.items()}}
    else:
        from dnakit.core import DNA
        from dnakit.visualization import plot_circular_map, plot_linear_map, plot_sequence

        wanted = settings.get("id")
        pick = next(((n, s) for n, s in records if wanted is None or n == wanted), None)
        if pick is None:
            raise InputFormatError(f"No sequence with id {wanted!r}.", code="INVALID_INPUT")
        name, sequence = pick
        features = settings.get("features")
        topology = "circular" if kind == "circular_map" else str(settings.get("topology", "linear"))
        dna = DNA(sequence, id=name, features=features if kind != "sequence" else None, topology=topology)  # type: ignore[arg-type]
        if kind == "sequence":
            artifact = plot_sequence(dna, verbose=False)
            data = {"kind": kind, "id": name, "length": len(sequence)}
        else:
            artifact = (plot_linear_map if kind == "linear_map" else plot_circular_map)(dna, title="Linear DNA map" if title is None and kind == "linear_map" else ("Circular DNA map" if title is None else str(title)), verbose=False)
            data = {"kind": kind, "id": name, "length": len(sequence), "feature_count": len(features or [])}  # type: ignore[arg-type]
    return candidate_ok(_MODULE, "plot", input, _artifact_result(artifact, data, settings), options=options)


def save_figure(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    # Not part of ``__all__``: the public ``dnakit.visualization.save_figure`` is the Matplotlib saver in ``figure.py``.
    return run_candidate_operation('visualization', 'save_figure', input, options=options)


__all__ = ['export_figure', 'generate_analysis_summary', 'generate_clinical_pdf_report', 'generate_html_genome_report', 'generate_qc_html_report', 'plot']

from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
