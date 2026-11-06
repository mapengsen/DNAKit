"""Biopython ``Bio.Graphics`` figures added by integration pack P19: genome diagrams with cross-links and
chromosome ideograms.

Both functions draw with ReportLab through Biopython and write one new file (``pdf``, ``svg`` or ``eps``;
``png`` needs the ReportLab Cairo backend ``rlPyCairo``). They return the file path together with the data
that was drawn. Colours are CSS/ReportLab colour names (``red``, ``lightgreen``) or ``#RRGGBB``. Biopython
and ReportLab are imported only when a function runs (extras ``bio-formats`` and ``reports``).
"""

from __future__ import annotations

import importlib.metadata
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-formats"
_REPORT_EXTRA = "reports"
_PACKAGE = "biopython"
_FORMATS = ("pdf", "svg", "eps", "png")
_NAME = re.compile(r"[^\x00-\x1f]{1,100}")
_HEX = re.compile(r"#[0-9A-Fa-f]{6}")
_MAX_TRACKS = 20
_MAX_FEATURES = 20_000
_MAX_LENGTH = 10**10
_PNG_HINT = "Install the ReportLab Cairo backend: python -m pip install rlPyCairo (or choose pdf, svg or eps)."


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _label(value: object, name: str) -> str:
    if not isinstance(value, str) or not _NAME.fullmatch(value):
        raise ConfigurationError(f"{name} must be a non-empty text without control characters.")
    return value


def _color(value: object, name: str, default: str) -> Any:
    colors = require("reportlab.lib.colors", _REPORT_EXTRA)
    text = default if value is None else value
    if not isinstance(text, str):
        raise ConfigurationError(f"{name} must be a colour name or #RRGGBB.")
    if _HEX.fullmatch(text):
        return colors.HexColor(text)
    named = colors.getAllNamedColors()
    if text.lower() in named:
        return named[text.lower()]
    raise ConfigurationError(f"{name}: {text!r} is not a known colour name or #RRGGBB.")


def _strand(value: object) -> int:
    if value is None or value in (0, "."):
        return 0
    if value in (1, "+", "1"):
        return 1
    if value in (-1, "-", "-1"):
        return -1
    raise ConfigurationError("strand must be +1, -1, 0 or one of '+', '-', '.'.")


def _format(value: object) -> str:
    if not isinstance(value, str) or value.lower() not in _FORMATS:
        raise ConfigurationError(f"image_format must be one of {list(_FORMATS)}.")
    return value.lower()


def _feature(item: object, length: int, where: str) -> dict[str, Any]:
    if not isinstance(item, Mapping):
        raise ConfigurationError(f"{where}: each feature must be a mapping with start and end.")
    unknown = set(item) - {"start", "end", "strand", "label", "color"}
    if unknown:
        raise ConfigurationError(f"{where}: unknown feature keys {sorted(unknown)}.")
    start = _integer(item.get("start"), f"{where}.start", minimum=0, maximum=length)
    end = _integer(item.get("end"), f"{where}.end", minimum=start + 1, maximum=length)
    label = item.get("label")
    return {
        "start": start,
        "end": end,
        "strand": _strand(item.get("strand")),
        "label": _label(label, f"{where}.label") if label is not None else None,
        "color": item.get("color"),
    }


def _write_error(error: Exception, fmt: str, target: Path) -> Exception:
    target.unlink(missing_ok=True)
    if fmt == "png":
        return BackendUnavailableError(f"PNG output is not available: {error}", hint=_PNG_HINT)
    return ConfigurationError(f"Biopython could not draw the figure: {error}")


def genome_diagram_crosslinks(
    tracks: Sequence[Mapping[str, Any]],
    links: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    image_format: str = "pdf",
    title: str | None = None,
    width_cm: float = 24.0,
    track_height_cm: float = 3.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw genome tracks as stacked bands joined by cross-links (``Bio.Graphics.GenomeDiagram.CrossLink``).

    ``tracks`` is a list of ``{"name", "length", "features": [{"start", "end", "strand", "label",
    "color"}]}`` (0-based half-open coordinates; strand +1/-1/0 draws arrow sigils). ``links`` is a list of
    ``{"track_a", "start_a", "end_a", "track_b", "start_b", "end_b", "color", "flip"}`` where ``track_a`` and
    ``track_b`` are track names or indices and must differ; ``flip=True`` draws an inverted ribbon (a
    reversed match). All tracks are drawn to the scale of the longest. The figure is written to the new
    file ``output_path``. Returns the path and the drawn tracks and links.
    """

    fmt = _format(image_format)
    if (
        isinstance(tracks, (str, bytes))
        or not isinstance(tracks, Sequence)
        or not 2 <= len(tracks) <= _MAX_TRACKS
    ):
        raise ConfigurationError(f"tracks must hold 2 to {_MAX_TRACKS} tracks.")
    if isinstance(links, (str, bytes)) or not isinstance(links, Sequence) or not links:
        raise ConfigurationError("links must be a non-empty list of link records.")
    for value, name in ((width_cm, "width_cm"), (track_height_cm, "track_height_cm")):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 3 <= value <= 200:
            raise ConfigurationError(f"{name} must be a number between 3 and 200.")
    target = cli.new_file(output_path, "output_path")
    parsed: list[dict[str, Any]] = []
    for index, track in enumerate(tracks):
        if not isinstance(track, Mapping):
            raise ConfigurationError("each track must be a mapping.")
        unknown = set(track) - {"name", "length", "features"}
        if unknown:
            raise ConfigurationError(f"tracks[{index}]: unknown keys {sorted(unknown)}.")
        length = _integer(
            track.get("length"), f"tracks[{index}].length", minimum=1, maximum=_MAX_LENGTH
        )
        features = track.get("features", [])
        if (
            isinstance(features, (str, bytes))
            or not isinstance(features, Sequence)
            or len(features) > _MAX_FEATURES
        ):
            raise ConfigurationError(f"tracks[{index}].features must be a list.")
        parsed.append(
            {
                "name": _label(track.get("name"), f"tracks[{index}].name"),
                "length": length,
                "features": [
                    _feature(item, length, f"tracks[{index}].features[{n}]")
                    for n, item in enumerate(features)
                ],
            }
        )
    if len({item["name"] for item in parsed}) != len(parsed):
        raise ConfigurationError("track names must be distinct.")

    def locate(value: object, where: str) -> int:
        if isinstance(value, int) and not isinstance(value, bool) and 0 <= value < len(parsed):
            return value
        for index, item in enumerate(parsed):
            if item["name"] == value:
                return index
        raise ConfigurationError(f"{where} must be a track name or index.")

    parsed_links = []
    for index, link in enumerate(links):
        if not isinstance(link, Mapping):
            raise ConfigurationError("each link must be a mapping.")
        unknown = set(link) - {
            "track_a",
            "start_a",
            "end_a",
            "track_b",
            "start_b",
            "end_b",
            "color",
            "flip",
        }
        if unknown:
            raise ConfigurationError(f"links[{index}]: unknown keys {sorted(unknown)}.")
        a, b = (
            locate(link.get("track_a"), f"links[{index}].track_a"),
            locate(link.get("track_b"), f"links[{index}].track_b"),
        )
        if a == b:
            raise ConfigurationError(f"links[{index}]: a link must join two different tracks.")
        flip = link.get("flip", False)
        if not isinstance(flip, bool):
            raise ConfigurationError(f"links[{index}].flip must be True or False.")
        start_a = _integer(
            link.get("start_a"), f"links[{index}].start_a", minimum=0, maximum=parsed[a]["length"]
        )
        end_a = _integer(
            link.get("end_a"),
            f"links[{index}].end_a",
            minimum=start_a + 1,
            maximum=parsed[a]["length"],
        )
        start_b = _integer(
            link.get("start_b"), f"links[{index}].start_b", minimum=0, maximum=parsed[b]["length"]
        )
        end_b = _integer(
            link.get("end_b"),
            f"links[{index}].end_b",
            minimum=start_b + 1,
            maximum=parsed[b]["length"],
        )
        parsed_links.append(
            {
                "track_a": parsed[a]["name"],
                "start_a": start_a,
                "end_a": end_a,
                "track_b": parsed[b]["name"],
                "start_b": start_b,
                "end_b": end_b,
                "color": link.get("color"),
                "flip": flip,
                "_a": a,
                "_b": b,
            }
        )
    diagram_title = _label(title, "title") if title is not None else None

    gd = require("Bio.Graphics.GenomeDiagram", _EXTRA)
    cm = require("reportlab.lib.units", _REPORT_EXTRA).cm
    seqfeature = require("Bio.SeqFeature", _EXTRA)
    longest = max(item["length"] for item in parsed)
    diagram = gd.Diagram(diagram_title or "genome diagram")
    objects = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for index, item in enumerate(parsed):
            # odd track levels leave a gap between bands for the cross-link ribbons
            track = diagram.new_track(
                2 * (len(parsed) - index) - 1,
                name=item["name"],
                greytrack=True,
                greytrack_labels=1,
                start=0,
                end=item["length"],
                height=1.0,
            )
            feature_set = track.new_set()
            for feature in item["features"]:
                feature_set.add_feature(
                    seqfeature.SeqFeature(
                        seqfeature.FeatureLocation(
                            feature["start"], feature["end"], strand=feature["strand"] or None
                        )
                    ),
                    color=_color(feature["color"], "feature color", "steelblue"),
                    label=feature["label"] is not None,
                    name=feature["label"] or "",
                    sigil="ARROW" if feature["strand"] else "BOX",
                )
            objects.append(track)
        for link in parsed_links:
            diagram.cross_track_links.append(
                gd.CrossLink(
                    (objects[link["_a"]], link["start_a"], link["end_a"]),
                    (objects[link["_b"]], link["start_b"], link["end_b"]),
                    _color(link["color"], "link color", "lightgreen"),
                    flip=link["flip"],
                )
            )
        try:
            diagram.draw(
                format="linear",
                pagesize=(width_cm * cm, (len(parsed) * 2 * track_height_cm / 2 + 2) * cm),
                fragments=1,
                start=0,
                end=longest,
            )
            diagram.write(str(target), fmt.upper())
        except Exception as exc:  # reportlab raises several exception types
            raise _write_error(exc, fmt, target) from exc
    return {
        "status": "completed",
        "operation": "genome_diagram_crosslinks",
        "output": str(target),
        "image_format": fmt,
        "bytes": target.stat().st_size,
        "tracks": [{key: item[key] for key in ("name", "length", "features")} for item in parsed],
        "links": [
            {key: value for key, value in link.items() if not key.startswith("_")}
            for link in parsed_links
        ],
        "provenance": _provenance(parser="Bio.Graphics.GenomeDiagram (Diagram, CrossLink)"),
    }


def chromosome_ideogram(
    chromosomes: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    title: str = "Chromosomes",
    image_format: str = "pdf",
    telomere_length: int | None = None,
    page_width_cm: float = 24.0,
    page_height_cm: float = 18.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw chromosomes side by side as an ideogram with telomere caps (``Bio.Graphics.BasicChromosome``).

    ``chromosomes`` is a list of ``{"name", "length", "features": [{"start", "end", "strand", "label",
    "color"}]}``; features are marked beside the chromosome body (strand +1 on the right, -1 on the left,
    0 across). All chromosomes use one scale (the longest chromosome). ``telomere_length`` (bases) sets
    the size of the end caps, default 2 percent of the longest chromosome. The figure is written to the new
    file ``output_path``; returns the path and the drawn data (``Organism.draw``, ``Chromosome``,
    ``TelomereSegment``, ``AnnotatedChromosomeSegment``).
    """

    fmt = _format(image_format)
    if (
        isinstance(chromosomes, (str, bytes))
        or not isinstance(chromosomes, Sequence)
        or not 1 <= len(chromosomes) <= 60
    ):
        raise ConfigurationError("chromosomes must hold 1 to 60 chromosomes.")
    for value, name in ((page_width_cm, "page_width_cm"), (page_height_cm, "page_height_cm")):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 5 <= value <= 200:
            raise ConfigurationError(f"{name} must be a number between 5 and 200.")
    heading = _label(title, "title")
    target = cli.new_file(output_path, "output_path")
    parsed: list[dict[str, Any]] = []
    for index, chromosome in enumerate(chromosomes):
        if not isinstance(chromosome, Mapping):
            raise ConfigurationError("each chromosome must be a mapping.")
        unknown = set(chromosome) - {"name", "length", "features"}
        if unknown:
            raise ConfigurationError(f"chromosomes[{index}]: unknown keys {sorted(unknown)}.")
        length = _integer(
            chromosome.get("length"), f"chromosomes[{index}].length", minimum=1, maximum=_MAX_LENGTH
        )
        features = chromosome.get("features", [])
        if (
            isinstance(features, (str, bytes))
            or not isinstance(features, Sequence)
            or len(features) > _MAX_FEATURES
        ):
            raise ConfigurationError(f"chromosomes[{index}].features must be a list.")
        parsed.append(
            {
                "name": _label(chromosome.get("name"), f"chromosomes[{index}].name"),
                "length": length,
                "features": [
                    _feature(item, length, f"chromosomes[{index}].features[{n}]")
                    for n, item in enumerate(features)
                ],
            }
        )
    longest = max(item["length"] for item in parsed)
    cap = (
        _integer(telomere_length, "telomere_length", minimum=1, maximum=_MAX_LENGTH)
        if telomere_length is not None
        else max(1, longest // 50)
    )
    basic = require("Bio.Graphics.BasicChromosome", _EXTRA)
    cm = require("reportlab.lib.units", _REPORT_EXTRA).cm
    organism = basic.Organism(output_format=fmt)
    organism.page_size = (page_width_cm * cm, page_height_cm * cm)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for item in parsed:
            chromosome = basic.Chromosome(item["name"])
            chromosome.scale_num = longest + 2 * cap
            start_cap = basic.TelomereSegment()
            start_cap.scale = cap
            chromosome.add(start_cap)
            if item["features"]:
                body = basic.AnnotatedChromosomeSegment(
                    item["length"],
                    [
                        (
                            feature["start"],
                            feature["end"],
                            feature["strand"],
                            feature["label"] or "",
                            _color(feature["color"], "feature color", "blue"),
                        )
                        for feature in item["features"]
                    ],
                )
            else:
                body = basic.ChromosomeSegment()
            body.scale = item["length"]
            chromosome.add(body)
            end_cap = basic.TelomereSegment(inverted=True)
            end_cap.scale = cap
            chromosome.add(end_cap)
            organism.add(chromosome)
        try:
            organism.draw(str(target), heading)
        except Exception as exc:  # reportlab raises several exception types
            raise _write_error(exc, fmt, target) from exc
    return {
        "status": "completed",
        "operation": "chromosome_ideogram",
        "output": str(target),
        "image_format": fmt,
        "bytes": target.stat().st_size,
        "telomere_length": cap,
        "chromosomes": parsed,
        "provenance": _provenance(
            parser="Bio.Graphics.BasicChromosome (Organism, Chromosome, segments)"
        ),
    }


__all__ = ["chromosome_ideogram", "genome_diagram_crosslinks"]

_register_public_verbose(globals(), tuple(__all__))
