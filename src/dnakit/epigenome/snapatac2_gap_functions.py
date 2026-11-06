"""snapatac2 registry downloads and plots as typed DNAKit functions.

The user installs snapatac2 (``pip install "dnakit[bio-snapatac2]"``). Download functions call the snapatac2
dataset and genome loaders through their module attributes and return the local paths, the parameters, a UTC
retrieval time and the parsed records. Plot functions run the matching snapatac2 plotting function with
``show=False`` and draw the plotted figure with Matplotlib into one new PNG file; they return the image path
and the plotted values. Plotly's PNG export needs a Chrome installation, so the figure is not exported by
Plotly itself. The coverage plot reads the same coverage tracks through snapatac2's internal ``get_coverage``
and draws them with Matplotlib, because the upstream function draws through pyplot's global state.

Inputs are checked before any upstream call. Missing optional packages raise ``BackendUnavailableError`` with
the install hint; download failures raise ``DownloadError``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError, DownloadError

_EXTRA = "bio-snapatac2"
_NUMPY_EXTRA = "bio-analysis"
_PLOT_EXTRA = "bio-plots"
_SERVICE = "snapatac2 dataset registry (pooch cache)"
_CACHE_HINT = (
    "Check the network connection. snapatac2 downloads from its registry and fallback URLs; "
    "set SNAP_DATA_DIR to choose the cache folder."
)
_MAX_DATA_VALUES = 20_000
_GENOMES: dict[str, tuple[str, str, str]] = {
    # canonical preset name: (registry dataset, FASTA file, annotation file), as in snapatac2.genome
    "GRCh37": ("gencode_GRCh37", "gencode_v41_GRCh37.fa.gz", "gencode_v41_GRCh37.gff3.gz"),
    "GRCh38": ("gencode_GRCh38", "gencode_v41_GRCh38.fa.gz", "gencode_v41_GRCh38.gff3.gz"),
    "GRCm39": ("gencode_GRCm39", "gencode_vM30_GRCm39.fa.gz", "gencode_vM30_GRCm39.gff3.gz"),
    "GRCm38": ("gencode_GRCm38", "gencode_vM25_GRCm38.fa.gz", "gencode_vM25_GRCm38.gff3.gz"),
}
_GENOME_ALIASES = {"hg19": "GRCh37", "hg38": "GRCh38", "mm39": "GRCm39", "mm10": "GRCm38"}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of: {', '.join(allowed)}.")
    return value


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _new_png(path: str | Path) -> Path:
    target = Path(path).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_png must end with .png.")
    if target.exists():
        raise ConfigurationError(f"output_png must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for output_png does not exist: {target.parent}")
    return target


def _existing_file(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _paths(value: object) -> list[str]:
    if isinstance(value, (str, Path)):
        return [str(value)]
    if isinstance(value, (list, tuple)) and all(isinstance(item, (str, Path)) for item in value):
        return [str(item) for item in value]
    return []


def _download(call: Callable[[], Any], *, dataset: str) -> Any:
    try:
        return call()
    except KeyError as exc:
        raise ConfigurationError(
            f"snapatac2 does not know the dataset or file {dataset!r}: {exc}"
        ) from exc
    except TypeError as exc:
        raise ConfigurationError(f"snapatac2 rejected the request for {dataset!r}: {exc}") from exc
    except (OSError, ValueError, RuntimeError) as exc:
        raise DownloadError(
            f"snapatac2 could not fetch {dataset!r}: {exc}", hint=_CACHE_HINT
        ) from exc


def _result(
    dataset: str,
    parameters: Mapping[str, Any],
    paths: Sequence[str],
    data: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "status": "completed",
        "service": _SERVICE,
        "dataset": dataset,
        "parameters": dict(parameters),
        "paths": list(paths),
        "retrieved_utc": _utc_now(),
        "data": dict(data or {}),
    }


def _motif_summary(motifs: Sequence[Any]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for motif in motifs:
        row: dict[str, Any] = {"id": str(motif.id), "name": str(motif.name)}
        family = getattr(motif, "family", None)
        if family is not None:
            row["family"] = str(family)
        rows.append(row)
    return {"motif_count": len(rows), "motifs": rows}


def _network_arrays(
    node_types: Sequence[str], edges: Sequence[Mapping[str, Any]]
) -> tuple[list[str], list[dict[str, Any]]]:
    """Validate a network given as node types and edges that refer to nodes by index."""
    if not node_types:
        raise ConfigurationError("node_types must contain at least one node.")
    types = [_text(item, "node type") for item in node_types]
    rows: list[dict[str, Any]] = []
    for position, edge in enumerate(edges):
        if not isinstance(edge, Mapping):
            raise ConfigurationError(f"edges[{position}] must be a mapping.")
        source, target = edge.get("source"), edge.get("target")
        for label, node in (("source", source), ("target", target)):
            if isinstance(node, bool) or not isinstance(node, int) or not 0 <= node < len(types):
                raise ConfigurationError(
                    f"edges[{position}].{label} must be a node index below {len(types)}."
                )
        row: dict[str, Any] = {"source": source, "target": target}
        for key in ("distance", "cor_score", "regr_score"):
            value = edge.get(key)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, (int, float))
            ):
                raise ConfigurationError(f"edges[{position}].{key} must be a number or absent.")
            row[key] = None if value is None else float(value)
        rows.append(row)
    if not rows:
        raise ConfigurationError("edges must contain at least one edge.")
    return types, rows


def _rustworkx_graph(types: Sequence[str], rows: Sequence[Mapping[str, Any]], rx: Any) -> Any:
    from types import SimpleNamespace

    graph = rx.PyDiGraph()
    indices = [graph.add_node(SimpleNamespace(type=node_type)) for node_type in types]
    for row in rows:
        payload = SimpleNamespace(
            distance=row["distance"], cor_score=row["cor_score"], regr_score=row["regr_score"]
        )
        graph.add_edge(indices[row["source"]], indices[row["target"]], payload)
    return graph


def load_snapatac2_dataset(
    dataset: str, *, file: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Fetch one file of a snapatac2 registry dataset into the local cache (snapatac2.datasets.load)."""
    name = _text(dataset, "dataset")
    if file is not None:
        file = _text(file, "file")
    datasets = require("snapatac2.datasets", _EXTRA)
    kwargs: dict[str, Any] = {} if file is None else {"file": file}
    result = _download(lambda: datasets.load(name, **kwargs), dataset=name)
    paths = _paths(result)
    data: dict[str, Any] = {}
    if not paths:
        data["object_type"] = type(result).__name__
        data["object_count"] = len(result) if isinstance(result, list) else None
    return _result(name, {"file": file}, paths, data)


def fetch_snapatac2_pbmc500(
    *, data_type: str = "fragment", downsample: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Download the 10x 500-cell PBMC scATAC-seq example as fragment, BAM or FASTQ files (snapatac2.datasets.pbmc500)."""
    kind = _choice(data_type, "data_type", ("fragment", "bam", "fastq"))
    if not isinstance(downsample, bool):
        raise ConfigurationError("downsample must be a boolean.")
    if downsample and kind != "fragment":
        raise ConfigurationError("downsample is only available with data_type='fragment'.")
    datasets = require("snapatac2.datasets", _EXTRA)
    result = _download(
        lambda: datasets.pbmc500(type=kind, downsample=downsample), dataset="atac_pbmc_500"
    )
    return _result("atac_pbmc_500", {"data_type": kind, "downsample": downsample}, _paths(result))


def fetch_snapatac2_pbmc5k(*, data_type: str = "fragment", verbose: bool = True) -> dict[str, Any]:
    """Download the 10x 5k PBMC scATAC-seq example as a fragment, h5ad or annotated h5ad file (snapatac2.datasets.pbmc5k)."""
    kind = _choice(data_type, "data_type", ("fragment", "h5ad", "annotated_h5ad"))
    datasets = require("snapatac2.datasets", _EXTRA)
    result = _download(lambda: datasets.pbmc5k(type=kind), dataset="atac_pbmc_5k")
    return _result("atac_pbmc_5k", {"data_type": kind}, _paths(result))


def fetch_snapatac2_pbmc10k_multiome(
    *,
    modality: str = "RNA",
    data_type: str = "h5ad",
    verbose: bool = True,
) -> dict[str, Any]:
    """Download the 10x 10k PBMC multiome RNA or ATAC example (snapatac2.datasets.pbmc10k_multiome)."""
    mode = _choice(modality, "modality", ("RNA", "ATAC"))
    kind = _choice(data_type, "data_type", ("fragment", "h5ad"))
    datasets = require("snapatac2.datasets", _EXTRA)
    result = _download(
        lambda: datasets.pbmc10k_multiome(modality=mode, type=kind), dataset="pbmc_10k_multiome"
    )
    return _result("pbmc_10k_multiome", {"modality": mode, "data_type": kind}, _paths(result))


def fetch_snapatac2_colon_transverse(*, verbose: bool = True) -> dict[str, Any]:
    """Download and extract the transverse colon scATAC-seq fragment files, one per sample (snapatac2.datasets.colon)."""
    datasets = require("snapatac2.datasets", _EXTRA)
    samples = _download(lambda: datasets.colon(), dataset="colon_transverse")
    rows = [{"sample": str(name), "fragment_file": str(path)} for name, path in samples]
    return _result(
        "colon_transverse", {}, [row["fragment_file"] for row in rows], {"samples": rows}
    )


def fetch_snapatac2_hea_ccre(*, verbose: bool = True) -> dict[str, Any]:
    """Download the curated human colon cis-regulatory element BED file (snapatac2.datasets.cre_HEA)."""
    datasets = require("snapatac2.datasets", _EXTRA)
    path = _download(lambda: datasets.cre_HEA(), dataset="HEA_cCRE")
    return _result("HEA_cCRE", {}, _paths(path))


def fetch_snapatac2_cis_bp_motifs(*, unique: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Download the CIS-BP transcription factor motifs and return their names and identifiers (snapatac2.datasets.cis_bp)."""
    if not isinstance(unique, bool):
        raise ConfigurationError("unique must be a boolean.")
    datasets = require("snapatac2.datasets", _EXTRA)
    motifs = _download(lambda: datasets.cis_bp(unique=unique), dataset="cisBP_human")
    data = {"registry_file": "cisBP_human.meme", **_motif_summary(motifs)}
    return _result("cisBP_human", {"unique": unique}, [], data)


def fetch_snapatac2_meuleman_2020_motifs(*, verbose: bool = True) -> dict[str, Any]:
    """Download the Meuleman 2020 grouped motifs and return their names, identifiers and families (snapatac2.datasets.Meuleman_2020)."""
    datasets = require("snapatac2.datasets", _EXTRA)
    motifs = _download(lambda: datasets.Meuleman_2020(), dataset="Meuleman_2020")
    data = {"registry_file": "Meuleman_2020.meme", **_motif_summary(motifs)}
    return _result("Meuleman_2020", {}, [], data)


def fetch_snapatac2_genome_files(
    genome: str, *, kind: str = "both", verbose: bool = True
) -> dict[str, Any]:
    """Download the FASTA and/or gene annotation of a genome preset through snapatac2.genome.load."""
    name = _text(genome, "genome")
    canonical = _GENOME_ALIASES.get(name, name)
    if canonical not in _GENOMES:
        raise ConfigurationError(
            f"genome must be one of: {', '.join(sorted(_GENOMES) + sorted(_GENOME_ALIASES))}."
        )
    choice = _choice(kind, "kind", ("fasta", "annotation", "both"))
    dataset, fasta_name, annotation_name = _GENOMES[canonical]
    genome_module = require("snapatac2.genome", _EXTRA)
    pooch = require("pooch", _EXTRA)
    paths: dict[str, str] = {}
    if choice in ("fasta", "both"):
        fasta = _download(
            lambda: genome_module.load(
                dataset, file=fasta_name, processor=pooch.Decompress(method="gzip")
            ),
            dataset=dataset,
        )
        paths["fasta"] = str(fasta)
    if choice in ("annotation", "both"):
        annotation = _download(
            lambda: genome_module.load(dataset, file=annotation_name), dataset=dataset
        )
        paths["annotation"] = str(annotation)
    return _result(
        dataset,
        {"genome": canonical, "requested_name": name, "kind": choice},
        list(paths.values()),
        paths,
    )


# ---------------------------------------------------------------------------
# Plots: upstream figures are drawn with Matplotlib (Plotly's PNG export needs Chrome, which is not installed).
# ---------------------------------------------------------------------------

_COLORS = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)
_LINKAGES = ("single", "complete", "average", "weighted", "centroid", "median", "ward")


def _numpy() -> Any:
    return require("numpy", _NUMPY_EXTRA)


def _plain(value: Any) -> Any:
    """Return a JSON-friendly copy: numpy values become Python values and non-finite floats become None."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        return _plain(value.item())
    return str(value)


def _decode(value: Any, np: Any) -> Any:
    """Decode Plotly's binary array encoding ({dtype, bdata, shape}) into numpy arrays."""
    if isinstance(value, Mapping):
        if "bdata" in value and "dtype" in value:
            import base64

            array = np.frombuffer(base64.b64decode(value["bdata"]), dtype=np.dtype(value["dtype"]))
            shape = value.get("shape")
            if shape:
                array = array.reshape(
                    tuple(int(part) for part in str(shape).split(",") if part.strip())
                )
            return array
        return {str(key): _decode(item, np) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_decode(item, np) for item in value]
    return value


def _spec_from_figure(figure: Any) -> dict[str, Any]:
    """Return a Plotly figure as {'data': [...], 'layout': {...}} with arrays decoded and the template removed."""
    np = _numpy()
    raw = figure.to_plotly_json()
    layout = {key: value for key, value in raw.get("layout", {}).items() if key != "template"}
    return {
        "data": [_decode(trace, np) for trace in raw.get("data", [])],
        "layout": _decode(layout, np),
    }


def _title(value: Any) -> str | None:
    if isinstance(value, Mapping):
        value = value.get("text")
    return value if isinstance(value, str) and value.strip() else None


def _axis_key(ref: str) -> str:
    return ref[0] + "axis" + ref[1:]


def _domain(layout: Mapping[str, Any], ref: str) -> tuple[float, float]:
    domain = (layout.get(_axis_key(ref)) or {}).get("domain") or [0.0, 1.0]
    return float(domain[0]), float(domain[1])


def _array(np: Any, values: Any, name: str) -> Any:
    array = np.asarray(values)
    if array.dtype.kind not in "biuf":
        raise ConfigurationError(
            f"The figure has non-numeric {name} values, which this renderer does not draw."
        )
    return array.astype(float)


def _coordinates(np: Any, values: Any, count: int) -> tuple[Any, list[str] | None]:
    """Numeric positions for a categorical or numeric axis, with category labels when the axis is categorical."""
    if values is None:
        return np.arange(count, dtype=float), None
    array = np.asarray(values)
    if array.shape[0] != count:
        raise ConfigurationError(
            "A figure coordinate does not match the number of matrix rows or columns."
        )
    if array.dtype.kind in "biuf":
        return array.astype(float), None
    return np.arange(count, dtype=float), [str(item) for item in array]


def _edges(np: Any, centers: Any) -> Any:
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5])
    middle = (centers[1:] + centers[:-1]) / 2
    first = centers[0] - (middle[0] - centers[0])
    last = centers[-1] + (centers[-1] - middle[-1])
    return np.concatenate([[first], middle, [last]])


def _css_color(text: str) -> Any:
    import re

    match = re.fullmatch(r"rgba?\(([\d.\s,]+)\)", text.strip())
    if match:
        parts = [float(part) for part in match.group(1).split(",")[:3]]
        return tuple(part / 255.0 for part in parts)
    return text


def _colormap(scale: Any) -> Any:
    colors = require("matplotlib.colors", _PLOT_EXTRA)
    maps = require("matplotlib", _PLOT_EXTRA).colormaps
    if isinstance(scale, str):
        for name in (scale, scale.lower()):
            if name in maps:
                return maps[name]
    if isinstance(scale, (list, tuple)) and scale:
        palette = [
            _css_color(str(pair[1]))
            for pair in scale
            if isinstance(pair, (list, tuple)) and len(pair) == 2
        ]
        if palette:
            return colors.LinearSegmentedColormap.from_list("plotly_scale", palette)
    return maps["viridis"]


def _legend_label(trace: Mapping[str, Any]) -> str | None:
    name = trace.get("name")
    if trace.get("showlegend") is False or not isinstance(name, str) or not name:
        return None
    return name


def _trace_color(trace: Mapping[str, Any], index: int) -> Any:
    marker = trace.get("marker") or {}
    if isinstance(marker.get("color"), str):
        return _css_color(marker["color"])
    line = trace.get("line") or {}
    if isinstance(line.get("color"), str):
        return _css_color(line["color"])
    return _COLORS[index % len(_COLORS)]


def _draw_scatter(
    ax: Any, trace: Mapping[str, Any], index: int, layout: Mapping[str, Any], np: Any
) -> tuple[Any, list[str] | None]:
    y = _array(np, trace.get("y"), "y")
    x, x_labels = _coordinates(np, trace.get("x"), y.size)
    mode = str(trace.get("mode") or "lines")
    marker = trace.get("marker") or {}
    size = float(marker.get("size") or 6)
    opacity = float(marker.get("opacity") if marker.get("opacity") is not None else 1.0)
    label = _legend_label(trace)
    if "lines" in mode:
        ax.plot(x, y, color=_trace_color(trace, index), linewidth=1.0, alpha=opacity, label=label)
    if "markers" not in mode:
        return None, x_labels
    color = marker.get("color")
    if isinstance(color, np.ndarray) and color.dtype.kind in "biuf":
        scale = (layout.get("coloraxis") or {}).get("colorscale", "viridis")
        mappable = ax.scatter(
            x,
            y,
            s=size**2,
            c=color.astype(float),
            cmap=_colormap(scale),
            alpha=opacity,
            linewidths=0,
        )
        return mappable, x_labels
    ax.scatter(
        x, y, s=size**2, color=_trace_color(trace, index), alpha=opacity, label=label, linewidths=0
    )
    return None, x_labels


def _draw_bar(ax: Any, trace: Mapping[str, Any], index: int, np: Any) -> list[str] | None:
    height = _array(np, trace.get("y"), "y")
    positions, labels = _coordinates(np, trace.get("x"), height.size)
    ax.bar(
        positions, height, width=0.8, color=_trace_color(trace, index), label=_legend_label(trace)
    )
    return labels


def _draw_violins(ax: Any, traces: Sequence[tuple[int, Mapping[str, Any]]], np: Any) -> None:
    positions: list[float] = []
    labels: list[str] = []
    datasets: list[Any] = []
    colors: list[str] = []
    for slot, (index, trace) in enumerate(traces):
        values = _array(np, trace.get("y"), "y").ravel()
        values = values[np.isfinite(values)]
        if values.size == 0:
            continue
        name = str(trace.get("name") or f"group {slot + 1}")
        if values.size == 1 or float(np.ptp(values)) == 0.0:
            ax.scatter([slot], [float(values[0])], color=_trace_color(trace, index), s=20)
        else:
            datasets.append(values)
            positions.append(slot)
            colors.append(_trace_color(trace, index))
        labels.append(name)
    if datasets:
        parts = ax.violinplot(datasets, positions=positions, showmedians=True, widths=0.8)
        for body, color in zip(parts["bodies"], colors, strict=False):
            body.set_facecolor(color)
            body.set_alpha(0.6)
    ax.set_xticks(list(range(len(traces))))
    ax.set_xticklabels(
        labels + [""] * (len(traces) - len(labels)), rotation=45, ha="right", fontsize=6
    )


def _draw_field(
    ax: Any, trace: Mapping[str, Any], kind: str, layout: Mapping[str, Any], np: Any
) -> tuple[Any, Any, Any]:
    z = np.asarray(trace.get("z"), dtype=float)
    if z.ndim != 2:
        raise ConfigurationError("A heatmap or contour trace needs a two-dimensional z matrix.")
    rows, columns = z.shape
    xs, x_labels = _coordinates(np, trace.get("x"), columns)
    ys, y_labels = _coordinates(np, trace.get("y"), rows)
    cmap = _colormap(trace.get("colorscale", "viridis"))
    if kind == "heatmap":
        mappable = ax.pcolormesh(_edges(np, xs), _edges(np, ys), z, cmap=cmap, shading="flat")
    else:
        mappable = ax.contourf(xs, ys, z, levels=12, cmap=cmap)
    return mappable, x_labels, y_labels


def _draw_shapes(ax: Any, layout: Mapping[str, Any], ref_x: str, np: Any) -> None:
    for shape in layout.get("shapes") or []:
        if shape.get("type") != "line" or shape.get("xref", "x") != ref_x:
            continue
        if shape.get("x0") is not None and shape.get("x0") == shape.get("x1"):
            ax.axvline(float(shape["x0"]), color="0.4", linestyle="--", linewidth=0.8)


def _style_axes(
    ax: Any, layout: Mapping[str, Any], ref_x: str, ref_y: str, labels: tuple[Any, Any], field: bool
) -> None:
    x_axis = layout.get(_axis_key(ref_x)) or {}
    y_axis = layout.get(_axis_key(ref_y)) or {}
    if x_axis.get("type") == "log":
        ax.set_xscale("log")
    if y_axis.get("type") == "log":
        ax.set_yscale("log")
    if _title(x_axis.get("title")):
        ax.set_xlabel(_title(x_axis.get("title")), fontsize=8)
    if _title(y_axis.get("title")):
        ax.set_ylabel(_title(y_axis.get("title")), fontsize=8)
    if y_axis.get("side") == "right":
        ax.yaxis.tick_right()
    if y_axis.get("autorange") == "reversed":
        ax.invert_yaxis()
    x_labels, y_labels = labels
    if x_axis.get("showticklabels") is False:
        ax.tick_params(axis="x", labelbottom=False)
    if y_axis.get("showticklabels") is False:
        ax.tick_params(axis="y", labelleft=False, labelright=False)
    for axis, axis_layout, categories in (
        (ax.xaxis, x_axis, x_labels),
        (ax.yaxis, y_axis, y_labels),
    ):
        if axis_layout.get("visible") is False:
            axis.set_ticks([])
            continue
        if (
            field
            and axis_layout.get("tickvals") is not None
            and axis_layout.get("ticktext") is not None
        ):
            values = [float(item) for item in axis_layout["tickvals"]]
            text = [str(item) for item in axis_layout["ticktext"]]
        elif categories is not None:
            values = list(range(len(categories)))
            text = list(categories)
        else:
            continue
        if axis_layout.get("showticklabels") is False:
            continue
        rotation = 45 if len(text) > 6 else 0
        if axis is ax.xaxis:
            ax.set_xticks(values)
            ax.set_xticklabels(
                text, rotation=rotation, ha="right" if rotation else "center", fontsize=7
            )
        else:
            ax.set_yticks(values)
            ax.set_yticklabels(text, fontsize=7)


def _draw_group(
    figure: Any,
    ax: Any,
    members: Sequence[tuple[int, Mapping[str, Any]]],
    layout: Mapping[str, Any],
    ref_x: str,
    ref_y: str,
    np: Any,
) -> list[tuple[Any, str | None]]:
    labels: tuple[Any, Any] = (None, None)
    colorbars: list[tuple[Any, str | None]] = []
    violins: list[tuple[int, Mapping[str, Any]]] = []
    field = False
    for index, trace in members:
        kind = str(trace.get("type", "scatter"))
        if kind == "scatter":
            mappable, x_labels = _draw_scatter(ax, trace, index, layout, np)
            if x_labels is not None:
                labels = (x_labels, labels[1])
            if mappable is not None:
                scale_title = _title(
                    ((layout.get("coloraxis") or {}).get("colorbar") or {}).get("title")
                )
                colorbars.append((mappable, scale_title))
        elif kind == "bar":
            labels = (_draw_bar(ax, trace, index, np), labels[1])
        elif kind == "violin":
            violins.append((index, trace))
        elif kind in ("heatmap", "contour"):
            field = True
            mappable, x_labels, y_labels = _draw_field(ax, trace, kind, layout, np)
            labels = (x_labels or labels[0], y_labels or labels[1])
            colorbar = trace.get("colorbar") or {}
            colorbars.append((mappable, _title(colorbar.get("title"))))
        else:
            raise ConfigurationError(f"The matplotlib renderer does not draw {kind} traces.")
    if violins:
        _draw_violins(ax, violins, np)
    _draw_shapes(ax, layout, ref_x, np)
    _style_axes(ax, layout, ref_x, ref_y, labels, field)
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=6, frameon=False)
    return colorbars


def _draw_3d(
    figure: Any, traces: Sequence[Mapping[str, Any]], layout: Mapping[str, Any], np: Any
) -> None:
    ax = figure.add_subplot(111, projection="3d")
    scene = layout.get("scene") or {}
    for index, trace in enumerate(traces):
        if trace.get("type") != "scatter3d":
            raise ConfigurationError(
                "A figure with scatter3d traces can only contain scatter3d traces."
            )
        marker = trace.get("marker") or {}
        size = float(marker.get("size") or 6)
        ax.scatter(
            _array(np, trace.get("x"), "x"),
            _array(np, trace.get("y"), "y"),
            _array(np, trace.get("z"), "z"),
            s=size**2,
            color=_trace_color(trace, index),
            alpha=float(marker.get("opacity") if marker.get("opacity") is not None else 1.0),
            label=_legend_label(trace),
        )
    for name, setter in (
        ("xaxis", ax.set_xlabel),
        ("yaxis", ax.set_ylabel),
        ("zaxis", ax.set_zlabel),
    ):
        title = _title((scene.get(name) or {}).get("title"))
        if title:
            setter(title, fontsize=8)
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=6, frameon=False)


def _draw_figure(spec: Mapping[str, Any], target: Path) -> None:
    """Draw a Plotly-style figure dictionary with Matplotlib and write it to a new PNG file."""
    np = _numpy()
    figure_module = require("matplotlib.figure", _PLOT_EXTRA)
    canvas_module = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    layout = spec.get("layout") or {}
    traces = list(spec.get("data") or [])
    if not traces:
        raise ConfigurationError("The figure has no traces to draw.")
    width = int(layout.get("width") or 600)
    height = int(layout.get("height") or 400)
    figure = figure_module.Figure(figsize=(width / 100.0, height / 100.0), dpi=100)
    canvas_module.FigureCanvasAgg(figure)
    if any(trace.get("type") == "scatter3d" for trace in traces):
        _draw_3d(figure, traces, layout, np)
    else:
        # A colour bar takes a narrow column on the right; the plot area stops before it.
        has_colorbar = "coloraxis" in layout or any(
            trace.get("type") in ("heatmap", "contour") for trace in traces
        )
        left, right = 0.10, (0.80 if has_colorbar else 0.95)
        groups: dict[tuple[str, str], list[tuple[int, Mapping[str, Any]]]] = {}
        for index, trace in enumerate(traces):
            key = (str(trace.get("xaxis") or "x"), str(trace.get("yaxis") or "y"))
            groups.setdefault(key, []).append((index, trace))
        colorbars: list[tuple[Any, str | None]] = []
        for (ref_x, ref_y), members in groups.items():
            x0, x1 = _domain(layout, ref_x)
            y0, y1 = _domain(layout, ref_y)
            ax = figure.add_axes(
                [
                    left + (right - left) * x0,
                    0.12 + 0.72 * y0,
                    (right - left) * (x1 - x0),
                    0.72 * (y1 - y0),
                ]
            )
            colorbars.extend(_draw_group(figure, ax, members, layout, ref_x, ref_y, np))
        if colorbars:
            colorbar_axis = figure.add_axes([0.925, 0.12, 0.018, 0.72]) if has_colorbar else None
            for mappable, title in colorbars:
                if colorbar_axis is None:
                    bar = figure.colorbar(mappable)
                else:
                    bar = figure.colorbar(mappable, cax=colorbar_axis)
                if title and colorbar_axis is not None:
                    colorbar_axis.set_title(title, fontsize=7)
                elif title:
                    bar.set_label(title, fontsize=7)
    title = _title(layout.get("title"))
    if title:
        figure.suptitle(title, fontsize=9)
    with require("matplotlib", _PLOT_EXTRA).rc_context({"savefig.bbox": "standard"}):
        figure.savefig(str(target), format="png", dpi=100)


def _bounded(np: Any, values: Any) -> Any:
    array = np.asarray(values)
    if array.size > _MAX_DATA_VALUES:
        return {"omitted": True, "length": int(array.size)}
    return _plain(array.tolist())


def _figure_data(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Summarise the values a figure plots: every trace array (up to a size limit) and the axis labels."""
    np = _numpy()
    layout = spec.get("layout") or {}
    traces: list[dict[str, Any]] = []
    for trace in spec.get("data") or []:
        entry: dict[str, Any] = {
            "type": str(trace.get("type", "scatter")),
            "name": _plain(trace.get("name")),
        }
        for key in ("x", "y", "z"):
            if trace.get(key) is not None:
                entry[key] = _bounded(np, trace[key])
        traces.append(entry)
    data: dict[str, Any] = {"traces": traces, "title": _title(layout.get("title"))}
    for ref in ("x", "y"):
        axis = layout.get(_axis_key(ref)) or {}
        data[f"{ref}axis_title"] = _title(axis.get("title"))
        data[f"{ref}axis_type"] = axis.get("type", "linear")
        if axis.get("ticktext") is not None:
            data[f"{ref}axis_ticktext"] = _plain(list(axis["ticktext"]))
    return data


def _plot_result(
    operation: str,
    spec: Mapping[str, Any],
    target: Path,
    parameters: Mapping[str, Any],
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    _draw_figure(spec, target)
    data = _figure_data(spec)
    if extra:
        data.update(_plain(dict(extra)))
    return {
        "status": "completed",
        "operation": operation,
        "image_path": str(target),
        "image_bytes": target.stat().st_size,
        "parameters": _plain(dict(parameters)),
        "data": data,
    }


def _values(values: Any, name: str, *, minimum: int = 1, positive: bool = False) -> Any:
    np = _numpy()
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a sequence of numbers.") from exc
    if array.ndim != 1 or array.size < minimum:
        raise ConfigurationError(
            f"{name} must be a one-dimensional sequence with at least {minimum} values."
        )
    if not np.all(np.isfinite(array)):
        raise ConfigurationError(f"{name} must not contain NaN or infinite values.")
    if positive and np.any(array <= 0):
        raise ConfigurationError(f"{name} must be positive.")
    return array


def _labels(values: Any, count: int, name: str) -> list[Any]:
    items = list(values)
    if len(items) != count:
        raise ConfigurationError(f"{name} must have one value per point ({count}).")
    if all(isinstance(item, (int, float)) and not isinstance(item, bool) for item in items):
        return [float(item) for item in items]
    return [str(item) for item in items]


def _names(values: Sequence[str] | None, count: int, name: str) -> list[str] | None:
    if values is None:
        return None
    items = [_text(item, name) for item in values]
    if len(items) != count:
        raise ConfigurationError(f"{name} must have {count} names.")
    return items


def _matrix(values: Any, name: str) -> Any:
    np = _numpy()
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular matrix of numbers.") from exc
    if array.ndim != 2 or min(array.shape) < 1:
        raise ConfigurationError(f"{name} must be a non-empty two-dimensional matrix.")
    if not np.all(np.isfinite(array)):
        raise ConfigurationError(f"{name} must not contain NaN or infinite values.")
    return array


def _network_graph(
    node_types: Sequence[str], edges: Sequence[Mapping[str, Any]]
) -> tuple[Any, list[dict[str, Any]]]:
    types, rows = _network_arrays(node_types, edges)
    rx = require("rustworkx", _EXTRA)
    return _rustworkx_graph(types, rows, rx), rows


def render_snapatac2_figure_png(
    figure: Mapping[str, Any],
    output_png: str,
    *,
    width: int = 600,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Render a Plotly figure dictionary with snapatac2's render_plot rules and write it to a new PNG file."""
    target = _new_png(output_png)
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    if not isinstance(figure, Mapping) or not isinstance(figure.get("data"), Sequence):
        raise ConfigurationError("figure must be a Plotly figure dictionary with a 'data' list.")
    go = require("plotly.graph_objects", _EXTRA)
    plotting = require("snapatac2.plotting", _EXTRA)
    try:
        plotly_figure = go.Figure(
            {"data": list(figure["data"]), "layout": dict(figure.get("layout") or {})}
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"figure is not a valid Plotly figure: {exc}") from exc
    rendered = plotting.render_plot(
        plotly_figure, width=width_px, height=height_px, interactive=True, show=False, out_file=None
    )
    return _plot_result(
        "render_plot", _spec_from_figure(rendered), target, {"width": width_px, "height": height_px}
    )


def plot_snapatac2_heatmap_png(
    matrix: Sequence[Sequence[float]],
    output_png: str,
    *,
    row_names: Sequence[str] | None = None,
    column_names: Sequence[str] | None = None,
    cluster_rows: bool = True,
    cluster_columns: bool = True,
    colorscale: str = "Blues",
    linkage: str = "ward",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a numeric matrix as a heatmap with optional hierarchical clustering dendrograms (snapatac2.plotting.heatmap)."""
    target = _new_png(output_png)
    values = _matrix(matrix, "matrix")
    rows, columns = values.shape
    if not isinstance(cluster_rows, bool) or not isinstance(cluster_columns, bool):
        raise ConfigurationError("cluster_rows and cluster_columns must be booleans.")
    if cluster_rows and rows < 2:
        raise ConfigurationError("Clustering rows needs at least two rows.")
    if cluster_columns and columns < 2:
        raise ConfigurationError("Clustering columns needs at least two columns.")
    names_rows = _names(row_names, rows, "row_names")
    names_columns = _names(column_names, columns, "column_names")
    method = _choice(linkage, "linkage", _LINKAGES)
    scale = _text(colorscale, "colorscale")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.heatmap(
        values,
        row_names=names_rows,
        column_names=names_columns,
        cluster_columns=cluster_columns,
        cluster_rows=cluster_rows,
        colorscale=scale,
        linkage=method,
        show=False,
    )
    if names_rows is not None and not cluster_rows:
        figure.update_layout(yaxis={"tickvals": list(range(rows)), "ticktext": names_rows})
    if names_columns is not None and not cluster_columns:
        figure.update_layout(xaxis={"tickvals": list(range(columns)), "ticktext": names_columns})
    return _plot_result(
        "heatmap",
        _spec_from_figure(figure),
        target,
        {
            "rows": rows,
            "columns": columns,
            "cluster_rows": cluster_rows,
            "cluster_columns": cluster_columns,
            "colorscale": scale,
            "linkage": method,
        },
    )


def plot_snapatac2_kde2d_png(
    x: Sequence[float],
    y: Sequence[float],
    output_png: str,
    *,
    log_x: bool = False,
    log_y: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a two-dimensional kernel density contour of paired observations (snapatac2.plotting.kde2d)."""
    target = _new_png(output_png)
    if not isinstance(log_x, bool) or not isinstance(log_y, bool):
        raise ConfigurationError("log_x and log_y must be booleans.")
    xs = _values(x, "x", minimum=3, positive=log_x)
    ys = _values(y, "y", minimum=3, positive=log_y)
    if xs.size != ys.size:
        raise ConfigurationError("x and y must have the same length.")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.kde2d(xs, ys, log_x=log_x, log_y=log_y)
    return _plot_result(
        "kde2d",
        _spec_from_figure(figure),
        target,
        {"points": int(xs.size), "log_x": log_x, "log_y": log_y},
    )


def plot_snapatac2_scatter_png(
    x: Sequence[float],
    y: Sequence[float],
    output_png: str,
    *,
    color: Sequence[Any] | None = None,
    x_label: str = "X",
    y_label: str = "Y",
    color_label: str = "Color",
    marker_size: float = 2.0,
    marker_opacity: float = 0.5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a two-dimensional scatter plot with optional categorical or numeric colours (snapatac2.plotting.scatter)."""
    target = _new_png(output_png)
    xs = _values(x, "x")
    ys = _values(y, "y")
    if xs.size != ys.size:
        raise ConfigurationError("x and y must have the same length.")
    colors = None if color is None else _labels(color, int(xs.size), "color")
    if marker_size <= 0 or not 0 <= marker_opacity <= 1:
        raise ConfigurationError(
            "marker_size must be positive and marker_opacity must be between 0 and 1."
        )
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.scatter(
        xs.tolist(),
        ys.tolist(),
        color=colors,
        x_label=_text(x_label, "x_label"),
        y_label=_text(y_label, "y_label"),
        color_label=_text(color_label, "color_label"),
        marker_size=float(marker_size),
        marker_opacity=float(marker_opacity),
        show=False,
    )
    return _plot_result(
        "scatter",
        _spec_from_figure(figure),
        target,
        {
            "points": int(xs.size),
            "marker_size": float(marker_size),
            "marker_opacity": float(marker_opacity),
        },
    )


def plot_snapatac2_scatter3d_png(
    x: Sequence[float],
    y: Sequence[float],
    z: Sequence[float],
    output_png: str,
    *,
    color: Sequence[Any] | None = None,
    x_label: str = "X",
    y_label: str = "Y",
    z_label: str = "Z",
    color_label: str = "Color",
    marker_size: float = 2.0,
    marker_opacity: float = 0.5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a three-dimensional scatter plot with optional categorical or numeric colours (snapatac2.plotting.scatter3d)."""
    target = _new_png(output_png)
    xs = _values(x, "x")
    ys = _values(y, "y")
    zs = _values(z, "z")
    if not xs.size == ys.size == zs.size:
        raise ConfigurationError("x, y and z must have the same length.")
    colors = None if color is None else _labels(color, int(xs.size), "color")
    if marker_size <= 0 or not 0 <= marker_opacity <= 1:
        raise ConfigurationError(
            "marker_size must be positive and marker_opacity must be between 0 and 1."
        )
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.scatter3d(
        xs.tolist(),
        ys.tolist(),
        zs.tolist(),
        color=colors,
        x_label=_text(x_label, "x_label"),
        y_label=_text(y_label, "y_label"),
        z_label=_text(z_label, "z_label"),
        color_label=_text(color_label, "color_label"),
        marker_size=float(marker_size),
        marker_opacity=float(marker_opacity),
        show=False,
    )
    return _plot_result(
        "scatter3d",
        _spec_from_figure(figure),
        target,
        {
            "points": int(xs.size),
            "marker_size": float(marker_size),
            "marker_opacity": float(marker_opacity),
        },
    )


def plot_snapatac2_network_scores_png(
    node_types: Sequence[str],
    edges: Sequence[Mapping[str, Any]],
    output_png: str,
    *,
    score_name: str,
    width: int = 800,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the mean edge score of a regulatory network in distance-to-TSS bins (snapatac2.plotting.network_scores)."""
    target = _new_png(output_png)
    score = _choice(score_name, "score_name", ("cor_score", "regr_score"))
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    graph, rows = _network_graph(node_types, edges)
    if any(row["distance"] is None for row in rows):
        raise ConfigurationError("Every edge needs a distance for network_scores.")
    if all(row[score] is None for row in rows):
        raise ConfigurationError(f"No edge has a {score} value.")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.network_scores(
        graph, score_name=score, width=width_px, height=height_px, show=False
    )
    return _plot_result(
        "network_scores",
        _spec_from_figure(figure),
        target,
        {"score_name": score, "edges": len(rows)},
    )


def plot_snapatac2_network_edge_stat_png(
    node_types: Sequence[str],
    edges: Sequence[Mapping[str, Any]],
    output_png: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw correlation-score violins for each source and target node type pair (snapatac2.plotting.network_edge_stat)."""
    target = _new_png(output_png)
    graph, rows = _network_graph(node_types, edges)
    if all(row["cor_score"] is None for row in rows):
        raise ConfigurationError("No edge has a cor_score value.")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.network_edge_stat(graph, show=False)
    return _plot_result(
        "network_edge_stat", _spec_from_figure(figure), target, {"edges": len(rows)}
    )


def plot_snapatac2_valid_cells_png(
    counts: Sequence[float],
    output_png: str,
    *,
    width: int = 500,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot ranked barcode counts on log-log axes to choose a cell-calling cutoff (snapatac2.plotting.valid_cells)."""
    target = _new_png(output_png)
    values = _values(counts, "counts", minimum=2, positive=True)
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.valid_cells(
        [float(item) for item in values], width=width_px, height=height_px, show=False
    )
    return _plot_result(
        "valid_cells", _spec_from_figure(figure), target, {"barcodes": int(values.size)}
    )


def plot_snapatac2_tsse_png(
    n_fragment: Sequence[float],
    tsse: Sequence[float],
    output_png: str,
    *,
    min_fragment: int = 500,
    width: int = 500,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot TSS enrichment against the number of unique fragments per cell (snapatac2.plotting.tsse)."""
    target = _new_png(output_png)
    fragments = _values(n_fragment, "n_fragment", minimum=2)
    scores = _values(tsse, "tsse", minimum=2)
    if fragments.size != scores.size:
        raise ConfigurationError("n_fragment and tsse must have one value per cell.")
    if isinstance(min_fragment, bool) or not isinstance(min_fragment, int) or min_fragment < 0:
        raise ConfigurationError("min_fragment must be a non-negative integer.")
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    anndata = require("anndata", _EXTRA)
    pandas = require("pandas", _EXTRA)
    obs = pandas.DataFrame(
        {"n_fragment": fragments, "tsse": scores},
        index=[f"cell{i}" for i in range(fragments.size)],
    )
    adata = anndata.AnnData(obs=obs)
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.tsse(
        adata, min_fragment=min_fragment, width=width_px, height=height_px, show=False
    )
    selected = int((fragments >= min_fragment).sum())
    return _plot_result(
        "tsse",
        _spec_from_figure(figure),
        target,
        {"cells": int(fragments.size), "min_fragment": min_fragment},
        {"cells_plotted": selected},
    )


def plot_snapatac2_frag_size_distr_png(
    distribution: Sequence[float],
    output_png: str,
    *,
    max_recorded_size: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the fragment size distribution from counts per fragment length (snapatac2.plotting.frag_size_distr)."""
    target = _new_png(output_png)
    counts = _values(distribution, "distribution", minimum=2)
    if bool(_numpy().any(counts < 0)):
        raise ConfigurationError("distribution counts must not be negative.")
    size = _positive_int(max_recorded_size, "max_recorded_size")
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.frag_size_distr(counts, max_recorded_size=size, show=False)
    return _plot_result(
        "frag_size_distr",
        _spec_from_figure(figure),
        target,
        {"max_recorded_size": size, "sizes": int(counts.size)},
    )


def plot_snapatac2_spectral_eigenvalues_png(
    eigenvalues: Sequence[float],
    output_png: str,
    *,
    width: int = 600,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot spectral embedding eigenvalues and the elbow that snapatac2 picks (snapatac2.plotting.spectral_eigenvalues)."""
    target = _new_png(output_png)
    values = _values(eigenvalues, "eigenvalues", minimum=2)
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    anndata = require("anndata", _EXTRA)
    adata = anndata.AnnData()
    adata.uns["spectral_eigenvalue"] = values
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.spectral_eigenvalues(adata, width=width_px, height=height_px, show=False)
    elbow = adata.uns.get("num_eigen")
    return _plot_result(
        "spectral_eigenvalues",
        _spec_from_figure(figure),
        target,
        {"eigenvalues": int(values.size)},
        {"num_eigen": None if elbow is None else int(elbow)},
    )


def plot_snapatac2_regions_png(
    matrix: Sequence[Sequence[float]],
    peak_names: Sequence[str],
    cell_groups: Sequence[str],
    peaks: Mapping[str, Sequence[str]],
    output_png: str,
    *,
    groupby: str = "group",
    width: int = 600,
    height: int = 400,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the grouped, RPKM-normalised accessibility of selected peaks (snapatac2.plotting.regions)."""
    import re

    target = _new_png(output_png)
    counts = _matrix(matrix, "matrix")
    cells, regions = counts.shape
    names = _names(peak_names, regions, "peak_names")
    if len(set(names or [])) != regions:
        raise ConfigurationError("peak_names must be unique.")
    if any(not re.fullmatch(r"[^:\s]+:\d+-\d+", name) for name in names or []):
        raise ConfigurationError("peak_names must be formatted as chrom:start-end.")
    groups = _names(cell_groups, cells, "cell_groups")
    label = _text(groupby, "groupby")
    if not isinstance(peaks, Mapping) or not peaks:
        raise ConfigurationError("peaks must map at least one group name to a list of peak names.")
    selection: dict[str, list[str]] = {}
    for group, members in peaks.items():
        chosen = _names(members, len(members), f"peaks[{group!r}]") or []
        unknown = [item for item in chosen if item not in set(names or [])]
        if not chosen or unknown:
            raise ConfigurationError(
                f"peaks[{group!r}] must list peaks from peak_names (unknown: {unknown[:3]})."
            )
        selection[str(group)] = chosen
    width_px = _positive_int(width, "width")
    height_px = _positive_int(height, "height")
    anndata = require("anndata", _EXTRA)
    pandas = require("pandas", _EXTRA)
    adata = anndata.AnnData(
        X=counts,
        obs=pandas.DataFrame({label: groups}, index=[f"cell{i}" for i in range(cells)]),
        var=pandas.DataFrame(index=names),
    )
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.regions(
        adata, groupby=label, peaks=selection, width=width_px, height=height_px, show=False
    )
    return _plot_result(
        "regions",
        _spec_from_figure(figure),
        target,
        {"cells": cells, "peaks": regions, "groupby": label},
        {"selected_peaks": sum(len(items) for items in selection.values())},
    )


def plot_snapatac2_umap_png(
    embedding: Sequence[Sequence[float]],
    output_png: str,
    *,
    color: Sequence[Any] | None = None,
    marker_size: float | None = None,
    marker_opacity: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a two- or three-dimensional embedding such as UMAP, coloured by group or value (snapatac2.plotting.umap)."""
    target = _new_png(output_png)
    coordinates = _matrix(embedding, "embedding")
    if coordinates.shape[1] not in (2, 3):
        raise ConfigurationError("embedding must have two or three columns.")
    if marker_size is not None and marker_size <= 0:
        raise ConfigurationError("marker_size must be positive.")
    if not 0 <= marker_opacity <= 1:
        raise ConfigurationError("marker_opacity must be between 0 and 1.")
    labels = (
        None if color is None else _numpy().asarray(_labels(color, coordinates.shape[0], "color"))
    )
    if labels is None and coordinates.shape[1] == 3:
        # snapatac2 2.10 fails for a 3D embedding without colours (scatter3d names a missing 'color' column),
        # so every point is put in one group; the legend then reads 'all cells'.
        labels = _numpy().asarray(["all cells"] * coordinates.shape[0])
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.umap(
        coordinates,
        color=labels,
        marker_size=None if marker_size is None else float(marker_size),
        marker_opacity=float(marker_opacity),
        show=False,
    )
    return _plot_result(
        "umap",
        _spec_from_figure(figure),
        target,
        {"points": int(coordinates.shape[0]), "dimensions": int(coordinates.shape[1])},
    )


def plot_snapatac2_motif_enrichment_png(
    enrichment: Mapping[str, Mapping[str, Sequence[Any]]],
    output_png: str,
    *,
    min_log_fc: float = 1.0,
    max_fdr: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw motif enrichment scores of several groups as a clustered heatmap (snapatac2.plotting.motif_enrichment)."""
    target = _new_png(output_png)
    if not isinstance(enrichment, Mapping) or len(enrichment) < 2:
        raise ConfigurationError("enrichment must map at least two group names to motif tables.")
    if min_log_fc < 0 or not 0 < max_fdr <= 1:
        raise ConfigurationError("min_log_fc must not be negative and max_fdr must be in (0, 1].")
    polars = require("polars", _EXTRA)
    required = ("id", "log2(fold change)", "adjusted p-value", "p-value")
    tables: dict[str, Any] = {}
    reference_ids: list[str] | None = None
    for group, columns in enrichment.items():
        if not isinstance(columns, Mapping) or any(name not in columns for name in required):
            raise ConfigurationError(
                f"enrichment[{group!r}] must be a table with the columns {', '.join(required)}."
            )
        if len({len(list(columns[name])) for name in required}) != 1:
            raise ConfigurationError(f"enrichment[{group!r}] columns must have the same length.")
        ids = [str(item) for item in columns["id"]]
        if reference_ids is None:
            reference_ids = ids
        elif ids != reference_ids:
            raise ConfigurationError("Every group must list the same motif ids in the same order.")
        tables[str(group)] = polars.DataFrame(
            {str(name): list(values) for name, values in columns.items()}
        )
    plotting = require("snapatac2.plotting", _EXTRA)
    figure = plotting.motif_enrichment(
        tables, min_log_fc=float(min_log_fc), max_fdr=float(max_fdr), show=False
    )
    return _plot_result(
        "motif_enrichment",
        _spec_from_figure(figure),
        target,
        {
            "groups": len(tables),
            "motifs": len(reference_ids or []),
            "min_log_fc": float(min_log_fc),
            "max_fdr": float(max_fdr),
        },
    )


def plot_snapatac2_coverage_png(
    input_h5ad: str,
    output_png: str,
    *,
    region: str,
    groupby: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw fragment coverage tracks of grouped cells across one region of an h5ad file (snapatac2 coverage data)."""
    import re

    target = _new_png(output_png)
    source = _existing_file(input_h5ad, "input_h5ad")
    text = _text(region, "region")
    match = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", text)
    if match is None or int(match.group(2)) >= int(match.group(3)):
        raise ConfigurationError(
            "region must be formatted as chrom:start-end with start below end."
        )
    start, end = int(match.group(2)), int(match.group(3))
    label = _text(groupby, "groupby")
    snapatac2 = require("snapatac2", _EXTRA)
    internal = require("snapatac2._snapatac2", _EXTRA)
    np = _numpy()
    figure_module = require("matplotlib.figure", _PLOT_EXTRA)
    canvas_module = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    colormap = require("matplotlib", _PLOT_EXTRA).colormaps["tab10"]
    adata = snapatac2.read(str(source), backed="r")
    try:
        try:
            groups = [item for item in adata.obs[label]]
        except (KeyError, RuntimeError) as exc:
            raise ConfigurationError(f"groupby column not found in obs: {label}") from exc
        coverage = internal.get_coverage(adata, text, groups)
    finally:
        close = getattr(adata, "close", None)
        if callable(close):
            close()
    tracks = sorted(coverage.items())
    if not tracks:
        raise ConfigurationError("The region has no coverage tracks.")
    signals = [np.asarray(signal, dtype=float) for _, signal in tracks]
    positions = int(signals[0].size)
    top = max(float(signal.max()) for signal in signals)
    figure = figure_module.Figure(figsize=(6.0, 1.2 * len(tracks)), dpi=100, layout="constrained")
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.subplots(len(tracks), 1, sharex=True, squeeze=False).ravel()
    for index, (axis, (name, signal)) in enumerate(
        zip(axes, zip([name for name, _ in tracks], signals, strict=False), strict=False)
    ):
        axis.fill_between(range(signal.size), 0, signal, color=colormap(index % 10))
        axis.set_title(str(name), fontsize=7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_ylim(0, top)
    axes[-1].set_xticks([0, positions - 1])
    axes[-1].set_xticklabels([str(start), str(end)])
    axes[-1].set_xlabel(text)
    figure.supylabel("RPM")
    with require("matplotlib", _PLOT_EXTRA).rc_context({"savefig.bbox": "standard"}):
        figure.savefig(str(target), format="png", dpi=100)
    data = {
        "region": text,
        "groupby": label,
        "positions": positions,
        "tracks": [
            {"group": str(name), "max": float(signal.max()), "signal": _bounded(np, signal)}
            for (name, _), signal in zip(tracks, signals, strict=False)
        ],
    }
    return {
        "status": "completed",
        "operation": "coverage",
        "image_path": str(target),
        "image_bytes": target.stat().st_size,
        "parameters": {"input_h5ad": str(source), "region": text, "groupby": label},
        "data": data,
    }


__all__ = [
    "fetch_snapatac2_cis_bp_motifs",
    "fetch_snapatac2_colon_transverse",
    "fetch_snapatac2_genome_files",
    "fetch_snapatac2_hea_ccre",
    "fetch_snapatac2_meuleman_2020_motifs",
    "fetch_snapatac2_pbmc5k",
    "fetch_snapatac2_pbmc10k_multiome",
    "fetch_snapatac2_pbmc500",
    "load_snapatac2_dataset",
    "plot_snapatac2_coverage_png",
    "plot_snapatac2_frag_size_distr_png",
    "plot_snapatac2_heatmap_png",
    "plot_snapatac2_kde2d_png",
    "plot_snapatac2_motif_enrichment_png",
    "plot_snapatac2_network_edge_stat_png",
    "plot_snapatac2_network_scores_png",
    "plot_snapatac2_regions_png",
    "plot_snapatac2_scatter3d_png",
    "plot_snapatac2_scatter_png",
    "plot_snapatac2_spectral_eigenvalues_png",
    "plot_snapatac2_tsse_png",
    "plot_snapatac2_umap_png",
    "plot_snapatac2_valid_cells_png",
    "render_snapatac2_figure_png",
]

_register_public_verbose(globals())
