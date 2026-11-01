"""cooler and cooltools operations that DNAKit did not yet expose (integration pack P38).

Covers region-to-bin look-ups, bin-table queries, sanitising of pixel and contact records, cooler file-type
detection, the building blocks of the cooltools dot caller (kernels, scoring stages, FDR thresholds, q-value
annotation, clustering, HiCCUPS filtering), pixel counting per block, combining and interpolating log-binned
expected, expected snippets around features and iterative correction of rectangular matrices.

Install with ``python -m pip install "dnakit[bio-cooltools]"`` (cooler, cooltools, pandas, numpy). Each function
imports its dependency only when it is called. Tables travel as tab-separated files with a header line; results
are JSON-friendly dictionaries. Output files are always new files (an existing path is refused). Tested with
cooler 0.10.4 and cooltools 0.7.1.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p38_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_COOLTOOLS_EXTRA = "bio-cooltools"
_MAX_ROWS = 5_000_000
_MAX_BINS = 200_000
_MAX_DENSE_CELLS = 25_000_000
_KERNEL_TYPES = ("donut", "vertical", "horizontal", "lowleft", "upright")
_DOT_STAGES = ("scored", "histogram", "thresholds", "pixels", "qvalues")


def _open_cooler(path: Path) -> Any:
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        clr = cooler.Cooler(str(path))
        int(clr.info["nbins"])
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {path} ({exc})") from exc
    return clr


def _read_table(path: Path, name: str, *, sep: str = "\t", header: Any = 0, **kw: Any) -> Any:
    pandas = require("pandas", _COOLER_EXTRA)
    try:
        frame = pandas.read_csv(path, sep=sep, header=header, nrows=_MAX_ROWS + 1, **kw)
    except (OSError, ValueError, pandas.errors.ParserError) as exc:
        raise ConfigurationError(f"{name} is not a readable table: {exc}") from exc
    if len(frame) > _MAX_ROWS:
        raise ConfigurationError(f"{name} has more than {_MAX_ROWS} rows.")
    return frame


def _write_table(frame: Any, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def _region(value: object, name: str = "region") -> str:
    return h.text(value, name)


def _chromsizes(value: Mapping[str, int]) -> Any:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError("chromsizes must be a non-empty mapping of chromosome to length.")
    pandas = require("pandas", _COOLER_EXTRA)
    sizes = {
        h.text(k, "chromosome name"): h.integer(v, "chromosome length", 1, 4_000_000_000)
        for k, v in value.items()
    }
    return pandas.Series(sizes, dtype="int64")


def _bins_from_chromsizes(chromsizes: Mapping[str, int], binsize: int) -> Any:
    cooler = require("cooler", _COOLER_EXTRA)
    size = h.integer(binsize, "binsize", 1, 100_000_000)
    series = _chromsizes(chromsizes)
    if int(sum(-(-int(v) // size) for v in series)) > _MAX_BINS:
        raise ConfigurationError(f"The bin layout has more than {_MAX_BINS} bins.")
    return cooler.binnify(series, size)


def cooler_region_bin_ids(
    cool_path: str, region: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Map a genomic region to cooler bin IDs, using ``cooler.Cooler.offset`` and ``Cooler.extent``.

    ``region`` is a UCSC-style string such as ``chr1`` or ``chr1:1,000,000-2,000,000``. Returns the ID of
    the first bin that overlaps the region (``offset``) and the half-open range of bin IDs ``[lo, hi)``
    (``extent``), plus the number of bins.
    """

    clr = _open_cooler(h.existing(cool_path, "cool_path"))
    text = _region(region)
    try:
        offset = int(clr.offset(text))
        lo, hi = (int(v) for v in clr.extent(text))
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"Cannot resolve region {text!r} in this cooler: {exc}") from exc
    return {
        "status": "completed",
        "operation": "cooler_region_bin_ids",
        "region": text,
        "offset": offset,
        "extent": [lo, hi],
        "bin_count": hi - lo,
        "provenance": h.provenance("cooler", "cooler.Cooler.offset and cooler.Cooler.extent"),
    }


def query_cooler_bins(
    cool_path: str,
    *,
    region: str | None = None,
    columns: Sequence[str] | None = None,
    limit: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query the bin table of a cooler (chrom, start, end and balancing weights), using ``Cooler.bins``.

    ``region`` restricts the table to a UCSC-style region (all bins when omitted); ``columns`` selects columns
    such as ``["chrom", "start", "end", "weight"]`` (all stored columns when omitted); ``limit`` caps the
    number of returned rows (the total is reported).
    """

    clr = _open_cooler(h.existing(cool_path, "cool_path"))
    cap = h.integer(limit, "limit", 1, 100_000)
    wanted = None if columns is None else [h.text(c, "column") for c in columns]
    selector = clr.bins()
    try:
        table = selector[:] if region is None else selector.fetch(_region(region))
        if wanted:
            table = table[wanted]
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Cannot query the bin table: {exc}") from exc
    return {
        "status": "completed",
        "operation": "query_cooler_bins",
        "region": region,
        "bin_count": int(len(table)),
        "returned": min(int(len(table)), cap),
        "columns": [str(c) for c in table.columns],
        "bins": h.records(table.reset_index(), cap),
        "provenance": h.provenance("cooler", "cooler.Cooler.bins().fetch"),
    }


def sanitize_cooler_pixels(
    pixels_path: str,
    chromsizes: Mapping[str, int],
    binsize: int,
    *,
    output_path: str | None = None,
    is_one_based: bool = False,
    tril_action: str = "reflect",
    sort: bool = True,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sanitise pre-binned pixel records (bin-ID base, upper-triangle orientation), using
    ``cooler.create.sanitize_pixels``.

    The input is a tab-separated file with a header containing ``bin1_id``, ``bin2_id`` and optional value
    columns such as ``count``. Bin IDs refer to fixed bins of ``binsize`` over ``chromsizes``.
    ``is_one_based`` shifts the IDs to zero-based, ``tril_action`` is ``reflect`` (flip lower-triangle
    records), ``drop`` or ``raise``. The sanitised table is written to ``output_path`` when given.
    """

    source = h.existing(pixels_path, "pixels_path")
    target = None if output_path is None else h.new_path(output_path, "output_path")
    action = h.choice(tril_action, "tril_action", ("reflect", "drop", "raise"))
    cooler = require("cooler", _COOLER_EXTRA)
    bins = _bins_from_chromsizes(chromsizes, binsize)
    frame = _read_table(source, "pixels_path")
    if not {"bin1_id", "bin2_id"}.issubset(frame.columns):
        raise ConfigurationError("pixels_path needs the columns bin1_id and bin2_id.")
    cap = h.integer(limit, "limit", 0, 100_000)
    try:
        pipeline = cooler.create.sanitize_pixels(
            bins,
            is_one_based=h.flag(is_one_based, "is_one_based"),
            tril_action=action,
            sort=h.flag(sort, "sort"),
        )
        clean = pipeline(frame)
    except Exception as exc:  # cooler raises BadInputError and ValueError for invalid records
        raise ConfigurationError(f"cooler could not sanitise the pixels: {exc}") from exc
    if target is not None:
        _write_table(clean, target)
    return {
        "status": "completed",
        "operation": "sanitize_cooler_pixels",
        "input_rows": int(len(frame)),
        "output_rows": int(len(clean)),
        "tril_action": action,
        "output_path": None if target is None else str(target),
        "pixels": h.records(clean, cap),
        "provenance": h.provenance("cooler", "cooler.create.sanitize_pixels"),
    }


def sanitize_cooler_records(
    records_path: str,
    chromsizes: Mapping[str, int],
    binsize: int,
    *,
    schema: str = "bg2",
    output_path: str | None = None,
    is_one_based: bool = False,
    tril_action: str = "reflect",
    sort: bool = True,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sanitise genomic contact records and assign bin IDs, using ``cooler.create.sanitize_records``.

    ``schema`` selects the column layout: ``bg2`` (chrom1, start1, end1, chrom2, start2, end2, count) or
    ``pairs`` (chrom1, pos1, chrom2, pos2 ...). The tab-separated input must have a header with those names.
    Records are binned with fixed ``binsize`` bins over ``chromsizes``; ``tril_action`` is ``reflect``,
    ``drop`` or ``raise``. The result has ``bin1_id`` and ``bin2_id`` added and is written to ``output_path``.
    """

    source = h.existing(records_path, "records_path")
    target = None if output_path is None else h.new_path(output_path, "output_path")
    layout = h.choice(schema, "schema", ("bg2", "pairs"))
    action = h.choice(tril_action, "tril_action", ("reflect", "drop", "raise"))
    cooler = require("cooler", _COOLER_EXTRA)
    bins = _bins_from_chromsizes(chromsizes, binsize)
    frame = _read_table(source, "records_path")
    cap = h.integer(limit, "limit", 0, 100_000)
    try:
        pipeline = cooler.create.sanitize_records(
            bins,
            schema=layout,
            is_one_based=h.flag(is_one_based, "is_one_based"),
            tril_action=action,
            sort=h.flag(sort, "sort"),
        )
        clean = pipeline(frame)
    except Exception as exc:  # cooler raises BadInputError, KeyError and ValueError
        raise ConfigurationError(f"cooler could not sanitise the records: {exc}") from exc
    if target is not None:
        _write_table(clean, target)
    return {
        "status": "completed",
        "operation": "sanitize_cooler_records",
        "schema": layout,
        "input_rows": int(len(frame)),
        "output_rows": int(len(clean)),
        "output_path": None if target is None else str(target),
        "records": h.records(clean, cap),
        "provenance": h.provenance("cooler", "cooler.create.sanitize_records"),
    }


def detect_cooler_file_type(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Detect whether a file or URI is a single-resolution cooler, a multi-resolution cooler or a
    single-cell cooler, using ``cooler.fileops.is_cooler``, ``is_multires_file`` and ``is_scool_file``.

    ``path`` is a file path, optionally with a ``::/group`` suffix. ``file_type`` is ``multires``, ``scool``,
    ``cooler`` or ``unknown``.
    """

    text = h.text(path, "path")
    fileops = require("cooler.fileops", _COOLER_EXTRA)
    base = text.split("::", 1)[0]
    if not Path(base).expanduser().is_file():
        raise ConfigurationError(f"Input file does not exist: {base}")
    try:
        is_cool = bool(fileops.is_cooler(text))
        is_multires = bool(fileops.is_multires_file(base))
        is_scool = bool(fileops.is_scool_file(base))
    except OSError as exc:
        raise ConfigurationError(f"Cannot open the file: {exc}") from exc
    if is_multires:
        kind = "multires"
    elif is_scool:
        kind = "scool"
    elif is_cool:
        kind = "cooler"
    else:
        kind = "unknown"
    return {
        "status": "completed",
        "operation": "detect_cooler_file_type",
        "path": text,
        "is_cooler": is_cool,
        "is_multires": is_multires,
        "is_scool": is_scool,
        "file_type": kind,
        "provenance": h.provenance("cooler", "cooler.fileops.is_cooler / is_multires_file / is_scool_file"),
    }


def build_dot_calling_kernels(
    *,
    binsize: int | None = None,
    kernel_half_width: int | None = None,
    donut_inner_half_width: int | None = None,
    kernel_types: Sequence[str] = ("donut", "vertical", "horizontal", "lowleft"),
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the convolution kernels of the dot caller (donut, vertical, horizontal, lowleft, upright).

    With ``binsize`` the cooltools recommendation (``dotfinder.recommend_kernels``) is returned for 4 to 28 kb
    resolutions (donut w, p = 7, 4 near 5 kb; 5, 2 near 10 kb; 3, 1 near 25 kb) using the four HiCCUPS kernel
    types. With ``kernel_half_width`` (w) and ``donut_inner_half_width`` (p) the kernels of the requested
    ``kernel_types`` are built directly (``cooltools.lib.numutils.get_kernel``).
    """

    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    if binsize is not None:
        size = h.integer(binsize, "binsize", 1, 100_000_000)
        try:
            kernels = dotfinder.recommend_kernels(size)
        except ValueError as exc:
            raise ConfigurationError(str(exc)) from exc
        source = "recommend_kernels"
    else:
        if kernel_half_width is None or donut_inner_half_width is None:
            raise ConfigurationError("Give binsize, or both kernel_half_width and donut_inner_half_width.")
        w = h.integer(kernel_half_width, "kernel_half_width", 1, 50)
        p = h.integer(donut_inner_half_width, "donut_inner_half_width", 0, w)
        kinds = [h.choice(k, "kernel type", _KERNEL_TYPES) for k in kernel_types]
        if not kinds:
            raise ConfigurationError("kernel_types must not be empty.")
        kernels = {k: numutils.get_kernel(w, p, k) for k in kinds}
        source = "get_kernel"
    first = next(iter(kernels.values()))
    return {
        "status": "completed",
        "operation": "build_dot_calling_kernels",
        "source": source,
        "width": int(first.shape[0]),
        "kernel_half_width": int((first.shape[0] - 1) // 2),
        "kernels": {name: [[int(v) for v in row] for row in k] for name, k in kernels.items()},
        "provenance": h.provenance("cooltools", "cooltools.api.dotfinder.recommend_kernels / lib.numutils.get_kernel"),
    }


def validate_dot_calling_kernels(
    kernels: Mapping[str, Sequence[Sequence[int]]],
    binsize: int,
    *,
    max_nans_tolerated: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate user-supplied dot-calling kernels, using ``cooltools.api.dotfinder.is_compatible_kernels``.

    ``kernels`` maps names to square integer matrices of the same odd size greater than 1; ``max_nans_tolerated``
    must not exceed the kernel width. Returns ``valid`` and, when invalid, the upstream ``reason``.
    """

    numpy = require("numpy", _COOLER_EXTRA)
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    if not isinstance(kernels, Mapping) or not kernels:
        raise ConfigurationError("kernels must be a non-empty mapping of name to matrix.")
    size = h.integer(binsize, "binsize", 1, 100_000_000)
    nans = h.integer(max_nans_tolerated, "max_nans_tolerated", 0, 1000)
    try:
        arrays = {str(k): numpy.asarray(v, dtype=float) for k, v in kernels.items()}
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"kernels must be numeric matrices: {exc}") from exc
    try:
        dotfinder.is_compatible_kernels(arrays, size, nans)
        valid, reason = True, None
    except ValueError as exc:
        valid, reason = False, str(exc)
    return {
        "status": "completed",
        "operation": "validate_dot_calling_kernels",
        "valid": valid,
        "reason": reason,
        "kernel_names": list(arrays),
        "provenance": h.provenance("cooltools", "cooltools.api.dotfinder.is_compatible_kernels"),
    }


def _dot_setup(
    cool_path: str,
    expected_path: str,
    weight: str,
    expected_value_col: str,
    kernel_half_width: int | None,
    donut_inner_half_width: int | None,
    max_loci_separation: int,
    tile_size: int,
    n_lambda_bins: int,
) -> dict[str, Any]:
    np = require("numpy", _COOLER_EXTRA)
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    clr = _open_cooler(h.existing(cool_path, "cool_path"))
    name = h.text(weight, "weight")
    if name not in clr.bins()[:1].columns:
        raise ConfigurationError(f"This cooler has no {name!r} column; balance it first.")
    value_col = h.text(expected_value_col, "expected_value_col")
    expected = _read_table(h.existing(expected_path, "expected_path"), "expected_path")
    if not {"region1", "region2", "dist", value_col}.issubset(expected.columns):
        raise ConfigurationError(f"expected_path needs region1, region2, dist and {value_col}.")
    expected = expected.set_index(["region1", "region2", "dist"]).sort_index()
    binsize = int(clr.binsize)
    if kernel_half_width is None:
        try:
            kernels = dotfinder.recommend_kernels(binsize)
        except ValueError as exc:
            raise ConfigurationError(str(exc)) from exc
    else:
        w = h.integer(kernel_half_width, "kernel_half_width", 2, 20)
        p = h.integer(donut_inner_half_width if donut_inner_half_width is not None else 1,
                      "donut_inner_half_width", 1, w - 1)
        kernels = {k: numutils.get_kernel(w, p, k) for k in ("donut", "vertical", "horizontal", "lowleft")}
    half = int((max(len(k) for k in kernels.values()) - 1) // 2)
    separation = dotfinder.bp_to_bins(h.integer(max_loci_separation, "max_loci_separation", 1, 1_000_000_000), binsize)
    tile_bins = dotfinder.bp_to_bins(h.integer(tile_size, "tile_size", 1, 1_000_000_000), binsize)
    bins = h.integer(n_lambda_bins, "n_lambda_bins", 40, 50)
    base = 2 ** (1 / 3)
    ledges = np.concatenate(([-np.inf], np.logspace(0, bins - 1, num=bins, base=base, dtype=np.float64), [np.inf]))
    view = common.make_cooler_view(clr)
    tiles = list(dotfinder.generate_tiles_diag_band(clr, view, half, tile_bins, separation))
    return {
        "clr": clr, "expected": expected, "weight": name, "value_col": value_col, "kernels": kernels,
        "separation": separation, "ledges": ledges, "tiles": tiles, "dotfinder": dotfinder,
    }


def calculate_dot_calling_stage(
    cool_path: str,
    expected_path: str,
    output_path: str,
    *,
    stage: str = "pixels",
    weight: str = "weight",
    expected_value_col: str = "balanced.avg",
    kernel_half_width: int | None = None,
    donut_inner_half_width: int | None = None,
    max_loci_separation: int = 10_000_000,
    tile_size: int = 5_000_000,
    n_lambda_bins: int = 40,
    lambda_bin_fdr: float = 0.1,
    max_nans_tolerated: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the lambda-binning dot-calling procedure of cooltools up to one stage and write that stage's table.

    ``stage`` selects the output: ``scored`` (locally adjusted expected per kernel and Poisson-ready counts for
    every pixel of every tile: ``score_tile`` with ``get_adjusted_expected_tile_some_nans``), ``histogram``
    (lambda-bin by observed-count histogram per kernel: ``scoring_and_histogramming_step``), ``thresholds`` (FDR
    threshold per lambda bin and kernel: ``determine_thresholds``), ``pixels`` (pixels passing the thresholds for
    every kernel: ``scoring_and_extraction_step`` / ``extract_scored_pixels``) or ``qvalues`` (those pixels
    annotated with q-values: ``annotate_pixels_with_qvalues``, then genomic coordinates chrom1 ... end2 as cooltools ``dots`` does). ``cool_path`` must be balanced; ``expected_path``
    is a tab-separated cis expected table (region1, region2, dist, balanced.avg). Kernels are the cooltools
    recommendation for the resolution unless ``kernel_half_width`` (and ``donut_inner_half_width``) are given.
    The table goes to ``output_path`` (histogram and thresholds in long format with a ``kernel`` column).
    """

    which = h.choice(stage, "stage", _DOT_STAGES)
    target = h.new_path(output_path, "output_path")
    fdr = h.number(lambda_bin_fdr, "lambda_bin_fdr", 1e-6, 1.0)
    nans = h.integer(max_nans_tolerated, "max_nans_tolerated", 0, 1000)
    pandas = require("pandas", _COOLER_EXTRA)
    s = _dot_setup(cool_path, expected_path, weight, expected_value_col, kernel_half_width,
                   donut_inner_half_width, max_loci_separation, tile_size, n_lambda_bins)
    df = s["dotfinder"]
    common = dict(clr=s["clr"], expected_indexed=s["expected"], expected_value_col=s["value_col"],
                  clr_weight_name=s["weight"], kernels=s["kernels"], max_nans_tolerated=nans)
    try:
        if which == "scored":
            parts = [df.score_tile(tile, band_to_cover=s["separation"], **common) for tile in s["tiles"]]
            table = pandas.concat(parts, ignore_index=True).sort_values(["bin1_id", "bin2_id"]) if parts else pandas.DataFrame()
            table = table.drop_duplicates().reset_index(drop=True)
        else:
            hist = df.scoring_and_histogramming_step(
                tiles=s["tiles"], ledges=s["ledges"], loci_separation_bins=s["separation"], nproc=1, **common
            )
            if which == "histogram":
                frames = []
                for kernel, matrix in hist.items():
                    long = matrix.copy()
                    long.columns = [str(c) for c in long.columns]
                    long = long.reset_index().melt(id_vars=long.index.name or "count",
                                                    var_name="lambda_bin", value_name="n_pixels")
                    long.insert(0, "kernel", kernel)
                    frames.append(long)
                table = pandas.concat(frames, ignore_index=True)
            else:
                thresholds, qvalues = df.determine_thresholds(hist, fdr)
                if which == "thresholds":
                    table = pandas.concat(
                        [pandas.DataFrame({"kernel": k, "lambda_bin": [str(i) for i in v.index],
                                           "threshold_count": v.to_numpy()}) for k, v in thresholds.items()],
                        ignore_index=True)
                else:
                    pixels = df.scoring_and_extraction_step(
                        tiles=s["tiles"], ledges=s["ledges"], thresholds=thresholds,
                        loci_separation_bins=s["separation"], nproc=1, **common)
                    table = pixels
                    if which == "qvalues":
                        table = df.annotate_pixels_with_qvalues(pixels, qvalues)
                        cooler = require("cooler", _COOLER_EXTRA)
                        table = cooler.annotate(table, s["clr"].bins()[["chrom", "start", "end"]], replace=True)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools dot calling failed at stage {which}: {exc}") from exc
    _write_table(table, target)
    return {
        "status": "completed",
        "operation": "calculate_dot_calling_stage",
        "stage": which,
        "kernels": list(s["kernels"]),
        "tiles": len(s["tiles"]),
        "lambda_bin_fdr": fdr,
        "rows": int(len(table)),
        "columns": [str(c) for c in table.columns],
        "output_path": str(target),
        "provenance": h.provenance("cooltools", "cooltools.api.dotfinder lambda-binning (HiCCUPS-like) stages"),
    }


def cluster_dot_pixels(
    pixels_path: str,
    output_path: str,
    *,
    clustering_radius: int = 20_000,
    count_column: str = "count",
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster enriched pixels into dot centroids with Birch clustering, using
    ``cooltools.api.dotfinder.clustering_step`` (which calls ``clust_2D_pixels``).

    ``pixels_path`` is a tab-separated table of enriched pixels annotated with chrom1, start1, chrom2, start2 and
    the observed count column (the ``qvalues`` output of ``calculate_dot_calling_stage`` joined to genomic
    coordinates). Clustering is done per chromosome (or per ``region`` column when present) on the start
    coordinates with ``clustering_radius`` in bp; the brightest pixel of every cluster is kept with its
    centroid (``cstart1``, ``cstart2``), ``c_label`` and ``c_size``.
    """

    source = h.existing(pixels_path, "pixels_path")
    target = h.new_path(output_path, "output_path")
    radius = h.integer(clustering_radius, "clustering_radius", 1, 100_000_000)
    count = h.text(count_column, "count_column")
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    frame = _read_table(source, "pixels_path")
    try:
        kwargs = {"assigned_regions_name": "region"} if "region" in frame.columns else {}
        centroids = dotfinder.clustering_step(frame, radius, obs_raw_name=count, **kwargs)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools could not cluster the pixels: {exc}") from exc
    centroids = centroids.reset_index(drop=True)
    _write_table(centroids, target)
    return {
        "status": "completed",
        "operation": "cluster_dot_pixels",
        "input_pixels": int(len(frame)),
        "clusters": int(len(centroids)),
        "clustering_radius": radius,
        "output_path": str(target),
        "provenance": h.provenance("cooltools", "cooltools.api.dotfinder.clustering_step and clust_2D_pixels"),
    }


def filter_dot_centroids(
    centroids_path: str,
    output_path: str,
    *,
    count_column: str = "count",
    enrichment_factor_vh: float = 1.5,
    enrichment_factor_d_and_ll: float = 1.75,
    enrichment_factor_d_or_ll: float = 2.0,
    fdr_orphan_threshold: float = 0.02,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter dot centroids by the HiCCUPS enrichment criteria, using
    ``cooltools.api.dotfinder.cluster_filtering_hiccups``.

    The tab-separated input is the output of ``cluster_dot_pixels`` scored with the four HiCCUPS kernels (donut,
    vertical, horizontal, lowleft; columns ``la_exp.<kernel>.value`` and ``la_exp.<kernel>.qval``). Defaults are
    those of cooltools (1.5, 1.75, 2.0 and 0.02).
    """

    source = h.existing(centroids_path, "centroids_path")
    target = h.new_path(output_path, "output_path")
    dotfinder = require("cooltools.api.dotfinder", _COOLTOOLS_EXTRA)
    frame = _read_table(source, "centroids_path")
    try:
        kept = dotfinder.cluster_filtering_hiccups(
            frame,
            obs_raw_name=h.text(count_column, "count_column"),
            enrichment_factor_vh=h.number(enrichment_factor_vh, "enrichment_factor_vh", 0, 1000),
            enrichment_factor_d_and_ll=h.number(enrichment_factor_d_and_ll, "enrichment_factor_d_and_ll", 0, 1000),
            enrichment_factor_d_or_ll=h.number(enrichment_factor_d_or_ll, "enrichment_factor_d_or_ll", 0, 1000),
            FDR_orphan_threshold=h.number(fdr_orphan_threshold, "fdr_orphan_threshold", 0, 4),
        )
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    _write_table(kept, target)
    return {
        "status": "completed",
        "operation": "filter_dot_centroids",
        "input_centroids": int(len(frame)),
        "kept_centroids": int(len(kept)),
        "output_path": str(target),
        "provenance": h.provenance("cooltools", "cooltools.api.dotfinder.cluster_filtering_hiccups"),
    }


def count_block_pixels(
    width: int, height: int, *, bad_bins_x: int = 0, bad_bins_y: int = 0, verbose: bool = True
) -> dict[str, Any]:
    """Count total and masked (bad) pixels in a rectangular block of a contact map, using
    ``cooltools.api.expected.count_all_pixels_per_block`` and ``count_bad_pixels_per_block``.

    ``width`` and ``height`` are the block sides in bins; ``bad_bins_x`` and ``bad_bins_y`` count masked bins on
    each side. Returns total, bad and valid pixel counts.
    """

    expected = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    x = h.integer(width, "width", 1, 10**9)
    y = h.integer(height, "height", 1, 10**9)
    bx = h.integer(bad_bins_x, "bad_bins_x", 0, x)
    by = h.integer(bad_bins_y, "bad_bins_y", 0, y)
    total = int(expected.count_all_pixels_per_block(x, y))
    bad = int(expected.count_bad_pixels_per_block(x, y, bx, by))
    return {
        "status": "completed",
        "operation": "count_block_pixels",
        "total_pixels": total,
        "bad_pixels": bad,
        "valid_pixels": total - bad,
        "provenance": h.provenance("cooltools", "cooltools.api.expected.count_all/bad_pixels_per_block"),
    }


def combine_logbinned_expected(
    binned_expected_path: str,
    output_path: str,
    *,
    binned_slope_path: str | None = None,
    pc_name: str = "balanced.avg",
    spread_funcs: str = "logstd",
    spread_funcs_slope: str = "std",
    minmax_drop_bins: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine by-region log-binned expected into a genome-wide P(s) curve and slope, using
    ``cooltools.api.expected.combine_binned_expected``.

    ``binned_expected_path`` is the tab-separated output of ``logbin_expected`` (columns region1, region2,
    dist_bin_start, dist_bin_end, n_valid, balanced.avg ...). With ``binned_slope_path`` (the slope table of
    ``logbin_expected`` as a tab-separated file without index) the spread of the slopes is estimated too. ``spread_funcs`` is ``logstd``, ``std`` or
    ``minmax``. The combined P(s) table is written to ``output_path``.
    """

    source = h.existing(binned_expected_path, "binned_expected_path")
    target = h.new_path(output_path, "output_path")
    expected = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    frame = _read_table(source, "binned_expected_path")
    slope = None
    if binned_slope_path is not None:
        slope = _read_table(h.existing(binned_slope_path, "binned_slope_path"), "binned_slope_path")
    try:
        scal, slope_df = expected.combine_binned_expected(
            frame,
            binned_exp_slope=slope,
            Pc_name=h.text(pc_name, "pc_name"),
            spread_funcs=h.choice(spread_funcs, "spread_funcs", ("logstd", "std", "minmax")),
            spread_funcs_slope=h.choice(spread_funcs_slope, "spread_funcs_slope", ("std", "minmax")),
            minmax_drop_bins=h.integer(minmax_drop_bins, "minmax_drop_bins", 0, 100),
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not combine the expected: {exc}") from exc
    scal = scal.reset_index(drop=True)
    _write_table(scal, target)
    return {
        "status": "completed",
        "operation": "combine_logbinned_expected",
        "rows": int(len(scal)),
        "columns": [str(c) for c in scal.columns],
        "slope_rows": None if slope_df is None else int(len(slope_df)),
        "output_path": str(target),
        "preview": h.records(scal, 5),
        "provenance": h.provenance("cooltools", "cooltools.api.expected.combine_binned_expected"),
    }


def interpolate_logbinned_expected(
    expected_path: str,
    binned_expected_path: str,
    output_path: str,
    *,
    columns: Sequence[str] = ("balanced.avg",),
    kind: str = "quadratic",
    by_region: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Interpolate log-binned expected back to per-diagonal values, using
    ``cooltools.api.expected.interpolate_expected``.

    ``expected_path`` is the per-diagonal expected (region1, region2, dist ...), ``binned_expected_path`` the
    log-binned or combined expected (the tab-separated ``logbin_expected`` or ``combine_logbinned_expected``
    table with ``dist.avg``). ``kind`` is a scipy ``interp1d`` kind. With
    ``by_region`` false one combined curve is used for every region.
    """

    exp_source = h.existing(expected_path, "expected_path")
    bin_source = h.existing(binned_expected_path, "binned_expected_path")
    target = h.new_path(output_path, "output_path")
    cols = [h.text(c, "column") for c in columns]
    expected = require("cooltools.api.expected", _COOLTOOLS_EXTRA)
    exp = _read_table(exp_source, "expected_path")
    binned = _read_table(bin_source, "binned_expected_path")
    try:
        out = expected.interpolate_expected(
            exp, binned, columns=cols,
            kind=h.choice(kind, "kind", ("linear", "nearest", "slinear", "quadratic", "cubic")),
            by_region=h.flag(by_region, "by_region"),
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not interpolate the expected: {exc}") from exc
    _write_table(out, target)
    return {
        "status": "completed",
        "operation": "interpolate_logbinned_expected",
        "rows": int(len(out)),
        "columns": [str(c) for c in out.columns],
        "output_path": str(target),
        "preview": h.records(out, 5),
        "provenance": h.provenance("cooltools", "cooltools.api.expected.interpolate_expected"),
    }


def snip_expected_around_features(
    cool_path: str,
    expected_path: str,
    features_path: str,
    *,
    flank: int = 100_000,
    expected_value_col: str = "balanced.avg",
    min_diag: int = 2,
    limit: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut expected-matrix snippets around features, using ``cooltools.api.snipping.ExpectedSnipper``
    (``select`` and ``snip``).

    ``features_path`` is a BED-like tab-separated file with a header (chrom, start, end) of features to centre
    on; every snippet spans the feature midpoint plus and minus ``flank`` bp, aligned to bins. Diagonals below
    ``min_diag`` are NaN (reported as null); windows leaving the chromosome are all null. Features on different
    chromosomes use cis expected of that chromosome (the whole chromosome is the region).
    """

    clr = _open_cooler(h.existing(cool_path, "cool_path"))
    exp_table = _read_table(h.existing(expected_path, "expected_path"), "expected_path")
    feats = _read_table(h.existing(features_path, "features_path"), "features_path")
    if not {"chrom", "start", "end"}.issubset(feats.columns):
        raise ConfigurationError("features_path needs the columns chrom, start and end.")
    pad = h.integer(flank, "flank", 1, 100_000_000)
    diag = h.integer(min_diag, "min_diag", 0, 1_000_000)
    cap = h.integer(limit, "limit", 1, 1000)
    snipping = require("cooltools.api.snipping", _COOLTOOLS_EXTRA)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    binsize = int(clr.binsize)
    view = common.make_cooler_view(clr)
    try:
        snipper = snipping.ExpectedSnipper(
            clr, exp_table, view_df=view, min_diag=diag, expected_value_col=h.text(expected_value_col, "expected_value_col")
        )
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    selected: dict[str, Any] = {}
    snippets = []
    for row in feats.head(cap).itertuples(index=False):
        chrom = str(row.chrom)
        if chrom not in set(view["name"]):
            raise ConfigurationError(f"Chromosome {chrom!r} is not in the cooler.")
        if chrom not in selected:
            selected[chrom] = snipper.select(chrom, chrom)
        mid = (int(row.start) + int(row.end)) // 2
        lo = (mid - pad) // binsize * binsize
        hi = (mid + pad) // binsize * binsize + binsize
        try:
            matrix = snipper.snip(selected[chrom], chrom, chrom, (lo, hi, lo, hi))
        except (AssertionError, ValueError, IndexError):
            matrix = None
        snippets.append({
            "chrom": chrom, "start": int(row.start), "end": int(row.end), "window": [lo, hi],
            "matrix": None if matrix is None else [[h.plain(float(v)) for v in r] for r in matrix],
        })
    return {
        "status": "completed",
        "operation": "snip_expected_around_features",
        "binsize": binsize,
        "flank": pad,
        "snippets": snippets,
        "provenance": h.provenance("cooltools", "cooltools.api.snipping.ExpectedSnipper"),
    }


def iterative_correction_rectangular(
    matrix: Sequence[Sequence[float]],
    *,
    max_iter: int = 1000,
    tol: float = 1e-5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Iterative correction (matrix balancing) of an asymmetric, possibly rectangular contact matrix, using
    ``cooltools.lib.numutils.iterative_correction_asymmetric``.

    The matrix is rows by columns (for example two chromosomes in trans). Returns the corrected matrix, the bias
    vector of the rows (``row_bias``, length = rows) and of the columns (``column_bias``) and the convergence
    report. Note: upstream names its two returned vectors in the order column, row.
    """

    numpy = require("numpy", _COOLER_EXTRA)
    numutils = require("cooltools.lib.numutils", _COOLTOOLS_EXTRA)
    try:
        array = numpy.asarray(matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"matrix must be numeric: {exc}") from exc
    if array.ndim != 2 or min(array.shape) < 2 or array.size > _MAX_DENSE_CELLS:
        raise ConfigurationError("matrix must be 2-D with at least 2 rows and 2 columns.")
    if not numpy.isfinite(array).all() or (array < 0).any():
        raise ConfigurationError("matrix must contain finite non-negative values.")
    corrected, col_bias, row_bias, report = numutils.iterative_correction_asymmetric(
        array, max_iter=h.integer(max_iter, "max_iter", 1, 100_000), tol=h.number(tol, "tol", 1e-12, 1.0), verbose=False
    )
    return {
        "status": "completed",
        "operation": "iterative_correction_rectangular",
        "shape": [int(array.shape[0]), int(array.shape[1])],
        "converged": bool(report[0]),
        "iterations": int(report[1]),
        "corrected": [[h.plain(float(v)) for v in r] for r in corrected],
        "row_bias": [h.plain(float(v)) for v in row_bias],
        "column_bias": [h.plain(float(v)) for v in col_bias],
        "provenance": h.provenance("cooltools", "cooltools.lib.numutils.iterative_correction_asymmetric"),
    }


__all__ = [
    "build_dot_calling_kernels",
    "calculate_dot_calling_stage",
    "cluster_dot_pixels",
    "combine_logbinned_expected",
    "cooler_region_bin_ids",
    "count_block_pixels",
    "detect_cooler_file_type",
    "filter_dot_centroids",
    "interpolate_logbinned_expected",
    "iterative_correction_rectangular",
    "query_cooler_bins",
    "sanitize_cooler_pixels",
    "sanitize_cooler_records",
    "snip_expected_around_features",
    "validate_dot_calling_kernels",
]

_register_public_verbose(globals(), tuple(__all__))
