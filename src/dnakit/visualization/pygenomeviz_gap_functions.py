"""pyGenomeViz example data, a genome feature figure and three alignment figures.

- ``load_example_dataset()`` and ``load_example_gff()`` download pyGenomeViz example genomes from
  the pygenomeviz-data repository on GitHub into a cache folder (network, via
  ``pygenomeviz.utils``).
- ``plot_genome_features()`` draws the CDS features of GenBank genomes as stacked tracks in one PNG.
- ``compare_genomes_mmseqs()``, ``compare_genomes_mummer()`` and
  ``compare_genomes_progressive_mauve()`` run MMseqs2 (``easy-rbh``), MUMmer (``nucmer`` or
  ``promer``, ``delta-filter``, ``show-coords``) or progressiveMauve with fixed argument lists and
  no shell, then draw the genomes and their alignment links with pyGenomeViz.

The aligner programs are installed by the user. They are looked up on PATH unless an explicit path
is given, and a missing program raises BackendUnavailableError with an install hint. Output folders
and PNG files must not exist yet (a folder may also be empty). Results are JSON-friendly
dictionaries with status, backend, commands, outputs, a bounded log tail, the image path and the
underlying data.

Upstream code that is not called here: ``pygenomeviz.scripts.simpleplot.run`` and
``Genbank.extract_features`` read ``SeqFeature.strand``, which Biopython 1.88 no longer defines, so
CDS selection reads ``location.strand`` instead. The upstream aligner scripts run their programs
through ``shell=True`` strings; the argument lists are built here.
"""

from __future__ import annotations

import csv
import datetime as dt
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


def _require_isolated(module: str, extra: str) -> Any:
    """Import a pygenomeviz module without leaking its global matplotlib and warnings changes.

    Importing pygenomeviz sets ``savefig.bbox`` and other matplotlib settings and turns warnings off for
    the whole process, which would change every other DNAKit plot. The settings are restored afterwards.
    """
    import warnings

    try:
        import matplotlib
    except ImportError:  # require() below raises the installation hint
        return require(module, extra)
    saved_rc = matplotlib.rcParams.copy()
    saved_filters = list(warnings.filters)
    try:
        return require(module, extra)
    finally:
        matplotlib.rcParams.update(saved_rc)
        warnings.filters[:] = saved_filters


_EXTRA = "bio-pygenomeviz"
_SERVICE = "pygenomeviz-data repository on GitHub"
_MAX_TAIL = 4000
_LARGE_OUTPUT = 100_000_000
_NAME = re.compile(r"[A-Za-z0-9._-]+")
_PLOTSTYLES = ("bigarrow", "arrow")
_ALIGN_TYPES = ("left", "center", "right")
_TICK_STYLES = ("bar", "axis")
_SEQTYPES = ("protein", "nucleotide")
_NORMAL_LINK = "grey"
_INVERTED_LINK = "red"
_MMSEQS_HINT = (
    "Install MMseqs2 (conda install -c bioconda mmseqs2, or a release binary from "
    "https://github.com/soedinglab/MMseqs2) and put mmseqs on PATH or pass mmseqs_path."
)
_MUMMER_HINT = (
    "Install MUMmer 4 (conda install -c bioconda mummer4, or apt install mummer) and put nucmer, "
    "promer, delta-filter and show-coords on PATH or pass the folder that contains them as "
    "mummer_dir."
)
_MAUVE_HINT = (
    "Install progressiveMauve from the Mauve release (http://darlinglab.org/mauve/) and put "
    "progressiveMauve on PATH or pass progressive_mauve_path."
)


@dataclass(frozen=True)
class _Cds:
    """One CDS feature inside the genome range (0-based half-open coordinates)."""

    start: int
    end: int
    strand: int
    qualifiers: Mapping[str, Sequence[str]]


@dataclass
class _Steps:
    """Records of the external program calls made by one facade run."""

    backend: str
    commands: list[list[str]] = field(default_factory=list)
    return_codes: list[int] = field(default_factory=list)
    texts: list[str] = field(default_factory=list)

    def run(
        self,
        program: str,
        args: Sequence[str],
        *,
        cwd: Path,
        timeout_seconds: float,
        max_output_bytes: int = 1_000_000,
    ) -> tuple[int, str]:
        result = execute_bounded_command(
            program,
            list(args),
            backend_id=self.backend,
            cwd=cwd,
            timeout_seconds=timeout_seconds,
            max_output_bytes=max_output_bytes,
            verbose=False,
        )
        text = result.output or ""
        self.commands.append([Path(program).name, *args])
        self.return_codes.append(result.return_code)
        self.texts.append(text)
        return result.return_code, text

    def tail(self) -> str:
        return "\n".join(self.texts)[-_MAX_TAIL:]


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _cache(cache_dir: str | None) -> str | None:
    return None if cache_dir is None else str(Path(cache_dir).expanduser())


def _choice(value: object, allowed: Sequence[object], name: str) -> None:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}.")


def _positive_number(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive number.")
    return float(value)


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _percent(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 100:
        raise ConfigurationError(f"{name} must be a percentage between 0 and 100.")
    return float(value)


def _position(value: object) -> int:
    return int(str(value).replace("<", "").replace(">", ""))


def _genomes(resources: Sequence[str], *, minimum: int) -> list[Any]:
    """GenBank genome objects from specs 'file', 'file:start-end' or 'file::-1' (CLI syntax)."""
    if isinstance(resources, (str, bytes)) or len(resources) < minimum:
        raise ConfigurationError(f"gbk_resources must list at least {minimum} GenBank file(s).")
    specs = [str(item) for item in resources]
    for spec in specs:
        input_path(spec.split(":")[0])
    scripts = _require_isolated("pygenomeviz.scripts", _EXTRA)
    genomes = list(scripts.gbk_files2gbk_objects(specs))
    _check_names([genome.name for genome in genomes])
    return genomes


def _check_names(names: Sequence[str]) -> None:
    for name in names:
        if not _NAME.fullmatch(name):
            raise ConfigurationError(
                f"Genome name {name!r} must use only letters, digits, '.', '_' and '-' "
                "(it names output files)."
            )
    if len(set(names)) != len(names):
        raise ConfigurationError("Genome names must be unique.")


def _sequence_files(values: Sequence[str], *, minimum: int) -> list[Path]:
    if isinstance(values, (str, bytes)) or len(values) < minimum:
        raise ConfigurationError(f"seq_files must list at least {minimum} sequence files.")
    return [input_path(str(value)) for value in values]


def _new_png(value: str) -> Path:
    target = Path(value).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    if target.exists():
        raise ConfigurationError(f"image_path must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for image_path does not exist: {target.parent}")
    return target


def _check_new_folder(value: str) -> Path:
    folder = Path(value).expanduser().resolve()
    if folder.exists() and (not folder.is_dir() or any(folder.iterdir())):
        raise ConfigurationError(f"output_dir must be a new or empty folder: {folder}")
    return folder


def _executable(name: str, explicit: str | None, hint: str) -> str:
    """An external program: the explicit file path, or the first match on PATH."""
    if explicit is not None:
        candidate = Path(explicit).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name} was not found at {candidate}.", hint=hint)
        return str(candidate.resolve())
    found = shutil.which(name)
    if found is None:
        raise BackendUnavailableError(f"{name} was not found on PATH.", hint=hint)
    return found


def _cds_spans(gbk: Any) -> list[_Cds]:
    """CDS features with a translation inside the genome range, selected as pygenomeviz does.

    The strand comes from ``feature.location.strand``. Upstream reads ``feature.strand``, which
    Biopython 1.88 no longer defines (AttributeError: 'SeqFeature' has no attribute 'strand').
    """
    min_range, max_range = gbk.min_range, gbk.max_range
    spans: list[_Cds] = []
    base_len = 0
    for record in gbk.records:
        for feature in record.features:
            if feature.type != "CDS" or feature.qualifiers.get("translation", [None])[0] is None:
                continue
            strand = -1 if feature.location.strand == -1 else 1
            if strand == -1:
                start = _position(feature.location.parts[-1].start) + base_len
                end = _position(feature.location.parts[0].end) + base_len
            else:
                start = _position(feature.location.parts[0].start) + base_len
                end = _position(feature.location.parts[-1].end) + base_len
            if not min_range <= start <= end <= max_range:
                continue
            spans.append(_Cds(start, end, strand, feature.qualifiers))
        base_len += len(record.seq)
    return spans


def _genome_summary(gbk: Any) -> dict[str, Any]:
    return {
        "name": gbk.name,
        "min_range": int(gbk.min_range),
        "max_range": int(gbk.max_range),
        "range_size": int(gbk.range_size),
        "reverse": bool(gbk.reverse),
        "full_genome_length": int(gbk.full_genome_length),
        "cds_count": len(_cds_spans(gbk)),
    }


def _write_cds_protein_fasta(gbk: Any, path: Path) -> None:
    """Protein FASTA of the CDS features, with the ID layout pygenomeviz parses from RBH hits."""
    spans = _cds_spans(gbk)
    if not spans:
        raise ConfigurationError(f"Genome {gbk.name} has no CDS with a translation in its range.")
    seqio = require("Bio.SeqIO", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    records = []
    for index, cds in enumerate(spans, 1):
        protein_id = cds.qualifiers.get("protein_id", [None])[0]
        product = cds.qualifiers.get("product", [""])[0]
        location_id = f"|{cds.start}_{cds.end}_{cds.strand}|"
        if protein_id is None:
            seq_id = f"GENE{index:06d}{location_id}"
        else:
            seq_id = f"GENE{index:06d}_{protein_id}{location_id}"
        translation = cds.qualifiers["translation"][0]
        records.append(
            record_module.SeqRecord(seq_module.Seq(translation), id=seq_id, description=product)
        )
    seqio.write(records, str(path), "fasta-2line")


def _genome_viz(
    pgv: Any,
    *,
    fig_width: float,
    fig_track_height: float,
    align_type: str,
    tick_style: str | None,
    link_track_ratio: float,
) -> Any:
    return pgv.GenomeViz(
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        feature_track_ratio=1.0,
        link_track_ratio=link_track_ratio,
        tick_track_ratio=1.0,
        align_type=align_type,
        tick_style=tick_style,
        tick_labelsize=15,
        plot_size_thr=0,
    )


def _draw_cds(track: Any, gbk: Any, *, plotstyle: str, facecolor: str, linewidth: float) -> None:
    for cds in _cds_spans(gbk):
        track.add_feature(
            cds.start,
            cds.end,
            cds.strand,
            plotstyle=plotstyle,
            arrow_shaft_ratio=0.5,
            facecolor=cds.qualifiers.get("facecolor", [facecolor])[0],
            linewidth=linewidth,
        )


def _save(figure: Any, image: Path, dpi: int) -> None:
    plt = require("matplotlib.pyplot", _EXTRA)
    figure.savefig(image, dpi=dpi, pad_inches=0.5, bbox_inches="tight")
    figure.clear()
    plt.close(figure)


def _draw_feature_alignment(
    pgv: Any,
    genomes: Sequence[Any],
    alignments: Sequence[Any],
    image: Path,
    *,
    fig_width: float,
    dpi: int,
) -> None:
    """Genome CDS tracks, alignment links and an identity colour bar (pygenomeviz run layout)."""
    gv = _genome_viz(
        pgv,
        fig_width=fig_width,
        fig_track_height=1.0,
        align_type="center",
        tick_style=None,
        link_track_ratio=5.0,
    )
    for gbk in genomes:
        track = gv.add_feature_track(
            gbk.name, size=gbk.range_size, start_pos=gbk.min_range, labelsize=20
        )
        _draw_cds(track, gbk, plotstyle="bigarrow", facecolor="orange", linewidth=0.0)
    min_identity = int(min(ac.identity for ac in alignments)) if alignments else 0
    for ac in alignments:
        gv.add_link(
            ac.ref_link,
            ac.query_link,
            _NORMAL_LINK,
            _INVERTED_LINK,
            curve=True,
            v=ac.identity,
            vmin=min_identity,
        )
    figure = gv.plotfig(dpi=dpi)
    if alignments:
        bar_colors = [_NORMAL_LINK]
        if any(ac.is_inverted for ac in alignments):
            bar_colors.append(_INVERTED_LINK)
        gv.set_colorbar(figure, bar_colors, vmin=min_identity, bar_height=0.2, bar_width=0.01)
    _save(figure, image, dpi)


def _spaced_colors(cmap: str, count: int) -> list[str]:
    """``count`` distinct colours taken evenly from [0, 1) of a colormap (no wrap-around)."""
    mpl = require("matplotlib", _EXTRA)
    colors = require("matplotlib.colors", _EXTRA)
    colormap = mpl.colormaps[cmap]
    return [colors.to_hex(colormap(index / count), keep_alpha=True) for index in range(count)]


def _draw_block_alignment(
    pgv: Any,
    names: Sequence[str],
    bbone: Path,
    alignments: Sequence[Any],
    image: Path,
    *,
    fig_width: float,
    dpi: int,
    cmap: str,
) -> None:
    """Conserved-block tracks with alignment links (the pygenomeviz progressiveMauve layout)."""
    pmauve = _require_isolated("pygenomeviz.scripts.pmauve", _EXTRA)
    gv = _genome_viz(
        pgv,
        fig_width=fig_width,
        fig_track_height=1.0,
        align_type="center",
        tick_style=None,
        link_track_ratio=5.0,
    )
    name2maxsize = pmauve.get_name2maxsize(list(names), bbone)
    name2blocks = pmauve.get_name2blocks(list(alignments))
    for name in names:
        blocks = name2blocks.get(name, [])
        track = gv.add_feature_track(name, name2maxsize[name], labelsize=20)
        for (start, end, strand), color in zip(
            blocks, _spaced_colors(cmap, len(blocks)), strict=True
        ):
            track.add_feature(start, end, strand, plotstyle="box", facecolor=color)
    for ac in alignments:
        gv.add_link(ac.ref_link, ac.query_link, _NORMAL_LINK, "tomato", curve=True)
    gv.savefig(image, dpi=dpi)


def _parse_rbh(path: Path, ref_name: str, query_name: str) -> list[Any]:
    """Alignment coordinates from an ``mmseqs easy-rbh`` result (identity in %, no duplicates)."""
    align = _require_isolated("pygenomeviz.align", _EXTRA)
    coords: list[Any] = []
    seen: list[str] = []
    with path.open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            ref_start, ref_end, ref_strand = (int(v) for v in row[0].split("|")[-2].split("_"))
            ref_length = ref_end - ref_start + 1
            if ref_strand == -1:
                ref_start, ref_end = ref_end, ref_start
            query_start, query_end, query_strand = (
                int(v) for v in row[1].split("|")[-2].split("_")
            )
            query_length = query_end - query_start + 1
            if query_strand == -1:
                query_start, query_end = query_end, query_start
            identity = int(float(row[2]) * 10000) / 100
            ref_key = f"{ref_name} {ref_start} {ref_end} {ref_strand}"
            query_key = f"{query_name} {query_start} {query_end} {query_strand}"
            if ref_key in seen or query_key in seen:
                continue
            seen.extend([ref_key, query_key])
            coords.append(
                align.AlignCoord(
                    ref_start,
                    ref_end,
                    query_start,
                    query_end,
                    ref_length,
                    query_length,
                    identity,
                    ref_name,
                    query_name,
                )
            )
    return coords


def _parse_bbone(path: Path, names: Sequence[str], refid: int) -> list[Any]:
    """Conserved blocks of a progressiveMauve backbone file (pygenomeviz parse as a function)."""
    align = _require_isolated("pygenomeviz.align", _EXTRA)
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genome_num = len(header) // 2
        if genome_num != len(names):
            raise BackendExecutionError(
                f"progressiveMauve backbone lists {genome_num} genomes, expected {len(names)}.",
                context={"backbone": str(path)},
            )
        ref_idx = refid * 2
        rows: list[list[int]] = []
        for raw in reader:
            if not raw:
                continue
            row = [int(col) for col in raw]
            if row[ref_idx] < 0:
                row = [col * -1 for col in row]
            if row.count(0) >= 2:
                continue
            if abs(row[ref_idx] - row[ref_idx + 1]) < 20:
                continue
            rows.append(row)
    rows.sort(key=lambda row: row[ref_idx])
    coords: list[Any] = []
    for row in rows:
        for index in range(genome_num - 1):
            idx = index * 2
            rstart, rend = row[idx], row[idx + 1]
            if rstart < 0 and rend < 0:
                rstart, rend = abs(rend), abs(rstart)
            qstart, qend = row[idx + 2], row[idx + 3]
            if qstart < 0 and qend < 0:
                qstart, qend = abs(qend), abs(qstart)
            coords.append(
                align.AlignCoord(
                    rstart,
                    rend,
                    qstart,
                    qend,
                    abs(rend - rstart) + 1,
                    abs(qend - qstart) + 1,
                    0,
                    names[index],
                    names[index + 1],
                )
            )
    return coords


def _finish_alignments(
    folder: Path, coords: Sequence[Any], *, min_length: int, min_identity: float
) -> list[Any]:
    """Write every alignment to align_coords.tsv; return those passing the drawing thresholds."""
    align = _require_isolated("pygenomeviz.align", _EXTRA)
    align.AlignCoord.write(list(coords), folder / "align_coords.tsv")
    filtered: list[Any] = align.AlignCoord.filter(list(coords), min_length, min_identity)
    return filtered


def _shift(coords: Sequence[Any], offsets: Mapping[str, int]) -> None:
    """Move range-relative coordinates to genome coordinates; unknown names are an output error."""
    for ac in coords:
        unknown = [name for name in (ac.ref_name, ac.query_name) if name not in offsets]
        if unknown:
            raise BackendExecutionError(
                f"show-coords reported sequence names that were not submitted: {unknown}.",
                context={"submitted": sorted(offsets)},
            )
        ac.ref_start += offsets[ac.ref_name]
        ac.ref_end += offsets[ac.ref_name]
        ac.query_start += offsets[ac.query_name]
        ac.query_end += offsets[ac.query_name]


def _failed(steps: _Steps) -> dict[str, Any]:
    return {
        "status": "failed",
        "backend": steps.backend,
        "commands": steps.commands,
        "return_codes": steps.return_codes,
        "outputs": [],
        "image_path": None,
        "log_tail": steps.tail(),
        "data": {},
    }


def _completed(
    steps: _Steps, outputs: Sequence[Path], image: Path, data: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "status": "completed",
        "backend": steps.backend,
        "commands": steps.commands,
        "return_codes": steps.return_codes,
        "outputs": [str(path) for path in outputs],
        "image_path": str(image),
        "log_tail": steps.tail(),
        "data": dict(data),
    }


def _alignment_data(genomes: Sequence[Any], alignments: Sequence[Any]) -> dict[str, Any]:
    return {
        "genomes": [_genome_summary(gbk) for gbk in genomes],
        "alignment_count": len(alignments),
        "alignments": [asdict(ac) for ac in alignments],
    }


def load_example_dataset(
    name: str,
    *,
    cache_dir: str | None = None,
    overwrite_cache: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a pyGenomeViz example GenBank dataset from GitHub and list its files and links."""
    utils = _require_isolated("pygenomeviz.utils", _EXTRA)
    config = _require_isolated("pygenomeviz.config", _EXTRA)
    if name not in config.DATASETS:
        raise ConfigurationError(
            f"Unknown pyGenomeViz example dataset {name!r}; "
            f"choose one of {sorted(config.DATASETS)}."
        )
    cache = _cache(cache_dir)
    gbk_files, links = utils.load_example_dataset(
        name, cache_dir=cache, overwrite_cache=overwrite_cache
    )
    return {
        "status": "completed",
        "service": _SERVICE,
        "source_url": f"{config.GITHUB_DATA_URL}{name}/",
        "dataset": name,
        "request": {"name": name, "cache_dir": cache, "overwrite_cache": overwrite_cache},
        "retrieved_at_utc": _utc_now(),
        "gbk_files": [str(path) for path in gbk_files],
        "link_count": len(links),
        "links": [asdict(link) for link in links],
    }


def load_example_gff(
    filename: str,
    *,
    cache_dir: str | None = None,
    overwrite_cache: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a pyGenomeViz example GFF file from GitHub into a cache folder; return its path."""
    utils = _require_isolated("pygenomeviz.utils", _EXTRA)
    config = _require_isolated("pygenomeviz.config", _EXTRA)
    if filename not in config.EXAMPLE_GFF_FILES:
        raise ConfigurationError(
            f"Unknown pyGenomeViz example GFF {filename!r}; "
            f"choose one of {list(config.EXAMPLE_GFF_FILES)}."
        )
    cache = _cache(cache_dir)
    gff_path = utils.load_example_gff(filename, cache_dir=cache, overwrite_cache=overwrite_cache)
    return {
        "status": "completed",
        "service": _SERVICE,
        "source_url": f"{config.GITHUB_DATA_URL}gff/{filename}",
        "filename": filename,
        "request": {"filename": filename, "cache_dir": cache, "overwrite_cache": overwrite_cache},
        "retrieved_at_utc": _utc_now(),
        "gff_path": str(Path(gff_path).expanduser().resolve()),
    }


def plot_genome_features(
    gbk_resources: Sequence[str],
    image_path: str,
    *,
    fig_width: float = 15.0,
    fig_track_height: float = 1.0,
    feature_color: str = "orange",
    feature_plotstyle: str = "bigarrow",
    align_type: str = "left",
    tick_style: str | None = None,
    dpi: int = 300,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the CDS features of GenBank genomes as stacked pyGenomeViz tracks in one PNG."""
    genomes = _genomes(gbk_resources, minimum=1)
    target = _new_png(image_path)
    _positive_number(fig_width, "fig_width")
    _positive_number(fig_track_height, "fig_track_height")
    _choice(feature_plotstyle, _PLOTSTYLES, "feature_plotstyle")
    _choice(align_type, _ALIGN_TYPES, "align_type")
    if tick_style is not None:
        _choice(tick_style, _TICK_STYLES, "tick_style")
    _positive_int(dpi, "dpi")
    pgv = _require_isolated("pygenomeviz", _EXTRA)
    gv = _genome_viz(
        pgv,
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        align_type=align_type,
        tick_style=tick_style,
        link_track_ratio=5.0,
    )
    for gbk in genomes:
        track = gv.add_feature_track(
            gbk.name, size=gbk.range_size, start_pos=gbk.min_range, labelsize=20
        )
        _draw_cds(track, gbk, plotstyle=feature_plotstyle, facecolor=feature_color, linewidth=0.0)
    gv.savefig(target, dpi=dpi)
    return {
        "status": "completed",
        "image_path": str(target),
        "data": {"genomes": [_genome_summary(gbk) for gbk in genomes]},
    }


def compare_genomes_mmseqs(
    gbk_resources: Sequence[str],
    output_dir: str,
    *,
    evalue: float = 1e-3,
    min_identity: float = 0.0,
    threads: int = 1,
    fig_width: float = 15.0,
    dpi: int = 300,
    timeout_seconds: float = 3600.0,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Link neighbouring genomes by MMseqs2 reciprocal best protein hits; draw result.png.

    Each pair of consecutive CDS protein sets goes through ``mmseqs easy-rbh`` with E-value
    ``evalue``. Hits with identity below ``min_identity`` (percent) are not drawn; all hits are kept
    in align_coords.tsv.
    """
    genomes = _genomes(gbk_resources, minimum=2)
    _positive_number(evalue, "evalue")
    _percent(min_identity, "min_identity")
    _positive_int(threads, "threads")
    _positive_number(fig_width, "fig_width")
    _positive_int(dpi, "dpi")
    folder = _check_new_folder(output_dir)
    program = _executable("mmseqs", mmseqs_path, _MMSEQS_HINT)
    folder.mkdir(parents=True, exist_ok=True)
    proteins: list[Path] = []
    for gbk in genomes:
        path = folder / f"{gbk.name}.faa"
        _write_cds_protein_fasta(gbk, path)
        proteins.append(path)
    steps = _Steps("mmseqs")
    rbh_files: list[Path] = []
    coords: list[Any] = []
    for index in range(len(genomes) - 1):
        ref_name, query_name = genomes[index].name, genomes[index + 1].name
        rbh = folder / f"{index + 1:02d}_{ref_name}-{query_name}_rbh.tsv"
        work = folder / f"mmseqs_work_{index + 1:02d}"
        work.mkdir()
        args = [
            "easy-rbh",
            str(proteins[index]),
            str(proteins[index + 1]),
            str(rbh),
            str(work),
            "--threads",
            str(threads),
            "-e",
            str(evalue),
            "-v",
            "0",
        ]
        code, _ = steps.run(program, args, cwd=folder, timeout_seconds=timeout_seconds)
        shutil.rmtree(work, ignore_errors=True)
        if code != 0 or not rbh.is_file():
            return _failed(steps)
        rbh_files.append(rbh)
        coords.extend(_parse_rbh(rbh, ref_name, query_name))
    alignments = _finish_alignments(folder, coords, min_length=0, min_identity=min_identity)
    image = folder / "result.png"
    _draw_feature_alignment(
        _require_isolated("pygenomeviz", _EXTRA),
        genomes,
        alignments,
        image,
        fig_width=fig_width,
        dpi=dpi,
    )
    outputs = [*rbh_files, folder / "align_coords.tsv", image]
    return _completed(steps, outputs, image, _alignment_data(genomes, alignments))


def compare_genomes_mummer(
    gbk_resources: Sequence[str],
    output_dir: str,
    *,
    seqtype: str = "protein",
    min_length: int = 0,
    min_identity: float = 0.0,
    one_to_one: bool = False,
    fig_width: float = 15.0,
    dpi: int = 300,
    timeout_seconds: float = 3600.0,
    mummer_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align neighbouring genomes with MUMmer and draw the links to output_dir/result.png.

    ``seqtype='protein'`` runs promer (translated comparison, the pyGenomeViz default); the
    ``'nucleotide'`` choice runs nucmer. The steps are ``--mum`` alignment, ``delta-filter -m``
    (``-1`` with ``one_to_one``) and ``show-coords -H -T -r -k``. Alignments shorter than
    ``min_length`` or with identity below ``min_identity`` (percent) are not drawn.
    """
    genomes = _genomes(gbk_resources, minimum=2)
    _choice(seqtype, _SEQTYPES, "seqtype")
    if isinstance(min_length, bool) or not isinstance(min_length, int) or min_length < 0:
        raise ConfigurationError("min_length must be a non-negative integer.")
    _percent(min_identity, "min_identity")
    _positive_number(fig_width, "fig_width")
    _positive_int(dpi, "dpi")
    folder = _check_new_folder(output_dir)
    aligner = "promer" if seqtype == "protein" else "nucmer"

    def located(name: str) -> str | None:
        return None if mummer_dir is None else str(Path(mummer_dir).expanduser() / name)

    programs = {
        name: _executable(name, located(name), _MUMMER_HINT)
        for name in (aligner, "delta-filter", "show-coords")
    }
    folder.mkdir(parents=True, exist_ok=True)
    fasta_files: list[Path] = []
    for gbk in genomes:
        path = folder / f"{gbk.name}.fna"
        gbk.write_genome_fasta(path)
        fasta_files.append(path)
    align = _require_isolated("pygenomeviz.align", _EXTRA)
    offsets = {gbk.name: int(gbk.min_range) for gbk in genomes}
    steps = _Steps("mummer")
    coords: list[Any] = []
    for index in range(len(genomes) - 1):
        prefix = folder / f"out{index}"
        delta = folder / f"out{index}.delta"
        code, _ = steps.run(
            programs[aligner],
            ["--mum", str(fasta_files[index]), str(fasta_files[index + 1]), f"--prefix={prefix}"],
            cwd=folder,
            timeout_seconds=timeout_seconds,
        )
        if code != 0 or not delta.is_file():
            return _failed(steps)
        code, filtered_text = steps.run(
            programs["delta-filter"],
            ["-1" if one_to_one else "-m", str(delta)],
            cwd=folder,
            timeout_seconds=timeout_seconds,
            max_output_bytes=_LARGE_OUTPUT,
        )
        if code != 0:
            return _failed(steps)
        filtered = folder / f"filter_out{index}.delta"
        filtered.write_text(filtered_text)
        code, coords_text = steps.run(
            programs["show-coords"],
            ["-H", "-T", "-r", "-k", str(filtered)],
            cwd=folder,
            timeout_seconds=timeout_seconds,
            max_output_bytes=_LARGE_OUTPUT,
        )
        if code != 0:
            return _failed(steps)
        coords_file = folder / f"coords{index}.tsv"
        rows = [
            line for line in coords_text.splitlines() if line.split("\t", 1)[0].strip().isdigit()
        ]
        coords_file.write_text("".join(f"{line}\n" for line in rows))
        pair = align.MUMmer.parse_coords_file(coords_file, seqtype)
        _shift(pair, offsets)
        coords.extend(pair)
        for work in (delta, filtered, coords_file):
            work.unlink(missing_ok=True)
    alignments = _finish_alignments(
        folder, coords, min_length=min_length, min_identity=min_identity
    )
    image = folder / "result.png"
    _draw_feature_alignment(
        _require_isolated("pygenomeviz", _EXTRA),
        genomes,
        alignments,
        image,
        fig_width=fig_width,
        dpi=dpi,
    )
    outputs = [folder / "align_coords.tsv", image]
    return _completed(steps, outputs, image, _alignment_data(genomes, alignments))


def compare_genomes_progressive_mauve(
    seq_files: Sequence[str],
    output_dir: str,
    *,
    refid: int = 0,
    cmap: str = "hsv",
    fig_width: float = 15.0,
    dpi: int = 300,
    timeout_seconds: float = 3600.0,
    progressive_mauve_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align sequence files with progressiveMauve; draw blocks and links to output_dir/result.png.

    Inputs are FASTA or GenBank files; a genome is named after its file without the suffix. Blocks
    shorter than 20 bp are dropped, as in pyGenomeViz. ``refid`` is the reference whose coordinates
    set the block order.
    """
    files = _sequence_files(seq_files, minimum=2)
    names = [path.with_suffix("").name for path in files]
    _check_names(names)
    refid = int(refid)
    if not 0 <= refid < len(files):
        raise ConfigurationError(f"refid must be between 0 and {len(files) - 1}.")
    mpl = require("matplotlib", _EXTRA)
    if cmap not in mpl.colormaps:
        raise ConfigurationError(f"cmap {cmap!r} is not a matplotlib colormap name.")
    _positive_number(fig_width, "fig_width")
    _positive_int(dpi, "dpi")
    folder = _check_new_folder(output_dir)
    program = _executable("progressiveMauve", progressive_mauve_path, _MAUVE_HINT)
    folder.mkdir(parents=True, exist_ok=True)
    seq_dir = folder / "seqfiles"
    seq_dir.mkdir()
    copies: list[Path] = []
    for source in files:
        target = seq_dir / source.name
        if target.suffix == ".gbff":
            target = target.with_suffix(".gbk")
        shutil.copy(source, target)
        copies.append(target)
    xmfa = folder / "mauve.xmfa"
    bbone = folder / "mauve_bbone.tsv"
    steps = _Steps("progressiveMauve")
    code, _ = steps.run(
        program,
        [f"--output={xmfa}", f"--backbone-output={bbone}", *[str(path) for path in copies]],
        cwd=folder,
        timeout_seconds=timeout_seconds,
    )
    if code != 0 or not bbone.is_file():
        return _failed(steps)
    alignments = _parse_bbone(bbone, names, refid)
    align = _require_isolated("pygenomeviz.align", _EXTRA)
    align.AlignCoord.write(alignments, folder / "align_coords.tsv")
    image = folder / "result.png"
    _draw_block_alignment(
        _require_isolated("pygenomeviz", _EXTRA),
        names,
        bbone,
        alignments,
        image,
        fig_width=fig_width,
        dpi=dpi,
        cmap=cmap,
    )
    data = {
        "genomes": [
            {"name": name, "sequence_file": str(path)}
            for name, path in zip(names, files, strict=True)
        ],
        "alignment_count": len(alignments),
        "alignments": [asdict(ac) for ac in alignments],
    }
    return _completed(steps, [xmfa, bbone, folder / "align_coords.tsv", image], image, data)


__all__ = [
    "compare_genomes_mmseqs",
    "compare_genomes_mummer",
    "compare_genomes_progressive_mauve",
    "load_example_dataset",
    "load_example_gff",
    "plot_genome_features",
]
_register_public_verbose(globals())
