"""Bakta command-line facades (pack P11): fixed argument lists, no shell, bounded execution.

``annotate_bakta_features`` runs ``bakta`` (genome annotation with a choice of the feature types to detect),
``export_bakta_annotation`` runs ``bakta_io`` (re-writes a Bakta JSON annotation in every output format) and
``plot_bakta_genome`` runs ``bakta_plot`` (circular genome plots from a Bakta JSON annotation). Bakta 1.12.2 is the
tested version. The user installs Bakta (``conda install -c bioconda bakta`` or ``python -m pip install bakta``), its
external tools (tRNAscan-SE, ARAGORN, Infernal, PILER-CR, DIAMOND, BLAST+, AMRFinderPlus) and downloads a database
(``bakta_db download``); programs are put on PATH or passed as ``bakta_path``, ``bakta_io_path`` and ``bakta_plot_path``.
Every run writes into a fresh ``output_directory`` (never overwritten); the Bakta result files are in its
``annotation``, ``export`` or ``plots`` subfolder and ``tool.log`` holds the program's standard output and error. Results
are JSON-serialisable dictionaries with status, version, command, output files, small parsed results and a bounded log
tail. Optional parameters default to ``None`` and are then not passed, so Bakta's own defaults apply.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p11_cli as _p
from dnakit.exceptions import ConfigurationError

_TESTED_VERSION = "1.12.2"
_CITATION = "https://doi.org/10.1099/mgen.0.000685"
_HINT = (
    "Install Bakta 1.12 or later (for example `conda install -c bioconda bakta`, or `python -m pip install bakta`) with its "
    "external tools, download a database with `bakta_db download`, then put bakta on PATH or pass bakta_path."
)
_IO_HINT = "bakta_io is installed with Bakta (`python -m pip install bakta`); put it on PATH or pass bakta_io_path."
_PLOT_HINT = "bakta_plot is installed with Bakta (`python -m pip install bakta`); put it on PATH or pass bakta_plot_path."
# Feature name used by annotate_bakta_features -> the Bakta workflow switch that skips it.
_FEATURE_SWITCHES = {
    "cds": "--skip-cds",
    "pseudogenes": "--skip-pseudo",
    "small_orfs": "--skip-sorf",
    "trna": "--skip-trna",
    "tmrna": "--skip-tmrna",
    "rrna": "--skip-rrna",
    "ncrna": "--skip-ncrna",
    "ncrna_regions": "--skip-ncrna-region",
    "crispr": "--skip-crispr",
    "gaps": "--skip-gap",
    "origins": "--skip-ori",
}
_MAX_JSON_BYTES = 500_000_000
_GRAM = ("+", "-", "?")


def _annotation_files(directory: Path, prefix: str) -> dict[str, Path]:
    """Every file ``<prefix>.<extension>`` that Bakta wrote, keyed by the extension."""
    found: dict[str, Path] = {}
    for path in sorted(directory.glob(f"{prefix}.*")):
        if path.is_file():
            found[path.name[len(prefix) + 1 :]] = path
    return found


def _feature_counts(json_path: Path) -> dict[str, int] | None:
    """Number of annotated features per Bakta feature type, from the JSON annotation."""
    if json_path.stat().st_size > _MAX_JSON_BYTES:
        return None
    with json_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    counts: dict[str, int] = {}
    for feature in data.get("features", []):
        kind = str(feature.get("type", "unknown"))
        counts[kind] = counts.get(kind, 0) + 1
        if kind == "cds" and "pseudogene" in feature:
            counts["pseudogene"] = counts.get("pseudogene", 0) + 1
        if kind == "cds" and "hypothetical" in feature:
            counts["hypothetical"] = counts.get("hypothetical", 0) + 1
    return counts


def _summary(path: Path) -> dict[str, Any]:
    """The ``key: value`` lines of the Bakta summary report (``<prefix>.txt``)."""
    values: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        key, separator, value = line.partition(": ")
        if separator and key and not key.startswith(("DOI", "URL")):
            values[key.strip()] = _p.number(value.strip())
    return values


def annotate_bakta_features(
    genome: str | Path,
    database: str | Path,
    output_directory: str | Path,
    *,
    features: Sequence[str] | None = None,
    prefix: str = "annotation",
    min_contig_length: int | None = None,
    genus: str | None = None,
    species: str | None = None,
    strain: str | None = None,
    plasmid: str | None = None,
    complete: bool = False,
    prodigal_training_file: str | Path | None = None,
    translation_table: int | None = None,
    gram: str | None = None,
    locus: str | None = None,
    locus_tag: str | None = None,
    locus_tag_increment: int | None = None,
    keep_contig_headers: bool = False,
    compliant: bool = False,
    replicons: str | Path | None = None,
    regions: str | Path | None = None,
    trusted_proteins: str | Path | None = None,
    trusted_hmms: str | Path | None = None,
    metagenome: bool = False,
    partial: bool = False,
    skip_overlap_filters: bool = False,
    plot: bool = False,
    threads: int | None = None,
    bakta_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate a bacterial genome, MAG or plasmid with Bakta, choosing which feature types are detected.

    ``features`` selects the detection steps to run (default: all, as Bakta does): ``cds`` (protein-coding genes with
    their functional annotation: UniRef100/UniParc/PSC/PSCC database look-ups, AMRFinderPlus resistance genes and the
    curated expert protein sequences matched with DIAMOND), ``pseudogenes``, ``small_orfs`` (sORFs below 30 aa),
    ``trna`` (tRNAscan-SE), ``tmrna`` (ARAGORN), ``rrna`` and ``ncrna`` and ``ncrna_regions`` (Infernal cmscan against
    Rfam), ``crispr`` (PILER-CR), ``gaps`` (runs of N in the assembly) and ``origins`` (oriC/oriV/oriT, BLASTn against
    the origin databases). Steps that are not listed are skipped (``--skip-*``); ``pseudogenes`` needs ``cds``.
    ``database`` is a Bakta database folder. Further options are Bakta's own: organism (``genus``, ``species``,
    ``strain``, ``plasmid``), ``complete`` replicons, ``translation_table`` (11, 4 or 25), ``gram`` (``+``, ``-`` or
    ``?`` for signal peptide prediction), locus naming (``locus``, ``locus_tag``, ``locus_tag_increment``),
    ``compliant`` INSDC output, ``replicons`` table, pre-annotated ``regions`` (GFF3 or GenBank), ``trusted_proteins``
    (FASTA) and ``trusted_hmms`` for CDS annotation, ``metagenome`` mode, ``partial`` genes at sequence ends and
    ``skip_overlap_filters`` (``--skip-filter``). ``plot`` also draws the circular genome plots (off by default; use
    ``plot_bakta_genome`` for them). Result keys: ``annotation_directory``, ``outputs`` (keyed by extension: ``json``,
    ``tsv``, ``gff3``, ``gbff``, ``embl``, ``fna``, ``ffn``, ``faa``, ``hypotheticals.tsv``, ``hypotheticals.faa``,
    ``inference.tsv``, ``txt``), ``features`` (the steps run), ``feature_counts`` (per type, from the JSON) and
    ``summary`` (the report in ``<prefix>.txt``).
    """
    genome_file = _p.file_path(genome, "genome")
    database_path = Path(database).expanduser().resolve()
    if not database_path.is_dir():
        raise ConfigurationError(f"database must be an existing Bakta database folder: {database_path}")
    if features is None:
        selected = list(_FEATURE_SWITCHES)
    elif isinstance(features, str) or not features:
        raise ConfigurationError("features must be a non-empty list of feature names.")
    else:
        selected = [_cli.choice(item, tuple(_FEATURE_SWITCHES), "features") for item in features]
    if "pseudogenes" in selected and "cds" not in selected:
        raise ConfigurationError("The pseudogenes step needs the cds step.")
    stem = _p.file_stem(prefix, "prefix")
    arguments = ["--db", str(database_path), "--prefix", stem]
    _p.opt_int(arguments, "--min-contig-length", min_contig_length, "min_contig_length", minimum=1)
    _p.opt_text(arguments, "--genus", genus, "genus")
    _p.opt_text(arguments, "--species", species, "species")
    _p.opt_text(arguments, "--strain", strain, "strain")
    _p.opt_text(arguments, "--plasmid", plasmid, "plasmid")
    _p.switch(arguments, "--complete", complete, "complete")
    _p.opt_file(arguments, "--prodigal-tf", prodigal_training_file, "prodigal_training_file")
    if translation_table is not None:
        arguments.extend(("--translation-table", str(_cli.choice(translation_table, (4, 11, 25), "translation_table"))))  # type: ignore[arg-type]
    if gram is not None:
        arguments.extend(("--gram", _cli.choice(gram, _GRAM, "gram")))
    _p.opt_text(arguments, "--locus", locus, "locus")
    _p.opt_text(arguments, "--locus-tag", locus_tag, "locus_tag")
    if locus_tag_increment is not None:
        arguments.extend(("--locus-tag-increment", str(_cli.choice(locus_tag_increment, (1, 5, 10), "locus_tag_increment"))))  # type: ignore[arg-type]
    _p.switch(arguments, "--keep-contig-headers", keep_contig_headers, "keep_contig_headers")
    _p.switch(arguments, "--compliant", compliant, "compliant")
    _p.opt_file(arguments, "--replicons", replicons, "replicons")
    _p.opt_file(arguments, "--regions", regions, "regions")
    _p.opt_file(arguments, "--proteins", trusted_proteins, "trusted_proteins")
    _p.opt_file(arguments, "--hmms", trusted_hmms, "trusted_hmms")
    _p.switch(arguments, "--meta", metagenome, "metagenome")
    _p.switch(arguments, "--partial", partial, "partial")
    _p.switch(arguments, "--skip-filter", skip_overlap_filters, "skip_overlap_filters")
    arguments.extend(switch for name, switch in _FEATURE_SWITCHES.items() if name not in selected)
    if not _p.flag(plot, "plot"):
        arguments.append("--skip-plot")
    _p.opt_int(arguments, "--threads", threads, "threads", minimum=1, maximum=1024)
    executable = _cli.program("bakta", explicit=bakta_path, hint=_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    annotation = directory / "annotation"
    command = [*arguments, "--output", str(annotation), genome_file]
    _p.run_checked("annotate_bakta_features", executable, command, directory, timeout_seconds=timeout_seconds)
    files = _annotation_files(annotation, stem)
    _p.require_outputs("annotate_bakta_features", directory, {name: annotation / f"{stem}.{name}" for name in ("json", "tsv", "gff3")})
    return _p.result(
        "annotate_bakta_features",
        "bakta",
        _p.version_of(executable),
        ["bakta", *command],
        directory,
        files,
        tool="Bakta",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        annotation_directory=str(annotation),
        features=selected,
        feature_counts=_feature_counts(files["json"]),
        summary=_summary(files["txt"]) if "txt" in files else {},
    )


def export_bakta_annotation(
    annotation_json: str | Path,
    output_directory: str | Path,
    *,
    prefix: str | None = None,
    bakta_io_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a Bakta JSON annotation in all output formats (``bakta_io``), without re-running the annotation.

    ``annotation_json`` is the ``.json`` file written by ``annotate_bakta_features`` (or by ``bakta``). The export
    contains the human readable feature table (``tsv``), the feature inference table (``inference.tsv``), the
    hypothetical CDS table (``hypotheticals.tsv``), GenBank (``gbff``), EMBL (``embl``) and GFF3 (``gff3``) flat files,
    the genome sequences (``fna``), the CDS/feature nucleotide sequences (``ffn``), the translated CDS sequences
    (``faa``), the hypothetical CDS proteins (``hypotheticals.faa``), the circular genome plot (``png`` and ``svg``)
    and the genome and annotation summary report (``txt``). ``prefix`` defaults to the stem of the input file. Result
    keys: ``export_directory``, ``outputs`` (keyed by extension) and ``summary`` (the report).
    """
    source = Path(_p.file_path(annotation_json, "annotation_json"))
    stem = _p.file_stem(prefix, "prefix") if prefix is not None else _p.file_stem(source.stem, "annotation_json")
    executable = _cli.program("bakta_io", explicit=bakta_io_path, hint=_IO_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    export = directory / "export"
    command = ["--prefix", stem, "--output", str(export), str(source)]
    _p.run_checked("export_bakta_annotation", executable, command, directory, timeout_seconds=timeout_seconds)
    files = _annotation_files(export, stem)
    _p.require_outputs("export_bakta_annotation", directory, {name: export / f"{stem}.{name}" for name in ("tsv", "gff3", "gbff", "embl", "txt")})
    return _p.result(
        "export_bakta_annotation",
        "bakta_io",
        _p.version_of(executable),
        ["bakta_io", *command],
        directory,
        files,
        tool="Bakta",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        export_directory=str(export),
        summary=_summary(files["txt"]),
    )


def plot_bakta_genome(
    annotation_json: str | Path,
    output_directory: str | Path,
    *,
    plot_type: str = "features",
    prefix: str | None = None,
    sequences: Sequence[str] | None = None,
    config: str | Path | None = None,
    label: str | None = None,
    size: int | None = None,
    dpi: int | None = None,
    bakta_plot_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw circular genome maps from a Bakta JSON annotation (``bakta_plot``).

    ``plot_type`` ``features`` (default) colours features by type; ``cog`` colours CDS by COG functional category.
    ``sequences`` selects the sequences to plot by one-based number or name (default all); ``config`` is a YAML plot
    configuration; ``label`` is the centre label (``|`` makes a line break); ``size`` is 4, 8 or 16 inches (Bakta
    default 8) and ``dpi`` 150, 300 or 600 (Bakta default 300). The images are ``outputs['png']`` and
    ``outputs['svg']``. Result keys: ``plot_directory`` and ``plot_type``.
    """
    source = Path(_p.file_path(annotation_json, "annotation_json"))
    kind = _cli.choice(plot_type, ("features", "cog"), "plot_type")
    stem = _p.file_stem(prefix, "prefix") if prefix is not None else _p.file_stem(source.stem, "annotation_json")
    arguments = ["--type", kind, "--prefix", stem]
    if sequences is not None:
        if isinstance(sequences, str) or not sequences:
            raise ConfigurationError("sequences must be a non-empty list of sequence numbers or names.")
        names = [_cli.text(str(item), "sequences") for item in sequences]
        if any("," in item for item in names):
            raise ConfigurationError("sequence names must not contain commas.")
        arguments.extend(("--sequences", ",".join(names)))
    _p.opt_file(arguments, "--config", config, "config")
    _p.opt_text(arguments, "--label", label, "label")
    if size is not None:
        arguments.extend(("--size", str(_cli.choice(size, (4, 8, 16), "size"))))  # type: ignore[arg-type]
    if dpi is not None:
        arguments.extend(("--dpi", str(_cli.choice(dpi, (150, 300, 600), "dpi"))))  # type: ignore[arg-type]
    executable = _cli.program("bakta_plot", explicit=bakta_plot_path, hint=_PLOT_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    plots = directory / "plots"
    command = [*arguments, "--output", str(plots), str(source)]
    _p.run_checked("plot_bakta_genome", executable, command, directory, timeout_seconds=timeout_seconds)
    files = _annotation_files(plots, stem)
    _p.require_outputs("plot_bakta_genome", directory, {name: plots / f"{stem}.{name}" for name in ("png", "svg")})
    return _p.result(
        "plot_bakta_genome",
        "bakta_plot",
        _p.version_of(executable),
        ["bakta_plot", *command],
        directory,
        files,
        tool="Bakta",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        plot_directory=str(plots),
        plot_type=kind,
    )


__all__ = ["annotate_bakta_features", "export_bakta_annotation", "plot_bakta_genome"]
_register_public_verbose(globals())
