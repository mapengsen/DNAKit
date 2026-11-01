"""GEMINI command-line facades: queries, exports, statistics, inheritance models and gene-level tools on a GEMINI database.

GEMINI (https://github.com/arq5x/gemini) stores the variants of a VCF in a SQLite database and offers sub-commands on it. The
user installs it (``conda install -c bioconda gemini``, 0.30.2 is the tested version), loads a database with
``dnakit.variants.load_gemini_database`` and puts ``gemini`` on PATH or passes ``gemini_path``. Each public function runs one
GEMINI sub-command with an explicit typed signature mapped to a fixed set of flags: no shell, no generic "run any command".
Sub-commands that print a table write it to a new TSV file (standard output streamed to the file, bounded size); a missing
program raises ``BackendUnavailableError`` with an install hint. Sub-commands that change the database (amend, annotate,
set_somatic, merge_chunks) work on a copy in a new path (``output_database``) unless ``in_place=True`` is given.
Results are JSON-serialisable dictionaries with status, backend, version, command, return code, existing outputs, a bounded
log tail and, for tables, a row count and a short preview. Options left at ``None`` or ``False`` are not passed, so GEMINI's
own defaults apply (named in the docstrings where they matter). The ``--annotation-dir`` option is exposed as
``annotation_dir`` and the annotation files themselves are the caller's responsibility (nothing is downloaded).
"""

from __future__ import annotations

import re
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p10_cli as _p10
from dnakit.exceptions import ConfigurationError

_PROGRAM = "gemini"
_TESTED_VERSION = "0.30.2"
_CITATION = "https://doi.org/10.1371/journal.pcbi.1003153"
_HINT = (
    "Install GEMINI 0.30.2 or later (`conda install -c bioconda gemini`), then put gemini on PATH or pass gemini_path. "
    "The annotation files (`gemini update --dataonly`) are installed separately."
)
_COLUMNS = re.compile(r"[A-Za-z0-9_.,*() ]+")
_FAMILIES = re.compile(r"[A-Za-z0-9_.,-]+")
_NAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*")


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


class _Args:
    """GEMINI arguments in order; ``values`` skips ``None`` and ``switch`` adds a flag only when it is on."""

    def __init__(self, *items: str) -> None:
        self.items: list[str] = list(items)

    def switch(self, flag: str, on: bool, name: str) -> None:
        if not isinstance(on, bool):
            raise ConfigurationError(f"{name} must be True or False.")
        if on:
            self.items.append(flag)

    def value(self, flag: str, value: object | None) -> None:
        if value is not None:
            self.items += [flag, str(value)]

    def many(self, flag: str, values: Sequence[str] | None) -> None:
        for value in values or ():
            self.items += [flag, value]


def _sql(value: str, name: str) -> str:
    """One SQL or filter expression (a single argument; never a shell string)."""
    if not isinstance(value, str) or not value.strip() or len(value) > 32_768 or "\x00" in value:
        raise ConfigurationError(
            f"{name} must be non-empty text of at most 32768 characters without NUL."
        )
    return value


def _optional_sql(value: str | None, name: str) -> str | None:
    return None if value is None else _sql(value, name)


def _list(values: Sequence[str] | None, name: str) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    return [_cli.text(value, name) for value in values]


def _columns(value: str | None, name: str = "columns") -> str | None:
    if value is None:
        return None
    text = _cli.text(value, name)
    if _COLUMNS.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} must be a comma-separated list of column names (letters, digits, '_')."
        )
    return text


def _optional_int(
    value: int | None, name: str, *, minimum: int = 0, maximum: int | None = None
) -> int | None:
    return None if value is None else _cli.integer(value, name, minimum=minimum, maximum=maximum)


def _optional_number(
    value: float | None, name: str, *, minimum: float | None = None, maximum: float | None = None
) -> float | None:
    return None if value is None else _cli.number(value, name, minimum=minimum, maximum=maximum)


def _database(path: str | Path) -> str:
    return str(_cli.existing_file(path, "database"))


def _global(annotation_dir: str | None) -> list[str]:
    if annotation_dir is None:
        return []
    directory = Path(annotation_dir).expanduser().resolve()
    if not directory.is_dir():
        raise ConfigurationError(f"annotation_dir must be an existing directory: {directory}")
    return ["--annotation-dir", str(directory)]


def _check_executable(gemini_path: str | None) -> str:
    """Resolve GEMINI before anything is written, so a missing program leaves no output behind."""
    return _cli.program(_PROGRAM, explicit=gemini_path, hint=_HINT)


def _finish(
    operation: str,
    rendered: dict[str, Any],
    parse: Callable[[], dict[str, Any]] | None = None,
    **extra: Any,
) -> dict[str, Any]:
    result = dict(rendered)
    result["operation"] = operation
    result["provenance"] = _p10.provenance(
        "GEMINI", _TESTED_VERSION, rendered.get("version"), _CITATION
    )
    if result["status"] == "completed":
        result.update(extra)
        if parse is not None:
            result.update(parse())
    return result


def _table_command(
    operation: str,
    arguments: Sequence[str],
    output: Path,
    *,
    gemini_path: str | None,
    header: bool,
    annotation_dir: str | None,
    environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Run ``gemini <annotation-dir> <arguments>`` with standard output streamed to ``output``; preview the table."""
    executable = _check_executable(gemini_path)
    rendered = _p10.stdout_to_file(
        _PROGRAM,
        [*_global(annotation_dir), *arguments],
        output=output,
        cwd=output.parent,
        explicit=gemini_path,
        hint=_HINT,
        program_path=executable,
        environment=environment,
    )
    return _finish(
        operation,
        rendered,
        lambda: {"table": str(output), **_p10.table_preview(output, header=header, comment=None)},
    )


def _new_output(path: str | Path, name: str = "output_path") -> Path:
    return _cli.new_file(path, name)


def _target_database(
    database: str, output_database: str | None, in_place: bool
) -> tuple[str, Path | None]:
    """The database to change: the original (``in_place``) or a copy made at ``output_database``."""
    if not isinstance(in_place, bool):
        raise ConfigurationError("in_place must be True or False.")
    if in_place == (output_database is not None):
        raise ConfigurationError(
            "Give exactly one of output_database (a copy is changed) or in_place=True."
        )
    source = _database(database)
    if in_place:
        return source, None
    return source, _cli.new_file(output_database, "output_database")  # type: ignore[arg-type]


def _mutating(
    operation: str,
    subcommand: str,
    database: str,
    output_database: str | None,
    in_place: bool,
    options: Sequence[str],
    *,
    gemini_path: str | None,
    annotation_dir: str | None,
    positional_after: Sequence[str] = (),
) -> dict[str, Any]:
    """Run a sub-command that updates a database, on a copy unless ``in_place``."""
    source, copy = _target_database(database, output_database, in_place)
    executable = _check_executable(gemini_path)
    if copy is not None:
        shutil.copyfile(source, copy)
        target = copy
    else:
        target = Path(source)
    arguments = [*_global(annotation_dir), subcommand, str(target), *options, *positional_after]
    rendered = _cli.run(
        _PROGRAM, arguments, cwd=target.parent, expected=[target], explicit=executable, hint=_HINT
    )
    rendered["command"] = [_PROGRAM, *arguments]
    return _finish(
        operation,
        rendered,
        database=str(target),
        changed_copy=copy is not None,
    )


# ---------------------------------------------------------------------------
# Database construction and changes
# ---------------------------------------------------------------------------


def gemini_merge_chunks(
    chunk_databases: Sequence[str],
    output_database: str,
    *,
    vcf_path: str | None = None,
    anno_type: str | None = None,
    index: bool = False,
    skip_pls: bool = False,
    tempdir: str | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the partial databases written by a chunked parallel load into one database (``gemini merge_chunks``).

    ``chunk_databases`` are the chunk databases (``gemini load_chunk`` output of one VCF). The merged database is the new
    file ``output_database``. ``vcf_path`` is the original VCF (GEMINI reads extra annotation fields from it) and
    ``anno_type`` one of ``snpEff``, ``VEP``, ``BCFT`` or ``all``. ``index=True`` builds all database indexes (do it
    only in the last of several merges); ``tempdir`` is a local, non-NFS directory for SQLite working files.
    """
    chunks = [_database(path) for path in _list(chunk_databases, "chunk_databases") or ()]
    out = _new_output(output_database, "output_database")
    args = _Args("merge_chunks")
    args.items += ["--chunkdb", *chunks]
    args.value("--db", out)
    args.value("--vcf", None if vcf_path is None else str(_cli.existing_file(vcf_path, "vcf_path")))
    if anno_type is not None:
        args.value("-t", _cli.choice(anno_type, ("snpEff", "VEP", "BCFT", "all"), "anno_type"))
    args.switch("--index", index, "index")
    args.switch("--skip-pls", skip_pls, "skip_pls")
    if tempdir is not None:
        directory = Path(tempdir).expanduser().resolve()
        if not directory.is_dir():
            raise ConfigurationError(f"tempdir must be an existing directory: {directory}")
        args.value("--tempdir", directory)
    executable = _check_executable(gemini_path)
    arguments = [*_global(annotation_dir), *args.items]
    rendered = _cli.run(
        _PROGRAM, arguments, cwd=out.parent, expected=[out], explicit=executable, hint=_HINT
    )
    return _finish("gemini_merge_chunks", rendered, database=str(out))


def gemini_amend(
    database: str,
    sample_file: str | None = None,
    *,
    output_database: str | None = None,
    in_place: bool = False,
    clear: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Amend the sample or pedigree information of a GEMINI database from a PED file (``gemini amend --sample``).

    ``sample_file`` is a PED(+) file with new columns or values. ``clear=True`` first sets every value of the columns that
    the file provides to NULL. The change is made on a copy at ``output_database`` or, with ``in_place=True``, on
    ``database`` itself. Result key: ``database``.
    """
    options: list[str] = []
    if sample_file is None and not clear:
        raise ConfigurationError("Give sample_file (and optionally clear=True).")
    if sample_file is not None:
        options += ["--sample", str(_cli.existing_file(sample_file, "sample_file"))]
    if not isinstance(clear, bool):
        raise ConfigurationError("clear must be True or False.")
    if clear:
        options.append("--clear")
    return _mutating(
        "gemini_amend",
        "amend",
        database,
        output_database,
        in_place,
        options,
        gemini_path=gemini_path,
        annotation_dir=annotation_dir,
    )


def gemini_annotate(
    database: str,
    annotation_file: str,
    column_names: str,
    *,
    output_database: str | None = None,
    in_place: bool = False,
    annotation_type: str = "extract",
    extract_columns: str | None = None,
    column_types: str | None = None,
    operations: str | None = None,
    region_only: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate the variants of a GEMINI database with a tabix-indexed BED or VCF file (``gemini annotate``).

    ``annotation_file`` is a bgzipped, tabix-indexed BED (or VCF) and ``column_names`` the name(s) of the new variant
    table column(s) (for a VCF: the INFO field(s) to pull). ``annotation_type`` is ``boolean`` (overlap flag),
    ``count`` (number of overlapping records) or ``extract`` (GEMINI's default; values taken from
    ``extract_columns`` and combined by ``operations`` such as ``sum``, ``mean``, ``list``, ``first``). ``column_types``
    gives ``integer``, ``float`` or ``text`` per column. For a VCF annotation, ``region_only`` ignores REF/ALT.
    Works on a copy at ``output_database`` or, with ``in_place=True``, on ``database`` itself.
    """
    anno = str(_cli.existing_file(annotation_file, "annotation_file"))
    names = _columns(column_names, "column_names")
    options = [
        "-f",
        anno,
        "-c",
        names or "",
        "-a",
        _cli.choice(annotation_type, ("boolean", "count", "extract"), "annotation_type"),
    ]
    if extract_columns is not None:
        options += ["-e", _columns(extract_columns, "extract_columns") or ""]
    if column_types is not None:
        text = _cli.text(column_types, "column_types")
        if not all(item in ("integer", "float", "text") for item in text.split(",")):
            raise ConfigurationError(
                "column_types must be integer, float or text (comma-separated per column)."
            )
        options += ["-t", text]
    if operations is not None:
        text = _cli.text(operations, "operations")
        allowed = {
            "sum",
            "mean",
            "median",
            "min",
            "max",
            "mode",
            "list",
            "uniq_list",
            "first",
            "last",
        }
        if not all(item in allowed for item in text.split(",")):
            raise ConfigurationError(
                f"operations must be comma-separated values of {sorted(allowed)}."
            )
        options += ["-o", text]
    if not isinstance(region_only, bool):
        raise ConfigurationError("region_only must be True or False.")
    if region_only:
        options.append("--region-only")
    return _mutating(
        "gemini_annotate",
        "annotate",
        database,
        output_database,
        in_place,
        options,
        gemini_path=gemini_path,
        annotation_dir=annotation_dir,
    )


def gemini_set_somatic(
    database: str,
    *,
    output_database: str | None = None,
    in_place: bool = False,
    min_depth: float | None = None,
    min_quality: float | None = None,
    min_somatic_score: float | None = None,
    max_normal_alt_freq: float | None = None,
    max_normal_alt_count: int | None = None,
    min_normal_depth: int | None = None,
    min_tumor_alt_freq: float | None = None,
    min_tumor_alt_count: int | None = None,
    min_tumor_depth: int | None = None,
    chrom: str | None = None,
    dry_run: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag somatic variants (``is_somatic``) by comparing tumor and normal samples of a GEMINI database (``gemini set_somatic``).

    The thresholds are GEMINI's options of the same name (all default to 0 in GEMINI): combined tumor+normal depth, VCF
    QUAL, somatic score, maximum alternate frequency/count in the normal, minimum depth in the normal and the minimum
    alternate frequency/count/depth in the tumor. ``chrom`` restricts the tagging to one chromosome; ``dry_run=True``
    only reports what would be set (the database is not changed, but is still copied unless ``in_place``).
    The tumor/normal pairs come from the phenotype or PED information stored in the database.
    """
    options: list[str] = []
    for flag, value, name, kind in (
        ("--min-depth", min_depth, "min_depth", "number"),
        ("--min-qual", min_quality, "min_quality", "number"),
        ("--min-somatic-score", min_somatic_score, "min_somatic_score", "number"),
        ("--max-norm-alt-freq", max_normal_alt_freq, "max_normal_alt_freq", "number"),
        ("--max-norm-alt-count", max_normal_alt_count, "max_normal_alt_count", "int"),
        ("--min-norm-depth", min_normal_depth, "min_normal_depth", "int"),
        ("--min-tumor-alt-freq", min_tumor_alt_freq, "min_tumor_alt_freq", "number"),
        ("--min-tumor-alt-count", min_tumor_alt_count, "min_tumor_alt_count", "int"),
        ("--min-tumor-depth", min_tumor_depth, "min_tumor_depth", "int"),
    ):
        if value is not None:
            checked = (
                _cli.integer(value, name, minimum=0)
                if kind == "int"
                else _cli.number(value, name, minimum=0.0)
            )
            options += [flag, str(checked)]
    if chrom is not None:
        options += ["--chrom", _cli.text(chrom, "chrom")]
    if not isinstance(dry_run, bool):
        raise ConfigurationError("dry_run must be True or False.")
    if dry_run:
        options.append("--dry-run")
    return _mutating(
        "gemini_set_somatic",
        "set_somatic",
        database,
        output_database,
        in_place,
        options,
        gemini_path=gemini_path,
        annotation_dir=annotation_dir,
    )


# ---------------------------------------------------------------------------
# Query, export and region
# ---------------------------------------------------------------------------


def gemini_carrier_summary_by_phenotype(
    database: str,
    phenotype_column: str,
    output_path: str,
    *,
    query: str = "SELECT chrom, start, end, ref, alt, gene FROM variants",
    gt_filter: str | None = None,
    region: str | None = None,
    header: bool = True,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Counts of carriers and non-carriers of each variant, stratified by a sample phenotype column (``gemini query --carrier-summary-by-phenotype``).

    ``phenotype_column`` is a column of the samples table (for example ``phenotype``); the query's variants get extra columns with the
    number of carriers and non-carriers per phenotype value. ``query`` is the SQL that selects the variants (GEMINI
    requires one; the default lists the basic columns), ``gt_filter`` restricts genotypes and ``region`` the
    positions (``chr1:10-20``). The TSV goes to ``output_path``. Result keys: ``columns``, ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("query")
    args.switch("--header", header, "header")
    args.value("--carrier-summary-by-phenotype", _cli.text(phenotype_column, "phenotype_column"))
    args.value("--gt-filter", _optional_sql(gt_filter, "gt_filter"))
    args.value("--region", None if region is None else _cli.text(region, "region"))
    args.items += ["-q", _sql(query, "query"), db]
    return _table_command(
        "gemini_carrier_summary_by_phenotype",
        args.items,
        out,
        gemini_path=gemini_path,
        header=header,
        annotation_dir=annotation_dir,
    )


_DUMP_TABLES = ("variants", "genotypes", "samples", "tfam")


def gemini_dump(
    database: str,
    output_path: str,
    *,
    table: str = "variants",
    header: bool = False,
    separator: str = "\t",
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a table of a GEMINI database as delimited text (``gemini dump``).

    ``table`` selects the export: ``"variants"`` (``--variants``, every column of the variants table), ``"genotypes"``
    (``--genotypes``, one line per variant and sample), ``"samples"`` (``--samples``, the samples/pedigree table) or
    ``"tfam"`` (``--tfam``, the sample information in TFAM format). ``header=True`` adds the column names;
    ``separator`` is the output column separator (GEMINI's default is a tab). Result keys: ``columns``, ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    kind = _cli.choice(table, _DUMP_TABLES, "table")
    if (
        not isinstance(separator, str)
        or not separator
        or "\x00" in separator
        or separator.startswith("-")
    ):
        raise ConfigurationError(
            "separator must be a non-empty string that does not start with '-'."
        )
    args = _Args("dump", f"--{kind}")
    args.switch("--header", header, "header")
    args.value("--sep", separator)
    args.items.append(db)
    result = _table_command(
        "gemini_dump",
        args.items,
        out,
        gemini_path=gemini_path,
        header=header and kind != "tfam",
        annotation_dir=annotation_dir,
    )
    result["table_kind"] = kind
    return result


def gemini_region(
    database: str,
    output_path: str,
    *,
    gene: str | None = None,
    region: str | None = None,
    columns: str | None = None,
    filter: str | None = None,
    header: bool = False,
    show_samples: bool = False,
    output_format: str | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the variants in a gene (by symbol) or a genomic region (``gemini region --gene`` / ``--reg``).

    Give exactly one of ``gene`` (a gene symbol such as ``TP53``) or ``region`` (``chr1:100-200``). ``columns`` is the
    comma-separated list of columns to return (GEMINI's default is ``*``), ``filter`` an SQL restriction on the
    variants, ``show_samples`` adds the names of the samples that carry each variant and ``output_format``
    is ``JSON``, ``TPED`` or ``default``. Result keys: ``columns``, ``row_count``, ``rows``.
    """
    if (gene is None) == (region is None):
        raise ConfigurationError("Give exactly one of gene or region.")
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("region")
    args.value("--gene", None if gene is None else _cli.text(gene, "gene"))
    args.value("--reg", None if region is None else _cli.text(region, "region"))
    args.switch("--header", header, "header")
    args.value("--columns", _columns(columns))
    args.value("--filter", _optional_sql(filter, "filter"))
    args.switch("--show-samples", show_samples, "show_samples")
    if output_format is not None:
        args.value(
            "--format", _cli.choice(output_format, ("default", "JSON", "TPED"), "output_format")
        )
    args.items.append(db)
    return _table_command(
        "gemini_region",
        args.items,
        out,
        gemini_path=gemini_path,
        header=header,
        annotation_dir=annotation_dir,
    )


# ---------------------------------------------------------------------------
# Statistics and windows
# ---------------------------------------------------------------------------

_STATISTICS = {
    "tstv": "--tstv",
    "tstv_coding": "--tstv-coding",
    "tstv_noncoding": "--tstv-noncoding",
    "snp_counts": "--snp-counts",
    "sfs": "--sfs",
    "mds": "--mds",
    "vars_by_sample": "--vars-by-sample",
    "gts_by_sample": "--gts-by-sample",
}


def gemini_stats(
    database: str,
    output_path: str,
    *,
    statistic: str = "snp_counts",
    query: str | None = None,
    gt_filter: str | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Basic statistics of a GEMINI database (``gemini stats``).

    ``statistic`` selects the report: ``"snp_counts"`` (``--snp-counts``, count of each substitution type A->G, G->T, ...),
    ``"sfs"`` (``--sfs``, site frequency spectrum), ``"mds"`` (``--mds``, pairwise genetic distance between samples,
    i.e. multidimensional-scaling input), ``"vars_by_sample"`` (``--vars-by-sample``, variants per sample),
    ``"gts_by_sample"`` (``--gts-by-sample``, genotype-class counts per sample), ``"tstv"``, ``"tstv_coding"`` and
    ``"tstv_noncoding"`` (transition/transversion ratios) and ``"summarize"`` (``--summarize``, summarise the variants of
    your SQL ``query`` by sample). ``gt_filter`` restricts genotypes. Result keys: ``columns`` (None), ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    kind = _cli.choice(statistic, (*_STATISTICS, "summarize"), "statistic")
    args = _Args("stats")
    if kind == "summarize":
        if query is None:
            raise ConfigurationError("statistic='summarize' needs a query.")
        args.value("--summarize", _sql(query, "query"))
    else:
        if query is not None:
            raise ConfigurationError("query only applies to statistic='summarize'.")
        args.items.append(_STATISTICS[kind])
    args.value("--gt-filter", _optional_sql(gt_filter, "gt_filter"))
    args.items.append(db)
    result = _table_command(
        "gemini_stats",
        args.items,
        out,
        gemini_path=gemini_path,
        header=False,
        annotation_dir=annotation_dir,
    )
    result["statistic"] = kind
    return result


def gemini_windower(
    database: str,
    output_path: str,
    *,
    analysis: str = "hwe",
    window_size: int = 1_000_000,
    step_size: int = 0,
    operation: str = "mean",
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Windowed population-genetic statistic along the genome: Hardy-Weinberg or nucleotide diversity (``gemini windower``).

    ``analysis`` is ``"hwe"`` (``-t hwe``, windowed Hardy-Weinberg equilibrium test; GEMINI's default) or ``"nucl_div"``
    (``-t nucl_div``, windowed nucleotide diversity). ``window_size`` is the window length in bp (default 1,000,000),
    ``step_size`` the step in bp (0 lets GEMINI use the window size) and ``operation`` how the per-variant values of a
    window are combined: ``mean`` (default), ``median``, ``min``, ``max`` or ``collapse``. Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("windower")
    args.value("-t", _cli.choice(analysis, ("hwe", "nucl_div"), "analysis"))
    args.value("-w", _cli.integer(window_size, "window_size", minimum=1))
    args.value("-s", _cli.integer(step_size, "step_size", minimum=0))
    args.value(
        "-o", _cli.choice(operation, ("mean", "median", "min", "max", "collapse"), "operation")
    )
    args.items.append(db)
    return _table_command(
        "gemini_windower",
        args.items,
        out,
        gemini_path=gemini_path,
        header=False,
        annotation_dir=annotation_dir,
    )


# ---------------------------------------------------------------------------
# Inheritance models and gene-level tools
# ---------------------------------------------------------------------------

_MODELS = (
    "autosomal_recessive",
    "autosomal_dominant",
    "x_linked_recessive",
    "x_linked_dominant",
    "x_linked_de_novo",
)


def gemini_inheritance_model(
    database: str,
    output_path: str,
    *,
    model: str = "autosomal_recessive",
    columns: str | None = None,
    filter: str | None = None,
    min_kindreds: int | None = None,
    families: str | None = None,
    allow_unaffected: bool = False,
    min_sample_depth: int | None = None,
    min_gq: int | None = None,
    lenient: bool = False,
    max_phred_likelihood: int | None = None,
    x_chromosomes: Sequence[str] | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find candidate variants that fit a Mendelian inheritance model in the families of a GEMINI database (``gemini autosomal_recessive`` and friends).

    ``model`` is ``"autosomal_recessive"``, ``"autosomal_dominant"``, ``"x_linked_recessive"``, ``"x_linked_dominant"`` or
    ``"x_linked_de_novo"`` (the GEMINI sub-command of the same name). The pedigree comes from the database. ``columns`` selects the
    returned columns (default ``*``), ``filter`` restricts the variants (SQL), ``min_kindreds`` is the minimum
    number of families with a candidate in a gene (default 1), ``families`` a comma-separated list of family
    ids, ``allow_unaffected`` reports candidates that also affect unaffected samples, ``min_sample_depth``
    (``-d``) and ``min_gq`` set per-sample depth and genotype-quality minima (default 0). Autosomal models also take
    ``lenient`` (looser family structure) and ``max_phred_likelihood`` (``--gt-pl-max``); X-linked models take
    ``x_chromosomes`` (``-X``, names of the X chromosome when not ``chrX``/``X``). Options of the other kind raise
    ConfigurationError. Result keys: ``model``, ``row_count``, ``rows``.
    """
    kind = _cli.choice(model, _MODELS, "model")
    db = _database(database)
    out = _new_output(output_path)
    autosomal = kind.startswith("autosomal")
    if not autosomal and (lenient or max_phred_likelihood is not None):
        raise ConfigurationError(
            "lenient and max_phred_likelihood only apply to the autosomal models."
        )
    if autosomal and x_chromosomes:
        raise ConfigurationError("x_chromosomes only applies to the X-linked models.")
    args = _Args(kind)
    args.value("--columns", _columns(columns))
    args.value("--filter", _optional_sql(filter, "filter"))
    args.value("--min-kindreds", _optional_int(min_kindreds, "min_kindreds", minimum=1))
    if families is not None:
        text = _cli.text(families, "families")
        if _FAMILIES.fullmatch(text) is None:
            raise ConfigurationError("families must be a comma-separated list of family ids.")
        args.value("--families", text)
    args.switch("--allow-unaffected", allow_unaffected, "allow_unaffected")
    args.value("-d", _optional_int(min_sample_depth, "min_sample_depth", minimum=0))
    args.value("--min-gq", _optional_int(min_gq, "min_gq", minimum=0))
    args.switch("--lenient", lenient, "lenient")
    args.value(
        "--gt-pl-max", _optional_int(max_phred_likelihood, "max_phred_likelihood", minimum=0)
    )
    for chromosome in _list(x_chromosomes, "x_chromosomes") or ():
        args.value("-X", chromosome)
    args.items.append(db)
    result = _table_command(
        "gemini_inheritance_model",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )
    result["model"] = kind
    return result


def gemini_pathways(
    database: str,
    output_path: str,
    *,
    ensembl_version: str = "68",
    lof: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map the genes and variants of a GEMINI database to KEGG pathways (``gemini pathways``).

    ``ensembl_version`` is the Ensembl gene version (GEMINI supports 66 to 71; default 68). ``lof=True``
    (``--lof``) reports the pathways of individuals, genes and sites with loss-of-function variants. Result keys:
    ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    version = _cli.text(ensembl_version, "ensembl_version")
    if not version.isdigit():
        raise ConfigurationError("ensembl_version must be a number such as '68'.")
    args = _Args("pathways", "-v", version)
    args.switch("--lof", lof, "lof")
    args.items.append(db)
    return _table_command(
        "gemini_pathways",
        args.items,
        out,
        gemini_path=gemini_path,
        header=False,
        annotation_dir=annotation_dir,
    )


def gemini_lof_sieve(
    database: str,
    output_path: str,
    *,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prioritise loss-of-function variants by the likelihood that they are real LoF (``gemini lof_sieve``).

    Reports each LoF variant with the transcript position, the fraction of the transcript affected and whether the
    transcript is canonical, which separate likely true LoF from variants near the 3' end. Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    return _table_command(
        "gemini_lof_sieve",
        ["lof_sieve", db],
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


def gemini_burden(
    database: str,
    output_path: str,
    *,
    nonsynonymous: bool = False,
    cases: Sequence[str] | None = None,
    controls: Sequence[str] | None = None,
    min_aaf: float | None = None,
    max_aaf: float | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene-level variant burden counts per sample (or case/control group) (``gemini burden``).

    By default every high-impact variant in coding regions counts towards a gene's burden; ``nonsynonymous=True``
    (``--nonsynonymous``) counts all nonsynonymous variants. ``cases`` and ``controls`` are sample names; with both,
    GEMINI adds the burden test of cases against controls. ``min_aaf``/``max_aaf`` bound the alternate allele
    frequency of the variants used (GEMINI's defaults 0.0 and 1.0). The C-alpha test is
    ``gemini_calpha_test``. Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("burden")
    args.switch("--nonsynonymous", nonsynonymous, "nonsynonymous")
    case_names = _list(cases, "cases")
    control_names = _list(controls, "controls")
    if case_names:
        args.items += ["--cases", *case_names]
    if control_names:
        args.items += ["--controls", *control_names]
    args.value("--min-aaf", _optional_number(min_aaf, "min_aaf", minimum=0.0, maximum=1.0))
    args.value("--max-aaf", _optional_number(max_aaf, "max_aaf", minimum=0.0, maximum=1.0))
    args.items.append(db)
    return _table_command(
        "gemini_burden",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


def gemini_calpha_test(
    database: str,
    output_path: str,
    cases: Sequence[str],
    controls: Sequence[str],
    *,
    permutations: int = 0,
    nonsynonymous: bool = False,
    min_aaf: float | None = None,
    max_aaf: float | None = None,
    save_tscores: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene-level C-alpha association test between cases and controls (``gemini burden --calpha``).

    ``cases`` and ``controls`` are sample names of the database (both required). ``permutations`` is the number of
    permutations for the permutation p-value (0: asymptotic only; try 1000 to start); ``save_tscores`` keeps the
    permuted T-scores in a file next to the working directory of GEMINI. The other options are those of
    ``gemini_burden``. Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    case_names = _list(cases, "cases")
    control_names = _list(controls, "controls")
    args = _Args("burden", "--calpha")
    args.switch("--nonsynonymous", nonsynonymous, "nonsynonymous")
    args.items += ["--cases", *(case_names or [])]
    args.items += ["--controls", *(control_names or [])]
    args.value("--permutations", _cli.integer(permutations, "permutations", minimum=0))
    args.value("--min-aaf", _optional_number(min_aaf, "min_aaf", minimum=0.0, maximum=1.0))
    args.value("--max-aaf", _optional_number(max_aaf, "max_aaf", minimum=0.0, maximum=1.0))
    args.switch("--save_tscores", save_tscores, "save_tscores")
    args.items.append(db)
    return _table_command(
        "gemini_calpha_test",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


def gemini_interactions(
    database: str,
    output_path: str,
    *,
    lof: bool = False,
    gene: str | None = None,
    radius: int | None = None,
    edges_file: str | None = None,
    variant_mode: bool = False,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene interaction partners of the genes that carry variants in the samples (``gemini interactions`` and ``gemini lof_interactions``).

    ``lof=False`` uses ``gemini interactions``: with ``gene`` (``-g``) as the root and ``radius`` (``-r``) the breadth-first
    search depth (0 and up). ``lof=True`` uses ``gemini lof_interactions`` for genes with loss-of-function variants
    (``gene`` does not apply). ``edges_file`` replaces the default HPRD network with ``geneA|geneB`` lines.
    ``variant_mode`` (``--var``) also returns variant information (impact, biotype) for the interacting genes. Result
    keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    if lof and gene is not None:
        raise ConfigurationError("gene does not apply to lof=True (gemini lof_interactions).")
    args = _Args("lof_interactions" if lof else "interactions")
    if gene is not None:
        args.value("-g", _cli.text(gene, "gene"))
    args.value("-r", _optional_int(radius, "radius", minimum=0))
    args.value(
        "--edges", None if edges_file is None else str(_cli.existing_file(edges_file, "edges_file"))
    )
    args.switch("--var", variant_mode, "variant_mode")
    args.items.append(db)
    return _table_command(
        "gemini_interactions",
        args.items,
        out,
        gemini_path=gemini_path,
        header=False,
        annotation_dir=annotation_dir,
    )


def gemini_fusions(
    database: str,
    output_path: str,
    *,
    in_cosmic_census: bool = False,
    min_quality: float | None = None,
    evidence_type: str | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect somatic gene fusions from structural-variant breakends in a GEMINI database (``gemini fusions``).

    ``in_cosmic_census=True`` keeps fusions where one or both genes are in the COSMIC cancer census, ``min_quality`` is the
    minimum VCF QUAL and ``evidence_type`` the supporting evidence (``PE``, ``SR`` or ``PE,SR``). Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("fusions")
    args.switch("--in_cosmic_census", in_cosmic_census, "in_cosmic_census")
    args.value("--min_qual", _optional_number(min_quality, "min_quality", minimum=0.0))
    if evidence_type is not None:
        args.value(
            "--evidence_type", _cli.choice(evidence_type, ("PE", "SR", "PE,SR"), "evidence_type")
        )
    args.items.append(db)
    return _table_command(
        "gemini_fusions",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


def gemini_gene_wise(
    database: str,
    output_path: str,
    *,
    genotype_filters: Sequence[str] | None = None,
    required_genotype_filters: Sequence[str] | None = None,
    min_filters: int | None = None,
    where: str | None = None,
    filter: str | None = None,
    columns: str | None = None,
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Aggregate variants per gene with genotype filters, reporting genes that pass enough filters (``gemini gene_wise``).

    ``genotype_filters`` (``--gt-filter``) are genotype expressions evaluated per variant; a gene is reported
    when at least ``min_filters`` of them (default 1) pass among its variants. ``required_genotype_filters``
    (``--gt-filter-required``) must pass for a variant to count at all. ``where`` is the SQL WHERE clause that selects
    the variants (GEMINI's default is ``is_exonic = 1 AND impact_severity != 'LOW'``), ``filter`` an extra SQL
    restriction and ``columns`` the columns to return (default ``chrom,start,end,gene,impact,impact_severity,max_aaf_all``).
    Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("gene_wise")
    args.value("--min-filters", _optional_int(min_filters, "min_filters", minimum=1))
    args.value("--where", _optional_sql(where, "where"))
    for expression in genotype_filters or ():
        args.value("--gt-filter", _sql(expression, "genotype_filters"))
    for expression in required_genotype_filters or ():
        args.value("--gt-filter-required", _sql(expression, "required_genotype_filters"))
    args.value("--filter", _optional_sql(filter, "filter"))
    args.value("--columns", _columns(columns))
    args.items.append(db)
    return _table_command(
        "gemini_gene_wise",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


def gemini_sex_check(
    database: str,
    output_path: str,
    *,
    chromosome: str = "chrX",
    annotation_dir: str | None = None,
    gemini_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check the recorded sex of each sample against its genotypes on the X chromosome (``gemini qc --mode sex``).

    ``chromosome`` is the name of the X chromosome in the database (GEMINI's default ``chrX``; use ``X`` for
    Ensembl-style names). The table has one row per sample with the sex in the pedigree, the number of homozygous/
    heterozygous sites and the inferred sex. Result keys: ``row_count``, ``rows``.
    """
    db = _database(database)
    out = _new_output(output_path)
    args = _Args("qc", "--mode", "sex", "--chrom", _cli.text(chromosome, "chromosome"), db)
    return _table_command(
        "gemini_sex_check",
        args.items,
        out,
        gemini_path=gemini_path,
        header=True,
        annotation_dir=annotation_dir,
    )


__all__ = [
    "gemini_amend",
    "gemini_annotate",
    "gemini_burden",
    "gemini_calpha_test",
    "gemini_carrier_summary_by_phenotype",
    "gemini_dump",
    "gemini_fusions",
    "gemini_gene_wise",
    "gemini_inheritance_model",
    "gemini_interactions",
    "gemini_lof_sieve",
    "gemini_merge_chunks",
    "gemini_pathways",
    "gemini_region",
    "gemini_set_somatic",
    "gemini_sex_check",
    "gemini_stats",
    "gemini_windower",
]

_register_public_verbose(globals(), tuple(__all__))
