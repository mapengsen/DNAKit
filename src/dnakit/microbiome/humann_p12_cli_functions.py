"""HUMAnN command-line facades: gene-family and pathway quantification and the HUMAnN table utilities.

Each public function runs one fixed HUMAnN program (https://github.com/biobakery/humann; 3.9 is the tested version):
``humann`` for the quantification of gene families and pathways from existing alignments or tables, and the
``humann_*`` table utilities. No shell is used and every argument comes from a validated parameter. The user installs
HUMAnN (``conda install -c bioconda humann`` or ``pip install humann``); the programs are looked up on PATH, then in the
folder of the running interpreter, or in ``humann_dir`` when that folder is given. A missing program raises
``BackendUnavailableError`` with an install hint. Results are JSON-serialisable dictionaries with the status, backend,
version, command, outputs, a small parsed view of the tables and a bounded log tail; outputs are written to new paths
and never overwrite existing files. Profiling reads from FASTQ is available as ``dnakit.microbiome.profile_humann``.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p12_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="humann",
    executable="humann",
    hint=(
        "Install HUMAnN 3.x (for example `conda install -c bioconda humann` or `pip install humann`) and put its "
        "programs on PATH or pass humann_dir. Some utilities need the HUMAnN utility mapping data "
        "(`humann_databases --download utility_mapping full <folder>`)."
    ),
    tested_version="3.9",
    citation="https://doi.org/10.7554/eLife.65088",
    version_arguments=("--version",),
    version_pattern=r"humann\S*\s+v?(\d[\w.]*)",
)
_INPUT_FORMATS = ("sam", "bam", "blastm8", "genetable", "biom")
_PREVIEW_ROWS = 20
_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_LEVELS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
_RENAME_CHOICES = (
    "kegg-orthology",
    "kegg-pathway",
    "kegg-module",
    "ec",
    "metacyc-rxn",
    "metacyc-pwy",
    "pfam",
    "eggnog",
    "go",
    "infogo1000",
    "uniref50",
)


def _explicit(humann_dir: str | None, name: str) -> str | None:
    if humann_dir is None:
        return None
    folder = Path(humann_dir).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError(f"humann_dir must be an existing folder: {folder}")
    return str(folder / name)


def _version(humann_dir: str | None) -> str | None:
    try:
        executable = _cli.program(
            "humann", explicit=_explicit(humann_dir, "humann"), hint=_TOOL.hint
        )
    except Exception:
        return None
    return _cli._version_of(executable, ("--version",), r"humann\S*\s+v?(\d[\w.]*)")


def _run(
    name: str,
    arguments: list[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    humann_dir: str | None,
    operation: str,
    timeout_seconds: float = 86_400.0,
    **extra: Any,
) -> dict[str, Any]:
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=cwd,
        expected=list(expected),
        explicit=_explicit(humann_dir, name),
        executable_name=name,
        timeout_seconds=timeout_seconds,
    )
    if record.get("version") is None:
        record["version"] = _version(humann_dir)
    return _cli.finish(record, operation, **extra)


def _table(path: Path, preview: int = _PREVIEW_ROWS) -> dict[str, Any]:
    """A small view of a HUMAnN TSV table: its columns, the number of feature rows and the first rows."""
    columns: list[str] = []
    rows: list[list[Any]] = []
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not columns:
                columns = [item.lstrip("# ").strip() for item in fields]
                continue
            if not line.strip() or line.startswith("#"):
                continue
            count += 1
            if len(rows) < preview:
                rows.append([fields[0], *(_cli.convert_number(item) for item in fields[1:])])
    return {"file": str(path), "columns": columns, "row_count": count, "rows": rows}


def _table_or_none(path: Path) -> dict[str, Any] | None:
    return _table(path) if path.suffix == ".tsv" else {"file": str(path)}


def _name(value: str, label: str) -> str:
    if not isinstance(value, str) or not _SAFE_NAME.match(value):
        raise ConfigurationError(
            f"{label} must start with a letter or digit and contain only letters, digits, '_', '.', '-'."
        )
    return value


def humann_quantify_alignments(
    input_path: str,
    output_dir: str,
    *,
    input_format: str = "blastm8",
    minpath: bool = True,
    gap_fill: bool = True,
    pathways: str = "metacyc",
    pathways_database: Sequence[str] | None = None,
    id_mapping: str | None = None,
    annotation_gene_index: int | None = None,
    output_basename: str | None = None,
    output_format: str = "tsv",
    output_max_decimals: int | None = None,
    remove_stratified_output: bool = False,
    remove_column_description_output: bool = False,
    remove_temp_output: bool = False,
    nucleotide_database: str | None = None,
    protein_database: str | None = None,
    threads: int = 1,
    timeout_seconds: float = 86_400.0,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Quantify gene families and pathways from alignments or from a gene family table with ``humann``.

    This function covers three operations of the work list through ``input_format``: ``sam`` or ``bam`` (read alignments to
    a UniRef/ChocoPhlAn reference, ``--input-format sam|bam``), ``blastm8`` (a BLAST-like tabular alignment file such as
    DIAMOND output, ``--input-format blastm8``) and ``genetable`` or ``biom`` (an existing gene family table, from which
    only pathway abundance and coverage are computed). The XIPE option (``--xipe on``) is not offered: in HUMAnN 3.9 it
    fails on Python 3 (``xipe.py`` uses the Python 2 name ``file``). The search steps of HUMAnN are not run. ``minpath`` and ``gap_fill`` switch the MinPath and
    gap-fill steps (both on by default as in HUMAnN; MinPath must be installed for ``minpath=True``); ``pathways`` is
    ``metacyc`` or ``unipathway``; ``pathways_database`` is one or two custom pathway mapping files (``--pathways-database``);
    ``id_mapping`` maps alignment reference ids to gene families (``--id-mapping``) and ``annotation_gene_index`` is the
    column of the gene id in the reference annotation (``--annotation-gene-index``). ``output_basename`` names the outputs
    (default: the input file name before the first dot); ``output_format`` is ``tsv`` or ``biom``;
    ``remove_stratified_output`` and ``remove_column_description_output`` simplify the tables.
    ``output_dir`` must not exist yet. Result keys: ``output_dir``, ``gene_families`` (absent for ``genetable``/``biom``
    input), ``pathway_abundance`` and ``pathway_coverage`` (each with ``file``, ``columns``, ``row_count`` and the first
    ``rows``; for ``biom`` output only ``file``).
    """
    source = _cli.existing_file(input_path, "input_path")
    fmt = _cli.choice(input_format, _INPUT_FORMATS, "input_format")
    directory = _cli.new_directory(output_dir, "output_dir")
    base = _name(
        output_basename if output_basename is not None else source.name.split(".", 1)[0],
        "output_basename",
    )
    out_fmt = _cli.choice(output_format, ("tsv", "biom"), "output_format")
    arguments = [
        "--input",
        str(source),
        "--output",
        str(directory),
        "--input-format",
        fmt,
        "--output-basename",
        base,
        "--output-format",
        out_fmt,
        "--threads",
        str(_cli.threads_value(threads)),
        "--minpath",
        "on" if _cli.flag(minpath, "minpath") else "off",
        "--gap-fill",
        "on" if _cli.flag(gap_fill, "gap_fill") else "off",
        "--pathways",
        _cli.choice(pathways, ("metacyc", "unipathway"), "pathways"),
    ]
    if pathways_database is not None:
        files = _cli.existing_files(list(pathways_database), "pathways_database")
        if len(files) > 2:
            raise ConfigurationError("pathways_database takes at most two files.")
        arguments += ["--pathways-database", ",".join(str(path) for path in files)]
    if id_mapping is not None:
        arguments += ["--id-mapping", str(_cli.existing_file(id_mapping, "id_mapping"))]
    if annotation_gene_index is not None:
        arguments += [
            "--annotation-gene-index",
            str(_cli.integer(annotation_gene_index, "annotation_gene_index", minimum=0)),
        ]
    if output_max_decimals is not None:
        arguments += [
            "--output-max-decimals",
            str(_cli.integer(output_max_decimals, "output_max_decimals", minimum=0, maximum=20)),
        ]
    if nucleotide_database is not None:
        arguments += [
            "--nucleotide-database",
            str(_cli.existing_directory(nucleotide_database, "nucleotide_database")),
        ]
    if protein_database is not None:
        arguments += [
            "--protein-database",
            str(_cli.existing_directory(protein_database, "protein_database")),
        ]
    for value, flag_name, option in (
        (remove_stratified_output, "remove_stratified_output", "--remove-stratified-output"),
        (
            remove_column_description_output,
            "remove_column_description_output",
            "--remove-column-description-output",
        ),
        (remove_temp_output, "remove_temp_output", "--remove-temp-output"),
    ):
        if _cli.flag(value, flag_name):
            arguments.append(option)
    suffix = ".tsv" if out_fmt == "tsv" else ".biom"
    expected = {
        "pathway_abundance": directory / f"{base}_pathabundance{suffix}",
        "pathway_coverage": directory / f"{base}_pathcoverage{suffix}",
    }
    if fmt not in ("genetable", "biom"):
        expected["gene_families"] = directory / f"{base}_genefamilies{suffix}"
    result = _run(
        "humann",
        arguments,
        cwd=directory,
        expected=list(expected.values()),
        humann_dir=humann_dir,
        operation="humann_quantify_alignments",
        timeout_seconds=timeout_seconds,
        input_format=fmt,
        output_dir=str(directory),
    )
    for key, path in expected.items():
        result[key] = _table_or_none(path)
    return result


def humann_join_tables(
    input_dir: str,
    output_path: str,
    *,
    file_name: str | None = None,
    search_subdirectories: bool = False,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join the per-sample gene, pathway or taxonomy tables of a folder into one multi-sample table (``humann_join_tables``).

    Every table in ``input_dir`` (or only those whose file name contains ``file_name``, ``--file_name``) becomes columns of
    ``output_path``; ``search_subdirectories`` also reads sub-folders (``-s``). Result keys: ``table`` (``file``,
    ``columns``, ``row_count`` and the first ``rows``).
    """
    folder = _cli.existing_directory(input_dir, "input_dir")
    target = _cli.new_file(output_path, "output_path")
    arguments = ["--input", str(folder), "--output", str(target)]
    if file_name is not None:
        arguments += ["--file_name", _cli.text(file_name, "file_name")]
    if _cli.flag(search_subdirectories, "search_subdirectories"):
        arguments.append("--search-subdirectories")
    result = _run(
        "humann_join_tables",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_join_tables",
    )
    result["table"] = _table(target)
    return result


def humann_split_table(
    input_path: str,
    output_dir: str,
    *,
    taxonomy_index: int | None = None,
    taxonomy_level: str | None = None,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a merged gene or pathway table into one table per sample (``humann_split_table``).

    ``output_dir`` must not exist; one table per sample column is written there. ``taxonomy_index`` is the index of the
    gene in the taxonomy data and ``taxonomy_level`` (Kingdom, Phylum, Class, Order, Family, Genus or Species) the level of
    taxonomy for the output when the input comes from PICRUSt ``metagenome_contributions.py``. Result keys: ``tables`` (the
    files written, sorted) and ``table_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    directory = _cli.new_directory(output_dir, "output_dir")
    arguments = ["--input", str(source), "--output", str(directory)]
    if taxonomy_index is not None:
        arguments += [
            "--taxonomy_index",
            str(_cli.integer(taxonomy_index, "taxonomy_index", minimum=0)),
        ]
    if taxonomy_level is not None:
        arguments += [
            "--taxonomy_level",
            _cli.choice(taxonomy_level, (*_LEVELS, "Species"), "taxonomy_level"),
        ]
    result = _run(
        "humann_split_table",
        arguments,
        cwd=directory,
        expected=[],
        humann_dir=humann_dir,
        operation="humann_split_table",
    )
    tables = sorted(str(path) for path in directory.iterdir() if path.is_file())
    if not tables:
        raise ConfigurationError("humann_split_table wrote no table.")
    result.update(tables=tables, table_count=len(tables), outputs=tables)
    return result


def humann_split_stratified_table(
    input_path: str,
    output_dir: str,
    *,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a stratified table into a stratified table and an unstratified (community totals) table (``humann_split_stratified_table``).

    ``input_path`` is a TSV (plain, gzip or bzip2) or BIOM HUMAnN table; ``output_dir`` must not exist and receives
    ``<name>_stratified.tsv`` and ``<name>_unstratified.tsv``. Result keys: ``stratified`` and ``unstratified`` (each with
    ``file``, ``columns``, ``row_count`` and the first ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    directory = _cli.new_directory(output_dir, "output_dir")
    result = _run(
        "humann_split_stratified_table",
        ["--input", str(source), "--output", str(directory)],
        cwd=directory,
        expected=[],
        humann_dir=humann_dir,
        operation="humann_split_stratified_table",
    )
    stratified = sorted(directory.glob("*_stratified.tsv"))
    unstratified = sorted(directory.glob("*_unstratified.tsv"))
    stratified = [path for path in stratified if path not in unstratified]
    if len(stratified) != 1 or len(unstratified) != 1:
        raise ConfigurationError(
            "humann_split_stratified_table did not write the expected two tables."
        )
    result.update(
        stratified=_table(stratified[0]),
        unstratified=_table(unstratified[0]),
        outputs=[str(stratified[0]), str(unstratified[0])],
    )
    return result


def humann_rename_table(
    input_path: str,
    output_path: str,
    *,
    names: str | None = None,
    custom_mapping: str | None = None,
    simplify: bool = False,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attach feature names (glosses) to the gene or pathway ids of a table (``humann_rename_table``).

    Give exactly one of ``names`` (a built-in mapping: kegg-orthology, kegg-pathway, kegg-module, ec, metacyc-rxn,
    metacyc-pwy, pfam, eggnog, go, infogo1000 or uniref50; these need the HUMAnN utility mapping data) and
    ``custom_mapping`` (a two-column TSV, optionally gzipped, from feature id to full name). ``simplify`` removes
    non-alphanumeric characters from the names. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first
    ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    if (names is None) == (custom_mapping is None):
        raise ConfigurationError("Give exactly one of names and custom_mapping.")
    arguments = ["--input", str(source)]
    if names is not None:
        arguments += ["--names", _cli.choice(names, _RENAME_CHOICES, "names")]
    else:
        arguments += ["--custom", str(_cli.existing_file(custom_mapping, "custom_mapping"))]  # type: ignore[arg-type]
    if _cli.flag(simplify, "simplify"):
        arguments.append("--simplify")
    arguments += ["--output", str(target)]
    result = _run(
        "humann_rename_table",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_rename_table",
    )
    result["table"] = _table(target)
    return result


def humann_renorm_table(
    input_path: str,
    output_path: str,
    *,
    units: str = "cpm",
    mode: str = "community",
    include_special: bool = True,
    update_sample_names: bool = False,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Renormalise gene or pathway abundances to copies per million or relative abundance (``humann_renorm_table``).

    ``units`` is ``cpm`` (copies per million, the HUMAnN default) or ``relab`` (relative abundance); ``mode`` normalises all
    levels of a stratified table by the ``community`` total (default) or by their own totals (``levelwise``);
    ``include_special`` keeps the special features UNMAPPED, UNINTEGRATED and UNGROUPED (``--special y``, the default);
    ``update_sample_names`` replaces ``-RPK`` in the sample names by the new suffix (``--update-snames``). Result keys:
    ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        "--input",
        str(source),
        "--units",
        _cli.choice(units, ("cpm", "relab"), "units"),
        "--mode",
        _cli.choice(mode, ("community", "levelwise"), "mode"),
        "--special",
        "y" if _cli.flag(include_special, "include_special") else "n",
    ]
    if _cli.flag(update_sample_names, "update_sample_names"):
        arguments.append("--update-snames")
    arguments += ["--output", str(target)]
    result = _run(
        "humann_renorm_table",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_renorm_table",
    )
    result["table"] = _table(target)
    return result


def humann_regroup_table(
    input_path: str,
    output_path: str,
    *,
    groups: str | None = None,
    custom_groups: str | None = None,
    reversed_mapping: bool = False,
    function: str = "sum",
    precision: int | None = None,
    include_ungrouped: bool = True,
    keep_protected: bool = True,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Regroup gene families into functional groups such as GO, KO, EC, Pfam or MetaCyc reactions (``humann_regroup_table``).

    Give exactly one of ``groups`` (a built-in mapping, ``uniref90_rxn`` or ``uniref50_rxn``, which needs the HUMAnN
    utility mapping data; other mappings such as GO, KO, EC and Pfam are supplied with ``custom_groups``) and
    ``custom_groups`` (a TSV or gzipped TSV with one group per line followed by its member features; with
    ``reversed_mapping=True`` one feature per line followed by its groups). ``function`` combines grouped features
    (``sum`` or ``mean``), ``precision`` rounds to that many decimals, ``include_ungrouped`` adds an UNGROUPED group
    (``--ungrouped Y``) and ``keep_protected`` carries UNMAPPED through (``--protected Y``). Result keys: ``table`` (``file``,
    ``columns``, ``row_count`` and the first ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    if (groups is None) == (custom_groups is None):
        raise ConfigurationError("Give exactly one of groups and custom_groups.")
    arguments = ["--input", str(source)]
    if groups is not None:
        arguments += ["--groups", _cli.choice(groups, ("uniref90_rxn", "uniref50_rxn"), "groups")]
        if reversed_mapping:
            raise ConfigurationError("reversed_mapping only applies to custom_groups.")
    else:
        arguments += ["--custom", str(_cli.existing_file(custom_groups, "custom_groups"))]  # type: ignore[arg-type]
        if _cli.flag(reversed_mapping, "reversed_mapping"):
            arguments.append("--reversed")
    arguments += ["--function", _cli.choice(function, ("sum", "mean"), "function")]
    if precision is not None:
        arguments += [
            "--precision",
            str(_cli.integer(precision, "precision", minimum=0, maximum=15)),
        ]
    arguments += [
        "--ungrouped",
        "Y" if _cli.flag(include_ungrouped, "include_ungrouped") else "N",
        "--protected",
        "Y" if _cli.flag(keep_protected, "keep_protected") else "N",
        "--output",
        str(target),
    ]
    result = _run(
        "humann_regroup_table",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_regroup_table",
    )
    result["table"] = _table(target)
    return result


def humann_reduce_table(
    input_path: str,
    output_path: str,
    *,
    function: str = "max",
    sort_by: str | None = None,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reduce a gene or pathway table to one value per feature with a summary function (``humann_reduce_table``).

    ``function`` is ``max`` (the HUMAnN default), ``sum``, ``min`` or ``mean`` and is applied across the sample columns
    of every feature; ``sort_by`` (``level``, ``name`` or ``value``) sorts the output. Result keys: ``table`` (``file``,
    ``columns``, ``row_count`` and the first ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        "--input",
        str(source),
        "--output",
        str(target),
        "--function",
        _cli.choice(function, ("sum", "min", "max", "mean"), "function"),
    ]
    if sort_by is not None:
        arguments += ["--sort-by", _cli.choice(sort_by, ("level", "name", "value"), "sort_by")]
    result = _run(
        "humann_reduce_table",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_reduce_table",
    )
    result["table"] = _table(target)
    return result


def humann_unpack_pathways(
    gene_table: str,
    pathway_table: str,
    output_path: str,
    *,
    gene_mapping: str | None = None,
    pathway_mapping: str | None = None,
    remove_taxonomy: bool = False,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unpack pathway abundances into the abundances of the genes that make up each pathway (``humann_unpack_pathways``).

    ``gene_table`` is the gene family (or EC) abundance table and ``pathway_table`` the pathway abundance table of the same
    samples. ``gene_mapping`` (gene family to reaction) and ``pathway_mapping`` (reaction to pathway) default to the MetaCyc
    files shipped with HUMAnN. ``remove_taxonomy`` drops the taxonomy strata (``--remove-taxonomy``). Result keys: ``table``
    (``file``, ``columns``, ``row_count`` and the first ``rows``).
    """
    genes = _cli.existing_file(gene_table, "gene_table")
    pathways = _cli.existing_file(pathway_table, "pathway_table")
    target = _cli.new_file(output_path, "output_path")
    arguments = ["--input-genes", str(genes), "--input-pathways", str(pathways)]
    if gene_mapping is not None:
        arguments += ["--gene-mapping", str(_cli.existing_file(gene_mapping, "gene_mapping"))]
    if pathway_mapping is not None:
        arguments += [
            "--pathway-mapping",
            str(_cli.existing_file(pathway_mapping, "pathway_mapping")),
        ]
    if _cli.flag(remove_taxonomy, "remove_taxonomy"):
        arguments.append("--remove-taxonomy")
    arguments += ["--output", str(target)]
    result = _run(
        "humann_unpack_pathways",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_unpack_pathways",
    )
    result["table"] = _table(target)
    return result


def humann_infer_taxonomy(
    input_path: str,
    output_path: str,
    *,
    level: str = "Family",
    database: str | None = None,
    database_file: str | None = None,
    mode: str = "totals",
    lca_choice: str = "humann_lca",
    threshold: float = 0.001,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer the taxonomy of unclassified gene families from the lowest common ancestor of their UniRef cluster (``humann_infer_taxonomy``).

    ``input_path`` is a HUMAnN gene family table. The LCA data come from ``database`` (the name of a tree-of-life/LCA file
    of the HUMAnN utility mapping data, for example ``uniref90`` or ``uniref50``) or from ``database_file`` (a path,
    ``--devdb``); give exactly one. ``level`` is the target taxonomic level (Kingdom, Phylum, Class, Order, Family or
    Genus), ``mode`` selects the rows used (``totals``, ``unclassified`` or ``stratified``), ``lca_choice`` the per-gene
    annotation (``source_tax``, ``uniref_lca`` or ``humann_lca``) and ``threshold`` the minimum frequency of a new taxon.
    Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``) and ``summary`` (the
    stratification summary line of the program, when present).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    if (database is None) == (database_file is None):
        raise ConfigurationError("Give exactly one of database and database_file.")
    arguments = [
        "--input",
        str(source),
        "--output",
        str(target),
        "--level",
        _cli.choice(level, _LEVELS, "level"),
        "--mode",
        _cli.choice(mode, ("totals", "unclassified", "stratified"), "mode"),
        "--lca-choice",
        _cli.choice(lca_choice, ("source_tax", "uniref_lca", "humann_lca"), "lca_choice"),
        "--threshold",
        format(_cli.number(threshold, "threshold", minimum=0.0, maximum=1.0), ".12g"),
    ]
    if database is not None:
        arguments += ["--database", _name(database, "database")]
    else:
        arguments += ["--devdb", str(_cli.existing_file(database_file, "database_file"))]  # type: ignore[arg-type]
    result = _run(
        "humann_infer_taxonomy",
        arguments,
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_infer_taxonomy",
    )
    match = re.search(r"Summary: .*", result["log_tail"])
    result.update(table=_table(target), summary=match.group(0) if match else None)
    return result


def humann_expand_cluster(
    input_path: str,
    uniref90: str,
    output_path: str,
    *,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand one UniRef90 gene family to all UniRef90 families in its UniRef50 cluster (``humann_expand_cluster``).

    ``input_path`` is a UniRef90 gene family table (TSV) and ``uniref90`` the family of interest (with or without the
    ``UniRef90_`` prefix). The output table keeps the input rows of every UniRef90 that clusters with it. The program
    reads the UniRef90-to-UniRef50 mapping of the HUMAnN utility mapping data. Result keys: ``table`` (``file``, ``columns``,
    ``row_count`` and the first ``rows``).
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    gene = _cli.text(uniref90, "uniref90")
    if not re.match(r"^(UniRef90_)?[A-Za-z0-9_]+$", gene):
        raise ConfigurationError("uniref90 must be a UniRef90 identifier.")
    result = _run(
        "humann_expand_cluster",
        ["--input", str(source), "--gene", gene, "--output", str(target)],
        cwd=target.parent,
        expected=[target],
        humann_dir=humann_dir,
        operation="humann_expand_cluster",
    )
    result["table"] = _table(target)
    return result


def humann_build_custom_database(
    fasta_path: str,
    output_dir: str,
    *,
    id_mapping: str | None = None,
    taxonomic_profile: str | None = None,
    database_format: str = "fasta",
    genus_abundance_threshold: float | None = None,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a custom nucleotide or protein search database from gene sequences, optionally limited to genera (``humann_build_custom_database``).

    ``fasta_path`` is a FASTA file of gene sequences. ``id_mapping`` maps FASTA ids to genus names and
    ``taxonomic_profile`` is a MetaPhlAn profile; with both, only the genera whose abundance reaches
    ``genus_abundance_threshold`` (``--genus-abundance-threshold``) are kept. ``database_format`` is ``fasta`` (a filtered
    FASTA file, the default) or ``diamond`` (also builds a DIAMOND index, which needs ``diamond`` on PATH). ``output_dir``
    must not exist. Result keys: ``output_dir`` and ``database_files`` (the files written, sorted).
    """
    source = _cli.existing_file(fasta_path, "fasta_path")
    directory = _cli.new_directory(output_dir, "output_dir")
    arguments = ["--input", str(source), "--output", str(directory)]
    if id_mapping is not None:
        arguments += ["--id-mapping", str(_cli.existing_file(id_mapping, "id_mapping"))]
    if taxonomic_profile is not None:
        arguments += [
            "--taxonomic-profile",
            str(_cli.existing_file(taxonomic_profile, "taxonomic_profile")),
        ]
    arguments += ["--format", _cli.choice(database_format, ("fasta", "diamond"), "database_format")]
    if genus_abundance_threshold is not None:
        arguments += [
            "--genus-abundance-threshold",
            format(
                _cli.number(genus_abundance_threshold, "genus_abundance_threshold", minimum=0.0),
                ".12g",
            ),
        ]
    result = _run(
        "humann_build_custom_database",
        arguments,
        cwd=directory,
        expected=[],
        humann_dir=humann_dir,
        operation="humann_build_custom_database",
    )
    files = sorted(str(path) for path in directory.iterdir())
    if not files:
        raise ConfigurationError("humann_build_custom_database wrote no file.")
    result.update(output_dir=str(directory), database_files=files, outputs=files)
    return result


def humann_barplot(
    input_path: str,
    focal_feature: str,
    output_image: str,
    *,
    top_taxa: int = 18,
    as_genera: bool = False,
    exclude_unclassified: bool = False,
    remove_zeros: bool = False,
    sort: Sequence[str] | None = None,
    focal_metadata: str | None = None,
    last_metadata: str | None = None,
    scaling: str = "original",
    units: str = "generic",
    dimensions: Sequence[float] = (11.0, 6.0),
    sample_order_output: str | None = None,
    humann_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the taxon-stratified contributions of one gene family or pathway across samples (``humann_barplot``).

    ``input_path`` is a stratified HUMAnN table (TSV or BIOM) and ``focal_feature`` the feature id to plot (the id, not the
    full name). The figure is written to ``output_image``; its extension selects the format (png, pdf, svg). ``top_taxa``
    is the largest number of taxa highlighted (default 18), ``as_genera`` collapses species to genera,
    ``exclude_unclassified`` drops the unclassified stratum, ``remove_zeros`` skips samples with a zero total for the
    feature. ``sort`` lists sample sorting methods in order (none, sum, dominant, braycurtis, braycurtis_w, metadata);
    ``focal_metadata`` and ``last_metadata`` use metadata rows of the table; ``scaling`` is ``original``, ``logstack`` or
    ``totalsum``; ``units`` names the y axis and ``dimensions`` is the width and height in inches. ``sample_order_output``
    saves the sample order. Result keys: ``image`` (the figure path, with its size in ``image_bytes``) and
    ``sample_order_file``. The plotted values are in the input table.
    """
    source = _cli.existing_file(input_path, "input_path")
    feature = _cli.text(focal_feature, "focal_feature")
    target = _cli.new_file(output_image, "output_image")
    if target.suffix.lower() not in (
        ".png",
        ".pdf",
        ".svg",
        ".jpg",
        ".jpeg",
        ".tif",
        ".tiff",
        ".eps",
    ):
        raise ConfigurationError(
            "output_image must end with an image extension such as .png, .pdf or .svg."
        )
    arguments = [
        "--input",
        str(source),
        "--focal-feature",
        feature,
        "--output",
        str(target),
        "--top-taxa",
        str(_cli.integer(top_taxa, "top_taxa", minimum=1, maximum=500)),
        "--scaling",
        _cli.choice(scaling, ("original", "logstack", "totalsum"), "scaling"),
        "--units",
        _cli.text(units, "units"),
    ]
    for value, flag_name, option in (
        (as_genera, "as_genera", "--as-genera"),
        (exclude_unclassified, "exclude_unclassified", "--exclude-unclassified"),
        (remove_zeros, "remove_zeros", "--remove-zeros"),
    ):
        if _cli.flag(value, flag_name):
            arguments.append(option)
    if sort:
        methods = list(sort)
        allowed = ("none", "sum", "dominant", "braycurtis", "braycurtis_w", "metadata")
        if any(item not in allowed for item in methods):
            raise ConfigurationError(f"sort must be a list of {list(allowed)}.")
        arguments += ["--sort", *methods]
    if focal_metadata is not None:
        arguments += ["--focal-metadata", _cli.text(focal_metadata, "focal_metadata")]
    if last_metadata is not None:
        arguments += ["--last-metadata", _cli.text(last_metadata, "last_metadata")]
    if len(list(dimensions)) != 2:
        raise ConfigurationError("dimensions must be two numbers: width and height in inches.")
    width, height = (
        _cli.number(item, "dimensions", minimum=1.0, maximum=100.0) for item in dimensions
    )
    arguments += ["--dimensions", format(width, ".12g"), format(height, ".12g")]
    expected = [target]
    order: Path | None = None
    if sample_order_output is not None:
        order = _cli.new_file(sample_order_output, "sample_order_output")
        arguments += ["--write-sample-order", str(order)]
        expected.append(order)
    result = _run(
        "humann_barplot",
        arguments,
        cwd=target.parent,
        expected=expected,
        humann_dir=humann_dir,
        operation="humann_barplot",
    )
    result.update(
        image=str(target),
        image_bytes=target.stat().st_size,
        sample_order_file=str(order) if order else None,
    )
    return result


__all__ = [
    "humann_barplot",
    "humann_build_custom_database",
    "humann_expand_cluster",
    "humann_infer_taxonomy",
    "humann_join_tables",
    "humann_quantify_alignments",
    "humann_reduce_table",
    "humann_regroup_table",
    "humann_rename_table",
    "humann_renorm_table",
    "humann_split_stratified_table",
    "humann_split_table",
    "humann_unpack_pathways",
]
_register_public_verbose(globals())
