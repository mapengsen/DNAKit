"""Kraken 2 command-line facades (tested with Kraken 2.17.2 and its ``k2`` wrapper): read classification with the
per-taxon report, database inspection, merging of classifications from several databases, and every step of database
construction (download a library, add your own sequences, build, standard database, special 16S/GTDB databases).

The user installs Kraken 2 (``conda install -c bioconda kraken2``, or from https://github.com/DerrickWood/kraken2)
and puts it on PATH, or passes ``kraken2_path`` (the folder that holds ``kraken2``, ``kraken2-build``,
``kraken2-inspect`` and ``k2``). Each public function runs one fixed Kraken 2 program with an argument list built
from typed, validated parameters: no shell, no generic runner, outputs are new files. Only the options that the
caller sets are passed, so every other option keeps Kraken 2's own default. Results are JSON-serialisable dictionaries
(status, backend, version, command, return code, outputs, log tail).

The download and build functions contact NCBI (and for GTDB the GTDB server) and can take hours and large disk space;
ordinary tests never call them.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = (
    "Install Kraken 2 (for example `conda install -c bioconda kraken2`, or from "
    "https://github.com/DerrickWood/kraken2), then put it on PATH or pass kraken2_path (the folder that "
    "contains kraken2, kraken2-build, kraken2-inspect and k2)."
)
_BACKENDS = ("kraken2-build", "k2")
_MERGE_POLICIES = ("first-classified", "last-classified", "resolve-lca")
_SPECIAL = ("greengenes", "silva", "rdp", "gtdb")
_LIBRARY_NAME = re.compile(r"^[A-Za-z0-9_.:-]+$")
_ASSEMBLY_LEVELS = ("chromosome", "complete_genome", "scaffold", "contig")
_ASSEMBLY_SOURCES = ("refseq", "genbank", "all")
_PROCESSED = re.compile(r"(\d+) sequences(?: \([\d.]+ Mbp\))? processed")
_CLASSIFIED = re.compile(r"(\d+) sequences classified \(([\d.]+)%\)")
_UNCLASSIFIED = re.compile(r"(\d+) sequences unclassified")
_HASH_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")
_MAX_ROWS = 100_000


def _executable(name: str, kraken2_path: str | None) -> str | None:
    if kraken2_path is None:
        return None
    folder = Path(kraken2_path).expanduser()
    if not folder.is_dir():
        raise ConfigurationError(f"kraken2_path must be a folder that contains the Kraken 2 programs: {folder}")
    return str(folder / name)


def _run(
    program: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    kraken2_path: str | None,
    timeout_seconds: float,
    stdout_file: Path | None = None,
) -> dict[str, Any]:
    return _p14_cli.run(
        program,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=_executable(program, kraken2_path),
        hint=_HINT,
        stdout_file=stdout_file,
        timeout_seconds=timeout_seconds,
    )


def _timeout(value: float) -> float:
    return _gap_cli.number(value, "timeout_seconds", minimum=1.0, maximum=86_400.0)


def _database_dir(value: str, name: str = "database") -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing Kraken 2 database folder: {path}")
    return path


def _built_database(value: str, name: str = "database") -> Path:
    path = _database_dir(value, name)
    missing = [item for item in _HASH_FILES if not (path / item).is_file()]
    if missing:
        raise ConfigurationError(f"{name} is not a built Kraken 2 database (missing {', '.join(missing)}): {path}")
    return path


def _database_list(value: str | Sequence[str], name: str = "databases") -> list[Path]:
    items = [value] if isinstance(value, str) else list(value)
    if not items:
        raise ConfigurationError(f"{name} must name at least one database.")
    paths = [_built_database(item, name) for item in items]
    if any("," in str(path) for path in paths):
        raise ConfigurationError("Database paths must not contain commas.")
    return paths


def _build_dir(value: str, name: str = "database") -> Path:
    """A database folder to build in: it is created when missing and must have an existing parent."""
    path = Path(value).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    path.mkdir(exist_ok=True)
    return path


def _report_rows(path: Path, *, minimizer_data: bool) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    columns = ["percent", "clade_reads", "taxon_reads"]
    if minimizer_data:
        columns += ["minimizers", "distinct_minimizers"]
    columns += ["rank", "taxid", "name"]
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns):
                continue
            row = dict(zip(columns, fields, strict=True))
            row["name"] = row["name"].strip()
            rows.append(row)
            if len(rows) >= _MAX_ROWS:
                break
    return rows


def _summary(result: dict[str, Any]) -> None:
    log = result["log_tail"]
    processed, classified, unclassified = _PROCESSED.search(log), _CLASSIFIED.search(log), _UNCLASSIFIED.search(log)
    if processed is not None:
        result["sequences_processed"] = int(processed.group(1))
    if classified is not None:
        result["sequences_classified"] = int(classified.group(1))
        result["percent_classified"] = float(classified.group(2))
    if unclassified is not None:
        result["sequences_unclassified"] = int(unclassified.group(1))


def kraken2_classify_reads(
    reads: Sequence[str],
    databases: str | Sequence[str],
    output: str,
    *,
    report: str | None = None,
    paired: bool = False,
    interleaved: bool = False,
    quick: bool = False,
    confidence: float | None = None,
    minimum_hit_groups: int | None = None,
    minimum_base_quality: int | None = None,
    merge_policy: str | None = None,
    threads: int | None = None,
    use_names: bool = False,
    memory_mapping: bool = False,
    report_minimizer_data: bool = False,
    report_zero_counts: bool = False,
    use_mpa_style: bool = False,
    classified_out: str | None = None,
    unclassified_out: str | None = None,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify reads against one or more Kraken 2 databases (``kraken2``; ``k2 classify`` for several databases).

    ``databases`` is a built database folder or a list of them; with several, ``k2 classify`` assigns each read the
    taxon chosen by ``merge_policy`` ('resolve-lca' is the k2 default and takes the lowest common ancestor of the
    per-database calls, 'first-classified' or 'last-classified'). A protein database (built with
    ``kraken2_build_database(protein=True)``) classifies DNA reads by six-frame translated search with the same
    call. ``quick`` stops at the first hit, ``confidence`` (0-1), ``minimum_hit_groups`` and ``minimum_base_quality``
    tune the call. ``paired`` means the files hold mates (two files per sample, or one with ``interleaved``, which
    needs ``k2``). ``report`` writes the per-taxon clade-count report (Kraken report format; ``use_mpa_style``,
    ``report_zero_counts`` and ``report_minimizer_data`` shape it) and is parsed into ``report_rows``. ``output``
    receives the per-read classification (taxid with names when ``use_names``). The counts of processed, classified
    and unclassified sequences are read from the program's log. Writes new files.
    """
    if isinstance(reads, str) or not reads:
        raise ConfigurationError("reads must be a non-empty list of FASTA/FASTQ files.")
    read_paths = [str(_gap_cli.existing_file(item, "reads")) for item in reads]
    database_paths = _database_list(databases)
    output_path = _gap_cli.new_file(output, "output")
    multiple = len(database_paths) > 1
    use_k2 = multiple or merge_policy is not None or interleaved
    program = "k2" if use_k2 else "kraken2"
    arguments = ["classify"] if use_k2 else []
    arguments += ["--db", ",".join(str(path) for path in database_paths), "--output", str(output_path)]
    expected = [output_path]
    if report is not None:
        report_path = _gap_cli.new_file(report, "report")
        arguments += ["--report", str(report_path)]
        expected.append(report_path)
    for value, flag_name, name in (
        (classified_out, "--classified-out", "classified_out"),
        (unclassified_out, "--unclassified-out", "unclassified_out"),
    ):
        if value is not None:
            if not Path(value).expanduser().resolve().parent.is_dir():
                raise ConfigurationError(f"The folder for {name} does not exist.")
            arguments += [flag_name, str(_gap_cli.new_file(value, name))]
    if merge_policy is not None:
        arguments += ["--merge-policy", _gap_cli.choice(merge_policy, _MERGE_POLICIES, "merge_policy")]
    if _p14_cli.flag(paired, "paired") or _p14_cli.flag(interleaved, "interleaved"):
        if paired and interleaved:
            raise ConfigurationError("Choose paired or interleaved, not both.")
        arguments.append("--paired" if paired else "--interleaved")
    if confidence is not None:
        arguments += ["--confidence", str(_gap_cli.number(confidence, "confidence", minimum=0.0, maximum=1.0))]
    if minimum_hit_groups is not None:
        arguments += ["--minimum-hit-groups", str(_gap_cli.integer(minimum_hit_groups, "minimum_hit_groups", minimum=1))]
    if minimum_base_quality is not None:
        arguments += ["--minimum-base-quality", str(_gap_cli.integer(minimum_base_quality, "minimum_base_quality", minimum=0, maximum=93))]
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    for value, flag_name, name in (
        (quick, "--quick", "quick"),
        (use_names, "--use-names", "use_names"),
        (memory_mapping, "--memory-mapping", "memory_mapping"),
        (report_minimizer_data, "--report-minimizer-data", "report_minimizer_data"),
        (report_zero_counts, "--report-zero-counts", "report_zero_counts"),
        (use_mpa_style, "--use-mpa-style", "use_mpa_style"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if report is None and (report_minimizer_data or report_zero_counts or use_mpa_style):
        raise ConfigurationError("report_minimizer_data, report_zero_counts and use_mpa_style need report.")
    arguments += read_paths
    result = _run(
        program, arguments, cwd=output_path.parent, expected=expected, kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )
    result["databases"] = [str(path) for path in database_paths]
    _summary(result)
    if report is not None and result["status"] == "completed" and not use_mpa_style:
        result["report_rows"] = _report_rows(Path(report).expanduser().resolve(), minimizer_data=report_minimizer_data)
    return result


def kraken2_inspect_database(
    database: str,
    output: str,
    *,
    skip_counts: bool = False,
    use_mpa_style: bool = False,
    report_zero_counts: bool = False,
    memory_mapping: bool = False,
    kraken2_path: str | None = None,
    timeout_seconds: float = 3_600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a Kraken 2 database and its per-taxon minimizer counts (``k2 inspect``, the same report as ``kraken2-inspect``).

    The output file starts with ``#`` lines (database path, k-mer and minimizer length, spaced mask, number of taxonomy
    nodes, hash table size) followed, unless ``skip_counts``, by a report with one line per taxon: percent of
    minimizers, minimizers in the clade, minimizers at the taxon, rank, taxid and name. ``use_mpa_style`` formats it
    like kraken-mpa-report, ``report_zero_counts`` also lists taxa without minimizers and ``memory_mapping`` avoids
    loading the database into RAM. ``header`` and ``rows`` hold the parsed text. Writes a new file.
    """
    db = _built_database(database)
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["inspect", "--db", str(db), "--output", str(output_path)]
    for value, flag_name, name in (
        (skip_counts, "--skip-counts", "skip_counts"),
        (use_mpa_style, "--use-mpa-style", "use_mpa_style"),
        (report_zero_counts, "--report-zero-counts", "report_zero_counts"),
        (memory_mapping, "--memory-mapping", "memory_mapping"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    result = _run(
        "k2", arguments, cwd=output_path.parent, expected=[output_path], kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )
    if result["status"] == "completed":
        lines = output_path.read_text(encoding="utf-8", errors="replace").splitlines()
        result["header"] = [line[1:].strip() for line in lines if line.startswith("#")]
        if not use_mpa_style:
            rows = []
            for line in lines:
                fields = line.split("\t")
                if not line.startswith("#") and len(fields) == 6:
                    rows.append(dict(zip(("percent", "clade_minimizers", "taxon_minimizers", "rank", "taxid", "name"), fields, strict=True)))
                if len(rows) >= _MAX_ROWS:
                    break
            for row in rows:
                row["name"] = row["name"].strip()
                row["percent"] = row["percent"].strip()
            result["rows"] = rows
    return result


def kraken2_merge_classifications(
    classification_files: Sequence[str],
    databases: str | Sequence[str],
    output: str,
    *,
    report: str | None = None,
    use_names: bool = False,
    use_mpa_style: bool = False,
    report_zero_counts: bool = False,
    threads: int | None = None,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the per-read classification files of the same reads made against different databases (``k2 merge-classifications``).

    Each read gets the lowest common ancestor of the taxa assigned in the input files. ``databases`` are the
    databases that produced the files (their taxonomy is needed for the LCA); ``report`` writes the clade-count
    report of the merged result. Writes new files. Requires Kraken 2 with the ``k2`` wrapper.
    """
    if isinstance(classification_files, str) or len(classification_files) < 2:
        raise ConfigurationError("classification_files must list at least two classification files.")
    files = [str(_gap_cli.existing_file(item, "classification_files")) for item in classification_files]
    database_paths = _database_list(databases)
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["merge-classifications", "--db", ",".join(str(path) for path in database_paths), "--output", str(output_path)]
    expected = [output_path]
    if report is not None:
        report_path = _gap_cli.new_file(report, "report")
        arguments += ["--report", str(report_path)]
        expected.append(report_path)
    elif use_mpa_style or report_zero_counts:
        raise ConfigurationError("use_mpa_style and report_zero_counts need report.")
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    for value, flag_name, name in (
        (use_names, "--use-names", "use_names"),
        (use_mpa_style, "--use-mpa-style", "use_mpa_style"),
        (report_zero_counts, "--report-zero-counts", "report_zero_counts"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    arguments += files
    return _run(
        "k2", arguments, cwd=output_path.parent, expected=expected, kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )


def _backend(value: str) -> str:
    return _gap_cli.choice(value, _BACKENDS, "backend")


def _library_outputs(result: dict[str, Any], folder: Path) -> dict[str, Any]:
    library = folder / "library"
    result["outputs"] = [str(library)] if library.is_dir() else []
    if result["status"] == "completed" and not library.is_dir():
        result["status"] = "failed"
    return result


def kraken2_download_library(
    database: str,
    library: str,
    *,
    backend: str = "kraken2-build",
    protein: bool = False,
    no_masking: bool = False,
    use_ftp: bool = False,
    threads: int | None = None,
    assembly_source: str | None = None,
    assembly_levels: Sequence[str] | None = None,
    has_annotation: bool = False,
    resume: bool = False,
    blast_volumes: str | None = None,
    masker_threads: int | None = None,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download NCBI genomic sequences into the ``library`` folder of a Kraken 2 database.

    With ``backend='kraken2-build'`` (``kraken2-build --download-library``) ``library`` is one of archaea, bacteria,
    plasmid, viral, human, fungi, plant, protozoa, nr, nt, UniVec, UniVec_Core; ``use_ftp`` downloads by FTP instead
    of rsync. With ``backend='k2'`` (``k2 download-library``) ``library`` may also be an NCBI taxonomy id, a
    BioProject or an assembly accession, and ``assembly_source`` ('refseq', 'genbank', 'all'), ``assembly_levels``
    (chromosome, complete_genome, scaffold, contig), ``has_annotation``, ``resume``, ``blast_volumes`` (nr/nt volumes
    such as ``'1,3,5-8'``) and ``masker_threads`` apply. Low-complexity masking is on unless ``no_masking``
    (it needs dustmasker/segmasker or k2mask). The taxonomy must be downloaded separately before building. Contacts
    NCBI over the network and can take hours.
    """
    folder = _build_dir(database)
    name = _gap_cli.text(library, "library")
    if not _LIBRARY_NAME.fullmatch(name):
        raise ConfigurationError("library may contain letters, digits and _ . : - only.")
    chosen = _backend(backend)
    k2_only = {
        "assembly_source": assembly_source, "assembly_levels": assembly_levels, "blast_volumes": blast_volumes,
        "masker_threads": masker_threads,
    }
    if chosen == "kraken2-build":
        if any(value is not None for value in k2_only.values()) or has_annotation or resume:
            raise ConfigurationError("assembly_source, assembly_levels, has_annotation, resume, blast_volumes and masker_threads need backend='k2'.")
        arguments = ["--download-library", name, "--db", str(folder)]
        if _p14_cli.flag(use_ftp, "use_ftp"):
            arguments.append("--use-ftp")
    else:
        if use_ftp:
            raise ConfigurationError("use_ftp applies to backend='kraken2-build' only.")
        arguments = ["download-library", "--db", str(folder), "--library", name]
        if assembly_source is not None:
            arguments += ["--assembly-source", _gap_cli.choice(assembly_source, _ASSEMBLY_SOURCES, "assembly_source")]
        if assembly_levels is not None:
            if isinstance(assembly_levels, str) or not assembly_levels:
                raise ConfigurationError("assembly_levels must be a non-empty list.")
            arguments += ["--assembly-levels", *(_gap_cli.choice(item, _ASSEMBLY_LEVELS, "assembly_levels") for item in assembly_levels)]
        if _p14_cli.flag(has_annotation, "has_annotation"):
            arguments.append("--has-annotation")
        if _p14_cli.flag(resume, "resume"):
            arguments.append("--resume")
        if blast_volumes is not None:
            if not re.fullmatch(r"[0-9,.:-]+", blast_volumes):
                raise ConfigurationError("blast_volumes must look like '1,3,5-8'.")
            arguments += ["--blast-volumes", blast_volumes]
        if masker_threads is not None:
            if no_masking:
                raise ConfigurationError("masker_threads cannot be combined with no_masking.")
            arguments += ["--masker-threads", str(_gap_cli.integer(masker_threads, "masker_threads", minimum=1, maximum=256))]
    if _p14_cli.flag(protein, "protein"):
        arguments.append("--protein")
    if _p14_cli.flag(no_masking, "no_masking"):
        arguments.append("--no-masking")
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    result = _run(
        chosen, arguments, cwd=folder, expected=[], kraken2_path=kraken2_path, timeout_seconds=_timeout(timeout_seconds)
    )
    return _library_outputs(result, folder)


def kraken2_add_to_library(
    database: str,
    fasta_files: Sequence[str],
    *,
    backend: str = "kraken2-build",
    protein: bool = False,
    no_masking: bool = False,
    threads: int | None = None,
    masker_threads: int | None = None,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add your own FASTA sequences to the library of a Kraken 2 database (``kraken2-build --add-to-library``; ``k2 add-to-library``).

    Every sequence needs a taxonomy id, either in the header as ``>name|kraken:taxid|<taxid>`` or through an NCBI
    accession listed in the database's ``seqid2taxid.map``. Low-complexity sequence is masked unless ``no_masking``
    (needs dustmasker/segmasker, or k2mask for ``backend='k2'``). With ``backend='kraken2-build'`` one command is run
    per file and the result lists them in ``runs``; with ``backend='k2'`` one command adds all files and
    ``threads``/``masker_threads`` apply. The database folder is created when missing.
    """
    folder = _build_dir(database)
    if isinstance(fasta_files, str) or not fasta_files:
        raise ConfigurationError("fasta_files must be a non-empty list of FASTA files.")
    files = [str(_gap_cli.existing_file(item, "fasta_files")) for item in fasta_files]
    chosen = _backend(backend)
    timeout = _timeout(timeout_seconds)
    if chosen == "kraken2-build":
        if threads is not None or masker_threads is not None:
            raise ConfigurationError("threads and masker_threads for add-to-library need backend='k2'.")
        runs = []
        for item in files:
            arguments = ["--add-to-library", item, "--db", str(folder)]
            if _p14_cli.flag(protein, "protein"):
                arguments.append("--protein")
            if _p14_cli.flag(no_masking, "no_masking"):
                arguments.append("--no-masking")
            run = _run(chosen, arguments, cwd=folder, expected=[], kraken2_path=kraken2_path, timeout_seconds=timeout)
            runs.append(run)
            if run["status"] != "completed":
                break
        result = dict(runs[-1])
        result["runs"] = runs
        result["status"] = "completed" if len(runs) == len(files) and all(r["status"] == "completed" for r in runs) else "failed"
    else:
        arguments = ["add-to-library", "--db", str(folder), "--file", *files]
        if _p14_cli.flag(protein, "protein"):
            arguments.append("--protein")
        if _p14_cli.flag(no_masking, "no_masking"):
            arguments.append("--no-masking")
        elif masker_threads is not None:
            arguments += ["--masker-threads", str(_gap_cli.integer(masker_threads, "masker_threads", minimum=1, maximum=256))]
        if threads is not None:
            arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
        result = _run(chosen, arguments, cwd=folder, expected=[], kraken2_path=kraken2_path, timeout_seconds=timeout)
    return _library_outputs(result, folder)


def _kmer_options(
    kmer_length: int | None,
    minimizer_length: int | None,
    minimizer_spaces: int | None,
    load_factor: float | None,
    max_db_size: int | None,
    threads: int | None,
    fast_build: bool,
    protein: bool,
) -> list[str]:
    arguments: list[str] = []
    for name, value, flag_name in (
        ("kmer_length", kmer_length, "--kmer-len"),
        ("minimizer_length", minimizer_length, "--minimizer-len"),
        ("max_db_size", max_db_size, "--max-db-size"),
    ):
        if value is not None:
            arguments += [flag_name, str(_gap_cli.integer(value, name, minimum=1))]
    if minimizer_spaces is not None:
        arguments += ["--minimizer-spaces", str(_gap_cli.integer(minimizer_spaces, "minimizer_spaces", minimum=0))]
    if load_factor is not None:
        value = _gap_cli.number(load_factor, "load_factor", minimum=0.0, maximum=1.0)
        if value <= 0.0:
            raise ConfigurationError("load_factor must be greater than 0.")
        arguments += ["--load-factor", str(value)]
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if _p14_cli.flag(fast_build, "fast_build"):
        arguments.append("--fast-build")
    if _p14_cli.flag(protein, "protein"):
        arguments.append("--protein")
    return arguments


def kraken2_build_database(
    database: str,
    *,
    backend: str = "kraken2-build",
    kmer_length: int | None = None,
    minimizer_length: int | None = None,
    minimizer_spaces: int | None = None,
    load_factor: float | None = None,
    max_db_size: int | None = None,
    threads: int | None = None,
    fast_build: bool = False,
    protein: bool = False,
    parts: int | None = None,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the minimizer hash and taxonomy index of a Kraken 2 database from its library (``kraken2-build --build``; ``k2 build``).

    The database folder must already hold ``taxonomy/`` (``names.dmp`` and ``nodes.dmp``) and a ``library/`` with
    sequences (see kraken2_download_library and kraken2_add_to_library). ``kmer_length`` and ``minimizer_length``
    default to Kraken 2's 35/31 bases (15/12 amino acids with ``protein=True``, which builds a database for
    six-frame translated search of DNA reads); ``minimizer_spaces`` defaults to 7 (0 for protein). ``load_factor`` is
    the hash table fill (default 0.7), ``max_db_size`` caps the table in bytes (the library is down-sampled to fit),
    ``fast_build`` skips deterministic multi-threaded building. ``parts`` (k2 only) builds a database that was split
    into parts. Writes ``hash.k2d``, ``opts.k2d`` and ``taxo.k2d`` in the database folder.
    """
    folder = _database_dir(database)
    chosen = _backend(backend)
    expected = [folder / item for item in _HASH_FILES]
    if any(path.exists() for path in expected):
        raise ConfigurationError(f"The database is already built: {folder}")
    options = _kmer_options(kmer_length, minimizer_length, minimizer_spaces, load_factor, max_db_size, threads, fast_build, protein)
    if chosen == "kraken2-build":
        if parts is not None:
            raise ConfigurationError("parts needs backend='k2'.")
        arguments = ["--build", "--db", str(folder), *options]
    else:
        arguments = ["build", "--db", str(folder), *options]
        if parts is not None:
            arguments += ["--parts", str(_gap_cli.integer(parts, "parts", minimum=1))]
    return _run(
        chosen, arguments, cwd=folder, expected=expected, kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )


def kraken2_build_standard_database(
    database: str,
    *,
    backend: str = "kraken2-build",
    no_masking: bool = False,
    use_ftp: bool = False,
    kmer_length: int | None = None,
    minimizer_length: int | None = None,
    minimizer_spaces: int | None = None,
    load_factor: float | None = None,
    max_db_size: int | None = None,
    threads: int | None = None,
    fast_build: bool = False,
    protein: bool = False,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and build the standard Kraken 2 database (``kraken2-build --standard``; ``k2 build --standard``).

    The standard database holds RefSeq archaea, bacteria, viral, plasmid, human and UniVec_Core sequences with the
    NCBI taxonomy; the download is tens of gigabytes and the build needs tens of gigabytes of memory and hours of
    time. The other parameters are as in kraken2_build_database; ``use_ftp`` (kraken2-build) downloads by FTP.
    Writes ``hash.k2d``, ``opts.k2d`` and ``taxo.k2d`` in the database folder.
    """
    folder = _build_dir(database)
    chosen = _backend(backend)
    expected = [folder / item for item in _HASH_FILES]
    if any(path.exists() for path in expected):
        raise ConfigurationError(f"The database is already built: {folder}")
    options = _kmer_options(kmer_length, minimizer_length, minimizer_spaces, load_factor, max_db_size, threads, fast_build, protein)
    if _p14_cli.flag(no_masking, "no_masking"):
        options.append("--no-masking")
    if chosen == "kraken2-build":
        arguments = ["--standard", "--db", str(folder), *options]
        if _p14_cli.flag(use_ftp, "use_ftp"):
            arguments.append("--use-ftp")
    else:
        if use_ftp:
            raise ConfigurationError("use_ftp applies to backend='kraken2-build' only.")
        arguments = ["build", "--standard", "--db", str(folder), *options]
    return _run(
        chosen, arguments, cwd=folder, expected=expected, kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )


def kraken2_build_special_database(
    database: str,
    special: str,
    *,
    backend: str = "kraken2-build",
    kmer_length: int | None = None,
    minimizer_length: int | None = None,
    minimizer_spaces: int | None = None,
    load_factor: float | None = None,
    max_db_size: int | None = None,
    threads: int | None = None,
    fast_build: bool = False,
    gtdb_files: str | None = None,
    gtdb_use_ncbi_taxonomy: bool = False,
    gtdb_representatives_only: bool = False,
    gtdb_split_database: bool = False,
    gtdb_server: str | None = None,
    no_masking: bool = False,
    kraken2_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and build a special reference database (``kraken2-build --special``; ``k2 build --special``).

    ``special`` is 'greengenes', 'silva' or 'rdp' (16S rRNA databases; RDP's download URLs no longer work upstream)
    or 'gtdb' (genomes of the Genome Taxonomy Database; ``backend='k2'`` only). For GTDB, ``gtdb_files`` names the
    genome files (or a pattern) to use instead of downloading, ``gtdb_use_ncbi_taxonomy`` uses NCBI taxonomy ids,
    ``gtdb_representatives_only`` keeps the representative genomes, ``gtdb_split_database`` splits the database into
    parts and ``gtdb_server`` selects the download server. The k-mer, minimizer, load-factor, size and thread options
    are as in kraken2_build_database. Writes ``hash.k2d``, ``opts.k2d`` and ``taxo.k2d`` in the database folder. Contacts
    the data servers over the network and can take hours.
    """
    folder = _build_dir(database)
    chosen = _backend(backend)
    kind = _gap_cli.choice(special, _SPECIAL, "special")
    expected = [folder / item for item in _HASH_FILES]
    if any(path.exists() for path in expected):
        raise ConfigurationError(f"The database is already built: {folder}")
    gtdb_options: list[str] = []
    if gtdb_files is not None:
        gtdb_options += ["--gtdb-files", _gap_cli.text(gtdb_files, "gtdb_files")]
    if gtdb_server is not None:
        gtdb_options += ["--gtdb-server", _gap_cli.text(gtdb_server, "gtdb_server")]
    for value, flag_name, name in (
        (gtdb_use_ncbi_taxonomy, "--gtdb-use-ncbi-taxonomy", "gtdb_use_ncbi_taxonomy"),
        (gtdb_representatives_only, "--gtdb-representatives-only", "gtdb_representatives_only"),
        (gtdb_split_database, "--gtdb-split-database", "gtdb_split_database"),
    ):
        if _p14_cli.flag(value, name):
            gtdb_options.append(flag_name)
    if (kind == "gtdb" or gtdb_options) and chosen != "k2":
        raise ConfigurationError("The GTDB database and the gtdb_* options need backend='k2'.")
    if gtdb_options and kind != "gtdb":
        raise ConfigurationError("The gtdb_* options apply to special='gtdb' only.")
    options = _kmer_options(kmer_length, minimizer_length, minimizer_spaces, load_factor, max_db_size, threads, fast_build, False)
    if _p14_cli.flag(no_masking, "no_masking"):
        options.append("--no-masking")
    if chosen == "kraken2-build":
        arguments = ["--special", kind, "--db", str(folder), *options]
    else:
        arguments = ["build", "--special", kind, "--db", str(folder), *gtdb_options, *options]
    return _run(
        chosen, arguments, cwd=folder, expected=expected, kraken2_path=kraken2_path,
        timeout_seconds=_timeout(timeout_seconds),
    )


__all__ = [
    "kraken2_add_to_library",
    "kraken2_build_database",
    "kraken2_build_special_database",
    "kraken2_build_standard_database",
    "kraken2_classify_reads",
    "kraken2_download_library",
    "kraken2_inspect_database",
    "kraken2_merge_classifications",
]
_register_public_verbose(globals())
