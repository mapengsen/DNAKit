"""Centrifuge command-line facades: index building, index inspection, reference downloads and reference preparation.

Each public function runs one fixed Centrifuge program (https://github.com/DaehwanKimLab/centrifuge; 1.0.4 is the
tested version) with an explicit, validated argument list: ``centrifuge-build``, ``centrifuge-inspect``,
``centrifuge-download``, ``centrifuge-kreport`` and the Perl helpers ``centrifuge-compress.pl``,
``centrifuge-BuildSharedSequence.pl``, ``centrifuge-sort-nt.pl``, ``centrifuge-RemoveN.pl`` and
``centrifuge-RemoveEmptySequence.pl``. No shell is used. The user installs Centrifuge
(``conda install -c bioconda centrifuge``); the programs are looked up on PATH, or in ``centrifuge_dir`` when that
folder is given. A missing program raises ``BackendUnavailableError`` with an install hint. Results are
JSON-serialisable dictionaries with the status, backend, version, command, outputs, parsed tables and a bounded log
tail; outputs are written to new explicit paths and never overwrite existing files. Read classification with
``dnakit.microbiome.classify_centrifuge``, which was already available.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p12_cli as _cli
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_TOOL = _cli.Tool(
    name="centrifuge",
    executable="centrifuge",
    hint=(
        "Install Centrifuge 1.0.4 or later (for example `conda install -c bioconda centrifuge`, or build it from "
        "https://github.com/DaehwanKimLab/centrifuge), then put its programs on PATH or pass centrifuge_dir."
    ),
    tested_version="1.0.4",
    citation="https://doi.org/10.1101/gr.210641.116",
    version_arguments=("--version",),
    version_pattern=r"version\s+(\S+)",
)
_INDEX_SUFFIXES = (".1.cf", ".2.cf", ".3.cf", ".4.cf")
_DOMAINS = (
    "bacteria",
    "viral",
    "archaea",
    "fungi",
    "protozoa",
    "invertebrate",
    "plant",
    "vertebrate_mammalian",
    "vertebrate_other",
)
_ASSEMBLY_LEVELS = ("Complete Genome", "Chromosome", "Scaffold", "Contig", "Any")
_DOWNLOAD_DATABASES = ("refseq", "genbank", "contaminants", "taxonomy")
_INSPECT_OUTPUTS = {
    "fasta": [],
    "names": ["-n"],
    "summary": ["-s"],
    "conversion_table": ["--conversion-table"],
    "taxonomy_tree": ["--taxonomy-tree"],
    "name_table": ["--name-table"],
    "size_table": ["--size-table"],
    "estimate_n_kmers": ["--estimate-n-kmers"],
}
_SAFE_WORD = re.compile(r"^[A-Za-z0-9_.-]+$")
_SAFE_CATEGORY = re.compile(r"^[A-Za-z][A-Za-z ]*$")
_MAX_ROWS = 10_000


def _explicit(centrifuge_dir: str | None, name: str) -> str | None:
    """The path of ``name`` inside ``centrifuge_dir`` (None means: look it up on PATH)."""
    if centrifuge_dir is None:
        return None
    folder = Path(centrifuge_dir).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError(f"centrifuge_dir must be an existing folder: {folder}")
    return str(folder / name)


def _version_for_scripts(centrifuge_dir: str | None) -> str | None:
    """The Centrifuge version, read from the main program (the Perl helpers have no version option)."""
    try:
        executable = _cli.program(
            "centrifuge", explicit=_explicit(centrifuge_dir, "centrifuge"), hint=_TOOL.hint
        )
    except Exception:
        return None
    return _cli._version_of(executable, ("--version",), r"version\s+(\S+)")


def _finish(
    record: dict[str, Any], operation: str, centrifuge_dir: str | None, **extra: Any
) -> dict[str, Any]:
    """Add the operation name and, when the program has no version option of its own, the Centrifuge version."""
    if record.get("version") is None:
        record["version"] = _version_for_scripts(centrifuge_dir)
    return _cli.finish(record, operation, **extra)


def _index_prefix(value: str, *, must_exist: bool) -> Path:
    prefix = Path(_cli.free_text(value, "index_prefix")).expanduser().resolve()
    if "," in str(prefix):
        raise ConfigurationError("index_prefix must not contain commas.")
    exists = [prefix.with_name(prefix.name + suffix).is_file() for suffix in _INDEX_SUFFIXES]
    if must_exist and not all(exists):
        raise ConfigurationError(
            f"No complete Centrifuge index (.1.cf to .4.cf) was found for prefix {prefix}."
        )
    if not must_exist:
        if any(exists):
            raise ConfigurationError(f"A Centrifuge index already exists for prefix {prefix}.")
        if not prefix.parent.is_dir():
            raise ConfigurationError(f"The folder for index_prefix does not exist: {prefix.parent}")
    return prefix


def _tsv_rows(
    path: Path, names: Sequence[str], text_columns: frozenset[str], limit: int
) -> tuple[list[dict[str, Any]], int]:
    rows: list[dict[str, Any]] = []
    total = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            total += 1
            if len(rows) < limit:
                fields = [item.strip() for item in line.rstrip("\n").split("\t")]
                fields = [item for item in fields if item != "|"]
                rows.append(
                    {
                        name: (value if name in text_columns else _cli.convert_number(value))
                        for name, value in zip(names, fields, strict=False)
                    }
                )
    return rows, total


def _fasta_stats(path: Path) -> tuple[int, int]:
    records = bases = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def centrifuge_build(
    reference_files: Sequence[str],
    index_prefix: str,
    conversion_table: str,
    taxonomy_tree: str,
    *,
    name_table: str | None = None,
    size_table: str | None = None,
    threads: int = 1,
    offrate: int | None = None,
    ftab_chars: int | None = None,
    bmax: int | None = None,
    bmax_divisor: int | None = None,
    diff_cover_period: int | None = None,
    seed: int | None = None,
    kmer_count: int | None = None,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a Centrifuge classification index from reference FASTA files and a taxonomy (``centrifuge-build``).

    ``reference_files`` are DNA FASTA files; ``conversion_table`` maps every sequence id to a taxonomy id (two tab
    separated columns), ``taxonomy_tree`` is the NCBI ``nodes.dmp`` file and ``name_table`` (optional) the NCBI
    ``names.dmp`` file; ``size_table`` is an optional table of genome sizes. ``index_prefix`` is the path prefix of the
    four index files ``<prefix>.1.cf`` to ``<prefix>.4.cf`` (none may exist yet). ``threads`` is ``-p``, ``offrate`` ``-o``
    (the suffix array is sampled every 2^offrate characters), ``ftab_chars`` ``-t``, ``bmax``, ``bmax_divisor``
    (``--bmaxdivn``), ``diff_cover_period`` (``--dcv``) and ``seed`` are passed on unchanged. ``kmer_count=K`` also
    counts the distinct k-mers of the reference for k from 1 to K (``--kmer-count K``, K at most 64) and returns them in
    ``distinct_kmers`` (a mapping from k, as text, to the count). Result keys: ``index_prefix`` and ``index_files``.
    """
    files = _cli.existing_files(list(reference_files), "reference_files")
    if any("," in str(path) for path in files):
        raise ConfigurationError("reference file paths must not contain commas.")
    prefix = _index_prefix(index_prefix, must_exist=False)
    arguments = [
        "-p",
        str(_cli.threads_value(threads)),
        "--conversion-table",
        str(_cli.existing_file(conversion_table, "conversion_table")),
        "--taxonomy-tree",
        str(_cli.existing_file(taxonomy_tree, "taxonomy_tree")),
    ]
    if name_table is not None:
        arguments += ["--name-table", str(_cli.existing_file(name_table, "name_table"))]
    if size_table is not None:
        arguments += ["--size-table", str(_cli.existing_file(size_table, "size_table"))]
    for value, flag_name, option in (
        (offrate, "offrate", "-o"),
        (ftab_chars, "ftab_chars", "-t"),
        (bmax, "bmax", "--bmax"),
        (bmax_divisor, "bmax_divisor", "--bmaxdivn"),
        (diff_cover_period, "diff_cover_period", "--dcv"),
        (seed, "seed", "--seed"),
        (kmer_count, "kmer_count", "--kmer-count"),
    ):
        if value is not None:
            maximum = 64 if flag_name == "kmer_count" else None
            arguments += [
                option,
                str(
                    _cli.integer(
                        value, flag_name, minimum=0 if flag_name == "seed" else 1, maximum=maximum
                    )
                ),
            ]
    arguments += [",".join(str(path) for path in files), str(prefix)]
    expected = [prefix.with_name(prefix.name + suffix) for suffix in _INDEX_SUFFIXES]
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=prefix.parent,
        expected=expected,
        explicit=_explicit(centrifuge_dir, "centrifuge-build"),
        executable_name="centrifuge-build",
    )
    kmers = dict(re.findall(r"Number of distinct (\d+)-mers is (\d+)", record["log_tail"]))
    return _finish(
        record,
        "centrifuge_build",
        centrifuge_dir,
        index_prefix=str(prefix),
        index_files=[str(path) for path in expected],
        distinct_kmers={key: int(value) for key, value in kmers.items()},
    )


def centrifuge_inspect(
    index_prefix: str,
    output_path: str,
    *,
    output: str = "fasta",
    line_width: int = 60,
    max_records: int = _MAX_ROWS,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract information from a Centrifuge index (``centrifuge-inspect``).

    ``output`` selects what is written to ``output_path`` and parsed: ``fasta`` the indexed reference sequences as
    FASTA with ``line_width`` bases per line (``-a``, default 60; ``record_count`` and ``base_count`` are returned),
    ``names`` the reference sequence names (``-n``), ``summary`` the index parameters and the sequences with their
    lengths (``-s``; returned as ``parameters`` and ``records``), ``conversion_table`` the sequence-to-taxon table
    (``--conversion-table``), ``taxonomy_tree`` (``--taxonomy-tree``), ``name_table`` (``--name-table``),
    ``size_table`` the summed reference length per taxon (``--size-table``) and ``estimate_n_kmers`` the approximate
    number of distinct k-mers of the reference (``--estimate-n-kmers``; returned as ``estimated_kmers``). Tables are
    returned as ``records`` (at most ``max_records``) with ``record_count``.
    """
    prefix = _index_prefix(index_prefix, must_exist=True)
    target = _cli.new_file(output_path, "output_path")
    kind = _cli.choice(output, tuple(_INSPECT_OUTPUTS), "output")
    limit = _cli.integer(max_records, "max_records", minimum=0, maximum=10_000_000)
    arguments = list(_INSPECT_OUTPUTS[kind])
    if kind == "fasta":
        arguments += ["-a", str(_cli.integer(line_width, "line_width", minimum=1))]
    arguments.append(str(prefix))
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=target.parent,
        stdout_to=target,
        explicit=_explicit(centrifuge_dir, "centrifuge-inspect"),
        executable_name="centrifuge-inspect",
    )
    extra: dict[str, Any] = {"output_file": str(target), "output": kind}
    text_columns = frozenset({"name", "sequence"})
    if kind == "fasta":
        records, bases = _fasta_stats(target)
        extra.update(record_count=records, base_count=bases)
    elif kind == "names":
        names = [
            line
            for line in target.read_text(encoding="utf-8", errors="replace").splitlines()
            if line
        ]
        extra.update(records=names[:limit], record_count=len(names))
    elif kind == "summary":
        parameters: dict[str, Any] = {}
        sequences: list[dict[str, Any]] = []
        for line in target.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split("\t")
            if fields[0].startswith("Sequence-") and len(fields) >= 3:
                sequences.append(
                    {
                        "index": int(fields[0].split("-", 1)[1]),
                        "name": fields[1],
                        "length": int(fields[2]),
                    }
                )
            elif len(fields) == 2:
                parameters[fields[0]] = _cli.convert_number(fields[1])
        extra.update(parameters=parameters, records=sequences[:limit], record_count=len(sequences))
    elif kind == "estimate_n_kmers":
        match = re.search(
            r"kmers in the reference sequence:\s*(\d+)",
            target.read_text(encoding="utf-8", errors="replace"),
        )
        if match is None:
            raise BackendExecutionError("centrifuge-inspect did not report a k-mer estimate.")
        extra["estimated_kmers"] = int(match.group(1))
    else:
        columns = {
            "conversion_table": ["sequence", "taxid"],
            "taxonomy_tree": ["taxid", "parent_taxid", "rank"],
            "name_table": ["taxid", "name"],
            "size_table": ["taxid", "total_length"],
        }[kind]
        rows, total = _tsv_rows(target, columns, text_columns | {"rank"}, limit)
        extra.update(records=rows, record_count=total)
    return _finish(record, "centrifuge_inspect", centrifuge_dir, **extra)


def centrifuge_download(
    output_dir: str,
    database: str,
    *,
    domains: Sequence[str] | None = None,
    assembly_level: str = "Complete Genome",
    refseq_category: str | None = None,
    taxids: Sequence[int] | None = None,
    include_rna: bool = False,
    filter_unplaced: bool = False,
    mask_low_complexity: bool = False,
    taxid_in_headers: bool = False,
    threads: int = 1,
    download_program: str | None = None,
    mapping_file: str | None = None,
    timeout_seconds: float = 86_400.0,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download reference sequences or the NCBI taxonomy for a Centrifuge index (``centrifuge-download``; needs the network).

    ``database`` selects the source: ``refseq`` or ``genbank`` download genome FASTA files of the listed ``domains``
    (one or more of bacteria, viral, archaea, fungi, protozoa, invertebrate, plant, vertebrate_mammalian and
    vertebrate_other); ``contaminants`` downloads the UniVec and EmVec vector sequences, and ``taxonomy`` the NCBI
    ``nodes.dmp`` and ``names.dmp``. The options of the genome databases are ``assembly_level`` (``-a``; one of
    "Complete Genome", "Chromosome", "Scaffold", "Contig" or "Any"), ``refseq_category`` (``-c``), ``taxids`` (``-t``),
    ``include_rna`` to also download RNA sequences (``-r``), ``filter_unplaced`` (``-u``), ``mask_low_complexity`` to
    mask low-complexity regions with ``dustmasker`` (``-m``, needs BLAST+), ``taxid_in_headers`` to put the taxonomy id
    into every header (``-l``), ``threads`` (``-P``) and ``download_program`` (``-g``: rsync, curl or wget).
    ``output_dir`` must not exist; the program writes into it. Centrifuge prints the sequence-id to taxonomy-id table
    to standard output; it is saved as ``seqid_to_taxid.map`` in ``output_dir`` (``mapping_file`` may name another new
    path). Result keys: ``output_dir``, ``mapping_file``, ``mapping_rows`` (number of rows) and ``files`` (the
    downloaded files, at most 1000, as paths relative to ``output_dir``).
    """
    kind = _cli.choice(database, _DOWNLOAD_DATABASES, "database")
    directory = _cli.new_directory(output_dir, "output_dir")
    mapping = (
        _cli.new_file(mapping_file, "mapping_file")
        if mapping_file is not None
        else directory / "seqid_to_taxid.map"
    )
    arguments = [
        "-o",
        str(directory),
        "-P",
        str(_cli.integer(threads, "threads", minimum=1, maximum=64)),
    ]
    if kind in ("refseq", "genbank"):
        chosen = list(domains) if domains is not None else ["bacteria"]
        if not chosen or any(item not in _DOMAINS for item in chosen):
            raise ConfigurationError(f"domains must be a non-empty list of {list(_DOMAINS)}.")
        arguments += ["-d", ",".join(chosen)]
        arguments += ["-a", _cli.choice(assembly_level, _ASSEMBLY_LEVELS, "assembly_level")]
        if refseq_category is not None:
            if not isinstance(refseq_category, str) or not _SAFE_CATEGORY.match(refseq_category):
                raise ConfigurationError(
                    "refseq_category must be words made of letters and spaces."
                )
            arguments += ["-c", refseq_category]
        if taxids:
            arguments += [
                "-t",
                ",".join(str(_cli.integer(item, "taxids", minimum=1)) for item in taxids),
            ]
        for value, flag_name, option in (
            (include_rna, "include_rna", "-r"),
            (filter_unplaced, "filter_unplaced", "-u"),
            (mask_low_complexity, "mask_low_complexity", "-m"),
            (taxid_in_headers, "taxid_in_headers", "-l"),
        ):
            if _cli.flag(value, flag_name):
                arguments.append(option)
    elif (
        include_rna
        or filter_unplaced
        or mask_low_complexity
        or taxids
        or domains
        or refseq_category
    ):
        raise ConfigurationError("genome options only apply to database='refseq' or 'genbank'.")
    elif taxid_in_headers and kind == "taxonomy":
        raise ConfigurationError("taxid_in_headers does not apply to database='taxonomy'.")
    elif taxid_in_headers:
        arguments.append("-l")
    if download_program is not None:
        arguments += [
            "-g",
            _cli.choice(download_program, ("rsync", "curl", "wget"), "download_program"),
        ]
    arguments.append(kind)
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=directory,
        stdout_to=mapping,
        explicit=_explicit(centrifuge_dir, "centrifuge-download"),
        executable_name="centrifuge-download",
        timeout_seconds=timeout_seconds,
    )
    files = sorted(
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file() and path != mapping
    )
    rows = sum(
        1
        for line in mapping.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    )
    return _finish(
        record,
        "centrifuge_download",
        centrifuge_dir,
        database=kind,
        output_dir=str(directory),
        mapping_file=str(mapping),
        mapping_rows=rows,
        files=files[:1000],
        file_count=len(files),
    )


def centrifuge_kreport(
    index_prefix: str,
    classification_files: Sequence[str],
    output_path: str,
    *,
    min_score: int | None = None,
    min_length: int | None = None,
    no_lca: bool = False,
    show_zeros: bool = False,
    is_count_table: bool = False,
    max_records: int = _MAX_ROWS,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a Kraken-style taxonomic report from Centrifuge classification output (``centrifuge-kreport``).

    ``index_prefix`` is the index used for the classification (its taxonomy is read), ``classification_files`` the
    Centrifuge ``-S`` output files. ``min_score`` and ``min_length`` keep only hits with at least that score or hit
    length; ``no_lca`` reports count fractions at the taxa instead of the lowest common ancestor of multiple
    assignments; ``show_zeros`` also lists clades without reads; ``is_count_table`` reads files with the two columns
    ``taxid`` and ``count`` instead. The report is saved in ``output_path``. Result keys: ``records`` (``percent``,
    ``clade_reads``, ``taxon_reads``, ``rank`` code, ``taxid``, ``name`` and ``depth``; at most ``max_records``) and
    ``record_count``.
    """
    prefix = _index_prefix(index_prefix, must_exist=True)
    inputs = _cli.existing_files(list(classification_files), "classification_files")
    target = _cli.new_file(output_path, "output_path")
    limit = _cli.integer(max_records, "max_records", minimum=0, maximum=10_000_000)
    arguments = ["-x", str(prefix)]
    if min_score is not None:
        arguments += ["--min-score", str(_cli.integer(min_score, "min_score"))]
    if min_length is not None:
        arguments += ["--min-length", str(_cli.integer(min_length, "min_length", minimum=0))]
    for value, flag_name, option in (
        (no_lca, "no_lca", "--no-lca"),
        (show_zeros, "show_zeros", "--show-zeros"),
        (is_count_table, "is_count_table", "--is-count-table"),
    ):
        if _cli.flag(value, flag_name):
            arguments.append(option)
    arguments += [str(path) for path in inputs]
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=target.parent,
        stdout_to=target,
        explicit=_explicit(centrifuge_dir, "centrifuge-kreport"),
        executable_name="centrifuge-kreport",
    )
    rows: list[dict[str, Any]] = []
    total = 0
    for line in target.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) != 6:
            continue
        total += 1
        if len(rows) < limit:
            name = fields[5]
            rows.append(
                {
                    "percent": float(fields[0]),
                    "clade_reads": int(fields[1]),
                    "taxon_reads": int(fields[2]),
                    "rank": fields[3],
                    "taxid": int(fields[4]),
                    "name": name.strip(),
                    "depth": (len(name) - len(name.lstrip(" "))) // 2,
                }
            )
    return _finish(
        record,
        "centrifuge_kreport",
        centrifuge_dir,
        output_file=str(target),
        records=rows,
        record_count=total,
    )


def centrifuge_compress_genomes(
    download_dir: str,
    taxonomy_dir: str,
    output_dir: str,
    *,
    output_name: str = "compressed",
    level: str = "species",
    mapping_file: str | None = None,
    threads: int = 1,
    max_genome_size: int = 50_000_000,
    no_compress: bool = False,
    no_dustmasker: bool = False,
    centrifuge_dir: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compress a reference genome library by merging near-identical genomes of a taxonomic level (``centrifuge-compress.pl``).

    ``download_dir`` holds the genome FASTA files (``.fna``, ``.fa``, ``.ffn``, ``.fasta``), ``taxonomy_dir`` the NCBI
    ``nodes.dmp`` and ``names.dmp`` (and ``gi_taxid_nucl.dmp`` unless ``mapping_file``, a two-column file mapping
    sequence ids to taxonomy ids, is given). ``level`` is the taxonomic level at which genomes are merged (``-l``, default
    species), ``max_genome_size`` the largest genome that is compressed (``--maxG``), ``no_compress`` only collects and
    cleans the sequences (``-noCompress``) and ``no_dustmasker`` skips low-complexity masking (``-noDustmasker``). The
    Perl program needs Jellyfish 2 and MUMmer (``nucmer``) on PATH (and ``dustmasker`` for masking). ``output_dir`` must
    not exist; the program runs inside it and writes ``<output_name>.fa``, ``<output_name>.map`` and
    ``<output_name>.size`` there. Result keys: ``output_dir``, ``fasta``, ``taxid_map``, ``size_table``, ``sequence_count``
    and ``base_count``.
    """
    genomes = _cli.existing_directory(download_dir, "download_dir")
    taxonomy = _cli.existing_directory(taxonomy_dir, "taxonomy_dir")
    name = _cli.free_text(output_name, "output_name")
    if not _SAFE_WORD.match(name):
        raise ConfigurationError("output_name must contain only letters, digits, '_', '.' and '-'.")
    rank = _cli.free_text(level, "level")
    if not _SAFE_WORD.match(rank):
        raise ConfigurationError("level must be a rank name such as species or genus.")
    directory = _cli.new_directory(output_dir, "output_dir")
    arguments = [
        "-l",
        rank,
        "-o",
        name,
        "-t",
        str(_cli.threads_value(threads)),
        "--maxG",
        str(_cli.integer(max_genome_size, "max_genome_size", minimum=1)),
    ]
    if mapping_file is not None:
        arguments += ["-map", str(_cli.existing_file(mapping_file, "mapping_file"))]
    if _cli.flag(no_compress, "no_compress"):
        arguments.append("-noCompress")
    if _cli.flag(no_dustmasker, "no_dustmasker"):
        arguments.append("-noDustmasker")
    arguments += [str(genomes), str(taxonomy)]
    expected = [directory / f"{name}.fa", directory / f"{name}.map", directory / f"{name}.size"]
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=directory,
        expected=expected,
        explicit=_explicit(centrifuge_dir, "centrifuge-compress.pl"),
        executable_name="centrifuge-compress.pl",
        timeout_seconds=timeout_seconds,
    )
    records, bases = _fasta_stats(expected[0])
    return _finish(
        record,
        "centrifuge_compress_genomes",
        centrifuge_dir,
        output_dir=str(directory),
        fasta=str(expected[0]),
        taxid_map=str(expected[1]),
        size_table=str(expected[2]),
        sequence_count=records,
        base_count=bases,
    )


def centrifuge_shared_sequence(
    fasta_files: Sequence[str],
    output_dir: str,
    *,
    prefix: str = "shared",
    kmer_size: int = 53,
    kmer_portion: float = 0.01,
    nucmer_identity: int = 99,
    overlap: int = 250,
    fragment: bool = False,
    centrifuge_dir: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect sequence segments shared by several genomes and merge them (``centrifuge-BuildSharedSequence.pl``).

    ``fasta_files`` are the genomes of one taxon (at least two). The program counts k-mers with Jellyfish, picks the
    pair of genomes that share most k-mers and aligns them with MUMmer ``nucmer``; both must be on PATH. ``kmer_size`` is
    ``-kmerSize`` (default 53), ``kmer_portion`` the fraction of k-mers sampled (``-kmerPortion``, default 0.01),
    ``nucmer_identity`` the minimum percent identity of shared segments (``-nucmerIdy``, default 99), ``overlap`` the bases
    kept at the ends of a shared segment (``-overlap``, default 250) and ``fragment`` keeps the shared segments as separate
    fragments instead of one chunk (``-fragment``). ``output_dir`` must not exist; the program runs inside it. Result
    keys: ``output_dir``, ``files`` (the ``<prefix>_<n>.fa`` FASTA files that were written, with their ``records`` and
    ``bases``) and ``program_output`` (the progress lines).
    """
    paths = _cli.existing_files(list(fasta_files), "fasta_files", minimum=2)
    name = _cli.free_text(prefix, "prefix")
    if not _SAFE_WORD.match(name):
        raise ConfigurationError("prefix must contain only letters, digits, '_', '.' and '-'.")
    directory = _cli.new_directory(output_dir, "output_dir")
    listing = directory / f"{name}.file_list.txt"
    listing.write_text("".join(f"{path}\n" for path in paths), encoding="utf-8")
    arguments = [
        str(listing),
        "-prefix",
        name,
        "-kmerSize",
        str(_cli.integer(kmer_size, "kmer_size", minimum=1)),
        "-kmerPortion",
        format(_cli.number(kmer_portion, "kmer_portion", minimum=0.0, maximum=1.0), ".12g"),
        "-nucmerIdy",
        str(_cli.integer(nucmer_identity, "nucmer_identity", minimum=1, maximum=100)),
        "-overlap",
        str(_cli.integer(overlap, "overlap", minimum=0)),
    ]
    if _cli.flag(fragment, "fragment"):
        arguments.append("-fragment")
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=directory,
        explicit=_explicit(centrifuge_dir, "centrifuge-BuildSharedSequence.pl"),
        executable_name="centrifuge-BuildSharedSequence.pl",
        timeout_seconds=timeout_seconds,
    )
    produced = []
    for path in sorted(directory.glob(f"{name}_*.fa")):
        count, bases = _fasta_stats(path)
        produced.append({"file": str(path), "records": count, "bases": bases})
    record["outputs"] = [item["file"] for item in produced]
    return _finish(
        record,
        "centrifuge_shared_sequence",
        centrifuge_dir,
        output_dir=str(directory),
        files=produced,
        program_output=record["log_tail"],
    )


def centrifuge_sort_nt(
    sequence_file: str,
    mapping_files: Sequence[str],
    output_path: str,
    *,
    taxid_map_output: str | None = None,
    unmapped_output: str | None = None,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort the records of a nucleotide FASTA by taxonomy id (``centrifuge-sort-nt.pl``).

    ``sequence_file`` is a nucleotide FASTA whose headers start with an accession with version (for example the NCBI nt
    database), ``mapping_files`` are one or more NCBI ``accession2taxid`` files (plain or ``.gz``; four columns:
    accession, accession.version, taxid, gi). The sorted FASTA goes to ``output_path``. ``taxid_map_output`` (``-m``)
    receives the accession-to-taxid table of the sequences that were found and ``unmapped_output`` (``-a``) the file the
    upstream script writes when some accession has no mapping (Centrifuge 1.0.4 writes the accessions that do have a
    mapping into it, a defect of the upstream script that this facade does not alter; the file is not required to exist
    afterwards). Result keys: ``output_file``, ``record_count`` and ``base_count``.
    """
    source = _cli.existing_file(sequence_file, "sequence_file")
    maps = _cli.existing_files(list(mapping_files), "mapping_files")
    target = _cli.new_file(output_path, "output_path")
    arguments: list[str] = []
    extra_outputs: list[Path] = []
    if taxid_map_output is not None:
        extra = _cli.new_file(taxid_map_output, "taxid_map_output")
        arguments += ["-m", str(extra)]
        extra_outputs.append(extra)
    if unmapped_output is not None:
        extra = _cli.new_file(unmapped_output, "unmapped_output")
        arguments += ["-a", str(extra)]
    arguments += [str(source), *(str(path) for path in maps)]
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=target.parent,
        stdout_to=target,
        expected=extra_outputs,
        explicit=_explicit(centrifuge_dir, "centrifuge-sort-nt.pl"),
        executable_name="centrifuge-sort-nt.pl",
    )
    records, bases = _fasta_stats(target)
    return _finish(
        record,
        "centrifuge_sort_nt",
        centrifuge_dir,
        output_file=str(target),
        record_count=records,
        base_count=bases,
    )


def centrifuge_remove_n(
    input_path: str,
    output_path: str,
    *,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete the N bases of every FASTA record and re-wrap the sequence to 80 columns (``centrifuge-RemoveN.pl``).

    Both ``N`` and ``n`` are removed; other bases keep their case. A record that consists only of Ns keeps its header
    with no sequence line. Result keys: ``output_file``, ``record_count`` and ``base_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    record = _cli.invoke(
        _TOOL,
        [str(source)],
        cwd=target.parent,
        stdout_to=target,
        explicit=_explicit(centrifuge_dir, "centrifuge-RemoveN.pl"),
        executable_name="centrifuge-RemoveN.pl",
    )
    records, bases = _fasta_stats(target)
    return _finish(
        record,
        "centrifuge_remove_n",
        centrifuge_dir,
        output_file=str(target),
        record_count=records,
        base_count=bases,
    )


def centrifuge_remove_empty_sequences(
    input_path: str,
    output_path: str,
    *,
    centrifuge_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the FASTA records that have an empty sequence, as left behind by masking (``centrifuge-RemoveEmptySequence.pl``).

    Result keys: ``output_file``, ``record_count`` and ``base_count`` of the cleaned file.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    record = _cli.invoke(
        _TOOL,
        [str(source)],
        cwd=target.parent,
        stdout_to=target,
        explicit=_explicit(centrifuge_dir, "centrifuge-RemoveEmptySequence.pl"),
        executable_name="centrifuge-RemoveEmptySequence.pl",
    )
    records, bases = _fasta_stats(target)
    return _finish(
        record,
        "centrifuge_remove_empty_sequences",
        centrifuge_dir,
        output_file=str(target),
        record_count=records,
        base_count=bases,
    )


__all__ = [
    "centrifuge_build",
    "centrifuge_compress_genomes",
    "centrifuge_download",
    "centrifuge_inspect",
    "centrifuge_kreport",
    "centrifuge_remove_empty_sequences",
    "centrifuge_remove_n",
    "centrifuge_shared_sequence",
    "centrifuge_sort_nt",
]
_register_public_verbose(globals())
