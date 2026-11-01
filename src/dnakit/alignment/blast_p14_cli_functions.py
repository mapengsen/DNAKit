"""BLAST+ command-line facades (tested with BLAST+ 2.17.0): nucleotide search with every ``blastn`` task, translated
search against a conserved-domain database, sequence retrieval from a BLAST database, MegaBLAST indexing, WindowMasker
repeat masking, lower-case-mask conversion and searches of SRA run databases.

The user installs BLAST+ (``conda install -c bioconda blast``, or the NCBI release archive) and puts it on PATH, or
passes ``blast_path`` (the folder that holds the BLAST+ programs). Each public function runs one fixed BLAST+ program
with an argument list built from typed, validated parameters: no shell, no generic runner, outputs are new files.
Only the options that the caller sets are passed, so every other option keeps BLAST+'s own default. Results are
JSON-serialisable dictionaries (status, backend, version, command, return code, outputs, log tail, small parsed rows).
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
    "Install BLAST+ 2.12 or later (for example `conda install -c bioconda blast`, or the archive from "
    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), then put it on PATH or pass blast_path "
    "(the folder that contains the BLAST+ programs)."
)
_VERSION_ARGS = ("-version",)
_NUCLEOTIDE_TASKS = ("blastn", "blastn-short", "dc-megablast", "megablast", "rmblastn")
_STRANDS = ("both", "plus", "minus")
_OUTFMT = re.compile(r"^\d{1,2}(?: [A-Za-z_]+)*$")
_STANDARD_COLUMNS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore",
)
_RANGE = re.compile(r"^\d+-(?:\d+)?$")
_MAX_ROWS = 10_000


def _executable(name: str, blast_path: str | None) -> str | None:
    """``blast_path`` is the folder that holds the BLAST+ programs; None means PATH."""
    if blast_path is None:
        return None
    folder = Path(blast_path).expanduser()
    if not folder.is_dir():
        raise ConfigurationError(f"blast_path must be a folder that contains the BLAST+ programs: {folder}")
    return str(folder / name)


def _outfmt(value: int | str | None) -> list[str]:
    if value is None:
        return []
    text = str(value)
    if isinstance(value, bool) or not _OUTFMT.fullmatch(text) or not 0 <= int(text.split()[0]) <= 18:
        raise ConfigurationError("outfmt must be an integer 0-18 or a string such as '6 qseqid sseqid pident'.")
    return ["-outfmt", text]


def _database(value: str | Path, name: str, *, suffixes: Sequence[str]) -> str:
    prefix = Path(value).expanduser().resolve()
    if not any(next(prefix.parent.glob(f"{prefix.name}.{suffix}"), None) is not None for suffix in suffixes):
        raise ConfigurationError(f"{name} is not a BLAST database (no {prefix.name}.{'/'.join(suffixes)} files): {prefix}")
    return str(prefix)


def _common_search_options(
    evalue: float | None,
    max_target_seqs: int | None,
    num_threads: int | None,
    outfmt: int | str | None,
    word_size: int | None,
    max_hsps: int | None,
) -> list[str]:
    arguments: list[str] = []
    if evalue is not None:
        arguments += ["-evalue", str(_gap_cli.number(evalue, "evalue", minimum=0.0))]
    if word_size is not None:
        arguments += ["-word_size", str(_gap_cli.integer(word_size, "word_size", minimum=2))]
    if max_target_seqs is not None:
        arguments += ["-max_target_seqs", str(_gap_cli.integer(max_target_seqs, "max_target_seqs", minimum=1))]
    if max_hsps is not None:
        arguments += ["-max_hsps", str(_gap_cli.integer(max_hsps, "max_hsps", minimum=1))]
    if num_threads is not None:
        arguments += ["-num_threads", str(_gap_cli.integer(num_threads, "num_threads", minimum=1, maximum=256))]
    return arguments + _outfmt(outfmt)


def _tabular_rows(result: dict[str, Any], output: Path, outfmt: int | str | None) -> None:
    """Add parsed rows for tabular output (outfmt 6, 7 or 10) when the run succeeded."""
    if outfmt is None or str(outfmt).split()[0] not in {"6", "7", "10"} or result["status"] != "completed":
        return
    parts = str(outfmt).split()
    columns = list(_STANDARD_COLUMNS) if len(parts) == 1 or parts[1] == "std" else parts[1:]
    separator = "," if parts[0] == "10" else "\t"
    rows: list[dict[str, str]] = []
    with output.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(dict(zip(columns, line.rstrip("\n").split(separator), strict=False)))
            if len(rows) >= _MAX_ROWS:
                break
    result["rows"] = rows
    result["row_count"] = len(rows)


def _run(
    program: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    blast_path: str | None,
    stdout_file: Path | None = None,
) -> dict[str, Any]:
    return _p14_cli.run(
        program,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=_executable(program, blast_path),
        hint=_HINT,
        version_args=_VERSION_ARGS,
        stdout_file=stdout_file,
    )


def blastn_search(
    query: str,
    output: str,
    *,
    database: str | None = None,
    subject: str | None = None,
    task: str = "megablast",
    outfmt: int | str | None = None,
    evalue: float | None = None,
    word_size: int | None = None,
    max_target_seqs: int | None = None,
    max_hsps: int | None = None,
    perc_identity: float | None = None,
    strand: str | None = None,
    dust: bool | None = None,
    gapopen: int | None = None,
    gapextend: int | None = None,
    penalty: int | None = None,
    reward: int | None = None,
    template_type: str | None = None,
    template_length: int | None = None,
    num_threads: int | None = None,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide-nucleotide BLAST (``blastn``) with any of its tasks, against a BLAST database or a subject FASTA.

    ``task`` selects the algorithm: ``'blastn'`` (classic blastn, word size 11), ``'blastn-short'`` (queries under 50
    bases, for example primers or probes), ``'dc-megablast'`` (discontiguous megablast, cross-species; set
    ``template_type`` 'coding', 'optimal' or 'coding_and_optimal' and ``template_length`` 16, 18 or 21),
    ``'megablast'`` (the BLAST+ default) and ``'rmblastn'`` (repeat-element scoring; needs a BLAST+ build that includes
    RMBlast). Give exactly one of ``database`` (the prefix of a database made by makeblastdb) or ``subject`` (a FASTA
    file). ``outfmt`` is the BLAST+ output format (0-18 or a string such as ``'6 qseqid sseqid pident evalue'``);
    formats 6, 7 and 10 are also parsed into ``rows``. ``dust`` True/False turns query low-complexity filtering on or
    off. Options left as None keep BLAST+'s defaults. Writes a new file.
    """
    query_path = _gap_cli.existing_file(query, "query")
    output_path = _gap_cli.new_file(output, "output")
    if (database is None) == (subject is None):
        raise ConfigurationError("Give exactly one of database and subject.")
    arguments = ["-query", str(query_path), "-task", _gap_cli.choice(task, _NUCLEOTIDE_TASKS, "task")]
    if database is not None:
        arguments += ["-db", _database(database, "database", suffixes=("nin", "ndb", "nal"))]
    else:
        assert subject is not None
        arguments += ["-subject", str(_gap_cli.existing_file(subject, "subject"))]
    arguments += ["-out", str(output_path)]
    arguments += _common_search_options(evalue, max_target_seqs, num_threads, outfmt, word_size, max_hsps)
    if perc_identity is not None:
        arguments += ["-perc_identity", str(_gap_cli.number(perc_identity, "perc_identity", minimum=0.0, maximum=100.0))]
    if strand is not None:
        arguments += ["-strand", _gap_cli.choice(strand, _STRANDS, "strand")]
    if dust is not None:
        arguments += ["-dust", "yes" if _p14_cli.flag(dust, "dust") else "no"]
    for flag_name, value in (("gapopen", gapopen), ("gapextend", gapextend)):
        if value is not None:
            arguments += [f"-{flag_name}", str(_gap_cli.integer(value, flag_name, minimum=0))]
    if penalty is not None:
        arguments += ["-penalty", str(_gap_cli.integer(penalty, "penalty", maximum=-1))]
    if reward is not None:
        arguments += ["-reward", str(_gap_cli.integer(reward, "reward", minimum=1))]
    if template_type is not None:
        arguments += ["-template_type", _gap_cli.choice(template_type, ("coding", "coding_and_optimal", "optimal"), "template_type")]
    if template_length is not None:
        if template_length not in (16, 18, 21):
            raise ConfigurationError("template_length must be 16, 18 or 21.")
        arguments += ["-template_length", str(template_length)]
    result = _run("blastn", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)
    result["task"] = task
    _tabular_rows(result, output_path, outfmt)
    return result


def rpstblastn_search(
    query: str,
    database: str,
    output: str,
    *,
    outfmt: int | str | None = None,
    evalue: float | None = None,
    max_target_seqs: int | None = None,
    query_gencode: int | None = None,
    strand: str | None = None,
    num_threads: int | None = None,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a DNA query, translated in six frames, against a conserved-domain PSSM database (``rpstblastn``).

    ``database`` is the prefix of an RPS-BLAST database such as NCBI's ``Cdd`` (made with makeprofiledb);
    ``query_gencode`` is the NCBI genetic-code number used to translate the query (BLAST+ default 1); ``strand`` is
    'both', 'plus' or 'minus'. Output format and parsed rows are as in blastn_search. Writes a new file.
    """
    query_path = _gap_cli.existing_file(query, "query")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-query", str(query_path),
        "-db", _database(database, "database", suffixes=("rps", "pin", "pal")),
        "-out", str(output_path),
    ]
    arguments += _common_search_options(evalue, max_target_seqs, num_threads, outfmt, None, None)
    if query_gencode is not None:
        arguments += ["-query_gencode", str(_gap_cli.integer(query_gencode, "query_gencode", minimum=1, maximum=33))]
    if strand is not None:
        arguments += ["-strand", _gap_cli.choice(strand, _STRANDS, "strand")]
    result = _run("rpstblastn", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)
    _tabular_rows(result, output_path, outfmt)
    return result


def blastdbcmd_fetch_sequences(
    database: str,
    output: str,
    *,
    entries: Sequence[str] | None = None,
    entry_batch: str | None = None,
    database_type: str = "nucl",
    sequence_range: str | None = None,
    strand: str | None = None,
    outfmt: str | None = None,
    target_only: bool = False,
    line_length: int | None = None,
    long_seqids: bool = False,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve sequences (or sub-ranges) from a BLAST database by identifier (``blastdbcmd``).

    Give ``entries`` (identifiers, or ``['all']`` for every sequence; ``sequence_range`` such as ``'10-50'`` or ``'10-'`` and
    ``strand`` 'plus' or 'minus' then apply to the single entry) or ``entry_batch`` (a file with one identifier per
    line, optionally followed by a range, strand or mask id). ``outfmt`` is the blastdbcmd format string (default
    ``%f``, FASTA; for example ``'%a %l'`` lists accession and length). ``database_type`` is 'nucl' or 'prot'.
    Writes a new file.
    """
    db = _database(database, "database", suffixes=("nin", "nsq", "ndb", "pin", "psq", "pdb", "nal", "pal"))
    output_path = _gap_cli.new_file(output, "output")
    if (entries is None) == (entry_batch is None):
        raise ConfigurationError("Give exactly one of entries and entry_batch.")
    arguments = ["-db", db, "-dbtype", _gap_cli.choice(database_type, ("nucl", "prot"), "database_type")]
    if entries is not None:
        if isinstance(entries, str) or not entries:
            raise ConfigurationError("entries must be a non-empty list of identifiers.")
        names = [_gap_cli.text(item, "entries") for item in entries]
        if any("," in item for item in names):
            raise ConfigurationError("An entry must not contain a comma.")
        arguments += ["-entry", ",".join(names)]
        if sequence_range is not None:
            if not _RANGE.fullmatch(sequence_range):
                raise ConfigurationError("sequence_range must look like '10-50' or '10-'.")
            arguments += ["-range", sequence_range]
        if strand is not None:
            arguments += ["-strand", _gap_cli.choice(strand, ("plus", "minus"), "strand")]
    else:
        assert entry_batch is not None
        if sequence_range is not None or strand is not None:
            raise ConfigurationError("sequence_range and strand cannot be combined with entry_batch; put them in the file.")
        arguments += ["-entry_batch", str(_gap_cli.existing_file(entry_batch, "entry_batch"))]
    if outfmt is not None:
        arguments += ["-outfmt", _p14_cli.line_text(outfmt, "outfmt")]
    if _p14_cli.flag(target_only, "target_only"):
        arguments.append("-target_only")
    if line_length is not None:
        arguments += ["-line_length", str(_gap_cli.integer(line_length, "line_length", minimum=1))]
    if _p14_cli.flag(long_seqids, "long_seqids"):
        arguments.append("-long_seqids")
    arguments += ["-out", str(output_path)]
    return _run("blastdbcmd", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


def makembindex_build_index(
    input_path: str,
    *,
    input_format: str = "blastdb",
    output_prefix: str | None = None,
    legacy: bool | None = None,
    idmap: bool = False,
    nmer: int | None = None,
    ws_hint: int | None = None,
    volsize: int | None = None,
    stride: int | None = None,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a MegaBLAST index of a nucleotide BLAST database or FASTA file (``makembindex``).

    ``input_format='blastdb'`` (the default) indexes the database whose prefix is ``input_path`` in the new index
    format; BLAST+ writes the index volumes ``<prefix>.NN.idx`` next to the database and ignores ``output_prefix``.
    ``input_format='fasta'`` indexes a FASTA file in the old (deprecated) index format and needs ``output_prefix``,
    a new prefix for the volumes ``<output_prefix>.NN.idx``. ``legacy`` selects the legacy offset-list layout
    (BLAST+ default True); ``nmer`` is the indexed word length (8-15), ``ws_hint`` the likely search word size,
    ``volsize`` the volume size in MB, ``stride`` the distance between stored positions. Search with the index
    through blastn's ``-use_index`` option.
    """
    fmt = _gap_cli.choice(input_format, ("blastdb", "fasta"), "input_format")
    if fmt == "blastdb":
        if output_prefix is not None:
            raise ConfigurationError("output_prefix is ignored by makembindex for blastdb input; leave it None.")
        source = _database(input_path, "input_path", suffixes=("nin", "ndb", "nal"))
        stem = Path(source)
        arguments = ["-input", source, "-iformat", "blastdb"]
    else:
        if output_prefix is None:
            raise ConfigurationError("output_prefix is required for fasta input.")
        source = str(_gap_cli.existing_file(input_path, "input_path"))
        stem = Path(output_prefix).expanduser().resolve()
        if not stem.parent.is_dir():
            raise ConfigurationError(f"The folder for output_prefix does not exist: {stem.parent}")
        if any(stem.parent.glob(f"{stem.name}*")):
            raise ConfigurationError(f"output_prefix must not exist yet: {stem}")
        arguments = ["-input", source, "-iformat", "fasta", "-output", str(stem), "-old_style_index", "true"]
    if any(stem.parent.glob(f"{stem.name}.*.idx")) and fmt == "blastdb":
        raise ConfigurationError(f"An index for this database already exists: {stem}.NN.idx")
    if legacy is not None:
        arguments += ["-legacy", "true" if _p14_cli.flag(legacy, "legacy") else "false"]
    if _p14_cli.flag(idmap, "idmap"):
        arguments += ["-idmap", "true"]
    for name, value, low, high in (
        ("nmer", nmer, 8, 15), ("ws_hint", ws_hint, 1, None), ("volsize", volsize, 1, None), ("stride", stride, 1, None)
    ):
        if value is not None:
            arguments += [f"-{name}", str(_gap_cli.integer(value, name, minimum=low, maximum=high))]
    result = _run("makembindex", arguments, cwd=stem.parent, expected=[], blast_path=blast_path)
    produced = sorted(str(path) for path in stem.parent.glob(f"{stem.name}.*.idx") if path.is_file())
    result["outputs"] = produced
    if result["status"] == "completed" and not produced:
        result["status"] = "failed"
    return result


def windowmasker_mask_sequences(
    input_path: str,
    unit_counts: str,
    output: str,
    *,
    input_format: str = "fasta",
    output_format: str = "interval",
    window: int | None = None,
    t_thres: int | None = None,
    t_extend: int | None = None,
    t_low: int | None = None,
    t_high: int | None = None,
    set_t_low: int | None = None,
    set_t_high: int | None = None,
    dust: bool = False,
    dust_level: int | None = None,
    parse_seqids: bool = False,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask repeats in genome sequences with WindowMasker, using a unit-count file from windowmasker_count_units.

    ``output_format`` is 'interval' (the BLAST+ default: ``>id`` lines followed by ``start - end`` intervals), 'fasta'
    (soft-masked, lower case) or one of the maskinfo/seqloc ASN.1 and XML formats; ``input_format`` is 'fasta',
    'blastdb' or 'seqids'. ``window`` is the window size, ``t_thres`` the score that triggers masking, ``t_extend``
    the score above which masking may extend, ``t_low`` and ``t_high`` the lowest and highest useful unit score,
    ``set_t_low`` and ``set_t_high`` the replacement scores; ``dust`` combines window masking with DUST. Options left
    as None keep the values stored in the count file. Writes a new file.
    """
    source = _gap_cli.existing_file(input_path, "input_path")
    counts = _gap_cli.existing_file(unit_counts, "unit_counts")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-ustat", str(counts), "-in", str(source),
        "-infmt", _gap_cli.choice(input_format, ("fasta", "blastdb", "seqids"), "input_format"),
        "-outfmt", _gap_cli.choice(
            output_format,
            ("interval", "fasta", "maskinfo_asn1_bin", "maskinfo_asn1_text", "maskinfo_xml", "seqloc_asn1_bin", "seqloc_asn1_text", "seqloc_xml"),
            "output_format",
        ),
        "-out", str(output_path),
    ]
    for name, value in (
        ("window", window), ("t_thres", t_thres), ("t_extend", t_extend), ("t_low", t_low), ("t_high", t_high),
        ("set_t_low", set_t_low), ("set_t_high", set_t_high),
    ):
        if value is not None:
            arguments += [f"-{name}", str(_gap_cli.integer(value, name, minimum=1 if name != "t_extend" else 0))]
    if _p14_cli.flag(dust, "dust"):
        arguments += ["-dust", "true"]
    if dust_level is not None:
        arguments += ["-dust_level", str(_gap_cli.integer(dust_level, "dust_level", minimum=0))]
    if _p14_cli.flag(parse_seqids, "parse_seqids"):
        arguments.append("-parse_seqids")
    return _run("windowmasker", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


def windowmasker_count_units(
    input_path: str,
    output: str,
    *,
    unit_length: int | None = None,
    genome_size: int | None = None,
    input_format: str = "fasta",
    input_is_list: bool = False,
    check_duplicates: bool = False,
    memory_mb: int | None = None,
    counts_format: str = "ascii",
    counts_size_mb: int | None = None,
    metadata: str | None = None,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count k-mer units of genome sequences for WindowMasker (``windowmasker -mk_counts``).

    ``unit_length`` is the unit (k-mer) length 1-16 (WindowMasker chooses it from the genome size when None),
    ``genome_size`` the total genome length when it is known in advance, ``input_is_list`` says the input file lists
    FASTA file names, ``memory_mb`` bounds the memory for counting, ``counts_format`` is 'ascii', 'binary', 'oascii'
    or 'obinary' (the optimised formats are smaller and faster to load) and ``metadata`` is free text stored in the
    count file. The result feeds windowmasker_mask_sequences. Writes a new file.
    """
    source = _gap_cli.existing_file(input_path, "input_path")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-mk_counts", "-in", str(source), "-out", str(output_path),
        "-infmt", _gap_cli.choice(input_format, ("fasta", "blastdb", "seqids"), "input_format"),
        "-sformat", _gap_cli.choice(counts_format, ("ascii", "binary", "oascii", "obinary"), "counts_format"),
    ]
    if unit_length is not None:
        arguments += ["-unit", str(_gap_cli.integer(unit_length, "unit_length", minimum=1, maximum=16))]
    if genome_size is not None:
        arguments += ["-genome_size", str(_gap_cli.integer(genome_size, "genome_size", minimum=1))]
    if _p14_cli.flag(input_is_list, "input_is_list"):
        arguments += ["-fa_list", "true"]
    if _p14_cli.flag(check_duplicates, "check_duplicates"):
        arguments += ["-checkdup", "true"]
    if memory_mb is not None:
        arguments += ["-mem", str(_gap_cli.integer(memory_mb, "memory_mb", minimum=1))]
    if counts_size_mb is not None:
        arguments += ["-smem", str(_gap_cli.integer(counts_size_mb, "counts_size_mb", minimum=1))]
    if metadata is not None:
        arguments += ["-meta", _p14_cli.line_text(metadata, "metadata")]
    return _run("windowmasker", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


def windowmasker_convert_counts(
    input_path: str,
    output: str,
    *,
    counts_format: str = "ascii",
    counts_size_mb: int | None = None,
    metadata: str | None = None,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a WindowMasker unit-count file to another count format (``windowmasker -convert``).

    ``counts_format`` is the target format: 'ascii', 'binary', 'oascii' or 'obinary'. Writes a new file.
    """
    source = _gap_cli.existing_file(input_path, "input_path")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-convert", "-in", str(source), "-out", str(output_path),
        "-sformat", _gap_cli.choice(counts_format, ("ascii", "binary", "oascii", "obinary"), "counts_format"),
    ]
    if counts_size_mb is not None:
        arguments += ["-smem", str(_gap_cli.integer(counts_size_mb, "counts_size_mb", minimum=1))]
    if metadata is not None:
        arguments += ["-meta", _p14_cli.line_text(metadata, "metadata")]
    return _run("windowmasker", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


def convert2blastmask_from_fasta(
    input_path: str,
    output: str,
    *,
    masking_algorithm: str,
    masking_options: str,
    output_format: str = "maskinfo_asn1_text",
    parse_seqids: bool = False,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert lower-case masking in a FASTA file into a mask file for makeblastdb (``convert2blastmask``).

    ``masking_algorithm`` names the program that produced the masking ('dust', 'seg', 'windowmasker', 'repeat' or
    'other'); ``masking_options`` is free text that records its options. ``output_format`` is 'maskinfo_asn1_text'
    (BLAST+ default), 'maskinfo_asn1_bin', 'maskinfo_xml' or 'interval'. Use the output with
    ``makeblastdb -mask_data``. Writes a new file.
    """
    source = _gap_cli.existing_file(input_path, "input_path")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-in", str(source), "-out", str(output_path),
        "-outfmt", _gap_cli.choice(output_format, ("interval", "maskinfo_asn1_bin", "maskinfo_asn1_text", "maskinfo_xml"), "output_format"),
        "-masking_algorithm", _gap_cli.text(masking_algorithm, "masking_algorithm"),
        "-masking_options", _p14_cli.line_text(masking_options, "masking_options"),
    ]
    if _p14_cli.flag(parse_seqids, "parse_seqids"):
        arguments.append("-parse_seqids")
    return _run("convert2blastmask", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


def blast_vdb_search(
    query: str,
    accessions: Sequence[str],
    output: str,
    *,
    program: str = "blastn_vdb",
    outfmt: int | str | None = None,
    evalue: float | None = None,
    max_target_seqs: int | None = None,
    num_threads: int | None = None,
    task: str | None = None,
    genetic_code: int | None = None,
    sra_mode: str | None = None,
    include_filtered_reads: bool = False,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search SRA run accessions with a nucleotide query (``blastn_vdb``) or a protein query (``tblastn_vdb``).

    ``program`` is 'blastn_vdb' (nucleotide query against the reads of the runs; ``task`` as in blastn_search) or
    'tblastn_vdb' (protein query against the six-frame translation of the reads; ``genetic_code`` is the NCBI
    database genetic code). ``accessions`` are SRA run accessions such as ``['SRR000001']``; the programs read them
    from NCBI over the network. ``sra_mode`` is the SRA search mode and ``include_filtered_reads`` also searches
    reads that SRA flags as filtered. Output format and parsed rows are as in blastn_search. Writes a new file.
    """
    chosen = _gap_cli.choice(program, ("blastn_vdb", "tblastn_vdb"), "program")
    query_path = _gap_cli.existing_file(query, "query")
    output_path = _gap_cli.new_file(output, "output")
    if isinstance(accessions, str) or not accessions:
        raise ConfigurationError("accessions must be a non-empty list of SRA run accessions.")
    names = [_gap_cli.text(item, "accessions") for item in accessions]
    if any(re.search(r"\s", item) for item in names):
        raise ConfigurationError("An accession must not contain white space.")
    arguments = ["-query", str(query_path), "-db", " ".join(names), "-out", str(output_path)]
    arguments += _common_search_options(evalue, max_target_seqs, num_threads, outfmt, None, None)
    if task is not None:
        if chosen != "blastn_vdb":
            raise ConfigurationError("task applies to blastn_vdb only.")
        arguments += ["-task", _gap_cli.choice(task, _NUCLEOTIDE_TASKS, "task")]
    if genetic_code is not None:
        if chosen != "tblastn_vdb":
            raise ConfigurationError("genetic_code applies to tblastn_vdb only.")
        arguments += ["-db_gencode", str(_gap_cli.integer(genetic_code, "genetic_code", minimum=1, maximum=33))]
    if sra_mode is not None:
        arguments += ["-sra_mode", _gap_cli.text(sra_mode, "sra_mode")]
    if _p14_cli.flag(include_filtered_reads, "include_filtered_reads"):
        arguments.append("-include_filtered_reads")
    result = _run(chosen, arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)
    _tabular_rows(result, output_path, outfmt)
    return result


def blast_vdb_cmd_fetch_reads(
    accessions: Sequence[str],
    output: str,
    *,
    view: str = "sequences",
    entries: Sequence[str] | None = None,
    entry_batch: str | None = None,
    outfmt: str | None = None,
    line_length: int | None = None,
    include_filtered_reads: bool = False,
    blast_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve sequence records from SRA run accessions (``blast_vdb_cmd``); the runs are read from NCBI over the network.

    ``view='sequences'`` writes the reads named by ``entries`` (identifiers such as ``'gnl|SRA|SRR000001.1.2'``, or
    ``['all']``) or listed in ``entry_batch`` (one identifier per line); ``outfmt`` is the format string (default
    ``%f``, FASTA; ``'%a %l'`` lists accession and length) and ``line_length`` the FASTA line width (default 80).
    ``include_filtered_reads`` also returns reads that SRA flags as filtered. ``view='info'`` prints the number of
    sequences and bases of each run, ``'references'`` dumps the reference sequences of aligned runs, and
    ``'paths'`` / ``'paths_all'`` print the local top-level or all paths of the run files. ``accessions`` are SRA run
    accessions such as ``['SRR000001']``. Writes a new file.
    """
    chosen = _gap_cli.choice(view, ("sequences", "info", "references", "paths", "paths_all"), "view")
    output_path = _gap_cli.new_file(output, "output")
    if isinstance(accessions, str) or not accessions:
        raise ConfigurationError("accessions must be a non-empty list of SRA run accessions.")
    names = [_gap_cli.text(item, "accessions") for item in accessions]
    if any(re.search(r"\s", item) for item in names):
        raise ConfigurationError("An accession must not contain white space.")
    arguments = ["-db", " ".join(names), "-out", str(output_path)]
    if chosen == "sequences":
        if (entries is None) == (entry_batch is None):
            raise ConfigurationError("Give exactly one of entries and entry_batch with view='sequences'.")
        if entries is not None:
            if isinstance(entries, str) or not entries:
                raise ConfigurationError("entries must be a non-empty list of identifiers.")
            identifiers = [_gap_cli.text(item, "entries") for item in entries]
            if any("," in item for item in identifiers):
                raise ConfigurationError("An entry must not contain a comma.")
            arguments += ["-entry", ",".join(identifiers)]
        else:
            assert entry_batch is not None
            arguments += ["-entry_batch", str(_gap_cli.existing_file(entry_batch, "entry_batch"))]
        if outfmt is not None:
            arguments += ["-outfmt", _p14_cli.line_text(outfmt, "outfmt")]
        if line_length is not None:
            arguments += ["-line_length", str(_gap_cli.integer(line_length, "line_length", minimum=1))]
        if _p14_cli.flag(include_filtered_reads, "include_filtered_reads"):
            arguments.append("-include_filtered_reads")
    else:
        if entries is not None or entry_batch is not None or outfmt is not None or line_length is not None or include_filtered_reads:
            raise ConfigurationError("entries, entry_batch, outfmt, line_length and include_filtered_reads apply to view='sequences' only.")
        arguments.append({"info": "-info", "references": "-ref", "paths": "-paths", "paths_all": "-paths_all"}[chosen])
    return _run("blast_vdb_cmd", arguments, cwd=output_path.parent, expected=[output_path], blast_path=blast_path)


__all__ = [
    "blast_vdb_cmd_fetch_reads",
    "blast_vdb_search",
    "blastdbcmd_fetch_sequences",
    "blastn_search",
    "convert2blastmask_from_fasta",
    "makembindex_build_index",
    "rpstblastn_search",
    "windowmasker_convert_counts",
    "windowmasker_count_units",
    "windowmasker_mask_sequences",
]
_register_public_verbose(globals())
