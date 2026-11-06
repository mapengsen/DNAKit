"""dendropy operations as DNAKit functions: sequence and tree format conversion, sumtrees and sumlabels summaries,
NCBI Entrez record retrieval, MUSCLE alignment and PAUP* tree estimation.

Install with ``python -m pip install "dnakit[bio-dendropy]"``. Format conversion, Entrez retrieval and the PAUP*
tree operations call the dendropy library itself. The dendropy-format, sumtrees and sumlabels functions run the
console applications with a fixed argument list (no shell). MUSCLE and PAUP* are external programs: they are
looked up before anything runs, and a missing program raises ``BackendUnavailableError`` with an install hint.
Outputs are always new files; an existing output path is refused, never overwritten.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError, DNAKitError

_EXTRA = "bio-dendropy"
_DENDROPY_HINT = (
    'Install dendropy 5.1 or later: python -m pip install "dnakit[bio-dendropy]" '
    "(dendropy also provides the dendropy-format, sumtrees and sumlabels console scripts)."
)
_SUMTREES_HINT = "sumtrees comes with dendropy 5.1: python -m pip install 'dendropy>=5.1,<6' (or dnakit[bio-dendropy])."
_SUMLABELS_HINT = "sumlabels comes with dendropy 5.1: python -m pip install 'dendropy>=5.1,<6' (or dnakit[bio-dendropy])."
_MUSCLE_HINT = (
    "Install MUSCLE (drive5.org/muscle5, or conda install -c bioconda muscle) and put it on PATH or pass muscle_path. "
    "The default call uses the stdin-to-stdout protocol of dendropy's muscle_align (MUSCLE 3.8 style); "
    "pass muscle_args for another command line."
)
_PAUP_HINT = (
    "Install PAUP* 4.0 (paup.sc.fsu.edu, or conda install -c bioconda paup), name the executable paup, "
    "put it on PATH or pass paup_path."
)
_MAX_STDOUT_BYTES = 100_000_000
_LOG_TAIL_CHARS = 4000
_ENTREZ_MAX_IDS = 200
_SHELL_SAFE_PATH = re.compile(r"[A-Za-z0-9_./+:@%-]+")
_PAUP_TOKEN = re.compile(r"[A-Za-z0-9_.=:+-]+")


def _dendropy() -> Any:
    return require("dendropy", _EXTRA)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _free_text(value: object, name: str) -> str:
    """Non-empty single-line text; it is always passed in the --name=value form, so it may start with '-'."""
    if not isinstance(value, str) or not value.strip() or "\x00" in value or "\n" in value:
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def _existing_files(values: Sequence[str], name: str) -> list[Path]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or len(values) == 0:
        raise ConfigurationError(f"{name} must be a non-empty list of existing files.")
    return [_gap_cli.existing_file(item, name) for item in values]


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


# ---------------------------------------------------------------------------
# Sequence and tree format conversion (dendropy.application.dendropy_format)
# ---------------------------------------------------------------------------

_INPUT_FORMATS = (
    "fasta",
    "newick",
    "nexml",
    "nexus",
    "phylip-strict",
    "phylip-relaxed-singlespace",
    "phylip-relaxed-multispace",
    "phylip-strict-interleaved",
    "phylip-relaxed-singlespace-interleaved",
    "phylip-relaxed-multispace-interleaved",
)
_OUTPUT_FORMATS = ("fasta", "newick", "nexml", "nexus", "phylip", "phylip-strict")
_DATA_TYPES = ("dna", "rna", "standard")
_RECODE_UNCERTAIN = ("missing", "gap")
_PHYLIP_INPUT_OPTIONS: Mapping[str, Mapping[str, bool]] = {
    "phylip-strict": {"strict": True, "multispace_delimiter": False, "interleaved": False},
    "phylip-relaxed-singlespace": {
        "strict": False,
        "multispace_delimiter": False,
        "interleaved": False,
    },
    "phylip-relaxed-multispace": {
        "strict": False,
        "multispace_delimiter": True,
        "interleaved": False,
    },
    "phylip-strict-interleaved": {
        "strict": True,
        "multispace_delimiter": False,
        "interleaved": True,
    },
    "phylip-relaxed-singlespace-interleaved": {
        "strict": False,
        "multispace_delimiter": False,
        "interleaved": True,
    },
    "phylip-relaxed-multispace-interleaved": {
        "strict": False,
        "multispace_delimiter": True,
        "interleaved": True,
    },
}


def _read_plan(input_format: str, data_type: str | None) -> tuple[str, dict[str, Any]]:
    """Reader schema and keyword arguments, built the way dendropy-format builds them for each --from value."""
    options: dict[str, Any] = {}
    schema = input_format
    if input_format.startswith("phylip"):
        schema = "phylip"
        options.update(_PHYLIP_INPUT_OPTIONS[input_format])
    if schema in ("fasta", "phylip"):
        if data_type is None:
            raise ConfigurationError(
                "data_type (dna, rna or standard) is required when reading FASTA or PHYLIP input."
            )
        options["data_type"] = data_type
    return schema, options


def _write_plan(output_format: str, unquoted_underscores: bool) -> tuple[str, dict[str, Any]]:
    options: dict[str, Any] = {}
    schema = output_format
    if output_format == "phylip-strict":
        schema = "phylip"
        options["strict"] = True
    if schema in ("nexus", "newick") and unquoted_underscores:
        options["unquoted_underscores"] = True
    return schema, options


def _recode_uncertain(dataset: Any, dendropy: Any, mode: str) -> None:
    """Replace ambiguous and polymorphic DNA states by gaps ('-') or missing data ('?').

    dendropy-format's own 'missing' branch reads alphabet.missing, which dendropy 5.1 does not define (AttributeError);
    the missing-data state of the alphabet (alphabet.no_data_state, symbol '?') is what the flag documents.
    """
    alphabet = dendropy.DNA_STATE_ALPHABET
    replacement = alphabet.gap if mode == "gap" else alphabet.no_data_state
    uncertain = set(alphabet._polymorphic_states) | set(alphabet._ambiguous_states)
    for char_matrix in dataset.char_matrices:
        for taxon in char_matrix:
            sequence = char_matrix[taxon]
            for index, state in enumerate(sequence):
                if state in uncertain:
                    sequence[index] = replacement


def convert_format(
    source: str,
    output: str,
    *,
    input_format: str,
    output_format: str | None = None,
    data_type: str | None = None,
    unquoted_underscores: bool = False,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a DNA sequence file or a tree file to another format with dendropy's reader and writer.

    The steps are those of dendropy-format's convert command: the phylip-* input names select the PHYLIP reader
    options; FASTA and PHYLIP input need data_type (dna, rna or standard); the output format defaults to the
    input format; unquoted_underscores applies to NEXUS and Newick output; recode_uncertain ('missing' or 'gap')
    replaces ambiguous and polymorphic DNA codes. This does not call the console application's convert()
    function, which closes sys.stdout when it returns and so breaks a running Python session. Writes a new file.
    """
    dendropy = _dendropy()
    source_path = _gap_cli.existing_file(source, "source")
    output_path = _gap_cli.new_file(output, "output")
    input_name = _gap_cli.choice(input_format, _INPUT_FORMATS, "input_format")
    checked_data_type = (
        None if data_type is None else _gap_cli.choice(data_type, _DATA_TYPES, "data_type")
    )
    read_schema, read_options = _read_plan(input_name, checked_data_type)
    target = (
        read_schema
        if output_format is None
        else _gap_cli.choice(output_format, _OUTPUT_FORMATS, "output_format")
    )
    write_schema, write_options = _write_plan(
        target, _flag(unquoted_underscores, "unquoted_underscores")
    )
    mode = (
        None
        if recode_uncertain is None
        else _gap_cli.choice(recode_uncertain, _RECODE_UNCERTAIN, "recode_uncertain")
    )
    with source_path.open() as handle:
        dataset = dendropy.DataSet.get(file=handle, schema=read_schema, **read_options)
    if mode is not None:
        _recode_uncertain(dataset, dendropy, mode)
    with output_path.open("x") as handle:
        dataset.write(file=handle, schema=write_schema, **write_options)
    return {
        "status": "completed",
        "backend": "dendropy",
        "version": dendropy.__version__,
        "operation": "convert",
        "input_format": input_name,
        "output_format": target,
        "data_type": checked_data_type,
        "tree_lists": len(dataset.tree_lists),
        "character_matrices": len(dataset.char_matrices),
        "outputs": [str(output_path)],
    }


def format_to_nexus(
    source: str,
    output: str,
    *,
    input_format: str,
    data_type: str | None = None,
    unquoted_underscores: bool = False,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence or tree file to NEXUS (dendropy-format --to nexus); see convert_format."""
    return convert_format(
        source,
        output,
        input_format=input_format,
        output_format="nexus",
        data_type=data_type,
        unquoted_underscores=unquoted_underscores,
        recode_uncertain=recode_uncertain,
    )


def format_to_newick(
    source: str,
    output: str,
    *,
    input_format: str,
    data_type: str | None = None,
    unquoted_underscores: bool = False,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence or tree file to Newick (dendropy-format --to newick); see convert_format."""
    return convert_format(
        source,
        output,
        input_format=input_format,
        output_format="newick",
        data_type=data_type,
        unquoted_underscores=unquoted_underscores,
        recode_uncertain=recode_uncertain,
    )


def format_to_nexml(
    source: str,
    output: str,
    *,
    input_format: str,
    data_type: str | None = None,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence or tree file to NeXML (dendropy-format --to nexml); see convert_format."""
    return convert_format(
        source,
        output,
        input_format=input_format,
        output_format="nexml",
        data_type=data_type,
        recode_uncertain=recode_uncertain,
    )


def format_to_fasta(
    source: str,
    output: str,
    *,
    input_format: str,
    data_type: str | None = None,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence or tree file to FASTA (dendropy-format --to fasta); see convert_format."""
    return convert_format(
        source,
        output,
        input_format=input_format,
        output_format="fasta",
        data_type=data_type,
        recode_uncertain=recode_uncertain,
    )


def format_to_phylip(
    source: str,
    output: str,
    *,
    input_format: str,
    data_type: str | None = None,
    recode_uncertain: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence file to relaxed PHYLIP (dendropy-format --to phylip); see convert_format."""
    return convert_format(
        source,
        output,
        input_format=input_format,
        output_format="phylip",
        data_type=data_type,
        recode_uncertain=recode_uncertain,
    )


def _run_stdout_program(
    name: str,
    arguments: Sequence[str],
    *,
    output: Path,
    explicit: str | None,
    hint: str,
    version: str,
) -> dict[str, Any]:
    """Run a fixed program and save its standard output as a new file (the bounded runner merges stderr into it)."""
    executable = _gap_cli.program(name, explicit=explicit, hint=hint)
    result = execute_bounded_command(
        executable,
        list(arguments),
        backend_id=name,
        cwd=output.parent,
        timeout_seconds=_gap_cli.DEFAULT_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_STDOUT_BYTES,
        verbose=False,
    )
    complete = result.return_code == 0
    if complete:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(result.output)
    return {
        "status": "completed" if complete else "failed",
        "backend": name,
        "version": version,
        "command": [name, *arguments],
        "return_code": result.return_code,
        "outputs": [str(output)] if output.is_file() else [],
        "log_tail": "" if complete else _tail(result.output),
    }


def run_dendropy_format(
    source: str,
    output: str,
    *,
    input_format: str,
    output_format: str | None = None,
    data_type: str | None = None,
    unquoted_underscores: bool = False,
    recode_uncertain: str | None = None,
    dendropy_format_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the dendropy-format console application on a file; its standard output becomes the new file ``output``.

    The flags are the application's own: --from, --to, --data-type, --unquoted-underscores and
    --recode-uncertain. The same conversion without a subprocess is convert_format.
    """
    dendropy = _dendropy()
    source_path = _gap_cli.existing_file(source, "source")
    output_path = _gap_cli.new_file(output, "output")
    input_name = _gap_cli.choice(input_format, _INPUT_FORMATS, "input_format")
    checked_data_type = (
        None if data_type is None else _gap_cli.choice(data_type, _DATA_TYPES, "data_type")
    )
    if checked_data_type is None and (input_name == "fasta" or input_name.startswith("phylip")):
        raise ConfigurationError(
            "data_type (dna, rna or standard) is required when reading FASTA or PHYLIP input."
        )
    arguments = ["--from", input_name]
    if output_format is not None:
        arguments += ["--to", _gap_cli.choice(output_format, _OUTPUT_FORMATS, "output_format")]
    if checked_data_type is not None:
        arguments += ["--data-type", checked_data_type]
    if _flag(unquoted_underscores, "unquoted_underscores"):
        arguments.append("--unquoted-underscores")
    if recode_uncertain is not None:
        mode = _gap_cli.choice(recode_uncertain, _RECODE_UNCERTAIN, "recode_uncertain")
        if mode == "missing":
            raise ConfigurationError(
                "dendropy-format 5.1 fails with --recode-uncertain missing (AttributeError: 'DnaStateAlphabet' object "
                "has no attribute 'missing'). Use recode_uncertain='gap' here, or convert_format, which writes '?'."
            )
        arguments += ["--recode-uncertain", mode]
    arguments.append(str(source_path))
    return _run_stdout_program(
        "dendropy-format",
        arguments,
        output=output_path,
        explicit=dendropy_format_path,
        hint=_DENDROPY_HINT,
        version=dendropy.__version__,
    )


# ---------------------------------------------------------------------------
# Tree summaries: sumtrees (dendropy.application.sumtrees) and sumlabels (dendropy.application.sumlabels)
# ---------------------------------------------------------------------------

_SUMTREES_INPUT_FORMATS = ("nexus/newick", "nexus", "newick", "phylip", "nexml")
_SUMTREES_OUTPUT_FORMATS = ("nexus", "newick", "phylip", "nexml")
_SUMTREES_TARGETS = ("consensus", "mcct", "msct")
_SUMTREES_EDGES = (
    "mean-length",
    "median-length",
    "mean-age",
    "median-age",
    "support",
    "keep",
    "clear",
)
_SUMTREES_LABELS = ("support", "keep", "clear")


def summarize_tree_samples(
    tree_files: Sequence[str],
    output: str,
    *,
    input_format: str | None = None,
    burnin: int = 0,
    rooted: bool | None = None,
    weighted_trees: bool = False,
    summary_target: str | None = None,
    target_tree: str | None = None,
    min_clade_freq: float | None = None,
    edges: str | None = None,
    labels: str = "support",
    percentages: bool = False,
    decimals: int = 8,
    output_format: str | None = None,
    sumtrees_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise tree samples (for example MCMC output) with the sumtrees application.

    summary_target chooses the summary tree: 'consensus' (the default, as in sumtrees), 'mcct' (maximum clade
    credibility tree) or 'msct'. target_tree instead scores the clades of one tree against the samples (sumtrees
    -t). The two cannot be combined: sumtrees reports the conflict and then crashes, so this refuses it first.
    burnin drops that many trees from the start of each file; edges sets the branch lengths (mean-length,
    median-length, mean-age, median-age); labels sets the node labels (support, keep, clear); percentages and
    decimals format the support values. Runs the sumtrees program with a fixed argument list and writes the
    summary to a new file.
    """
    sources = _existing_files(tree_files, "tree_files")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["--quiet"]
    if input_format is not None:
        arguments += [
            "--input-format",
            _gap_cli.choice(input_format, _SUMTREES_INPUT_FORMATS, "input_format"),
        ]
    arguments += ["--burnin", str(_gap_cli.integer(burnin, "burnin", minimum=0))]
    if rooted is True:
        arguments.append("--force-rooted")
    elif rooted is False:
        arguments.append("--force-unrooted")
    if _flag(weighted_trees, "weighted_trees"):
        arguments.append("--weighted-trees")
    if target_tree is not None and summary_target is not None:
        raise ConfigurationError(
            "Give either summary_target or target_tree, not both (sumtrees crashes on that combination)."
        )
    if target_tree is not None:
        arguments += [
            "--target-tree-filepath",
            str(_gap_cli.existing_file(target_tree, "target_tree")),
        ]
    else:
        target = _gap_cli.choice(summary_target or "consensus", _SUMTREES_TARGETS, "summary_target")
        arguments += ["--summary-target", target]
    if min_clade_freq is not None:
        value = _gap_cli.number(min_clade_freq, "min_clade_freq", minimum=0.0, maximum=1.0)
        arguments += ["--min-clade-freq", repr(value)]
    if edges is not None:
        arguments += ["--edges", _gap_cli.choice(edges, _SUMTREES_EDGES, "edges")]
    arguments += ["--labels", _gap_cli.choice(labels, _SUMTREES_LABELS, "labels")]
    if _flag(percentages, "percentages"):
        arguments.append("--percentages")
    arguments += ["--decimals", str(_gap_cli.integer(decimals, "decimals", minimum=0))]
    if output_format is not None:
        arguments += [
            "--output-tree-format",
            _gap_cli.choice(output_format, _SUMTREES_OUTPUT_FORMATS, "output_format"),
        ]
    arguments += ["--output", str(output_path)]
    arguments += [str(path) for path in sources]
    return _gap_cli.run(
        "sumtrees",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        explicit=sumtrees_path,
        hint=_SUMTREES_HINT,
    )


def transfer_tree_labels(
    tree_files: Sequence[str],
    target_tree: str,
    output: str,
    *,
    separator: str = "/",
    preserve_target_labels: bool = False,
    rooted: bool | None = None,
    include_taxa_block: bool = True,
    additional_comments: str | None = None,
    to_newick: bool = False,
    sumlabels_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy the labels that source trees give to matching splits onto one target tree (the sumlabels application).

    Each split of the target tree receives the labels found on the same split in the source trees, joined by
    separator. preserve_target_labels keeps the target's own labels; include_taxa_block and to_newick control the
    output format (NEXUS with a taxa block by default). Runs the sumlabels program with a fixed argument list and
    writes a new file.
    """
    sources = _existing_files(tree_files, "tree_files")
    target = _gap_cli.existing_file(target_tree, "target_tree")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "--quiet",
        "--target",
        str(target),
        f"--separator={_free_text(separator, 'separator')}",
    ]
    if _flag(preserve_target_labels, "preserve_target_labels"):
        arguments.append("--preserve-target-labels")
    if rooted is True:
        arguments.append("--rooted")
    elif rooted is False:
        arguments.append("--unrooted")
    if not _flag(include_taxa_block, "include_taxa_block"):
        arguments.append("--no-taxa-block")
    if additional_comments is not None:
        arguments.append(
            f"--additional-comments={_free_text(additional_comments, 'additional_comments')}"
        )
    if _flag(to_newick, "to_newick"):
        arguments.append("--to-newick")
    arguments += ["--output", str(output_path)]
    arguments += [str(path) for path in sources]
    return _gap_cli.run(
        "sumlabels",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        explicit=sumlabels_path,
        hint=_SUMLABELS_HINT,
    )


# ---------------------------------------------------------------------------
# NCBI Entrez efetch (dendropy.interop.entrez)
# ---------------------------------------------------------------------------


def _entrez_identifiers(ids: object) -> list[str]:
    if (
        isinstance(ids, (str, bytes))
        or not isinstance(ids, Sequence)
        or not 0 < len(ids) <= _ENTREZ_MAX_IDS
    ):
        raise ConfigurationError(f"ids must be a list of 1 to {_ENTREZ_MAX_IDS} Entrez UIDs.")
    identifiers: list[str] = []
    for value in ids:
        if isinstance(value, bool) or not isinstance(value, (str, int)):
            raise ConfigurationError(
                "Each Entrez identifier must be an integer or a string of digits."
            )
        digits = str(value).strip()
        if not re.fullmatch(r"[0-9]+", digits):
            raise ConfigurationError(
                f"Entrez identifiers are numeric UIDs; got {digits!r}. Look up accessions with ESearch first."
            )
        identifiers.append(digits)
    return identifiers


def _record_count(records: str, rettype: str, retmode: str) -> int | None:
    if rettype == "fasta" and retmode == "text":
        return sum(1 for line in records.splitlines() if line.startswith(">"))
    return None


def fetch_entrez_records(
    database: str,
    ids: Sequence[str | int],
    *,
    rettype: str = "fasta",
    retmode: str = "text",
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch sequence records by numeric Entrez UID from an NCBI database (E-utilities efetch).

    database is the Entrez database (for example 'nuccore' for nucleotide records); rettype and retmode are the
    efetch return type and mode (rettype 'fasta' with retmode 'text' gives FASTA, and record_count counts its
    headers). The request is made by dendropy's efetch, which sends it over HTTP to NCBI. Returns the records as
    text together with the service name, the access date (UTC) and the final request URL. Nothing is cached.
    """
    _dendropy()
    from dendropy.interop import entrez as _entrez

    db = _gap_cli.text(database, "database")
    identifiers = _entrez_identifiers(ids)
    rettype_name = _gap_cli.text(rettype, "rettype")
    mode = _gap_cli.choice(retmode, ("text", "xml"), "retmode")
    mail = None
    if email is not None:
        mail = _free_text(email, "email")
        if "@" not in mail:
            raise ConfigurationError("email must be an address containing '@'.")
    access_date = datetime.now(timezone.utc).date().isoformat()
    response = _entrez.efetch(db, ",".join(identifiers), rettype_name, retmode=mode, email=mail)
    try:
        body = response.read()
        final_url = response.geturl() if hasattr(response, "geturl") else None
    finally:
        response.close()
    records = body.decode("utf-8") if isinstance(body, bytes) else str(body)
    return {
        "status": "completed",
        "backend": "NCBI Entrez E-utilities",
        "service": "efetch",
        "database": db,
        "rettype": rettype_name,
        "retmode": mode,
        "identifier_count": len(identifiers),
        "record_count": _record_count(records, rettype_name, mode),
        "access_date": access_date,
        "url": final_url,
        "records": records,
    }


# ---------------------------------------------------------------------------
# MUSCLE (dendropy.interop.muscle): external program, looked up before running
# ---------------------------------------------------------------------------


def align_with_muscle(
    fasta: str,
    output: str,
    *,
    muscle_args: Sequence[str] | None = None,
    muscle_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align DNA sequences with the MUSCLE program (dendropy's muscle_align) and save the alignment as FASTA.

    Without muscle_args the program is started with no arguments and the sequences go in on standard input,
    which is how dendropy calls it (MUSCLE 3.8 style). MUSCLE 5 needs its own arguments, for example
    ['-align', input, '-output', output]; pass them in muscle_args. A non-zero exit raises BackendExecutionError
    with MUSCLE's message. Writes a new file.
    """
    dendropy = _dendropy()
    from dendropy.interop import muscle as _muscle

    source = _gap_cli.existing_file(fasta, "fasta")
    output_path = _gap_cli.new_file(output, "output")
    arguments = _program_arguments(muscle_args, "muscle_args")
    executable = _gap_cli.program("muscle", explicit=muscle_path, hint=_MUSCLE_HINT)
    sequences = dendropy.DnaCharacterMatrix.get(path=str(source), schema="fasta")
    try:
        aligned = _muscle.muscle_align(
            sequences, muscle_args=arguments or None, muscle_path=executable
        )
    except Exception as exc:  # dendropy raises a bare Exception that carries MUSCLE's stderr
        raise BackendExecutionError(
            f"MUSCLE failed: {exc}", context={"backend_id": "muscle"}
        ) from exc
    with output_path.open("x") as handle:
        aligned.write_to_stream(handle, schema="fasta")
    return {
        "status": "completed",
        "backend": "muscle",
        "version": _gap_cli._version_of(executable),
        "command": [executable, *arguments],
        "sequence_count": len(aligned),
        "alignment_length": aligned.max_sequence_size,
        "outputs": [str(output_path)],
    }


def _program_arguments(values: Sequence[str] | None, name: str) -> list[str]:
    if values is None:
        return []
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of strings.")
    arguments: list[str] = []
    for value in values:
        if not isinstance(value, str) or not value or "\x00" in value:
            raise ConfigurationError(
                f"Each item of {name} must be a non-empty string without NUL characters."
            )
        arguments.append(value)
    return arguments


# ---------------------------------------------------------------------------
# PAUP* (dendropy.interop.paup): external program, looked up before running
# ---------------------------------------------------------------------------

_PAUP_CRITERIA = ("likelihood", "parsimony", "nj", "upgma")
_PAUP_SUBSTITUTION_TYPES = (1, 2, 6)


def _paup_executable(paup_path: str | None) -> str:
    executable = _gap_cli.program("paup", explicit=paup_path, hint=_PAUP_HINT)
    if not _SHELL_SAFE_PATH.fullmatch(executable):
        raise ConfigurationError(
            "The PAUP* path may contain only letters, digits and ._/+:@%- : dendropy starts PAUP* through a shell "
            f"in some of its functions. Got {executable!r}; put a link with a plain name on PATH."
        )
    return executable


def _paup_call(operation: str, call: Any) -> Any:
    try:
        return call()
    except DNAKitError:
        raise
    except (
        Exception
    ) as exc:  # dendropy reports PAUP* failures as ExternalServiceError or bare exceptions
        raise BackendExecutionError(
            f"PAUP* {operation} failed: {exc}", context={"backend_id": "paup"}
        ) from exc


def _paup_result(
    operation: str, executable: str, outputs: Sequence[Path], **extra: Any
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": "PAUP*",
        "operation": operation,
        "version": _gap_cli._version_of(executable),
        "outputs": [str(path) for path in outputs],
    }
    result.update(extra)
    return result


def _read_alignment(dendropy: Any, path: Path) -> Any:
    return dendropy.DnaCharacterMatrix.get(path=str(path), schema="fasta")


def _tree_schema(path: Path) -> str:
    with path.open() as handle:
        for line in handle:
            if line.strip():
                return "nexus" if line.strip().upper().startswith("#NEXUS") else "newick"
    return "newick"


def _read_trees(dendropy: Any, path: Path, namespace: Any) -> Any:
    return dendropy.TreeList.get(
        path=str(path), schema=_tree_schema(path), taxon_namespace=namespace
    )


def _write_tree(tree: Any, output: Path) -> None:
    with output.open("x") as handle:
        tree.write_to_stream(handle, schema="newick")


def _check_estimation(criterion: str, substitution_types: object) -> tuple[str, int]:
    name = _gap_cli.choice(criterion, _PAUP_CRITERIA, "criterion")
    if isinstance(substitution_types, bool) or substitution_types not in _PAUP_SUBSTITUTION_TYPES:
        raise ConfigurationError("substitution_types must be 1, 2 or 6.")
    return name, int(substitution_types)


def paup_estimate_tree(
    alignment: str,
    output_tree: str,
    *,
    criterion: str = "likelihood",
    substitution_types: int = 6,
    unequal_base_freqs: bool = True,
    gamma_rates: bool = True,
    prop_invar: bool = True,
    paup_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate a tree from a DNA alignment with PAUP* (dendropy's estimate_tree) and save it as Newick.

    criterion is likelihood (heuristic search under the substitution model), parsimony, nj or upgma;
    substitution_types (1, 2 or 6) sets the model as in PAUP*'s lset nst. The model options apply to likelihood
    estimation. The program is looked up first; a missing PAUP* raises BackendUnavailableError.
    """
    dendropy = _dendropy()
    from dendropy.interop import paup as _paup

    source = _gap_cli.existing_file(alignment, "alignment")
    output_path = _gap_cli.new_file(output_tree, "output_tree")
    criterion_name, types = _check_estimation(criterion, substitution_types)
    options = {
        "unequal_base_freqs": _flag(unequal_base_freqs, "unequal_base_freqs"),
        "gamma_rates": _flag(gamma_rates, "gamma_rates"),
        "prop_invar": _flag(prop_invar, "prop_invar"),
    }
    executable = _paup_executable(paup_path)
    matrix = _read_alignment(dendropy, source)
    tree = _paup_call(
        "estimate_tree",
        lambda: _paup.estimate_tree(
            matrix,
            tree_est_criterion=criterion_name,
            num_subst=types,
            paup_path=executable,
            **options,
        ),
    )
    _write_tree(tree, output_path)
    return _paup_result(
        "estimate_tree",
        executable,
        [output_path],
        criterion=criterion_name,
        substitution_types=types,
        taxon_count=len(matrix),
    )


def paup_estimate_ultrametric_tree(
    alignment: str,
    output_tree: str,
    *,
    topology_tree: str | None = None,
    paup_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate an ultrametric (clock-constrained) tree with PAUP* (dendropy's estimate_ultrametric_tree).

    Without topology_tree the topology is a neighbour-joining tree, then branch lengths are estimated under a
    Jukes-Cantor-style model with a clock. With topology_tree (a Newick or NEXUS file) that topology is kept and
    only its branch lengths are estimated under the clock. Writes the tree as Newick with branch lengths.
    """
    dendropy = _dendropy()
    from dendropy.interop import paup as _paup

    source = _gap_cli.existing_file(alignment, "alignment")
    output_path = _gap_cli.new_file(output_tree, "output_tree")
    executable = _paup_executable(paup_path)
    matrix = _read_alignment(dendropy, source)
    topology = None
    if topology_tree is not None:
        topology_path = _gap_cli.existing_file(topology_tree, "topology_tree")
        topology = _read_trees(dendropy, topology_path, matrix.taxon_namespace)[0]
    tree = _paup_call(
        "estimate_ultrametric_tree",
        lambda: _paup.estimate_ultrametric_tree(
            matrix, topology_tree=topology, paup_path=executable
        ),
    )
    _write_tree(tree, output_path)
    return _paup_result(
        "estimate_ultrametric_tree",
        executable,
        [output_path],
        taxon_count=len(matrix),
        topology_given=topology is not None,
    )


def paup_estimate_model(
    alignment: str,
    output_tree: str,
    *,
    tree_file: str | None = None,
    criterion: str = "likelihood",
    substitution_types: int = 6,
    unequal_base_freqs: bool = True,
    gamma_rates: bool = True,
    prop_invar: bool = True,
    user_branch_lengths: bool = True,
    paup_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the substitution model parameters for a DNA alignment with PAUP* (dendropy's estimate_model).

    With tree_file (the first tree of a Newick or NEXUS file is used) the model is fitted on that tree; otherwise
    a tree is first estimated with criterion. Returns the fitted values from PAUP*'s report (log_likelihood,
    kappa, rAC...rGT, base frequencies pA...pT, prop_invar, alpha) as floats or None when not reported, and
    writes the tree that PAUP* used as Newick.
    """
    dendropy = _dendropy()
    from dendropy.interop import paup as _paup

    source = _gap_cli.existing_file(alignment, "alignment")
    output_path = _gap_cli.new_file(output_tree, "output_tree")
    criterion_name, types = _check_estimation(criterion, substitution_types)
    executable = _paup_executable(paup_path)
    matrix = _read_alignment(dendropy, source)
    model_tree = None
    if tree_file is not None:
        tree_path = _gap_cli.existing_file(tree_file, "tree_file")
        model_tree = _read_trees(dendropy, tree_path, matrix.taxon_namespace)[0]
    tree, estimates = _paup_call(
        "estimate_model",
        lambda: _paup.estimate_model(
            matrix,
            tree_model=model_tree,
            num_subst=types,
            unequal_base_freqs=_flag(unequal_base_freqs, "unequal_base_freqs"),
            gamma_rates=_flag(gamma_rates, "gamma_rates"),
            prop_invar=_flag(prop_invar, "prop_invar"),
            tree_est_criterion=criterion_name,
            tree_user_brlens=_flag(user_branch_lengths, "user_branch_lengths"),
            paup_path=executable,
        ),
    )
    _write_tree(tree, output_path)
    return _paup_result(
        "estimate_model",
        executable,
        [output_path],
        substitution_types=types,
        estimates=dict(estimates),
    )


def paup_pscore_trees(
    trees_file: str,
    alignment: str,
    *,
    pset_options: Sequence[str] | None = None,
    paup_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parsimony score (PAUP* pscore) of every tree in a file against a DNA alignment.

    pset_options are PAUP* pset settings written as ``pset <options>`` before pscore, each a plain token such as
    'gapmode=missing'. The options are checked against a conservative token pattern because dendropy writes them
    into the command text unescaped. dendropy's pscore_option_list is not exposed: it joins its options without a
    separator to the scorefile path, so they cannot be passed correctly. Returns one integer score per tree.
    """
    dendropy = _dendropy()
    from dendropy.interop import paup as _paup

    tree_path = _gap_cli.existing_file(trees_file, "trees_file")
    source = _gap_cli.existing_file(alignment, "alignment")
    tokens = [] if pset_options is None else list(pset_options)
    for token in tokens:
        if not isinstance(token, str) or not _PAUP_TOKEN.fullmatch(token):
            raise ConfigurationError(
                f"pset_options must be plain PAUP* tokens such as gapmode=missing; got {token!r}."
            )
    executable = _paup_executable(paup_path)
    matrix = _read_alignment(dendropy, source)
    trees = _read_trees(dendropy, tree_path, matrix.taxon_namespace)
    scores = _paup_call(
        "pscore",
        lambda: _paup.pscore_trees(
            trees,
            matrix,
            pset_option_list=tokens or None,
            paup_path=executable,
        ),
    )
    return _paup_result(
        "pscore", executable, [], scores=[int(value) for value in scores], tree_count=len(trees)
    )


def paup_prune_taxa(
    trees_file: str,
    taxa: Sequence[str],
    output_trees: str,
    *,
    paup_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete the named taxa from every tree in a file with PAUP* (dendropy's prune_taxa_from_trees).

    Each name must occur in the trees. The pruned trees are saved as a new NEXUS file with their taxa block.
    """
    dendropy = _dendropy()
    from dendropy.interop import paup as _paup

    tree_path = _gap_cli.existing_file(trees_file, "trees_file")
    output_path = _gap_cli.new_file(output_trees, "output_trees")
    if isinstance(taxa, (str, bytes)) or not isinstance(taxa, Sequence) or len(taxa) == 0:
        raise ConfigurationError("taxa must be a non-empty list of taxon labels.")
    namespace = dendropy.TaxonNamespace()
    trees = _read_trees(dendropy, tree_path, namespace)
    selected = []
    for label in taxa:
        taxon = namespace.get_taxon(_free_text(label, "taxa item"))
        if taxon is None:
            raise ConfigurationError(f"Taxon not found in the trees: {label}")
        selected.append(taxon)
    executable = _paup_executable(paup_path)
    pruned = _paup_call(
        "prune_taxa_from_trees",
        lambda: _paup.prune_taxa_from_trees(trees, selected, paup_path=executable),
    )
    with output_path.open("x") as handle:
        pruned.write_to_stream(handle, schema="nexus")
    return _paup_result(
        "prune_taxa_from_trees",
        executable,
        [output_path],
        removed=[str(t.label) for t in selected],
        tree_count=len(pruned),
    )


__all__ = [
    "align_with_muscle",
    "convert_format",
    "fetch_entrez_records",
    "format_to_fasta",
    "format_to_newick",
    "format_to_nexml",
    "format_to_nexus",
    "format_to_phylip",
    "paup_estimate_model",
    "paup_estimate_tree",
    "paup_estimate_ultrametric_tree",
    "paup_prune_taxa",
    "paup_pscore_trees",
    "run_dendropy_format",
    "summarize_tree_samples",
    "transfer_tree_labels",
]
_register_public_verbose(globals())
