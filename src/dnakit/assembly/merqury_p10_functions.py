"""Merqury command-line facades: k-mer based assembly evaluation, trio hap-mer analysis and k-mer database utilities.

Merqury (https://github.com/marbl/merqury) is a collection of shell scripts around Meryl, with R plots, Java helpers and,
for the tracks, samtools and the UCSC ``wigToBigWig``. The user installs it (``conda install -c bioconda merqury``, which
brings Meryl, R, Java, bedtools and samtools) and either sets the ``MERQURY`` environment variable to the Merqury folder,
puts ``merqury.sh`` on PATH, or passes ``merqury_home``. 1.4.1 is the tested version. Each public function runs one fixed
Merqury script (``bash <home>/<script> <arguments>``, no shell string, ``MERQURY`` set for the child) in a fresh output
directory. Merqury scripts read and write next to their inputs, so inputs (FASTA files, Meryl databases, tables) are
linked into the output directory first and the scripts see only those links; the input folders are never written to.
Results are JSON-serialisable dictionaries with status, backend, version (of Meryl, when it is needed), command, return
code, existing outputs, a bounded log tail and small parsed values. A missing script, program (meryl, java, Rscript,
samtools, wigToBigWig, bedtools, mash, pigz) raises ``BackendUnavailableError`` with an install hint.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p10_cli as _p10
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError

_TESTED_VERSION = "1.4.1"
_CITATION = "https://doi.org/10.1186/s13059-020-02134-9"
_HINT = (
    "Install Merqury 1.4.1 or later with its dependencies (`conda install -c bioconda merqury`), then set the MERQURY "
    "environment variable to the Merqury folder, put merqury.sh on PATH, or pass merqury_home."
)
_PROGRAM_HINTS = {
    "meryl": "conda install -c bioconda meryl",
    "java": "conda install -c conda-forge openjdk",
    "Rscript": "conda install -c conda-forge r-base r-ggplot2 r-argparse r-scales",
    "samtools": "conda install -c bioconda samtools",
    "wigToBigWig": "conda install -c bioconda ucsc-wigtobigwig",
    "bedGraphToBigWig": "conda install -c bioconda ucsc-bedgraphtobigwig",
    "bedtools": "conda install -c bioconda bedtools",
    "mash": "conda install -c bioconda mash",
    "pigz": "conda install -c conda-forge pigz",
    "bash": "bash is part of every Linux and macOS system",
}
_NAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*")
_FASTA_SUFFIXES = (".fasta", ".fa", ".fasta.gz", ".fa.gz")


# ---------------------------------------------------------------------------
# Locating Merqury and its programs
# ---------------------------------------------------------------------------


def _is_home(path: Path) -> bool:
    return (path / "util" / "util.sh").is_file() and (path / "eval").is_dir()


def _home(merqury_home: str | None) -> Path:
    """The Merqury folder: explicit path, then $MERQURY, then next to ``merqury.sh`` on PATH."""
    candidates: list[Path] = []
    if merqury_home is not None:
        candidates.append(Path(merqury_home).expanduser())
    elif os.environ.get("MERQURY"):
        candidates.append(Path(os.environ["MERQURY"]).expanduser())
    else:
        found = shutil.which("merqury.sh") or str(Path(sys.executable).parent / "merqury.sh")
        located = Path(found)
        if located.exists():
            resolved = located.resolve()
            candidates += [resolved.parent, located.parent.parent / "share" / "merqury"]
            candidates += sorted(located.parent.parent.glob("share/merqury*"))
    for candidate in candidates:
        if _is_home(candidate):
            return candidate.resolve()
    raise BackendUnavailableError("The Merqury scripts were not found.", hint=_HINT)


def _path_with(home: Path) -> str:
    """PATH for the child: the Merqury folder's own bin, the interpreter's bin and the caller's PATH."""
    directories = [str(Path(sys.executable).parent)]
    for parent in (home.parent.parent, home.parent.parent.parent):
        if (parent / "bin").is_dir():
            directories.append(str(parent / "bin"))
    return os.pathsep.join([*directories, os.environ.get("PATH", "")])


def _check_programs(home: Path, names: Sequence[str]) -> dict[str, str]:
    search = _path_with(home)
    found: dict[str, str] = {}
    for name in names:
        location = shutil.which(name, path=search)
        if location is None:
            raise BackendUnavailableError(
                f"{name} was not found on PATH.", hint=_PROGRAM_HINTS.get(name, _HINT)
            )
        found[name] = location
    return found


def _script(home: Path, relative: str) -> Path:
    path = home / relative
    if not path.is_file():
        raise BackendUnavailableError(f"This Merqury installation has no {relative}.", hint=_HINT)
    return path


def _meryl_version(programs: Mapping[str, str]) -> str | None:
    executable = programs.get("meryl")
    return None if executable is None else _cli._version_of(executable)


def _run(
    relative: str,
    arguments: Sequence[str],
    *,
    home: Path,
    cwd: Path,
    expected: Sequence[Path],
    needs: Sequence[str],
    interpreter: str = "bash",
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    extra_environment: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Run ``<interpreter> <home>/<relative> <arguments>`` in ``cwd``; completed when the exit code is 0 and ``expected`` exist."""
    script = _script(home, relative)
    programs = _check_programs(home, [interpreter, *needs])
    environment = {"MERQURY": str(home), "PATH": _path_with(home)}
    if extra_environment:
        environment.update(extra_environment)
    result = execute_bounded_command(
        programs[interpreter],
        [str(script), *arguments],
        backend_id="merqury",
        cwd=cwd,
        timeout_seconds=timeout_seconds,
        max_output_bytes=_cli.MAX_OUTPUT_BYTES,
        monitored_output_paths=list(expected),
        max_monitored_output_bytes=_cli.MAX_ARTIFACT_BYTES,
        environment_overrides=environment,
        verbose=False,
    )
    produced = [str(path) for path in expected if path.exists()]
    complete = result.return_code == 0 and len(produced) == len(expected)
    return {
        "status": "completed" if complete else "failed",
        "backend": "merqury",
        "version": _meryl_version(programs),
        "command": ["merqury", relative, *arguments],
        "return_code": result.return_code,
        "outputs": produced,
        "log_tail": (result.output or "")[-_cli.LOG_TAIL_CHARS :],
        "_text": result.output or "",
    }


def _finish(
    operation: str,
    rendered: dict[str, Any],
    home: Path,
    parse: Any = None,
    **extra: Any,
) -> dict[str, Any]:
    result = {key: value for key, value in rendered.items() if key != "_text"}
    result["operation"] = operation
    result["merqury_home"] = str(home)
    result["provenance"] = _p10.provenance("Merqury", _TESTED_VERSION, None, _CITATION) | {
        "meryl_version": rendered.get("version")
    }
    if result["status"] == "completed":
        result.update(extra)
        if parse is not None:
            result.update(parse(rendered["_text"]))
    return result


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def _file(path: str | Path, name: str) -> Path:
    return _cli.existing_file(path, name)


def _meryl_db(path: str | Path, name: str) -> Path:
    """An existing Meryl database (a directory)."""
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise ConfigurationError(f"{name} must be an existing Meryl database directory: {resolved}")
    return resolved


def _fasta(path: str | Path, name: str) -> Path:
    resolved = _file(path, name)
    if not resolved.name.endswith(_FASTA_SUFFIXES):
        raise ConfigurationError(
            f"{name} must be named *.fa, *.fasta, *.fa.gz or *.fasta.gz: {resolved.name}"
        )
    return resolved


def _meryl_named(path: str | Path, name: str) -> Path:
    resolved = _meryl_db(path, name)
    if not resolved.name.endswith(".meryl"):
        raise ConfigurationError(f"{name} must be a directory named *.meryl: {resolved.name}")
    return resolved


def _stage(directory: Path, *paths: Path) -> list[str]:
    """Link inputs into ``directory`` under their own names and return the names (no two inputs may share a name)."""
    names: list[str] = []
    for path in paths:
        if path.name in names:
            raise ConfigurationError(
                f"Two inputs have the same file name {path.name!r}; rename one of them."
            )
        link = directory / path.name
        link.symlink_to(path)
        names.append(path.name)
    return names


def _name(value: str, name: str) -> str:
    text = _cli.text(value, name)
    if _NAME.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} may contain only letters, digits, '_', '.' and '-': {text!r}"
        )
    return text


def _stem(file_name: str) -> str:
    for suffix in (".gz", ".fasta", ".fa"):
        if file_name.endswith(suffix):
            file_name = file_name[: -len(suffix)]
    return file_name


def _new_directory(path: str | Path, home: Path | None = None) -> Path:
    return _cli.new_directory(path, "output_dir")


def _read_table(path: Path, *, header: bool = True, limit: int = 20) -> dict[str, Any]:
    return _p10.table_preview(path, header=header, comment=None, limit=limit)


def _number(text: str) -> int | float | str:
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


# ---------------------------------------------------------------------------
# eval/
# ---------------------------------------------------------------------------


def merqury_per_sequence_qv(
    assembly_fasta: str,
    read_meryl: str,
    output_dir: str,
    *,
    prefix: str = "per_seq",
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the consensus quality value (QV) of every sequence of an assembly from a read k-mer database (``eval/per_seq_qv.sh``).

    ``read_meryl`` is the Meryl database of the reads. The assembly is counted into k-mers of the same size, the
    k-mers found only in the assembly are looked up per sequence, and QV = -10 log10(error rate) with the
    error rate 1 - (1 - e/n)^(1/k). Output ``<prefix>.qv`` in the new ``output_dir`` has the columns sequence,
    assembly-only k-mers, k-mers, QV and error rate. Result keys: ``qv_file``, ``columns``, ``row_count``, ``rows``.
    """
    fasta = _fasta(assembly_fasta, "assembly_fasta")
    reads = _meryl_db(read_meryl, "read_meryl")
    base = _name(prefix, "prefix")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "bash"])
    directory = _new_directory(output_dir)
    fasta_name, reads_name = _stage(directory, fasta, reads)
    qv = directory / f"{base}.qv"
    rendered = _run(
        "eval/per_seq_qv.sh",
        [fasta_name, reads_name, base],
        home=home,
        cwd=directory,
        expected=[qv],
        needs=["meryl", "meryl-lookup", "awk"],
    )

    def parse(_: str) -> dict[str, Any]:
        table = _p10.table_preview(qv, header=False, comment=None)
        table["columns"] = ["sequence", "error_kmers", "total_kmers", "qv", "error_rate"]
        return {"qv_file": str(qv), **table}

    return _finish("merqury_per_sequence_qv", rendered, home, parse)


def merqury_kmer_multiplicity_track(
    assembly_fasta: str,
    meryl_db: str,
    output_dir: str,
    *,
    source: str = "assembly",
    prefix: str = "track",
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a k-mer count track over the assembly sequences as WIG and bigWig (``eval/asm_multiplicity.sh`` and ``eval/read_multiplicity.sh``).

    ``source="assembly"`` (``asm_multiplicity.sh``) counts how often each assembly k-mer occurs in the assembly
    itself (copy number; ``meryl_db`` is the Meryl database of the assembly); ``source="reads"``
    (``read_multiplicity.sh``) reports the read multiplicity of each assembly k-mer (``meryl_db`` is the read
    database). The files are ``<prefix>.copies.wig`` / ``<prefix>.asm_multiplicity.bigWig`` or
    ``<prefix>.read.wig`` / ``<prefix>.read_multiplicity.bigwig``. Needs samtools and wigToBigWig.
    """
    kind = _cli.choice(source, ("assembly", "reads"), "source")
    fasta = _fasta(assembly_fasta, "assembly_fasta")
    database = _meryl_db(meryl_db, "meryl_db")
    base = _name(prefix, "prefix")
    home = _home(merqury_home)
    _check_programs(home, ["meryl-lookup", "samtools", "wigToBigWig", "bash"])
    directory = _new_directory(output_dir)
    fasta_name, db_name = _stage(directory, fasta, database)
    # asm_multiplicity.sh tests for the FASTA (not the index) before indexing, so the index is made here
    import subprocess

    subprocess.run(  # noqa: S603 - fixed argument vector, no shell
        [_check_programs(home, ["samtools"])["samtools"], "faidx", fasta_name],
        cwd=directory, check=False, capture_output=True, timeout=600,
    )
    if kind == "assembly":
        script, wig, big = (
            "eval/asm_multiplicity.sh",
            f"{base}.copies.wig",
            f"{base}.asm_multiplicity.bigWig",
        )
    else:
        script, wig, big = (
            "eval/read_multiplicity.sh",
            f"{base}.read.wig",
            f"{base}.read_multiplicity.bigwig",
        )
    rendered = _run(
        script,
        [fasta_name, db_name, base],
        home=home,
        cwd=directory,
        expected=[directory / wig, directory / big],
        needs=["meryl-lookup", "samtools", "wigToBigWig"],
    )
    return _finish(
        "merqury_kmer_multiplicity_track",
        rendered,
        home,
        source=kind,
        wig_file=str(directory / wig),
        bigwig_file=str(directory / big),
    )


def merqury_false_duplications(
    spectra_cn_hist: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count assembly k-mers that are extra copies in the single- and two-copy read peaks (``eval/false_duplications.sh``).

    ``spectra_cn_hist`` is the ``*.spectra-cn.hist`` table of a pseudo-haplotype or mixed-haplotype assembly
    (columns Copies, kmer_multiplicity, Count). The cutoff is 1.5 times the read multiplicity at the 1-copy
    peak. Nothing is written. Result keys: ``cutoff``, ``copies`` (k-mer counts for 1, 2, 3, 4 and >4 copies),
    ``duplicated_kmers`` (copies above 1), ``all_kmers`` and ``duplicated_percent``.
    """
    hist = _file(spectra_cn_hist, "spectra_cn_hist")
    home = _home(merqury_home)
    rendered = _run(
        "eval/false_duplications.sh",
        [str(hist)],
        home=home,
        cwd=hist.parent,
        expected=[],
        needs=["awk"],
    )

    def parse(text: str) -> dict[str, Any]:
        lines = [line for line in text.splitlines() if line.strip()]
        row = next(
            (line.split("\t") for line in reversed(lines) if line.startswith(str(hist))), None
        )
        if row is None or len(row) < 10:
            raise DNAKitError(f"Unexpected output of false_duplications.sh: {text[-300:]!r}")
        values = [_number(item) for item in row[1:]]
        return {
            "cutoff": values[0],
            "copies": dict(zip(("1", "2", "3", "4", ">4"), values[1:6], strict=True)),
            "duplicated_kmers": values[6],
            "all_kmers": values[7],
            "duplicated_percent": values[8],
        }

    return _finish("merqury_false_duplications", rendered, home, parse)


def merqury_false_duplications_track(
    spectra_cn_hist: str,
    assembly_fasta: str,
    assembly_meryl: str,
    reads_meryl: str,
    output_dir: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate the extra-copy k-mers of ``merqury_false_duplications`` along the assembly as BED and WIG tracks (``eval/false_duplications_track.sh``).

    ``assembly_meryl`` is the assembly's Meryl database (the ``<assembly>.meryl`` that Merqury builds) and
    ``reads_meryl`` the read database. Outputs in the new ``output_dir``: ``<assembly>.fd.bed`` (merged
    intervals) and ``<assembly>.fd.wig``. Result keys: ``bed_file``, ``wig_file``, ``interval_count``.
    """
    hist = _file(spectra_cn_hist, "spectra_cn_hist")
    fasta = _fasta(assembly_fasta, "assembly_fasta")
    asm_db = _meryl_named(assembly_meryl, "assembly_meryl")
    reads = _meryl_db(reads_meryl, "reads_meryl")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "meryl-lookup", "bedtools", "bash"])
    directory = _new_directory(output_dir)
    expected_db = f"{_stem(fasta.name)}.meryl"
    if asm_db.name != expected_db:
        raise ConfigurationError(f"assembly_meryl must be named {expected_db} for {fasta.name}.")
    hist_name, fasta_name, _asm_name, reads_name = _stage(directory, hist, fasta, asm_db, reads)
    stem = _stem(fasta_name)
    bed, wig = directory / f"{stem}.fd.bed", directory / f"{stem}.fd.wig"
    rendered = _run(
        "eval/false_duplications_track.sh",
        [hist_name, fasta_name, reads_name],
        home=home,
        cwd=directory,
        expected=[bed, wig],
        needs=["meryl", "meryl-lookup", "bedtools"],
    )

    def parse(_: str) -> dict[str, Any]:
        return {
            "bed_file": str(bed),
            "wig_file": str(wig),
            "interval_count": _p10.count_records(bed, comment=None),
        }

    return _finish("merqury_false_duplications_track", rendered, home, parse)


def merqury_mash_qv(
    assembly_fasta: str,
    read_paths: Sequence[str],
    output_dir: str,
    *,
    k: int = 21,
    cpus: int = 1,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Quick assembly QV estimate from a Mash sketch of the assembly screened against the reads (``merqury-mash.sh``).

    ``k`` is the k-mer size and ``read_paths`` the read files (FASTA/FASTQ, possibly gzipped); they are listed in a
    file of file names as the script expects. Mash screen identities are converted to QV; the values do not match
    the full Merqury run exactly (the script says so). Output files ``<assembly>.msh.idy`` and ``<assembly>.msh.qv``
    are written next to the linked assembly. Result keys: ``qv_file``, ``columns``, ``rows``.
    """
    fasta = _fasta(assembly_fasta, "assembly_fasta")
    if isinstance(read_paths, str) or not read_paths:
        raise ConfigurationError("read_paths must be a non-empty list of read files.")
    reads = [_file(path, "read_paths") for path in read_paths]
    size = _cli.integer(k, "k", minimum=3, maximum=63)
    threads = _cli.integer(cpus, "cpus", minimum=1, maximum=256)
    home = _home(merqury_home)
    _check_programs(home, ["mash", "bash"])
    directory = _new_directory(output_dir)
    (fasta_name,) = _stage(directory, fasta)
    fofn = directory / "reads.fofn"
    fofn.write_text("".join(f"{path}\n" for path in reads), encoding="utf-8")
    stem = _stem(fasta_name)
    qv = directory / f"{stem}.msh.qv"
    rendered = _run(
        "merqury-mash.sh",
        [str(size), fasta_name, fofn.name, str(threads)],
        home=home,
        cwd=directory,
        expected=[qv, directory / f"{stem}.msh.idy"],
        needs=["mash", "awk"],
    )

    def parse(_: str) -> dict[str, Any]:
        table = _p10.table_preview(qv, header=False, comment=None)
        table["columns"] = ["assembly", "reads", "qv", "error_rate"]
        return {"qv_file": str(qv), **table}

    return _finish("merqury_mash_qv", rendered, home, parse)


def merqury_best_k(
    genome_size: int,
    *,
    tolerable_collision_rate: float = 0.001,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Choose the k-mer size for a genome size and a tolerable k-mer collision rate (``best_k.sh``; Fofanov et al. 2004).

    ``genome_size`` is the haploid (or diploid, depending on what is evaluated) genome size in bp.
    ``tolerable_collision_rate`` is the error rate of the read set (Merqury's default 0.001 suits Illumina WGS).
    k = log4(genome_size * (1 - rate) / rate). Result keys: ``k`` (the real value printed by Merqury) and
    ``k_rounded_up``.
    """
    size = _cli.integer(genome_size, "genome_size", minimum=1)
    rate = _cli.number(
        tolerable_collision_rate, "tolerable_collision_rate", minimum=1e-12, maximum=0.5
    )
    home = _home(merqury_home)
    with_cwd = Path(os.getcwd()) if Path(os.getcwd()).is_dir() else Path("/")
    rendered = _run(
        "best_k.sh",
        [str(size), format(rate, ".12g")],
        home=home,
        cwd=with_cwd,
        expected=[],
        needs=["awk"],
    )

    def parse(text: str) -> dict[str, Any]:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        value = float(lines[-1])
        return {"k": value, "k_rounded_up": int(-(-value // 1))}

    return _finish(
        "merqury_best_k", rendered, home, parse, genome_size=size, tolerable_collision_rate=rate
    )


# ---------------------------------------------------------------------------
# build/
# ---------------------------------------------------------------------------


def _fofn(directory: Path, paths: Sequence[Path], name: str = "input.fofn") -> Path:
    fofn = directory / name
    fofn.write_text("".join(f"{path}\n" for path in paths), encoding="utf-8")
    return fofn


def merqury_split_10x_fastq(
    read_paths: Sequence[str],
    output_dir: str,
    *,
    threads: int = 8,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim the first 23 bases of 10x R1 reads and split each file into pieces of 300 million lines, gzip-compressed (``build/split_10x.sh``).

    ``read_paths`` are FASTQ R1 files (gzipped or not). The pieces go to ``<output_dir>/split/<name>.NNNN.fq.gz``;
    ``<fofn>.<n>`` lists the pieces of the n-th input. Needs pigz. Result keys: ``pieces`` (per input file).
    """
    reads = (
        [_file(path, "read_paths") for path in read_paths]
        if not isinstance(read_paths, str)
        else []
    )
    if not reads:
        raise ConfigurationError("read_paths must be a non-empty list of FASTQ files.")
    cpu = _cli.integer(threads, "threads", minimum=1, maximum=256)
    home = _home(merqury_home)
    _check_programs(home, ["pigz", "awk", "split", "bash"])
    directory = _new_directory(output_dir)
    fofn = _fofn(directory, reads)
    pieces: dict[str, list[str]] = {}
    last: dict[str, Any] | None = None
    text = ""
    for index, path in enumerate(reads, start=1):
        listing = directory / f"{fofn.name}.{index}"
        step = _run(
            "build/split_10x.sh",
            [fofn.name, str(index)],
            home=home,
            cwd=directory,
            expected=[listing],
            needs=["pigz", "awk", "split"],
            extra_environment={"SLURM_CPUS_PER_TASK": str(cpu)},
        )
        text += step["_text"]
        if last is None:
            last = step
        else:
            last["outputs"] += step["outputs"]
            last["return_code"] = last["return_code"] or step["return_code"]
            last["status"] = (
                "completed" if last["status"] == step["status"] == "completed" else "failed"
            )
            last["log_tail"] = (last["log_tail"] + step["log_tail"])[-_cli.LOG_TAIL_CHARS :]
        if listing.is_file():
            pieces[str(path)] = [
                line for line in listing.read_text(encoding="utf-8").splitlines() if line
            ]
    assert last is not None
    last["_text"] = text
    return _finish("merqury_split_10x_fastq", last, home, pieces=pieces)


def merqury_filter_solid_kmers(
    meryl_db: str,
    output_dir: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop the erroneous low-count k-mers of a read Meryl database using the ploidy-derived threshold (``build/filt.sh``).

    The k-mer histogram is written (``<db>.hist``), ``kmerHistToPloidyDepth.jar`` derives the ploidy and the
    count threshold (``<db>.hist.ploidy``) and ``meryl greater-than`` keeps the k-mers above it as
    ``<db>.gt<threshold>.meryl``. Result keys: ``threshold``, ``solid_database``, ``histogram_file``, ``ploidy_file``.
    """
    database = _meryl_named(meryl_db, "meryl_db")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "java", "bash"])
    directory = _new_directory(output_dir)
    (name,) = _stage(directory, database)
    stem = name[: -len(".meryl")]
    filt = directory / f"{stem}.filt"
    rendered = _run(
        "build/filt.sh",
        [name],
        home=home,
        cwd=directory,
        expected=[filt, directory / f"{stem}.hist", directory / f"{stem}.hist.ploidy"],
        needs=["meryl", "java"],
    )

    def parse(_: str) -> dict[str, Any]:
        threshold = int(filt.read_text(encoding="utf-8").split()[0])
        return {
            "threshold": threshold,
            "solid_database": str(directory / f"{stem}.gt{threshold}.meryl"),
            "histogram_file": str(directory / f"{stem}.hist"),
            "ploidy_file": str(directory / f"{stem}.hist.ploidy"),
        }

    return _finish("merqury_filter_solid_kmers", rendered, home, parse)


def merqury_combine_kmer_dbs(
    operation: str,
    meryl_dbs: Sequence[str],
    output_dir: str,
    output_name: str,
    *,
    k: int | None = None,
    filter_solid: bool = True,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine Meryl k-mer databases (``build/diff.sh``, ``build/intersect.sh`` and ``build/union_sum.sh``).

    ``operation`` selects the script. ``"difference"`` (``diff.sh``) takes exactly two databases and keeps the k-mers of
    the first that are not in the second (``<db1>_not_<db2>.meryl``); with ``filter_solid=True`` (Merqury's
    default) the result is then filtered with ``merqury_filter_solid_kmers``' rule and linked as
    ``<output_name>.meryl``, otherwise (``-no-filt``) only linked. ``"intersect"`` (``intersect.sh``) takes two
    databases and keeps the k-mers of both with the counts of the first (``<db1>_and_<db2>.meryl``, linked as
    ``<output_name>.meryl``). ``"union_sum"`` (``union_sum.sh``) takes two or more databases and the k-mer size ``k`` and
    writes the union with summed counts as ``<output_name>.k<k>.meryl``. Result key: ``database``.
    """
    kind = _cli.choice(operation, ("difference", "intersect", "union_sum"), "operation")
    if isinstance(meryl_dbs, str) or not meryl_dbs:
        raise ConfigurationError("meryl_dbs must be a list of Meryl database directories.")
    databases = [_meryl_named(path, "meryl_dbs") for path in meryl_dbs]
    label = _name(output_name, "output_name")
    if kind in ("difference", "intersect") and len(databases) != 2:
        raise ConfigurationError(f"operation={kind!r} needs exactly two databases.")
    if kind == "union_sum":
        if len(databases) < 2:
            raise ConfigurationError("operation='union_sum' needs at least two databases.")
        if k is None:
            raise ConfigurationError("operation='union_sum' needs the k-mer size k.")
        size = _cli.integer(k, "k", minimum=3, maximum=63)
    home = _home(merqury_home)
    _check_programs(
        home, ["meryl", "bash"] + (["java"] if kind == "difference" and filter_solid else [])
    )
    directory = _new_directory(output_dir)
    names = _stage(directory, *databases)
    if kind == "difference":
        arguments = [names[0], names[1], label] + ([] if filter_solid else ["-no-filt"])
        script, result_name = "build/diff.sh", f"{label}.meryl"
        needs = ["meryl"] + (["java"] if filter_solid else [])
    elif kind == "intersect":
        arguments = [names[0], names[1], label]
        script, result_name = "build/intersect.sh", f"{label}.meryl"
        needs = ["meryl"]
    else:
        fofn = _fofn(directory, [directory / name for name in names], "meryl.list")
        arguments = [str(size), fofn.name, label]
        script, result_name = "build/union_sum.sh", f"{label}.k{size}.meryl"
        needs = ["meryl"]
    target = directory / result_name
    rendered = _run(script, arguments, home=home, cwd=directory, expected=[target], needs=needs)
    return _finish(
        "merqury_combine_kmer_dbs", rendered, home, combine_operation=kind, database=str(target.resolve())
    )


def merqury_kmer_histogram_table(
    meryl_dbs: Mapping[str, str],
    output_path: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export the k-mer multiplicity histograms of Meryl databases as the table that Merqury's plot script reads (``plot/to_hist_for_plotting.sh``).

    ``meryl_dbs`` maps a label to a Meryl database. The TSV ``output_path`` has the columns ``kmer`` (the label),
    ``kmer_multiplicity`` and ``Count``. Result keys: ``columns``, ``row_count``, ``rows``.
    """
    if not isinstance(meryl_dbs, Mapping) or not meryl_dbs:
        raise ConfigurationError("meryl_dbs must be a non-empty {label: meryl database} mapping.")
    pairs = [
        (_name(label, "label"), _meryl_db(path, "meryl_dbs")) for label, path in meryl_dbs.items()
    ]
    out = _cli.new_file(output_path, "output_path")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "awk", "bash"])
    arguments: list[str] = []
    for label, path in pairs:
        arguments += [str(path), label]
    program = _check_programs(home, ["bash"])["bash"]
    rendered = _p10.stdout_to_file(
        "merqury",
        [str(_script(home, "plot/to_hist_for_plotting.sh")), *arguments],
        output=out,
        cwd=out.parent,
        explicit=None,
        hint=_HINT,
        program_path=program,
        environment={"MERQURY": str(home), "PATH": _path_with(home)},
    )
    rendered["_text"] = rendered["log_tail"]
    rendered["version"] = None
    return _finish(
        "merqury_kmer_histogram_table",
        rendered,
        home,
        lambda _: _p10.table_preview(out, header=True, comment=None),
    )


# ---------------------------------------------------------------------------
# trio/
# ---------------------------------------------------------------------------


def merqury_hapmers(
    hap1_meryl: str,
    hap2_meryl: str,
    output_dir: str,
    *,
    child_meryl: str | None = None,
    no_filt: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Derive haplotype-specific k-mers (hap-mers) from parental k-mer databases, and the inherited ones from a child (``trio/hapmers.sh``).

    ``hap1_meryl`` and ``hap2_meryl`` hold the k-mers of both parents (all built with the same k). The script
    writes ``<hap1>.only.meryl``, ``<hap2>.only.meryl`` and ``shrd.meryl``; with ``child_meryl`` it adds the
    inherited hap-mers ``<hap>.hapmer.meryl``, ``shrd.inherited.meryl``, ``cutoffs.txt`` and the
    ``inherited_hapmers`` spectra plot. ``no_filt=True`` skips the solid k-mer filtering of the parental
    databases (use only when they do not come from regular sequencing). Needs Java and, with a child, R.
    Result keys: ``hapmer_databases`` (the .meryl directories that exist), ``cutoffs`` (parsed ``cutoffs.txt``).
    """
    first = _meryl_named(hap1_meryl, "hap1_meryl")
    second = _meryl_named(hap2_meryl, "hap2_meryl")
    child = None if child_meryl is None else _meryl_named(child_meryl, "child_meryl")
    home = _home(merqury_home)
    needs = ["meryl", "java"] + (["Rscript"] if child is not None else [])
    _check_programs(home, [*needs, "bash"])
    directory = _new_directory(output_dir)
    names = _stage(directory, first, second, *([child] if child is not None else []))
    first_stem, second_stem = (name[: -len(".meryl")] for name in names[:2])
    expected = [
        directory / f"{first_stem}.only.meryl",
        directory / f"{second_stem}.only.meryl",
        directory / "shrd.meryl",
    ]
    if child is not None:
        expected += [
            directory / f"{first_stem}.hapmer.meryl",
            directory / f"{second_stem}.hapmer.meryl",
            directory / "cutoffs.txt",
        ]
    arguments = list(names) + (["-no-filt"] if no_filt else [])
    rendered = _run(
        "trio/hapmers.sh", arguments, home=home, cwd=directory, expected=expected, needs=needs
    )

    def parse(_: str) -> dict[str, Any]:
        databases = sorted(str(path) for path in directory.glob("*.meryl") if path.name not in names)
        cutoffs: dict[str, int] = {}
        cut = directory / "cutoffs.txt"
        if cut.is_file():
            for line in cut.read_text(encoding="utf-8").splitlines():
                parts = line.split()
                if len(parts) == 2 and parts[1].isdigit():
                    cutoffs[parts[0]] = int(parts[1])
        return {"hapmer_databases": databases, "cutoffs": cutoffs}

    return _finish("merqury_hapmers", rendered, home, parse)


def merqury_hapmer_spectra(
    reads_meryl: str,
    hap1_meryl: str,
    hap2_meryl: str,
    assemblies: Sequence[str],
    output_dir: str,
    output_prefix: str,
    *,
    homopolymer_compressed: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hap-mer copy-number spectra and hap-mer completeness of one or two assemblies (``trio/spectra-hap.sh``).

    ``reads_meryl`` holds all read k-mers; ``hap1_meryl`` and ``hap2_meryl`` the haplotype-specific k-mers with their
    read counts (for example ``mat.hapmer.meryl`` and ``pat.hapmer.meryl``); ``assemblies`` is one or two FASTA
    files. Per assembly and haplotype it writes ``<prefix>.<asm>.<hap>.spectra-cn.hist`` and plots (R),
    and ``<prefix>.completeness.stats`` (assembly, haplotype, hap-mers found, hap-mers total, percent).
    ``homopolymer_compressed`` is the script's ``-c`` flag. Result keys: ``completeness`` (rows) and ``histograms``.
    """
    reads = _meryl_named(reads_meryl, "reads_meryl")
    first = _meryl_named(hap1_meryl, "hap1_meryl")
    second = _meryl_named(hap2_meryl, "hap2_meryl")
    if isinstance(assemblies, str) or not 1 <= len(assemblies) <= 2:
        raise ConfigurationError("assemblies must list one or two FASTA files.")
    fastas = [_fasta(path, "assemblies") for path in assemblies]
    prefix = _name(output_prefix, "output_prefix")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "Rscript", "bash"])
    directory = _new_directory(output_dir)
    names = _stage(directory, reads, first, second, *fastas)
    stats = directory / f"{prefix}.completeness.stats"
    arguments = (["-c"] if homopolymer_compressed else []) + [*names, prefix]
    rendered = _run(
        "trio/spectra-hap.sh",
        arguments,
        home=home,
        cwd=directory,
        expected=[stats],
        needs=["meryl", "Rscript"],
    )

    def parse(_: str) -> dict[str, Any]:
        rows = []
        for line in stats.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) == 5:
                rows.append(
                    {
                        "assembly": parts[0],
                        "haplotype": parts[1],
                        "found": _number(parts[2]),
                        "total": _number(parts[3]),
                        "percent": _number(parts[4]),
                    }
                )
        return {
            "completeness": rows,
            "histograms": sorted(
                str(path) for path in directory.glob(f"{prefix}.*.spectra-hap-cn.hist")
            ),
        }

    return _finish("merqury_hapmer_spectra", rendered, home, parse)


def merqury_hapmer_counts(
    hap1_meryl: str,
    hap2_meryl: str,
    assemblies: Sequence[str],
    output_dir: str,
    output_prefix: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the hap-mers of each haplotype in every assembly contig and draw the blob plot (``trio/hap_blob.sh``).

    ``hap1_meryl`` and ``hap2_meryl`` are the hap-mer databases; ``assemblies`` one or two FASTA files. Writes the
    table ``<prefix>.hapmers.count`` (Assembly, Contig, hap1 count, hap2 count, Size) and the plot
    ``<prefix>.hapmers.blob.png`` (R). Result keys: ``count_file``, ``columns``, ``row_count``, ``rows``.
    """
    first = _meryl_named(hap1_meryl, "hap1_meryl")
    second = _meryl_named(hap2_meryl, "hap2_meryl")
    if isinstance(assemblies, str) or not 1 <= len(assemblies) <= 2:
        raise ConfigurationError("assemblies must list one or two FASTA files.")
    fastas = [_fasta(path, "assemblies") for path in assemblies]
    prefix = _name(output_prefix, "output_prefix")
    home = _home(merqury_home)
    _check_programs(home, ["meryl-lookup", "Rscript", "bash"])
    directory = _new_directory(output_dir)
    names = _stage(directory, first, second, *fastas)
    count = directory / f"{prefix}.hapmers.count"
    rendered = _run(
        "trio/hap_blob.sh",
        [*names, prefix],
        home=home,
        cwd=directory,
        expected=[count, directory / f"{prefix}.hapmers.blob.png"],
        needs=["meryl-lookup", "Rscript", "awk"],
    )
    return _finish(
        "merqury_hapmer_counts",
        rendered,
        home,
        lambda _: {
            "count_file": str(count),
            **_p10.table_preview(count, header=True, comment=None),
        },
    )


def merqury_hamming_error(
    hapmer_count_table: str,
    assemblies: Sequence[str],
    output_dir: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hamming (phasing) error of assemblies from the hap-mer counts per contig (``trio/hamming_error.sh``).

    ``hapmer_count_table`` is the ``*.hapmers.count`` file of ``merqury_hapmer_counts`` and ``assemblies`` the one or
    two FASTA files (their names, without suffix, match the first column of the table). Writes
    ``<table>.ham`` (per contig) and ``phase_hammingerror.stats`` in the new ``output_dir``. Result keys:
    ``rows`` (Assembly, Hap1, Hap2, Hap1Error(%), Hap2Error(%), MajorHapError(%)).
    """
    table = _file(hapmer_count_table, "hapmer_count_table")
    if isinstance(assemblies, str) or not 1 <= len(assemblies) <= 2:
        raise ConfigurationError("assemblies must list one or two FASTA files.")
    names = [_fasta(path, "assemblies").name for path in assemblies]
    home = _home(merqury_home)
    _check_programs(home, ["awk", "bash"])
    directory = _new_directory(output_dir)
    (table_name,) = _stage(directory, table)
    stats = directory / "phase_hammingerror.stats"
    rendered = _run(
        "trio/hamming_error.sh",
        [table_name, *names],
        home=home,
        cwd=directory,
        expected=[stats, directory / f"{table_name}.ham"],
        needs=["awk"],
    )
    return _finish(
        "merqury_hamming_error",
        rendered,
        home,
        lambda _: {
            "stats_file": str(stats),
            **_p10.table_preview(stats, header=True, comment=None),
        },
    )


def merqury_assign_reads(
    read_path: str,
    hap1_name: str,
    hap2_name: str,
    hap1_hapmer_meryl: str,
    hap2_hapmer_meryl: str,
    output_dir: str,
    *,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign reads to haplotype 1, haplotype 2 or unknown by their hap-mer content (``trio/assign_reads.sh``).

    One read file (FASTQ/FASTA, gzipped or not). For each read the output ``<read file>.assigned`` has the read
    id, k-mer count, hap-1 k-mers and how many were found, the same for hap 2, and the assignment (``hap1_name``,
    ``hap2_name`` or ``unk``). Result keys: ``assignments`` (file), ``counts`` (reads per assignment).
    """
    reads = _file(read_path, "read_path")
    one = _name(hap1_name, "hap1_name")
    two = _name(hap2_name, "hap2_name")
    first = _meryl_named(hap1_hapmer_meryl, "hap1_hapmer_meryl")
    second = _meryl_named(hap2_hapmer_meryl, "hap2_hapmer_meryl")
    home = _home(merqury_home)
    _check_programs(home, ["meryl-lookup", "awk", "bash"])
    directory = _new_directory(output_dir)
    names = _stage(directory, reads, first, second)
    out = directory / f"{reads.name}.assigned"
    rendered = _run(
        "trio/assign_reads.sh",
        [names[0], one, two, names[1], names[2]],
        home=home,
        cwd=directory,
        expected=[out],
        needs=["meryl-lookup", "awk"],
    )

    def parse(_: str) -> dict[str, Any]:
        counts: dict[str, int] = {}
        with out.open(encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[-1]:
                    counts[parts[-1]] = counts.get(parts[-1], 0) + 1
        return {"assignments": str(out), "counts": dict(sorted(counts.items()))}

    return _finish("merqury_assign_reads", rendered, home, parse)


def merqury_exclude_reads(
    meryl_db: str,
    read_pairs: Sequence[Sequence[str]],
    output_dir: str,
    *,
    tenx: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove every read (and its mate) that contains a k-mer of a Meryl database (``trio/exclude_reads.sh``).

    ``read_pairs`` is a list of ``[R1]`` or ``[R1, R2]`` FASTQ paths; the filtered reads are written under the new
    ``output_dir`` with the input file names. ``tenx=True`` ignores the first 23 bases of each R1 (10x barcodes).
    Result key: ``filtered_reads``.
    """
    database = _meryl_named(meryl_db, "meryl_db")
    if isinstance(read_pairs, str) or not read_pairs:
        raise ConfigurationError("read_pairs must be a non-empty list of [R1] or [R1, R2] paths.")
    pairs: list[list[Path]] = []
    for pair in read_pairs:
        if isinstance(pair, str) or not 1 <= len(pair) <= 2:
            raise ConfigurationError("each entry of read_pairs must be [R1] or [R1, R2].")
        pairs.append([_file(path, "read_pairs") for path in pair])
    home = _home(merqury_home)
    _check_programs(home, ["meryl-lookup", "awk", "bash"])
    directory = _new_directory(output_dir)
    (db_name,) = _stage(directory, database)
    map_file = directory / "input.map"
    map_file.write_text(
        "".join("\t".join(str(path) for path in pair) + "\n" for pair in pairs), encoding="utf-8"
    )
    out_directory = directory / "filtered"
    out_directory.mkdir()
    produced: list[str] = []
    final: dict[str, Any] | None = None
    text = ""
    for index, pair in enumerate(pairs, start=1):
        targets = [out_directory / str(path).lstrip("/") for path in pair]
        for target in targets:
            target.parent.mkdir(parents=True, exist_ok=True)
        step = _run(
            "trio/exclude_reads.sh",
            [db_name, map_file.name, "filtered", *(["-10x"] if tenx else []), str(index)],
            home=home,
            cwd=directory,
            expected=targets,
            needs=["meryl-lookup", "awk"],
        )
        text += step["_text"]
        if final is None:
            final = step
        else:
            final["outputs"] += step["outputs"]
            final["return_code"] = final["return_code"] or step["return_code"]
            final["status"] = (
                "completed" if final["status"] == step["status"] == "completed" else "failed"
            )
        produced += step["outputs"]
    assert final is not None
    final["_text"] = text
    return _finish("merqury_exclude_reads", final, home, filtered_reads=produced)


def merqury_phase_blocks(
    assembly_fasta: str,
    hap1_meryl: str,
    hap2_meryl: str,
    output_dir: str,
    output_prefix: str,
    *,
    num_switch: int = 100,
    short_range: int = 20000,
    include_gaps: bool = True,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build phase blocks from hap-mer markers along an assembly and count switch errors (``trio/phase_block.sh`` and ``trio/switch_error.sh``).

    ``phase_block.sh`` finds the gaps, looks up the hap-mers of both haplotypes along the assembly (BED and WIG
    tracks) and builds blocks with its fixed default (100 switches within 20,000 bp). When ``num_switch`` and
    ``short_range`` differ from that, ``switch_error.sh`` is run again on ``<prefix>.sort.bed`` with your
    values (``include_gaps=False`` restricts blocks to within contigs). Needs Java. Result keys:
    ``block_files`` (``*.phased_block.bed`` per setting), ``switch_error`` (rows of ``*.switches.txt``), ``stats_files``.
    """
    fasta = _fasta(assembly_fasta, "assembly_fasta")
    first = _meryl_named(hap1_meryl, "hap1_meryl")
    second = _meryl_named(hap2_meryl, "hap2_meryl")
    prefix = _name(output_prefix, "output_prefix")
    switches = _cli.integer(num_switch, "num_switch", minimum=1)
    ranges = _cli.integer(short_range, "short_range", minimum=1)
    if not isinstance(include_gaps, bool):
        raise ConfigurationError("include_gaps must be True or False.")
    home = _home(merqury_home)
    _check_programs(home, ["meryl", "meryl-lookup", "java", "bedtools", "awk", "bash"])
    directory = _new_directory(output_dir)
    names = _stage(directory, fasta, first, second)
    default_blocks = directory / f"{prefix}.100_20000.phased_block.bed"
    rendered = _run(
        "trio/phase_block.sh",
        [*names, prefix],
        home=home,
        cwd=directory,
        expected=[directory / f"{prefix}.sort.bed", default_blocks],
        needs=["meryl", "meryl-lookup", "java", "bedtools", "awk"],
    )
    text = rendered["_text"]
    if rendered["status"] == "completed" and (switches, ranges, include_gaps) != (100, 20000, True):
        custom = directory / f"{prefix}.{switches}_{ranges}.phased_block.bed"
        step = _run(
            "trio/switch_error.sh",
            [
                f"{prefix}.sort.bed",
                prefix,
                str(switches),
                str(ranges),
                *([] if include_gaps else ["F"]),
            ],
            home=home,
            cwd=directory,
            expected=[custom],
            needs=["java", "awk"],
        )
        rendered["outputs"] += step["outputs"]
        # switch_error.sh ends with a blob plot (R) that can fail on very small inputs; the tables are what counts
        switches_file = directory / f"{prefix}.{switches}_{ranges}.switches.txt"
        done = custom.is_file() and switches_file.is_file()
        rendered["status"] = "completed" if done else "failed"
        rendered["return_code"] = step["return_code"]
        if done and step["return_code"] != 0:
            rendered["plot_warning"] = "the trailing blob plot of switch_error.sh failed; see log_tail"
        rendered["log_tail"] = (rendered["log_tail"] + step["log_tail"])[-_cli.LOG_TAIL_CHARS :]
        text += step["_text"]
    rendered["_text"] = text

    def parse(_: str) -> dict[str, Any]:
        blocks = sorted(str(path) for path in directory.glob("*.phased_block.bed"))
        switch_rows = []
        for path in sorted(directory.glob("*.switches.txt")):
            switch_rows.append(path.read_text(encoding="utf-8").strip())
        stats = sorted(str(path) for path in directory.glob("*.phased_block.stats"))
        return {"block_files": blocks, "switch_error": switch_rows, "stats_files": stats}

    return _finish("merqury_phase_blocks", rendered, home, parse)


def merqury_block_n_stats(
    assemblies: Sequence[str],
    block_beds: Sequence[str],
    output_dir: str,
    output_prefix: str,
    *,
    genome_size: int | None = None,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Contig and phase-block N* / NG* statistics and the N* curves (``trio/block_n_stats.sh``).

    ``assemblies`` (one or two FASTA files) and ``block_beds`` (the matching ``*.phased_block.bed`` files from
    ``merqury_phase_blocks``, in the same order). ``genome_size`` (bp) adds NG* values. Writes the size tables and
    N* plots (R) with the prefix. Result key: ``files``.
    """
    if isinstance(assemblies, str) or isinstance(block_beds, str) or not 1 <= len(assemblies) <= 2:
        raise ConfigurationError("assemblies must list one or two FASTA files.")
    if len(block_beds) != len(assemblies):
        raise ConfigurationError("block_beds needs one phased-block BED per assembly.")
    fastas = [_fasta(path, "assemblies") for path in assemblies]
    beds = [_file(path, "block_beds") for path in block_beds]
    prefix = _name(output_prefix, "output_prefix")
    size = None if genome_size is None else _cli.integer(genome_size, "genome_size", minimum=1)
    home = _home(merqury_home)
    _check_programs(home, ["java", "bedtools", "Rscript", "awk", "bash"])
    directory = _new_directory(output_dir)
    staged = _stage(directory, *[item for pair in zip(fastas, beds, strict=True) for item in pair])
    arguments = [*staged, prefix] + ([] if size is None else [str(size)])
    rendered = _run(
        "trio/block_n_stats.sh",
        arguments,
        home=home,
        cwd=directory,
        expected=[],
        needs=["java", "bedtools", "Rscript", "awk"],
    )
    files = sorted(
        str(path)
        for path in directory.iterdir()
        if not path.is_symlink() and path.name != "input.fofn"
    )
    rendered["outputs"] = files
    if rendered["return_code"] == 0 and not any(
        name.startswith(prefix) for name in (Path(f).name for f in files)
    ):
        rendered["status"] = "failed"
    return _finish("merqury_block_n_stats", rendered, home, files=files)


# ---------------------------------------------------------------------------
# util/ and plot/
# ---------------------------------------------------------------------------


def merqury_bed_to_bigwig(
    bed_path: str,
    fasta_index: str,
    output_dir: str,
    *,
    convert_pipes: bool = True,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Turn a BED file of intervals into a bigWig coverage track (``util/bed_to_bigwig.sh``).

    ``fasta_index`` is the ``.fai`` of the assembly the BED refers to. Intervals are converted to a bedGraph with
    ``bedtools genomecov``, sorted per sequence and converted with ``bedGraphToBigWig``. ``convert_pipes=True``
    (Merqury's default) replaces ``|`` by ``_`` in sequence names. The result is ``<bed>.bigwig`` in the new
    ``output_dir``. Needs bedtools and bedGraphToBigWig.
    """
    bed = _file(bed_path, "bed_path")
    fai = _file(fasta_index, "fasta_index")
    if not bed.name.endswith(".bed"):
        raise ConfigurationError("bed_path must be named *.bed.")
    if not fai.name.endswith(".fai"):
        raise ConfigurationError("fasta_index must be named *.fai.")
    home = _home(merqury_home)
    _check_programs(home, ["bedtools", "bedGraphToBigWig", "bash"])
    directory = _new_directory(output_dir)
    bed_name, fai_name = _stage(directory, bed, fai)
    target = directory / (bed_name[: -len(".bed")] + ".bigwig")
    rendered = _run(
        "util/bed_to_bigwig.sh",
        [bed_name, fai_name, *([] if convert_pipes else ["F"])],
        home=home,
        cwd=directory,
        expected=[target],
        needs=["bedtools", "bedGraphToBigWig"],
    )
    return _finish("merqury_bed_to_bigwig", rendered, home, bigwig_file=str(target))


def _plot_result(directory: Path, prefix: str, pdf: bool) -> Path:
    return directory / f"{prefix}.{'pdf' if pdf else 'png'}"


def merqury_plot_spectra_cn(
    histogram: str,
    output_dir: str,
    output_prefix: str,
    *,
    only_histogram: str | None = None,
    cutoff_file: str | None = None,
    plot_type: str = "all",
    width: float = 6.0,
    height: float = 5.0,
    xmax: int = 0,
    ymax: int = 0,
    pdf: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot k-mer copy-number spectra from a ``*.spectra-cn.hist`` table (``plot/plot_spectra_cn.R``).

    ``histogram`` has the columns ``Assembly``/``Copies``, ``kmer_multiplicity`` and ``Count`` (see
    ``merqury_kmer_histogram_table``). ``only_histogram`` is the optional ``*.only.hist`` (assembly-only k-mers) and
    ``cutoff_file`` the optional solid k-mer cutoffs. ``plot_type`` is ``line``, ``fill``, ``stack`` or ``all`` (one
    file each: ``<prefix>.ln|fl|st.png``); ``xmax``/``ymax`` of 0 let the script choose. ``pdf=True`` writes PDF.
    Needs Rscript with ggplot2, argparse and scales. Result key: ``plots``.
    """
    table = _file(histogram, "histogram")
    zero = None if only_histogram is None else _file(only_histogram, "only_histogram")
    cutoff = None if cutoff_file is None else _file(cutoff_file, "cutoff_file")
    kind = _cli.choice(plot_type, ("line", "fill", "stack", "all"), "plot_type")
    prefix = _name(output_prefix, "output_prefix")
    home = _home(merqury_home)
    _check_programs(home, ["Rscript"])
    directory = _new_directory(output_dir)
    staged = _stage(directory, table, *(item for item in (zero, cutoff) if item is not None))
    arguments = [
        "-f",
        staged[0],
        "-o",
        prefix,
        "-t",
        kind,
        "-x",
        format(_cli.number(width, "width", minimum=0.1), ".6g"),
        "-y",
        format(_cli.number(height, "height", minimum=0.1), ".6g"),
        "-m",
        str(_cli.integer(xmax, "xmax", minimum=0)),
        "-n",
        str(_cli.integer(ymax, "ymax", minimum=0)),
    ]
    index = 1
    if zero is not None:
        arguments += ["-z", staged[index]]
        index += 1
    if cutoff is not None:
        arguments += ["-l", staged[index]]
    if pdf:
        arguments.append("-p")
    suffix = "pdf" if pdf else "png"
    codes = {"line": ["ln"], "fill": ["fl"], "stack": ["st"], "all": ["ln", "fl", "st"]}[kind]
    expected = [directory / f"{prefix}.{code}.{suffix}" for code in codes]
    rendered = _run(
        "plot/plot_spectra_cn.R",
        arguments,
        home=home,
        cwd=directory,
        expected=expected,
        needs=[],
        interpreter="Rscript",
    )
    return _finish(
        "merqury_plot_spectra_cn", rendered, home, plots=[str(path) for path in expected]
    )


def merqury_plot_blob(
    count_table: str,
    output_dir: str,
    output_prefix: str,
    *,
    width: float = 6.5,
    height: float = 6.0,
    pdf: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot hap-mer (or phased-block) counts per contig as a blob plot (``plot/plot_blob.R``).

    ``count_table`` is a tab-separated table with a header and the columns category, sequence id, hap1 count,
    hap2 count and sequence size (the ``*.hapmers.count`` of ``merqury_hapmer_counts``). Writes ``<prefix>.png`` (or
    ``.pdf``). Needs Rscript with ggplot2, argparse and scales. Result key: ``plot``.
    """
    table = _file(count_table, "count_table")
    prefix = _name(output_prefix, "output_prefix")
    home = _home(merqury_home)
    _check_programs(home, ["Rscript"])
    directory = _new_directory(output_dir)
    (staged,) = _stage(directory, table)
    arguments = [
        "-f",
        staged,
        "-o",
        prefix,
        "-x",
        format(_cli.number(width, "width", minimum=0.1), ".6g"),
        "-y",
        format(_cli.number(height, "height", minimum=0.1), ".6g"),
    ] + (["-p"] if pdf else [])
    target = _plot_result(directory, prefix, pdf)
    rendered = _run(
        "plot/plot_blob.R",
        arguments,
        home=home,
        cwd=directory,
        expected=[target],
        needs=[],
        interpreter="Rscript",
    )
    return _finish("merqury_plot_blob", rendered, home, plot=str(target))


def merqury_plot_block_n(
    output_dir: str,
    output_prefix: str,
    *,
    block_sizes: str | None = None,
    scaffold_sizes: str | None = None,
    contig_sizes: str | None = None,
    genome_size: int = 0,
    width: float = 6.0,
    height: float = 5.0,
    pdf: bool = False,
    merqury_home: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot N* (or NG* with ``genome_size``) continuity curves of phase blocks, scaffolds and contigs (``plot/plot_block_N.R``).

    Give at least one of ``block_sizes``, ``scaffold_sizes`` and ``contig_sizes``: ``.sizes`` tables with three tab-separated columns
    (type, group, size; longest first). Writes ``<prefix>.block.N.png`` and ``<prefix>.continuity.N.png`` (``NG`` instead of ``N`` with
    ``genome_size``; ``.pdf`` with ``pdf=True``) in the new ``output_dir``. Needs Rscript with ggplot2, argparse and scales.
    Result key: ``plots``.
    """
    sources = {"-b": block_sizes, "-s": scaffold_sizes, "-c": contig_sizes}
    if all(value is None for value in sources.values()):
        raise ConfigurationError(
            "Give at least one of block_sizes, scaffold_sizes and contig_sizes."
        )
    paths = {flag: _file(value, flag) for flag, value in sources.items() if value is not None}
    prefix = _name(output_prefix, "output_prefix")
    home = _home(merqury_home)
    _check_programs(home, ["Rscript"])
    directory = _new_directory(output_dir)
    staged = _stage(directory, *paths.values())
    arguments: list[str] = []
    for flag, name in zip(paths, staged, strict=True):
        arguments += [flag, name]
    arguments += [
        "-o",
        prefix,
        "-g",
        str(_cli.integer(genome_size, "genome_size", minimum=0)),
        "-x",
        format(_cli.number(width, "width", minimum=0.1), ".6g"),
        "-y",
        format(_cli.number(height, "height", minimum=0.1), ".6g"),
    ]
    if pdf:
        arguments.append("-p")
    suffix = "pdf" if pdf else "png"
    statistic = "NG" if genome_size else "N"
    plots = [directory / f"{prefix}.{kind}.{statistic}.{suffix}" for kind in ("block", "continuity")]
    rendered = _run(
        "plot/plot_block_N.R", arguments, home=home, cwd=directory, expected=plots, needs=[], interpreter="Rscript",
    )
    return _finish("merqury_plot_block_n", rendered, home, plots=[str(path) for path in plots])


__all__ = [
    "merqury_assign_reads",
    "merqury_bed_to_bigwig",
    "merqury_best_k",
    "merqury_block_n_stats",
    "merqury_combine_kmer_dbs",
    "merqury_exclude_reads",
    "merqury_false_duplications",
    "merqury_false_duplications_track",
    "merqury_filter_solid_kmers",
    "merqury_hamming_error",
    "merqury_hapmer_counts",
    "merqury_hapmer_spectra",
    "merqury_hapmers",
    "merqury_kmer_histogram_table",
    "merqury_kmer_multiplicity_track",
    "merqury_mash_qv",
    "merqury_per_sequence_qv",
    "merqury_phase_blocks",
    "merqury_plot_blob",
    "merqury_plot_block_n",
    "merqury_plot_spectra_cn",
    "merqury_split_10x_fastq",
]

_register_public_verbose(globals(), tuple(__all__))
