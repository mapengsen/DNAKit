"""MEGAHIT helper-program facades: read statistics, FASTG conversion and local assembly of contigs with reads.

MEGAHIT (1.2.x, tested against 1.2.9) ships ``megahit_core`` (also installed as ``megahit_toolkit``), the program behind the
``megahit`` driver. The user installs MEGAHIT (``conda install -c bioconda megahit`` or the release archive) and puts
``megahit_core`` on PATH or passes ``megahit_core_path``. Each function runs it with a fixed argument whitelist and no shell
and returns a JSON-serialisable dictionary (status, backend, version, command, outputs, parsed result, log tail). Output files
must not exist yet.

The full assembler is ``dnakit.assembly.external.assemble_megahit``. The remaining sub-programs of ``megahit_core`` (``count``,
``read2sdbg``, ``seq2sdbg``, ``iterate``, ``assemble``) are stages of that pipeline that exchange MEGAHIT-private binary files and
are therefore not exposed on their own.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="megahit_core",
    executable="megahit_core",
    hint=(
        "Install MEGAHIT 1.2.9 or later (for example `conda install -c bioconda megahit`, or the release archive from "
        "https://github.com/voutcn/megahit/releases) and put megahit_core on PATH or pass megahit_core_path."
    ),
    tested_version="1.2.9",
    citation="https://doi.org/10.1093/bioinformatics/btv033",
    version_arguments=("dumpversion",),
    version_pattern=r"(\d+\.\d+\.\d+)",
)
_STATS = re.compile(r"^(number reads|total size|longest|shortest|avg):\s*(\S+)\s*$", re.MULTILINE)


def megahit_read_statistics(
    reads: str,
    stats_path: str,
    *,
    megahit_core_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read count, total bases, longest, shortest and mean read length of a FASTA/FASTQ file (``megahit_core readstat``).

    ``reads`` may be gzip-compressed or not (MEGAHIT's parser reads the stream). The table printed by MEGAHIT is written to
    the new file ``stats_path`` and returned as ``read_count``, ``total_bases``, ``longest``, ``shortest`` and ``mean_length``.
    """
    source = _cli.existing_file(reads, "reads")
    target = _cli.new_file(stats_path, "stats_path")
    record = _cli.invoke(
        _TOOL,
        ["readstat"],
        cwd=target.parent,
        stdout_to=target,
        stdin_from=source,
        explicit=megahit_core_path,
        timeout_seconds=timeout_seconds,
    )
    found = dict(_STATS.findall(target.read_text(encoding="utf-8", errors="replace")))
    if len(found) != 5:
        raise ConfigurationError("megahit_core readstat did not print the expected table.")
    return _cli.finish(
        record,
        "read_set_statistics",
        stats_file=str(target),
        read_count=int(found["number reads"]),
        total_bases=int(found["total size"]),
        longest=int(found["longest"]),
        shortest=int(found["shortest"]),
        mean_length=float(found["avg"]),
    )


def megahit_contigs_to_fastg(
    contigs: str,
    k: int,
    fastg_path: str,
    *,
    megahit_core_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert MEGAHIT contigs (``k<k>.contigs.fa``) into a FASTG assembly graph (``megahit_core contig2fastg``).

    ``contigs`` is an intermediate contig file of a MEGAHIT run (for example ``intermediate_contigs/k21.contigs.fa``, whose
    headers carry ``flag=`` and ``multi=``) and ``k`` the k-mer size that produced it. The graph, in which the nodes are the
    contigs and the edges the k-1 overlaps between them (``NODE_1_length_..._cov_..._ID_1;`` records), is written to the new
    file ``fastg_path``. Result keys: ``node_count`` (records of the FASTG file).
    """
    source = _cli.existing_file(contigs, "contigs")
    size = _cli.integer(k, "k", minimum=3, maximum=255)
    target = _cli.new_file(fastg_path, "fastg_path")
    record = _cli.invoke(
        _TOOL,
        ["contig2fastg", str(size), str(source)],
        cwd=target.parent,
        stdout_to=target,
        explicit=megahit_core_path,
        timeout_seconds=timeout_seconds,
    )
    nodes = sum(
        1 for line in target.open(encoding="utf-8", errors="replace") if line.startswith(">")
    )
    return _cli.finish(record, "contigs_to_fastg", fastg_file=str(target), k=size, node_count=nodes)


def _prepared_contigs(source: Path, destination: Path) -> tuple[int, int]:
    """Copy contigs with the ``flag=`` and ``multi=`` header fields and write the ``.info`` file that MEGAHIT reads."""
    count = bases = 0
    name: str | None = None
    parts: list[str] = []
    entries: list[tuple[str, str]] = []
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(parts)))
                name, parts = line[1:], []
            elif name is not None:
                parts.append(line.strip().upper())
    if name is not None:
        entries.append((name, "".join(parts)))
    with destination.open("x", encoding="utf-8") as out:
        for header, sequence in entries:
            if not sequence:
                continue
            if "flag=" not in header or "multi=" not in header:
                header = f"{header.split()[0] if header.split() else 'contig'} flag=1 multi=1.0"
            out.write(f">{header}\n{sequence}\n")
            count += 1
            bases += len(sequence)
    destination.with_name(destination.name + ".info").write_text(
        f"{count} {bases}\n", encoding="utf-8"
    )
    return count, bases


def megahit_local_assembly(
    contigs: str,
    output_directory: str,
    *,
    single_reads: Sequence[str] = (),
    paired_reads: Sequence[Sequence[str]] = (),
    interleaved_reads: Sequence[str] = (),
    kmin: int = 11,
    kmax: int = 41,
    step: int = 6,
    seed_kmer: int = 31,
    min_contig_len: int = 200,
    min_mapping_len: int = 75,
    sparsity: int = 8,
    similarity: float = 0.8,
    threads: int = 1,
    megahit_core_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extend contigs or fill gaps between them by local assembly of the reads that map to their ends (``megahit_core local``).

    The reads are first converted into MEGAHIT's binary read library (``megahit_core buildlib``, which this call includes) and
    mapped to the ends of ``contigs``; the reads that map are assembled locally with iterating k from ``kmin`` to ``kmax`` in
    steps of ``step`` (``seed_kmer`` seeds the alignments, ``min_mapping_len`` and ``similarity`` filter them, ``sparsity`` thins
    the hash mapper, ``min_contig_len`` is the shortest contig that is used or reported). Give reads as ``single_reads``
    (files), ``paired_reads`` (pairs ``[read1, read2]``; paired reads are what lets the assembler bridge gaps) or
    ``interleaved_reads``. Contigs without MEGAHIT's ``flag=`` and ``multi=`` header fields get ``flag=1 multi=1.0``. The new
    contigs go to ``<output_directory>/local_contigs.fa``. Result keys: ``local_contig_file``, ``local_contig_count``,
    ``local_contig_bases``, ``input_contig_count`` and ``read_library`` (MEGAHIT's per-library summary).
    """
    source = _cli.existing_file(contigs, "contigs")
    singles = _cli.existing_files(single_reads, "single_reads", minimum=0) if single_reads else []
    pairs: list[tuple[Path, Path]] = []
    for pair in paired_reads:
        if isinstance(pair, str) or len(pair) != 2:
            raise ConfigurationError("paired_reads must be a list of [read1, read2] pairs.")
        pairs.append(
            (
                _cli.existing_file(pair[0], "paired_reads"),
                _cli.existing_file(pair[1], "paired_reads"),
            )
        )
    interleaved = (
        _cli.existing_files(interleaved_reads, "interleaved_reads", minimum=0)
        if interleaved_reads
        else []
    )
    if not (singles or pairs or interleaved):
        raise ConfigurationError("give reads: single_reads, paired_reads or interleaved_reads.")
    out = _cli.new_directory(output_directory, "output_directory")
    arguments_local = [
        "--kmin",
        str(_cli.integer(kmin, "kmin", minimum=3)),
        "--kmax",
        str(_cli.integer(kmax, "kmax", minimum=5)),
        "--step",
        str(_cli.integer(step, "step", minimum=1)),
        "--seed_kmer",
        str(_cli.integer(seed_kmer, "seed_kmer", minimum=5)),
        "--min_contig_len",
        str(_cli.integer(min_contig_len, "min_contig_len", minimum=1)),
        "--min_mapping_len",
        str(_cli.integer(min_mapping_len, "min_mapping_len", minimum=1)),
        "--sparsity",
        str(_cli.integer(sparsity, "sparsity", minimum=1)),
        "--similarity",
        format(_cli.number(similarity, "similarity", minimum=0.0, maximum=1.0), ".12g"),
        "-t",
        str(_cli.integer(threads, "threads", minimum=1, maximum=256)),
    ]
    if kmin > kmax:
        raise ConfigurationError("kmin must not exceed kmax.")
    prepared = out / "input_contigs.fa"
    input_count, _ = _prepared_contigs(source, prepared)
    library = out / "reads.lib"
    lines: list[str] = []
    for path in interleaved:
        lines += [str(path), f"interleaved {path}"]
    for first, second in pairs:
        lines += [f"{first},{second}", f"pe {first} {second}"]
    for path in singles:
        lines += [str(path), f"se {path}"]
    library.write_text("\n".join(lines) + "\n", encoding="utf-8")
    built = _cli.invoke(
        _TOOL,
        ["buildlib", str(library), str(library)],
        cwd=out,
        expected=[out / "reads.lib.bin", out / "reads.lib.lib_info"],
        explicit=megahit_core_path,
        timeout_seconds=timeout_seconds,
    )
    target = out / "local_contigs.fa"
    record = _cli.invoke(
        _TOOL,
        ["local", "-c", str(prepared), "-l", str(library), "-o", str(target), *arguments_local],
        cwd=out,
        expected=[target],
        explicit=megahit_core_path,
        timeout_seconds=timeout_seconds,
    )
    count = bases = 0
    for line in target.open(encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            count += 1
        else:
            bases += len(line.strip())
    summary = [
        line.split(" - ", 1)[-1] for line in built["log_tail"].splitlines() if "Lib " in line
    ]
    record["command"] = built["command"] + ["&&"] + record["command"]
    record["outputs"] = [str(target)]
    return _cli.finish(
        record,
        "local_assembly_of_contigs_with_reads",
        local_contig_file=str(target),
        local_contig_count=count,
        local_contig_bases=bases,
        input_contig_count=input_count,
        read_library=summary,
    )


__all__ = ["megahit_contigs_to_fastg", "megahit_local_assembly", "megahit_read_statistics"]
_register_public_verbose(globals())
