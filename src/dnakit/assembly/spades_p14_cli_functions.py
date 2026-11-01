"""SPAdes toolkit command-line facades (tested with SPAdes 4.3.0): read error correction (BayesHammer and
IonHammer), k-mer coverage read filtering, k-mer count estimation, assembly-graph building, simplification and
splitting, read and long-read alignment to a graph (gmapper, SPAligner), profile-HMM search in a graph (PathRacer),
repeat resolution with synthetic long reads (SPlitteR) and metagenome bin refinement (BinSPreader).

The user installs SPAdes (``conda install -c bioconda spades``, or the release archive from
https://github.com/ablab/spades) and puts it on PATH, or passes ``spades_path`` (the folder that holds the SPAdes
programs). Each public function runs one fixed SPAdes program with an argument list built from typed, validated
parameters: no shell, no generic runner, outputs are new files or folders. Read sets are given as file lists; the
YAML dataset description the programs need is written for you in a temporary folder. Only the options that the
caller sets are passed, so every other option keeps SPAdes' own default. Results are JSON-serialisable dictionaries
(status, backend, version, command, return code, outputs, log tail).
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = (
    "Install SPAdes 4.0 or later (for example `conda install -c bioconda spades`, or the archive from "
    "https://github.com/ablab/spades), then put it on PATH or pass spades_path (the folder that contains the "
    "SPAdes programs)."
)
_KMER_COUNT = re.compile(r"Kmer number estimation:\s*(\d+)")
_ORIENTATIONS = ("fr", "rf", "ff")
_GRAPH_FORMATS = {"unitigs": "--unitigs", "fastg": "--fastg", "gfa": "--gfa", "spades": "--spades"}


def _executable(name: str, spades_path: str | None) -> str | None:
    if spades_path is None:
        return None
    folder = Path(spades_path).expanduser()
    if not folder.is_dir():
        raise ConfigurationError(f"spades_path must be a folder that contains the SPAdes programs: {folder}")
    return str(folder / name)


def _run(
    program: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    spades_path: str | None,
    timeout_seconds: float = 86_400.0,
    version_args: Sequence[str] = ("--version",),
) -> dict[str, Any]:
    result = _p14_cli.run(
        program,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=_executable(program, spades_path),
        hint=_HINT,
        version_args=("--version",) if program == "spades.py" else ("--no-such-flag",),
        timeout_seconds=_gap_cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0),
    )
    if result["version"] is None:
        try:
            main = _gap_cli.program("spades.py", explicit=_executable("spades.py", spades_path), hint=_HINT)
            result["version"] = _p14_cli.version_of(main, ("--version",))
        except Exception:  # noqa: BLE001 - the version is informational only
            pass
    return result


def _reads(reads_1: Sequence[str] | None, reads_2: Sequence[str] | None, unpaired: Sequence[str] | None = None) -> tuple[list[str], list[str], list[str]]:
    def files(values: Sequence[str] | None, name: str) -> list[str]:
        if values is None:
            return []
        if isinstance(values, str) or not values:
            raise ConfigurationError(f"{name} must be a non-empty list of read files.")
        paths = [str(_gap_cli.existing_file(item, name)) for item in values]
        if any(re.search(r"[\"\n]", item) for item in paths):
            raise ConfigurationError(f"{name} paths must not contain quotes or newlines.")
        return paths

    first, second, single = files(reads_1, "reads_1"), files(reads_2, "reads_2"), files(unpaired, "unpaired_reads")
    if len(first) != len(second):
        raise ConfigurationError("reads_1 and reads_2 must list the same number of files.")
    if not first and not single:
        raise ConfigurationError("Give reads_1 and reads_2, or unpaired_reads.")
    return first, second, single


def _dataset_yaml(folder: Path, first: list[str], second: list[str], single: list[str], orientation: str = "fr") -> Path:
    """Write the SPAdes dataset description (YAML) for the given read files."""
    quoted = lambda items: ", ".join(f'"{item}"' for item in items)  # noqa: E731
    libraries = []
    if first:
        libraries.append(
            f'  {{\n    orientation: "{_gap_cli.choice(orientation, _ORIENTATIONS, "orientation")}",\n    type: "paired-end",\n'
            f"    right reads: [{quoted(second)}],\n    left reads: [{quoted(first)}]\n  }}"
        )
    if single:
        libraries.append(f'  {{\n    type: "single",\n    single reads: [{quoted(single)}]\n  }}')
    path = folder / "dataset.yaml"
    path.write_text("[\n" + ",\n".join(libraries) + "\n]\n", encoding="utf-8")
    return path


def _threads(value: int | None) -> list[str]:
    return [] if value is None else ["-t", str(_gap_cli.integer(value, "threads", minimum=1, maximum=256))]


def spades_correct_reads(
    output_dir: str,
    *,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    iontorrent: bool = False,
    threads: int | None = None,
    memory_gb: int | None = None,
    phred_offset: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Correct sequencing errors in short reads without assembling (``spades.py --only-error-correction``).

    Illumina-type reads are corrected with BayesHammer (``spades-hammer``); with ``iontorrent=True`` the Ion Torrent
    corrector IonHammer (``spades-ionhammer``) is used. Give paired files ``reads_1``/``reads_2`` and/or
    ``unpaired_reads``. ``memory_gb`` is the RAM limit (SPAdes default 250) and ``phred_offset`` 33 or 64. Writes the
    corrected FASTQ files to ``<output_dir>/corrected`` (``output_dir`` must be new).
    """
    first, second, single = _reads(reads_1, reads_2, unpaired_reads)
    folder = _gap_cli.new_directory(output_dir, "output_dir")
    arguments = ["--only-error-correction", "-o", str(folder)]
    if len(first) > 1:
        raise ConfigurationError("Give one pair of paired-read files; concatenate lanes first.")
    if first:
        arguments += ["--pe1-1", first[0], "--pe1-2", second[0]]
    for item in single:
        arguments += ["-s", item]
    if _p14_cli.flag(iontorrent, "iontorrent"):
        arguments.append("--iontorrent")
    if threads is not None:
        arguments += ["-t", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if memory_gb is not None:
        arguments += ["-m", str(_gap_cli.integer(memory_gb, "memory_gb", minimum=1))]
    if phred_offset is not None:
        arguments += ["--phred-offset", str(_gap_cli.choice(str(phred_offset), ("33", "64"), "phred_offset"))]
    result = _run("spades.py", arguments, cwd=folder, expected=[], spades_path=spades_path)
    corrected = sorted(str(path) for path in (folder / "corrected").glob("*") if path.is_file()) if (folder / "corrected").is_dir() else []
    result["corrected_files"] = corrected
    result["corrector"] = "ionhammer" if iontorrent else "bayeshammer"
    result["outputs"] = corrected
    if result["status"] == "completed" and not any(item.endswith((".fastq.gz", ".fastq")) for item in corrected):
        result["status"] = "failed"
    return result


def spades_filter_reads_by_kmer_coverage(
    output_dir: str,
    *,
    kmer: int,
    coverage_threshold: int,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    threads: int | None = None,
    drop_names: bool = False,
    drop_quality: bool = False,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove reads whose median k-mer count is low (``spades-read-filter``).

    For read pairs both reads must have a median ``kmer``-count above ``coverage_threshold`` to be kept; for
    unpaired reads the read itself. ``drop_names`` and ``drop_quality`` write FASTQ without read names/qualities
    (faster). ``output_dir`` must be new; the filtered FASTQ files and ``dataset.yaml`` are written there.
    """
    first, second, single = _reads(reads_1, reads_2, unpaired_reads)
    folder_out = Path(output_dir).expanduser().resolve()
    if folder_out.exists() or not folder_out.parent.is_dir():
        raise ConfigurationError(f"output_dir must be a new folder inside an existing folder: {folder_out}")
    with tempfile.TemporaryDirectory(prefix="dnakit_spades_") as temporary:
        yaml = _dataset_yaml(Path(temporary), first, second, single)
        arguments = ["-k", str(_gap_cli.integer(kmer, "kmer", minimum=1, maximum=255)), "-c", str(_gap_cli.integer(coverage_threshold, "coverage_threshold", minimum=0)), "-d", str(yaml), "-o", str(folder_out)]
        arguments += _threads(threads)
        for value, flag_name, name in ((drop_names, "--drop-names", "drop_names"), (drop_quality, "--drop-quality", "drop_quality")):
            if _p14_cli.flag(value, name):
                arguments.append(flag_name)
        result = _run("spades-read-filter", arguments, cwd=Path(temporary), expected=[], spades_path=spades_path, version_args=("--help",))
    files = sorted(str(path) for path in folder_out.glob("*") if path.is_file()) if folder_out.is_dir() else []
    result["outputs"] = files
    if result["status"] == "completed" and not any(item.endswith("dataset.yaml") for item in files):
        result["status"] = "failed"
    return result


def spades_estimate_kmer_count(
    kmer: int,
    *,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the number of distinct ``kmer``-mers in reads (``spades-kmer-estimating``).

    K-mers from reverse-complementary reads are not counted separately. ``distinct_kmers`` holds the estimate.
    """
    first, second, single = _reads(reads_1, reads_2, unpaired_reads)
    with tempfile.TemporaryDirectory(prefix="dnakit_spades_") as temporary:
        yaml = _dataset_yaml(Path(temporary), first, second, single)
        arguments = ["-k", str(_gap_cli.integer(kmer, "kmer", minimum=1, maximum=255)), "-d", str(yaml)] + _threads(threads)
        result = _run("spades-kmer-estimating", arguments, cwd=Path(temporary), expected=[], spades_path=spades_path, version_args=("--help",))
    match = _KMER_COUNT.search(result["log_tail"])
    result["distinct_kmers"] = int(match.group(1)) if match is not None else None
    result["outputs"] = []
    if result["status"] == "completed" and match is None:
        result["status"] = "failed"
    return result


def spades_build_assembly_graph(
    output: str,
    *,
    kmer: int,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    fasta: str | None = None,
    graph_format: str = "unitigs",
    infer_coverage: bool = False,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a de Bruijn assembly graph from reads or a FASTA file (``spades-gbuilder``).

    Give reads (``reads_1``/``reads_2``/``unpaired_reads``) or ``fasta`` (sequences to build the graph from).
    ``graph_format`` is 'unitigs' (FASTA of unitigs, SPAdes default), 'gfa' (GFA1), 'fastg' or 'spades' (SPAdes
    internal format); ``infer_coverage`` estimates edge coverage. ``kmer`` is the k-mer length (odd). Writes a new
    file.
    """
    output_path = _gap_cli.new_file(output, "output")
    chosen = _gap_cli.choice(graph_format, tuple(_GRAPH_FORMATS), "graph_format")
    with tempfile.TemporaryDirectory(prefix="dnakit_spades_") as temporary:
        if fasta is not None:
            if reads_1 is not None or reads_2 is not None or unpaired_reads is not None:
                raise ConfigurationError("Give reads or fasta, not both.")
            source = str(_gap_cli.existing_file(fasta, "fasta"))
        else:
            first, second, single = _reads(reads_1, reads_2, unpaired_reads)
            source = str(_dataset_yaml(Path(temporary), first, second, single))
        arguments = [source, str(output_path), "-k", str(_gap_cli.integer(kmer, "kmer", minimum=3, maximum=255)), _GRAPH_FORMATS[chosen]]
        if _p14_cli.flag(infer_coverage, "infer_coverage"):
            arguments.append("-c")
        arguments += _threads(threads)
        arguments += ["-tmp-dir", temporary]
        return _run("spades-gbuilder", arguments, cwd=output_path.parent, expected=[output_path], spades_path=spades_path, version_args=("--help",))


def spades_simplify_assembly_graph(
    graph: str,
    output_prefix: str,
    *,
    kmer: int,
    read_length: int,
    coverage: str | float | None = None,
    use_coverage_ratios: bool = False,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simplify an assembly graph (remove tips, bulges and low-coverage edges) with ``spades-gsimplifier``.

    ``graph`` is a GFA file; ``kmer`` and ``read_length`` describe the data it came from. ``coverage`` is the
    estimated (k+1)-mer coverage (a number, or 'auto'; default 0). ``use_coverage_ratios`` enables the procedures that
    use unitig coverage ratios. Writes ``<output_prefix>.gfa`` (the prefix must be new).
    """
    source = _gap_cli.existing_file(graph, "graph")
    prefix = Path(output_prefix).expanduser().resolve()
    target = Path(str(prefix) + ".gfa")
    if not prefix.parent.is_dir() or target.exists():
        raise ConfigurationError(f"output_prefix must be a new prefix inside an existing folder: {prefix}")
    with tempfile.TemporaryDirectory(prefix="dnakit_spades_", dir=prefix.parent) as temporary:
        arguments = [str(source), str(prefix), "-k", str(_gap_cli.integer(kmer, "kmer", minimum=3, maximum=255)), "--read-length", str(_gap_cli.integer(read_length, "read_length", minimum=1)), "--tmpdir", temporary]
        if coverage is not None:
            arguments += ["-c", "auto" if coverage == "auto" else str(_gap_cli.number(coverage, "coverage", minimum=0.0))]
        if _p14_cli.flag(use_coverage_ratios, "use_coverage_ratios"):
            arguments.append("--use-cov-ratios")
        arguments += _threads(threads)
        return _run("spades-gsimplifier", arguments, cwd=prefix.parent, expected=[target], spades_path=spades_path, version_args=("--help",))


def spades_split_assembly_graph(
    graph: str,
    output_dir: str,
    *,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split an assembly graph (GFA) into one GFA file per connected component (``spades-gfa-split``).

    Writes ``subgraph_<n>.gfa`` files into ``output_dir`` (a new folder that the program creates); ``components`` is
    their number.
    """
    source = _gap_cli.existing_file(graph, "graph")
    folder = Path(output_dir).expanduser().resolve()
    if folder.exists() or not folder.parent.is_dir():
        raise ConfigurationError(f"output_dir must be a new folder inside an existing folder: {folder}")
    result = _run("spades-gfa-split", [str(source), str(folder)], cwd=folder.parent, expected=[], spades_path=spades_path, version_args=("--help",))
    produced = sorted(str(path) for path in folder.glob("*.gfa")) if folder.is_dir() else []
    result["outputs"] = produced
    result["components"] = len(produced)
    if result["status"] == "completed" and not produced:
        result["status"] = "failed"
    return result


def spades_map_reads_to_graph(
    graph: str,
    output: str,
    *,
    kmer: int,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    library_index: int | None = None,
    mode: str | None = None,
    retain: str | None = None,
    hic: bool = False,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align reads to an assembly graph and write the paired-end distance histogram (``spades-gmapper``).

    ``graph`` is a GFA file and ``kmer`` its k-mer size. ``library_index`` selects the library of the dataset (0-based,
    default 0; paired libraries come first). ``mode`` is the inner alignment mode 'illumina', 'pacbio', 'intractg' or
    'hic'; ``retain`` keeps 'default', 'all', 'primary' or 'quality' alignments; ``hic`` enables Hi-C-aware pair
    handling. Writes a new file.
    """
    source = _gap_cli.existing_file(graph, "graph")
    output_path = _gap_cli.new_file(output, "output")
    first, second, single = _reads(reads_1, reads_2, unpaired_reads)
    with tempfile.TemporaryDirectory(prefix="dnakit_spades_") as temporary:
        yaml = _dataset_yaml(Path(temporary), first, second, single)
        arguments = [str(yaml), str(source), str(output_path), "-k", str(_gap_cli.integer(kmer, "kmer", minimum=3, maximum=255))]
        if library_index is not None:
            arguments += ["-l", str(_gap_cli.integer(library_index, "library_index", minimum=0))]
        if mode is not None:
            arguments.append("-X" + _gap_cli.choice(mode, ("illumina", "pacbio", "intractg", "hic"), "mode"))
        if retain is not None:
            arguments.append("-R" + _gap_cli.choice(retain, ("default", "all", "primary", "quality"), "retain"))
        if _p14_cli.flag(hic, "hic"):
            arguments.append("--hic")
        arguments += _threads(threads) + ["--tmp-dir", temporary]
        return _run("spades-gmapper", arguments, cwd=output_path.parent, expected=[output_path], spades_path=spades_path, version_args=("--help",))


def spaligner_align_long_reads(
    sequences: str,
    graph: str,
    output_dir: str,
    *,
    datatype: str,
    kmer: int,
    config: str | None = None,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align long reads or contigs to an assembly graph with SPAligner (``spaligner``).

    ``datatype`` is 'nanopore' or 'pacbio'; ``graph`` a GFA file (or SPAdes saves folder) built with k-mer size
    ``kmer`` (odd); ``config`` is the SPAligner YAML parameter file (default: the ``spaligner_config.yaml`` installed
    with SPAligner). Writes ``<output_dir>/alignment.tsv`` (columns: read, read start/end/length, path, ...);
    ``output_dir`` must be new.
    """
    reads = _gap_cli.existing_file(sequences, "sequences")
    graph_path = _gap_cli.existing_file(graph, "graph")
    folder = Path(output_dir).expanduser().resolve()
    if folder.exists() or not folder.parent.is_dir():
        raise ConfigurationError(f"output_dir must be a new folder inside an existing folder: {folder}")
    if config is not None:
        config_path = _gap_cli.existing_file(config, "config")
    else:
        executable = Path(_gap_cli.program("spaligner", explicit=_executable("spaligner", spades_path), hint=_HINT))
        candidates = [executable.parent.parent / "share" / "spaligner" / "spaligner_config.yaml", executable.parent / "spaligner_config.yaml"]
        found = next((item for item in candidates if item.is_file()), None)
        if found is None:
            raise ConfigurationError("spaligner_config.yaml was not found next to the installation; pass config.")
        config_path = found
    arguments = [
        str(config_path), "-d", _gap_cli.choice(datatype, ("nanopore", "pacbio"), "datatype"), "-s", str(reads),
        "-g", str(graph_path), "-k", str(_gap_cli.integer(kmer, "kmer", minimum=3, maximum=255)), "-o", str(folder),
    ] + _threads(threads)
    parent = folder.parent
    result = _run("spaligner", arguments, cwd=parent, expected=[], spades_path=spades_path, version_args=("--help",))
    table = folder / "alignment.tsv"
    result["outputs"] = [str(table)] if table.is_file() else []
    if result["status"] == "completed" and not table.is_file():
        result["status"] = "failed"
    return result


def pathracer_search_graph(
    query: str,
    graph: str,
    output_dir: str,
    *,
    query_type: str = "hmm",
    mode: str = "global",
    minimum_length: float | None = None,
    indel_rate: float | None = None,
    top: int | None = None,
    threads: int | None = None,
    memory_gb: int | None = None,
    seeding: str | None = None,
    evalue: float | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find paths in an assembly graph that match a profile HMM or a sequence (``pathracer``).

    ``query`` is an HMMER3 profile file (``query_type='hmm'``), or a FASTA file of nucleotide (``'nt'``) or amino-acid
    (``'aa'``) strings. ``graph`` is a GFA file or SPAdes saves folder. ``mode`` is 'global' (glocal, default) or
    'local'; ``minimum_length`` is the minimal matched length (a fraction of the HMM length if <= 1; default 0.9);
    ``indel_rate`` the expected nucleotide indel rate on edges (for amino-acid HMMs with frame shifts); ``top`` the
    number of paths extracted (default 10000); ``seeding`` one of 'edges', 'scaffolds', 'edges-scaffolds',
    'exhaustive', 'edges-1-by-1', 'scaffolds-1-by-1'; ``evalue`` the HMMER reporting E-value. Writes the result files
    (matched sequences, annotated graph, log) into the new folder ``output_dir``.
    """
    query_path = _gap_cli.existing_file(query, "query")
    graph_path = _gap_cli.existing_file(graph, "graph")
    folder = _gap_cli.new_directory(output_dir, "output_dir")
    arguments = [str(query_path), str(graph_path), "--output", str(folder)]
    arguments.append("--" + _gap_cli.choice(query_type, ("hmm", "nt", "aa"), "query_type"))
    arguments.append("--" + _gap_cli.choice(mode, ("global", "local"), "mode"))
    if minimum_length is not None:
        arguments += ["--length", str(_gap_cli.number(minimum_length, "minimum_length", minimum=0.0))]
    if indel_rate is not None:
        arguments += ["--indel-rate", str(_gap_cli.number(indel_rate, "indel_rate", minimum=0.0, maximum=1.0))]
    if top is not None:
        arguments += ["--top", str(_gap_cli.integer(top, "top", minimum=1))]
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if memory_gb is not None:
        arguments += ["--memory", str(_gap_cli.integer(memory_gb, "memory_gb", minimum=1))]
    if seeding is not None:
        arguments.append("--seed-" + _gap_cli.choice(seeding, ("edges", "scaffolds", "edges-scaffolds", "exhaustive", "edges-1-by-1", "scaffolds-1-by-1"), "seeding"))
    if evalue is not None:
        arguments += ["-E", str(_gap_cli.number(evalue, "evalue", minimum=0.0))]
    result = _run("pathracer", arguments, cwd=folder, expected=[folder / "pathracer.log"], spades_path=spades_path, version_args=("--help",))
    result["outputs"] = sorted(str(path) for path in folder.glob("*") if path.is_file())
    return result


def pathracer_align_hmm_to_sequences(
    hmm: str,
    sequences: str,
    output_dir: str,
    *,
    mode: str = "global",
    indel_rate: float | None = None,
    top: int | None = None,
    max_frameshifts: int | None = None,
    cutoff: float | None = None,
    exhaustive: bool = False,
    threads: int | None = None,
    memory_gb: int | None = None,
    evalue: float | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align an amino-acid profile HMM to nucleotide sequences, allowing frame shifts (``pathracer-seq-fs``).

    ``hmm`` is a HMMER3 amino-acid profile and ``sequences`` a nucleotide FASTA file. ``max_frameshifts`` is the
    largest number of frame shifts in a reported sequence (default 10), ``indel_rate`` the expected nucleotide indel
    rate (default 0.05), ``top`` the number of paths (default 100), ``cutoff`` the bit-score cutoff (a fraction of the
    HMM's GA cutoff if <= 1; default 0.7) and ``exhaustive`` disables the HMM pre-filter. Writes the result files
    into the new folder ``output_dir``.
    """
    hmm_path = _gap_cli.existing_file(hmm, "hmm")
    sequence_path = _gap_cli.existing_file(sequences, "sequences")
    folder = _gap_cli.new_directory(output_dir, "output_dir")
    arguments = [str(hmm_path), str(sequence_path), "--output", str(folder), "--" + _gap_cli.choice(mode, ("global", "local"), "mode")]
    if indel_rate is not None:
        arguments += ["--indel-rate", str(_gap_cli.number(indel_rate, "indel_rate", minimum=0.0, maximum=1.0))]
    if top is not None:
        arguments += ["--top", str(_gap_cli.integer(top, "top", minimum=1))]
    if max_frameshifts is not None:
        arguments += ["--max-fs", str(_gap_cli.integer(max_frameshifts, "max_frameshifts", minimum=0))]
    if cutoff is not None:
        arguments += ["--cutoff", str(_gap_cli.number(cutoff, "cutoff", minimum=0.0))]
    if _p14_cli.flag(exhaustive, "exhaustive"):
        arguments.append("--exhaustive")
    if threads is not None:
        arguments += ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if memory_gb is not None:
        arguments += ["--memory", str(_gap_cli.integer(memory_gb, "memory_gb", minimum=1))]
    if evalue is not None:
        arguments += ["-E", str(_gap_cli.number(evalue, "evalue", minimum=0.0))]
    result = _run("pathracer-seq-fs", arguments, cwd=folder, expected=[], spades_path=spades_path, version_args=("--help",))
    result["outputs"] = sorted(str(path) for path in folder.glob("*") if path.is_file())
    return result


def splitter_resolve_repeats(
    graph: str,
    library_yaml: str,
    output_dir: str,
    *,
    library_index: int | None = None,
    graph_type: str | None = None,
    resolution_mode: str | None = None,
    score_threshold: float | None = None,
    relative_threshold: float | None = None,
    scaffold_links: bool = False,
    mapping_kmer: int | None = None,
    sampling_factor: float | None = None,
    seed: int | None = None,
    assembly_info: str | None = None,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resolve repeats in an assembly graph with synthetic long reads / linked reads (``splitter``, SPlitteR).

    ``graph`` is a GFA (or binary) graph; ``library_yaml`` is the SPAdes YAML description of the synthetic-long-read
    library (barcoded reads). ``graph_type`` is 'mdbg' or 'blunt'; ``resolution_mode`` is 'diploid' or 'meta'
    (``assembly_info`` is the metaFlye ``assembly_info.txt`` for meta mode); ``score_threshold`` and
    ``relative_threshold`` set the link-index score thresholds, ``scaffold_links`` adds scaffold links,
    ``mapping_kmer`` the k for read mapping, ``sampling_factor`` and ``seed`` the barcode down-sampling.
    Writes the resolved graph into the new folder ``output_dir``.
    """
    graph_path = _gap_cli.existing_file(graph, "graph")
    yaml_path = _gap_cli.existing_file(library_yaml, "library_yaml")
    folder = _gap_cli.new_directory(output_dir, "output_dir")
    arguments = [str(graph_path), str(yaml_path), str(folder)]
    if library_index is not None:
        arguments += ["-l", str(_gap_cli.integer(library_index, "library_index", minimum=0))]
    arguments += _threads(threads)
    if graph_type is not None:
        arguments.append("-G" + _gap_cli.choice(graph_type, ("mdbg", "blunt"), "graph_type"))
    if resolution_mode is not None:
        arguments.append("-M" + _gap_cli.choice(resolution_mode, ("diploid", "meta"), "resolution_mode"))
    if score_threshold is not None:
        arguments += ["--score", str(_gap_cli.number(score_threshold, "score_threshold"))]
    if relative_threshold is not None:
        arguments += ["--rel-threshold", str(_gap_cli.number(relative_threshold, "relative_threshold", minimum=0.0))]
    if _p14_cli.flag(scaffold_links, "scaffold_links"):
        arguments.append("--scaffold-links")
    if mapping_kmer is not None:
        arguments += ["--mapping-k", str(_gap_cli.integer(mapping_kmer, "mapping_kmer", minimum=3))]
    if sampling_factor is not None:
        arguments += ["--sampling-factor", str(_gap_cli.number(sampling_factor, "sampling_factor", minimum=0.0))]
    if seed is not None:
        arguments += ["--seed", str(_gap_cli.integer(seed, "seed", minimum=0))]
    if assembly_info is not None:
        arguments += ["--assembly-info", str(_gap_cli.existing_file(assembly_info, "assembly_info"))]
    result = _run("splitter", arguments, cwd=folder, expected=[], spades_path=spades_path, version_args=("--help",))
    result["outputs"] = sorted(str(path) for path in folder.rglob("*") if path.is_file())
    return result


def binspreader_refine_bins(
    graph: str,
    binning: str,
    output: str,
    *,
    contig_paths: str | None = None,
    dataset_yaml: str | None = None,
    multiple_assignment: bool = False,
    assignment_strategy: str | None = None,
    refiner: str | None = None,
    cami_format: bool = False,
    zero_bin: bool = False,
    labeled_alpha: float | None = None,
    sparse_propagation: bool = False,
    split_reads: bool = False,
    read_bin_weight_threshold: float | None = None,
    threads: int | None = None,
    spades_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Refine and propagate metagenome bins over the assembly graph (``binspreader``, BinSPreader).

    ``graph`` is the assembly graph (GFA or binary), ``binning`` the initial contig-to-bin TSV from a binner, and
    ``contig_paths`` the ``contigs.paths`` file that maps contigs to graph edges. ``assignment_strategy`` is 'max' or
    'mle', ``refiner`` 'corr' (correction) or 'prop' (propagation); ``multiple_assignment`` allows a contig in
    several bins; ``cami_format`` reads/writes CAMI bioboxes binning; ``zero_bin`` writes a bin 0 for unbinned
    contigs; ``labeled_alpha`` is the label-correction alpha; ``sparse_propagation`` gradually reduces alpha; with
    ``split_reads`` (needs ``dataset_yaml``) reads are split by bin using ``read_bin_weight_threshold``. Writes the
    refined binning TSV to a new file.
    """
    graph_path = _gap_cli.existing_file(graph, "graph")
    binning_path = _gap_cli.existing_file(binning, "binning")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [str(graph_path), str(binning_path), str(output_path)]
    if contig_paths is not None:
        arguments += ["--paths", str(_gap_cli.existing_file(contig_paths, "contig_paths"))]
    if dataset_yaml is not None:
        arguments += ["--dataset", str(_gap_cli.existing_file(dataset_yaml, "dataset_yaml"))]
    arguments += _threads(threads)
    if assignment_strategy is not None:
        arguments.append("-S" + _gap_cli.choice(assignment_strategy, ("max", "mle"), "assignment_strategy"))
    if refiner is not None:
        arguments.append("-R" + _gap_cli.choice(refiner, ("corr", "prop"), "refiner"))
    if labeled_alpha is not None:
        arguments += ["-la", str(_gap_cli.number(labeled_alpha, "labeled_alpha", minimum=0.0, maximum=1.0))]
    for value, flag_name, name in (
        (multiple_assignment, "-m", "multiple_assignment"), (cami_format, "--cami", "cami_format"),
        (zero_bin, "--zero-bin", "zero_bin"), (sparse_propagation, "--sparse-propagation", "sparse_propagation"),
        (split_reads, "-r", "split_reads"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if split_reads and dataset_yaml is None:
        raise ConfigurationError("split_reads needs dataset_yaml.")
    if read_bin_weight_threshold is not None:
        arguments += ["-b", str(_gap_cli.number(read_bin_weight_threshold, "read_bin_weight_threshold", minimum=0.0, maximum=1.0))]
    return _run("binspreader", arguments, cwd=output_path.parent, expected=[output_path], spades_path=spades_path, version_args=("--help",))


__all__ = [
    "binspreader_refine_bins",
    "pathracer_align_hmm_to_sequences",
    "pathracer_search_graph",
    "spades_build_assembly_graph",
    "spades_correct_reads",
    "spades_estimate_kmer_count",
    "spades_filter_reads_by_kmer_coverage",
    "spades_map_reads_to_graph",
    "spades_simplify_assembly_graph",
    "spades_split_assembly_graph",
    "spaligner_align_long_reads",
    "splitter_resolve_repeats",
]
_register_public_verbose(globals())
