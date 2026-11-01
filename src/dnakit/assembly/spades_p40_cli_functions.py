"""SPAdes assembly modes as one DNAKit function (pack P40).

``spades_assemble`` runs ``spades.py`` with the mode flag that selects the assembler variant: single-cell
(``--sc``), metagenome (``--meta``), plasmid (``--plasmid``), metagenomic plasmid (``--metaplasmid``), metagenomic
viral (``--metaviral``), transcriptome (``--rna``), RNA viral (``--rnaviral``), coronavirus (``--corona``),
biosynthetic gene cluster (``--bio``), wastewater strain deconvolution (``--sewage``), Ion Torrent (``--iontorrent``),
or the default isolate mode. The program is installed by the user (``conda install -c bioconda spades``) and is found
on PATH, or in ``spades_path`` (the folder with the SPAdes programs). Arguments are a fixed whitelist; no shell.
Only options the caller sets are passed, so everything else keeps SPAdes' own default.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.assembly.spades_p14_cli_functions import _run
from dnakit.bioinformatics import _gap_cli
from dnakit.exceptions import ConfigurationError

_MODES = {
    "isolate": None,
    "single_cell": "--sc",
    "metagenome": "--meta",
    "plasmid": "--plasmid",
    "metaplasmid": "--metaplasmid",
    "metaviral": "--metaviral",
    "rna": "--rna",
    "rnaviral": "--rnaviral",
    "corona": "--corona",
    "biosynthetic": "--bio",
    "sewage": "--sewage",
    "iontorrent": "--iontorrent",
}
_PRIMARY = ("contigs.fasta", "scaffolds.fasta", "transcripts.fasta", "hard_filtered_transcripts.fasta",
            "soft_filtered_transcripts.fasta", "gene_clusters.fasta", "assembly_graph.fastg",
            "assembly_graph_with_scaffolds.gfa", "spades.log")


def _files(values: Sequence[str] | None, name: str, *, limit: int = 1) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of files.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} accepts at most {limit} file(s) here; concatenate lanes first.")
    return [str(_gap_cli.existing_file(item, name)) for item in values]


def spades_assemble(
    output_dir: str,
    *,
    mode: str = "isolate",
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    interleaved_reads: Sequence[str] | None = None,
    unpaired_reads: Sequence[str] | None = None,
    nanopore_reads: Sequence[str] | None = None,
    pacbio_reads: Sequence[str] | None = None,
    trusted_contigs: Sequence[str] | None = None,
    kmers: Sequence[int] | None = None,
    careful: bool = False,
    only_assembler: bool = False,
    threads: int | None = None,
    memory_gb: int | None = None,
    phred_offset: int | None = None,
    spades_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assemble reads with SPAdes in one of its modes (``spades.py`` and the ``*spades.py`` wrappers).

    ``mode`` selects the variant: ``isolate`` (default), ``single_cell`` (--sc, MDA data), ``metagenome``
    (--meta, metaSPAdes), ``plasmid`` (--plasmid, plasmidSPAdes), ``metaplasmid`` (--metaplasmid), ``metaviral``
    (--metaviral), ``rna`` (--rna, rnaSPAdes), ``rnaviral`` (--rnaviral), ``corona`` (--corona, HMM-guided
    coronaSPAdes), ``biosynthetic`` (--bio), ``sewage`` (--sewage, SARS-CoV-2 strain deconvolution) and
    ``iontorrent`` (--iontorrent). Give paired files ``reads_1``/``reads_2`` (one pair), ``interleaved_reads``
    and/or ``unpaired_reads``; ``nanopore_reads``/``pacbio_reads`` are long reads, ``trusted_contigs`` guides the
    assembly. ``careful`` adds mismatch correction (not allowed for metagenome, rna, rnaviral and biosynthetic
    modes). ``output_dir`` must be new. The result lists the assembly files found in it (contigs, scaffolds,
    transcripts, gene clusters, graphs, log).
    """

    selected = _gap_cli.choice(mode, tuple(_MODES), "mode")
    folder_arg = output_dir
    first, second = _files(reads_1, "reads_1"), _files(reads_2, "reads_2")
    inter = _files(interleaved_reads, "interleaved_reads")
    single = _files(unpaired_reads, "unpaired_reads", limit=100)
    nano = _files(nanopore_reads, "nanopore_reads", limit=100)
    pacbio = _files(pacbio_reads, "pacbio_reads", limit=100)
    trusted = _files(trusted_contigs, "trusted_contigs", limit=100)
    if len(first) != len(second):
        raise ConfigurationError("reads_1 and reads_2 must be given together.")
    if not (first or inter or single or nano or pacbio):
        raise ConfigurationError("Give at least one read file.")
    if careful and selected in ("metagenome", "rna", "rnaviral", "biosynthetic"):
        raise ConfigurationError(f"careful cannot be combined with mode {selected!r}.")
    if selected in ("rna", "rnaviral") and (nano or pacbio):
        raise ConfigurationError("rna modes take short reads only.")
    folder = _gap_cli.new_directory(folder_arg, "output_dir")
    arguments: list[str] = ["-o", str(folder)]
    flag = _MODES[selected]
    if flag is not None:
        arguments.append(flag)
    if first:
        arguments += ["--pe1-1", first[0], "--pe1-2", second[0]]
    for item in inter:
        arguments += ["--pe1-12", item]
    for item in single:
        arguments += ["-s", item]
    for item in nano:
        arguments += ["--nanopore", item]
    for item in pacbio:
        arguments += ["--pacbio", item]
    for item in trusted:
        arguments += ["--trusted-contigs", item]
    if kmers is not None:
        values = [_gap_cli.integer(k, "kmers", minimum=1, maximum=127) for k in kmers]
        if not values or any(k % 2 == 0 for k in values):
            raise ConfigurationError("kmers must be a non-empty list of odd integers.")
        arguments += ["-k", ",".join(str(k) for k in values)]
    if careful:
        arguments.append("--careful")
    if only_assembler:
        arguments.append("--only-assembler")
    if threads is not None:
        arguments += ["-t", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if memory_gb is not None:
        arguments += ["-m", str(_gap_cli.integer(memory_gb, "memory_gb", minimum=1))]
    if phred_offset is not None:
        arguments += ["--phred-offset", _gap_cli.choice(str(phred_offset), ("33", "64"), "phred_offset")]
    result = _run("spades.py", arguments, cwd=folder, expected=[], spades_path=spades_path,
                  timeout_seconds=timeout_seconds)
    found = [str(folder / name) for name in _PRIMARY if (folder / name).is_file()]
    result["mode"] = selected
    result["outputs"] = found
    if result["status"] == "completed" and not any(Path(p).name != "spades.log" for p in found):
        result["status"] = "failed"
    return result


__all__ = ["spades_assemble"]

_register_public_verbose(globals(), tuple(__all__))
