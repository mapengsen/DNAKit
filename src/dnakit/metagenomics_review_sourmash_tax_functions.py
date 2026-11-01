"""Taxonomic profile files from sourmash gather results (sourmash tax metagenome), one function per output format.

The inputs are gather CSV files written by ``sourmash gather`` and taxonomy lineage CSV files (an ``ident`` column plus
one column per rank). Each function writes one new file next to ``output_base``: ``<base>.krona.tsv``,
``<base>.summarized.csv``, ``<base>.human.txt``, ``<base>.lineage_summary.tsv`` or ``<base>.bioboxes.profile``.
Nothing is downloaded; the user supplies every file. Results are JSON-serialisable dictionaries with status, backend,
command, outputs and a bounded log tail.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.metagenomics_review_sourmash_functions import _existing, _nonempty, _run

_RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain")


def _tax_metagenome(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    output_format: str,
    suffix: str,
    *,
    rank: str | None,
    sourmash_path: str | None,
    verbose: bool,
    extra: Sequence[str] = (),
) -> dict[str, Any]:
    gathers = [str(_existing(p, "gather_csv")) for p in _nonempty(gather_csvs, "gather_csvs")]
    taxonomies = [str(_existing(p, "taxonomy_csv")) for p in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    base = Path(output_base).expanduser().resolve()
    if not base.parent.is_dir():
        raise ConfigurationError(f"The folder for output_base does not exist: {base.parent}")
    outfile = base.parent / (base.name + suffix)
    if outfile.exists():
        raise ConfigurationError(f"Output file must not exist yet: {outfile}")
    args = ["tax", "metagenome", "-g", *gathers, "-t", *taxonomies, "-F", output_format, "-o", str(base), *extra]
    if rank is not None:
        if rank not in _RANKS:
            raise ConfigurationError(f"rank must be one of {list(_RANKS)}.")
        args += ["--rank", rank]
    return _run(args, [outfile], sourmash_path, verbose)


def krona_taxonomy_profile(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    rank: str,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a Krona table of the fraction of each query assigned to each lineage at one rank (sourmash tax metagenome
    -F krona). The Krona formatting (aggregation by lineage at the rank) is done by sourmash; ``rank`` is required."""
    if rank not in _RANKS:
        raise ConfigurationError(f"rank must be one of {list(_RANKS)} for Krona output.")
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "krona", ".krona.tsv",
        rank=rank, sourmash_path=sourmash_path, verbose=verbose,
    )


def csv_taxonomy_summary(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the per-rank summary of gather results as CSV with columns query_name, rank, fraction, lineage and the
    matching statistics (sourmash tax metagenome -F csv_summary)."""
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "csv_summary", ".summarized.csv",
        rank=None, sourmash_path=sourmash_path, verbose=verbose,
    )


def human_taxonomy_summary(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    rank: str | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a human-readable table of the taxonomic lineages found in each query (sourmash tax metagenome -F human).
    ``rank`` limits the table to one rank; without it sourmash reports the lineages it finds."""
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "human", ".human.txt",
        rank=rank, sourmash_path=sourmash_path, verbose=verbose,
    )


def lineage_taxonomy_summary(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    rank: str,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the fraction of the query assigned to each lineage at one rank as a two-column table (sourmash tax
    metagenome -F lineage_summary). ``rank`` is required by sourmash for this format."""
    if rank not in _RANKS:
        raise ConfigurationError(f"rank must be one of {list(_RANKS)} for lineage_summary output.")
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "lineage_summary", ".lineage_summary.tsv",
        rank=rank, sourmash_path=sourmash_path, verbose=verbose,
    )


def bioboxes_taxonomy_profile(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the CAMI bioboxes profile of gather results (sourmash tax metagenome -F bioboxes). The lineage CSV must
    carry NCBI taxids in a ``taxpath`` column ('|'-separated, one per rank). Upstream sourmash 4.9.4 raises a
    TypeError when a reported lineage has no taxid, so such queries end with status ``failed`` and the log tail."""
    taxonomies = [_existing(p, "taxonomy_csv") for p in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    for path in taxonomies:
        with path.open(encoding="utf-8", newline="") as handle:
            header = handle.readline().strip().split(",")
        if "taxpath" not in header:
            raise ConfigurationError(
                f"bioboxes output needs a 'taxpath' column with NCBI taxids in the taxonomy CSV: {path}"
            )
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "bioboxes", ".bioboxes.profile",
        rank=None, sourmash_path=sourmash_path, verbose=verbose,
    )


def kreport_taxonomy_report(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a Kraken-style report of the gather results: percentage, reads, direct reads, rank code, taxid and name
    for each lineage (sourmash tax metagenome -F kreport). Unclassified query fraction is reported as ``U``."""
    return _tax_metagenome(
        gather_csvs, taxonomy_csvs, output_base, "kreport", ".kreport.txt",
        rank=None, sourmash_path=sourmash_path, verbose=verbose,
    )


def lingroup_taxonomy_report(
    gather_csvs: Sequence[str],
    lin_taxonomy_csvs: Sequence[str],
    lingroup_csv: str,
    output_base: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the gather results by LIN lingroup (sourmash tax metagenome -F lingroup with --lins and --lingroup).
    ``lin_taxonomy_csvs`` must carry a ``lin`` column; ``lingroup_csv`` has columns ``name`` and ``lin`` giving the
    prefix of each group."""
    lineages = [_existing(p, "lin_taxonomy_csv") for p in _nonempty(lin_taxonomy_csvs, "lin_taxonomy_csvs")]
    for path in lineages:
        with path.open(encoding="utf-8", newline="") as handle:
            header = handle.readline().strip().split(",")
        if "lin" not in header:
            raise ConfigurationError(f"lingroup output needs a 'lin' column in the taxonomy CSV: {path}")
    groups = _existing(lingroup_csv, "lingroup_csv")
    return _tax_metagenome(
        gather_csvs, [str(p) for p in lineages], output_base, "lingroup", ".lingroup.tsv",
        rank=None, sourmash_path=sourmash_path, verbose=verbose,
        extra=("--lins", "--lingroup", str(groups)),
    )


_register_public_verbose(globals())
