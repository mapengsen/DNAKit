"""Catalogue transforms from genomepy's UCSC, Ensembl and GENCODE genome listings.

genomepy builds a genome catalogue (a dictionary keyed by assembly name) from the provider web services
and from the UCSC MySQL server. The functions below perform the record transforms of that build on
explicit inputs that the caller supplies: assembly accessions from UCSC equivalence rows or from
sourceName fields, annotation table names from UCSC information_schema rows, and the GRCh37 records that
Ensembl and GENCODE add to their listings. They open no network or database connection. The downloads
themselves stay in the existing genomepy_* listing wrappers.

Every function returns a new catalogue and leaves its input unchanged. Where the upstream step mutates
its input in place, this is the only difference in calling convention; the differences in results are
listed in each docstring.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError

_ACCESSION = re.compile(r"(GC[AF]_\d{9}\.\d+)")
_ENSEMBL_GRCH37 = {
    "genome": "http://ftp.ensembl.org/pub/grch37/release-{}/fasta/homo_sapiens/dna",
    "annotation": (
        "http://ftp.ensembl.org/pub/grch37/release-{}/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
    ),
    "assembly_accession": "GCA_000001405.14",
    "taxonomy_id": 9606,
    "name": "Homo_sapiens",
    "scientific_name": "Homo sapiens",
    "url_name": "human",
    "assembly_name": "GRCh37",
    "division": "vertebrates",
    "base_count": "3137144693",
    "display_name": "",
    "genebuild": "",
}


def _catalogue(genomes: object, name: str) -> dict[str, dict[str, Any]]:
    if not isinstance(genomes, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of assembly name to record.")
    out: dict[str, dict[str, Any]] = {}
    for key, record in genomes.items():
        if not isinstance(key, str):
            raise ConfigurationError(f"{name} keys must be assembly names (str), not {key!r}.")
        if not isinstance(record, Mapping):
            raise ConfigurationError(f"{name}[{key!r}] must be a mapping.")
        copied = dict(record)
        for field, value in record.items():
            if isinstance(value, list):
                copied[field] = list(value)
        out[key] = copied
    return out


def _rows(rows: object, width: int, name: str) -> list[tuple[Any, ...]]:
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise ConfigurationError(f"{name} must be a sequence of {width}-field rows.")
    out: list[tuple[Any, ...]] = []
    for row in rows:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != width:
            raise ConfigurationError(f"each entry of {name} must have {width} fields, got {row!r}.")
        out.append(tuple(row))
    return out


def genomepy_ucsc_assembly_accessions(
    genomes: Mapping[str, Mapping[str, Any]],
    equivalences: Sequence[Sequence[Any]],
    *,
    verbose: bool = True,
) -> dict[str, dict[str, Any]]:
    """Set ``assembly_accession`` from UCSC assembly-equivalence rows.

    Upstream step: ``genomepy.providers.ucsc.add_accessions1`` (ucsc.py:439-485), which reads the
    ``asmEquivalent`` table through ``query_ucsc``. Here the rows are passed in.

    Each row is ``(source, destination, matchCount)``. For every source assembly the rows with the
    highest match count are kept. The GCA or GCF accession is taken from the destination name. When an
    assembly keeps more than one row, only the GCA accessions are used, as upstream does. Assemblies
    that are not in ``genomes`` are ignored.

    Difference from upstream: a kept row whose destination has no accession is not assigned. Upstream
    stores NaN in that case.

    Parameters
    ----------
    genomes : Mapping[str, Mapping[str, Any]]
        Catalogue keyed by UCSC assembly name.
    equivalences : Sequence[Sequence[Any]]
        Rows ``(source, destination, matchCount)``.
    verbose : bool
        Accepted for the verbosity contract; the transform prints nothing.

    Returns
    -------
    dict[str, dict[str, Any]]
        A new catalogue with ``assembly_accession`` set where an accession was found.
    """
    out = _catalogue(genomes, "genomes")
    rows = _rows(equivalences, 3, "equivalences")
    if not rows:
        return out
    best: dict[str, Any] = {}
    for source, _destination, match in rows:
        best[source] = match if source not in best else max(best[source], match)
    kept = [(source, destination) for source, destination, match in rows if match == best[source]]
    counts = Counter(source for source, _destination in kept)
    for source, destination in kept:
        found = _ACCESSION.search(str(destination))
        if found is None:
            continue
        accession = found.group(1)
        if counts[source] > 1 and accession[2] != "A":
            continue
        if source in out:
            out[source]["assembly_accession"] = accession
    return out


def genomepy_ucsc_sourcename_accessions(
    genomes: Mapping[str, Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, dict[str, Any]]:
    """Set ``assembly_accession`` from the ``sourceName`` field of each UCSC record.

    Upstream step: ``genomepy.providers.ucsc.add_accessions2`` (ucsc.py:488-499). The first GCA or GCF
    accession in ``sourceName`` is stored.

    Differences from upstream: records without ``sourceName`` raise ``ConfigurationError`` (upstream
    raises KeyError); a ``None`` sourceName is skipped (upstream raises TypeError); the input dict is
    not modified.

    Parameters
    ----------
    genomes : Mapping[str, Mapping[str, Any]]
        Catalogue whose records carry a ``sourceName`` string.
    verbose : bool
        Accepted for the verbosity contract; the transform prints nothing.

    Returns
    -------
    dict[str, dict[str, Any]]
        A new catalogue with ``assembly_accession`` set where ``sourceName`` contains an accession.
    """
    out = _catalogue(genomes, "genomes")
    for name, record in out.items():
        if "sourceName" not in record:
            raise ConfigurationError(f"genome record {name!r} has no sourceName field.")
        source_name = record["sourceName"]
        if source_name is None:
            continue
        if not isinstance(source_name, str):
            raise ConfigurationError(f"sourceName of {name!r} must be a string.")
        hit = _ACCESSION.search(source_name)
        if hit:
            record["assembly_accession"] = hit.group(1)
    return out


def genomepy_ucsc_annotation_links(
    genomes: Mapping[str, Mapping[str, Any]],
    table_rows: Sequence[Sequence[Any]],
    *,
    verbose: bool = True,
) -> dict[str, dict[str, Any]]:
    """Append the annotation tables that UCSC lists for each genome.

    Upstream step: ``genomepy.providers.ucsc.add_annotation_links`` (ucsc.py:502-520). The upstream
    query selects the tables ``ensGene``, ``ncbiRefSeq``, ``knownGene`` and ``refGene`` from
    ``information_schema``. Here the ``(database, table)`` rows of that result are passed in, and each
    table name is appended to the ``annotations`` list of the matching genome.

    Difference from upstream: a genome record without an ``annotations`` list raises
    ``ConfigurationError`` (upstream raises KeyError). Rows for genomes that are not in ``genomes`` are
    ignored, as upstream does.

    Parameters
    ----------
    genomes : Mapping[str, Mapping[str, Any]]
        Catalogue keyed by UCSC database name, whose records have an ``annotations`` list.
    table_rows : Sequence[Sequence[Any]]
        Rows ``(database, table)`` in the order the query returns them.
    verbose : bool
        Accepted for the verbosity contract; the transform prints nothing.

    Returns
    -------
    dict[str, dict[str, Any]]
        A new catalogue whose ``annotations`` lists include the table names.
    """
    out = _catalogue(genomes, "genomes")
    rows = _rows(table_rows, 2, "table_rows")
    for database, table in rows:
        if database not in out:
            continue
        if "annotations" not in out[database]:
            raise ConfigurationError(f"genome record {database!r} has no annotations list.")
        out[database]["annotations"].append(table)
    return out


def genomepy_ensembl_add_grch37(
    genomes: Mapping[str, Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, dict[str, Any]]:
    """Add the GRCh37 record that genomepy lists separately from the Ensembl release catalogue.

    Upstream step: ``genomepy.providers.ensembl.add_grch37`` (ensembl.py:278-298). The record is a fixed
    literal (Homo sapiens, GCA_000001405.14, Ensembl GRCh37 release 87 GTF); it is inserted under the
    key ``GRCh37``. Any existing ``GRCh37`` entry is replaced, as upstream does.

    Parameters
    ----------
    genomes : Mapping[str, Mapping[str, Any]]
        Ensembl catalogue keyed by assembly name.
    verbose : bool
        Accepted for the verbosity contract; the transform prints nothing.

    Returns
    -------
    dict[str, dict[str, Any]]
        A new catalogue with the GRCh37 record added.
    """
    out = _catalogue(genomes, "genomes")
    out["GRCh37"] = {key: value for key, value in _ENSEMBL_GRCH37.items()}
    return out


def genomepy_gencode_add_grch37(
    genomes: Mapping[str, Mapping[str, Any]],
    ftp_link: str,
    *,
    verbose: bool = True,
) -> dict[str, dict[str, Any]]:
    """Add the GRCh37 annotation record that GENCODE files under its own mapping folder.

    Upstream step: ``genomepy.providers.gencode.add_grch37`` (gencode.py:228-243). The latest human
    GRCh assembly in the catalogue gives the GENCODE release number, taken from the ``release_`` path
    segment of its first annotation link. The GRCh37 lift-over annotation URL is built from that
    release under ``ftp_link``.

    Difference from upstream: only keys that start with ``GRCh`` are considered when finding the latest
    human assembly (upstream takes the maximum over every key), and a catalogue without such a key or
    without a ``release_`` segment raises ``ConfigurationError``. The input dict is not modified.

    Parameters
    ----------
    genomes : Mapping[str, Mapping[str, Any]]
        GENCODE catalogue keyed by assembly name, with ``annotations`` lists.
    ftp_link : str
        Base FTP link of the GENCODE listing, for example ``ftp://ftp.ebi.ac.uk/pub/databases/gencode``.
    verbose : bool
        Accepted for the verbosity contract; the transform prints nothing.

    Returns
    -------
    dict[str, dict[str, Any]]
        A new catalogue with the GRCh37 record added.
    """
    out = _catalogue(genomes, "genomes")
    if not isinstance(ftp_link, str) or not ftp_link.strip():
        raise ConfigurationError("ftp_link must be a non-empty string.")
    grch = [key for key in out if key.startswith("GRCh") and key[4:].isdigit()]
    if not grch:
        raise ConfigurationError(
            "the GENCODE catalogue has no GRCh assembly to take the release from."
        )
    latest = f"GRCh{max(int(key[4:]) for key in grch)}"
    annotations = out[latest].get("annotations") or []
    if not annotations:
        raise ConfigurationError(f"{latest} has no annotation link to read the release from.")
    segments = [part for part in str(annotations[0]).split("/") if "release_" in part]
    if not segments:
        raise ConfigurationError(f"annotation link of {latest} has no release_ segment.")
    release = segments[0][-2:]
    out["GRCh37"] = {
        "annotations": [
            f"{ftp_link}/Gencode_human/release_{release}/GRCh37_mapping/"
            f"gencode.v{release}lift37.annotation.gtf.gz"
        ],
        "taxonomy_id": 9606,
        "species": "Homo sapiens",
        "text_search": "human",
    }
    return out


_register_public_verbose(globals())
