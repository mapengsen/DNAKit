"""varcode command-line operations that DNAKit did not wrap (Fig. 2 coverage gap, batch 0 review).

Two public functions, both built on the installed varcode and pyensembl packages:

- ``varcode_install_reference_data`` downloads and indexes the Ensembl reference data (GTF and cDNA
  FASTA) for reference names such as ``GRCh38``. This is the network step of the varcode command line
  (``download_and_install_reference_data``); files that are already present are not downloaded again.
- ``varcode_variant_gene_overlaps`` lists the genes that overlap explicit variants in a local GTF and
  genome FASTA (the varcode ``genes`` command), using the same local reference builder as the effect
  predictions in ``varcode_gap_functions``.

Missing packages raise ``BackendUnavailableError`` with the ``dnakit[bio-varcode]`` or
``dnakit[bio-pyensembl]`` install hint. The download function needs network access; its tests use a
fake pyensembl module and never reach the network.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError, DownloadError
from dnakit.molecular_lab.varcode_gap_functions import (
    _PYENSEMBL_EXTRA,
    _VARCODE_EXTRA,
    _document,
    _genome,
)

_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.\-]{1,64}$")
_BASES_PATTERN = re.compile(r"^[ACGTN]{1,256}$")
_CONTIG_PATTERN = re.compile(r"^[A-Za-z0-9_.|:\-]{1,128}$")
_MAX_NAMES = 16
_MAX_VARIANTS = 10000


def _names(value: object) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(
            "reference_names must be a list of reference names such as GRCh38."
        )
    if not 1 <= len(value) <= _MAX_NAMES:
        raise ConfigurationError(f"reference_names must hold 1 to {_MAX_NAMES} names.")
    names: list[str] = []
    for item in value:
        if not isinstance(item, str) or not _NAME_PATTERN.match(item):
            raise ConfigurationError(f"Invalid reference name {item!r}.")
        if item not in names:
            names.append(item)
    return names


def _integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer (1-based).")
    return value


def _variant(item: object, index: int) -> dict[str, Any]:
    if not isinstance(item, Mapping):
        raise ConfigurationError(
            f"variants[{index}] must be a mapping with contig, position, ref and alt."
        )
    missing = [key for key in ("contig", "position", "ref", "alt") if key not in item]
    if missing:
        raise ConfigurationError(f"variants[{index}] is missing {', '.join(missing)}.")
    contig = item["contig"]
    if not isinstance(contig, str) or not _CONTIG_PATTERN.match(contig):
        raise ConfigurationError(
            f"variants[{index}].contig must be a chromosome name such as chr1."
        )
    ref = item["ref"]
    alt = item["alt"]
    for label, allele in (("ref", ref), ("alt", alt)):
        if not isinstance(allele, str) or not _BASES_PATTERN.match(allele.upper()):
            raise ConfigurationError(
                f"variants[{index}].{label} must be a DNA allele (A, C, G, T or N)."
            )
    return {
        "contig": contig,
        "position": _integer(item["position"], f"variants[{index}].position"),
        "ref": ref.upper(),
        "alt": alt.upper(),
    }


def varcode_install_reference_data(
    reference_names: Sequence[str],
    *,
    allow_older_downloaded_release: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and index the Ensembl reference data that varcode needs, for each reference name.

    For every distinct name (for example ``GRCh38``) pyensembl picks the release: with
    ``allow_older_downloaded_release`` true the newest release that is already downloaded is used,
    otherwise the newest release. When that release's GTF and FASTA files are missing, they are
    downloaded from Ensembl and the pyensembl index is built, as the varcode command line does for the
    variant collections it parsed. Present files are never downloaded again. The result lists, per
    reference, the release, species, cache directory and whether files were downloaded and indexed.
    Network access happens only for missing files; a download failure raises DownloadError.
    """
    names = _names(reference_names)
    allow_older = bool(allow_older_downloaded_release)
    pyensembl = require("pyensembl", _PYENSEMBL_EXTRA)
    rows: list[dict[str, Any]] = []
    for name in names:
        try:
            genome = pyensembl.genome_for_reference_name(
                name, allow_older_downloaded_release=allow_older
            )
        except ValueError as exc:
            raise ConfigurationError(f"Unknown reference name {name!r}: {exc}") from exc
        present = bool(genome.required_local_files_exist())
        downloaded = False
        indexed = False
        if not present:
            try:
                genome.download()
                downloaded = True
            except (OSError, ValueError, KeyError) as exc:
                raise DownloadError(
                    f"Could not download the Ensembl files for {name}: {exc}",
                    hint="Check the network connection to the Ensembl FTP server and retry.",
                ) from exc
            try:
                genome.index()
                indexed = True
            except (OSError, ValueError, KeyError) as exc:
                raise BackendExecutionError(
                    f"pyensembl could not index the downloaded files for {name}: {exc}"
                ) from exc
        rows.append(
            {
                "reference_name": name,
                "release": int(genome.release),
                "species": str(genome.species.latin_name),
                "cache_directory": str(genome.download_cache.cache_directory_path),
                "files_present_before": present,
                "downloaded": downloaded,
                "indexed": indexed,
            }
        )
    return _document(
        "install_reference_data",
        references=rows,
        allow_older_downloaded_release=allow_older,
    )


def varcode_variant_gene_overlaps(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    cache_dir: str | None = None,
    raise_on_error: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the genes that overlap each variant, using a local GTF and genome FASTA (varcode ``genes``).

    Each variant is a mapping with ``contig`` (a chromosome name as in the GTF), ``position`` (1-based),
    ``ref`` and ``alt`` (DNA alleles). The reference index is built once per content digest of the GTF
    and FASTA, as in the varcode effect predictions. For each variant the result holds the overlapping
    genes (gene_id, gene_name, biotype, contig, start, end, strand) and, when a variant cannot be
    annotated (for example its contig is absent from the GTF), an error message. With
    ``raise_on_error`` true the first such failure raises BackendExecutionError instead.
    """
    if not isinstance(variants, Sequence) or isinstance(variants, (str, bytes)):
        raise ConfigurationError("variants must be a list of variant mappings.")
    if not 1 <= len(variants) <= _MAX_VARIANTS:
        raise ConfigurationError(f"variants must hold 1 to {_MAX_VARIANTS} entries.")
    parsed = [_variant(item, index) for index, item in enumerate(variants)]
    genome, reference = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    results: list[dict[str, Any]] = []
    failures = 0
    for index, item in enumerate(parsed):
        label = f"{item['contig']}:{item['position']}:{item['ref']}>{item['alt']}"
        entry: dict[str, Any] = {
            "variant_index": index,
            "variant": label,
            **item,
            "genes": [],
            "error": None,
        }
        try:
            variant = varcode.Variant(
                item["contig"],
                start=item["position"],
                ref=item["ref"],
                alt=item["alt"],
                genome=genome,
            )
            genes = variant.genes  # a property in varcode: the overlapping pyensembl Gene objects
        except (ValueError, KeyError) as exc:
            if raise_on_error:
                raise BackendExecutionError(f"varcode could not annotate {label}: {exc}") from exc
            failures += 1
            entry["error"] = str(exc)[:300]
            results.append(entry)
            continue
        entry["genes"] = [
            {
                "gene_id": gene.gene_id,
                "gene_name": gene.gene_name,
                "biotype": gene.biotype,
                "contig": gene.contig,
                "start": int(gene.start),
                "end": int(gene.end),
                "strand": gene.strand,
            }
            for gene in genes
        ]
        results.append(entry)
    return _document(
        "variant_gene_overlaps",
        variants=results,
        variant_count=len(results),
        annotation_failures=failures,
        reference=reference,
    )


__all__ = [
    "varcode_install_reference_data",
    "varcode_variant_gene_overlaps",
]

_register_public_verbose(globals())
