"""Explicit gget Ensembl downloads and gene-isoform maps (P41; network)."""

from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def gget_download_reference(
    species: str,
    output_dir: str,
    *,
    which: str = "dna",
    release: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resolve Ensembl files with gget.ref then download them through DNAKit's bounded downloader.

    Network operation. which=dna/gtf/cdna/cds/cdrna/all; files are never overwritten.
    The gget CLI uses the same ref(ftp=True) URLs and curl; this interface uses
    the existing atomic, checksummed DNAKit HTTPS downloader instead.
    """
    kind = h.choice(which, ("dna", "gtf", "cdna", "cds", "cdrna", "all"), "which")
    if release is not None:
        h.integer(release, "release", minimum=1)
    folder = Path(output_dir).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError("output_dir must be an existing directory.")
    module = require("gget", "bio-gget")
    urls = module.ref(
        h.text(species, "species"),
        which=kind,
        release=release,
        ftp=True,
        save=False,
        verbose=verbose,
    )
    if not isinstance(urls, list) or not urls or len(urls) > 100:
        raise ConfigurationError("gget did not return a bounded nonempty list of URLs.")
    targets = []
    for url in urls:
        parsed = urlparse(url)
        if parsed.hostname is None or not (
            parsed.hostname == "ensembl.org" or parsed.hostname.endswith(".ensembl.org")
        ):
            raise ConfigurationError("Reference URLs must be hosted by Ensembl.")
        name = Path(unquote(parsed.path)).name
        if not name or name in (".", ".."):
            raise ConfigurationError("Invalid reference filename.")
        targets.append(h.new_file(folder / name, "reference output"))
    if len(set(targets)) != len(targets):
        raise ConfigurationError("Reference URLs have conflicting filenames.")
    from dnakit.download.files import download_url

    for url, target in zip(urls, targets, strict=False):
        download_url(
            url.replace("ftp://", "https://", 1),
            target.parent,
            filename=target.name,
            verbose=verbose,
        )
    return {
        "status": "completed",
        "backend": "gget",
        "release": release,
        "retrieved_on": date.today().isoformat(),
        "urls": urls,
        "outputs": [str(x) for x in targets],
    }


def gget_isoform_structure_map(
    gene: str, *, isoform: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Query gget.g2p(resource='map') for an isoform's genomic structure (network)."""
    module = require("gget", "bio-gget")
    table = module.g2p(
        gene=h.text(gene, "gene"), resource="map", isoform=isoform, save=False, verbose=verbose
    )
    import json

    rows = None if table is None else json.loads(table.to_json(orient="records"))
    return {
        "status": "completed",
        "backend": "gget",
        "retrieved_on": date.today().isoformat(),
        "records": rows,
    }


__all__ = ["gget_download_reference", "gget_isoform_structure_map"]
_register_public_verbose(globals(), tuple(__all__))
