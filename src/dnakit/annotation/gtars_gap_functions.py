"""gtars VRS operation on local stores: HGVS to a GA4GH VRS allele identifier.

``hgvs_to_vrs_allele_id`` maps an HGVS expression on a transcript (c., n.) or genome (g.) to a VRS allele
identifier with gtars. It needs three local inputs that the user builds once: a transcript store file
(``gtars.reftx.TxStoreBuilder`` output), a refget sequence store directory (``gtars.refget.RefgetStore``), and
the collection digest of the reference sequences that the store holds. Nothing is downloaded.
"""

from __future__ import annotations

import importlib.metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-gtars"


def _version() -> str:
    try:
        return importlib.metadata.version("gtars")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def hgvs_to_vrs_allele_id(
    hgvs: str,
    transcript_store: str,
    refget_store: str,
    collection_digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map an HGVS variant to its GA4GH VRS allele identifier with gtars (gtars.vrs.hgvs.hgvs_to_vrs_id).

    ``transcript_store`` is a built transcript store file, ``refget_store`` a local refget store directory,
    and ``collection_digest`` the digest of the reference collection in that store. A reference mismatch (the
    HGVS reference base does not match the sequence) is reported as a ConfigurationError.
    """
    if not isinstance(hgvs, str) or not hgvs.strip():
        raise ConfigurationError("hgvs must be a non-empty HGVS string.")
    if not isinstance(collection_digest, str) or not collection_digest.strip():
        raise ConfigurationError("collection_digest must be a non-empty digest string.")
    store_file = Path(transcript_store).expanduser().resolve()
    if not store_file.is_file():
        raise ConfigurationError(f"transcript_store must be an existing file: {store_file}")
    store_dir = Path(refget_store).expanduser().resolve()
    if not store_dir.is_dir():
        raise ConfigurationError(f"refget_store must be an existing directory: {store_dir}")
    hgvs_module = require("gtars.vrs.hgvs", _EXTRA)
    reftx = require("gtars.reftx", _EXTRA)
    refget = require("gtars.refget", _EXTRA)
    try:
        provider = reftx.ReftxProvider(reftx.ReadonlyTxStore.open(str(store_file)))
    except reftx.TxStoreError as exc:
        raise ConfigurationError(f"transcript_store could not be opened: {exc}") from exc
    sequences = refget.RefgetStore.open_local(str(store_dir))
    try:
        vrs_id = hgvs_module.hgvs_to_vrs_id(hgvs.strip(), provider, sequences, collection_digest.strip())
    except hgvs_module.HgvsError as exc:
        raise ConfigurationError(f"{hgvs!r} could not be mapped: {exc}") from exc
    return {
        "status": "completed",
        "backend": "gtars",
        "version": _version(),
        "hgvs": hgvs.strip(),
        "vrs_allele_id": str(vrs_id),
        "transcript_store": str(store_file),
        "refget_store": str(store_dir),
        "collection_digest": collection_digest.strip(),
    }


__all__ = ["hgvs_to_vrs_allele_id"]
_register_public_verbose(globals())
