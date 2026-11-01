"""gtars refget store operations (integration pack P24, part 1).

A refget store is a local directory that holds sequences and sequence collections under their GA4GH
digests (sha512t24u), with optional aliases. These functions create, fill, query and clean such a
store through ``gtars.refget``. Every function takes the store directory as its first argument.
Reading functions need an existing store; writing functions create it when it is missing. Digests
may be given with or without the ``SQ.`` prefix. gtars is imported only when a function runs."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.gtars_functions import (
    _choice,
    _flag,
    _integer,
    _plain,
    _provenance,
    _text,
)
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-gtars"
_KINDS = ("sequence", "collection")
_ATTRIBUTES = (
    "names",
    "lengths",
    "sequences",
    "name_length_pairs",
    "sorted_name_length_pairs",
    "sorted_sequences",
)
_MAX_BASES = 100_000_000


def _refget() -> Any:
    return require("gtars.refget", _EXTRA)


def _open(store: str | Path, *, create: bool = False) -> Any:
    refget = _refget()
    path = Path(store).expanduser().resolve()
    try:
        if create:
            handle = refget.RefgetStore.on_disk(str(path))
        else:
            if not (path / "rgstore.json").is_file():
                raise ConfigurationError(f"store is not an existing refget store directory: {path}")
            handle = refget.RefgetStore.open_local(str(path))
    except ConfigurationError:
        raise
    except Exception as exc:  # gtars reports store problems with runtime errors
        raise ConfigurationError(f"gtars could not open the refget store {path}: {exc}") from exc
    handle.set_quiet(True)
    return handle


def _call(call: Any, what: str) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:  # gtars raises IOError, KeyError and RuntimeError for bad input
        raise ConfigurationError(f"gtars could not {what}: {exc}") from exc


def _digest(value: object, name: str = "digest") -> str:
    text = _text(value, name, maximum=128).strip()
    return text[3:] if text[:3].upper() == "SQ." else text


def _seq_meta(meta: Any) -> dict[str, Any]:
    description = getattr(meta, "description", None)
    return {
        "name": str(meta.name),
        "description": None if description in (None, "") else str(description),
        "length": int(meta.length),
        "alphabet": str(meta.alphabet),
        "sha512t24u": str(meta.sha512t24u),
        "md5": str(meta.md5),
    }


def _col_meta(meta: Any) -> dict[str, Any]:
    row = {
        "digest": str(meta.digest),
        "n_sequences": int(meta.n_sequences),
    }
    for field in (
        "names_digest",
        "lengths_digest",
        "sequences_digest",
        "name_length_pairs_digest",
        "sorted_name_length_pairs_digest",
        "sorted_sequences_digest",
    ):
        value = getattr(meta, field, None)
        if value is not None:
            row[field] = str(value)
    return row


def _store_text(store: str | Path) -> str:
    return str(Path(store).expanduser().resolve())


def _provenance_note(note: str) -> dict[str, str]:
    return _provenance("GA4GH refget sequences and seqcol collections; " + note)


def _resolve_collection(
    handle: Any, digest: str | None, namespace: str | None, alias: str | None
) -> str:
    if digest is not None:
        return _digest(digest)
    if namespace is None or alias is None:
        raise ConfigurationError("give either digest or both namespace and alias.")
    ns, name = _text(namespace, "namespace", maximum=128), _text(alias, "alias", maximum=512)
    meta = _call(
        lambda: handle.get_collection_metadata_by_alias(ns, name), "resolve the collection alias"
    )
    if meta is None:
        raise ConfigurationError(f"collection alias {ns}/{name} was not found.")
    return str(meta.digest)


def _resolve_sequence(
    handle: Any, digest: str | None, namespace: str | None, alias: str | None
) -> str:
    if digest is not None:
        return _digest(digest)
    if namespace is None or alias is None:
        raise ConfigurationError("give either digest or both namespace and alias.")
    ns, name = _text(namespace, "namespace", maximum=128), _text(alias, "alias", maximum=512)
    meta = _call(
        lambda: handle.get_sequence_metadata_by_alias(ns, name), "resolve the sequence alias"
    )
    if meta is None:
        raise ConfigurationError(f"sequence alias {ns}/{name} was not found.")
    return str(meta.sha512t24u)


# Filling a store.


def refget_store_add_fasta(
    store: str | Path,
    fasta_files: Sequence[str | Path],
    *,
    force: bool = False,
    use_collection_object: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Digest FASTA files and add each as a sequence collection (with its sequences) to a refget
    store.

    ``fasta_files`` may hold one file (``gtars RefgetStore.add_sequence_collection_from_fasta``)
    or many (``add_sequence_collections_from_fastas``). With ``use_collection_object=True`` each
    file is first digested into a ``SequenceCollection`` and then stored with
    ``add_sequence_collection``; the stored result is identical. Existing collections are
    skipped unless ``force`` is true. The store is created when missing."""

    if isinstance(fasta_files, (str, Path)) or not fasta_files:
        raise ConfigurationError("fasta_files must be a non-empty list of FASTA paths.")
    files = [str(input_path(item)) for item in fasta_files]
    overwrite = _flag(force, "force")
    via_object = _flag(use_collection_object, "use_collection_object")
    refget = _refget()
    handle = _open(store, create=True)
    rows = []
    if via_object:
        for path in files:
            collection = _call(lambda p=path: refget.digest_fasta(p), f"digest {path}")
            existed = _call(
                lambda c=collection: handle.get_collection_metadata(str(c.digest)),
                "check the store",
            )
            _call(
                lambda c=collection: handle.add_sequence_collection(c, force=overwrite),
                "add the collection",
            )
            rows.append(
                {
                    "file": path,
                    "digest": str(collection.digest),
                    "added": existed is None or overwrite,
                }
            )
    elif len(files) == 1:
        meta, added = _call(
            lambda: handle.add_sequence_collection_from_fasta(files[0], force=overwrite),
            "add the FASTA file",
        )
        rows.append({"file": files[0], "digest": str(meta.digest), "added": bool(added)})
    else:
        report = _call(
            lambda: handle.add_sequence_collections_from_fastas(files, jobs=1, force=overwrite),
            "add the FASTA files",
        )
        for path, (meta, added) in zip(files, report.collections, strict=True):
            rows.append({"file": path, "digest": str(meta.digest), "added": bool(added)})
    return {
        "status": "completed",
        "operation": "refget_store_add_fasta",
        "store": _store_text(store),
        "collections": rows,
        "collection_count": len(handle.list_collections()["results"]),
        "sequence_count": len(handle.list_sequences()),
        "provenance": _provenance_note("existing collections are skipped unless force=True"),
    }


def refget_store_add_sequence(
    store: str | Path,
    sequence: str,
    *,
    name: str | None = None,
    description: str | None = None,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Digest one sequence (``gtars.refget.digest_sequence``) and store it without any collection.

    The sequence is uppercased by gtars; the result holds its digests and the detected alphabet. The
    store is created when missing."""

    text = sequence_text(sequence)
    label = "" if name is None else _text(name, "name", maximum=256)
    note = None if description is None else _text(description, "description", maximum=1024)
    overwrite = _flag(force, "force")
    refget = _refget()
    record = _call(
        lambda: refget.digest_sequence(text.encode("ascii"), name=label, description=note), "digest"
    )
    handle = _open(store, create=True)
    _call(lambda: handle.add_sequence(record, force=overwrite), "add the sequence")
    _call(handle.write, "write the store index")
    return {
        "status": "completed",
        "operation": "refget_store_add_sequence",
        "store": _store_text(store),
        **_seq_meta(record.metadata),
        "provenance": _provenance_note("sequence is uppercased before digesting"),
    }


# Aliases.


def refget_alias_register(
    store: str | Path,
    kind: str,
    namespace: str,
    *,
    alias: str | None = None,
    digest: str | None = None,
    alias_file: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Register aliases for sequences or collections in a refget store.

    ``kind`` is ``sequence`` or ``collection``. Give one alias with ``alias`` and ``digest``
    (``add_sequence_alias`` / ``add_collection_alias``), or a two-column TSV file (alias, digest)
    with ``alias_file`` (``load_sequence_aliases`` / ``load_collection_aliases``). Digests in an
    alias file must be written without the ``SQ.`` prefix."""

    what = _choice(kind, "kind", _KINDS)
    ns = _text(namespace, "namespace", maximum=128)
    single = alias is not None or digest is not None
    if single == (alias_file is not None):
        raise ConfigurationError(
            "give either alias and digest, or alias_file, not both or neither."
        )
    handle = _open(store)
    if single:
        name = _text(alias, "alias", maximum=512)
        target = _digest(digest)
        if what == "sequence":
            _call(lambda: handle.add_sequence_alias(ns, name, target), "add the sequence alias")
        else:
            _call(lambda: handle.add_collection_alias(ns, name, target), "add the collection alias")
        count = 1
    else:
        path = str(input_path(alias_file))  # type: ignore[arg-type]
        loader = (
            handle.load_sequence_aliases if what == "sequence" else handle.load_collection_aliases
        )
        count = int(_call(lambda: loader(ns, path), "load the alias file"))
    return {
        "status": "completed",
        "operation": "refget_alias_register",
        "store": _store_text(store),
        "kind": what,
        "namespace": ns,
        "aliases_registered": count,
        "provenance": _provenance_note("aliases are stored in the store directory"),
    }


def refget_alias_remove(
    store: str | Path,
    kind: str,
    namespace: str,
    alias: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove one sequence or collection alias (``kind`` = ``sequence`` or ``collection``) from a
    store."""

    what = _choice(kind, "kind", _KINDS)
    ns, name = _text(namespace, "namespace", maximum=128), _text(alias, "alias", maximum=512)
    handle = _open(store)
    remover = handle.remove_sequence_alias if what == "sequence" else handle.remove_collection_alias
    removed = bool(_call(lambda: remover(ns, name), "remove the alias"))
    return {
        "status": "completed",
        "operation": "refget_alias_remove",
        "store": _store_text(store),
        "kind": what,
        "namespace": ns,
        "alias": name,
        "removed": removed,
        "provenance": _provenance_note("removed is false when the alias did not exist"),
    }


def refget_alias_list(
    store: str | Path,
    kind: str,
    *,
    namespace: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List alias namespaces of a store, or the aliases of one namespace (``kind`` = sequence or
    collection)."""

    what = _choice(kind, "kind", _KINDS)
    handle = _open(store)
    if namespace is None:
        names = (
            handle.list_sequence_alias_namespaces()
            if what == "sequence"
            else handle.list_collection_alias_namespaces()
        )
        return {
            "status": "completed",
            "operation": "refget_alias_list",
            "store": _store_text(store),
            "kind": what,
            "namespaces": sorted(str(item) for item in names),
            "provenance": _provenance_note("namespaces sorted by name"),
        }
    ns = _text(namespace, "namespace", maximum=128)
    lister = handle.list_sequence_aliases if what == "sequence" else handle.list_collection_aliases
    aliases = _call(lambda: lister(ns), "list the aliases")
    return {
        "status": "completed",
        "operation": "refget_alias_list",
        "store": _store_text(store),
        "kind": what,
        "namespace": ns,
        "aliases": None if aliases is None else sorted(str(item) for item in aliases),
        "provenance": _provenance_note("aliases is null when the namespace does not exist"),
    }


def refget_alias_reverse_lookup(
    store: str | Path,
    kind: str,
    digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return every (namespace, alias) pair that points at a sequence or collection digest."""

    what = _choice(kind, "kind", _KINDS)
    target = _digest(digest)
    handle = _open(store)
    getter = (
        handle.get_aliases_for_sequence if what == "sequence" else handle.get_aliases_for_collection
    )
    pairs = _call(lambda: getter(target), "look up the aliases")
    rows = sorted(
        ({"namespace": str(ns), "alias": str(al)} for ns, al in pairs),
        key=lambda row: (row["namespace"], row["alias"]),
    )
    return {
        "status": "completed",
        "operation": "refget_alias_reverse_lookup",
        "store": _store_text(store),
        "kind": what,
        "digest": target,
        "aliases": rows,
        "provenance": _provenance_note("aliases sorted by namespace then alias"),
    }


def refget_alias_pull(
    remote_url: str,
    cache_dir: str | Path,
    *,
    namespace: str | None = None,
    strategy: str = "keep-ours",
    verbose: bool = True,
) -> dict[str, Any]:
    """Pull alias tables from a remote refget store into a local cache store (network).

    ``strategy`` is ``keep-ours`` (default), ``keep-theirs`` or ``notify`` for conflicting aliases.
    With ``namespace=None`` every namespace the remote store advertises is pulled."""

    url = _text(remote_url, "remote_url", maximum=2048)
    if not url.startswith(("http://", "https://", "s3://")):
        raise ConfigurationError("remote_url must start with http://, https:// or s3://.")
    mode = _choice(strategy, "strategy", ("keep-ours", "keep-theirs", "notify"))
    ns = None if namespace is None else _text(namespace, "namespace", maximum=128)
    refget = _refget()
    cache = Path(cache_dir).expanduser().resolve()
    handle = _call(lambda: refget.RefgetStore.open_remote(str(cache), url), "open the remote store")
    handle.set_quiet(True)
    report = _call(lambda: handle.pull_aliases(namespace=ns, strategy=mode), "pull the aliases")
    return {
        "status": "completed",
        "operation": "refget_alias_pull",
        "remote_url": url,
        "cache_dir": str(cache),
        "strategy": mode,
        "report": _plain(report),
        "provenance": _provenance_note(
            "network operation; result lists pulled, skipped, not_found, conflicts"
        ),
    }


# Querying collections.


def refget_collection_compare(
    store: str | Path,
    digest_a: str,
    digest_b: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two collections in a store with the seqcol comparison (shared and unique attributes,
    order)."""

    first, second = _digest(digest_a, "digest_a"), _digest(digest_b, "digest_b")
    handle = _open(store)
    result = _call(lambda: handle.compare(first, second), "compare the collections")
    return {
        "status": "completed",
        "operation": "refget_collection_compare",
        "store": _store_text(store),
        "comparison": _plain(result),
        "provenance": _provenance_note("GA4GH seqcol comparison format"),
    }


def refget_collections_by_attribute(
    store: str | Path,
    attribute: str,
    attribute_digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the collections of a store that carry a given attribute digest (for example the same
    names)."""

    name = _choice(attribute, "attribute", _ATTRIBUTES)
    target = _text(attribute_digest, "attribute_digest", maximum=128)
    handle = _open(store)
    found = _call(lambda: handle.find_collections_by_attribute(name, target), "search by attribute")
    return {
        "status": "completed",
        "operation": "refget_collections_by_attribute",
        "store": _store_text(store),
        "attribute": name,
        "attribute_digest": target,
        "collections": sorted(str(item) for item in found),
        "provenance": _provenance_note("collection digests sorted"),
    }


def refget_collection_attribute_array(
    store: str | Path,
    attribute: str,
    attribute_digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the array (names, lengths or sequence digests) that stands behind an attribute
    digest."""

    name = _choice(attribute, "attribute", _ATTRIBUTES)
    target = _text(attribute_digest, "attribute_digest", maximum=128)
    handle = _open(store)
    values = _call(lambda: handle.get_attribute(name, target), "read the attribute")
    return {
        "status": "completed",
        "operation": "refget_collection_attribute_array",
        "store": _store_text(store),
        "attribute": name,
        "attribute_digest": target,
        "values": None if values is None else _plain(list(values)),
        "provenance": _provenance_note("values is null when the digest is unknown"),
    }


def refget_collection_digests(
    store: str | Path,
    digest: str,
    *,
    level: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a collection as seqcol level 1 (attribute digests) or level 2 (full arrays)
    (``level`` 1 or 2)."""

    depth = _integer(level, "level", 1, 2)
    target = _digest(digest)
    handle = _open(store)
    getter = handle.get_collection_level1 if depth == 1 else handle.get_collection_level2
    data = _call(lambda: getter(target), "read the collection")
    return {
        "status": "completed",
        "operation": "refget_collection_digests",
        "store": _store_text(store),
        "collection": target,
        "level": depth,
        "representation": _plain(data),
        "provenance": _provenance_note("level 1 holds attribute digests; level 2 holds the arrays"),
    }


def refget_collection_retrieve(
    store: str | Path,
    *,
    digest: str | None = None,
    namespace: str | None = None,
    alias: str | None = None,
    include_sequences: bool = False,
    max_bases: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve one collection by digest or by (namespace, alias): its metadata, and optionally its
    sequences.

    ``include_sequences=False`` returns the collection metadata and per-sequence metadata only
    (``get_collection_metadata``); ``True`` also loads and returns the residues
    (``get_collection``)."""

    withseq = _flag(include_sequences, "include_sequences")
    limit = _integer(max_bases, "max_bases", 1, _MAX_BASES)
    handle = _open(store)
    target = _resolve_collection(handle, digest, namespace, alias)
    meta = _call(lambda: handle.get_collection_metadata(target), "read the collection metadata")
    if meta is None:
        raise ConfigurationError(f"collection {target} is not in the store.")
    rows = []
    total = 0
    arrays = _call(lambda: handle.get_collection_level2(target), "read the collection")
    for key in arrays["sequences"]:
        info = _call(
            lambda k=key: handle.get_sequence_metadata(_digest(k)), "read the sequence metadata"
        )
        if info is None:
            raise ConfigurationError(f"sequence {key} of the collection is not in the store.")
        row = _seq_meta(info)
        total += row["length"]
        if withseq:
            if total > limit:
                raise ConfigurationError(f"collection has more than max_bases={limit} bases.")
            name = str(info.sha512t24u)
            _call(lambda n=name: handle.load_sequence(n), "load the sequence")
            row["sequence"] = _call(
                lambda n=name: handle.get_sequence(n).decode(), "read the sequence"
            )
        rows.append(row)
    return {
        "status": "completed",
        "operation": "refget_collection_retrieve",
        "store": _store_text(store),
        **_col_meta(meta),
        "total_bases": total,
        "sequences": rows,
        "provenance": _provenance_note("sequences in collection order"),
    }


def refget_collections_list(
    store: str | Path,
    *,
    attribute: str | None = None,
    attribute_digest: str | None = None,
    page: int = 0,
    page_size: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the collections of a store, optionally only those with one attribute digest (paged)."""

    index = _integer(page, "page", 0, 10**9)
    size = _integer(page_size, "page_size", 1, 10_000)
    if (attribute is None) != (attribute_digest is None):
        raise ConfigurationError("give attribute and attribute_digest together, or neither.")
    filters = None
    if attribute is not None:
        filters = {
            _choice(attribute, "attribute", _ATTRIBUTES): _text(
                attribute_digest, "attribute_digest", maximum=128
            )
        }
    handle = _open(store)
    result = _call(
        lambda: handle.list_collections(page=index, page_size=size, filters=filters),
        "list the collections",
    )
    rows = [_col_meta(item) for item in result["results"]]
    return {
        "status": "completed",
        "operation": "refget_collections_list",
        "store": _store_text(store),
        "collections": rows,
        "pagination": _plain(result["pagination"]),
        "provenance": _provenance_note("page numbers start at 0"),
    }


# Querying sequences.


def refget_sequence_retrieve(
    store: str | Path,
    *,
    digest: str | None = None,
    namespace: str | None = None,
    alias: str | None = None,
    include_sequence: bool = False,
    max_bases: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve one sequence by digest (sha512t24u or MD5) or by (namespace, alias): metadata and
    optionally residues.

    ``include_sequence=False`` returns only the sequence metadata (``get_sequence_metadata``);
    ``True`` also returns the residues (``get_sequence``)."""

    withseq = _flag(include_sequence, "include_sequence")
    limit = _integer(max_bases, "max_bases", 1, _MAX_BASES)
    handle = _open(store)
    target = _resolve_sequence(handle, digest, namespace, alias)
    meta = _call(lambda: handle.get_sequence_metadata(target), "read the sequence metadata")
    if meta is None and len(target) == 32 and all(c in "0123456789abcdefABCDEF" for c in target):
        record = _call(lambda: handle.get_sequence(target), "find the sequence by MD5")
        meta = record.metadata
    if meta is None:
        raise ConfigurationError(f"sequence {target} is not in the store.")
    result = {
        "status": "completed",
        "operation": "refget_sequence_retrieve",
        "store": _store_text(store),
        **_seq_meta(meta),
    }
    if withseq:
        if int(meta.length) > limit:
            raise ConfigurationError(
                f"sequence is {int(meta.length)} bases; max_bases={limit} is the limit."
            )
        key = str(meta.sha512t24u)
        _call(lambda: handle.load_sequence(key), "load the sequence")
        record = _call(lambda: handle.get_sequence(key), "read the sequence")
        result["sequence"] = record.decode()
    result["provenance"] = _provenance_note("residues are uppercase")
    return result


def refget_sequence_substring(
    store: str | Path,
    ranges: Sequence[Sequence[int]],
    *,
    digest: str | None = None,
    namespace: str | None = None,
    alias: str | None = None,
    collection: str | None = None,
    name: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract substrings of one stored sequence for 0-based half-open ``[start, end)`` ranges.

    The sequence is given by ``digest``, by (``namespace``, ``alias``), or by ``collection`` digest
    plus the sequence ``name`` inside it. Several ranges are served in one call (``get_substring`` /
    ``get_substrings``)."""

    if isinstance(ranges, (str, bytes)) or not ranges or len(ranges) > 100_000:
        raise ConfigurationError(
            "ranges must be a non-empty list of [start, end] pairs (at most 100000)."
        )
    spans = []
    for item in ranges:
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)) or len(item) != 2:
            raise ConfigurationError("each range must be a [start, end] pair.")
        start = _integer(item[0], "range start", 0, 10**12)
        end = _integer(item[1], "range end", 0, 10**12)
        if end < start:
            raise ConfigurationError("each range needs start <= end.")
        spans.append((start, end))
    handle = _open(store)
    if collection is not None:
        if digest is not None or alias is not None:
            raise ConfigurationError(
                "use either collection and name, or digest, or namespace and alias."
            )
        coll = _digest(collection, "collection")
        label = _text(name, "name", maximum=512)
        arrays = _call(lambda: handle.get_collection_level2(coll), "read the collection")
        if label not in arrays["names"]:
            raise ConfigurationError(f"sequence name {label!r} is not in collection {coll}.")
        target = _digest(arrays["sequences"][arrays["names"].index(label)])
    else:
        target = _resolve_sequence(handle, digest, namespace, alias)
    pieces = _call(lambda: handle.get_substrings(target, spans), "extract the substrings")
    return {
        "status": "completed",
        "operation": "refget_sequence_substring",
        "store": _store_text(store),
        "sequence": target,
        "substrings": [
            {"start": start, "end": end, "sequence": str(piece)}
            for (start, end), piece in zip(spans, pieces, strict=True)
        ],
        "provenance": _provenance_note("0-based half-open coordinates"),
    }


# Moving and removing collections.


def refget_collection_import(
    store: str | Path,
    source_store: str | Path,
    digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy one collection, with its sequences, aliases and FHR metadata, from one local store into
    another."""

    target = _digest(digest)
    source = _open(source_store)
    destination = _open(store, create=True)
    _call(lambda: destination.import_collection(source, target), "import the collection")
    meta = destination.get_collection_metadata(target)
    return {
        "status": "completed",
        "operation": "refget_collection_import",
        "store": _store_text(store),
        "source_store": _store_text(source_store),
        "imported": None if meta is None else _col_meta(meta),
        "provenance": _provenance_note("both stores must be disk-backed"),
    }


def refget_orphan_removal_plan(
    store: str | Path,
    digest: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Preview which sequences would be deleted if a collection were removed with orphan cleanup
    (changes nothing)."""

    target = _digest(digest)
    handle = _open(store)
    doomed = _call(lambda: handle.plan_orphan_removal(target), "plan the removal")
    return {
        "status": "completed",
        "operation": "refget_orphan_removal_plan",
        "store": _store_text(store),
        "collection": target,
        "orphan_sequences": sorted(str(item) for item in doomed),
        "provenance": _provenance_note("advisory scan without a lock"),
    }


def refget_collection_remove(
    store: str | Path,
    digest: str,
    *,
    remove_orphan_sequences: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove a collection (and its aliases and metadata) from a store, optionally with sequences no
    other collection uses."""

    target = _digest(digest)
    orphans = _flag(remove_orphan_sequences, "remove_orphan_sequences")
    handle = _open(store)
    removed = bool(
        _call(
            lambda: handle.remove_collection(target, remove_orphan_sequences=orphans),
            "remove the collection",
        )
    )
    return {
        "status": "completed",
        "operation": "refget_collection_remove",
        "store": _store_text(store),
        "collection": target,
        "removed": removed,
        "remove_orphan_sequences": orphans,
        "provenance": _provenance_note("removed is false when the collection was not in the store"),
    }


def refget_write_rgsi(
    fasta: str | Path,
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Digest a FASTA file and write its RGSI sequence index (collection digests plus per-sequence
    digests).

    The output file must not exist yet."""

    source = input_path(fasta)
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(
            f"Output file already exists and will not be overwritten: {target}"
        )
    refget = _refget()
    collection = _call(lambda: refget.digest_fasta(str(source)), "digest the FASTA file")
    _call(lambda: collection.write_rgsi(str(target)), "write the RGSI file")
    return {
        "status": "completed",
        "operation": "refget_write_rgsi",
        "fasta": str(source),
        "output": str(target),
        "collection_digest": str(collection.digest),
        "sequence_count": len(collection),
        "provenance": _provenance_note("RGSI is a TSV index"),
    }


__all__ = [
    "refget_alias_list",
    "refget_alias_pull",
    "refget_alias_register",
    "refget_alias_remove",
    "refget_alias_reverse_lookup",
    "refget_collection_attribute_array",
    "refget_collection_compare",
    "refget_collection_digests",
    "refget_collection_import",
    "refget_collection_remove",
    "refget_collection_retrieve",
    "refget_collections_by_attribute",
    "refget_collections_list",
    "refget_orphan_removal_plan",
    "refget_sequence_retrieve",
    "refget_sequence_substring",
    "refget_store_add_fasta",
    "refget_store_add_sequence",
    "refget_write_rgsi",
]

_register_public_verbose(globals(), tuple(__all__))
