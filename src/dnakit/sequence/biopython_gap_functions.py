"""Biopython network look-ups and DNA plots that DNAKit did not provide before (package ``Bio``).

Network functions call one upstream Biopython function through its module attribute at call
time, so the remote service is contacted only when a function is called and tests can replace
the upstream callable. They cover NCBI BLAST (``Bio.Blast.NCBIWWW.qblast`` and
``Bio.Blast.qblast``),
NCBI E-utilities (``Bio.Entrez``), KEGG REST (``Bio.KEGG.REST``) and TogoWS (``Bio.TogoWS``).
Each returns a JSON-serialisable dictionary with the service, the upstream function, the request
parameters, the UTC retrieval time, the Biopython version, the output path (or None) and the
parsed data.

Plot functions write one PNG to a new explicit path with the matplotlib object API (no pyplot
global state) and return the path together with the plotted numbers. ``phylo_draw_ascii_text``
writes a text file instead of an image, because an ASCII tree has no image.

Rules that apply to all functions here:

- An NCBI e-mail address is optional. It is set on the upstream module only while the call runs
  and is not written to the result. DNAKit never supplies an address.
- Sizes are bounded: at most 200 NCBI identifiers per request, 10 KEGG entries per
  ``kegg_get_entries`` request, 100 KEGG identifiers per conversion or link request and 1000
  TogoWS identifiers per iterator call.
- Output files must be new paths; an existing file is refused.
- Upstream exceptions (``ValueError`` for unknown identifiers, ``urllib`` errors, NCBI error
  responses) propagate unchanged.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import io
import itertools
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from itertools import accumulate
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_BIO_EXTRA = "bio-formats"
_PLOT_EXTRA = "bio-plots"
_EUTILS = "NCBI Entrez E-utilities (eutils.ncbi.nlm.nih.gov)"
_BLAST = "NCBI BLAST QBLAST (blast.ncbi.nlm.nih.gov)"
_KEGG = "KEGG REST (rest.kegg.jp)"
_TOGOWS = "TogoWS (togows.dbcls.jp)"
_BLAST_PROGRAMS = ("blastn", "blastp", "blastx", "tblastn", "tblastx")
_MAX_QUERY = 100_000
_MAX_TEXT = 4096
_MAX_IDS = 200
_MAX_ESEARCH = 10_000
_MAX_BLAST_HITS = 500
_MAX_KEGG_GET = 10
_MAX_KEGG_IDS = 100
_MAX_TOGOWS_ITER = 1000
_MAX_TOGOWS_IDS = 1_000_000
_MAX_SEQUENCE = 10_000_000
_MAX_TREE_TEXT = 5_000_000
_TOKEN = re.compile(r"[A-Za-z][A-Za-z0-9_.-]{0,63}")
_IDENTIFIER = re.compile(r"[A-Za-z0-9_.:|+\-]{1,256}")
_EMAIL = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")
_DNA_LETTERS = frozenset("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_KEGG_GET_OPTIONS = (None, "aaseq", "ntseq", "mol", "kcf", "kgml", "json")
_KEGG_FIND_OPTIONS = (None, "formula", "exact_mass", "mol_weight")
_KEGG_CONVERT_OPTIONS = (None, "turtle", "n-triple")


# --- shared helpers -------------------------------------------------------------------------------


def _bio(module: str) -> Any:
    return require(module, _BIO_EXTRA)


def _clean_text(value: object, name: str, *, max_length: int = _MAX_TEXT) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    text = value.strip()
    if len(text) > max_length:
        raise ConfigurationError(f"{name} is longer than {max_length} characters.")
    return text


def _token(value: object, name: str) -> str:
    if not isinstance(value, str) or not _TOKEN.fullmatch(value):
        raise ConfigurationError(f"{name} must be a short name such as 'pubmed', 'nt' or 'hsa'.")
    return value


def _identifier(value: object, name: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise ConfigurationError(f"{name} must be an identifier without spaces, commas or slashes.")
    return value


def _identifiers(values: object, name: str, *, limit: int = _MAX_IDS) -> list[str]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of identifier strings.")
    items = [_identifier(item, name) for item in values]
    if not items or len(items) > limit:
        raise ConfigurationError(f"{name} must hold 1 to {limit} identifiers.")
    return items


def _bounded_int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _positive_float(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive number.")
    return float(value)


def _email(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not _EMAIL.fullmatch(value.strip()):
        raise ConfigurationError("email must be an e-mail address when it is given.")
    return value.strip()


def _new_output(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _plain(value: Any) -> Any:
    """Convert Bio.Entrez parser elements (dict, list and str subclasses) to plain JSON types."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value)
    return str(value)


def _read_text(handle: Any) -> str:
    """Read a text or binary handle (or a string) completely and close it."""
    if isinstance(handle, str):
        return handle
    if isinstance(handle, bytes):
        return handle.decode("utf-8", errors="replace")
    try:
        raw = handle.read()
    finally:
        close = getattr(handle, "close", None)
        if close is not None:
            close()
    return raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else str(raw)


@contextlib.contextmanager
def _temporary_attribute(module: Any, name: str, value: str | None) -> Iterator[None]:
    """Set a module setting such as ``Bio.Entrez.email`` only while one upstream call runs."""
    if value is None:
        yield
        return
    previous = getattr(module, name, None)
    setattr(module, name, value)
    try:
        yield
    finally:
        setattr(module, name, previous)


def _envelope(
    service: str,
    upstream: str,
    parameters: Mapping[str, Any],
    data: Any,
    *,
    output_path: str | None = None,
) -> dict[str, Any]:
    return {
        "status": "completed",
        "service": service,
        "upstream": upstream,
        "parameters": _plain(dict(parameters)),
        "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "biopython_version": str(_bio("Bio").__version__),
        "output_path": output_path,
        "data": data,
    }


# --- NCBI Entrez E-utilities ----------------------------------------------------------------------


def _entrez(
    upstream: str,
    params: Mapping[str, Any],
    email: str | None,
    consume: Callable[[Any, Any], Any],
) -> Any:
    entrez = _bio("Bio.Entrez")
    call = getattr(entrez, upstream)
    with _temporary_attribute(entrez, "email", email):
        handle = call(**params)
        return consume(entrez, handle)


def _parsed(entrez: Any, handle: Any) -> Any:
    try:
        return _plain(entrez.read(handle))
    finally:
        handle.close()


def _text_of(entrez: Any, handle: Any) -> str:
    return _read_text(handle)


def ncbi_einfo(
    *, db: str | None = None, email: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List the NCBI Entrez databases, or describe one database, with Bio.Entrez.einfo."""
    params: dict[str, Any] = {}
    if db is not None:
        params["db"] = _token(db, "db")
    data = _entrez("einfo", params, _email(email), _parsed)
    return _envelope(_EUTILS, "Bio.Entrez.einfo", params, data)


def ncbi_esearch(
    db: str,
    term: str,
    *,
    retmax: int = 20,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search an NCBI Entrez database; return the total count and the UIDs (Bio.Entrez.esearch)."""
    params = {
        "db": _token(db, "db"),
        "term": _clean_text(term, "term"),
        "retmax": _bounded_int(retmax, "retmax", 0, _MAX_ESEARCH),
    }
    parsed = _entrez("esearch", params, _email(email), _parsed)
    ids = [str(item) for item in parsed.get("IdList", [])]
    data = {**parsed, "count": int(parsed.get("Count", 0)), "ids": ids}
    return _envelope(_EUTILS, "Bio.Entrez.esearch", params, data)


def ncbi_efetch(
    db: str,
    ids: Sequence[str],
    *,
    rettype: str = "fasta",
    retmode: str = "text",
    output_path: str | None = None,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch NCBI Entrez records as text with Bio.Entrez.efetch, optionally into a new file."""
    params = {
        "db": _token(db, "db"),
        "id": ",".join(_identifiers(ids, "ids")),
        "rettype": _token(rettype, "rettype"),
        "retmode": _token(retmode, "retmode"),
    }
    target = _new_output(output_path, "output_path") if output_path is not None else None
    text = _entrez("efetch", params, _email(email), _text_of)
    written = None
    if target is not None:
        target.write_text(text, encoding="utf-8")
        written = str(target)
    data = {"text": text, "characters": len(text)}
    return _envelope(_EUTILS, "Bio.Entrez.efetch", params, data, output_path=written)


def ncbi_esummary(
    db: str, ids: Sequence[str], *, email: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Retrieve NCBI Entrez document summaries for identifiers with Bio.Entrez.esummary."""
    params = {"db": _token(db, "db"), "id": ",".join(_identifiers(ids, "ids"))}
    data = _entrez("esummary", params, _email(email), _parsed)
    return _envelope(_EUTILS, "Bio.Entrez.esummary", params, data)


def ncbi_elink(
    dbfrom: str,
    db: str,
    ids: Sequence[str],
    *,
    linkname: str | None = None,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find NCBI Entrez records linked to identifiers in another database (Bio.Entrez.elink)."""
    params: dict[str, Any] = {
        "dbfrom": _token(dbfrom, "dbfrom"),
        "db": _token(db, "db"),
        "id": ",".join(_identifiers(ids, "ids")),
    }
    if linkname is not None:
        params["linkname"] = _token(linkname, "linkname")
    data = _entrez("elink", params, _email(email), _parsed)
    return _envelope(_EUTILS, "Bio.Entrez.elink", params, data)


def ncbi_epost(
    db: str, ids: Sequence[str], *, email: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Post identifiers to the NCBI history server (Bio.Entrez.epost); return the WebEnv key."""
    params = {"db": _token(db, "db"), "id": ",".join(_identifiers(ids, "ids"))}
    data = _entrez("epost", params, _email(email), _parsed)
    return _envelope(_EUTILS, "Bio.Entrez.epost", params, data)


# --- NCBI BLAST QBLAST ----------------------------------------------------------------------------


def _blast_summary(record: Any) -> dict[str, Any]:
    alignments = []
    for alignment in record.alignments:
        hsps = [
            {
                "score": int(hsp.score),
                "bits": float(hsp.bits),
                "expect": float(hsp.expect),
                "identities": int(hsp.identities),
                "align_length": int(hsp.align_length),
                "query_start": int(hsp.query_start),
                "query_end": int(hsp.query_end),
                "subject_start": int(hsp.sbjct_start),
                "subject_end": int(hsp.sbjct_end),
            }
            for hsp in alignment.hsps
        ]
        alignments.append(
            {
                "hit_id": str(alignment.hit_id),
                "hit_def": str(alignment.hit_def),
                "length": int(alignment.length),
                "hsps": hsps,
            }
        )
    query_id = getattr(record, "query_id", None)
    return {
        "query_id": None if query_id is None else str(query_id),
        "query": str(record.query),
        "query_letters": None if record.query_letters is None else int(record.query_letters),
        "database": str(record.database),
        "alignment_count": len(alignments),
        "alignments": alignments,
    }


def _run_qblast(
    module_name: str,
    function_name: str,
    program: str,
    database: str,
    sequence: str,
    *,
    hitlist_size: int,
    expect: float,
    entrez_query: str,
    output_xml: str | None,
    email: str | None,
) -> dict[str, Any]:
    if program not in _BLAST_PROGRAMS:
        raise ConfigurationError(f"program must be one of {', '.join(_BLAST_PROGRAMS)}.")
    query = _clean_text(sequence, "sequence", max_length=_MAX_QUERY)
    params = {
        "program": program,
        "database": _token(database, "database"),
        "sequence_length": len(query),
        "hitlist_size": _bounded_int(hitlist_size, "hitlist_size", 1, _MAX_BLAST_HITS),
        "expect": _positive_float(expect, "expect"),
        "entrez_query": _clean_text(entrez_query, "entrez_query", max_length=512),
        "format_type": "XML",
    }
    target = _new_output(output_xml, "output_xml") if output_xml is not None else None
    module = _bio(module_name)
    upstream = getattr(module, function_name)
    with _temporary_attribute(module, "email", _email(email)):
        handle = upstream(
            program=program,
            database=params["database"],
            sequence=query,
            hitlist_size=params["hitlist_size"],
            expect=params["expect"],
            entrez_query=params["entrez_query"],
            format_type="XML",
        )
        text = _read_text(handle)
    written = None
    if target is not None:
        target.write_text(text, encoding="utf-8")
        written = str(target)
    record = _bio("Bio.Blast.NCBIXML").read(io.StringIO(text))
    return _envelope(
        _BLAST,
        f"{module_name}.{function_name}",
        params,
        _blast_summary(record),
        output_path=written,
    )


def ncbiwww_qblast_search(
    program: str,
    database: str,
    sequence: str,
    *,
    hitlist_size: int = 50,
    expect: float = 10.0,
    entrez_query: str = "(none)",
    output_xml: str | None = None,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run NCBI BLAST QBLAST (Bio.Blast.NCBIWWW.qblast) and summarise the hits."""
    return _run_qblast(
        "Bio.Blast.NCBIWWW",
        "qblast",
        program,
        database,
        sequence,
        hitlist_size=hitlist_size,
        expect=expect,
        entrez_query=entrez_query,
        output_xml=output_xml,
        email=email,
    )


def blast_qblast_search(
    program: str,
    database: str,
    sequence: str,
    *,
    hitlist_size: int = 50,
    expect: float = 10.0,
    entrez_query: str = "(none)",
    output_xml: str | None = None,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run NCBI BLAST QBLAST (Bio.Blast.qblast) with a DNA query and summarise the hits."""
    return _run_qblast(
        "Bio.Blast",
        "qblast",
        program,
        database,
        sequence,
        hitlist_size=hitlist_size,
        expect=expect,
        entrez_query=entrez_query,
        output_xml=output_xml,
        email=email,
    )


# --- KEGG REST ------------------------------------------------------------------------------------


def _kegg_argument(value: object, name: str, *, limit: int = _MAX_KEGG_IDS) -> str | list[str]:
    """Accept a KEGG database name or entry string, or a list of identifiers (bounded)."""
    if isinstance(value, str):
        return _identifier(value, name)
    return _identifiers(value, name, limit=limit)


def _kegg_call(function_name: str, *args: Any) -> str:
    rest = _bio("Bio.KEGG.REST")
    return _read_text(getattr(rest, function_name)(*args))


def _kegg_rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line.strip()]


def kegg_database_info(database: str, *, verbose: bool = True) -> dict[str, Any]:
    """Describe a KEGG database (release and statistics) with Bio.KEGG.REST.kegg_info."""
    name = _identifier(database, "database")
    text = _kegg_call("kegg_info", name)
    data = {"text": text, "line_count": len(text.splitlines())}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_info", {"database": name}, data)


def kegg_list_entries(
    database: str | Sequence[str],
    *,
    organism: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List KEGG entries of a database or an organism's pathways (Bio.KEGG.REST.kegg_list)."""
    name = _kegg_argument(database, "database")
    org = None if organism is None else _token(organism, "organism")
    rows = _kegg_rows(_kegg_call("kegg_list", name, org))
    data = {"rows": rows, "row_count": len(rows)}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_list", {"database": name, "organism": org}, data)


def kegg_find_entries(
    database: str,
    query: str,
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a KEGG database by keyword, formula or mass (Bio.KEGG.REST.kegg_find)."""
    name = _token(database, "database")
    text_query = _clean_text(query, "query", max_length=512)
    if option not in _KEGG_FIND_OPTIONS:
        raise ConfigurationError(
            f"option must be one of {', '.join(str(o) for o in _KEGG_FIND_OPTIONS)}."
        )
    rows = _kegg_rows(_kegg_call("kegg_find", name, text_query, option))
    data = {"rows": rows, "row_count": len(rows)}
    params = {"database": name, "query": text_query, "option": option}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_find", params, data)


def kegg_convert_identifiers(
    target_db: str,
    source: str | Sequence[str],
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert KEGG identifiers to or from outside identifiers with Bio.KEGG.REST.kegg_conv."""
    target = _token(target_db, "target_db")
    source_value = _kegg_argument(source, "source")
    if option not in _KEGG_CONVERT_OPTIONS:
        raise ConfigurationError(
            f"option must be one of {', '.join(str(o) for o in _KEGG_CONVERT_OPTIONS)}."
        )
    rows = _kegg_rows(_kegg_call("kegg_conv", target, source_value, option))
    data = {"rows": rows, "row_count": len(rows)}
    params = {"target_db": target, "source": source_value, "option": option}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_conv", params, data)


def kegg_link_entries(
    target_db: str,
    source: str | Sequence[str],
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find KEGG entries linked to a database or to given entries (Bio.KEGG.REST.kegg_link)."""
    target = _token(target_db, "target_db")
    source_value = _kegg_argument(source, "source")
    if option not in _KEGG_CONVERT_OPTIONS:
        raise ConfigurationError(
            f"option must be one of {', '.join(str(o) for o in _KEGG_CONVERT_OPTIONS)}."
        )
    rows = _kegg_rows(_kegg_call("kegg_link", target, source_value, option))
    data = {"rows": rows, "row_count": len(rows)}
    params = {"target_db": target, "source": source_value, "option": option}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_link", params, data)


def kegg_get_entries(
    entries: str | Sequence[str],
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve KEGG entries as text (Bio.KEGG.REST.kegg_get); option selects the format."""
    if isinstance(entries, str):
        entry_value: str | list[str] = _identifier(entries, "entries")
    else:
        entry_value = _identifiers(entries, "entries", limit=_MAX_KEGG_GET)
    if option not in _KEGG_GET_OPTIONS:
        raise ConfigurationError(
            f"option must be one of {', '.join(str(o) for o in _KEGG_GET_OPTIONS)}; "
            "image is not supported."
        )
    text = _kegg_call("kegg_get", entry_value, option)
    data = {"text": text, "line_count": len(text.splitlines())}
    params = {"entries": entry_value, "option": option}
    return _envelope(_KEGG, "Bio.KEGG.REST.kegg_get", params, data)


# --- TogoWS ---------------------------------------------------------------------------------------


def _togows() -> Any:
    return _bio("Bio.TogoWS")


def togows_fetch_entry(
    database: str,
    identifier: str,
    *,
    output_format: str | None = None,
    field: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch one TogoWS record, optionally in one format or field (Bio.TogoWS.entry)."""
    db = _token(database, "database")
    ident = _identifier(identifier, "identifier")
    fmt = None if output_format is None else _token(output_format, "output_format")
    fld = None if field is None else _token(field, "field")
    handle = _togows().entry(db=db, id=ident, format=fmt, field=fld)
    text = _read_text(handle)
    params: dict[str, Any] = {"database": db, "identifier": ident}
    if fmt is not None:
        params["output_format"] = fmt
    if fld is not None:
        params["field"] = fld
    return _envelope(_TOGOWS, "Bio.TogoWS.entry", params, {"text": text, "characters": len(text)})


def togows_search_ids(
    database: str,
    query: str,
    *,
    offset: int | None = None,
    limit: int | None = None,
    output_format: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a TogoWS database; return one page of identifiers, or raw text in a format."""
    db = _token(database, "database")
    text_query = _clean_text(query, "query", max_length=1024)
    start = None if offset is None else _bounded_int(offset, "offset", 1, 1_000_000_000)
    count = None if limit is None else _bounded_int(limit, "limit", 1, _MAX_TOGOWS_IDS)
    fmt = None if output_format is None else _token(output_format, "output_format")
    text = _read_text(
        _togows().search(db=db, query=text_query, offset=start, limit=count, format=fmt)
    )
    params: dict[str, Any] = {"database": db, "query": text_query, "offset": start, "limit": count}
    if fmt is not None:
        params["output_format"] = fmt
        data: dict[str, Any] = {"text": text}
    else:
        ids = text.split()
        data = {"ids": ids, "count": len(ids)}
    return _envelope(_TOGOWS, "Bio.TogoWS.search", params, data)


def togows_count_matches(database: str, query: str, *, verbose: bool = True) -> dict[str, Any]:
    """Count the TogoWS matches of a query with Bio.TogoWS.search_count."""
    db = _token(database, "database")
    text_query = _clean_text(query, "query", max_length=1024)
    count = int(_togows().search_count(db=db, query=text_query))
    return _envelope(
        _TOGOWS, "Bio.TogoWS.search_count", {"database": db, "query": text_query}, {"count": count}
    )


def togows_iterate_ids(
    database: str,
    query: str,
    *,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect up to a bounded number of TogoWS identifiers via Bio.TogoWS.search_iter."""
    db = _token(database, "database")
    text_query = _clean_text(query, "query", max_length=1024)
    bound = _bounded_int(limit, "limit", 1, _MAX_TOGOWS_ITER)
    iterator = _togows().search_iter(db=db, query=text_query, limit=bound)
    ids = [str(identifier) for identifier in itertools.islice(iterator, bound)]
    data = {"ids": ids, "returned_count": len(ids)}
    return _envelope(
        _TOGOWS,
        "Bio.TogoWS.search_iter",
        {"database": db, "query": text_query, "limit": bound},
        data,
    )


def togows_convert_data(
    data: str,
    in_format: str,
    out_format: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a record between file formats with TogoWS (Bio.TogoWS.convert)."""
    payload = _clean_text(data, "data", max_length=_MAX_SEQUENCE)
    source_format = _token(in_format, "in_format")
    target_format = _token(out_format, "out_format")
    text = _read_text(_togows().convert(payload, source_format, target_format))
    params = {"in_format": source_format, "out_format": target_format, "characters": len(payload)}
    return _envelope(_TOGOWS, "Bio.TogoWS.convert", params, {"text": text, "characters": len(text)})


# --- Phylo trees and GC skew (matplotlib PNG output) ----------------------------------------------


def _plot(module: str) -> Any:
    return require(module, _PLOT_EXTRA)


def _new_png(path: str | Path, name: str) -> Path:
    target = _new_output(path, name)
    if target.suffix.lower() != ".png":
        raise ConfigurationError(f"{name} must end with .png.")
    return target


def phylo_draw_newick_png(
    newick: str,
    image_path: str,
    *,
    show_confidence: bool = True,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a Newick tree with Bio.Phylo.draw into a new PNG; return the tree data."""
    text = _clean_text(newick, "newick", max_length=_MAX_TREE_TEXT)
    resolution = _bounded_int(dpi, "dpi", 50, 600)
    target = _new_png(image_path, "image_path")
    phylo = _bio("Bio.Phylo")
    figure_module = _plot("matplotlib.figure")
    canvas_module = _plot("matplotlib.backends.backend_agg")
    tree = phylo.read(io.StringIO(text), "newick")
    figure = figure_module.Figure(figsize=(6, 4))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(1, 1, 1)
    phylo.draw(tree, do_show=False, axes=axes, show_confidence=bool(show_confidence))
    figure.savefig(target, format="png", dpi=resolution)
    normalized = io.StringIO()
    phylo.write(tree, normalized, "newick")
    terminals = tree.get_terminals()
    data = {
        "newick": normalized.getvalue().strip(),
        "tip_names": [str(tip.name) for tip in terminals],
        "terminal_count": len(terminals),
        "clade_count": sum(1 for _ in tree.find_clades()),
        "total_branch_length": float(tree.total_branch_length()),
    }
    return {"status": "completed", "image_path": str(target), "format": "png", "data": data}


def phylo_draw_ascii_text(
    newick: str,
    output_path: str,
    *,
    column_width: int = 80,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a Newick tree as ASCII art (Bio.Phylo.draw_ascii); write and return the text."""
    text = _clean_text(newick, "newick", max_length=_MAX_TREE_TEXT)
    width = _bounded_int(column_width, "column_width", 20, 1000)
    target = _new_output(output_path, "output_path")
    phylo = _bio("Bio.Phylo")
    tree = phylo.read(io.StringIO(text), "newick")
    buffer = io.StringIO()
    phylo.draw_ascii(tree, file=buffer, column_width=width)
    drawing = buffer.getvalue()
    target.write_text(drawing, encoding="utf-8")
    data = {
        "text": drawing,
        "line_count": len(drawing.splitlines()),
        "tip_names": [str(tip.name) for tip in tree.get_terminals()],
        "column_width": width,
    }
    return {"status": "completed", "output_path": str(target), "data": data}


def gc_skew_window_plot(
    sequence: str,
    image_path: str,
    *,
    window: int = 1000,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the GC skew per window and its running sum; save a line plot as a new PNG.

    The skew is (G - C) / (G + C) per window. The window values come from Bio.SeqUtils.GC_skew,
    the numeric part of Bio.SeqUtils.xGC_skew. xGC_skew draws on a Tk canvas and returns nothing,
    so it needs a display; this function writes a PNG and returns the numbers.
    """
    text = _clean_text(sequence, "sequence", max_length=_MAX_SEQUENCE)
    if set(text) - _DNA_LETTERS:
        raise ConfigurationError("sequence must contain only IUPAC DNA letters.")
    size = _bounded_int(window, "window", 1, _MAX_SEQUENCE)
    resolution = _bounded_int(dpi, "dpi", 50, 600)
    target = _new_png(image_path, "image_path")
    skew = [float(value) for value in _bio("Bio.SeqUtils").GC_skew(text, size)]
    accumulated = list(accumulate(skew))
    starts = [index * size for index in range(len(skew))]
    figure_module = _plot("matplotlib.figure")
    canvas_module = _plot("matplotlib.backends.backend_agg")
    figure = figure_module.Figure(figsize=(7, 3.5))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(1, 1, 1)
    axes.axhline(0.0, color="grey", linewidth=0.8)
    axes.plot(starts, skew, color="tab:blue", label="GC skew")
    axes.plot(starts, accumulated, color="tab:purple", label="accumulated GC skew")
    axes.set_xlabel("window start (nt)")
    axes.set_ylabel("(G - C) / (G + C)")
    axes.set_title(f"GC skew, window {size} nt")
    axes.legend()
    figure.tight_layout()
    figure.savefig(target, format="png", dpi=resolution)
    data = {
        "sequence_length": len(text),
        "window": size,
        "window_count": len(skew),
        "window_starts": starts,
        "gc_skew": skew,
        "accumulated_gc_skew": accumulated,
    }
    return {"status": "completed", "image_path": str(target), "format": "png", "data": data}


__all__ = [
    "blast_qblast_search",
    "gc_skew_window_plot",
    "kegg_convert_identifiers",
    "kegg_database_info",
    "kegg_find_entries",
    "kegg_get_entries",
    "kegg_link_entries",
    "kegg_list_entries",
    "ncbi_efetch",
    "ncbi_einfo",
    "ncbi_elink",
    "ncbi_epost",
    "ncbi_esearch",
    "ncbi_esummary",
    "ncbiwww_qblast_search",
    "phylo_draw_ascii_text",
    "phylo_draw_newick_png",
    "togows_convert_data",
    "togows_count_matches",
    "togows_fetch_entry",
    "togows_iterate_ids",
    "togows_search_ids",
]

_register_public_verbose(globals())
