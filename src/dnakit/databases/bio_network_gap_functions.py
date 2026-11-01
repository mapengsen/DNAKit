"""Biopython remote-service look-ups for NCBI, KEGG, TogoWS, UniProt and ExPASy.

Each public function calls one Biopython network function (``Bio.Entrez``, ``Bio.KEGG.REST``,
``Bio.TogoWS``, ``Bio.UniProt``, ``Bio.ExPASy``, ``Bio.Blast.NCBIWWW``) with the same arguments the
original takes, and returns a JSON-serialisable dictionary. The remote service is contacted at call
time; nothing is cached or stored. Every result records the service name and the access date
(``accessed``, ISO 8601, local date of the call) and the Biopython version used.

Rules that apply to all functions here:

- NCBI E-utilities (``entrez_*``) require the caller's e-mail address as the keyword-only ``email``
  argument. There is no default and DNAKit never supplies an address. The address is sent to NCBI as
  the ``email`` parameter and is not echoed in results.
- Biopython exceptions (``ValueError`` for unknown identifiers, ``urllib`` errors, Entrez error
  responses) propagate unchanged, so callers see the same failure the original would raise.
- ``ncbi_qblast_search`` is a long-running remote job; it has a documented ``timeout`` and raises
  ``BackendTimeoutError`` when the job has not returned in time. NCBI keeps running the job after
  the timeout; the caller cannot cancel it from here.
- Network use happens only when a function is called. Importing this module makes no requests.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import io
import re
import threading
from collections.abc import Callable, Iterator, Mapping, Sequence
from itertools import islice
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendTimeoutError, ConfigurationError

_BIO_EXTRA = "bio-formats"
_TOOL_NAME = "dnakit"
_EUTILS_SERVICE = "NCBI Entrez E-utilities (eutils.ncbi.nlm.nih.gov)"
_BLAST_SERVICE = "NCBI BLAST QBLAST URL API (blast.ncbi.nlm.nih.gov)"
_KEGG_SERVICE = "KEGG REST API (rest.kegg.jp)"
_TOGOWS_SERVICE = "TogoWS REST (togows.dbcls.jp)"
_UNIPROT_SERVICE = "UniProt REST search (rest.uniprot.org)"
_EXPASY_SERVICE = "UniProtKB flat-file retrieval (www.uniprot.org; Bio.ExPASy.get_sprot_raw)"
_KEGG_CONV_OPTIONS = frozenset({None, "turtle", "n-triple"})
_KEGG_FIND_OPTIONS = frozenset({None, "formula", "exact_mass", "mol_weight"})
_KEGG_GET_OPTIONS = frozenset({None, "aaseq", "ntseq", "mol", "kcf", "kgml", "json"})
_MAX_IDS = 200
_MAX_KEGG_GET = 10
_MAX_KEGG_LIST = 100
_EMAIL_PATTERN = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")
_ENTREZ_LOCK = threading.RLock()
_HSP_FIELDS = (
    "score",
    "bits",
    "expect",
    "num_alignments",
    "identities",
    "positives",
    "gaps",
    "align_length",
    "frame",
    "query_start",
    "query_end",
    "sbjct_start",
    "sbjct_end",
    "query",
    "match",
    "sbjct",
)


def _bio(module: str) -> Any:
    """Import a Biopython module, raising BackendUnavailableError with the extra's hint."""
    return require(module, _BIO_EXTRA)


def _envelope(
    operation: str,
    service: str,
    call: str,
    parameters: Mapping[str, Any],
    result: Any,
) -> dict[str, Any]:
    bio = _bio("Bio")
    return {
        "status": "ok",
        "operation": operation,
        "service": service,
        "accessed": _dt.date.today().isoformat(),
        "backend": {"name": "Biopython", "version": str(bio.__version__), "call": call},
        "parameters": dict(parameters),
        "result": result,
    }


def _text(value: object, name: str, *, max_length: int = 4096) -> str:
    """Require non-empty text without NUL; the value itself is passed on unchanged."""
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be non-empty text.", code="INVALID_QUERY")
    if len(value) > max_length or "\x00" in value:
        raise ConfigurationError(
            f"{name} exceeds {max_length} characters or contains NUL.", code="INVALID_QUERY"
        )
    return value


def _bounded_int(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.", code="INVALID_QUERY")
    if not minimum <= value <= maximum:
        raise ConfigurationError(
            f"{name} must be between {minimum} and {maximum}.", code="INVALID_QUERY"
        )
    return value


def _ids(value: object, name: str) -> str:
    """Return E-utilities/TogoWS identifiers as one comma-separated string."""
    items = [value] if isinstance(value, str) else value
    if not isinstance(items, (list, tuple)) or not items:
        raise ConfigurationError(
            f"{name} must be a string or a non-empty list of strings.", code="INVALID_QUERY"
        )
    if len(items) > _MAX_IDS:
        raise ConfigurationError(
            f"{name} accepts at most {_MAX_IDS} identifiers per call.", code="INVALID_QUERY"
        )
    cleaned: list[str] = []
    for item in items:
        identifier = _text(item, f"each item of {name}", max_length=256)
        if "," in identifier or any(ch.isspace() for ch in identifier):
            raise ConfigurationError(
                f"identifiers in {name} must not contain commas or whitespace.",
                code="INVALID_QUERY",
            )
        cleaned.append(identifier)
    return ",".join(cleaned)


def _check_email(email: object) -> str:
    """Require the caller's e-mail address for NCBI E-utilities; no default is supplied."""
    if not isinstance(email, str) or _EMAIL_PATTERN.fullmatch(email.strip()) is None:
        raise ConfigurationError(
            "NCBI E-utilities require the caller's e-mail address as email=...; "
            "DNAKit does not supply one.",
            code="INVALID_QUERY",
            hint="Pass email='you@your-institution.example' to every entrez_* call.",
        )
    return email.strip()


def _option(value: str | None, allowed: frozenset[str | None], name: str) -> str | None:
    if value not in allowed:
        choices = ", ".join(sorted(item for item in allowed if item is not None))
        raise ConfigurationError(f"{name} must be None or one of: {choices}.", code="INVALID_QUERY")
    return value


def _plain(value: Any) -> Any:
    """Convert Biopython parser objects (dict/list/str subclasses) to builtin JSON types."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return str(value)
    if value is None or isinstance(value, (bool, int, float)):
        return value
    return str(value)


def _read_text(handle: Any) -> str:
    try:
        data = handle.read()
    finally:
        with contextlib.suppress(Exception):
            handle.close()
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return str(data)


def _rows(text: str) -> list[list[str]]:
    """Split KEGG tab-separated output into rows of fields; blank lines are dropped."""
    return [line.split("\t") for line in text.splitlines() if line.strip()]


@contextlib.contextmanager
def _entrez(email: str, api_key: str | None) -> Iterator[Any]:
    """Set Bio.Entrez's e-mail, tool and API key for one call, then restore the previous values."""
    entrez = _bio("Bio.Entrez")
    with _ENTREZ_LOCK:
        saved = (entrez.email, entrez.tool, entrez.api_key)
        entrez.email, entrez.tool, entrez.api_key = email, _TOOL_NAME, api_key
        try:
            yield entrez
        finally:
            entrez.email, entrez.tool, entrez.api_key = saved


def entrez_search(
    database: str,
    term: str,
    *,
    email: str,
    retmax: int = 20,
    retstart: int = 0,
    sort: str | None = None,
    field: str | None = None,
    usehistory: bool = False,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run an NCBI ESearch query in any Entrez database and return the ID list and count.

    Calls ``Bio.Entrez.esearch(db=..., term=...)``. ``email`` is required. ``usehistory=True`` adds
    the WebEnv/QueryKey pair to the result. The result is the parsed ESearch record (``Count``,
    ``RetMax``, ``IdList``, ``WebEnv``, ``QueryKey``, ``QueryTranslation`` and so on).
    """
    db = _text(database, "database", max_length=64)
    query = _text(term, "term")
    mail = _check_email(email)
    options: dict[str, Any] = {
        "retmax": _bounded_int(retmax, "retmax", minimum=0, maximum=10_000),
        "retstart": _bounded_int(retstart, "retstart", minimum=0, maximum=10_000_000),
    }
    if sort is not None:
        options["sort"] = _text(sort, "sort", max_length=64)
    if field is not None:
        options["field"] = _text(field, "field", max_length=64)
    if usehistory:
        options["usehistory"] = "y"
    with _entrez(mail, api_key) as entrez:
        handle = entrez.esearch(db=db, term=query, **options)
        record = entrez.read(handle)
    return _envelope(
        "entrez_search",
        _EUTILS_SERVICE,
        "Bio.Entrez.esearch",
        {"database": db, "term": query, **options},
        _plain(record),
    )


def entrez_info(
    *,
    email: str,
    database: str | None = None,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the Entrez databases (no database) or describe one database (record counts, fields).

    Calls ``Bio.Entrez.einfo()`` or ``Bio.Entrez.einfo(db=...)``. ``email`` is required.
    """
    mail = _check_email(email)
    options: dict[str, Any] = {}
    if database is not None:
        options["db"] = _text(database, "database", max_length=64)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.einfo(**options)
        record = entrez.read(handle)
    return _envelope(
        "entrez_info",
        _EUTILS_SERVICE,
        "Bio.Entrez.einfo",
        dict(options),
        _plain(record),
    )


def entrez_fetch(
    database: str,
    ids: str | Sequence[str],
    *,
    email: str,
    rettype: str | None = None,
    retmode: str | None = None,
    retstart: int | None = None,
    retmax: int | None = None,
    strand: int | None = None,
    seq_start: int | None = None,
    seq_stop: int | None = None,
    parse_xml: bool = False,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch Entrez records by ID in any database (EFetch), for example FASTA or GenBank text.

    Calls ``Bio.Entrez.efetch(db=..., id=...)`` with the given ``rettype``/``retmode``/range
    arguments. ``email`` is required. The result holds the decoded response text, its size, and,
    when ``parse_xml=True`` and the response is XML, the Biopython-parsed tree (``parsed``).
    """
    db = _text(database, "database", max_length=64)
    id_string = _ids(ids, "ids")
    mail = _check_email(email)
    options: dict[str, Any] = {}
    if rettype is not None:
        options["rettype"] = _text(rettype, "rettype", max_length=32)
    if retmode is not None:
        options["retmode"] = _text(retmode, "retmode", max_length=32)
    if retstart is not None:
        options["retstart"] = _bounded_int(retstart, "retstart", minimum=0, maximum=10_000_000)
    if retmax is not None:
        options["retmax"] = _bounded_int(retmax, "retmax", minimum=0, maximum=10_000)
    if strand is not None:
        options["strand"] = _bounded_int(strand, "strand", minimum=1, maximum=2)
    if seq_start is not None:
        options["seq_start"] = _bounded_int(seq_start, "seq_start", minimum=1, maximum=2**40)
    if seq_stop is not None:
        options["seq_stop"] = _bounded_int(seq_stop, "seq_stop", minimum=1, maximum=2**40)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.efetch(db=db, id=id_string, **options)
        data = handle.read()
        with contextlib.suppress(Exception):
            handle.close()
        raw = data if isinstance(data, bytes) else str(data).encode("utf-8")
        result: dict[str, Any] = {
            "text": raw.decode("utf-8", errors="replace"),
            "content_bytes": len(raw),
            "rettype": rettype,
            "retmode": retmode,
        }
        if parse_xml:
            result["parsed"] = _plain(entrez.read(io.BytesIO(raw)))
    return _envelope(
        "entrez_fetch",
        _EUTILS_SERVICE,
        "Bio.Entrez.efetch",
        {"database": db, "ids": id_string.split(","), **options},
        result,
    )


def entrez_summary(
    database: str,
    ids: str | Sequence[str],
    *,
    email: str,
    version: str | None = None,
    retstart: int | None = None,
    retmax: int | None = None,
    webenv: str | None = None,
    query_key: str | None = None,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return ESummary document summaries for Entrez IDs (any database).

    Calls ``Bio.Entrez.esummary(db=..., id=...)``. Pass ``version="2.0"`` for the XML-style summary
    layout that the original library defaults to in its examples. ``email`` is required.
    """
    db = _text(database, "database", max_length=64)
    id_string = _ids(ids, "ids")
    mail = _check_email(email)
    options: dict[str, Any] = {}
    if version is not None:
        options["version"] = _text(version, "version", max_length=8)
    if retstart is not None:
        options["retstart"] = _bounded_int(retstart, "retstart", minimum=0, maximum=10_000_000)
    if retmax is not None:
        options["retmax"] = _bounded_int(retmax, "retmax", minimum=0, maximum=10_000)
    if webenv is not None:
        options["webenv"] = _text(webenv, "webenv", max_length=512)
    if query_key is not None:
        options["query_key"] = _text(query_key, "query_key", max_length=64)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.esummary(db=db, id=id_string, **options)
        record = entrez.read(handle)
    return _envelope(
        "entrez_summary",
        _EUTILS_SERVICE,
        "Bio.Entrez.esummary",
        {"database": db, "ids": id_string.split(","), **options},
        _plain(record),
    )


def entrez_link(
    dbfrom: str,
    db: str,
    ids: str | Sequence[str],
    *,
    email: str,
    cmd: str | None = None,
    linkname: str | None = None,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find linked records across Entrez databases (ELink), for example nuccore to gene.

    Calls ``Bio.Entrez.elink(dbfrom=..., db=..., id=...)``. ``cmd`` (for example ``"neighbor"``) and
    ``linkname`` (for example ``"nuccore_gene"``) are passed through when given. ``email`` is
    required.
    """
    source = _text(dbfrom, "dbfrom", max_length=64)
    target = _text(db, "db", max_length=64)
    id_string = _ids(ids, "ids")
    mail = _check_email(email)
    options: dict[str, Any] = {}
    if cmd is not None:
        options["cmd"] = _text(cmd, "cmd", max_length=64)
    if linkname is not None:
        options["linkname"] = _text(linkname, "linkname", max_length=128)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.elink(dbfrom=source, db=target, id=id_string, **options)
        record = entrez.read(handle)
    return _envelope(
        "entrez_link",
        _EUTILS_SERVICE,
        "Bio.Entrez.elink",
        {"dbfrom": source, "db": target, "ids": id_string.split(","), **options},
        _plain(record),
    )


def entrez_post(
    database: str,
    ids: str | Sequence[str],
    *,
    email: str,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Upload a list of IDs to the Entrez History server (EPost) and return WebEnv and QueryKey.

    Calls ``Bio.Entrez.epost(db=..., id=...)``. The returned pair can be passed to
    ``entrez_summary``, ``entrez_fetch`` or ``entrez_link`` in later calls. ``email`` is required.
    """
    db = _text(database, "database", max_length=64)
    id_string = _ids(ids, "ids")
    mail = _check_email(email)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.epost(db=db, id=id_string)
        record = entrez.read(handle)
    return _envelope(
        "entrez_post",
        _EUTILS_SERVICE,
        "Bio.Entrez.epost",
        {"database": db, "ids": id_string.split(",")},
        _plain(record),
    )


def entrez_spell(
    database: str,
    term: str,
    *,
    email: str,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return NCBI's spelling suggestion for an Entrez query term (ESpell).

    Calls ``Bio.Entrez.espell(db=..., term=...)``. ``CorrectedQuery`` is empty when the service
    suggests no change. ``email`` is required.
    """
    db = _text(database, "database", max_length=64)
    query = _text(term, "term")
    mail = _check_email(email)
    with _entrez(mail, api_key) as entrez:
        handle = entrez.espell(db=db, term=query)
        record = entrez.read(handle)
    return _envelope(
        "entrez_spell",
        _EUTILS_SERVICE,
        "Bio.Entrez.espell",
        {"database": db, "term": query},
        _plain(record),
    )


def kegg_rest_info(database: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return KEGG's release and statistics text for one database or organism (``/info``).

    Calls ``Bio.KEGG.REST.kegg_info(database)``, for example ``"kegg"``, ``"pathway"`` or ``"hsa"``.
    """
    db = _text(database, "database", max_length=64)
    kegg = _bio("Bio.KEGG.REST")

    return _envelope(
        "kegg_rest_info",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_info",
        {"database": db},
        {"text": _read_text(kegg.kegg_info(db))},
    )


def kegg_rest_list(
    database: str | Sequence[str],
    *,
    organism: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List KEGG entries of a database (``/list``), as tab-separated rows.

    Calls ``Bio.KEGG.REST.kegg_list(database, org=organism)``. ``database`` is a database name
    (``"pathway"``, ``"hsa"``, ``"compound"`` and so on) or a list of up to 100
    entries. ``organism``
    is only used with ``"pathway"`` or ``"module"``.
    """
    kegg = _bio("Bio.KEGG.REST")

    if isinstance(database, str):
        target: Any = _text(database, "database", max_length=64)
    elif isinstance(database, (list, tuple)) and 0 < len(database) <= _MAX_KEGG_LIST:
        target = [_text(item, "database entry", max_length=128) for item in database]
    else:
        raise ConfigurationError(
            f"database must be a name or a list of 1 to {_MAX_KEGG_LIST} entries.",
            code="INVALID_QUERY",
        )
    org = None if organism is None else _text(organism, "organism", max_length=16)
    text = _read_text(kegg.kegg_list(target, org=org) if org else kegg.kegg_list(target))
    return _envelope(
        "kegg_rest_list",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_list",
        {"database": target, "organism": org},
        {"text": text, "rows": _rows(text)},
    )


def kegg_rest_find(
    database: str,
    query: str | Sequence[str],
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a KEGG database by keyword or by molecular property (``/find``), as rows.

    Calls ``Bio.KEGG.REST.kegg_find(database, query, option)``. ``option`` (``"formula"``,
    ``"exact_mass"`` or ``"mol_weight"``) is only valid for the compound and drug databases.
    """
    kegg = _bio("Bio.KEGG.REST")

    db = _text(database, "database", max_length=64)
    if isinstance(query, str):
        search: Any = _text(query, "query", max_length=1024)
    elif isinstance(query, (list, tuple)) and 0 < len(query) <= _MAX_KEGG_GET:
        search = [_text(item, "query", max_length=1024) for item in query]
    else:
        raise ConfigurationError(
            "query must be text or a list of 1 to 10 search terms.", code="INVALID_QUERY"
        )
    chosen = _option(option, _KEGG_FIND_OPTIONS, "option")
    args: tuple[Any, ...] = (db, search) if chosen is None else (db, search, chosen)
    text = _read_text(kegg.kegg_find(*args))
    return _envelope(
        "kegg_rest_find",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_find",
        {"database": db, "query": search, "option": chosen},
        {"text": text, "rows": _rows(text)},
    )


def kegg_rest_get(
    entries: str | Sequence[str],
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve KEGG entries (``/get``): flat-file text, amino-acid or nucleotide sequences.

    Calls ``Bio.KEGG.REST.kegg_get(dbentries, option)``. Up to 10 entries per call, for example
    ``["hsa:672", "hsa:675"]``. ``option`` is ``None`` (flat file), ``"aaseq"``, ``"ntseq"``,
    ``"mol"``, ``"kcf"``, ``"kgml"`` or ``"json"``. The KEGG ``image`` option returns binary PNG
    data and is not offered here.
    """
    kegg = _bio("Bio.KEGG.REST")

    if isinstance(entries, str):
        joined = _text(entries, "entries", max_length=1024)
        listed = [joined]
    elif isinstance(entries, (list, tuple)) and 0 < len(entries) <= _MAX_KEGG_GET:
        listed = [_text(item, "entries", max_length=128) for item in entries]
        joined = listed if len(listed) > 1 else listed[0]
    else:
        raise ConfigurationError(
            "entries must be an entry name or a list of 1 to 10 entries.", code="INVALID_QUERY"
        )
    chosen = _option(option, _KEGG_GET_OPTIONS, "option")
    args: tuple[Any, ...] = (joined,) if chosen is None else (joined, chosen)
    text = _read_text(kegg.kegg_get(*args))
    return _envelope(
        "kegg_rest_get",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_get",
        {"entries": listed, "option": chosen},
        {"text": text},
    )


def kegg_rest_conv(
    target_db: str,
    source: str,
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert identifiers between KEGG and outside databases (``/conv``), as tab-separated rows.

    Calls ``Bio.KEGG.REST.kegg_conv(target_db, source)``. ``source`` is a database name (for example
    ``"hsa"``) or a list of entries (for example ``"hsa:672"``). Pairs in the output are
    ``source<TAB>target``. ``option`` is ``None``, ``"turtle"`` or ``"n-triple"``.
    """
    kegg = _bio("Bio.KEGG.REST")

    target = _text(target_db, "target_db", max_length=64)
    origin = _text(source, "source", max_length=1024)
    chosen = _option(option, _KEGG_CONV_OPTIONS, "option")
    args: tuple[Any, ...] = (target, origin) if chosen is None else (target, origin, chosen)
    text = _read_text(kegg.kegg_conv(*args))
    return _envelope(
        "kegg_rest_conv",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_conv",
        {"target_db": target, "source": origin, "option": chosen},
        {"text": text, "rows": _rows(text)},
    )


def kegg_rest_link(
    target_db: str,
    source: str,
    *,
    option: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find KEGG entries linked to a source database or entry (``/link``), as tab-separated rows.

    Calls ``Bio.KEGG.REST.kegg_link(target_db, source)``. For example ``kegg_rest_link("pathway",
    "hsa:672")`` returns the pathways of BRCA1. Pairs are ``source<TAB>target``. ``option`` is
    ``None``, ``"turtle"`` or ``"n-triple"``.
    """
    kegg = _bio("Bio.KEGG.REST")

    target = _text(target_db, "target_db", max_length=64)
    origin = _text(source, "source", max_length=1024)
    chosen = _option(option, _KEGG_CONV_OPTIONS, "option")
    args: tuple[Any, ...] = (target, origin) if chosen is None else (target, origin, chosen)
    text = _read_text(kegg.kegg_link(*args))
    return _envelope(
        "kegg_rest_link",
        _KEGG_SERVICE,
        "Bio.KEGG.REST.kegg_link",
        {"target_db": target, "source": origin, "option": chosen},
        {"text": text, "rows": _rows(text)},
    )


def togows_entry(
    database: str,
    identifiers: str | Sequence[str],
    *,
    format: str | None = None,
    field: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch database entries through TogoWS (for example NCBI nuccore or DDBJ records).

    Calls ``Bio.TogoWS.entry(db, id, format=..., field=...)``. ``format`` (for example ``"fasta"``,
    ``"gff"``, ``"json"``) and ``field`` must be ones TogoWS lists for the database; Biopython
    raises ``ValueError`` otherwise.
    """
    togows = _bio("Bio.TogoWS")
    db = _text(database, "database", max_length=64)
    id_string = _ids(identifiers, "identifiers")
    options: dict[str, Any] = {}
    if format is not None:
        options["format"] = _text(format, "format", max_length=32)
    if field is not None:
        options["field"] = _text(field, "field", max_length=64)
    text = _read_text(togows.entry(db, id_string, **options))
    return _envelope(
        "togows_entry",
        _TOGOWS_SERVICE,
        "Bio.TogoWS.entry",
        {"database": db, "identifiers": id_string.split(","), **options},
        {"text": text},
    )


def togows_search(
    database: str,
    query: str,
    *,
    offset: int | None = None,
    limit: int | None = None,
    format: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a TogoWS search (``/search``) and return the hits, one per line in plain-text mode.

    Calls ``Bio.TogoWS.search(db, query, offset=..., limit=..., format=...)``. ``offset`` and
    ``limit`` must be given together (1-based offset). Without ``format`` the ``ids`` list holds
    one hit per line; with ``format`` the raw response text is returned and ``ids`` is ``None``.
    """
    togows = _bio("Bio.TogoWS")
    db = _text(database, "database", max_length=64)
    search = _text(query, "query", max_length=1024)
    options: dict[str, Any] = {}
    if offset is not None:
        options["offset"] = _bounded_int(offset, "offset", minimum=1, maximum=10_000_000)
    if limit is not None:
        options["limit"] = _bounded_int(limit, "limit", minimum=1, maximum=10_000)
    if format is not None:
        options["format"] = _text(format, "format", max_length=32)
    text = _read_text(togows.search(db, search, **options))
    ids = None if format is not None else [line for line in text.splitlines() if line.strip()]
    return _envelope(
        "togows_search",
        _TOGOWS_SERVICE,
        "Bio.TogoWS.search",
        {"database": db, "query": search, **options},
        {"text": text, "ids": ids},
    )


def togows_search_count(database: str, query: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return how many TogoWS hits a search has (``/search/<db>/<query>/count``).

    Calls ``Bio.TogoWS.search_count(db, query)``. The count is an integer.
    """
    togows = _bio("Bio.TogoWS")
    db = _text(database, "database", max_length=64)
    search = _text(query, "query", max_length=1024)
    count = togows.search_count(db, search)
    return _envelope(
        "togows_search_count",
        _TOGOWS_SERVICE,
        "Bio.TogoWS.search_count",
        {"database": db, "query": search},
        {"count": int(count)},
    )


def togows_search_iter(
    database: str,
    query: str,
    *,
    limit: int = 100,
    batch: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect up to ``limit`` TogoWS hit identifiers, fetching ``batch`` per request.

    Calls ``Bio.TogoWS.search_iter(db, query, limit=..., batch=...)`` and stops after ``limit``
    identifiers, so the number of requests is bounded. ``limit`` is required to be a positive
    integer (no unbounded download).
    """
    togows = _bio("Bio.TogoWS")
    db = _text(database, "database", max_length=64)
    search = _text(query, "query", max_length=1024)
    bounded = _bounded_int(limit, "limit", minimum=1, maximum=100_000)
    step = _bounded_int(batch, "batch", minimum=1, maximum=1_000)
    ids = list(islice(togows.search_iter(db, search, limit=bounded, batch=step), bounded))
    return _envelope(
        "togows_search_iter",
        _TOGOWS_SERVICE,
        "Bio.TogoWS.search_iter",
        {"database": db, "query": search, "limit": bounded, "batch": step},
        {"ids": [str(item) for item in ids], "count": len(ids)},
    )


def togows_convert(
    data: str,
    in_format: str,
    out_format: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a record between formats with TogoWS (for example GenBank to FASTA).

    Calls ``Bio.TogoWS.convert(data, in_format, out_format)``. Biopython accepts only the pairs
    TogoWS lists; others raise ``ValueError``. For local conversions, use Biopython's SeqIO.
    """
    togows = _bio("Bio.TogoWS")
    payload = _text(data, "data", max_length=50_000_000)
    source = _text(in_format, "in_format", max_length=32)
    target = _text(out_format, "out_format", max_length=32)
    text = _read_text(togows.convert(payload, source, target))
    return _envelope(
        "togows_convert",
        _TOGOWS_SERVICE,
        "Bio.TogoWS.convert",
        {"in_format": source, "out_format": target},
        {"text": text},
    )


def uniprot_search_records(
    query: str,
    *,
    fields: Sequence[str] | None = None,
    limit: int = 100,
    batch_size: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search UniProtKB (for example ``"gene:BRCA1 AND organism_id:9606 AND reviewed:true"``).

    Calls ``Bio.UniProt.search(query, fields=..., batch_size=...)`` and returns at most ``limit``
    records as JSON objects. ``fields`` are UniProt return-field names (for example
    ``["accession", "gene_names"]``); ``None`` returns the service's default record layout, which
    is large, so set ``limit`` and ``batch_size`` accordingly.
    """
    uniprot = _bio("Bio.UniProt")
    search = _text(query, "query", max_length=2048)
    bounded = _bounded_int(limit, "limit", minimum=1, maximum=10_000)
    step = _bounded_int(batch_size, "batch_size", minimum=1, maximum=500)
    names: list[str] | None = None
    if fields is not None:
        if not isinstance(fields, (list, tuple)) or not fields:
            raise ConfigurationError(
                "fields must be a non-empty list of names.", code="INVALID_QUERY"
            )
        names = [_text(item, "field", max_length=64) for item in fields]
    results = uniprot.search(search, fields=names, batch_size=step)
    records = [_plain(item) for item in islice(results, bounded)]
    return _envelope(
        "uniprot_search_records",
        _UNIPROT_SERVICE,
        "Bio.UniProt.search",
        {"query": search, "fields": names, "limit": bounded, "batch_size": step},
        {"records": records, "returned": len(records)},
    )


def expasy_sprot_raw_entry(identifier: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the raw UniProtKB/Swiss-Prot flat-file text for an accession or entry name.

    Calls ``Bio.ExPASy.get_sprot_raw(identifier)``. Biopython raises ``ValueError`` when the entry
    does not exist.
    """
    expasy = _bio("Bio.ExPASy")
    entry = _text(identifier, "identifier", max_length=128)
    return _envelope(
        "expasy_sprot_raw_entry",
        _EXPASY_SERVICE,
        "Bio.ExPASy.get_sprot_raw",
        {"identifier": entry},
        {"text": _read_text(expasy.get_sprot_raw(entry))},
    )


def _blast_result(record: Any) -> dict[str, Any]:
    """Plain summary of one parsed BLAST XML record (hits and their HSPs)."""
    hits: list[dict[str, Any]] = []
    for alignment in getattr(record, "alignments", []):
        hsps = []
        for hsp in alignment.hsps:
            row = {name: _plain(getattr(hsp, name, None)) for name in _HSP_FIELDS}
            hsps.append(row)
        hits.append(
            {
                "hit_id": _plain(getattr(alignment, "hit_id", None)),
                "hit_def": _plain(getattr(alignment, "hit_def", None)),
                "title": _plain(getattr(alignment, "title", None)),
                "accession": _plain(getattr(alignment, "accession", None)),
                "length": _plain(getattr(alignment, "length", None)),
                "hsps": hsps,
            }
        )
    return {
        "program": _plain(getattr(record, "application", None)),
        "database": _plain(getattr(record, "database", None)),
        "query_id": _plain(getattr(record, "query_id", None)),
        "query_letters": _plain(getattr(record, "query_letters", None)),
        "hit_count": len(hits),
        "hits": hits,
    }


def _run_with_timeout(task: Callable[[], Any], timeout: float | None) -> Any:
    """Run a blocking call in a daemon thread and give up after ``timeout`` seconds."""
    if timeout is None:
        return task()
    outcome: dict[str, Any] = {}

    def target() -> None:
        try:
            outcome["value"] = task()
        except BaseException as exc:  # re-raised in the calling thread below
            outcome["error"] = exc

    worker = threading.Thread(target=target, name="dnakit-ncbi-qblast", daemon=True)
    worker.start()
    worker.join(timeout)
    if worker.is_alive():
        raise BackendTimeoutError(
            f"NCBI BLAST did not return within {timeout:g} seconds.",
            hint=(
                "The remote job may still be running on NCBI; qblast does not expose its RID "
                "while it waits, so the job cannot be cancelled from DNAKit. Re-run with a "
                "larger timeout, or use a shorter query or database."
            ),
        )
    if "error" in outcome:
        raise outcome["error"]
    return outcome.get("value")


def ncbi_qblast_search(
    program: str,
    database: str,
    sequence: str,
    *,
    timeout: float | None = 1800.0,
    expect: float = 10.0,
    hitlist_size: int = 50,
    entrez_query: str = "(none)",
    word_size: int | None = None,
    matrix_name: str | None = None,
    megablast: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Submit one BLAST search to NCBI QBLAST and return the parsed hits.

    Calls ``Bio.Blast.NCBIWWW.qblast(program, database, sequence, ...)`` and parses the XML with
    ``Bio.Blast.NCBIXML``. ``program`` is for example ``"blastn"``, ``"blastp"`` or ``"blastx"``;
    ``database`` is an NCBI BLAST database (for example ``"nt"`` or ``"refseq_protein"``);
    ``sequence`` is FASTA text or a bare sequence. This is a long-running remote job: NCBI polls
    every 60 seconds or longer. ``timeout`` (seconds, default 1800) bounds how long the call waits;
    ``None`` waits without limit. On timeout the function raises ``BackendTimeoutError``; the NCBI
    job itself continues and is not cancelled.
    """
    qblast_program = _text(program, "program", max_length=32)
    qblast_database = _text(database, "database", max_length=128)
    qblast_sequence = _text(sequence, "sequence", max_length=1_000_000)
    if timeout is not None and (not isinstance(timeout, (int, float)) or timeout <= 0):
        raise ConfigurationError(
            "timeout must be a positive number of seconds or None.", code="INVALID_QUERY"
        )
    options: dict[str, Any] = {
        "expect": float(expect),
        "hitlist_size": _bounded_int(hitlist_size, "hitlist_size", minimum=1, maximum=5000),
        "entrez_query": _text(entrez_query, "entrez_query", max_length=1024),
    }
    if word_size is not None:
        options["word_size"] = _bounded_int(word_size, "word_size", minimum=2, maximum=512)
    if matrix_name is not None:
        options["matrix_name"] = _text(matrix_name, "matrix_name", max_length=32)
    if megablast is not None:
        options["megablast"] = bool(megablast)
    ncbiwww = _bio("Bio.Blast.NCBIWWW")
    ncbixml = _bio("Bio.Blast.NCBIXML")

    def job() -> dict[str, Any]:
        handle = ncbiwww.qblast(qblast_program, qblast_database, qblast_sequence, **options)
        try:
            record = ncbixml.read(handle)
        finally:
            with contextlib.suppress(Exception):
                handle.close()
        return _blast_result(record)

    result = _run_with_timeout(job, None if timeout is None else float(timeout))
    return _envelope(
        "ncbi_qblast_search",
        _BLAST_SERVICE,
        "Bio.Blast.NCBIWWW.qblast",
        {
            "program": qblast_program,
            "database": qblast_database,
            "timeout_seconds": timeout,
            **options,
        },
        result,
    )


__all__ = [
    "entrez_fetch",
    "entrez_info",
    "entrez_link",
    "entrez_post",
    "entrez_search",
    "entrez_spell",
    "entrez_summary",
    "expasy_sprot_raw_entry",
    "kegg_rest_conv",
    "kegg_rest_find",
    "kegg_rest_get",
    "kegg_rest_info",
    "kegg_rest_link",
    "kegg_rest_list",
    "ncbi_qblast_search",
    "togows_convert",
    "togows_entry",
    "togows_search",
    "togows_search_count",
    "togows_search_iter",
    "uniprot_search_records",
]

_register_public_verbose(globals())
