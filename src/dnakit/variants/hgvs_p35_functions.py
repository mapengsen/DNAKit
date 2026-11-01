"""hgvs wrappers, pack P35: parsing, validation, coordinate mapping, projection and data-provider queries.

hgvs is imported only when a function runs. A missing package raises ``BackendUnavailableError`` with the
``dnakit[bio-hgvs]`` hint.

Functions that need transcript alignments or sequences (mapping, extrinsic validation, reference filling,
projection, protein consequence, the UTA queries) take a data source:

- ``data_source="local"`` (default) reads a ``local_data`` snapshot that you pass: a dictionary with the keys
  ``sequences`` (accession to DNA/protein string), ``transcripts`` (transcript accession to ``hgnc``,
  ``cds_start_i``, ``cds_end_i``, ``pro_ac``), ``alignments`` (a list of ``{"tx_ac", "alt_ac",
  "alt_aln_method", "alt_strand", "exons"}`` whose exons are ``{"ord", "tx_start_i", "tx_end_i",
  "alt_start_i", "alt_end_i", "cigar"}`` in transcript order), ``genes`` (HGNC symbol to a record) and
  optionally ``assemblies`` (assembly name to ``{accession: chromosome name}``). It implements the hgvs
  data-provider interface, so every mapper behaves as with UTA, without any network access.
- ``data_source="uta"`` connects to a UTA PostgreSQL database (the public ``uta.biocommons.org`` instance
  or ``uta_db_url``). It needs the ``psycopg2`` package (extra ``bio-hgvs-uta``) and a network connection;
  sequences come from SeqRepo (``HGVS_SEQREPO_DIR``) or NCBI. Nothing connects until a function runs.

The NCBI-mirror queries of hgvs (HGNC to gene ID and similar) are not provided: they need an ``ncbi``
PostgreSQL database that hgvs does not publish and the project does not host.
"""

from __future__ import annotations

import functools
import importlib.metadata
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-hgvs"
_UTA_EXTRA = "bio-hgvs-uta"
_SOURCES = ("local", "uta")
_TYPES = ("c", "g", "m", "n", "p", "r")
_QUERIES = (
    "gene_info",
    "tx_exons",
    "tx_for_gene",
    "tx_identity_info",
    "tx_info",
    "tx_mapping_options",
    "similar_transcripts",
    "protein_accession",
    "assembly_map",
)
_PREVALIDATION = ("NONE", "INTRINSIC", "EXTRINSIC")
_MAX_ROWS = 500
_MAX_LOCAL_SEQUENCE = 5_000_000


def _text(value: object, name: str, *, maximum: int = 500) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {maximum} characters.")
    return value.strip()


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of: {', '.join(allowed)}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version("hgvs")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        installed = None
    return {"package": "hgvs", "version": installed, **notes}


def _parser() -> Any:
    return require("hgvs.parser", _EXTRA).Parser()


def _parse_variant(text: object, name: str = "variant") -> Any:
    value = _text(text, name)
    try:
        return _parser().parse_hgvs_variant(value)
    except Exception as exc:  # hgvs raises HGVSParseError and grammar errors
        raise ConfigurationError(f"{name} is not a valid HGVS variant: {value!r} ({exc})") from exc


def _hgvs_call(action: Any, message: str) -> Any:
    exceptions = require("hgvs.exceptions", _EXTRA)
    try:
        return action()
    except ConfigurationError:
        raise
    except exceptions.HGVSError as exc:
        raise ConfigurationError(f"{message}: {type(exc).__name__}: {exc}") from exc
    except (AssertionError, AttributeError, IndexError, KeyError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"{message}: {type(exc).__name__}: {exc}") from exc


def _variant_info(variant: Any, *, show_reference: bool = False) -> dict[str, Any]:
    posedit = getattr(variant, "posedit", None)
    return {
        "hgvs": variant.format(conf={"max_ref_length": None}) if show_reference else str(variant),
        "accession": str(variant.ac),
        "type": str(variant.type),
        "posedit": None if posedit is None else str(posedit),
    }


# ---------------------------------------------------------------------------------------------
# Local data provider
# ---------------------------------------------------------------------------------------------


def _check_local_data(data: object) -> dict[str, Any]:
    if not isinstance(data, Mapping):
        raise ConfigurationError("local_data must be a dictionary (see the module description).")
    unknown = set(data) - {"sequences", "transcripts", "alignments", "genes", "assemblies"}
    if unknown:
        raise ConfigurationError(f"local_data has unknown keys: {sorted(unknown)}.")
    result: dict[str, Any] = {"sequences": {}, "transcripts": {}, "alignments": [], "genes": {}, "assemblies": {}}
    sequences = data.get("sequences") or {}
    if not isinstance(sequences, Mapping):
        raise ConfigurationError("local_data['sequences'] must map accessions to sequences.")
    for key, value in sequences.items():
        if not isinstance(key, str) or not isinstance(value, str) or len(value) > _MAX_LOCAL_SEQUENCE:
            raise ConfigurationError("local_data['sequences'] must map accession strings to sequence strings.")
        result["sequences"][key] = value.upper()
    transcripts = data.get("transcripts") or {}
    if not isinstance(transcripts, Mapping):
        raise ConfigurationError("local_data['transcripts'] must map transcript accessions to records.")
    for key, record in transcripts.items():
        if not isinstance(record, Mapping) or set(record) - {"hgnc", "cds_start_i", "cds_end_i", "pro_ac"}:
            raise ConfigurationError(
                f"local_data['transcripts'][{key!r}] must hold hgnc, cds_start_i, cds_end_i and pro_ac only."
            )
        start, end = record.get("cds_start_i"), record.get("cds_end_i")
        if (start is None) != (end is None):
            raise ConfigurationError(f"transcript {key!r}: cds_start_i and cds_end_i must both be set or both null.")
        result["transcripts"][key] = dict(record)
    alignments = data.get("alignments") or []
    if isinstance(alignments, (str, bytes)) or not isinstance(alignments, Sequence):
        raise ConfigurationError("local_data['alignments'] must be a list.")
    required = {"ord", "tx_start_i", "tx_end_i", "alt_start_i", "alt_end_i", "cigar"}
    for index, alignment in enumerate(alignments):
        if not isinstance(alignment, Mapping) or set(alignment) - {"tx_ac", "alt_ac", "alt_aln_method", "alt_strand", "exons"}:
            raise ConfigurationError(f"local_data['alignments'][{index}] has unknown or missing structure.")
        if alignment.get("tx_ac") not in result["transcripts"]:
            raise ConfigurationError(f"alignments[{index}].tx_ac is not in local_data['transcripts'].")
        if alignment.get("alt_strand") not in (1, -1):
            raise ConfigurationError(f"alignments[{index}].alt_strand must be 1 or -1.")
        exons = alignment.get("exons")
        if not isinstance(exons, Sequence) or not exons:
            raise ConfigurationError(f"alignments[{index}].exons must be a non-empty list.")
        for exon in exons:
            if not isinstance(exon, Mapping) or set(exon) != required:
                raise ConfigurationError(f"alignments[{index}] exons need exactly the keys {sorted(required)}.")
        result["alignments"].append(
            {
                "tx_ac": alignment["tx_ac"],
                "alt_ac": _text(alignment.get("alt_ac"), "alt_ac"),
                "alt_aln_method": _text(alignment.get("alt_aln_method"), "alt_aln_method"),
                "alt_strand": alignment["alt_strand"],
                "exons": [dict(exon) for exon in exons],
            }
        )
    genes = data.get("genes") or {}
    if not isinstance(genes, Mapping):
        raise ConfigurationError("local_data['genes'] must map HGNC symbols to records.")
    result["genes"] = {str(k): dict(v) for k, v in genes.items()}
    assemblies = data.get("assemblies") or {}
    if not isinstance(assemblies, Mapping):
        raise ConfigurationError("local_data['assemblies'] must map assembly names to accession maps.")
    result["assemblies"] = {str(k): dict(v) for k, v in assemblies.items()}
    return result


@functools.lru_cache(maxsize=1)
def _provider_class() -> Any:
    interface = require("hgvs.dataproviders.interface", _EXTRA)
    exceptions = require("hgvs.exceptions", _EXTRA)
    digests = require("bioutils.digests", _EXTRA)
    assemblies = require("bioutils.assemblies", _EXTRA)

    class _SeqFetcher:
        def __init__(self, provider: Any) -> None:
            self._provider = provider

        def fetch_seq(self, ac: str, start_i: int | None = None, end_i: int | None = None) -> str:
            return self._provider.get_seq(ac, start_i, end_i)

    class LocalDataProvider(interface.Interface):  # type: ignore[misc]
        """hgvs data provider that answers from an in-memory snapshot."""

        required_version = "1.1"

        def __init__(self, data: dict[str, Any]) -> None:
            self._data = data
            self.seqfetcher = _SeqFetcher(self)
            super().__init__()

        def __str__(self) -> str:
            return "LocalDataProvider <data_version:local; schema_version:1.1>"

        def data_version(self) -> str:
            return "local"

        def schema_version(self) -> str:
            return "1.1"

        def get_seq(self, ac: str, start_i: int | None = None, end_i: int | None = None) -> str:
            try:
                sequence = self._data["sequences"][ac]
            except KeyError:
                raise exceptions.HGVSDataNotAvailableError(f"No sequence for {ac} in the local data") from None
            return sequence[start_i:end_i]

        def get_acs_for_protein_seq(self, seq: str) -> list[str]:
            md5 = digests.seq_md5(seq)
            return [ac for ac, value in self._data["sequences"].items() if value == seq.upper()] + ["MD5_" + md5]

        def get_gene_info(self, gene: str) -> dict[str, Any] | None:
            record = self._data["genes"].get(gene)
            return None if record is None else {"hgnc": gene, **record}

        def _alignments(self, tx_ac: str | None = None) -> list[dict[str, Any]]:
            return [a for a in self._data["alignments"] if tx_ac is None or a["tx_ac"] == tx_ac]

        def _exon_rows(self, alignment: dict[str, Any]) -> list[dict[str, Any]]:
            rows = []
            for exon in alignment["exons"]:
                rows.append(
                    {
                        "tx_ac": alignment["tx_ac"],
                        "alt_ac": alignment["alt_ac"],
                        "alt_strand": alignment["alt_strand"],
                        "alt_aln_method": alignment["alt_aln_method"],
                        **exon,
                    }
                )
            return sorted(rows, key=lambda row: row["alt_start_i"])

        def get_tx_exons(self, tx_ac: str, alt_ac: str, alt_aln_method: str) -> list[dict[str, Any]]:
            for alignment in self._alignments(tx_ac):
                if alignment["alt_ac"] == alt_ac and alignment["alt_aln_method"] == alt_aln_method:
                    rows = self._exon_rows(alignment)
                    first = 0 if rows[0]["alt_strand"] == 1 else -1
                    if rows[first]["tx_start_i"] != 0:
                        raise exceptions.HGVSDataNotAvailableError(
                            f"Alignment is incomplete; cannot use transcript for mapping (tx_ac={tx_ac},"
                            f"alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
                        )
                    return rows
            raise exceptions.HGVSDataNotAvailableError(
                f"No tx_exons for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
            )

        def get_tx_for_gene(self, gene: str) -> list[dict[str, Any]]:
            rows = []
            for alignment in self._alignments():
                info = self._data["transcripts"][alignment["tx_ac"]]
                if info.get("hgnc") == gene and alignment["alt_aln_method"] != "transcript":
                    rows.append(self._info_row(alignment))
            return rows

        def _info_row(self, alignment: dict[str, Any]) -> dict[str, Any]:
            info = self._data["transcripts"][alignment["tx_ac"]]
            return {
                "hgnc": info.get("hgnc"),
                "cds_start_i": info.get("cds_start_i"),
                "cds_end_i": info.get("cds_end_i"),
                "tx_ac": alignment["tx_ac"],
                "alt_ac": alignment["alt_ac"],
                "alt_aln_method": alignment["alt_aln_method"],
            }

        def get_tx_for_region(self, alt_ac: str, alt_aln_method: str, start_i: int, end_i: int) -> list[dict[str, Any]]:
            return self.get_alignments_for_region(alt_ac, start_i, end_i, alt_aln_method)

        def get_alignments_for_region(
            self, alt_ac: str, start_i: int, end_i: int, alt_aln_method: str | None = None
        ) -> list[dict[str, Any]]:
            # UTA's query: transcripts whose aligned span starts before start_i and ends at or after end_i.
            rows = []
            for alignment in self._alignments():
                if alignment["alt_ac"] != alt_ac:
                    continue
                if alt_aln_method is not None and alignment["alt_aln_method"] != alt_aln_method:
                    continue
                low = min(e["alt_start_i"] for e in alignment["exons"])
                high = max(e["alt_end_i"] for e in alignment["exons"])
                if low < start_i and end_i <= high:
                    rows.append(
                        {
                            "tx_ac": alignment["tx_ac"],
                            "alt_ac": alt_ac,
                            "alt_strand": alignment["alt_strand"],
                            "alt_aln_method": alignment["alt_aln_method"],
                            "start_i": low,
                            "end_i": high,
                        }
                    )
            return rows

        def get_tx_identity_info(self, tx_ac: str) -> dict[str, Any]:
            info = self._data["transcripts"].get(tx_ac)
            lengths = None
            for alignment in self._alignments(tx_ac):
                ordered = sorted(alignment["exons"], key=lambda e: e["ord"])
                lengths = [e["tx_end_i"] - e["tx_start_i"] for e in ordered]
                break
            if info is None or lengths is None:
                if info is not None and tx_ac in self._data["sequences"]:
                    lengths = [len(self._data["sequences"][tx_ac])]
                else:
                    raise exceptions.HGVSDataNotAvailableError(f"No transcript definition for (tx_ac={tx_ac})")
            return {
                "tx_ac": tx_ac,
                "alt_ac": tx_ac,
                "alt_aln_method": "transcript",
                "cds_start_i": info.get("cds_start_i"),
                "cds_end_i": info.get("cds_end_i"),
                "lengths": lengths,
                "hgnc": info.get("hgnc"),
            }

        def get_tx_info(self, tx_ac: str, alt_ac: str, alt_aln_method: str) -> dict[str, Any]:
            rows = [
                self._info_row(a)
                for a in self._alignments(tx_ac)
                if a["alt_ac"] == alt_ac and a["alt_aln_method"] == alt_aln_method
            ]
            if not rows:
                raise exceptions.HGVSDataNotAvailableError(
                    f"No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
                )
            return rows[0]

        def get_tx_mapping_options(self, tx_ac: str) -> list[dict[str, Any]]:
            return [
                {"tx_ac": a["tx_ac"], "alt_ac": a["alt_ac"], "alt_aln_method": a["alt_aln_method"]}
                for a in self._alignments(tx_ac)
                if a["alt_aln_method"] != "transcript"
            ]

        def _cds_sequence(self, tx_ac: str) -> str | None:
            info = self._data["transcripts"].get(tx_ac)
            sequence = self._data["sequences"].get(tx_ac)
            if info is None or sequence is None or info.get("cds_start_i") is None:
                return None
            return sequence[info["cds_start_i"] : info["cds_end_i"]]

        @staticmethod
        def _gap_free(alignment: dict[str, Any]) -> bool:
            return all(set(e["cigar"]) <= set("0123456789=XM") for e in alignment["exons"])

        def _structures(self, tx_ac: str) -> dict[tuple[str, str], tuple[Any, Any]]:
            info = self._data["transcripts"].get(tx_ac, {})
            result = {}
            for alignment in self._alignments(tx_ac):
                if not self._gap_free(alignment):
                    continue
                full = tuple(sorted((e["alt_start_i"], e["alt_end_i"]) for e in alignment["exons"]))
                clipped: tuple[Any, ...] | None = None
                if info.get("cds_start_i") is not None:
                    pieces = []
                    for e in alignment["exons"]:
                        lo, hi = max(e["tx_start_i"], info["cds_start_i"]), min(e["tx_end_i"], info["cds_end_i"])
                        if lo >= hi:
                            continue
                        if alignment["alt_strand"] == 1:
                            pieces.append((e["alt_start_i"] + lo - e["tx_start_i"], e["alt_start_i"] + hi - e["tx_start_i"]))
                        else:
                            pieces.append((e["alt_end_i"] - (hi - e["tx_start_i"]), e["alt_end_i"] - (lo - e["tx_start_i"])))
                    clipped = tuple(sorted(pieces))
                result[(alignment["alt_ac"], alignment["alt_aln_method"])] = (full, clipped)
            return result

        def get_similar_transcripts(self, tx_ac: str) -> list[dict[str, Any]]:
            mine = self._structures(tx_ac)
            my_cds = self._cds_sequence(tx_ac)
            rows = []
            for other in sorted(self._data["transcripts"]):
                if other == tx_ac:
                    continue
                theirs = self._structures(other)
                shared = set(mine) & set(theirs)
                cds = self._cds_sequence(other)
                row = {
                    "tx_ac1": tx_ac,
                    "tx_ac2": other,
                    "cds_eq": my_cds is not None and my_cds == cds,
                    "es_fp_eq": any(mine[k][0] == theirs[k][0] for k in shared),
                    "cds_es_fp_eq": any(mine[k][1] is not None and mine[k][1] == theirs[k][1] for k in shared),
                }
                if row["cds_eq"] or row["es_fp_eq"] or row["cds_es_fp_eq"]:
                    rows.append(row)
            return rows

        def get_pro_ac_for_tx_ac(self, tx_ac: str) -> str | None:
            info = self._data["transcripts"].get(tx_ac)
            return None if info is None else info.get("pro_ac")

        def get_assembly_map(self, assembly_name: str) -> dict[str, str]:
            if assembly_name in self._data["assemblies"]:
                return dict(self._data["assemblies"][assembly_name])
            return assemblies.make_ac_name_map(assembly_name)

    return LocalDataProvider


def _provider(data_source: str, local_data: object, uta_db_url: str | None) -> Any:
    kind = _choice(data_source, "data_source", _SOURCES)
    if kind == "local":
        if local_data is None:
            raise ConfigurationError(
                "data_source='local' needs local_data (a snapshot of sequences, transcripts and alignments)."
            )
        return _provider_class()(_check_local_data(local_data))
    require("psycopg2", _UTA_EXTRA)
    uta = require("hgvs.dataproviders.uta", _EXTRA)
    url = None if uta_db_url is None else _text(uta_db_url, "uta_db_url")
    return _hgvs_call(lambda: uta.connect(url), "hgvs could not connect to the UTA database")


def _rows(rows: Any) -> list[dict[str, Any]]:
    return [{str(k): _plain(v) for k, v in dict(row).items()} for row in list(rows)[:_MAX_ROWS]]


def _plain(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


# ---------------------------------------------------------------------------------------------
# Parsing and validation
# ---------------------------------------------------------------------------------------------


def parse_hgvs_posedit(
    text: str,
    *,
    variant_type: str = "c",
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse the position-plus-edit part of an HGVS string, for example ``22+1A>T`` (``c``) or ``Ala22Trp`` (``p``).

    ``variant_type`` selects the grammar rule: ``c`` (coding DNA, ``Parser.parse_c_posedit``), ``g``
    (genomic), ``m`` (mitochondrial), ``n`` (non-coding), ``p`` (protein) or ``r`` (RNA). The result gives
    the normalized text, the position (start, end) and the edit (type, reference, alternate).
    """

    kind = _choice(variant_type, "variant_type", _TYPES)
    value = _text(text, "text")
    parser = _parser()
    posedit = _hgvs_call(lambda: getattr(parser, f"parse_{kind}_posedit")(value), f"{value!r} is not a valid {kind}. posedit")
    edit = getattr(posedit, "edit", None)
    pos = getattr(posedit, "pos", None)
    return {
        "status": "completed",
        "operation": "parse_hgvs_posedit",
        "variant_type": kind,
        "posedit": str(posedit),
        "position": None if pos is None else {"start": str(pos.start), "end": str(pos.end)},
        "edit": None
        if edit is None
        else {
            "type": str(getattr(edit, "type", "")),
            "ref": _plain(getattr(edit, "ref", None)),
            "alt": _plain(getattr(edit, "alt", None)),
        },
        "uncertain": bool(getattr(posedit, "uncertain", False)),
        "provenance": _provenance(method=f"hgvs.parser.Parser.parse_{kind}_posedit"),
    }


def parse_hgvs_position(
    text: str,
    *,
    variant_type: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse an HGVS position or interval string, for example ``NM_000001.1:c.26+2_57-3`` or ``26+2_57-3``.

    A string with an accession and type prefix (``NM_000001.1:c.26+2_57-3``) is parsed by
    ``parse_hgvs_position``, or by ``parse_c_hgvs_position`` and its siblings when ``variant_type`` (``c``,
    ``g``, ``m``, ``n``, ``p`` or ``r``) is given. A bare interval (``26+2_57-3``) needs ``variant_type`` and
    uses ``parse_c_interval`` and its siblings. The result gives the start and end positions, with base,
    offset and datum for coding positions.
    """

    parser = _parser()
    value = _text(text, "text")
    if variant_type is None:
        method, label = parser.parse_hgvs_position, "parse_hgvs_position"
    else:
        kind = _choice(variant_type, "variant_type", _TYPES)
        label = f"parse_{kind}_hgvs_position" if ":" in value else f"parse_{kind}_interval"
        method = getattr(parser, label)
    parsed = _hgvs_call(lambda: method(value), f"{value!r} is not a valid HGVS position")
    interval = getattr(parsed, "pos", parsed)

    def describe(position: Any) -> dict[str, Any]:
        info: dict[str, Any] = {"text": str(position), "base": _plain(getattr(position, "base", None))}
        if hasattr(position, "offset"):
            info["offset"] = int(position.offset)
        if hasattr(position, "datum"):
            info["datum"] = str(position.datum.name if hasattr(position.datum, "name") else position.datum)
        return info

    return {
        "status": "completed",
        "operation": "parse_hgvs_position",
        "position": str(parsed),
        "accession": _plain(getattr(parsed, "ac", None)),
        "type": _plain(getattr(parsed, "type", variant_type)),
        "start": describe(interval.start),
        "end": describe(interval.end),
        "provenance": _provenance(method=f"hgvs.parser.Parser.{label}"),
    }


def validate_hgvs_position(
    text: str,
    *,
    variant_type: str = "c",
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate an HGVS position or interval intrinsically (no reference data needed).

    ``text`` is an interval such as ``26+2_57-3`` of the given ``variant_type``. The check is hgvs'
    ``Interval.validate``: base zero, bases below one for sequence-start positions and start after end are
    errors. Returns ``valid`` with the validation level and message.
    """

    enums = require("hgvs.enums", _EXTRA)
    kind = _choice(variant_type, "variant_type", _TYPES)
    value = _text(text, "text")
    parser = _parser()
    interval = _hgvs_call(lambda: getattr(parser, f"parse_{kind}_interval")(value), f"{value!r} is not a valid {kind}. interval")
    level, message = interval.validate()
    return {
        "status": "completed",
        "operation": "validate_hgvs_position",
        "position": str(interval),
        "valid": level == enums.ValidationLevel.VALID,
        "level": str(level.name),
        "message": message,
        "provenance": _provenance(method="hgvs.location.Interval.validate"),
    }


def validate_hgvs_variant(
    variant: str,
    *,
    level: str = "intrinsic",
    strict: bool = True,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate an HGVS variant: intrinsic, extrinsic or both.

    ``level="intrinsic"`` checks the string itself (accession and type agree, insertion and deletion
    lengths, position order; ``IntrinsicValidator``) and needs no data. ``"extrinsic"`` checks the variant
    against the reference sequence and transcript bounds (``ExtrinsicValidator``; needs ``data_source``).
    ``"full"`` does both (``Validator``, as ``hgvs.easy.validate``). ``strict`` treats warnings as failures,
    as hgvs does by default. Returns ``valid`` and, if invalid, the hgvs message.
    """

    validator = require("hgvs.validator", _EXTRA)
    exceptions = require("hgvs.exceptions", _EXTRA)
    kind = _choice(level, "level", ("intrinsic", "extrinsic", "full"))
    parsed = _parse_variant(variant)
    strictness = _flag(strict, "strict")
    if kind == "intrinsic":
        engine = validator.IntrinsicValidator(strict=strictness)
    else:
        provider = _provider(data_source, local_data, uta_db_url)
        engine = (
            validator.ExtrinsicValidator(provider, strict=strictness)
            if kind == "extrinsic"
            else validator.Validator(provider, strict=strictness)
        )
    try:
        engine.validate(parsed)
        valid, message = True, None
    except (exceptions.HGVSInvalidVariantError, exceptions.HGVSInvalidIntervalError) as exc:
        valid, message = False, str(exc)
    except exceptions.HGVSError as exc:
        raise ConfigurationError(f"hgvs could not validate the variant: {type(exc).__name__}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "validate_hgvs_variant",
        "variant": str(parsed),
        "level": kind,
        "valid": valid,
        "message": message,
        "provenance": _provenance(method=f"hgvs.validator ({kind})"),
    }


def fill_hgvs_reference(
    variant: str,
    *,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    genomic_accession: str | None = None,
    alt_aln_method: str = "splign",
    verbose: bool = True,
) -> dict[str, Any]:
    """Fill the reference bases of a deletion, duplication, delins or repeat variant from the reference sequence.

    ``NM_x:c.5_7del`` becomes ``NM_x:c.5_7delACG`` (``SequenceVariant.fill_ref``; the text is written with the
    reference bases shown, and ``reference`` holds them). For transcript variants
    the reference is read through the transcript's genomic alignment (``genomic_accession`` selects it);
    variants that carry their reference already are returned unchanged.
    """

    parsed = _parse_variant(variant)
    provider = _provider(data_source, local_data, uta_db_url)
    alt_ac = None if genomic_accession is None else _text(genomic_accession, "genomic_accession")
    filled = _hgvs_call(
        lambda: parsed.fill_ref(provider, alt_ac=alt_ac, alt_aln_method=_text(alt_aln_method, "alt_aln_method")),
        "hgvs could not fill the reference",
    )
    return {
        "status": "completed",
        "operation": "fill_hgvs_reference",
        "input": variant.strip(),
        "output": _variant_info(filled, show_reference=True),
        "reference": _plain(getattr(filled.posedit.edit, "ref", None)),
        "provenance": _provenance(method="hgvs.sequencevariant.SequenceVariant.fill_ref"),
    }


def hgvs_length_change(
    variant: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Get the net length change of an HGVS variant (``PosEdit.length_change``).

    Negative values are net deletions, zero is a balanced change such as an SNV, positive values are net
    insertions. Deletion lengths can be inferred from exact interval coordinates even
    when deleted bases are omitted. Undefined changes raise ``ConfigurationError``.
    """

    parsed = _parse_variant(variant)
    if parsed.posedit is None:
        raise ConfigurationError("variant has no position and edit.")
    change = _hgvs_call(lambda: parsed.posedit.length_change(on_error_raise=True), "hgvs could not compute the length change")
    return {
        "status": "completed",
        "operation": "hgvs_length_change",
        "variant": str(parsed),
        "edit_type": str(parsed.posedit.edit.type),
        "length_change": int(change),
        "provenance": _provenance(method="hgvs.posedit.PosEdit.length_change"),
    }


# ---------------------------------------------------------------------------------------------
# Mapping and projection
# ---------------------------------------------------------------------------------------------


def _mapper_options(prevalidation_level: str | None, alt_aln_method: str, normalize: bool, replace_reference: bool) -> dict[str, Any]:
    level = None if prevalidation_level is None else _choice(prevalidation_level.upper(), "prevalidation_level", _PREVALIDATION)
    return {
        "prevalidation_level": None if level == "NONE" else level,
        "alt_aln_method": _text(alt_aln_method, "alt_aln_method"),
        "normalize": _flag(normalize, "normalize"),
        "replace_reference": _flag(replace_reference, "replace_reference"),
    }


def map_hgvs_variant(
    variant: str,
    *,
    to_type: str,
    transcript_accession: str | None = None,
    genomic_accession: str | None = None,
    assembly_name: str | None = None,
    alt_aln_method: str = "splign",
    prevalidation_level: str | None = "EXTRINSIC",
    normalize: bool = True,
    replace_reference: bool = True,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map an HGVS variant between genomic (g.), coding (c.) and non-coding (n.) coordinates.

    One function for the hgvs mappers ``g_to_c``, ``g_to_n``, ``g_to_t``, ``n_to_g``, ``c_to_g``,
    ``t_to_g``, ``c_to_n`` and ``n_to_c``: the input type is read from the variant and ``to_type`` is ``g``,
    ``c``, ``n`` or ``t`` (``t`` picks c. or n. by transcript type). To a transcript type from g. give
    ``transcript_accession``; to g. from a transcript give ``genomic_accession`` or ``assembly_name`` (then
    hgvs' ``AssemblyMapper`` finds the chromosome and normalizes, as ``hgvs.easy`` does); c. <-> n. needs
    neither. ``prevalidation_level`` (``EXTRINSIC`` by default, as in hgvs), ``normalize`` and
    ``replace_reference`` follow hgvs' defaults. ``data_source`` is ``local`` (``local_data``) or ``uta``.
    """

    mapping = require("hgvs.variantmapper", _EXTRA)
    assembly = require("hgvs.assemblymapper", _EXTRA)
    target = _choice(to_type, "to_type", ("g", "c", "n", "t"))
    parsed = _parse_variant(variant)
    source = parsed.type
    allowed = {"g": ("c", "n", "t"), "m": ("c", "n", "t"), "c": ("g", "n", "t"), "n": ("g", "c", "t")}
    if source not in allowed or target not in allowed[source] and not (source in "cn" and target == "g"):
        raise ConfigurationError(f"cannot map a {source}. variant to {target}.; use one of {allowed.get(source, ())}.")
    options = _mapper_options(prevalidation_level, alt_aln_method, normalize, replace_reference)
    provider = _provider(data_source, local_data, uta_db_url)
    tx_ac = None if transcript_accession is None else _text(transcript_accession, "transcript_accession")
    alt_ac = None if genomic_accession is None else _text(genomic_accession, "genomic_accession")
    method = f"{source}_to_{target}" if source != "m" else f"g_to_{target}"
    if target == "t":
        method = "g_to_t" if source in "gm" else "t_to_g"
    if source in "gm" and tx_ac is None:
        raise ConfigurationError("mapping from g. needs transcript_accession.")
    if target in ("g", "t") and source in "cn" and alt_ac is None and assembly_name is None:
        raise ConfigurationError("mapping to g. needs genomic_accession or assembly_name.")
    if assembly_name is not None:
        mapper = _hgvs_call(
            lambda: assembly.AssemblyMapper(
                provider,
                assembly_name=_text(assembly_name, "assembly_name"),
                alt_aln_method=options["alt_aln_method"],
                normalize=options["normalize"],
                prevalidation_level=options["prevalidation_level"],
                replace_reference=options["replace_reference"],
            ),
            "hgvs could not create the AssemblyMapper",
        )
        args: tuple[Any, ...] = (parsed, tx_ac) if source in "gm" else (parsed,)
        used = "AssemblyMapper"
    else:
        mapper = mapping.VariantMapper(
            provider,
            replace_reference=options["replace_reference"],
            prevalidation_level=options["prevalidation_level"],
        )
        if source in "gm":
            args = (parsed, tx_ac)
        elif target in ("g", "t"):
            args = (parsed, alt_ac)
        else:
            args = (parsed,)
        used = "VariantMapper"
    result = _hgvs_call(lambda: getattr(mapper, method)(*args, alt_aln_method=options["alt_aln_method"]) if used == "VariantMapper" else getattr(mapper, method)(*args), f"hgvs could not map the variant ({method})")
    return {
        "status": "completed",
        "operation": "map_hgvs_variant",
        "input": variant.strip(),
        "mapper_method": method,
        "mapper": used,
        "output": _variant_info(result),
        "provenance": _provenance(method=f"hgvs.{'assemblymapper' if used == 'AssemblyMapper' else 'variantmapper'}.{used}.{method}"),
    }


def map_hgvs_interval(
    interval: str,
    *,
    from_type: str,
    to_type: str,
    transcript_accession: str,
    genomic_accession: str | None = None,
    alt_aln_method: str = "splign",
    strict_bounds: bool | None = None,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map a position interval between genomic, coding and non-coding coordinates with an ``AlignmentMapper``.

    ``interval`` is a bare interval such as ``26+2_57-3`` of ``from_type``; ``from_type`` and ``to_type``
    are two of ``g``, ``c`` and ``n``. The mapper is built for ``transcript_accession`` and its
    ``genomic_accession`` (not needed for c. <-> n.). ``strict_bounds`` overrides hgvs' default of requiring
    transcript positions inside the transcript. This is the interval form of ``map_hgvs_variant``
    (``TranscriptMapper``/``AlignmentMapper`` ``g_to_c`` and siblings).
    """

    alignment = require("hgvs.alignmentmapper", _EXTRA)
    source = _choice(from_type, "from_type", ("g", "c", "n"))
    target = _choice(to_type, "to_type", ("g", "c", "n"))
    if source == target:
        raise ConfigurationError("from_type and to_type must differ.")
    value = _text(interval, "interval")
    parsed = _hgvs_call(lambda: getattr(_parser(), f"parse_{source}_interval")(value), f"{value!r} is not a valid {source}. interval")
    provider = _provider(data_source, local_data, uta_db_url)
    tx_ac = _text(transcript_accession, "transcript_accession")
    if {source, target} == {"c", "n"}:
        alt_ac, method = tx_ac, "transcript"
    else:
        if genomic_accession is None:
            raise ConfigurationError("genomic_accession is required when mapping to or from g.")
        alt_ac, method = _text(genomic_accession, "genomic_accession"), _text(alt_aln_method, "alt_aln_method")
    mapper = _hgvs_call(lambda: alignment.AlignmentMapper(provider, tx_ac, alt_ac, method), "hgvs could not create the AlignmentMapper")
    bounds = None if strict_bounds is None else _flag(strict_bounds, "strict_bounds")
    out = _hgvs_call(lambda: getattr(mapper, f"{source}_to_{target}")(parsed, strict_bounds=bounds), "hgvs could not map the interval")
    return {
        "status": "completed",
        "operation": "map_hgvs_interval",
        "input": value,
        "mapper_method": f"{source}_to_{target}",
        "output": str(out),
        "start": str(out.start),
        "end": str(out.end),
        "provenance": _provenance(method=f"hgvs.alignmentmapper.AlignmentMapper.{source}_to_{target}"),
    }


def predict_hgvs_protein_consequence(
    variant: str,
    *,
    protein_accession: str | None = None,
    genomic_accession: str | None = None,
    assembly_name: str | None = None,
    alt_aln_method: str = "splign",
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the protein consequence (p.) of a coding (c.) variant.

    This is ``VariantMapper.c_to_p`` (``AssemblyMapper.c_to_p`` when ``assembly_name`` is given and
    ``AssemblyMapper.t_to_p``-style handling: a non-coding n. variant returns ``non-coding``). The protein
    accession defaults to the transcript's own; hgvs uses the standard translation table, adds
    ``p.(...)`` uncertainty brackets for inferred changes (``inferred_p_is_uncertain``) and reports
    ``p.?``, ``p.=`` or frameshift/extension notation as appropriate.
    """

    mapping = require("hgvs.variantmapper", _EXTRA)
    assembly = require("hgvs.assemblymapper", _EXTRA)
    parsed = _parse_variant(variant)
    provider = _provider(data_source, local_data, uta_db_url)
    method = _text(alt_aln_method, "alt_aln_method")
    if parsed.type == "n":
        return {
            "status": "completed",
            "operation": "predict_hgvs_protein_consequence",
            "input": variant.strip(),
            "output": "non-coding",
            "provenance": _provenance(method="hgvs.assemblymapper.AssemblyMapper.t_to_p"),
        }
    if parsed.type != "c":
        raise ConfigurationError("variant must be a coding (c.) variant.")
    if assembly_name is not None:
        mapper = assembly.AssemblyMapper(provider, assembly_name=_text(assembly_name, "assembly_name"), alt_aln_method=method)
        result = _hgvs_call(lambda: mapper.c_to_p(parsed), "hgvs could not predict the protein consequence")
    else:
        vm = mapping.VariantMapper(provider)
        result = _hgvs_call(
            lambda: vm.c_to_p(
                parsed,
                pro_ac=None if protein_accession is None else _text(protein_accession, "protein_accession"),
                alt_ac=None if genomic_accession is None else _text(genomic_accession, "genomic_accession"),
                alt_aln_method=method,
            ),
            "hgvs could not predict the protein consequence",
        )
    return {
        "status": "completed",
        "operation": "predict_hgvs_protein_consequence",
        "input": variant.strip(),
        "output": _variant_info(result),
        "provenance": _provenance(method="hgvs.variantmapper.VariantMapper.c_to_p"),
    }


def project_hgvs_variant(
    variant: str,
    *,
    source_accession: str,
    destination_accession: str,
    genomic_accession: str,
    direction: str = "forward",
    source_alt_aln_method: str = "splign",
    destination_alt_aln_method: str = "splign",
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project a coding variant or interval between two transcripts through a common genomic sequence.

    ``direction="forward"`` projects from ``source_accession`` to ``destination_accession``
    (``Projector.project_variant_forward`` / ``project_interval_forward``); ``"backward"`` projects a variant
    on the destination back to the source. ``variant`` is a full c. variant on the transcript of the
    starting side, or a bare c. interval such as ``26+2_57-3``. The transcripts must be on the same strand
    of ``genomic_accession``.
    """

    projector_module = require("hgvs.projector", _EXTRA)
    way = _choice(direction, "direction", ("forward", "backward"))
    provider = _provider(data_source, local_data, uta_db_url)
    projector = _hgvs_call(
        lambda: projector_module.Projector(
            provider,
            _text(genomic_accession, "genomic_accession"),
            _text(source_accession, "source_accession"),
            _text(destination_accession, "destination_accession"),
            src_alt_aln_method=_text(source_alt_aln_method, "source_alt_aln_method"),
            dst_alt_aln_method=_text(destination_alt_aln_method, "destination_alt_aln_method"),
        ),
        "hgvs could not create the Projector",
    )
    value = _text(variant, "variant")
    if ":" in value:
        parsed = _parse_variant(value)
        result = _hgvs_call(
            lambda: getattr(projector, f"project_variant_{way}")(parsed), "hgvs could not project the variant"
        )
        output: Any = _variant_info(result)
        kind = "variant"
    else:
        interval = _hgvs_call(lambda: _parser().parse_c_interval(value), f"{value!r} is not a valid c. interval")
        result = _hgvs_call(
            lambda: getattr(projector, f"project_interval_{way}")(interval), "hgvs could not project the interval"
        )
        output = str(result)
        kind = "interval"
    return {
        "status": "completed",
        "operation": "project_hgvs_variant",
        "input": value,
        "kind": kind,
        "direction": way,
        "output": output,
        "provenance": _provenance(method=f"hgvs.projector.Projector.project_{kind}_{way}"),
    }


def list_hgvs_relevant_transcripts(
    variant: str,
    *,
    assembly_name: str = "GRCh38",
    alt_aln_method: str = "splign",
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the transcripts whose genomic alignment covers a genomic (g.) variant.

    This is ``AssemblyMapper.relevant_transcripts`` (UTA's ``get_tx_for_region`` query): the transcripts
    whose aligned span starts before the variant start and reaches its end, on the variant's chromosome.
    """

    assembly = require("hgvs.assemblymapper", _EXTRA)
    parsed = _parse_variant(variant)
    if parsed.type not in "gm":
        raise ConfigurationError("variant must be a genomic (g.) variant.")
    provider = _provider(data_source, local_data, uta_db_url)
    mapper = _hgvs_call(
        lambda: assembly.AssemblyMapper(
            provider, assembly_name=_text(assembly_name, "assembly_name"), alt_aln_method=_text(alt_aln_method, "alt_aln_method")
        ),
        "hgvs could not create the AssemblyMapper",
    )
    accessions = _hgvs_call(lambda: mapper.relevant_transcripts(parsed), "hgvs could not list the transcripts")
    return {
        "status": "completed",
        "operation": "list_hgvs_relevant_transcripts",
        "variant": str(parsed),
        "transcript_count": len(accessions),
        "transcripts": [str(item) for item in accessions],
        "provenance": _provenance(method="hgvs.assemblymapper.AssemblyMapper.relevant_transcripts"),
    }


def hgvs_to_vcf_allele(
    variant: str,
    *,
    assembly_name: str = "GRCh38",
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a genomic (g.) HGVS variant to a VCF-style allele (chromosome, position, ref, alt).

    This is hgvs' experimental ``Babelfish.hgvs_to_vcf``: the variant is left-normalized, deletions and
    insertions are anchored on the preceding base, duplications become an insertion, and the chromosome
    uses the official name (no ``chr`` prefix) of ``assembly_name``. Reference bases are read through the
    data source.
    """

    babelfish = require("hgvs.extras.babelfish", _EXTRA)
    parsed = _parse_variant(variant)
    provider = _provider(data_source, local_data, uta_db_url)
    fish = _hgvs_call(lambda: babelfish.Babelfish(provider, _text(assembly_name, "assembly_name")), "hgvs could not create the Babelfish")
    chrom, position, ref, alt, kind = _hgvs_call(lambda: fish.hgvs_to_vcf(parsed), "hgvs could not convert the variant")
    return {
        "status": "completed",
        "operation": "hgvs_to_vcf_allele",
        "variant": str(parsed),
        "assembly": assembly_name,
        "chrom": str(chrom),
        "pos": int(position),
        "ref": str(ref),
        "alt": str(alt),
        "edit_type": str(kind),
        "provenance": _provenance(method="hgvs.extras.babelfish.Babelfish.hgvs_to_vcf"),
    }


# ---------------------------------------------------------------------------------------------
# Data provider queries
# ---------------------------------------------------------------------------------------------


def fetch_hgvs_reference_sequence(
    accession: str,
    *,
    start_i: int | None = None,
    end_i: int | None = None,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch a reference sequence by accession and zero-based interval ``[start_i, end_i)``.

    ``data_source="local"`` reads ``local_data['sequences']``; ``"uta"`` goes through the UTA data
    provider's sequence fetcher, which reads SeqRepo (``HGVS_SEQREPO_DIR`` or ``HGVS_SEQREPO_URL``) or
    falls back to NCBI/Ensembl over the network. Intervals are as in ``SeqFetcher.fetch_seq``.
    """

    ac = _text(accession, "accession")
    for label, value in (("start_i", start_i), ("end_i", end_i)):
        if value is not None:
            _integer(value, label, minimum=0, maximum=10**10)
    provider = _provider(data_source, local_data, uta_db_url)
    sequence = _hgvs_call(lambda: provider.get_seq(ac, start_i, end_i), "hgvs could not fetch the sequence")
    return {
        "status": "completed",
        "operation": "fetch_hgvs_reference_sequence",
        "accession": ac,
        "start_i": start_i,
        "end_i": end_i,
        "length": len(sequence),
        "sequence": sequence,
        "source": data_source,
        "provenance": _provenance(method="hgvs.dataproviders.seqfetcher.SeqFetcher.fetch_seq / provider.get_seq"),
    }


def query_hgvs_transcript_database(
    query: str,
    *,
    gene: str | None = None,
    transcript_accession: str | None = None,
    genomic_accession: str | None = None,
    alt_aln_method: str | None = None,
    assembly_name: str | None = None,
    data_source: str = "local",
    local_data: dict[str, Any] | None = None,
    uta_db_url: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query the transcript database behind hgvs (UTA schema) for genes, transcripts and alignments.

    ``query`` selects the data-provider call: ``gene_info`` (HGNC symbol ``gene``; ``get_gene_info``),
    ``tx_exons`` (``transcript_accession``, ``genomic_accession``, ``alt_aln_method``; exon structure on the
    genome alignment; ``get_tx_exons``), ``tx_for_gene`` (transcripts of ``gene``; ``get_tx_for_gene``),
    ``tx_identity_info`` (transcript lengths and CDS bounds; ``get_tx_identity_info``), ``tx_info``
    (alignment info for transcript, genomic accession and method; ``get_tx_info``), ``tx_mapping_options``
    (alignment options of ``transcript_accession``; ``get_tx_mapping_options``), ``similar_transcripts``
    (transcripts with identical CDS or exon structure; ``get_similar_transcripts``), ``protein_accession``
    (the protein of a transcript; ``get_pro_ac_for_tx_ac``) and ``assembly_map`` (the chromosome accessions
    of ``assembly_name`` and their names; ``get_assembly_map``, which reads hgvs' bundled assembly data).
    Rows are returned as dictionaries.
    """

    kind = _choice(query, "query", _QUERIES)

    def needs(value: str | None, name: str) -> str:
        if value is None:
            raise ConfigurationError(f"{name} is required for query={kind!r}.")
        return _text(value, name)

    provider = _provider(data_source, local_data, uta_db_url)
    if kind == "gene_info":
        symbol = needs(gene, "gene")
        row = _hgvs_call(lambda: provider.get_gene_info(symbol), "hgvs could not read the gene")
        rows = [] if row is None else _rows([row])
    elif kind == "tx_exons":
        tx = needs(transcript_accession, "transcript_accession")
        alt = needs(genomic_accession, "genomic_accession")
        method = needs(alt_aln_method, "alt_aln_method")
        rows = _rows(_hgvs_call(lambda: provider.get_tx_exons(tx, alt, method), "hgvs could not read the exons"))
    elif kind == "tx_for_gene":
        symbol = needs(gene, "gene")
        rows = _rows(_hgvs_call(lambda: provider.get_tx_for_gene(symbol), "hgvs could not read the transcripts"))
    elif kind == "tx_identity_info":
        tx = needs(transcript_accession, "transcript_accession")
        rows = _rows([_hgvs_call(lambda: provider.get_tx_identity_info(tx), "hgvs could not read the transcript")])
    elif kind == "tx_info":
        tx = needs(transcript_accession, "transcript_accession")
        alt = needs(genomic_accession, "genomic_accession")
        method = needs(alt_aln_method, "alt_aln_method")
        rows = _rows([_hgvs_call(lambda: provider.get_tx_info(tx, alt, method), "hgvs could not read the transcript info")])
    elif kind == "tx_mapping_options":
        tx = needs(transcript_accession, "transcript_accession")
        rows = _rows(_hgvs_call(lambda: provider.get_tx_mapping_options(tx), "hgvs could not read the mapping options"))
    elif kind == "similar_transcripts":
        tx = needs(transcript_accession, "transcript_accession")
        rows = _rows(_hgvs_call(lambda: provider.get_similar_transcripts(tx), "hgvs could not read the similar transcripts"))
    elif kind == "protein_accession":
        tx = needs(transcript_accession, "transcript_accession")
        value = _hgvs_call(lambda: provider.get_pro_ac_for_tx_ac(tx), "hgvs could not read the protein")
        rows = [] if value is None else [{"tx_ac": tx, "pro_ac": str(value)}]
    else:
        name = needs(assembly_name, "assembly_name")
        mapping = _hgvs_call(lambda: provider.get_assembly_map(name), "hgvs could not read the assembly map")
        rows = [{"accession": str(k), "name": str(v)} for k, v in sorted(mapping.items())][:_MAX_ROWS]
    return {
        "status": "completed",
        "operation": "query_hgvs_transcript_database",
        "query": kind,
        "source": data_source,
        "row_count": len(rows),
        "rows": rows,
        "provenance": _provenance(method=f"hgvs.dataproviders.uta.UTABase ({kind})"),
    }


__all__ = [
    "fetch_hgvs_reference_sequence",
    "fill_hgvs_reference",
    "hgvs_length_change",
    "hgvs_to_vcf_allele",
    "list_hgvs_relevant_transcripts",
    "map_hgvs_interval",
    "map_hgvs_variant",
    "parse_hgvs_position",
    "parse_hgvs_posedit",
    "predict_hgvs_protein_consequence",
    "project_hgvs_variant",
    "query_hgvs_transcript_database",
    "validate_hgvs_position",
    "validate_hgvs_variant",
]

_register_public_verbose(globals(), tuple(__all__))
