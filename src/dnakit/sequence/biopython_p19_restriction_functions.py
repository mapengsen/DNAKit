"""Biopython ``Bio.Restriction`` and ``Bio.CAPS`` operations added by integration pack P19.

Enzymes are named as in the REBASE data shipped with Biopython (``EcoRI``, ``BamHI``, ...). Positions
reported by ``Bio.Restriction`` are 1-based: the number is the position of the first base after the cut on
the top strand. A sequence is linear by default; pass ``linear=False`` for a circular molecule.

``RestrictionBatch.lambdasplit`` and ``RestrictionBatch.split`` return an empty batch in Biopython 1.88
(the result is stored in a private attribute instead of the set), so ``filter_restriction_enzymes``
filters the enzyme properties itself with the same predicates (``is_blunt``, ``is_5overhang``, ``size`` ...).
Biopython is imported only when a function runs.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import io
import re
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_MAX_SEQUENCE = 50_000_000
_MAX_ENZYMES = 2_000
_ENZYME_NAME = re.compile(r"[A-Za-z0-9]{2,20}")
_END_TYPES = ("blunt", "5_overhang", "3_overhang")
_RELATIONS = ("compatible_end", "neoschizomers", "equischizomers", "isoschizomers")
_VIEWS = ("list", "number", "map")
_POSITION_FILTERS = (
    "only_between",
    "between",
    "show_only_between",
    "only_outside",
    "outside",
    "do_not_cut",
)


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _bio(module: str) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return require(module, _EXTRA)


def _restriction() -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _bio("Bio.Restriction")


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _enzyme(name: object) -> Any:
    restriction = _restriction()
    if not isinstance(name, str) or not _ENZYME_NAME.fullmatch(name):
        raise ConfigurationError(f"{name!r} is not a restriction enzyme name.")
    candidate = getattr(restriction, name, None)
    if candidate is None or not isinstance(candidate, _bio("Bio.Restriction.Restriction").RestrictionType):
        raise ConfigurationError(f"{name!r} is not a restriction enzyme in the Biopython REBASE data.")
    return candidate


def _enzymes(names: object, label: str = "enzymes") -> list[Any]:
    if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or not names:
        raise ConfigurationError(f"{label} must be a non-empty list of enzyme names.")
    if len(names) > _MAX_ENZYMES:
        raise ConfigurationError(f"{label} may hold at most {_MAX_ENZYMES} enzymes.")
    seen: dict[str, Any] = {}
    for name in names:
        seen.setdefault(str(name), _enzyme(name))
    return list(seen.values())


def _dna(sequence: object) -> Any:
    seq = _bio("Bio.Seq")
    if not isinstance(sequence, str) or not sequence or len(sequence) > _MAX_SEQUENCE:
        raise ConfigurationError("sequence must be a non-empty DNA string.")
    if re.search(r"[^ACGTRYMKSWBDHVNacgtrymkswbdhvn]", sequence):
        raise ConfigurationError("sequence may contain only IUPAC nucleotide letters.")
    return seq.Seq(sequence.upper())


def _end_type(enzyme: Any) -> str:
    if enzyme.is_blunt():
        return "blunt"
    if enzyme.is_5overhang():
        return "5_overhang"
    if enzyme.is_3overhang():
        return "3_overhang"
    return "undefined"


def _batch(restriction: Any, enzymes: Sequence[Any]) -> Any:
    return restriction.RestrictionBatch(list(enzymes))


def restriction_enzyme_properties(
    enzymes: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Describe restriction enzymes: site, cut pattern, end type, overhang, frequency and suppliers.

    For each enzyme (Biopython ``Restriction.<Enzyme>``) returns the recognition ``site`` and ``size``, the
    cut pattern from ``elucidate()`` (for example ``G^AATT_C``; ``^`` is the top-strand cut, ``_`` the
    bottom-strand cut), the end type (``blunt``, ``5_overhang`` or ``3_overhang`` from ``is_blunt``,
    ``is_5overhang``, ``is_3overhang``), the overhang length and sequence, whether the site is
    palindromic, ``frequency`` (the site occurs once per this many random bases), whether the enzyme is
    commercially available and its suppliers (codes and names).
    """

    records = []
    for enzyme in _enzymes(enzymes):
        suppliers = list(enzyme.supplier_list()) if enzyme.is_comm() else []
        records.append(
            {
                "name": str(enzyme),
                "site": enzyme.site,
                "size": int(enzyme.size),
                "elucidate": enzyme.elucidate(),
                "end_type": _end_type(enzyme),
                "overhang": int(enzyme.ovhg),
                "overhang_sequence": str(enzyme.ovhgseq),
                "palindromic": bool(enzyme.is_palindromic()),
                "frequency": float(enzyme.frequency()),
                "commercial": bool(enzyme.is_comm()),
                "supplier_codes": list(enzyme.suppl) if enzyme.is_comm() else [],
                "suppliers": suppliers,
            }
        )
    return {
        "status": "completed",
        "operation": "restriction_enzyme_properties",
        "enzyme_count": len(records),
        "enzymes": records,
        "provenance": _provenance(
            parser="Bio.Restriction enzyme classes (elucidate, is_blunt, is_5overhang, is_3overhang, "
            "frequency, suppl, supplier_list, is_comm)"
        ),
    }


def restriction_related_enzymes(
    enzyme: str,
    relation: str = "isoschizomers",
    *,
    among: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find enzymes related to ``enzyme`` (``Bio.Restriction.<Enzyme>.<relation>``).

    ``relation`` is ``compatible_end`` (enzymes producing ends that can be ligated to this enzyme's ends),
    ``neoschizomers`` (same site, different cut position), ``equischizomers`` (same site and cut position)
    or ``isoschizomers`` (both). The search runs over every REBASE enzyme in Biopython unless ``among``
    lists the enzyme names to consider.
    """

    chosen = _choice(relation, _RELATIONS, "relation")
    subject = _enzyme(enzyme)
    restriction = _restriction()
    batch = _batch(restriction, _enzymes(among, "among")) if among is not None else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        found = getattr(subject, chosen)(batch) if batch is not None else getattr(subject, chosen)()
    names = sorted(str(item) for item in found)
    return {
        "status": "completed",
        "operation": "restriction_related_enzymes",
        "enzyme": str(subject),
        "relation": chosen,
        "site": subject.site,
        "count": len(names),
        "enzymes": names,
        "provenance": _provenance(parser=f"Bio.Restriction.<Enzyme>.{chosen}"),
    }


def restriction_enzymes_by_supplier(
    supplier: str | None = None, *, verbose: bool = True
) -> dict[str, Any]:
    """List the enzymes sold by a supplier (``RestrictionBatch.add_supplier``) or list the suppliers.

    ``supplier`` is a one-letter REBASE supplier code (``N`` = New England Biolabs) or the full supplier
    name. Without ``supplier`` the function returns the table of supplier codes and names with the number of
    enzymes each sells.
    """

    restriction = _restriction()
    table = restriction.RestrictionBatch.suppl_codes()
    if supplier is None:
        counts = {}
        for code in table:
            batch = restriction.RestrictionBatch([])
            batch.add_supplier(code)
            counts[code] = len(batch)
        return {
            "status": "completed",
            "operation": "restriction_enzymes_by_supplier",
            "suppliers": [
                {"code": code, "name": name, "enzyme_count": counts[code]} for code, name in table.items()
            ],
            "provenance": _provenance(parser="Bio.Restriction.RestrictionBatch.suppl_codes"),
        }
    if not isinstance(supplier, str) or not supplier:
        raise ConfigurationError("supplier must be a supplier code or name.")
    code = supplier if supplier in table else next((k for k, v in table.items() if v == supplier), None)
    if code is None:
        raise ConfigurationError(
            f"{supplier!r} is not a supplier; codes are {sorted(table)} (names: {sorted(table.values())})."
        )
    batch = restriction.RestrictionBatch([])
    batch.add_supplier(code)
    return {
        "status": "completed",
        "operation": "restriction_enzymes_by_supplier",
        "supplier_code": code,
        "supplier": table[code],
        "enzyme_count": len(batch),
        "enzymes": sorted(str(item) for item in batch),
        "provenance": _provenance(parser="Bio.Restriction.RestrictionBatch.add_supplier"),
    }


def filter_restriction_enzymes(
    enzymes: Sequence[str] | None = None,
    *,
    end_type: str | None = None,
    commercial: bool | None = None,
    supplier: str | None = None,
    palindromic: bool | None = None,
    min_site_size: int | None = None,
    max_site_size: int | None = None,
    defined_site: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select enzymes by property (the job of ``RestrictionBatch.lambdasplit`` / ``split``).

    Starts from ``enzymes`` (default: every enzyme in Biopython's REBASE data) and keeps those that match all
    given criteria: ``end_type`` (``blunt``, ``5_overhang`` or ``3_overhang``), ``commercial`` availability,
    a one-letter ``supplier`` code, a ``palindromic`` site, the recognition-site size range
    (``min_site_size`` / ``max_site_size``) and ``defined_site`` (site without ambiguity codes).
    Biopython 1.88's own ``lambdasplit`` and ``split`` return empty batches, so the predicates are applied
    here.
    """

    restriction = _restriction()
    pool = _enzymes(enzymes) if enzymes is not None else list(restriction.AllEnzymes)
    if end_type is not None:
        _choice(end_type, _END_TYPES, "end_type")
    for value, label in ((commercial, "commercial"), (palindromic, "palindromic"), (defined_site, "defined_site")):
        if value is not None:
            _flag(value, label)
    low = None if min_site_size is None else _integer(min_site_size, "min_site_size", minimum=1, maximum=100)
    high = None if max_site_size is None else _integer(max_site_size, "max_site_size", minimum=1, maximum=100)
    if supplier is not None and supplier not in restriction.RestrictionBatch.suppl_codes():
        raise ConfigurationError(
            f"supplier must be one of the REBASE codes {sorted(restriction.RestrictionBatch.suppl_codes())}."
        )

    def keep(enzyme: Any) -> bool:
        return (
            (end_type is None or _end_type(enzyme) == end_type)
            and (commercial is None or bool(enzyme.is_comm()) == commercial)
            and (supplier is None or (enzyme.is_comm() and supplier in enzyme.suppl))
            and (palindromic is None or bool(enzyme.is_palindromic()) == palindromic)
            and (low is None or enzyme.size >= low)
            and (high is None or enzyme.size <= high)
            and (defined_site is None or bool(enzyme.is_defined()) == defined_site)
        )

    chosen = sorted(str(item) for item in pool if keep(item))
    return {
        "status": "completed",
        "operation": "filter_restriction_enzymes",
        "candidate_count": len(pool),
        "count": len(chosen),
        "enzymes": chosen,
        "provenance": _provenance(
            parser="Bio.Restriction enzyme predicates (is_blunt, is_5overhang, is_3overhang, is_comm, "
            "is_palindromic, size, is_defined)",
            reason="RestrictionBatch.lambdasplit/split are empty in Biopython 1.88",
        ),
    }


def _analysis(sequence: str, enzymes: Sequence[str], linear: bool) -> Any:
    restriction = _restriction()
    dna = _dna(sequence)
    batch = _batch(restriction, _enzymes(enzymes))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return restriction.Analysis(batch, dna, linear=_flag(linear, "linear"))


def _cuts(result: Mapping[Any, Sequence[int]]) -> dict[str, list[int]]:
    return {str(enzyme): [int(position) for position in positions] for enzyme, positions in sorted(result.items(), key=lambda item: str(item[0]))}


def restriction_analysis_report(
    sequence: str,
    enzymes: Sequence[str],
    *,
    linear: bool = True,
    view: str = "list",
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a sequence for restriction sites and print the result as text (``Analysis.print_as``).

    ``view`` is ``list`` (one line per enzyme, alphabetical, with the cut positions), ``number`` (number of
    sites per enzyme) or ``map`` (the sequence with the enzyme names written at their cut positions).
    Returns the text report, the cut positions (1-based) of every enzyme that cuts and the names of the
    enzymes that do not cut (``no_cut``).
    """

    chosen = _choice(view, _VIEWS, "view")
    analysis = _analysis(sequence, enzymes, linear)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        found = analysis.full()
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            analysis.print_as(chosen)
            analysis.print_that()
    text = buffer.getvalue()
    return {
        "status": "completed",
        "operation": "restriction_analysis_report",
        "view": chosen,
        "linear": linear,
        "sequence_length": len(sequence),
        "cuts": _cuts({enzyme: positions for enzyme, positions in found.items() if positions}),
        "no_cut": sorted(str(enzyme) for enzyme, positions in found.items() if not positions),
        "report": text,
        "provenance": _provenance(parser="Bio.Restriction.Restriction.Analysis.print_as / print_that"),
    }


def filter_restriction_results(
    sequence: str,
    enzymes: Sequence[str],
    *,
    linear: bool = True,
    end_type: str | None = None,
    cut_count: int | Sequence[int] | None = None,
    cuts: str | None = None,
    names: Sequence[str] | None = None,
    site_sizes: Sequence[int] | None = None,
    position_filter: str | None = None,
    start: int | None = None,
    end: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a restriction search and keep only the results that match the criteria.

    All given criteria are combined (each narrows the previous result) using the ``Bio.Restriction.Analysis``
    methods of the same names:

    * ``end_type``: ``blunt`` / ``5_overhang`` / ``3_overhang`` (``blunt``, ``overhang5``, ``overhang3``);
      as in Biopython, enzymes of that end type that do not cut are listed with an empty position list
      (add ``cuts="cutting"`` to drop them);
    * ``cut_count``: exact number of sites, or a list of allowed numbers (``with_N_sites``,
      ``with_number_list``); ``cuts="cutting"`` keeps enzymes that cut (``with_sites``) and
      ``cuts="non_cutting"`` those that do not (``without_site``);
    * ``names``: enzyme names (``with_name``); ``site_sizes``: recognition-site lengths (``with_site_size``;
      these two are applied directly because the two Biopython 1.88 methods do not accept names or sizes
      together with an earlier result);
    * ``position_filter`` with 1-based ``start`` and ``end``: ``only_between`` (all cuts inside the range),
      ``between`` (at least one cut inside), ``show_only_between`` (keep only the cuts inside),
      ``only_outside`` / ``outside`` (the same outside the range) or ``do_not_cut`` (enzymes with no cut
      in the range).
    """

    analysis = _analysis(sequence, enzymes, linear)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        current = dict(analysis.full())
        applied = []
        # Analysis methods treat an empty dictionary as "no dictionary given" and search again, so a
        # step that leaves nothing ends the chain.
        if end_type is not None:
            method = {"blunt": "blunt", "5_overhang": "overhang5", "3_overhang": "overhang3"}[
                _choice(end_type, _END_TYPES, "end_type")
            ]
            current = getattr(analysis, method)(current) if current else {}
            applied.append(f"end_type={end_type}")
        if cut_count is not None:
            if isinstance(cut_count, int) and not isinstance(cut_count, bool):
                number = _integer(cut_count, "cut_count", minimum=0, maximum=10**6)
                current = analysis.with_N_sites(number, current) if current else {}
            elif isinstance(cut_count, Sequence) and not isinstance(cut_count, (str, bytes)) and cut_count:
                values = [_integer(item, "cut_count", minimum=0, maximum=10**6) for item in cut_count]
                current = analysis.with_number_list(values, current) if current else {}
            else:
                raise ConfigurationError("cut_count must be an integer or a non-empty list of integers.")
            applied.append("cut_count")
        if cuts is not None:
            kind = _choice(cuts, ("cutting", "non_cutting"), "cuts")
            if current:
                current = analysis.with_sites(current) if kind == "cutting" else analysis.without_site(current)
            applied.append(f"cuts={kind}")
        if names is not None:
            if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or not names:
                raise ConfigurationError("names must be a non-empty list of enzyme names.")
            wanted = {str(_enzyme(name)) for name in names}
            current = {key: value for key, value in current.items() if str(key) in wanted}
            applied.append("names")
        if site_sizes is not None:
            if isinstance(site_sizes, (str, bytes)) or not isinstance(site_sizes, Sequence) or not site_sizes:
                raise ConfigurationError("site_sizes must be a non-empty list of integers.")
            sizes = {_integer(item, "site_sizes", minimum=1, maximum=100) for item in site_sizes}
            current = {key: value for key, value in current.items() if key.size in sizes}
            applied.append("site_sizes")
        if position_filter is not None:
            chosen = _choice(position_filter, _POSITION_FILTERS, "position_filter")
            first = _integer(start, "start", minimum=1, maximum=len(sequence))
            last = _integer(end, "end", minimum=first, maximum=len(sequence))
            current = getattr(analysis, chosen)(first, last, current) if current else {}
            applied.append(f"{chosen}({first},{last})")
    if not applied:
        raise ConfigurationError("Give at least one filter criterion.")
    return {
        "status": "completed",
        "operation": "filter_restriction_results",
        "filters": applied,
        "match_count": len(current),
        "results": _cuts(current),
        "provenance": _provenance(
            parser="Bio.Restriction.Restriction.Analysis (blunt, overhang5, overhang3, with_N_sites, "
            "with_number_list, with_sites, without_site, with_name, with_site_size, only_between, between, ...)"
        ),
    }


def caps_differential_cut_sites(
    alignment: Mapping[str, str],
    enzymes: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find CAPS marker candidates: restriction sites that cut some aligned sequences but not others.

    This is ``Bio.CAPS.CAPSMap``. For every enzyme, a site is reported when at least one aligned sequence
    is cut at that alignment position and at least one is not (cleaved amplified polymorphic sequence
    markers, Konieczny and Ausubel 1993). ``alignment`` maps sequence name to aligned DNA of equal length
    (IUPAC letters; gap characters are not accepted because ``Bio.Restriction`` cannot search them: use
    ``N`` or remove gap columns first). ``start`` is the 0-based alignment position of the cut; ``cuts_in`` and ``blocked_in``
    name the sequences.
    """

    if not isinstance(alignment, Mapping) or len(alignment) < 2:
        raise ConfigurationError("alignment must map at least two sequence names to aligned DNA.")
    names = list(alignment)
    rows = []
    for name in names:
        text = alignment[name]
        if not isinstance(text, str) or not text or re.search(r"[^ACGTRYMKSWBDHVNacgtrymkswbdhvn]", text):
            raise ConfigurationError(f"sequence {name!r} must be DNA text of IUPAC letters (no gap characters).")
        if len(text) > _MAX_SEQUENCE:
            raise ConfigurationError("alignment sequences are too long.")
        rows.append(text.upper())
    if len({len(text) for text in rows}) != 1:
        raise ConfigurationError("all aligned sequences must have the same length.")
    objects = _enzymes(enzymes)
    caps = _bio("Bio.CAPS")
    align = _bio("Bio.Align")
    seq = _bio("Bio.Seq")
    record = _bio("Bio.SeqRecord")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        msa = align.MultipleSeqAlignment(
            [record.SeqRecord(seq.Seq(text), id=name) for name, text in zip(names, rows)]
        )
        marker_map = caps.CAPSMap(msa, objects)
    sites = [
        {
            "enzyme": str(site.enzyme),
            "start": int(site.start),
            "cuts_in": [names[i] for i in site.cuts_in],
            "blocked_in": [names[i] for i in site.blocked_in],
        }
        for site in marker_map.dcuts
    ]
    return {
        "status": "completed",
        "operation": "caps_differential_cut_sites",
        "sequence_count": len(names),
        "alignment_length": len(rows[0]),
        "site_count": len(sites),
        "sites": sites,
        "provenance": _provenance(parser="Bio.CAPS.CAPSMap"),
    }


__all__ = [
    "caps_differential_cut_sites",
    "filter_restriction_enzymes",
    "filter_restriction_results",
    "restriction_analysis_report",
    "restriction_enzyme_properties",
    "restriction_enzymes_by_supplier",
    "restriction_related_enzymes",
]

_register_public_verbose(globals(), tuple(__all__))
