"""pydna wrappers added by integration pack P34 (molecular_lab area).

Enzymatic end treatments of double-stranded DNA with overhangs, restriction cutter analysis, strand
melting, assembly and PCR figures, restriction-ligation assembly, OpenCloning strategies and cloning
history, primer screening, recombinase sites, codon-usage values and small sequence utilities.

Sequences are plain JSON values. A double-stranded molecule is either a string (a blunt linear
molecule) or an object ``{"watson": ..., "crick": ..., "ovhg": int, "circular": bool}`` written as
pydna's ``Dseq`` takes it: ``watson`` is the top strand 5'->3', ``crick`` the bottom strand 5'->3',
``ovhg`` the overhang of the crick strand (pydna convention). A cloning fragment is a string or
``{"sequence": ..., "circular": bool, "name": str}``. pydna is imported only when a function runs.
Nothing here uses the network or the shell.
"""

from __future__ import annotations

import datetime
import importlib.metadata
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-assembly"
_PACKAGE = "pydna"
_MAX_BP = 2_000_000
_MAX_ITEMS = 500
_ERRORS = (ValueError, TypeError, AssertionError, KeyError, IndexError, AttributeError)
_DNA_TEXT = re.compile(r"^[A-Za-z]+$")
_ORGANISMS = ("sce",)
_CUTTER_CLASSES = ("none", "once", "twice", "n", "any")
_ASSEMBLY_FIGURES = ("compact", "detailed", "matplotlib")
_END_SIDES = ("left", "right")
_STRATEGY_OPERATIONS = ("none", "pcr", "gibson", "restriction_ligation")


# Validation helpers.


def _provenance() -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed}


def _integer(value: object, name: str, minimum: int = 0, maximum: int = _MAX_BP) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return value


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _text(value: object, name: str, *, maximum: int = _MAX_BP) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum} characters.")
    return value.strip()


def _letters(value: object, name: str, *, maximum: int = _MAX_BP) -> str:
    text = _text(value, name, maximum=maximum)
    if not _DNA_TEXT.match(text):
        raise ConfigurationError(f"{name} must contain letters only.")
    return text


def _nucleotides(value: object, name: str, default: str | None = None) -> str | None:
    if value is None:
        return default
    text = _text(value, name, maximum=8).lower()
    if set(text) - set("acgt"):
        raise ConfigurationError(f"{name} may contain only A, C, G and T.")
    return text


def _strand(value: object, name: str) -> str:
    """A strand written as pydna prints it: letters, with a space for a gap in the strand."""

    if not isinstance(value, str) or not value.strip(" "):
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > _MAX_BP:
        raise ConfigurationError(f"{name} must be at most {_MAX_BP} characters.")
    if not re.match(r"^[A-Za-z][A-Za-z ]*$", value) or value.endswith(" "):
        raise ConfigurationError(f"{name} must be letters, with spaces only inside for gaps.")
    return value


def _dseq(value: object, name: str) -> Any:
    """A pydna ``Dseq`` from a string or a ``{watson, crick, ovhg, circular}`` object."""

    dseq_module = require("pydna.dseq", _EXTRA)
    if isinstance(value, str):
        return dseq_module.Dseq(_letters(value, name))
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a DNA string or an object with a watson strand.")
    unknown = set(value) - {"watson", "crick", "ovhg", "circular"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    watson = _strand(value.get("watson"), f"{name}.watson")
    crick_value = value.get("crick")
    crick = None if crick_value is None else _strand(crick_value, f"{name}.crick")
    ovhg_value = value.get("ovhg")
    ovhg = None if ovhg_value is None else _integer(ovhg_value, f"{name}.ovhg", -_MAX_BP)
    circular = _flag(value.get("circular", False), f"{name}.circular")
    try:
        return dseq_module.Dseq(watson, crick, ovhg=ovhg, circular=circular)
    except _ERRORS as exc:
        raise ConfigurationError(f"pydna rejected {name}: {exc}") from exc


def _dseq_dict(dseq: Any) -> dict[str, Any]:
    return {
        "watson": str(dseq.watson),
        "crick": str(dseq.crick),
        "ovhg": None if dseq.ovhg is None else int(dseq.ovhg),
        "circular": bool(dseq.circular),
        "length": len(dseq),
        "is_blunt": bool(dseq.isblunt()),
        "representation": repr(dseq),
    }


def _record(value: object, name: str, index: int = 0) -> Any:
    """A pydna ``Dseqrecord`` from a string or ``{sequence|watson..., circular, name}``."""

    record_module = require("pydna.dseqrecord", _EXTRA)
    label = f"{name}[{index}]"
    if isinstance(value, str):
        rec = record_module.Dseqrecord(_letters(value, label))
        rec.name = f"fragment{index + 1}"
        return rec
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{label} must be a DNA string or an object.")
    unknown = set(value) - {"sequence", "watson", "crick", "ovhg", "circular", "name"}
    if unknown:
        raise ConfigurationError(f"{label} has unsupported keys: {sorted(unknown)}.")
    if "sequence" in value:
        circular = _flag(value.get("circular", False), f"{label}.circular")
        rec = record_module.Dseqrecord(_letters(value["sequence"], f"{label}.sequence"))
        if circular:
            rec = record_module.Dseqrecord(
                _letters(value["sequence"], f"{label}.sequence"), circular=True
            )
    else:
        rec = record_module.Dseqrecord(_dseq({k: v for k, v in value.items() if k != "name"}, label))
    rec.name = str(value.get("name") or f"fragment{index + 1}")[:40]
    return rec


def _records(values: object, name: str, minimum: int = 1) -> list[Any]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if not minimum <= len(values) <= _MAX_ITEMS:
        raise ConfigurationError(f"{name} must hold {minimum} to {_MAX_ITEMS} items.")
    return [_record(item, name, index) for index, item in enumerate(values)]


def _primers(values: object, name: str) -> list[Any]:
    primer_module = require("pydna.primer", _EXTRA)
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(values) > _MAX_ITEMS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_ITEMS} primers.")
    result = []
    for index, item in enumerate(values):
        label = f"{name}[{index}]"
        if isinstance(item, str):
            sequence, primer_name = _letters(item, label, maximum=500), f"primer{index + 1}"
        elif isinstance(item, Mapping) and "sequence" in item:
            sequence = _letters(item["sequence"], f"{label}.sequence", maximum=500)
            primer_name = str(item.get("name") or f"primer{index + 1}")[:40]
        else:
            raise ConfigurationError(f"{label} must be a DNA string or {{name, sequence}}.")
        result.append(primer_module.Primer(sequence, id=primer_name, name=primer_name))
    return result


def _enzymes(names: object, name: str) -> list[Any]:
    restriction = require("Bio.Restriction", "bio-core")
    if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or not names:
        raise ConfigurationError(f"{name} must be a non-empty list of enzyme names.")
    result = []
    for item in names:
        enzyme = None
        if isinstance(item, str) and re.fullmatch(r"[A-Za-z0-9]{2,20}", item):
            try:
                enzyme = restriction.AllEnzymes.get(item)
            except ValueError:
                enzyme = None
        if enzyme is None:
            raise ConfigurationError(f"{name} contains an unknown restriction enzyme: {item!r}.")
        result.append(enzyme)
    return result


def _now() -> str:
    return datetime.datetime.now().replace(microsecond=0).isoformat()


def _guard(operation: str, call: Any) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except _ERRORS as exc:
        raise ConfigurationError(f"pydna could not run {operation}: {exc}") from exc


def _done(operation: str, **fields: Any) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": _provenance()}


# End treatments.


def pydna_fill_in_ends(
    dna: str | dict[str, Any], *, nucleotides: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Fill in recessed 3' ends (5' overhangs) with a polymerase that has no exonuclease activity.

    This is the exo-minus Klenow fill-in of pydna ``Dseq.fill_in``. ``nucleotides`` limits the
    dNTPs that are present (any combination of A, C, G, T); None means all four. ``dna`` is a string
    or a ``{watson, crick, ovhg, circular}`` object. The result holds the new molecule.
    """

    molecule = _dseq(dna, "dna")
    letters = _nucleotides(nucleotides, "nucleotides")
    new = _guard("fill_in", lambda: molecule.fill_in(letters))
    return _done("fill_in", input=_dseq_dict(molecule), result=_dseq_dict(new), nucleotides=letters)


def pydna_mung_bean_blunt(dna: str | dict[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Blunt a molecule by removing single-stranded overhangs, as mung bean nuclease does.

    Covers pydna ``Dseq.nibble_to_blunt`` and its alias ``Dseq.mung``: 5' and 3' overhangs are
    both degraded, so the result is blunt whatever kind of overhang the input had.
    """

    molecule = _dseq(dna, "dna")
    new = _guard("nibble_to_blunt", lambda: molecule.nibble_to_blunt())
    return _done("mung_bean_blunt", input=_dseq_dict(molecule), result=_dseq_dict(new))


def pydna_t4_polymerase_polish(
    dna: str | dict[str, Any], *, nucleotides: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Polish ends with T4 DNA polymerase: fill in 5' overhangs and chew back 3' overhangs.

    This is pydna ``Dseq.T4``. ``nucleotides`` restricts the dNTPs present (A, C, G, T in any
    combination); None means all four. With fewer dNTPs the exonuclease activity stops at the first
    base that could be refilled, so partial sticky ends can be produced.
    """

    molecule = _dseq(dna, "dna")
    letters = _nucleotides(nucleotides, "nucleotides")
    new = _guard("T4", lambda: molecule.T4(letters))
    return _done("t4_polymerase_polish", input=_dseq_dict(molecule), result=_dseq_dict(new),
                 nucleotides=letters)


def pydna_exonuclease_five_prime(
    dna: str | dict[str, Any], *, end: str = "left", bases: int = 1, verbose: bool = True
) -> dict[str, Any]:
    """Resect one end of a molecule with a 5'->3' exonuclease, removing ``bases`` nucleotides.

    ``end`` is ``left`` (pydna ``Dseq.nibble_five_prime_left``) or ``right``
    (``Dseq.nibble_five_prime_right``). The 5' strand of that end is shortened, which exposes a 3'
    overhang.
    """

    molecule = _dseq(dna, "dna")
    side = _choice(end, _END_SIDES, "end")
    count = _integer(bases, "bases", 0, 100_000)
    method = f"nibble_five_prime_{side}"
    new = _guard(method, lambda: getattr(molecule, method)(count))
    return _done("exonuclease_five_prime", end=side, bases=count, input=_dseq_dict(molecule),
                 result=_dseq_dict(new))


def pydna_exonuclease_three_prime(
    dna: str | dict[str, Any], *, end: str = "left", bases: int = 1, verbose: bool = True
) -> dict[str, Any]:
    """Resect one end of a molecule with a 3'->5' exonuclease, removing ``bases`` nucleotides.

    ``end`` is ``left`` (pydna ``Dseq.nibble_three_prime_left``) or ``right``
    (``Dseq.nibble_three_prime_right``). The 3' strand of that end is shortened, which exposes a 5'
    overhang.
    """

    molecule = _dseq(dna, "dna")
    side = _choice(end, _END_SIDES, "end")
    count = _integer(bases, "bases", 0, 100_000)
    method = f"nibble_three_prime_{side}"
    new = _guard(method, lambda: getattr(molecule, method)(count))
    return _done("exonuclease_three_prime", end=side, bases=count, input=_dseq_dict(molecule),
                 result=_dseq_dict(new))


def pydna_terminal_transferase_tail(
    dna: str | dict[str, Any], *, nucleotides: str = "a", verbose: bool = True
) -> dict[str, Any]:
    """Add a tail of nucleotides to the 3' ends, as terminal transferase (TdT, Taq tailing) does.

    This is pydna ``Dseq.terminal_transferase`` (and ``Dseqrecord.terminal_transferase``).
    ``nucleotides`` is the string of bases that is added to each 3' end; the default ``a`` gives the
    single A overhang used in T/A cloning.
    """

    molecule = _dseq(dna, "dna")
    letters = _nucleotides(nucleotides, "nucleotides", "a")
    new = _guard("terminal_transferase", lambda: molecule.terminal_transferase(letters))
    return _done("terminal_transferase_tail", nucleotides=letters, input=_dseq_dict(molecule),
                 result=_dseq_dict(new))


def pydna_user_enzyme_treatment(dna: str | dict[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Simulate USER enzyme treatment: each uracil is excised and leaves a gap in its strand.

    This is pydna ``Dseq.user``. Write the uracil as ``U`` in the strand; short stretches
    between the gap and the end are released, which creates a 3' overhang.
    """

    molecule = _dseq(dna, "dna")
    new = _guard("user", lambda: molecule.user())
    return _done("user_enzyme_treatment", input=_dseq_dict(molecule), result=_dseq_dict(new))


def pydna_cast_to_double_strand(
    dna: str | dict[str, Any], *, end: str = "both", verbose: bool = True
) -> dict[str, Any]:
    """Convert single-stranded ends into double-stranded sequence by adding the complement.

    ``end`` is ``both`` (pydna ``Dseq.cast_to_ds``), ``left`` (``Dseq.cast_to_ds_left``) or
    ``right`` (``Dseq.cast_to_ds_right``). The result is blunt on the converted ends.
    """

    molecule = _dseq(dna, "dna")
    side = _choice(end, ("both", "left", "right"), "end")
    method = {"both": "cast_to_ds", "left": "cast_to_ds_left", "right": "cast_to_ds_right"}[side]
    new = _guard(method, lambda: getattr(molecule, method)())
    return _done("cast_to_double_strand", end=side, input=_dseq_dict(molecule),
                 result=_dseq_dict(new))


def pydna_melt_short_strands(
    dna: str | dict[str, Any], *, length: int, verbose: bool = True
) -> dict[str, Any]:
    """Release short single-stranded stretches from a duplex and report where they are.

    Single-stranded pieces up to ``length`` nucleotides that are separated by nicks or gaps are shed
    without changing the length of the remaining molecule. This covers pydna
    ``Dseq.get_ss_meltsites`` (the ``sites``), ``Dseq.melt`` (``melt_fragments``: the shed strands and the duplex pieces left when short double-
    stranded stretches also melt) and ``Dseq.melt_ss_dna`` (``remaining`` and ``shed``).
    """

    molecule = _dseq(dna, "dna")
    limit = _integer(length, "length", 1, 100_000)

    def run() -> dict[str, Any]:
        sites = molecule.get_ss_meltsites(limit)
        remaining, strands = molecule.melt_ss_dna(limit)
        return {
            "watson_sites": [list(map(int, pair)) for pair in sites[0]],
            "crick_sites": [list(map(int, pair)) for pair in sites[1]],
            "remaining": _dseq_dict(remaining),
            "shed": [_dseq_dict(item) for item in strands],
            "melt_fragments": [_dseq_dict(item) for item in molecule.melt(limit)],
        }

    found = _guard("melt_ss_dna", run)
    return _done("melt_short_strands", length=limit, input=_dseq_dict(molecule), **found)


def pydna_shed_strand_segments(
    dna: str | dict[str, Any],
    *,
    watson_cutpairs: list[list[int]] | None = None,
    crick_cutpairs: list[list[int]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Separate chosen segments of the top (watson) and bottom (crick) strand from a molecule.

    Each cut pair is ``[start, end]`` in the coordinates pydna ``Dseq.shed_ss_dna`` uses. The result
    has the new molecule and the single-stranded pieces that were shed.
    """

    molecule = _dseq(dna, "dna")

    def pairs(value: object, name: str) -> list[tuple[int, int]]:
        if value is None:
            return []
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise ConfigurationError(f"{name} must be a list of [start, end] pairs.")
        out = []
        for index, item in enumerate(value):
            if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 2:
                raise ConfigurationError(f"{name}[{index}] must be [start, end].")
            out.append((_integer(item[0], f"{name}[{index}][0]"), _integer(item[1], f"{name}[{index}][1]")))
        return out

    watson_pairs = pairs(watson_cutpairs, "watson_cutpairs")
    crick_pairs = pairs(crick_cutpairs, "crick_cutpairs")
    new, strands = _guard("shed_ss_dna", lambda: molecule.shed_ss_dna(watson_pairs, crick_pairs))
    return _done("shed_strand_segments", input=_dseq_dict(molecule), result=_dseq_dict(new),
                 shed=[_dseq_dict(item) for item in strands])


# Restriction analysis.


def pydna_restriction_cutters(
    dna: str | dict[str, Any],
    *,
    cutter_class: str = "once",
    cuts: int = 1,
    enzymes: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the restriction enzymes that cut a molecule a given number of times.

    ``cutter_class`` selects the group: ``none`` (pydna ``Dseq.no_cutters``), ``once``
    (``unique_cutters`` / ``once_cutters``), ``twice`` (``twice_cutters``), ``n`` (``n_cutters``
    with ``cuts`` sites) or ``any`` (``cutters``, at least one site). ``enzymes`` is a list of
    Biopython enzyme names; the default is pydna's batch of commercially available enzymes
    (``CommOnly``). The result also holds the number of cuts of every selected enzyme.
    """

    molecule = _dseq(dna, "dna")
    chosen = _choice(cutter_class, _CUTTER_CLASSES, "cutter_class")
    count = _integer(cuts, "cuts", 0, 10_000)
    restriction = require("Bio.Restriction", "bio-core")
    batch = restriction.CommOnly if enzymes is None else restriction.RestrictionBatch(
        _enzymes(enzymes, "enzymes")
    )
    method = {
        "none": "no_cutters", "once": "unique_cutters", "twice": "twice_cutters",
        "n": "n_cutters", "any": "cutters",
    }[chosen]

    def run() -> dict[str, Any]:
        selected = (
            molecule.n_cutters(count, batch) if chosen == "n" else getattr(molecule, method)(batch)
        )
        analysis = batch.search(molecule)
        counts = {str(e): len(analysis[e]) for e in sorted(selected, key=str)}
        return {"enzymes": sorted(counts), "cut_counts": counts}

    found = _guard(method, run)
    return _done("restriction_cutters", cutter_class=chosen,
                 cuts=count if chosen == "n" else {"none": 0, "once": 1, "twice": 2}.get(chosen),
                 batch_size=len(batch), **found)


def pydna_restriction_ligation_assembly(
    fragments: list[Any],
    enzymes: list[str],
    *,
    allow_blunt: bool = True,
    circular_only: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Digest fragments with restriction enzymes, ligate the compatible ends, and flag partial digests.

    Runs pydna ``assembly2.restriction_ligation_assembly`` and, for every product, pydna
    ``assembly2.partially_digested_seqs_in_restriction_ligation_assembly``: the input fragments
    that took part in a product without being cut between two sites of the enzymes. Fragments are
    strings or ``{sequence, circular, name}`` objects; ``enzymes`` are Biopython names. With
    ``circular_only`` only circular products are returned.
    """

    frags = _records(fragments, "fragments", 1)
    enzyme_list = _enzymes(enzymes, "enzymes")
    blunt = _flag(allow_blunt, "allow_blunt")
    only_circular = _flag(circular_only, "circular_only")
    assembly = require("pydna.assembly2", _EXTRA)

    def run() -> list[dict[str, Any]]:
        products = assembly.restriction_ligation_assembly(frags, enzyme_list, blunt, only_circular)
        rows = []
        for product in products:
            partial = assembly.partially_digested_seqs_in_restriction_ligation_assembly(product)
            rows.append({
                "sequence": str(product.seq.watson),
                "length": len(product),
                "circular": bool(product.circular),
                "partially_digested_inputs": [str(item.name) for item in partial],
            })
        return rows

    rows = _guard("restriction_ligation_assembly", run)
    return _done("restriction_ligation_assembly", enzymes=[str(e) for e in enzyme_list],
                 product_count=len(rows), products=rows)


# Figures.


def pydna_pcr_alignment_figure(
    template: str | dict[str, Any],
    forward_primer: str,
    reverse_primer: str,
    *,
    limit: int = 13,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run an in-silico PCR and draw the primer-template alignment as text.

    Uses pydna ``amplify.pcr`` and ``Amplicon.figure``. The figure shows the forward primer on the
    top strand and the reverse primer on the bottom strand, with tails marked. ``limit`` is the
    number of 3' bases that must anneal (pydna default 13). Exactly one product must form.
    """

    amplify = require("pydna.amplify", _EXTRA)
    tmpl = _record(template, "template")
    primers = _primers([forward_primer, reverse_primer], "primers")
    min_anneal = _integer(limit, "limit", 4, 100)

    def run() -> dict[str, Any]:
        amplicon = amplify.pcr(primers[0], primers[1], tmpl, limit=min_anneal)
        return {
            "figure": amplicon.figure(),
            "product": str(amplicon.seq.watson),
            "product_length": len(amplicon),
            "forward_position": int(amplicon.forward_primer.position),
            "reverse_position": int(amplicon.reverse_primer.position),
        }

    found = _guard("PCR", run)
    return _done("pcr_alignment_figure", limit=min_anneal, **found)


def pydna_assembly_figure(
    fragments: list[Any],
    *,
    limit: int = 25,
    figure: str = "compact",
    product_index: int = 0,
    circular: bool = False,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assemble fragments by their shared end sequences and draw the junctions.

    Uses pydna ``assembly.Assembly`` and its contigs. ``figure`` is ``compact`` (``Contig.figure``,
    fragment names with the overlap lengths), ``detailed`` (``Contig.detailed_figure``, the aligned
    overlap sequences) or ``matplotlib`` (``Contig.figure_mpl``, written to ``output_path`` as
    a PNG). ``limit`` is the minimum overlap length. ``circular`` chooses circular instead of linear
    products; ``product_index`` selects one of several products (longest first).
    """

    frags = _records(fragments, "fragments", 2)
    kind = _choice(figure, _ASSEMBLY_FIGURES, "figure")
    overlap = _integer(limit, "limit", 4, 10_000)
    index = _integer(product_index, "product_index", 0, 10_000)
    use_circular = _flag(circular, "circular")
    target = None
    if kind == "matplotlib":
        if not isinstance(output_path, str) or not output_path.strip():
            raise ConfigurationError("output_path (a .png path) is required for the matplotlib figure.")
        target = Path(output_path).expanduser().resolve()
        if target.suffix.lower() != ".png":
            raise ConfigurationError("output_path must end with .png.")

    def run() -> dict[str, Any]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            assembly = require("pydna.assembly", _EXTRA).Assembly(frags, limit=overlap)
            contigs = sorted(
                assembly.assemble_circular() if use_circular else assembly.assemble_linear(),
                key=len, reverse=True,
            )
        if not contigs:
            raise ValueError("no assembly product was found for these fragments")
        if index >= len(contigs):
            raise ValueError(f"product_index {index} is out of range for {len(contigs)} products")
        contig = contigs[index]
        out: dict[str, Any] = {
            "product_count": len(contigs),
            "product": str(contig.seq.watson),
            "product_length": len(contig),
            "circular": bool(contig.circular),
        }
        if kind == "compact":
            out["figure"] = contig.figure()
        elif kind == "detailed":
            out["figure"] = contig.detailed_figure()
        else:
            plt = require("matplotlib.pyplot", "bio-plots")
            fig = contig.figure_mpl()
            assert target is not None
            target.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(target, format="png", dpi=120)
            plt.close(fig)
            out["image_path"] = str(target)
        return out

    found = _guard("assembly figure", run)
    return _done("assembly_figure", limit=overlap, figure_type=kind, **found)


# OpenCloning strategies and cloning history.


def _strategy_module() -> Any:
    return require("pydna.opencloning_models", _EXTRA)


def pydna_build_cloning_strategy(
    fragments: list[Any],
    *,
    operation: str = "none",
    primers: list[Any] | None = None,
    enzymes: list[str] | None = None,
    limit: int = 25,
    product_index: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a cloning step with pydna and export the product and its history as OpenCloning JSON.

    ``operation`` is ``none`` (export the fragments as they are), ``pcr`` (``fragments[0]`` is the
    template; ``primers`` are the forward and reverse primer), ``gibson`` (Gibson assembly of the
    fragments with minimum overlap ``limit``; for ``pcr`` the primer anneal limit is ``limit`` capped at 14) or ``restriction_ligation`` (``enzymes`` required). The
    export is pydna ``CloningStrategy.from_dseqrecords`` followed by ``model_dump_json``; the JSON can
    be read back with ``pydna_cloning_strategy_to_records``.
    """

    chosen = _choice(operation, _STRATEGY_OPERATIONS, "operation")
    frags = _records(fragments, "fragments", 1)
    overlap = _integer(limit, "limit", 4, 10_000)
    index = _integer(product_index, "product_index", 0, 10_000)
    models = _strategy_module()

    def run() -> dict[str, Any]:
        if chosen == "none":
            products = frags
        elif chosen == "pcr":
            pair = _primers(primers, "primers")
            if len(pair) != 2:
                raise ValueError("primers must hold the forward and the reverse primer")
            products = require("pydna.assembly2", _EXTRA).pcr_assembly(
                frags[0], pair[0], pair[1], limit=min(overlap, 14)
            )
        elif chosen == "gibson":
            products = require("pydna.assembly2", _EXTRA).gibson_assembly(frags, limit=overlap)
        else:
            if enzymes is None:
                raise ValueError("enzymes are required for restriction_ligation")
            products = require("pydna.assembly2", _EXTRA).restriction_ligation_assembly(
                frags, _enzymes(enzymes, "enzymes")
            )
        if not products:
            raise ValueError("the operation gave no product")
        if chosen != "none":
            if index >= len(products):
                raise ValueError(f"product_index {index} is out of range for {len(products)} products")
            products = [products[index]]
        for number, product in enumerate(products):
            product.name = f"product{number + 1}"
        strategy = models.CloningStrategy.from_dseqrecords(products)
        return {
            "strategy_json": strategy.model_dump_json(),
            "product_count": len(products),
            "product_names": [str(p.name) for p in products],
            "history": [p.history() for p in products],
        }

    found = _guard("cloning strategy export", run)
    return _done("build_cloning_strategy", cloning_operation=chosen, **found)


def pydna_cloning_strategy_to_records(
    strategy_json: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Read an OpenCloning cloning-strategy JSON text and return the sequences it contains.

    Uses pydna ``opencloning_models.CloningStrategy.to_dseqrecords`` (which also reads the
    ``TextFileSequence`` entries through ``read_dseqrecord_from_text_file_sequence``). Each record
    reports its name, sequence, topology and the cloning history as text.
    """

    text = _text(strategy_json, "strategy_json", maximum=50_000_000)
    models = _strategy_module()

    def run() -> list[dict[str, Any]]:
        records = models.CloningStrategy.model_validate_json(text).to_dseqrecords()
        return [
            {
                "name": str(r.name),
                "sequence": str(r.seq.watson),
                "length": len(r),
                "circular": bool(r.circular),
                "history": r.history(),
            }
            for r in records
        ]

    rows = _guard("cloning strategy import", run)
    return _done("cloning_strategy_to_records", record_count=len(rows), records=rows)


def pydna_cloning_history(
    strategy_json: str, *, action: str = "display", verbose: bool = True
) -> dict[str, Any]:
    """Show, validate or normalise the cloning history of the sequences in an OpenCloning strategy.

    ``action`` is ``display`` (the history tree as text, pydna ``Dseqrecord.history``),
    ``validate`` (replay every recorded operation and compare the result with the stored sequence,
    ``Dseqrecord.validate_history``) or ``normalize`` (replay the operations and return the rebuilt
    sequences, ``Dseqrecord.normalize_history``). Records without history are reported as such.
    """

    mode = _choice(action, ("display", "validate", "normalize"), "action")
    text = _text(strategy_json, "strategy_json", maximum=50_000_000)
    models = _strategy_module()

    def run() -> list[dict[str, Any]]:
        records = models.CloningStrategy.model_validate_json(text).to_dseqrecords()
        rows = []
        for r in records:
            row: dict[str, Any] = {"name": str(r.name), "has_history": r.source is not None
                                   and bool(r.history())}
            if mode == "display":
                row["history"] = r.history()
            elif mode == "validate":
                try:
                    r.validate_history()
                    row["valid"] = True
                except _ERRORS + (Exception,) as exc:  # a failed replay is a result, not an error
                    row["valid"] = False
                    row["error"] = str(exc)[:500]
            else:
                new = r.normalize_history()
                row["sequence"] = str(new.seq.watson)
                row["unchanged"] = str(new.seq.watson) == str(r.seq.watson)
            rows.append(row)
        return rows

    rows = _guard("cloning history", run)
    return _done("cloning_history", action=mode, record_count=len(rows), records=rows)


# Records.


def pydna_dna_mass_from_amount(
    dna: str | dict[str, Any], *, amount_mol: float = 5e-14, verbose: bool = True
) -> dict[str, Any]:
    """Mass in grams of an amount of a double-stranded molecule (pydna ``Dseqrecord.m``).

    The mass is the molecular weight of the molecule (``Dseq.mw``, g/mol) times ``amount_mol``
    (mol). pydna's default amount is 5e-14 mol (0.05 pmol).
    """

    molecule = _dseq(dna, "dna")
    moles = _number(amount_mol, "amount_mol")
    if moles < 0:
        raise ConfigurationError("amount_mol must not be negative.")
    weight = float(_guard("mw", lambda: molecule.mw()))
    return _done("dna_mass_from_amount", amount_mol=moles, molecular_weight=weight, mass_g=weight * moles)


def pydna_find_peptide(
    dna: str, peptide: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Find a peptide in the three forward frames of a DNA sequence (pydna ``find_aminoacids``).

    Returns the DNA slice (0-based start, end) of the first match, or None. ``find_aa`` is the
    alias. The sequence is searched in its three forward frames (the reverse strand is not).
    """

    record = _record(dna, "dna")
    letters = _letters(peptide, "peptide", maximum=10_000)
    match = _guard("find_aminoacids", lambda: record.find_aminoacids(letters))
    start = None if match is None else int(match.start)
    end = None if match is None else int(match.stop)
    return _done("find_peptide", peptide=letters.upper(), found=match is not None, start=start,
                 end=end, dna=None if match is None else str(record.seq.watson)[start:end])


def pydna_map_sanger_traces(
    sequence: str | dict[str, Any],
    trace_glob: str,
    *,
    limit: int = 25,
    target_start: int | None = None,
    target_end: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map ABI Sanger trace files onto a sequence (pydna ``Dseqrecord.map_trace_files``).

    ``trace_glob`` is a file pattern such as ``/data/reads/*.ab1``. A trace maps when it shares a
    stretch of at least ``limit`` bases with the sequence (either strand). With ``target_start`` and
    ``target_end`` only traces that cover that region count as matching. The result lists the matching
    files with the mapped coordinates.
    """

    record = _record(sequence, "sequence")
    pattern = _text(trace_glob, "trace_glob", maximum=1000)
    minimum = _integer(limit, "limit", 4, 10_000)
    if (target_start is None) != (target_end is None):
        raise ConfigurationError("target_start and target_end must be given together.")
    if target_start is not None:
        low = _integer(target_start, "target_start")
        high = _integer(target_end, "target_end")
        if high <= low:
            raise ConfigurationError("target_end must be greater than target_start.")
        record.map_target = slice(low, high)

    def run() -> dict[str, Any]:
        files = record.map_trace_files(pattern, limit=minimum)
        mapped = [
            {"file": str(f.qualifiers["label"][0]),
             "start": int(f.location.start), "end": int(f.location.end)}
            for f in record.features if f.type == "trace"
        ]
        out: dict[str, Any] = {"matching_files": [str(f) for f in files], "mapped": mapped}
        if getattr(record, "matching_reads", None) is not None:
            out["not_matching_files"] = [
                str(t.annotations["filename"]) for t in record.not_matching_reads
            ]
        return out

    found = _guard("map_trace_files", run)
    return _done("map_sanger_traces", limit=minimum, **found)


def pydna_sync_circular_origin(
    sequence: str, reference: str, *, limit: int = 25, verbose: bool = True
) -> dict[str, Any]:
    """Rotate (and if needed reverse-complement) a circular sequence to start where a reference starts.

    This is pydna ``Dseqrecord.synced``. The two sequences must share a stretch of at least
    ``limit`` bases (``limit`` is lowered for short sequences, as pydna does).
    """

    record = require("pydna.dseqrecord", _EXTRA).Dseqrecord(_letters(sequence, "sequence"), circular=True)
    ref = _letters(reference, "reference")
    minimum = _integer(limit, "limit", 1, 10_000)
    new = _guard("synced", lambda: record.synced(ref, limit=minimum))
    return _done("sync_circular_origin", limit=minimum, input=str(record.seq.watson),
                 result=str(new.seq.watson), circular=bool(new.circular))


def pydna_orfs_to_features(
    sequence: str | dict[str, Any], *, min_size: int = 300, verbose: bool = True
) -> dict[str, Any]:
    """Annotate the open reading frames of both strands as CDS features with their translation.

    This is pydna ``Dseqrecord.orfs_to_features``. ``min_size`` is the minimum ORF length in bases
    (pydna default 300). Each feature has the start, end and strand, the translation and the SEGUID
    checksums of the DNA and the protein.
    """

    record = _record(sequence, "sequence")
    minimum = _integer(min_size, "min_size", 3, _MAX_BP)
    features = _guard("orfs_to_features", lambda: record.orfs_to_features(minimum))
    rows = [
        {
            "start": int(f.location.start),
            "end": int(f.location.end),
            "strand": int(f.location.strand),
            "length": int(f.location.end - f.location.start),
            "translation": f.qualifiers["translation"],
            "checksum": list(f.qualifiers["checksum"]),
        }
        for f in features
    ]
    return _done("orfs_to_features", min_size=minimum, feature_count=len(rows), features=rows)


def pydna_clean_genbank_text(genbank_text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Repair malformed GenBank text so that it can be parsed (pydna ``genbankfixer.gbtext_clean``).

    Returns the repaired GenBank text and the parsed JSON representation (name, size, topology,
    sequence, features) that pydna builds on the way.
    """

    text = _text(genbank_text, "genbank_text", maximum=50_000_000)
    fixer = require("pydna.genbankfixer", _EXTRA)
    cleaned = _guard("gbtext_clean", lambda: fixer.gbtext_clean(text))
    jseq = cleaned.jseq
    return _done("clean_genbank_text", genbank_text=cleaned.gbtext, name=jseq.get("name"),
                 length=len(jseq.get("sequence", "")), topology=jseq.get("topology"),
                 feature_count=len(jseq.get("features", [])))


# Primer screening.


def _screen_module() -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return require("pydna.primer_screen", _EXTRA)


def _screen(call: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return _guard("primer screen", call)


def _amplicon_rows(rows: Sequence[Any], primers: list[Any]) -> list[dict[str, Any]]:
    return [
        {
            "forward_primer": str(primers[int(r.fp)].name),
            "reverse_primer": str(primers[int(r.rp)].name),
            "forward_position": int(r.fposition),
            "reverse_position": int(r.rposition),
            "size": int(r.size),
        }
        for r in rows
    ]


def pydna_primer_pairs_by_size(
    template: str | dict[str, Any],
    primers: list[Any],
    *,
    short: int = 500,
    long: int = 2000,
    limit: int = 16,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find primer pairs from a primer list whose PCR product size lies in a range.

    This is pydna ``primer_screen.primer_pairs``. Product size includes the primers; only primers
    that bind in one place each are used. ``limit`` is the length of the 3' end used to find the
    binding sites. A circular template also yields products across the origin.
    """

    record = _record(template, "template")
    plist = _primers(primers, "primers")
    low, high = _integer(short, "short"), _integer(long, "long")
    seed = _integer(limit, "limit", 4, 100)
    screen = _screen_module()
    rows = _screen(lambda: screen.primer_pairs(record, plist, low, high, seed))
    return _done("primer_pairs_by_size", short=low, long=high, pair_count=len(rows),
                 pairs=_amplicon_rows(rows, plist))


def pydna_flanking_primer_pairs(
    template: str | dict[str, Any],
    primers: list[Any],
    *,
    target_start: int,
    target_end: int,
    limit: int = 16,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find primer pairs that flank a target region of the template.

    This is pydna ``primer_screen.flanking_primer_pairs``: the forward primer binds at or before
    ``target_start`` and the reverse primer at or after ``target_end`` (0-based template
    coordinates).
    """

    record = _record(template, "template")
    plist = _primers(primers, "primers")
    start, end = _integer(target_start, "target_start"), _integer(target_end, "target_end")
    seed = _integer(limit, "limit", 4, 100)
    screen = _screen_module()
    rows = _screen(lambda: screen.flanking_primer_pairs(record, plist, (start, end), seed))
    return _done("flanking_primer_pairs", target=[start, end], pair_count=len(rows),
                 pairs=_amplicon_rows(rows, plist))


def pydna_diagnostic_primer_sets(
    templates: list[Any],
    primers: list[Any],
    *,
    primers_per_set: int = 2,
    min_relative_size_difference: float = 0.2,
    short: int = 500,
    long: int = 1500,
    limit: int = 16,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select primer pairs or triplets whose PCR products have different sizes on each template.

    ``primers_per_set`` is 2 (pydna ``primer_screen.diff_primer_pairs``, a primer pair that gives a
    distinct product size from every template) or 3 (``primer_screen.diff_primer_triplets``, a shared
    primer plus one primer specific to each template). Such sets verify clones or genetic
    modifications by product size. Needs at least two templates. Two product sizes count as distinct
    when they differ by at least ``min_relative_size_difference`` (default 0.2, pydna's rule of 20 %
    of the larger product) of the larger one.
    """

    records = _records(templates, "templates", 2)
    plist = _primers(primers, "primers")
    size = _integer(primers_per_set, "primers_per_set", 2, 3)
    if size not in (2, 3):
        raise ConfigurationError("primers_per_set must be 2 or 3.")
    low, high = _integer(short, "short"), _integer(long, "long")
    seed = _integer(limit, "limit", 4, 100)
    fraction = _number(min_relative_size_difference, "min_relative_size_difference")
    if not 0 <= fraction < 1:
        raise ConfigurationError("min_relative_size_difference must be at least 0 and below 1.")
    screen = _screen_module()

    def distinct(a: int, b: int) -> bool:
        return abs(a - b) >= fraction * max(a, b)

    def convert(item: Any) -> dict[str, Any]:
        record, fp, rp, length = item
        return {"template": str(record.name), "forward_primer": str(plist[int(fp)].name),
                "reverse_primer": str(plist[int(rp)].name), "product_size": int(length)}

    if size == 2:
        sets = _screen(lambda: screen.diff_primer_pairs(
            records, plist, low, high, seed, callback=distinct))
        rows = [[convert(item) for item in group] for group in sets]
    else:
        sets = _screen(lambda: screen.diff_primer_triplets(
            records, plist, seed, low, high, callback=distinct))
        rows = [[convert(item) for item in group] for group in sets]
    return _done("diagnostic_primer_sets", primers_per_set=size, set_count=len(rows), sets=rows)


# Recombinase sites.


def pydna_find_recombinase_sites(
    sequence: str | dict[str, Any],
    site1: str,
    site2: str,
    *,
    site1_name: str = "site1",
    site2_name: str = "site2",
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the two recognition sites of a site-specific recombinase on both strands.

    This is pydna ``recombinase.Recombinase.find`` and ``Recombinase.annotate``
    (and the same methods of ``RecombinaseCollection``). The homology core of each site is written in
    lowercase, as in ``ATGCCCTAAaaTT`` / ``AAaaTTTTTTTCCCT``, and the two cores must match. The result
    lists the 0-based locations per site and the annotated GenBank features.
    """

    record = _record(sequence, "sequence")
    first = _letters(site1, "site1", maximum=200)
    second = _letters(site2, "site2", maximum=200)
    names = (_text(site1_name, "site1_name", maximum=60), _text(site2_name, "site2_name", maximum=60))
    module = require("pydna.recombinase", _EXTRA)

    def run() -> dict[str, Any]:
        recombinase = module.Recombinase(first, second, names[0], names[1])
        found = recombinase.find(record)
        annotated = recombinase.annotate(record)
        return {
            "sites": {
                name: [{"start": int(loc.start), "end": int(loc.end), "strand": int(loc.strand)}
                       for loc in locations]
                for name, locations in found.items()
            },
            "annotated_features": [
                {"type": str(f.type), "start": int(f.location.start), "end": int(f.location.end),
                 "strand": int(f.location.strand)}
                for f in annotated.features
            ],
        }

    found = _guard("Recombinase.find", run)
    return _done("find_recombinase_sites", **found)


# Codon usage values and sequence utilities.


def _codon_tables(organism: str) -> Any:
    codon = require("pydna.codon", _EXTRA)
    key = _choice(organism, tuple(codon.weights), "organism")
    return key, codon


def pydna_terminal_codon_usage(
    cds: str, *, organism: str = "sce", codon: str = "both", verbose: bool = True
) -> dict[str, Any]:
    """Look up the usage value of the first and last codon of a coding sequence.

    ``codon`` is ``start`` (pydna ``Seq.startcodon``), ``stop`` (``Seq.stopcodon``) or ``both``. The
    value is the entry of pydna's start-codon or stop-codon usage table for ``organism`` (only
    ``sce``, Saccharomyces cerevisiae, is available) and is None when the codon is not in the table.
    """

    text = _letters(cds, "cds").upper()
    which = _choice(codon, ("start", "stop", "both"), "codon")
    key, tables = _codon_tables(organism)
    result: dict[str, Any] = {"organism": key}
    if which in ("start", "both"):
        result["start_codon"] = text[:3]
        result["start_value"] = tables.start[key].get(text[:3])
    if which in ("stop", "both"):
        result["stop_codon"] = text[-3:]
        result["stop_value"] = tables.stop[key].get(text[-3:])
    return _done("terminal_codon_usage", codon=which, **result)


def pydna_cds_expression_report(
    cds: str, *, organism: str = "sce", verbose: bool = True
) -> dict[str, Any]:
    """Summarise a coding sequence for expression: CAI, GC, start and stop codon, N-end, rare codons.

    This is the report of pydna ``Seq.express``: length in codons, codon adaptation index (geometric
    mean of the organism's codon weights, as ``pydna_codon_adaptation_index``), GC fraction, start and
    stop codon values, the N-end rule value of the second amino acid, the count of each rare codon
    and the fraction of rare codons. ``organism`` is ``sce``.
    """

    from dnakit.molecular_lab.pydna_gap_functions import pydna_codon_adaptation_index

    text = _letters(cds, "cds").upper()
    if len(text) < 6 or len(text) % 3:
        raise ConfigurationError("cds must be a whole number of codons (at least two).")
    key, tables = _codon_tables(organism)
    biopython = require("Bio.SeqUtils", "bio-core")
    seq_module = require("Bio.Seq", "bio-core")
    codons = [text[i : i + 3] for i in range(0, len(text), 3)]
    rare = list(tables.rare_codons[key])
    counts = {codon: codons.count(codon) for codon in rare}
    second = biopython.seq3(str(seq_module.Seq(codons[1]).translate()))
    cai = pydna_codon_adaptation_index(text, weights=tables.weights[key], verbose=False)["cai"]
    return _done(
        "cds_expression_report",
        organism=key,
        codon_count=len(codons),
        label=f"{text[:3]}...{text[-3:]}",
        cai=cai,
        gc=round(biopython.gc_fraction(text), 3),
        start_value=tables.start[key].get(text[:3]),
        stop_value=tables.stop[key].get(text[-3:]),
        n_end=tables.n_end[key].get(second),
        rare_codon_counts=counts,
        rare_fraction=round(sum(counts.values()) / len(codons), 3),
    )


def pydna_is_open_reading_frame(
    sequence: str, *, table: int = 1, verbose: bool = True
) -> dict[str, Any]:
    """Check whether a sequence is an open reading frame on its 5'->3' strand (pydna ``isorf``).

    True when it starts with a start codon, ends with a stop codon, is a whole number of codons and
    has no internal stop codon, under the NCBI translation ``table`` (default 1, the standard code).
    """

    text = _letters(sequence, "sequence")
    number = _integer(table, "table", 1, 50)
    seqrecord = require("pydna.seqrecord", _EXTRA)
    answer = _guard("isorf", lambda: seqrecord.SeqRecord(text).isorf(number))
    return _done("is_open_reading_frame", table=number, is_orf=bool(answer), length=len(text))


def pydna_seguid_stamp(
    sequence: str | dict[str, Any],
    *,
    timestamp: str | None = None,
    tool: str = "pydna",
    comment: str = "",
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the SEGUID checksum stamp of a record and write it into the record's COMMENT.

    This is pydna ``SeqRecord.stamp``. The stamp is ``cdseguid=...`` for a circular molecule and
    ``ldseguid=...`` for a linear one; the comment reads ``<tool> <stamp> <timestamp> <comment>``.
    ``timestamp`` defaults to the current time (ISO format); pass a fixed text for reproducible
    output.
    """

    record = _record(sequence, "sequence")
    stamp_tool = _text(tool, "tool", maximum=60)
    note = comment if isinstance(comment, str) else ""
    moment = None if timestamp is None else _text(timestamp, "timestamp", maximum=60)
    stamp = _guard("stamp", lambda: record.stamp(
        now=(lambda: moment) if moment is not None else _now,
        tool=stamp_tool, comment=note,
    ))
    return _done("seguid_stamp", stamp=str(stamp), comment=record.annotations.get("comment", ""),
                 circular=bool(record.circular))


def pydna_iupac_site_to_regex(site: str, *, verbose: bool = True) -> dict[str, Any]:
    """Convert a degenerate (IUPAC) DNA site into a regular expression (pydna ``compute_regex_site``).

    Ambiguity codes become character classes (``R`` becomes ``[AG]``) and the pattern is case
    insensitive. The result also reports whether the pattern is the same for the reverse
    complement.
    """

    text = _letters(site, "site", maximum=1000)
    module = require("pydna.sequence_regex", _EXTRA)
    pattern = _guard("compute_regex_site", lambda: module.compute_regex_site(text))
    return _done("iupac_site_to_regex", site=text, regex=pattern)


__all__ = [
    "pydna_assembly_figure",
    "pydna_build_cloning_strategy",
    "pydna_cast_to_double_strand",
    "pydna_cds_expression_report",
    "pydna_clean_genbank_text",
    "pydna_cloning_history",
    "pydna_cloning_strategy_to_records",
    "pydna_diagnostic_primer_sets",
    "pydna_dna_mass_from_amount",
    "pydna_exonuclease_five_prime",
    "pydna_exonuclease_three_prime",
    "pydna_fill_in_ends",
    "pydna_find_peptide",
    "pydna_find_recombinase_sites",
    "pydna_flanking_primer_pairs",
    "pydna_iupac_site_to_regex",
    "pydna_is_open_reading_frame",
    "pydna_map_sanger_traces",
    "pydna_melt_short_strands",
    "pydna_mung_bean_blunt",
    "pydna_orfs_to_features",
    "pydna_pcr_alignment_figure",
    "pydna_primer_pairs_by_size",
    "pydna_restriction_cutters",
    "pydna_restriction_ligation_assembly",
    "pydna_seguid_stamp",
    "pydna_shed_strand_segments",
    "pydna_sync_circular_origin",
    "pydna_t4_polymerase_polish",
    "pydna_terminal_codon_usage",
    "pydna_terminal_transferase_tail",
    "pydna_user_enzyme_treatment",
]

_register_public_verbose(globals(), tuple(__all__))
