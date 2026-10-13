"""pydna functions that the family check found to have no existing dnakit wrapper (H2).

Each function imports pydna only when it runs. A missing package raises
``BackendUnavailableError`` with the ``dnakit[bio-assembly]`` install hint. Positions are
0-based. Every function returns a JSON-safe dictionary with ``status``, ``operation``, the
result fields and ``provenance``. Inputs are DNA strings of A, C, G and T; they are upper-cased.
"""

from __future__ import annotations

import importlib.metadata
import math
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-assembly"
_MAX_PAIR = 20_000
_MAX_PRIMER = 200
_MAX_PRODUCTS = 200
_MAX_MATCHES = 1000
_ANNEALING_POLYMERASES = ("taq", "dbd")


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _positive(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a positive number.")
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ConfigurationError(f"{name} must be a positive number.")
    return number


def _dna(value: object, name: str, *, maximum: int = _MAX_PAIR) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _primer(value: object, name: str) -> str:
    text = _dna(value, name, maximum=_MAX_PRIMER)
    if len(text) < 6:
        raise ConfigurationError(f"{name} must be at least 6 bases long.")
    return text


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("pydna")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "pydna", "version": installed, **notes}


def _run(action: Any, message: str) -> Any:
    try:
        return action()
    except (AssertionError, AttributeError, IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{message}: {exc}") from exc


def find_common_substrings_with_topology(
    sequence_x: str,
    sequence_y: str,
    *,
    min_length: int = 25,
    circular_x: bool = False,
    circular_y: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find common substrings of two DNA molecules, wrapping around circular ones.

    A circular molecule is searched across its origin, so a match may start near the end
    and continue at the beginning; the match is never longer than the shorter molecule.
    For two linear molecules the matches equal those of ``find_common_substrings``. Each
    match is ``(x_start, y_start, length)`` with 0-based starts. Matches are sorted and
    capped at 1000 rows; ``match_count`` gives the total.
    """

    assembly = require("pydna.assembly2", _EXTRA)
    record_class = require("pydna.dseqrecord", _EXTRA).Dseqrecord
    first = _dna(sequence_x, "sequence_x")
    second = _dna(sequence_y, "sequence_y")
    limit = _integer(min_length, "min_length", minimum=1, maximum=_MAX_PAIR)
    circular_first = _flag(circular_x, "circular_x")
    circular_second = _flag(circular_y, "circular_y")
    raw = _run(
        lambda: assembly.common_sub_strings(
            record_class(first, circular=circular_first),
            record_class(second, circular=circular_second),
            limit,
        ),
        "pydna could not compare the sequences",
    )
    rows = sorted(
        ({"x_start": int(a), "y_start": int(b), "length": int(c)} for a, b, c in raw),
        key=lambda row: (row["x_start"], row["y_start"], row["length"]),
    )
    return {
        "status": "completed",
        "operation": "find_common_substrings_with_topology",
        "min_length": limit,
        "circular_x": circular_first,
        "circular_y": circular_second,
        "match_count": len(rows),
        "matches": rows[:_MAX_MATCHES],
        "provenance": _provenance(method="pydna.assembly2.common_sub_strings"),
    }


def simulate_pcr_assembly(
    template: str,
    forward_primer: str,
    reverse_primer: str,
    *,
    min_annealing: int = 14,
    mismatches: int = 0,
    template_circular: bool = False,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate every PCR product of a primer pair, allowing primer mismatches.

    This is pydna's assembly-based PCR. Each primer must anneal over at least
    ``min_annealing`` matching bases; up to ``mismatches`` mismatched bases are tolerated
    in the annealing region, which ``simulate_pcr_amplicon`` does not allow. A circular
    template can give a product across the origin, and more than one site gives more than
    one product. Products are sorted by length and sequence; at most ``max_products`` are
    listed and ``product_count`` gives the total. No product gives an empty list.
    """

    assembly = require("pydna.assembly2", _EXTRA)
    record_class = require("pydna.dseqrecord", _EXTRA).Dseqrecord
    primer_class = require("pydna.primer", _EXTRA).Primer
    text = _dna(template, "template")
    forward = _dna(forward_primer, "forward_primer", maximum=_MAX_PRIMER)
    reverse = _dna(reverse_primer, "reverse_primer", maximum=_MAX_PRIMER)
    limit = _integer(min_annealing, "min_annealing", minimum=8, maximum=60)
    allowed = _integer(mismatches, "mismatches", minimum=0, maximum=5)
    circular = _flag(template_circular, "template_circular")
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    products = _run(
        lambda: assembly.pcr_assembly(
            record_class(text, circular=circular),
            primer_class(forward),
            primer_class(reverse),
            limit=limit,
            mismatches=allowed,
        ),
        "pydna could not simulate the PCR",
    )
    rows = sorted(
        (
            {
                "sequence": str(item.seq).upper(),
                "length": len(item),
                "circular": bool(item.circular),
            }
            for item in products
        ),
        key=lambda row: (row["length"], row["sequence"], row["circular"]),
    )
    return {
        "status": "completed",
        "operation": "simulate_pcr_assembly",
        "template_length": len(text),
        "template_circular": circular,
        "min_annealing": limit,
        "mismatches": allowed,
        "product_count": len(rows),
        "products": rows[:listed],
        "provenance": _provenance(method="pydna.assembly2.pcr_assembly"),
    }


def primer_melting_temperature_dbd(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Compute a primer melting temperature with pydna's rule for DNA-binding-domain polymerases.

    This is the nearest-neighbour Tm with the SantaLucia 1998 table, 50 mM Na+, 1.5 mM Mg2+
    and 250 nM of each strand that pydna uses for Phusion-like and Pfu-Sso7d enzymes. It
    differs from ``primer_melting_temperature``, which uses the Taq buffer. The value is
    rounded to three decimals.
    """

    module = require("pydna.tm", _EXTRA)
    text = _primer(sequence, "sequence")
    value = _run(lambda: float(module.tm_dbd(text)), "pydna could not compute Tm")
    return {
        "status": "completed",
        "operation": "primer_melting_temperature_dbd",
        "sequence": text,
        "length": len(text),
        "tm_celsius": round(value, 3),
        "provenance": _provenance(method="pydna.tm.tm_dbd", rounding="3 decimals"),
    }


def primer_melting_temperature_breslauer(
    sequence: str,
    *,
    primer_nm: float = 500.0,
    salt_mm: float = 50.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a primer melting temperature with the Breslauer parameters in pydna.

    This is pydna's ``tmbresluc`` formula for polymerases with a DNA binding domain. It uses
    Breslauer nearest-neighbour enthalpy and entropy values, the primer concentration
    ``primer_nm`` (nM) and the monovalent cation concentration ``salt_mm`` (mM). The value
    is rounded to three decimals.
    """

    module = require("pydna.tm", _EXTRA)
    text = _primer(sequence, "sequence")
    primer_conc = _positive(primer_nm, "primer_nm")
    salt_conc = _positive(salt_mm, "salt_mm")
    value = _run(
        lambda: float(module.tmbresluc(text, primerc=primer_conc, saltc=salt_conc)),
        "pydna could not compute Tm",
    )
    return {
        "status": "completed",
        "operation": "primer_melting_temperature_breslauer",
        "sequence": text,
        "length": len(text),
        "primer_nm": primer_conc,
        "salt_mm": salt_conc,
        "tm_celsius": round(value, 3),
        "provenance": _provenance(method="pydna.tm.tmbresluc", rounding="3 decimals"),
    }


def amplicon_melting_temperature(
    sequence: str,
    *,
    potassium_m: float = 0.05,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the melting temperature of a PCR product from its length and GC content.

    This is the Rychlik, Spencer and Rhoads (1990) formula that pydna uses for the amplicon:
    ``81.5 + 0.41 * GC% + 16.6 * log10(K) - 675 / length``, where ``potassium_m`` is the
    monovalent cation concentration in mol/L. The value is rounded to three decimals.
    """

    module = require("pydna.tm", _EXTRA)
    text = _dna(sequence, "sequence")
    if len(text) < 20:
        raise ConfigurationError("sequence must be at least 20 bases long.")
    cation = _positive(potassium_m, "potassium_m")
    value = _run(lambda: float(module.tm_product(text, K=cation)), "pydna could not compute Tm")
    gc = sum(base in "GC" for base in text) / len(text)
    return {
        "status": "completed",
        "operation": "amplicon_melting_temperature",
        "length": len(text),
        "gc_fraction": round(gc, 6),
        "potassium_m": cation,
        "tm_celsius": round(value, 3),
        "provenance": _provenance(method="pydna.tm.tm_product", rounding="3 decimals"),
    }


def pcr_annealing_temperature(
    forward_primer: str,
    reverse_primer: str,
    product: str,
    *,
    polymerase: str = "taq",
    verbose: bool = True,
) -> dict[str, Any]:
    """Suggest a PCR annealing temperature for a primer pair and its product with pydna.

    ``polymerase`` is ``"taq"`` (default; the Rychlik rule, ``0.3 * lower primer Tm +
    0.7 * product Tm - 14.9``) or ``"dbd"`` (the lower primer Tm plus 3 degrees, with the
    DNA-binding-domain Tm rule). ``forward_primer`` and ``reverse_primer`` are the primer
    annealing parts and ``product`` is the amplicon sequence. Values are rounded to three
    decimals.
    """

    module = require("pydna.tm", _EXTRA)
    if polymerase not in _ANNEALING_POLYMERASES:
        allowed = ", ".join(_ANNEALING_POLYMERASES)
        raise ConfigurationError(f"polymerase must be one of: {allowed}.")
    forward = _primer(forward_primer, "forward_primer")
    reverse = _primer(reverse_primer, "reverse_primer")
    amplicon = _dna(product, "product")
    if len(amplicon) < 20:
        raise ConfigurationError("product must be at least 20 bases long.")
    if polymerase == "taq":
        function, tm_function = module.ta_default, module.tm_default
        function_name = "ta_default"
    else:
        function, tm_function = module.ta_dbd, module.tm_dbd
        function_name = "ta_dbd"
    value = _run(lambda: float(function(forward, reverse, amplicon)), "pydna could not compute Ta")
    forward_tm = _run(lambda: float(tm_function(forward)), "pydna could not compute Tm")
    reverse_tm = _run(lambda: float(tm_function(reverse)), "pydna could not compute Tm")
    return {
        "status": "completed",
        "operation": "pcr_annealing_temperature",
        "polymerase": polymerase,
        "forward_tm_celsius": round(forward_tm, 3),
        "reverse_tm_celsius": round(reverse_tm, 3),
        "product_length": len(amplicon),
        "annealing_celsius": round(value, 3),
        "provenance": _provenance(method=f"pydna.tm.{function_name}", rounding="3 decimals"),
    }


__all__ = [
    "amplicon_melting_temperature",
    "find_common_substrings_with_topology",
    "pcr_annealing_temperature",
    "primer_melting_temperature_breslauer",
    "primer_melting_temperature_dbd",
    "simulate_pcr_assembly",
]

_register_public_verbose(globals(), tuple(__all__))
