"""primer3-py end-stability, primer-pair design and sanitizing helpers (optional).

primer3-py is GPL-2.0. It is imported only inside each function, so dnakit does not
depend on it at import time, and no primer3 code is copied here. The Tm and
hairpin/dimer calculations are exposed by :mod:`dnakit.molecular_lab.backends` as
``calculate_primer3_tm``, ``calculate_primer3_hairpin`` and ``calculate_primer3_dimer``.
"""

from __future__ import annotations

import importlib.metadata
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require, sequence_text
from dnakit.exceptions import ConfigurationError

_MAX_OLIGO = 200
_MAX_TEMPLATE = 10_000
_MAX_SANITIZE = 100_000


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version("primer3-py")
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": "primer3-py", "version": version, "license": "GPL-2.0", **notes}


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return number


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _acgtn(value: object, name: str, *, maximum: int, minimum: int = 1) -> str:
    text = sequence_text(value)
    if not minimum <= len(text) <= maximum:
        raise ConfigurationError(f"{name} must have {minimum} to {maximum} bases.")
    if set(text) - set("ACGTN"):
        raise ConfigurationError(f"{name} must contain only A, C, G, T and N.")
    return text


def _float_or_none(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def calculate_primer3_end_stability(
    primer: str,
    template: str,
    *,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """3' end stability of ``primer`` annealed to ``template`` (primer3-py ``calc_end_stability``).

    Concentrations are mM for ``mv_conc``, ``dv_conc`` and ``dntp_conc``, and nM
    for ``dna_conc``. ``dg``, ``dh`` and ``ds`` are returned for the end-stability
    structure. When no structure is found, the thermodynamic values are ``None``.
    """

    primer3 = require("primer3", "bio-primer3")
    primer_text = _acgtn(primer, "primer", maximum=_MAX_OLIGO)
    template_text = _acgtn(template, "template", maximum=_MAX_OLIGO)
    mv = _number(mv_conc, "mv_conc", minimum=0.0, maximum=1000.0)
    dv = _number(dv_conc, "dv_conc", minimum=0.0, maximum=100.0)
    dntp = _number(dntp_conc, "dntp_conc", minimum=0.0, maximum=100.0)
    dna = _number(dna_conc, "dna_conc", minimum=0.0, maximum=10_000.0)
    temp = _number(temp_c, "temp_c", minimum=0.0, maximum=100.0)
    loop = _integer(max_loop, "max_loop", minimum=1, maximum=100)
    try:
        result = primer3.calc_end_stability(
            primer_text,
            template_text,
            mv_conc=mv,
            dv_conc=dv,
            dntp_conc=dntp,
            dna_conc=dna,
            temp_c=temp,
            max_loop=loop,
        )
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 could not compute end stability: {exc}") from exc
    found = bool(result.structure_found)
    return {
        "status": "completed",
        "operation": "calculate_primer3_end_stability",
        "primer_length": len(primer_text),
        "template_length": len(template_text),
        "structure_found": found,
        "tm_c": _float_or_none(result.tm) if found else None,
        "dg_cal_per_mol": _float_or_none(result.dg) if found else None,
        "dh_cal_per_mol": _float_or_none(result.dh) if found else None,
        "ds_cal_per_mol_k": _float_or_none(result.ds) if found else None,
        "conditions": {
            "mv_mM": mv,
            "dv_mM": dv,
            "dntp_mM": dntp,
            "dna_nM": dna,
            "temp_c": temp,
            "max_loop": loop,
        },
        "provenance": _provenance(method="calc_end_stability"),
    }


def _target(value: object, length: int) -> list[int] | None:
    if value is None:
        return None
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ConfigurationError("target must be (start, length) in template coordinates.")
    start = _integer(value[0], "target start", minimum=0, maximum=length - 1)
    span = _integer(value[1], "target length", minimum=1, maximum=length)
    if start + span > length:
        raise ConfigurationError("target extends past the end of the template.")
    return [start, span]


def _pair(raw: dict[str, Any], index: int) -> dict[str, Any]:
    left_start, left_length = raw[f"PRIMER_LEFT_{index}"]
    right_end, right_length = raw[f"PRIMER_RIGHT_{index}"]

    def oligo(sequence: str, start: int, length: int, prefix: str) -> dict[str, Any]:
        return {
            "sequence": str(sequence),
            "start": int(start),
            "length": int(length),
            "tm_c": _float_or_none(raw[f"PRIMER_{prefix}_{index}_TM"]),
            "gc_percent": _float_or_none(raw[f"PRIMER_{prefix}_{index}_GC_PERCENT"]),
        }

    return {
        "rank": index + 1,
        "product_size": int(raw[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"]),
        "penalty": _float_or_none(raw[f"PRIMER_PAIR_{index}_PENALTY"]),
        "left": oligo(raw[f"PRIMER_LEFT_{index}_SEQUENCE"], left_start, left_length, "LEFT"),
        "right": oligo(
            raw[f"PRIMER_RIGHT_{index}_SEQUENCE"],
            int(right_end) - int(right_length) + 1,
            right_length,
            "RIGHT",
        ),
    }


def design_primer3_pairs(
    template: str,
    *,
    target: tuple[int, int] | None = None,
    min_product: int = 100,
    max_product: int = 300,
    num_return: int = 5,
    min_size: int = 18,
    opt_size: int = 20,
    max_size: int = 25,
    min_tm: float = 57.0,
    opt_tm: float = 60.0,
    max_tm: float = 63.0,
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design PCR primer pairs for a DNA template with primer3-py ``design_primers``.

    ``target`` is an optional ``(start, length)`` region that the amplicon must
    cover, in 0-based template coordinates. Each returned pair lists the left and
    right oligos with their 0-based 5' start, and the product size. Pairs are in
    the order primer3 ranks them, from lowest to highest penalty.
    """

    primer3 = require("primer3", "bio-primer3")
    text = _acgtn(template, "template", maximum=_MAX_TEMPLATE, minimum=50)
    region = _target(target, len(text))
    lower = _integer(min_product, "min_product", minimum=50, maximum=len(text))
    upper = _integer(max_product, "max_product", minimum=lower, maximum=len(text))
    count = _integer(num_return, "num_return", minimum=1, maximum=50)
    min_len = _integer(min_size, "min_size", minimum=10, maximum=35)
    opt_len = _integer(opt_size, "opt_size", minimum=min_len, maximum=35)
    max_len = _integer(max_size, "max_size", minimum=opt_len, maximum=35)
    low_tm = _number(min_tm, "min_tm", minimum=30.0, maximum=90.0)
    mid_tm = _number(opt_tm, "opt_tm", minimum=low_tm, maximum=90.0)
    high_tm = _number(max_tm, "max_tm", minimum=mid_tm, maximum=90.0)
    low_gc = _number(min_gc, "min_gc", minimum=0.0, maximum=100.0)
    high_gc = _number(max_gc, "max_gc", minimum=low_gc, maximum=100.0)
    seq_args: dict[str, Any] = {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": text}
    if region is not None:
        seq_args["SEQUENCE_TARGET"] = region
    global_args: dict[str, Any] = {
        "PRIMER_PRODUCT_SIZE_RANGE": [[lower, upper]],
        "PRIMER_NUM_RETURN": count,
        "PRIMER_MIN_SIZE": min_len,
        "PRIMER_OPT_SIZE": opt_len,
        "PRIMER_MAX_SIZE": max_len,
        "PRIMER_MIN_TM": low_tm,
        "PRIMER_OPT_TM": mid_tm,
        "PRIMER_MAX_TM": high_tm,
        "PRIMER_MIN_GC": low_gc,
        "PRIMER_MAX_GC": high_gc,
    }
    try:
        raw = primer3.design_primers(seq_args, global_args)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 could not design primers: {exc}") from exc
    if raw.get("PRIMER_ERROR"):
        raise ConfigurationError(f"primer3 rejected the design: {raw['PRIMER_ERROR']}")
    found = int(raw.get("PRIMER_PAIR_NUM_RETURNED", 0))
    return {
        "status": "completed",
        "operation": "design_primer3_pairs",
        "template_length": len(text),
        "target": region,
        "product_size_range": [lower, upper],
        "pair_count": found,
        "pairs": [_pair(raw, index) for index in range(found)],
        "provenance": _provenance(method="design_primers"),
    }


def sanitize_primer3_sequence(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Replace IUPAC ambiguity codes with N, as primer3's sequence sanitizer does.

    Standard A, C, G, T, N letters are kept in their original case. Characters
    that are not IUPAC codes, such as X, Z or gaps, raise ``ConfigurationError``.
    ``changed_count`` counts the positions that were replaced.
    """

    primer3 = require("primer3", "bio-primer3")
    if not isinstance(sequence, str) or not sequence or len(sequence) > _MAX_SANITIZE:
        raise ConfigurationError(
            f"sequence must be a non-empty string of at most {_MAX_SANITIZE} characters."
        )
    try:
        cleaned = str(primer3.p3helpers.sanitize_sequence(sequence))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"primer3 rejected the sequence: {exc}") from exc
    return {
        "status": "completed",
        "operation": "sanitize_primer3_sequence",
        "input_length": len(sequence),
        "sanitized": cleaned,
        "changed_count": sum(
            1 for before, after in zip(sequence, cleaned, strict=True) if before != after
        ),
        "provenance": _provenance(method="p3helpers.sanitize_sequence"),
    }


__all__ = [
    "calculate_primer3_end_stability",
    "design_primer3_pairs",
    "sanitize_primer3_sequence",
]

_register_public_verbose(globals(), tuple(__all__))
