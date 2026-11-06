"""Native bench calculations, protocol planning, Sanger trace handling and small-read assembly.

Pure standard library except for nearest-neighbour Tm, molecular weight and extinction coefficients,
which are delegated to :mod:`dnakit.thermodynamics`. Every function is deterministic and bounded.
"""

from __future__ import annotations

import csv
import html
import io
import json
import math
import re
import struct
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._candidate_runtime import candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError

AVOGADRO = 6.02214076e23
BP_MASS = 660.0  # average g/mol per base pair of dsDNA
NT_MASS = 330.0  # average g/mol per nucleotide of ssDNA
A260_FACTOR = {"dsDNA": 50.0, "ssDNA": 33.0, "RNA": 40.0}  # ng/uL per A260 unit (1 cm path)
_COMP = str.maketrans("ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb")


def run_op(
    operation: str, fn: Any, value: object, options: Mapping[str, object] | None
) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError(
            "options must be a mapping or None.",
            code="INVALID_OPTION",
            hint="Pass options={'key': value}.",
        )
    opts = dict(options or {})
    warnings: list[str] = []
    result = fn(value, opts, warnings)
    return candidate_ok("molecular_lab", operation, value, result, options=opts, warnings=warnings)


# --------------------------------------------------------------------------- helpers


def _num(
    source: Mapping[str, object],
    key: str,
    *,
    positive: bool = True,
    required: bool = True,
    default: float | None = None,
) -> float | None:
    value = source.get(key, default)
    if value is None:
        if required:
            raise ConfigurationError(
                f"parameter {key!r} is required.",
                code="MISSING_OPTION",
                hint=f"Pass {key} in the input mapping or in options.",
            )
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"parameter {key!r} must be a finite number; got {value!r}.", code="INVALID_OPTION"
        )
    if positive and value <= 0:
        raise ConfigurationError(
            f"parameter {key!r} must be > 0; got {value!r}.", code="INVALID_OPTION"
        )
    if not positive and value < 0:
        raise ConfigurationError(
            f"parameter {key!r} must be >= 0; got {value!r}.", code="INVALID_OPTION"
        )
    return float(value)


def _int(
    source: Mapping[str, object], key: str, default: int | None, *, minimum: int, maximum: int
) -> int:
    value = source.get(key, default)
    if value is None:
        raise ConfigurationError(f"parameter {key!r} is required.", code="MISSING_OPTION")
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"parameter {key!r} must be an integer in [{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
        )
    return value


def _choice(
    source: Mapping[str, object], key: str, default: str | None, allowed: Sequence[str]
) -> str:
    value = source.get(key, default)
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(
            f"parameter {key!r} must be one of {', '.join(allowed)}; got {value!r}.",
            code="INVALID_OPTION",
        )
    return value


def _flag(source: Mapping[str, object], key: str, default: bool) -> bool:
    value = source.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"parameter {key!r} must be a boolean; got {value!r}.", code="INVALID_OPTION"
        )
    return value


def _params(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Parameters from the input mapping, overridden by ``options``."""

    merged: dict[str, object] = dict(value) if isinstance(value, Mapping) else {}
    merged.update(options)
    return merged


def _r(x: float | None, digits: int = 6) -> float | None:
    if x is None:
        return None
    if x == 0 or not math.isfinite(x):
        return x
    return float(f"{x:.{digits}g}")


def dna(
    value: object, name: str = "sequence", *, strict: bool = True, maximum: int = 5_000_000
) -> str:
    if isinstance(value, Mapping):
        for key in ("sequence", "seq"):
            if key in value:
                value = value[key]
                break
    if not isinstance(value, str):
        raise InputFormatError(
            f"{name} must be a DNA string; got {type(value).__name__}.", code="INVALID_INPUT"
        )
    text = "".join(value.split()).upper().replace("U", "T")
    if text.startswith(">"):
        raise InputFormatError(
            f"{name} looks like FASTA; pass bare sequences or parse it first.", code="INVALID_INPUT"
        )
    if not text:
        raise InputFormatError(f"{name} is empty.", code="EMPTY_INPUT")
    pattern = r"[ACGT]+" if strict else r"[ACGTRYKMSWBDHVN]+"
    if not re.fullmatch(pattern, text):
        raise InputFormatError(
            f"{name} must contain only {'A, C, G and T' if strict else 'IUPAC DNA letters'}.",
            code="INVALID_SEQUENCE",
            hint="Remove gaps, numbers and ambiguity codes." if strict else None,
        )
    if len(text) > maximum:
        raise InputFormatError(f"{name} is longer than {maximum} bases.", code="INPUT_TOO_LARGE")
    return text


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _sequences(value: object, limit: int = 1000) -> list[tuple[str, str]]:
    """(name, sequence) pairs from a string, FASTA text, a list, or a mapping name -> sequence."""

    items: list[tuple[str, str]] = []
    if isinstance(value, str):
        if value.lstrip().startswith(">"):
            name, parts = None, []
            for line in value.splitlines():
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        items.append((name, "".join(parts)))
                    name, parts = (line[1:].split() or [f"seq{len(items) + 1}"])[0], []
                else:
                    parts.append(line)
            if name is not None:
                items.append((name, "".join(parts)))
        else:
            items = [("sequence", value)]
    elif isinstance(value, Mapping):
        if "sequence" in value or "seq" in value:
            items = [
                (
                    str(value.get("id", value.get("name", "sequence"))),
                    str(value.get("sequence", value.get("seq"))),
                )
            ]
        else:
            items = [(str(k), v) for k, v in value.items()]  # type: ignore[misc]
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        for i, item in enumerate(value):
            if isinstance(item, str):
                items.append((f"seq{i + 1}", item))
            elif isinstance(item, Mapping):
                items.append(
                    (
                        str(item.get("id", item.get("name", f"seq{i + 1}"))),
                        item.get("sequence", item.get("seq")),
                    )
                )  # type: ignore[arg-type]
            elif isinstance(item, Sequence) and len(item) == 2:
                items.append((str(item[0]), item[1]))
            else:
                raise InputFormatError(
                    f"item {i + 1} is not a sequence string, mapping or (id, sequence) pair.",
                    code="INVALID_INPUT",
                )
    else:
        raise InputFormatError(
            f"unsupported input type {type(value).__name__}.",
            code="INVALID_INPUT",
            hint="Pass a sequence, FASTA text, a list or a mapping.",
        )
    if not items:
        raise InputFormatError("no sequences found.", code="EMPTY_INPUT")
    if len(items) > limit:
        raise InputFormatError(f"more than {limit} sequences.", code="INPUT_TOO_LARGE")
    for name, seq in items:
        if not isinstance(seq, str):
            raise InputFormatError(f"sequence {name!r} must be a string.", code="INVALID_INPUT")
    return items


# --------------------------------------------------------------------------- thermodynamics wrappers


def _conditions(params: Mapping[str, object]) -> Any:
    from dnakit.thermodynamics import ThermodynamicConditions

    sodium = _num(params, "sodium_mm", default=50.0)
    strand = _num(params, "primer_nm", default=250.0)
    assert sodium is not None and strand is not None
    return ThermodynamicConditions(
        sodium_molar=sodium / 1000.0, strand_concentration_molar=strand * 1e-9
    )


def nn_tm(seq: str, params: Mapping[str, object]) -> float:
    """Nearest-neighbour Tm (SantaLucia 1998 unified parameters, monovalent salt correction) via dnakit.thermodynamics."""

    from dnakit.core import DNASequence
    from dnakit.thermodynamics import melting_temperature

    return float(
        melting_temperature(
            DNASequence(seq), conditions=_conditions(params), verbose=False
        ).tm_celsius
    )


def basic_tm(seq: str) -> float:
    """Basic GC rule for oligonucleotides of 14 nt or more: Tm = 64.9 + 41 (G + C - 16.4) / N."""

    if len(seq) < 14:
        raise ConfigurationError(
            "the 'basic' Tm formula is defined for oligos of 14 nt or more.",
            code="INVALID_OPTION",
            hint="Use method='wallace' (2-13 nt) or 'nearest_neighbor'.",
        )
    gc = seq.count("G") + seq.count("C")
    return 64.9 + 41.0 * (gc - 16.4) / len(seq)


def wallace_tm(seq: str) -> float:
    if not 2 <= len(seq) <= 13:
        raise ConfigurationError(
            "the Wallace rule is defined for 2-13 nt oligos.", code="INVALID_OPTION"
        )
    return 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (seq.count("G") + seq.count("C"))


def _tm(seq: str, method: str, params: Mapping[str, object]) -> float:
    if method == "wallace":
        return wallace_tm(seq)
    if method == "basic":
        return basic_tm(seq)
    if len(seq) < 2:
        raise ConfigurationError("nearest-neighbour Tm needs at least 2 nt.", code="INVALID_OPTION")
    return nn_tm(seq, params)


# --------------------------------------------------------------------------- calculate

_CALCULATIONS = (
    "dilution",
    "molarity",
    "copy_number",
    "ligation",
    "a260",
    "tm",
    "pcr_yield",
    "gel_size",
)


def _dilution(p: Mapping[str, object]) -> dict[str, Any]:
    c1 = _num(p, "c1", required=False)
    v1 = _num(p, "v1", required=False)
    c2 = _num(p, "c2", required=False)
    v2 = _num(p, "v2", required=False)
    known = [x is not None for x in (c1, v1, c2, v2)]
    if sum(known) != 3:
        raise ConfigurationError(
            "dilution needs exactly three of c1, v1, c2, v2.",
            code="INVALID_OPTION",
            hint="C1*V1 = C2*V2: supply three values to solve the fourth.",
        )
    if c1 is None:
        assert v1 and c2 and v2
        c1 = c2 * v2 / v1
        solved = "c1"
    elif v1 is None:
        assert c1 and c2 and v2
        v1 = c2 * v2 / c1
        solved = "v1"
    elif c2 is None:
        assert c1 and v1 and v2
        c2 = c1 * v1 / v2
        solved = "c2"
    else:
        assert c1 and v1 and c2
        v2 = c1 * v1 / c2
        solved = "v2"
    assert c1 is not None and v1 is not None and c2 is not None and v2 is not None
    if c2 > c1 * (1 + 1e-12):
        raise ConfigurationError(
            "the final concentration exceeds the stock concentration (c2 > c1); dilution is impossible.",
            code="INVALID_OPTION",
        )
    return {
        "solved": solved,
        "c1": _r(c1),
        "v1": _r(v1),
        "c2": _r(c2),
        "v2": _r(v2),
        "diluent_volume": _r(v2 - v1),
        "dilution_factor": _r(c1 / c2),
        "formula": "C1 * V1 = C2 * V2 (concentration and volume units must be consistent)",
    }


def _mw_per_unit(p: Mapping[str, object], strand: str) -> tuple[float, str]:
    if p.get("sequence") is not None:
        from dnakit.core import DNASequence
        from dnakit.thermodynamics import molecular_weight

        seq = dna(p["sequence"])
        mw = float(
            molecular_weight(
                DNASequence(seq), strand="double" if strand == "ds" else "single", verbose=False
            ).value_dalton
        )
        return mw, f"exact from sequence ({len(seq)} {'bp' if strand == 'ds' else 'nt'})"
    length = _num(p, "length")
    assert length is not None
    return (
        length * (BP_MASS if strand == "ds" else NT_MASS),
        f"{BP_MASS if strand == 'ds' else NT_MASS:g} g/mol per {'bp' if strand == 'ds' else 'nt'} x length",
    )


def _molarity(p: Mapping[str, object]) -> dict[str, Any]:
    strand = _choice(p, "strand", "ds", ("ds", "ss"))
    mw, basis = _mw_per_unit(p, strand)
    ng_ul = _num(p, "mass_concentration_ng_per_ul", required=False)
    nm = _num(p, "molar_concentration_nm", required=False)
    if (ng_ul is None) == (nm is None):
        raise ConfigurationError(
            "give exactly one of mass_concentration_ng_per_ul or molar_concentration_nm.",
            code="INVALID_OPTION",
        )
    if ng_ul is not None:
        # ng/uL = 1e-3 g/L ; mol/L = (ng/uL * 1e-3) / MW ; nM = mol/L * 1e9
        nm_out = ng_ul * 1e-3 / mw * 1e9
        return {
            "molecular_weight": _r(mw),
            "mw_basis": basis,
            "mass_concentration_ng_per_ul": ng_ul,
            "molar_concentration_nm": _r(nm_out),
            "pmol_per_ul": _r(nm_out / 1000.0),
            "formula": "nM = (ng/uL * 1e-3 / MW) * 1e9",
        }
    assert nm is not None
    return {
        "molecular_weight": _r(mw),
        "mw_basis": basis,
        "molar_concentration_nm": nm,
        "mass_concentration_ng_per_ul": _r(nm * 1e-9 * mw / 1e-3),
        "formula": "ng/uL = nM * 1e-9 * MW / 1e-3",
    }


def _copy_number(p: Mapping[str, object]) -> dict[str, Any]:
    strand = _choice(p, "strand", "ds", ("ds", "ss"))
    mw, basis = _mw_per_unit(p, strand)
    ng = _num(p, "mass_ng")
    assert ng is not None
    copies = ng * 1e-9 / mw * AVOGADRO
    return {
        "molecular_weight": _r(mw),
        "mw_basis": basis,
        "mass_ng": ng,
        "copies": _r(copies, 9),
        "formula": "copies = mass_g / MW * 6.02214076e23",
    }


def _ligation(p: Mapping[str, object]) -> dict[str, Any]:
    vec_bp = _num(p, "vector_bp")
    ins_bp = _num(p, "insert_bp")
    vec_ng = _num(p, "vector_ng")
    ratio = _num(p, "insert_to_vector_ratio", default=3.0)
    assert vec_bp and ins_bp and vec_ng and ratio
    ins_ng = ratio * vec_ng * ins_bp / vec_bp
    out: dict[str, Any] = {
        "insert_ng": _r(ins_ng),
        "vector_ng": vec_ng,
        "insert_to_vector_ratio": ratio,
        "vector_pmol": _r(vec_ng * 1e-9 / (vec_bp * BP_MASS) * 1e12),
        "insert_pmol": _r(ins_ng * 1e-9 / (ins_bp * BP_MASS) * 1e12),
        "formula": "insert_ng = ratio * vector_ng * insert_bp / vector_bp (equal ratio of ends, 660 g/mol per bp)",
    }
    vc = _num(p, "vector_conc_ng_per_ul", required=False)
    ic = _num(p, "insert_conc_ng_per_ul", required=False)
    if vc:
        out["vector_volume_ul"] = _r(vec_ng / vc)
    if ic:
        out["insert_volume_ul"] = _r(ins_ng / ic)
    return out


def _a260(p: Mapping[str, object]) -> dict[str, Any]:
    a = _num(p, "a260", positive=False)
    dil = _num(p, "dilution_factor", default=1.0)
    path = _num(p, "path_length_cm", default=1.0)
    assert a is not None and dil and path
    if p.get("sequence") is not None:
        from dnakit.core import DNASequence
        from dnakit.thermodynamics import extinction_coefficient_260nm, molecular_weight

        seq = dna(p["sequence"])
        eps = float(
            extinction_coefficient_260nm(DNASequence(seq), verbose=False).value_m_inverse_cm_inverse
        )
        mw = float(molecular_weight(DNASequence(seq), verbose=False).value_dalton)
        molar = a * dil / (eps * path)
        return {
            "method": "Beer-Lambert with the theoretical oligo extinction coefficient",
            "extinction_coefficient": eps,
            "molecular_weight": _r(mw),
            "molar_concentration_um": _r(molar * 1e6),
            "mass_concentration_ng_per_ul": _r(molar * mw * 1e3),
            "formula": "c = A260 * dilution / (epsilon * path)",
        }
    kind = _choice(p, "nucleic_acid", "dsDNA", tuple(A260_FACTOR))
    conc = a * A260_FACTOR[kind] * dil / path
    return {
        "method": f"{kind}: {A260_FACTOR[kind]:g} ng/uL per A260 unit at 1 cm",
        "mass_concentration_ng_per_ul": _r(conc),
        "formula": "c = A260 * factor * dilution / path_length",
    }


def _tm_calc(p: Mapping[str, object]) -> dict[str, Any]:
    seq = dna(p.get("sequence"), "sequence")
    method = _choice(p, "method", "nearest_neighbor", ("nearest_neighbor", "wallace", "basic"))
    tm = _tm(seq, method, p)
    out: dict[str, Any] = {
        "sequence": seq,
        "length": len(seq),
        "method": method,
        "tm_celsius": _r(tm, 5),
        "gc_fraction": _r((seq.count("G") + seq.count("C")) / len(seq)),
    }
    if method == "nearest_neighbor":
        out["conditions"] = {
            "sodium_mm": _num(p, "sodium_mm", default=50.0),
            "primer_nm": _num(p, "primer_nm", default=250.0),
        }
        out["model"] = (
            "SantaLucia 1998 unified nearest-neighbour parameters with a monovalent-salt entropy correction (dnakit.thermodynamics)"
        )
    elif method == "basic":
        out["formula"] = "Tm = 64.9 + 41 * (G + C - 16.4) / N"
    else:
        out["formula"] = "Tm = 2 (A + T) + 4 (G + C)"
    return out


def _pcr_yield(p: Mapping[str, object]) -> dict[str, Any]:
    n0 = _num(p, "template_copies")
    cycles = _int(p, "cycles", None, minimum=1, maximum=100)
    eff = _num(p, "efficiency", default=1.0)
    assert n0 and eff is not None
    if not 0 < eff <= 1:
        raise ConfigurationError("efficiency must be in (0, 1].", code="INVALID_OPTION")
    out = n0 * (1 + eff) ** cycles
    return {
        "template_copies": n0,
        "cycles": cycles,
        "efficiency": eff,
        "product_copies": _r(out, 9),
        "fold_amplification": _r((1 + eff) ** cycles, 9),
        "formula": "N = N0 * (1 + E)^cycles (exponential phase only; no plateau)",
    }


def _gel_size(p: Mapping[str, object]) -> dict[str, Any]:
    ladder = p.get("ladder")
    if not isinstance(ladder, Sequence) or isinstance(ladder, (str, bytes)) or len(ladder) < 2:
        raise ConfigurationError(
            "ladder must be a list of at least two {'bp', 'distance'} points.",
            code="MISSING_OPTION",
        )
    xs, ys = [], []
    for point in ladder:
        if not isinstance(point, Mapping):
            raise ConfigurationError(
                "ladder points must be {'bp': size, 'distance': migration} mappings.",
                code="INVALID_OPTION",
            )
        bp, dist = _num(point, "bp"), _num(point, "distance", positive=False)
        assert bp is not None and dist is not None
        xs.append(dist)
        ys.append(math.log10(bp))
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        raise ConfigurationError("ladder distances must not all be equal.", code="INVALID_OPTION")
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=False)) / sxx
    intercept = my - slope * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys, strict=False))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    raw = p.get("distances", p.get("distance"))
    if raw is None:
        raise ConfigurationError(
            "distances (list of migration distances) is required.", code="MISSING_OPTION"
        )
    distances = [raw] if isinstance(raw, (int, float)) and not isinstance(raw, bool) else list(raw)  # type: ignore[arg-type]
    sizes = []
    lo, hi = min(xs), max(xs)
    outside = 0
    for d in distances:
        if isinstance(d, bool) or not isinstance(d, (int, float)):
            raise ConfigurationError("distances must be numbers.", code="INVALID_OPTION")
        if not lo <= d <= hi:
            outside += 1
        sizes.append(_r(10 ** (intercept + slope * d), 6))
    return {
        "slope_log10bp_per_distance": _r(slope),
        "intercept_log10bp": _r(intercept),
        "r_squared": _r(r2),
        "estimated_bp": sizes,
        "extrapolated_points": outside,
        "formula": "least-squares line log10(bp) = a + b * distance through the ladder bands; sizes outside the ladder range are extrapolated",
    }


def op_calculate(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    params = _params(value, options)
    which = params.get("calculation")
    if which not in _CALCULATIONS:
        raise ConfigurationError(
            f"options['calculation'] must be one of {', '.join(_CALCULATIONS)}; got {which!r}.",
            code="MISSING_OPTION",
            hint="Example: calculate({'c1': 10, 'c2': 1, 'v2': 100}, options={'calculation': 'dilution'}).",
        )
    fn = {
        "dilution": _dilution,
        "molarity": _molarity,
        "copy_number": _copy_number,
        "ligation": _ligation,
        "a260": _a260,
        "tm": _tm_calc,
        "pcr_yield": _pcr_yield,
        "gel_size": _gel_size,
    }[which]
    return {"calculation": which, **fn(params)}


# --------------------------------------------------------------------------- analyze / validate


def _longest_run(seq: str) -> int:
    return max(len(m.group(0)) for m in re.finditer(r"(.)\1*", seq))


def _analyze_one(seq: str, params: Mapping[str, object]) -> dict[str, Any]:
    from dnakit.core import DNASequence
    from dnakit.thermodynamics import extinction_coefficient_260nm, molecular_weight

    method = _choice(
        params, "tm_method", "nearest_neighbor", ("nearest_neighbor", "wallace", "basic")
    )
    n = len(seq)
    counts = {b: seq.count(b) for b in "ACGT"}
    gc = counts["G"] + counts["C"]
    out: dict[str, Any] = {
        "length": n,
        "composition": counts,
        "gc_fraction": _r(gc / n),
        "reverse_complement": revcomp(seq),
        "is_self_complementary": seq == revcomp(seq),
        "longest_homopolymer": _longest_run(seq),
        "ends_with_gc_clamp": seq[-1] in "GC",
    }
    obj = DNASequence(seq)
    out["molecular_weight_da"] = _r(float(molecular_weight(obj, verbose=False).value_dalton), 8)
    if n >= 2:
        eps = float(extinction_coefficient_260nm(obj, verbose=False).value_m_inverse_cm_inverse)
        out["extinction_coefficient_260nm"] = eps
        out["nmol_per_od260"] = _r(1e6 / eps)
        out["ug_per_od260"] = _r(1e3 * out["molecular_weight_da"] / eps)
    try:
        out["tm_celsius"] = _r(_tm(seq, method, params), 5)
        out["tm_method"] = method
    except ConfigurationError as exc:
        out["tm_celsius"] = None
        out["tm_note"] = exc.message
    return out


def op_analyze(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    items = _sequences(value)
    rows = {}
    for name, seq in items:
        rows[name] = _analyze_one(dna(seq, name, maximum=100_000), options)
    return {
        "sequences": rows,
        "count": len(rows),
        "definition": "Composition, GC fraction, reverse complement, anhydrous single-strand molecular weight, theoretical 260 nm extinction "
        "coefficient (nearest-neighbour hypochromicity), nmol and ug per OD260 and melting temperature. Tm defaults to the SantaLucia 1998 "
        "nearest-neighbour model with sodium_mm=50 and primer_nm=250 (options tm_method, sodium_mm, primer_nm).",
    }


def _check(name: str, ok: bool, detail: str) -> dict[str, Any]:
    return {"check": name, "passed": ok, "detail": detail}


def _validate_primer(seq: str, p: Mapping[str, object]) -> dict[str, Any]:
    min_len = _int(p, "min_length", 18, minimum=1, maximum=200)
    max_len = _int(p, "max_length", 30, minimum=min_len, maximum=200)
    gc_lo = _num(p, "min_gc", positive=False, default=0.40)
    gc_hi = _num(p, "max_gc", positive=False, default=0.60)
    tm_lo = _num(p, "min_tm", positive=False, default=52.0)
    tm_hi = _num(p, "max_tm", positive=False, default=65.0)
    max_run = _int(p, "max_homopolymer", 4, minimum=1, maximum=50)
    assert gc_lo is not None and gc_hi is not None and tm_lo is not None and tm_hi is not None
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / n
    tail = seq[-5:]
    clamp = sum(1 for b in tail if b in "GC")
    method = _choice(p, "tm_method", "nearest_neighbor", ("nearest_neighbor", "wallace", "basic"))
    checks = [
        _check("length", min_len <= n <= max_len, f"{n} nt (allowed {min_len}-{max_len})"),
        _check(
            "gc_content",
            gc_lo <= gc <= gc_hi,
            f"{gc * 100:.1f}% (allowed {gc_lo * 100:g}-{gc_hi * 100:g}%)",
        ),
        _check(
            "three_prime_gc_clamp",
            seq[-1] in "GC" and clamp <= 3,
            f"last 5 nt {tail} has {clamp} G/C; 3' base {seq[-1]} (want 3' G/C and <= 3 G/C in the last 5)",
        ),
        _check(
            "homopolymer_runs",
            _longest_run(seq) <= max_run,
            f"longest run {_longest_run(seq)} (allowed <= {max_run})",
        ),
        _check(
            "dinucleotide_repeats",
            re.search(r"(..)\1{3,}", seq) is None,
            "no dinucleotide unit repeated 4+ times"
            if not re.search(r"(..)\1{3,}", seq)
            else "contains a dinucleotide unit repeated 4+ times",
        ),
    ]
    tm: float | None
    try:
        tm = _tm(seq, method, p)
        checks.append(
            _check(
                "melting_temperature",
                tm_lo <= tm <= tm_hi,
                f"{tm:.1f} C ({method}; allowed {tm_lo:g}-{tm_hi:g} C)",
            )
        )
    except ConfigurationError as exc:
        tm = None
        checks.append(_check("melting_temperature", False, exc.message))
    return {
        "sequence": seq,
        "tm_celsius": _r(tm, 5),
        "passed": all(c["passed"] for c in checks),
        "checks": checks,
    }


def op_validate(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    params = _params(None, options)
    if isinstance(value, Mapping) and "forward" in value and "reverse" in value:
        fwd = _validate_primer(dna(value["forward"], "forward", maximum=200), params)
        rev = _validate_primer(dna(value["reverse"], "reverse", maximum=200), params)
        max_dtm = _num(params, "max_tm_difference", positive=False, default=5.0)
        assert max_dtm is not None
        pair_checks = []
        if fwd["tm_celsius"] is not None and rev["tm_celsius"] is not None:
            d = abs(fwd["tm_celsius"] - rev["tm_celsius"])
            pair_checks.append(
                _check(
                    "tm_difference",
                    d <= max_dtm,
                    f"|Tm1 - Tm2| = {d:.1f} C (allowed <= {max_dtm:g})",
                )
            )
        ok = fwd["passed"] and rev["passed"] and all(c["passed"] for c in pair_checks)
        return {
            "kind": "primer_pair",
            "passed": ok,
            "forward": fwd,
            "reverse": rev,
            "pair_checks": pair_checks,
            "rules": _RULES,
        }
    rows = {
        name: _validate_primer(dna(seq, name, maximum=200), params)
        for name, seq in _sequences(value, 500)
    }
    return {
        "kind": "primers",
        "passed": all(r["passed"] for r in rows.values()),
        "primers": rows,
        "rules": _RULES,
    }


_RULES = (
    "Rule-of-thumb PCR primer criteria (adjustable): length 18-30 nt, GC 40-60%, Tm 52-65 C (nearest-neighbour, 50 mM Na+, 250 nM), "
    "3' terminal G/C with at most 3 G/C in the last 5 bases, homopolymer runs <= 4, no dinucleotide repeated 4+ times, pair Tm "
    "difference <= 5 C. Passing does not predict amplification; secondary structure and off-target binding are not assessed."
)


# --------------------------------------------------------------------------- parse

_PREFIX = {
    "": 1.0,
    "k": 1e3,
    "m": 1e-3,
    "u": 1e-6,
    "µ": 1e-6,
    "μ": 1e-6,
    "n": 1e-9,
    "p": 1e-12,
    "f": 1e-15,
}
_QTY = re.compile(r"^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*(.+?)\s*$")


def _unit(text: str) -> tuple[float, str]:
    for base in ("mol", "g", "L", "l", "M"):
        if text.endswith(base):
            prefix = text[: -len(base)]
            if prefix in _PREFIX:
                return _PREFIX[prefix], "L" if base == "l" else base
    raise InputFormatError(
        f"unrecognised unit {text!r}.",
        code="INVALID_UNIT",
        hint="Use g, L, mol, M with prefixes k, m, u, n, p, f (e.g. ng/uL, mM, pmol).",
    )


def parse_quantity(text: str) -> dict[str, Any]:
    if not isinstance(text, str):
        raise InputFormatError(
            "quantity must be a string such as '50 ng/uL'.", code="INVALID_INPUT"
        )
    m = _QTY.match(text)
    if not m:
        raise InputFormatError(
            f"cannot parse quantity {text!r}.",
            code="INVALID_QUANTITY",
            hint="Use '<number> <unit>', e.g. '50 ng/uL' or '10 mM'.",
        )
    value = float(m.group(1))
    unit = m.group(2).replace(" ", "")
    if unit.upper() == "X":
        return {
            "text": text,
            "value": value,
            "unit": "X",
            "dimension": "fold_concentration",
            "si_value": value,
            "si_unit": "x",
        }
    if "/" in unit:
        top, bottom = unit.split("/", 1)
        f1, b1 = _unit(top)
        f2, b2 = _unit(bottom)
        if (b1, b2) == ("g", "L"):
            dim, si_unit = "mass_concentration", "g/L"
        elif (b1, b2) == ("mol", "L"):
            dim, si_unit = "molar_concentration", "mol/L"
        else:
            raise InputFormatError(
                f"unsupported compound unit {unit!r}.",
                code="INVALID_UNIT",
                hint="Supported ratios: mass/volume and amount/volume.",
            )
        si = value * f1 / f2
    else:
        factor, base = _unit(unit)
        dim, si_unit = {
            "g": ("mass", "g"),
            "L": ("volume", "L"),
            "mol": ("amount", "mol"),
            "M": ("molar_concentration", "mol/L"),
        }[base]
        si = value * factor
    return {
        "text": text,
        "value": value,
        "unit": unit,
        "dimension": dim,
        "si_value": _r(si, 12),
        "si_unit": si_unit,
    }


_PLATES = {96: (8, 12), 384: (16, 24)}


def _well(label: str, rows: int, cols: int) -> tuple[int, int]:
    m = re.fullmatch(r"([A-Za-z]{1,2})0*(\d{1,2})", label.strip())
    if not m:
        raise InputFormatError(
            f"invalid well {label!r}.", code="INVALID_WELL", hint="Wells look like A1 or H12."
        )
    letters = m.group(1).upper()
    row = 0
    for ch in letters:
        row = row * 26 + ord(ch) - 64
    row -= 1
    col = int(m.group(2)) - 1
    if not (0 <= row < rows and 0 <= col < cols):
        raise InputFormatError(
            f"well {label!r} is outside a {rows}x{cols} plate.", code="INVALID_WELL"
        )
    return row, col


def _row_label(row: int) -> str:
    return chr(65 + row) if row < 26 else chr(64 + row // 26) + chr(65 + row % 26)


def parse_wells(text: str, plate: int, order: str) -> list[str]:
    rows, cols = _PLATES[plate]
    out: list[str] = []
    for part in re.split(r"[,\s;]+", text.strip()):
        if not part:
            continue
        if ":" in part:
            a, b = part.split(":", 1)
            r0, c0 = _well(a, rows, cols)
            r1, c1 = _well(b, rows, cols)
            rlo, rhi, clo, chi = min(r0, r1), max(r0, r1), min(c0, c1), max(c0, c1)
            cells = (
                [(r, c) for r in range(rlo, rhi + 1) for c in range(clo, chi + 1)]
                if order == "row"
                else [(r, c) for c in range(clo, chi + 1) for r in range(rlo, rhi + 1)]
            )
        else:
            cells = [_well(part, rows, cols)]
        out.extend(f"{_row_label(r)}{c + 1}" for r, c in cells)
    if not out:
        raise InputFormatError("no wells found.", code="EMPTY_INPUT")
    return out


def op_parse(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    kind = _choice(options, "kind", "quantity", ("quantity", "wells"))
    if kind == "quantity":
        items = (
            [value]
            if isinstance(value, str)
            else list(value)
            if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray))
            else None
        )
        if items is None or not items:
            raise InputFormatError(
                "input must be a quantity string or a list of them.", code="INVALID_INPUT"
            )
        if len(items) > 10_000:
            raise InputFormatError("more than 10000 quantities.", code="INPUT_TOO_LARGE")
        return {
            "kind": kind,
            "quantities": [parse_quantity(x) for x in items],
            "definition": "'<number> <unit>' with units g, L, mol, M and prefixes k, m, u (or the micro sign), n, p, f; ratios mass/volume and amount/volume; 'X' for fold concentration. SI values are g, L, mol, g/L and mol/L.",
        }
    if not isinstance(value, str):
        raise InputFormatError(
            "wells input must be text such as 'A1:B3, D5'.", code="INVALID_INPUT"
        )
    plate = _int(options, "plate", 96, minimum=96, maximum=384)
    if plate not in _PLATES:
        raise ConfigurationError("plate must be 96 or 384.", code="INVALID_OPTION")
    order = _choice(options, "order", "row", ("row", "column"))
    wells = parse_wells(value, plate, order)
    return {
        "kind": kind,
        "plate": plate,
        "order": order,
        "wells": wells,
        "count": len(wells),
        "definition": f"Well lists and rectangular ranges ('A1:B3') on a {_PLATES[plate][0]}x{_PLATES[plate][1]} plate, expanded {order}-major.",
    }


# --------------------------------------------------------------------------- plan


_PCR_DEFAULT = [
    {"name": "reaction buffer", "stock": 5.0, "final": 1.0, "unit": "X"},
    {"name": "dNTP mix", "stock": 10.0, "final": 0.2, "unit": "mM each"},
    {"name": "forward primer", "stock": 10.0, "final": 0.5, "unit": "uM"},
    {"name": "reverse primer", "stock": 10.0, "final": 0.5, "unit": "uM"},
    {"name": "DNA polymerase", "stock": 5.0, "final": 0.025, "unit": "U/uL"},
]


def _plan_pcr(p: Mapping[str, object]) -> dict[str, Any]:
    volume = _num(p, "reaction_volume_ul", default=25.0)
    n = _int(p, "n_reactions", 1, minimum=1, maximum=100_000)
    overage = _num(p, "overage_fraction", positive=False, default=0.1)
    template = _num(p, "template_volume_ul", positive=False, default=1.0)
    assert volume and overage is not None and template is not None
    comps = p.get("components") or _PCR_DEFAULT
    if not isinstance(comps, Sequence) or isinstance(comps, (str, bytes)):
        raise ConfigurationError(
            "components must be a list of {name, stock, final, unit} mappings.",
            code="INVALID_OPTION",
        )
    rows = []
    used = template
    scale = n * (1 + overage)
    for c in comps:
        if not isinstance(c, Mapping) or "name" not in c:
            raise ConfigurationError(
                "each component needs name, stock and final.", code="INVALID_OPTION"
            )
        stock, final = _num(c, "stock"), _num(c, "final")
        assert stock and final
        if final > stock:
            raise ConfigurationError(
                f"component {c['name']!r}: final concentration exceeds the stock.",
                code="INVALID_OPTION",
            )
        v = final / stock * volume
        used += v
        rows.append(
            {
                "component": str(c["name"]),
                "stock": stock,
                "final": final,
                "unit": c.get("unit", ""),
                "volume_per_reaction_ul": _r(v),
                "master_mix_volume_ul": _r(v * scale),
            }
        )
    water = volume - used
    if water < -1e-9:
        raise ConfigurationError(
            f"components and template ({used:g} uL) exceed the reaction volume ({volume:g} uL).",
            code="INVALID_OPTION",
            hint="Use more concentrated stocks.",
        )
    rows.append(
        {
            "component": "water",
            "volume_per_reaction_ul": _r(max(water, 0.0)),
            "master_mix_volume_ul": _r(max(water, 0.0) * scale),
        }
    )
    mm_per_rxn = volume - template
    return {
        "protocol": "pcr",
        "reaction_volume_ul": volume,
        "n_reactions": n,
        "overage_fraction": overage,
        "master_mix_for_reactions": _r(scale),
        "components": rows,
        "master_mix_total_ul": _r(mm_per_rxn * scale),
        "master_mix_aliquot_per_tube_ul": _r(mm_per_rxn),
        "template_per_tube_ul": template,
        "definition": "volume = final / stock * reaction volume per component; water makes up the rest; the master mix (everything but the template) is scaled by n_reactions * (1 + overage_fraction); the template is added per tube.",
    }


def _plan_serial(p: Mapping[str, object]) -> dict[str, Any]:
    c0 = _num(p, "start_concentration")
    factor = _num(p, "dilution_factor", default=10.0)
    steps = _int(p, "steps", None, minimum=1, maximum=100)
    final_v = _num(p, "final_volume_ul", default=100.0)
    assert c0 and factor and final_v
    if factor <= 1:
        raise ConfigurationError("dilution_factor must be > 1.", code="INVALID_OPTION")
    transfer = final_v / factor
    rows = [{"tube": 0, "concentration": _r(c0), "note": "stock"}]
    for k in range(1, steps + 1):
        rows.append(
            {
                "tube": k,
                "concentration": _r(c0 / factor**k),
                "transfer_from_previous_ul": _r(transfer),
                "diluent_ul": _r(final_v - transfer),
                "total_volume_ul": final_v,
            }
        )
    return {
        "protocol": "serial_dilution",
        "dilution_factor": factor,
        "steps": steps,
        "tubes": rows,
        "total_diluent_ul": _r((final_v - transfer) * steps),
        "definition": "Each tube receives final_volume / factor of the previous tube plus diluent; concentration_k = start / factor^k.",
    }


def _plan_ligation(p: Mapping[str, object]) -> dict[str, Any]:
    out = _ligation(p)
    total = _num(p, "reaction_volume_ul", default=20.0)
    assert total
    used = (out.get("vector_volume_ul") or 0) + (out.get("insert_volume_ul") or 0)
    ligase_buffer = total / 10.0
    ligase = 1.0
    out.update(
        {
            "protocol": "ligation",
            "reaction_volume_ul": total,
            "ligase_buffer_10x_ul": _r(ligase_buffer),
            "ligase_ul": ligase,
        }
    )
    if "vector_volume_ul" in out and "insert_volume_ul" in out:
        water = total - used - ligase_buffer - ligase
        if water < -1e-9:
            raise ConfigurationError(
                "DNA, buffer and ligase exceed the reaction volume.",
                code="INVALID_OPTION",
                hint="Use more concentrated DNA or a larger volume.",
            )
        out["water_ul"] = _r(max(water, 0.0))
    out["definition"] = (
        "Insert mass from the insert:vector molar ratio (660 g/mol per bp); buffer is 1/10 of the volume (10x stock) and 1 uL ligase; water to volume."
    )
    return out


def op_plan(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    params = _params(value, options)
    protocol = _choice(params, "protocol", "pcr", ("pcr", "serial_dilution", "ligation"))
    return {"pcr": _plan_pcr, "serial_dilution": _plan_serial, "ligation": _plan_ligation}[
        protocol
    ](params)


# --------------------------------------------------------------------------- export


def op_export(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    fmt = _choice(options, "format", "csv", ("csv", "tsv", "json", "fasta"))
    records: list[dict[str, Any]]
    if isinstance(value, Mapping) and not any(k in value for k in ("records", "sequence")):
        records = [{"id": str(k), "sequence": v} for k, v in value.items()]
    elif isinstance(value, Mapping) and "records" in value:
        return op_export(value["records"], options, warnings)
    elif isinstance(value, Mapping):
        records = [dict(value)]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value:
        records = []
        for i, item in enumerate(value):
            if isinstance(item, Mapping):
                records.append(dict(item))
            elif isinstance(item, str):
                records.append({"id": f"seq{i + 1}", "sequence": item})
            else:
                raise InputFormatError(
                    f"record {i + 1} must be a mapping or a sequence string.", code="INVALID_INPUT"
                )
    else:
        raise InputFormatError(
            "input must be a list of records, a mapping id -> sequence, or one record.",
            code="INVALID_INPUT",
        )
    if len(records) > 100_000:
        raise InputFormatError("more than 100000 records.", code="INPUT_TOO_LARGE")
    columns = options.get("columns")
    if columns is None:
        columns = []
        for r in records:
            for k in r:
                if k not in columns:
                    columns.append(k)
    elif not isinstance(columns, Sequence) or isinstance(columns, (str, bytes)) or not columns:
        raise ConfigurationError(
            "columns must be a non-empty list of field names.", code="INVALID_OPTION"
        )
    if fmt in ("csv", "tsv"):
        buf = io.StringIO()
        writer = csv.writer(buf, delimiter="," if fmt == "csv" else "\t", lineterminator="\n")
        writer.writerow(columns)
        for r in records:
            writer.writerow([_cell(r.get(c)) for c in columns])
        text = buf.getvalue()
    elif fmt == "json":
        text = (
            json.dumps(
                [{c: r.get(c) for c in columns} for r in records],
                indent=2,
                ensure_ascii=False,
                default=str,
            )
            + "\n"
        )
    else:
        width = _int(options, "line_width", 80, minimum=10, maximum=1000)
        lines = []
        for i, r in enumerate(records):
            seq = r.get("sequence", r.get("seq"))
            if not isinstance(seq, str) or not seq:
                raise InputFormatError(
                    f"record {i + 1} has no 'sequence' for FASTA output.", code="INVALID_INPUT"
                )
            lines.append(f">{r.get('id', r.get('name', f'seq{i + 1}'))}")
            lines.extend(seq[j : j + width] for j in range(0, len(seq), width))
        text = "\n".join(lines) + "\n"
    target = options.get("output_path")
    written = None
    if target is not None:
        path = Path(str(target))
        if not path.parent.is_dir():
            raise ConfigurationError(
                f"directory {path.parent} does not exist.", code="INVALID_OPTION"
            )
        if path.exists() and not _flag(options, "overwrite", False):
            raise ConfigurationError(
                f"{path} already exists.",
                code="OUTPUT_EXISTS",
                hint="Pass overwrite=True to replace it.",
            )
        path.write_text(text, encoding="utf-8")
        written = str(path)
    return {
        "format": fmt,
        "records": len(records),
        "columns": list(columns),
        "text": text,
        "output_path": written,
    }


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


# --------------------------------------------------------------------------- plot (virtual gel)


def op_plot(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    if isinstance(value, Mapping) and "lanes" in value:
        value = value["lanes"]
    if isinstance(value, Mapping):
        lanes = [(str(k), v) for k, v in value.items()]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        lanes = [(f"lane {i + 1}", v) for i, v in enumerate(value)]
    else:
        raise InputFormatError(
            "input must be a mapping lane name -> fragment sizes (bp) or a list of such lists.",
            code="INVALID_INPUT",
            hint="Example: {'ladder': [1000, 500, 100], 'PCR': [650]}",
        )
    if not lanes or len(lanes) > 40:
        raise InputFormatError("1 to 40 lanes are supported.", code="INVALID_INPUT")
    lo = _num(options, "min_bp", default=100.0)
    hi = _num(options, "max_bp", default=10000.0)
    assert lo and hi
    if hi <= lo:
        raise ConfigurationError("max_bp must exceed min_bp.", code="INVALID_OPTION")
    sizes_by_lane = []
    for name, sizes in lanes:
        if not isinstance(sizes, Sequence) or isinstance(sizes, (str, bytes)):
            raise InputFormatError(
                f"lane {name!r} must be a list of sizes in bp.", code="INVALID_INPUT"
            )
        clean = []
        for s in sizes:
            if isinstance(s, bool) or not isinstance(s, (int, float)) or s <= 0:
                raise InputFormatError(
                    f"lane {name!r}: fragment sizes must be positive numbers.", code="INVALID_INPUT"
                )
            clean.append(float(s))
        sizes_by_lane.append((name, clean))
    lane_w, top, gel_h, left = 54, 50, 320, 60
    width = left + lane_w * len(lanes) + 20
    height = top + gel_h + 30
    span = math.log10(hi) - math.log10(lo)

    def y_of(bp: float) -> float:
        frac = (math.log10(hi) - math.log10(min(max(bp, lo), hi))) / span
        return top + 15 + (gel_h - 30) * frac

    out_of_range = sum(1 for _, ss in sizes_by_lane for s in ss if not lo <= s <= hi)
    if out_of_range:
        warnings.append(
            f"{out_of_range} fragments lie outside [{lo:g}, {hi:g}] bp and are drawn at the gel edge."
        )
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="10">',
        "<title>Virtual agarose gel</title>",
        f'<rect width="{width}" height="{height}" fill="white"/>',
        f'<rect x="{left - 8}" y="{top}" width="{lane_w * len(lanes) + 16}" height="{gel_h}" fill="#1c1c1c"/>',
    ]
    top_mass = max((s for _, ss in sizes_by_lane for s in ss), default=1.0)
    for i, (name, sizes) in enumerate(sizes_by_lane):
        x = left + i * lane_w
        svg.append(
            f'<text x="{x + lane_w / 2 - 4}" y="{top - 8}" text-anchor="middle" fill="#222">{html.escape(name[:12])}</text>'
        )
        svg.append(
            f'<rect x="{x + 4}" y="{top + 4}" width="{lane_w - 14}" height="6" fill="#444"/>'
        )
        for s in sizes:
            opacity = 0.35 + 0.65 * math.sqrt(s / top_mass)
            svg.append(
                f'<rect x="{x + 4}" y="{y_of(s) - 1.5:.1f}" width="{lane_w - 14}" height="3" fill="#f2f2f2" fill-opacity="{opacity:.2f}"><title>{s:g} bp</title></rect>'
            )
    marks = [m for m in (100, 200, 500, 1000, 2000, 5000, 10000) if lo <= m <= hi]
    for m in marks:
        svg.append(
            f'<text x="{left - 12}" y="{y_of(m) + 3:.1f}" text-anchor="end" fill="#222">{m}</text>'
        )
    svg.append("</svg>")
    text = "\n".join(svg)
    written = None
    target = options.get("output_path")
    if target is not None:
        path = Path(str(target))
        if path.suffix.lower() != ".svg" or not path.parent.is_dir():
            raise ConfigurationError(
                "output_path must be an .svg file in an existing directory.", code="INVALID_OPTION"
            )
        if path.exists() and not _flag(options, "overwrite", False):
            raise ConfigurationError(
                f"{path} already exists.",
                code="OUTPUT_EXISTS",
                hint="Pass overwrite=True to replace it.",
            )
        path.write_text(text, encoding="utf-8")
        written = str(path)
    return {
        "svg": text,
        "output_path": written,
        "lanes": len(lanes),
        "min_bp": lo,
        "max_bp": hi,
        "model": "Schematic gel: migration distance is linear in log10(size) between min_bp and max_bp (the usual semilog behaviour of agarose); band brightness grows with fragment mass for equimolar fragments. It is a drawing, not a simulation of a specific gel percentage or buffer.",
    }


# --------------------------------------------------------------------------- optimize (primer length for target Tm)


def op_optimize(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    template = dna(value, "template", maximum=100_000)
    target = _num(options, "target_tm", positive=False, default=60.0)
    min_len = _int(options, "min_length", 18, minimum=6, maximum=100)
    max_len = _int(options, "max_length", 30, minimum=min_len, maximum=100)
    direction = _choice(options, "direction", "forward", ("forward", "reverse"))
    clamp = _flag(options, "require_gc_clamp", False)
    method = _choice(options, "tm_method", "nearest_neighbor", ("nearest_neighbor", "basic"))
    assert target is not None
    if len(template) < min_len:
        raise InputFormatError(
            f"template ({len(template)} nt) is shorter than min_length={min_len}.",
            code="INPUT_TOO_SHORT",
        )
    candidates = []
    for n in range(min_len, min(max_len, len(template)) + 1):
        seq = template[:n] if direction == "forward" else revcomp(template[-n:])
        tm = _tm(seq, method, options)
        candidates.append(
            {
                "length": n,
                "primer": seq,
                "tm_celsius": _r(tm, 5),
                "tm_error": _r(abs(tm - target), 5),
                "gc_fraction": _r((seq.count("G") + seq.count("C")) / n),
                "gc_clamp": seq[-1] in "GC",
            }
        )
    pool = [c for c in candidates if c["gc_clamp"]] if clamp else candidates
    if not pool:
        raise InputFormatError(
            "no candidate length ends in G/C (require_gc_clamp).",
            code="NO_CANDIDATE",
            hint="Allow other lengths or disable require_gc_clamp.",
        )
    best = min(pool, key=lambda c: (c["tm_error"], c["length"]))
    return {
        "direction": direction,
        "target_tm": target,
        "tm_method": method,
        "best": best,
        "candidates": candidates,
        "definition": "Primers are the first (forward) or reverse complement of the last (reverse) min_length..max_length bases of the template; "
        "the length whose Tm is closest to target_tm wins (ties: shorter). Tm = nearest-neighbour (SantaLucia 1998, sodium_mm=50, primer_nm=250 by default) "
        "or the basic GC formula; no secondary-structure or specificity check.",
    }


# --------------------------------------------------------------------------- Sanger traces


_IUPAC2 = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
}
_ABIF_TYPES = {
    2: ("c", 1),
    4: ("h", 2),
    5: ("i", 4),
    7: ("f", 4),
    18: ("c", 1),
    19: ("c", 1),
    3: ("H", 2),
    1: ("B", 1),
    13: ("B", 1),
}


def _abif_directory(data: bytes) -> dict[str, Any]:
    if data[:4] != b"ABIF":
        raise InputFormatError("not an ABIF file (missing 'ABIF' marker).", code="INVALID_AB1")
    try:
        _name, _number, _etype, esize, nelem, _dsize, doffset, _handle = struct.unpack(
            ">4sI2H4I", data[6:34]
        )
    except struct.error as exc:
        raise InputFormatError("truncated ABIF header.", code="INVALID_AB1") from exc
    tags: dict[str, Any] = {}
    if esize != 28 or nelem > 10_000:
        raise InputFormatError("unexpected ABIF directory layout.", code="INVALID_AB1")
    for i in range(nelem):
        start = doffset + i * esize
        entry = data[start : start + esize]
        if len(entry) < esize:
            raise InputFormatError("truncated ABIF directory.", code="INVALID_AB1")
        tname, tnum, ttype, _tsize, tcount, tdatasize, tdataoffset, _h = struct.unpack(
            ">4sI2H4I", entry[:28]
        )
        key = f"{tname.decode('ascii', 'replace')}{tnum}"
        if key not in _WANTED_TAGS:
            continue
        raw = (
            entry[20 : 20 + tdatasize]
            if tdatasize <= 4
            else data[tdataoffset : tdataoffset + tdatasize]
        )
        if len(raw) < tdatasize:
            raise InputFormatError(f"truncated ABIF data for tag {key}.", code="INVALID_AB1")
        if ttype not in _ABIF_TYPES:
            continue
        code, width = _ABIF_TYPES[ttype]
        if code == "c":
            tags[key] = raw[:tcount].decode("latin-1")
        else:
            tags[key] = list(struct.unpack(f">{tcount}{code}", raw[: tcount * width]))
    return tags


_WANTED_TAGS = {
    "DATA9",
    "DATA10",
    "DATA11",
    "DATA12",
    "FWO_1",
    "PLOC1",
    "PLOC2",
    "PBAS1",
    "PBAS2",
    "PCON1",
    "PCON2",
}


def _trace_input(
    value: object, params: Mapping[str, object]
) -> tuple[dict[str, list[float]], list[int] | None, dict[str, Any]]:
    meta: dict[str, Any] = {}
    if isinstance(value, (str, Path)) and not (isinstance(value, str) and "\n" in value):
        path = Path(str(value))
        tags = _abif_directory(_read_bytes(path))
        order = tags.get("FWO_1")
        if not isinstance(order, str) or sorted(order) != list("ACGT"):
            raise InputFormatError(
                "ABIF file has no valid FWO_1 channel order.", code="INVALID_AB1"
            )
        traces = {}
        for k, base in enumerate(order):
            data = tags.get(f"DATA{9 + k}")
            if not data:
                raise InputFormatError(f"ABIF file lacks channel DATA{9 + k}.", code="INVALID_AB1")
            traces[base] = [float(x) for x in data]
        peaks = tags.get("PLOC2") or tags.get("PLOC1")
        meta["source"] = "ABIF file"
        meta["channel_order"] = order
        called = tags.get("PBAS2") or tags.get("PBAS1")
        if called:
            meta["instrument_basecalls"] = called
        return traces, [int(x) for x in peaks] if peaks else None, meta
    if not isinstance(value, Mapping) or "traces" not in value:
        raise InputFormatError(
            "input must be an .ab1 path or {'traces': {'A': [...], 'C': [...], 'G': [...], 'T': [...]}, 'peak_locations': [...] (optional)}.",
            code="INVALID_INPUT",
        )
    raw = value["traces"]
    if not isinstance(raw, Mapping) or set(raw) != set("ACGT"):
        raise InputFormatError(
            "traces must have exactly the channels A, C, G and T.", code="INVALID_INPUT"
        )
    traces = {}
    for base in "ACGT":
        series = raw[base]
        if (
            not isinstance(series, Sequence)
            or isinstance(series, (str, bytes))
            or any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in series)
        ):
            raise InputFormatError(
                f"channel {base} must be a list of numbers.", code="INVALID_INPUT"
            )
        traces[base] = [float(x) for x in series]
    lengths = {len(v) for v in traces.values()}
    if len(lengths) != 1 or lengths == {0}:
        raise InputFormatError(
            "all four channels must have the same non-zero length.", code="INVALID_INPUT"
        )
    if lengths.pop() > 2_000_000:
        raise InputFormatError("trace longer than 2,000,000 samples.", code="INPUT_TOO_LARGE")
    peaks_raw = value.get("peak_locations")
    peaks = None
    if peaks_raw is not None:
        if (
            not isinstance(peaks_raw, Sequence)
            or isinstance(peaks_raw, (str, bytes))
            or any(isinstance(x, bool) or not isinstance(x, int) for x in peaks_raw)
        ):
            raise InputFormatError(
                "peak_locations must be a list of integer sample indexes.", code="INVALID_INPUT"
            )
        peaks = list(peaks_raw)
    meta["source"] = "mapping"
    return traces, peaks, meta


def _read_bytes(path: Path) -> bytes:
    if not path.is_file():
        raise InputFormatError(f"file not found: {path}", code="INPUT_FILE_NOT_FOUND")
    if path.stat().st_size > 50_000_000:
        raise InputFormatError("file larger than 50 MB.", code="INPUT_TOO_LARGE")
    return path.read_bytes()


def op_call_bases(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    traces, peaks, meta = _trace_input(value, options)
    n = len(traces["A"])
    mixed_ratio = _num(options, "mixed_peak_ratio", default=0.5)
    min_height = _num(options, "min_peak_height", default=0.1)
    min_spacing = _int(options, "min_peak_spacing", 3, minimum=1, maximum=1000)
    assert mixed_ratio and min_height
    if not 0 < mixed_ratio <= 1 or not 0 < min_height < 1:
        raise ConfigurationError(
            "mixed_peak_ratio must be in (0, 1] and min_peak_height in (0, 1).",
            code="INVALID_OPTION",
        )
    detected = False
    if peaks is None:
        env = [max(traces[b][i] for b in "ACGT") for i in range(n)]
        top = max(env)
        if top <= 0:
            raise InputFormatError("trace has no positive signal.", code="INVALID_INPUT")
        peaks = []
        for i in range(1, n - 1):
            if env[i] >= min_height * top and env[i] > env[i - 1] and env[i] >= env[i + 1]:
                if peaks and i - peaks[-1] < min_spacing:
                    if env[i] > env[peaks[-1]]:
                        peaks[-1] = i
                    continue
                peaks.append(i)
        detected = True
    if any(not 0 <= p < n for p in peaks):
        raise InputFormatError("a peak location lies outside the trace.", code="INVALID_INPUT")
    bases, ratios, mixed = [], [], []
    for index, p in enumerate(peaks):
        heights = sorted(((traces[b][p], b) for b in "ACGT"), key=lambda hb: (-hb[0], hb[1]))
        (h1, b1), (h2, b2) = heights[0], heights[1]
        if h1 <= 0:
            bases.append("N")
            ratios.append(None)
            continue
        ratio = h2 / h1
        ratios.append(round(ratio, 4))
        if ratio >= mixed_ratio:
            code = _IUPAC2[frozenset((b1, b2))]
            bases.append(code)
            mixed.append(
                {
                    "index": index,
                    "sample": p,
                    "bases": b1 + b2,
                    "code": code,
                    "secondary_to_primary": round(ratio, 4),
                }
            )
        else:
            bases.append(b1)
    seq = "".join(bases)
    result = {
        **meta,
        "samples": n,
        "peaks": len(peaks),
        "peaks_detected_automatically": detected,
        "sequence": seq,
        "secondary_to_primary_peak_ratio": ratios,
        "mixed_positions": mixed,
        "peak_locations": list(peaks)[:100_000],
        "definition": "At each peak location the base is the channel with the highest signal; if the second-highest channel is at least "
        "mixed_peak_ratio (default 0.5) of it, the IUPAC two-base code is reported instead. Without peak_locations, peaks are local maxima of the "
        "per-sample maximum across channels above min_peak_height (default 10%) of the global maximum, at least min_peak_spacing samples apart. "
        "No Phred-style quality values are estimated.",
    }
    if "instrument_basecalls" in meta:
        called = meta["instrument_basecalls"]
        same = sum(1 for a, b in zip(seq, called, strict=False) if a == b)
        result["agreement_with_instrument_basecalls"] = round(
            same / max(1, max(len(seq), len(called))), 4
        )
    return result


# --------------------------------------------------------------------------- Sanger read assembly


def _overlaps(
    a: str, b: str, min_overlap: int, min_identity: float, seed: int
) -> tuple[int, float] | None:
    """Best (overlap length, identity) with ``b`` starting inside ``a`` and running to the end of ``a`` (or ending inside it).

    Candidate start offsets come from exact seeds of length ``seed`` taken at every non-overlapping offset of ``b``; each candidate is
    then verified by an ungapped comparison of the whole overlap.
    """

    if len(b) < min_overlap:
        return None
    starts: set[int] = set()
    for offset in range(0, len(b) - seed + 1, seed):
        anchor = b[offset : offset + seed]
        i = a.find(anchor)
        while i != -1:
            if i - offset >= 0:
                starts.add(i - offset)
            i = a.find(anchor, i + 1)
    best: tuple[int, float] | None = None
    for start in starts:
        ov = min(len(a) - start, len(b))
        if ov < min_overlap:
            continue
        mism = sum(1 for x, y in zip(a[start : start + ov], b[:ov], strict=False) if x != y)
        ident = 1.0 - mism / ov
        if ident >= min_identity and (best is None or (ov, ident) > best):
            best = (ov, ident)
    return best


def op_assemble_sanger(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    items = [(n, dna(s, n, strict=False, maximum=20_000)) for n, s in _sequences(value, 200)]
    min_overlap = _int(options, "min_overlap", 20, minimum=4, maximum=5000)
    identity = _num(options, "min_identity", default=0.95)
    both = _flag(options, "both_strands", True)
    assert identity is not None
    if not 0.5 <= identity <= 1.0:
        raise ConfigurationError("min_identity must be in [0.5, 1].", code="INVALID_OPTION")
    if len(items) > 100:
        raise InputFormatError("at most 100 reads can be assembled.", code="INPUT_TOO_LARGE")
    seed = max(4, min(12, min_overlap // 2))
    # contig = {"seq": str, "reads": [{id, strand, offset, length}]}
    contigs = [
        {"seq": s, "reads": [{"id": n, "strand": "+", "offset": 0, "length": len(s)}]}
        for n, s in items
    ]

    def flip(c: dict[str, Any]) -> dict[str, Any]:
        total = len(c["seq"])
        return {
            "seq": revcomp(c["seq"]),
            "reads": [
                {
                    "id": r["id"],
                    "strand": "-" if r["strand"] == "+" else "+",
                    "offset": total - r["offset"] - r["length"],
                    "length": r["length"],
                }
                for r in c["reads"]
            ],
        }

    while len(contigs) > 1:
        best = None  # (overlap, identity, i, j, flip_j)
        for i, x in enumerate(contigs):
            for j, y in enumerate(contigs):
                if i == j:
                    continue
                for flip_y in (False, True) if both else (False,):
                    yseq = revcomp(y["seq"]) if flip_y else y["seq"]
                    found = _overlaps(x["seq"], yseq, min_overlap, identity, seed)
                    if found and (
                        best is None
                        or (found[0], found[1], -i, -j) > (best[0], best[1], -best[2], -best[3])
                    ):
                        best = (found[0], found[1], i, j, flip_y)
        if best is None:
            break
        ov, _ident, i, j, flip_y = best
        x = contigs[i]
        y = flip(contigs[j]) if flip_y else contigs[j]
        start = len(x["seq"]) - ov
        # y begins at x-offset `start` (it may end inside x when y is contained)
        merged_seq = x["seq"] + y["seq"][ov:] if len(y["seq"]) > ov else x["seq"]
        reads = list(x["reads"]) + [{**r, "offset": r["offset"] + start} for r in y["reads"]]
        contigs = [c for k, c in enumerate(contigs) if k not in (i, j)] + [
            {"seq": merged_seq, "reads": reads}
        ]
    contigs.sort(key=lambda c: (-len(c["seq"]), c["seq"]))
    out = []
    singlets = []
    for _k, c in enumerate(contigs, 1):
        if len(c["reads"]) == 1:
            singlets.append({"id": c["reads"][0]["id"], "sequence": c["seq"]})
        else:
            out.append(
                {
                    "id": f"contig{len(out) + 1}",
                    "length": len(c["seq"]),
                    "sequence": c["seq"],
                    "reads": sorted(c["reads"], key=lambda r: r["offset"]),
                }
            )
    return {
        "reads": len(items),
        "contigs": out,
        "singlets": singlets,
        "min_overlap": min_overlap,
        "min_identity": identity,
        "both_strands": both,
        "definition": "Greedy overlap-layout: repeatedly merge the pair of contigs with the longest suffix/prefix overlap (>= min_overlap bases, "
        "identity >= min_identity, reads or their reverse complements) until none remains; the consensus keeps the earlier contig's base in "
        "overlaps (no quality weighting, no gap-handling), so reads should be quality-trimmed first. Repeats longer than the overlap can mis-join reads.",
    }
