"""Summary-statistics tables for ``dnakit.gwas``: parsing, standardisation, QC plots data, loci and
meta-analysis.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

from . import _impl_stats as st
from ._impl_common import (
    dumps_json,
    load_table,
    opt_bool,
    opt_choice,
    opt_number,
    require_option,
    serialize_records,
    to_number,
    write_text_file,
)

CANONICAL = (
    "variant_id",
    "chrom",
    "pos",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    "eaf",
    "n",
    "z",
)
GENOME_WIDE = 5e-8
SUGGESTIVE = 1e-5
NEGLOG10_P_FLOOR = -math.log10(5e-324)  # 323.3, the largest finite -log10(p) of a double

_ALIASES: dict[str, tuple[str, ...]] = {
    "variant_id": (
        "variantid",
        "snp",
        "snpid",
        "rsid",
        "rs",
        "rsids",
        "variant",
        "id",
        "markername",
        "marker",
    ),
    "chrom": ("chrom", "chr", "chromosome", "chromosomenumber"),
    "pos": ("pos", "bp", "position", "basepairlocation", "bppos", "basepair", "chrpos"),
    "effect_allele": (
        "effectallele",
        "ea",
        "a1",
        "allele1",
        "alt",
        "testedallele",
        "effectalleles",
    ),
    "other_allele": (
        "otherallele",
        "oa",
        "nea",
        "a2",
        "allele2",
        "ref",
        "noneffectallele",
        "referenceallele",
    ),
    "beta": ("beta", "b", "effect", "effectsize", "logor", "betahat"),
    "odds_ratio": ("or", "oddsratio"),
    "se": ("se", "stderr", "standarderror", "sebeta", "stderror"),
    "p": ("p", "pval", "pvalue", "pvalues", "pvalmlog"),
    "eaf": ("eaf", "freq", "freq1", "af", "effectallelefrequency", "a1freq", "afreq", "allelefreq"),
    "n": ("n", "ntotal", "samplesize", "neff", "nsamples"),
    "z": ("z", "zscore", "zstat", "zstatistic"),
    "maf": ("maf", "minorallelefrequency"),
    "info": ("info", "imputationinfo", "rsq"),
    "ld_score": ("ldscore", "ld", "l2", "ldscores"),
}
_ALIAS_LOOKUP = {alias: canon for canon, aliases in _ALIASES.items() for alias in aliases}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _norm_key(key: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key).lower())


def normalize_chrom(value: object) -> str | None:
    if value is None:
        return None
    text = (
        str(int(value)) if isinstance(value, float) and value == int(value) else str(value).strip()
    )
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.upper()
    return {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(text, text) or None


def chrom_order(chrom: str) -> int:
    if chrom.isdigit():
        return int(chrom)
    return {"X": 23, "Y": 24, "XY": 25, "MT": 26}.get(chrom, 100 + sum(map(ord, chrom)))


def revcomp(allele: str) -> str:
    return allele.translate(_COMPLEMENT)[::-1]


def is_palindromic(a: str | None, b: str | None) -> bool:
    return a is not None and b is not None and len(a) == 1 and len(b) == 1 and a == revcomp(b)


def _number(record: Mapping[str, object], key: str, row: int) -> float | None:
    try:
        return to_number(record.get(key))
    except ValueError as exc:
        raise InputFormatError(
            f"row {row}: column '{key}' has the non-numeric value {record.get(key)!r}.",
            code="NON_NUMERIC_VALUE",
            hint="Numeric columns must hold numbers or NA.",
        ) from exc


def standardize_rows(
    records: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, Any]], dict[str, str], list[str]]:
    """Rename recognised columns to the canonical names and derive missing z / p / se /
    variant_id.
    """

    if not records:
        raise InputFormatError(
            "the table has no rows.", code="EMPTY_TABLE", hint="Provide at least one variant."
        )
    column_map: dict[str, str] = {}
    all_keys = list(dict.fromkeys(str(key) for rec in records for key in rec))
    for key in all_keys:
        canon = _ALIAS_LOOKUP.get(_norm_key(key))
        if canon is not None and canon not in column_map.values():
            column_map[key] = canon
    if not column_map:
        raise InputFormatError(
            "no recognisable GWAS columns found.",
            code="UNKNOWN_COLUMNS",
            hint="Use columns such as SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N.",
        )
    rows: list[dict[str, Any]] = []
    derived: set[str] = set()
    for i, rec in enumerate(records, start=1):
        row: dict[str, Any] = {name: None for name in CANONICAL}
        odds_ratio = None
        extra_values: dict[str, object] = {}
        for key, canon in column_map.items():
            raw = rec.get(key)
            if canon in ("variant_id",):
                row[canon] = None if raw is None or str(raw).strip() == "" else str(raw)
            elif canon == "chrom":
                row[canon] = normalize_chrom(raw)
            elif canon == "pos":
                value = _number(rec, key, i)
                row[canon] = None if value is None else int(value)
            elif canon in ("effect_allele", "other_allele"):
                row[canon] = (
                    None if raw is None or str(raw).strip() == "" else str(raw).strip().upper()
                )
            elif canon == "odds_ratio":
                odds_ratio = _number(rec, key, i)
            elif canon in ("beta", "se", "p", "eaf", "n", "z"):
                row[canon] = _number(rec, key, i)
                if canon == "n" and row[canon] is not None and row[canon] == int(row[canon]):
                    row[canon] = int(row[canon])
            elif canon in ("maf", "info", "ld_score"):
                extra_values[canon] = _number(rec, key, i)
        if row["beta"] is None and odds_ratio is not None:
            if odds_ratio <= 0:
                raise InputFormatError(
                    f"row {i}: odds ratio must be positive.",
                    code="INVALID_VALUE",
                    hint="Odds ratios are converted with beta = ln(OR).",
                )
            row["beta"] = math.log(odds_ratio)
            derived.add("beta")
        if row["z"] is None and row["beta"] is not None and row["se"]:
            row["z"] = row["beta"] / row["se"]
            derived.add("z")
        if row["p"] is None and row["z"] is not None:
            row["p"] = st.two_sided_p_from_z(row["z"])
            derived.add("p")
        if (
            row["se"] is None
            and row["beta"] is not None
            and row["p"] is not None
            and 0 < row["p"] < 1
            and row["beta"] != 0
        ):
            row["se"] = abs(row["beta"]) / st.z_from_two_sided_p(row["p"])
            derived.add("se")
        if row["z"] is None and row["beta"] is not None and row["se"]:
            row["z"] = row["beta"] / row["se"]
            derived.add("z")
        if row["variant_id"] is None and row["chrom"] is not None and row["pos"] is not None:
            alleles = [a for a in (row["other_allele"], row["effect_allele"]) if a]
            row["variant_id"] = ":".join([row["chrom"], str(row["pos"]), *alleles])
            derived.add("variant_id")
        row.update({k: v for k, v in extra_values.items() if v is not None})
        rows.append(row)
    warnings = []
    ignored = [k for k in all_keys if k not in column_map]
    if ignored:
        warnings.append("ignored unrecognised columns: " + ", ".join(ignored))
    return rows, column_map, warnings


def load_sumstats(
    source: object, opts: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, str], list[str]]:
    max_rows = opt_number(opts, "max_rows", 500_000, minimum=1, integer=True)
    records = load_table(
        source, what="summary statistics", max_rows=max_rows, delimiter=opts.get("delimiter")
    )
    return standardize_rows(records)


def _need(rows: Sequence[Mapping[str, Any]], *fields: str) -> None:
    missing = [f for f in fields if not any(r.get(f) is not None for r in rows)]
    if missing:
        raise InputFormatError(
            "required column(s) missing or empty: " + ", ".join(missing),
            code="MISSING_COLUMNS",
            hint="Add the columns (aliases such as P/PVAL, BETA, SE, BP/POS are recognised).",
        )


def neglog10p(p: float | None) -> float | None:
    if p is None:
        return None
    if p <= 0:
        return NEGLOG10_P_FLOOR  # p underflowed to 0: report the floating-point floor, not infinity
    return -math.log10(p)


def chi2_from_p(p: float) -> float:
    z = st.z_from_two_sided_p(p)
    return z * z


def lambda_gc_of(rows: Sequence[Mapping[str, Any]]) -> float | None:
    chi = [chi2_from_p(r["p"]) for r in rows if r.get("p") is not None and 0 < r["p"] <= 1]
    if not chi:
        return None
    return st.median(chi) / st.MEDIAN_CHI2_1DF


# -------------------------------------------------------------------- read / standardize / validate
def read_table(source: object, opts: Mapping[str, Any]) -> dict[str, Any]:
    max_rows = opt_number(opts, "max_rows", 500_000, minimum=1, integer=True)
    records = load_table(
        source, what="GWAS table", max_rows=max_rows, delimiter=opts.get("delimiter")
    )
    if opt_bool(opts, "standardize", False):
        rows, cmap, _ = standardize_rows(records)
        return {
            "columns": list(CANONICAL),
            "column_map": cmap,
            "n_records": len(rows),
            "records": rows,
        }
    columns = list(records[0].keys()) if records else []
    return {"columns": columns, "n_records": len(records), "records": records}


def standardize(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, cmap, warnings = load_sumstats(source, opts)
    return {
        "records": rows,
        "column_map": cmap,
        "columns": list(CANONICAL),
        "n_records": len(rows),
    }, warnings


def validate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, cmap, notes = load_sumstats(source, opts)
    errors: list[str] = []
    warnings: list[str] = list(notes)
    n = len(rows)

    def count(pred: Any) -> int:
        return sum(1 for r in rows if pred(r))

    has_id = any(r["variant_id"] for r in rows)
    has_pos = any(r["chrom"] is not None and r["pos"] is not None for r in rows)
    if not has_id and not has_pos:
        errors.append("no variant identifier: need variant_id (SNP) or chrom + pos")
    if not (
        any(r["p"] is not None for r in rows)
        or any(r["beta"] is not None and r["se"] for r in rows)
    ):
        errors.append("no association statistic: need p, z or beta with se")
    bad_p = count(lambda r: r["p"] is not None and not 0.0 <= r["p"] <= 1.0)
    if bad_p:
        errors.append(f"{bad_p} p-values outside [0, 1]")
    bad_se = count(lambda r: r["se"] is not None and not r["se"] > 0)
    if bad_se:
        errors.append(f"{bad_se} standard errors that are not > 0")
    bad_eaf = count(lambda r: r["eaf"] is not None and not 0.0 <= r["eaf"] <= 1.0)
    if bad_eaf:
        errors.append(f"{bad_eaf} effect-allele frequencies outside [0, 1]")
    bad_n = count(lambda r: r["n"] is not None and r["n"] <= 0)
    if bad_n:
        errors.append(f"{bad_n} sample sizes that are not > 0")
    nonfinite = count(
        lambda r: any(r[k] is not None and not math.isfinite(r[k]) for k in ("beta", "se", "z"))
    )
    if nonfinite:
        errors.append(f"{nonfinite} non-finite beta/se/z values")
    acgt = re.compile(r"^[ACGT]+$")
    bad_alleles = count(
        lambda r: any(
            a is not None and not acgt.match(a) for a in (r["effect_allele"], r["other_allele"])
        )
    )
    if bad_alleles:
        warnings.append(f"{bad_alleles} variants have alleles other than A/C/G/T strings")
    same_alleles = count(
        lambda r: r["effect_allele"] is not None and r["effect_allele"] == r["other_allele"]
    )
    if same_alleles:
        errors.append(f"{same_alleles} variants have identical effect and other alleles")
    ids = [r["variant_id"] for r in rows if r["variant_id"]]
    duplicated = len(ids) - len(set(ids))
    if duplicated:
        warnings.append(f"{duplicated} duplicated variant identifiers")
    if not any(r["effect_allele"] for r in rows):
        warnings.append("no effect allele column; signed effects cannot be harmonised")
    if not any(r["eaf"] is not None for r in rows):
        warnings.append("no effect-allele frequency column")
    return {
        "valid": not errors,
        "n_records": n,
        "errors": errors,
        "warnings": warnings,
        "column_map": cmap,
    }, []


# --------------------------------------------------------------------------- filter
def _parse_region(region: object) -> tuple[str, int, int]:
    if isinstance(region, Mapping):
        chrom = normalize_chrom(region.get("chrom"))
        start, end = region.get("start"), region.get("end")
    elif isinstance(region, str):
        match = re.fullmatch(r"\s*([\w.]+)\s*:\s*(\d+)\s*-\s*(\d+)\s*", region)
        if not match:
            raise ConfigurationError(
                f"region {region!r} is not 'chrom:start-end'.",
                code="INVALID_OPTION",
                hint="Use e.g. '7:1000000-2000000'.",
            )
        chrom, start, end = (
            normalize_chrom(match.group(1)),
            int(match.group(2)),
            int(match.group(3)),
        )
    else:
        raise ConfigurationError(
            "region must be 'chrom:start-end' or a mapping with chrom/start/end.",
            code="INVALID_OPTION",
            hint="Use e.g. '7:1000000-2000000'.",
        )
    if (
        chrom is None
        or not isinstance(start, (int, float))
        or not isinstance(end, (int, float))
        or start > end
    ):
        raise ConfigurationError(
            "region needs chrom and start <= end.",
            code="INVALID_OPTION",
            hint="Use e.g. '7:1000000-2000000'.",
        )
    return chrom, int(start), int(end)


def filter_rows(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    max_p = opt_number(opts, "max_p", None, minimum=0, maximum=1)
    min_maf = opt_number(opts, "min_maf", None, minimum=0, maximum=0.5)
    min_n = opt_number(opts, "min_n", None, minimum=0)
    min_info = opt_number(opts, "min_info", None, minimum=0)
    snp_only = opt_bool(opts, "snp_only", False)
    drop_palindromic = opt_bool(opts, "remove_palindromic", False)
    require_p = opt_bool(opts, "require_p", False)
    chroms = opts.get("chromosomes")
    chrom_set = None
    if chroms is not None:
        if isinstance(chroms, str) or not isinstance(chroms, Sequence):
            raise ConfigurationError(
                "chromosomes must be a list such as ['1', 'X'].",
                code="INVALID_OPTION",
                hint="Pass a list of chromosome names.",
            )
        chrom_set = {normalize_chrom(c) for c in chroms}
    region = _parse_region(opts["region"]) if opts.get("region") is not None else None
    removed: dict[str, int] = defaultdict(int)
    kept = []
    for r in rows:
        reason = None
        if require_p and r["p"] is None:
            reason = "missing_p"
        elif max_p is not None and (r["p"] is None or r["p"] > max_p):
            reason = "p_above_threshold"
        elif min_maf is not None and (r["eaf"] is None or min(r["eaf"], 1 - r["eaf"]) < min_maf):
            reason = "maf_below_threshold"
        elif min_n is not None and (r["n"] is None or r["n"] < min_n):
            reason = "n_below_threshold"
        elif min_info is not None and (r.get("info") is None or r["info"] < min_info):
            reason = "info_below_threshold"
        elif chrom_set is not None and r["chrom"] not in chrom_set:
            reason = "chromosome_excluded"
        elif region is not None and not (
            r["chrom"] == region[0] and r["pos"] is not None and region[1] <= r["pos"] <= region[2]
        ):
            reason = "outside_region"
        elif snp_only and not (
            r["effect_allele"]
            and r["other_allele"]
            and len(r["effect_allele"]) == 1
            and len(r["other_allele"]) == 1
        ):
            reason = "not_a_snp"
        elif drop_palindromic and is_palindromic(r["effect_allele"], r["other_allele"]):
            reason = "palindromic"
        if reason:
            removed[reason] += 1
        else:
            kept.append(r)
    return {
        "records": kept,
        "n_input": len(rows),
        "n_kept": len(kept),
        "removed": dict(removed),
    }, notes


# --------------------------------------------------------------------------- harmonize
def _align_alleles(
    a1: str,
    a2: str,
    b1: str,
    b2: str,
    e_t: float | None,
    e_r: float | None,
    palindromes: str,
    tolerance: float,
) -> tuple[str, bool] | str:
    """Return (action, flip_sign) that expresses (a1 effect, a2 other) in the reference's (b1,
    b2), or a drop reason.
    """

    if is_palindromic(a1, a2):
        if (a1, a2) not in ((b1, b2), (b2, b1)):
            return "allele_mismatch"
        if palindromes == "assume_forward_strand":
            return ("none", False) if (a1, a2) == (b1, b2) else ("swap", True)
        if e_t is None or e_r is None or abs(e_t - 0.5) <= tolerance or abs(e_r - 0.5) <= tolerance:
            return "ambiguous_palindrome"
        same_side = (e_t > 0.5) == (e_r > 0.5)
        if (a1, a2) == (b1, b2):
            return ("none", False) if same_side else ("strand_flip_swap", True)
        return ("strand_flip", False) if same_side else ("swap", True)
    c1, c2 = revcomp(a1), revcomp(a2)
    if (a1, a2) == (b1, b2):
        return "none", False
    if (a1, a2) == (b2, b1):
        return "swap", True
    if (c1, c2) == (b1, b2):
        return "strand_flip", False
    if (c1, c2) == (b2, b1):
        return "strand_flip_swap", True
    return "allele_mismatch"


def harmonize(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    if isinstance(source, Mapping) and "target" in source and "reference" in source:
        target_src, ref_src = source["target"], source["reference"]
    else:
        target_src, ref_src = (
            source,
            require_option(
                opts, "reference", "harmonisation aligns the input to a reference table"
            ),
        )
    target, _, notes = load_sumstats(target_src, opts)
    reference, _, _ = load_sumstats(ref_src, opts)
    match_on = opt_choice(opts, "match_on", "auto", ("auto", "id", "position"))
    palindromes = opt_choice(
        opts, "palindromes", "drop_ambiguous", ("drop_ambiguous", "assume_forward_strand")
    )
    tolerance = opt_number(opts, "palindrome_eaf_tolerance", 0.08, minimum=0, maximum=0.5)
    _need(target, "effect_allele", "other_allele")
    _need(reference, "effect_allele", "other_allele")

    def key(r: Mapping[str, Any], mode: str) -> str | None:
        if mode == "position":
            return None if r["chrom"] is None or r["pos"] is None else f"{r['chrom']}:{r['pos']}"
        return r["variant_id"]

    mode = match_on
    if mode == "auto":
        mode = (
            "position"
            if all(r["chrom"] is not None and r["pos"] is not None for r in target + reference)
            else "id"
        )
    index: dict[str, dict[str, Any]] = {}
    for r in reference:
        k = key(r, mode)
        if k is not None and k not in index:
            index[k] = r
    out: list[dict[str, Any]] = []
    counts: dict[str, int] = defaultdict(int)
    for r in target:
        k = key(r, mode)
        ref = index.get(k) if k is not None else None
        if ref is None:
            counts["not_in_reference"] += 1
            continue
        a1, a2 = r["effect_allele"], r["other_allele"]
        b1, b2 = ref["effect_allele"], ref["other_allele"]
        if a1 is None or a2 is None or b1 is None or b2 is None:
            counts["missing_alleles"] += 1
            continue
        verdict = _align_alleles(a1, a2, b1, b2, r["eaf"], ref["eaf"], palindromes, tolerance)
        if isinstance(verdict, str):
            counts[verdict] += 1
            continue
        action, flip_sign = verdict
        new = dict(r)
        new["effect_allele"], new["other_allele"] = b1, b2
        if flip_sign:
            for field in ("beta", "z"):
                if new[field] is not None:
                    new[field] = -new[field]
            if new["eaf"] is not None:
                new["eaf"] = 1.0 - new["eaf"]
        new["harmonize_action"] = action
        counts[action] += 1
        out.append(new)
    return {
        "records": out,
        "n_target": len(target),
        "n_harmonized": len(out),
        "actions": dict(counts),
        "matched_on": mode,
    }, notes


# --------------------------------------------------------------------------- write / export
def _table_columns(rows: Sequence[Mapping[str, Any]], opts: Mapping[str, Any]) -> list[str]:
    columns = opts.get("columns")
    if columns is not None:
        if (
            isinstance(columns, str)
            or not isinstance(columns, Sequence)
            or not all(isinstance(c, str) for c in columns)
        ):
            raise ConfigurationError(
                "columns must be a list of column names.",
                code="INVALID_OPTION",
                hint="Pass e.g. ['variant_id', 'beta', 'se', 'p'].",
            )
        unknown = [c for c in columns if c not in rows[0]]
        if unknown:
            raise ConfigurationError(
                f"unknown column(s): {unknown}",
                code="INVALID_OPTION",
                hint=f"Available columns: {list(rows[0])}",
            )
        return list(columns)
    return [c for c in rows[0] if any(r.get(c) is not None for r in rows)]


def write_table(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    columns = _table_columns(rows, opts)
    path = opts.get("path", opts.get("output"))
    fmt = opt_choice(
        opts,
        "format",
        "csv" if path and str(path).endswith((".csv", ".csv.gz")) else "tsv",
        ("tsv", "csv"),
    )
    text = serialize_records(rows, columns, "\t" if fmt == "tsv" else ",")
    written = None
    if path is not None:
        written = write_text_file(path, text, overwrite=opt_bool(opts, "overwrite", False))
    return {
        "format": fmt,
        "columns": columns,
        "n_records": len(rows),
        "text": text,
        "path": written,
    }, notes


def export_table(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    fmt = opt_choice(opts, "format", "tsv", ("tsv", "csv", "json", "jsonl", "bed", "ldsc", "metal"))
    if fmt in ("tsv", "csv"):
        columns = _table_columns(rows, opts)
        text = serialize_records(rows, columns, "\t" if fmt == "tsv" else ",")
    elif fmt in ("json", "jsonl"):
        text = dumps_json(rows) if fmt == "json" else "\n".join(dumps_json(r) for r in rows) + "\n"
    elif fmt == "bed":
        _need(rows, "chrom", "pos")
        lines = []
        for r in rows:
            if r["chrom"] is None or r["pos"] is None:
                continue
            score = neglog10p(r["p"])
            lines.append(
                "\t".join(
                    [
                        "chr" + r["chrom"],
                        str(r["pos"] - 1),
                        str(r["pos"]),
                        r["variant_id"] or ".",
                        "." if score is None else repr(score),
                    ]
                )
            )
        text = "\n".join(lines) + "\n"
    elif fmt == "ldsc":
        _need(rows, "variant_id", "effect_allele", "other_allele", "z")
        sel = [
            {
                "SNP": r["variant_id"],
                "A1": r["effect_allele"],
                "A2": r["other_allele"],
                "N": r["n"],
                "Z": r["z"],
            }
            for r in rows
            if r["z"] is not None
        ]
        text = serialize_records(sel, ["SNP", "A1", "A2", "N", "Z"], "\t")
    else:  # metal
        _need(rows, "variant_id", "effect_allele", "other_allele", "beta", "se")
        sel = [
            {
                "MarkerName": r["variant_id"],
                "Allele1": r["effect_allele"],
                "Allele2": r["other_allele"],
                "Freq1": r["eaf"],
                "Effect": r["beta"],
                "StdErr": r["se"],
                "P-value": r["p"],
            }
            for r in rows
            if r["beta"] is not None and r["se"]
        ]
        text = serialize_records(
            sel, ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"], "\t"
        )
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, text, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"format": fmt, "n_records": len(rows), "text": text, "path": written}, notes


# ---------------------------------------------------------------------- calculate / genomic control
def calculate(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    statistic = opt_choice(
        opts,
        "statistic",
        "summary",
        ("summary", "lambda_gc", "z_scores", "neglog10p", "maf", "n_significant", "adjusted_p"),
    )
    threshold = opt_number(opts, "significance", GENOME_WIDE, minimum=0, maximum=1, exclusive=False)
    if statistic == "lambda_gc":
        _need(rows, "p")
        return {
            "statistic": "lambda_gc",
            "value": lambda_gc_of(rows),
            "median_chi2": lambda_gc_of(rows) * st.MEDIAN_CHI2_1DF
            if lambda_gc_of(rows) is not None
            else None,
        }, notes
    if statistic == "z_scores":
        _need(rows, "z")
        return {
            "statistic": "z_scores",
            "values": [{"variant_id": r["variant_id"], "z": r["z"]} for r in rows],
        }, notes
    if statistic == "neglog10p":
        _need(rows, "p")
        return {
            "statistic": "neglog10p",
            "values": [
                {"variant_id": r["variant_id"], "neglog10p": neglog10p(r["p"])} for r in rows
            ],
        }, notes
    if statistic == "maf":
        _need(rows, "eaf")
        return {
            "statistic": "maf",
            "values": [
                {
                    "variant_id": r["variant_id"],
                    "maf": None if r["eaf"] is None else min(r["eaf"], 1 - r["eaf"]),
                }
                for r in rows
            ],
        }, notes
    _need(rows, "p")
    valid = [r for r in rows if r["p"] is not None]
    if statistic == "adjusted_p":
        bonf, bh = (
            st.bonferroni([r["p"] for r in valid]),
            st.benjamini_hochberg([r["p"] for r in valid]),
        )
        return {
            "statistic": "adjusted_p",
            "n_tests": len(valid),
            "values": [
                {"variant_id": r["variant_id"], "p": r["p"], "p_bonferroni": a, "q_bh": b}
                for r, a, b in zip(valid, bonf, bh, strict=False)
            ],
        }, notes
    significant = [r for r in valid if r["p"] <= threshold]
    if statistic == "n_significant":
        return {
            "statistic": "n_significant",
            "threshold": threshold,
            "value": len(significant),
        }, notes
    lead = min(valid, key=lambda r: r["p"])
    chi = [chi2_from_p(r["p"]) for r in valid if 0 < r["p"] <= 1]
    return {
        "statistic": "summary",
        "n_variants": len(rows),
        "n_with_p": len(valid),
        "min_p": lead["p"],
        "lead_variant": lead["variant_id"],
        "lambda_gc": lambda_gc_of(rows),
        "mean_chi2": st.mean(chi) if chi else None,
        "median_chi2": st.median(chi) if chi else None,
        "significance_threshold": threshold,
        "n_significant": len(significant),
        "n_chromosomes": len({r["chrom"] for r in rows if r["chrom"]}),
    }, notes


def apply_genomic_control(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "p")
    lam = lambda_gc_of(rows)
    given = opt_number(opts, "lambda_gc", None, minimum=0, exclusive=True)
    lam_est = given if given is not None else lam
    only_inflated = opt_bool(opts, "only_inflated", True)
    used = max(1.0, lam_est) if only_inflated else lam_est
    out = []
    for r in rows:
        new = dict(r)
        if r["p"] is not None:
            chi = chi2_from_p(r["p"]) if r["p"] > 0 else math.inf
            new["p_gc"] = st.chi2_sf(chi / used, 1) if math.isfinite(chi) else 0.0
        else:
            new["p_gc"] = None
        new["se_gc"] = None if r["se"] is None else r["se"] * math.sqrt(used)
        new["z_gc"] = None if r["z"] is None else r["z"] / math.sqrt(used)
        out.append(new)
    return {"lambda_gc": lam, "lambda_used": used, "records": out, "n_records": len(out)}, notes


# --------------------------------------------------------------------------- plot data
def _downsample(values: list[Any], max_points: int, key: Any) -> list[Any]:
    if len(values) <= max_points:
        return values
    ordered = sorted(values, key=key, reverse=True)
    head = ordered[: max_points // 2]
    rest = ordered[max_points // 2 :]
    step = len(rest) / (max_points - len(head))
    tail = [rest[int(i * step)] for i in range(max_points - len(head))]
    return head + tail


def compute_qq_data(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "p")
    max_points = opt_number(opts, "max_points", 50_000, minimum=10, integer=True)
    ps = sorted(r["p"] for r in rows if r["p"] is not None and 0 <= r["p"] <= 1)
    if not ps:
        raise InputFormatError(
            "no valid p-values in [0, 1].", code="NO_PVALUES", hint="Provide a p column."
        )
    n = len(ps)
    points = [
        {"rank": i + 1, "expected": -math.log10((i + 0.5) / n), "observed": neglog10p(p)}
        for i, p in enumerate(ps)
    ]
    kept = (
        _downsample(points, max_points, key=lambda pt: pt["observed"]) if n > max_points else points
    )
    kept = sorted(kept, key=lambda pt: pt["rank"], reverse=True)
    return {
        "n": n,
        "points": kept,
        "lambda_gc": lambda_gc_of(rows),
        "downsampled": len(kept) < n,
        "expected_definition": "-log10((rank - 0.5) / n) with rank 1 = smallest p",
    }, notes


def compute_manhattan_data(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    min_nlp = opt_number(opts, "min_neglog10p", 0.0, minimum=0)
    max_points = opt_number(opts, "max_points", 200_000, minimum=10, integer=True)
    usable = [r for r in rows if r["chrom"] and r["pos"] is not None and r["p"] is not None]
    lengths: dict[str, int] = {}
    for r in usable:
        lengths[r["chrom"]] = max(lengths.get(r["chrom"], 0), r["pos"])
    offsets: dict[str, int] = {}
    chroms = []
    running = 0
    for chrom in sorted(lengths, key=chrom_order):
        offsets[chrom] = running
        chroms.append(
            {
                "chrom": chrom,
                "offset": running,
                "length": lengths[chrom],
                "center": running + lengths[chrom] / 2.0,
            }
        )
        running += lengths[chrom]
    points = []
    for r in usable:
        v = neglog10p(r["p"])
        if v is not None and v >= min_nlp:
            points.append(
                {
                    "variant_id": r["variant_id"],
                    "chrom": r["chrom"],
                    "pos": r["pos"],
                    "x": offsets[r["chrom"]] + r["pos"],
                    "neglog10p": v,
                }
            )
    truncated = len(points) > max_points
    if truncated:
        points = _downsample(points, max_points, key=lambda pt: pt["neglog10p"])
    points.sort(key=lambda pt: pt["x"])
    return {
        "points": points,
        "chromosomes": chroms,
        "total_length": running,
        "truncated": truncated,
        "thresholds": {
            "genome_wide": {"p": GENOME_WIDE, "neglog10p": -math.log10(GENOME_WIDE)},
            "suggestive": {"p": SUGGESTIVE, "neglog10p": -math.log10(SUGGESTIVE)},
        },
        "x_definition": (
            "cumulative base-pair position: chromosome offset (sum of preceding chromosomes' "
            "maximum observed positions) + pos"
        ),
    }, notes


def _ld_r2_to_lead(opts: Mapping[str, Any], lead_id: str, ids: Sequence[str]) -> dict[str, float]:
    from ._impl_models import correlation_lookup  # local import to avoid a cycle

    return correlation_lookup(opts, lead_id, ids)


def compute_locus_plot_data(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    flank = opt_number(opts, "flank", 250_000, minimum=0, integer=True)
    usable = [r for r in rows if r["chrom"] and r["pos"] is not None and r["p"] is not None]
    lead_id = opts.get("lead")
    if lead_id is None:
        lead = min(usable, key=lambda r: r["p"])
    else:
        matches = [r for r in usable if r["variant_id"] == lead_id]
        if not matches:
            raise ConfigurationError(
                f"lead variant {lead_id!r} not found in the table.",
                code="UNKNOWN_LEAD",
                hint="Use a variant_id present in the input.",
            )
        lead = matches[0]
    window = [
        r for r in usable if r["chrom"] == lead["chrom"] and abs(r["pos"] - lead["pos"]) <= flank
    ]
    r2_map: dict[str, float] = {}
    if opts.get("ld") is not None or opts.get("genotypes") is not None:
        r2_map = _ld_r2_to_lead(opts, lead["variant_id"], [r["variant_id"] for r in window])
    points = []
    for r in sorted(window, key=lambda r: r["pos"]):
        r2 = r2_map.get(r["variant_id"])
        if r["variant_id"] == lead["variant_id"]:
            r2, band = 1.0, "lead"
        elif r2 is None:
            band = "unknown"
        elif r2 >= 0.8:
            band = "0.8-1.0"
        elif r2 >= 0.6:
            band = "0.6-0.8"
        elif r2 >= 0.4:
            band = "0.4-0.6"
        elif r2 >= 0.2:
            band = "0.2-0.4"
        else:
            band = "0.0-0.2"
        points.append(
            {
                "variant_id": r["variant_id"],
                "pos": r["pos"],
                "neglog10p": neglog10p(r["p"]),
                "r2": r2,
                "ld_band": band,
            }
        )
    return {
        "lead": {
            "variant_id": lead["variant_id"],
            "chrom": lead["chrom"],
            "pos": lead["pos"],
            "p": lead["p"],
        },
        "region": {
            "chrom": lead["chrom"],
            "start": max(1, lead["pos"] - flank),
            "end": lead["pos"] + flank,
        },
        "points": points,
    }, notes


# --------------------------------------------------------------------------- lead SNPs, loci
def _significant(rows: Sequence[Mapping[str, Any]], threshold: float) -> list[dict[str, Any]]:
    return [
        dict(r)
        for r in rows
        if r["p"] is not None and r["p"] <= threshold and r["chrom"] and r["pos"] is not None
    ]


def leads_from_rows(
    rows: Sequence[Mapping[str, Any]], threshold: float, window: int, max_leads: int | None = None
) -> list[dict[str, Any]]:
    candidates = sorted(
        _significant(rows, threshold), key=lambda r: (r["p"], chrom_order(r["chrom"]), r["pos"])
    )
    leads: list[dict[str, Any]] = []
    used: set[int] = set()
    for i, cand in enumerate(candidates):
        if i in used:
            continue
        absorbed = 0
        for j, other in enumerate(candidates):
            if (
                j != i
                and j not in used
                and other["chrom"] == cand["chrom"]
                and abs(other["pos"] - cand["pos"]) <= window
            ):
                used.add(j)
                absorbed += 1
        used.add(i)
        lead = {
            k: cand[k]
            for k in (
                "variant_id",
                "chrom",
                "pos",
                "effect_allele",
                "other_allele",
                "beta",
                "se",
                "p",
                "eaf",
            )
        }
        lead["neglog10p"] = neglog10p(cand["p"])
        lead["n_nearby_significant"] = absorbed
        leads.append(lead)
        if max_leads is not None and len(leads) >= max_leads:
            break
    return leads


def select_lead_snps(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    threshold = opt_number(opts, "p_threshold", GENOME_WIDE, minimum=0, maximum=1)
    window = opt_number(opts, "window", 500_000, minimum=0, integer=True)
    max_leads = opt_number(opts, "max_leads", None, minimum=1, integer=True)
    leads = leads_from_rows(rows, threshold, window, max_leads)
    leads.sort(key=lambda r: (chrom_order(r["chrom"]), r["pos"]))
    return {
        "leads": leads,
        "n_leads": len(leads),
        "p_threshold": threshold,
        "window": window,
    }, notes


def merge_nearby_significant_hits(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    threshold = opt_number(opts, "p_threshold", GENOME_WIDE, minimum=0, maximum=1)
    distance = opt_number(opts, "distance", 500_000, minimum=0, integer=True)
    hits = sorted(_significant(rows, threshold), key=lambda r: (chrom_order(r["chrom"]), r["pos"]))
    groups: list[list[dict[str, Any]]] = []
    for hit in hits:
        if (
            groups
            and groups[-1][-1]["chrom"] == hit["chrom"]
            and hit["pos"] - groups[-1][-1]["pos"] <= distance
        ):
            groups[-1].append(hit)
        else:
            groups.append([hit])
    merged = []
    for g in groups:
        lead = min(g, key=lambda r: r["p"])
        merged.append(
            {
                "chrom": g[0]["chrom"],
                "start": g[0]["pos"],
                "end": g[-1]["pos"],
                "n_hits": len(g),
                "lead_variant": lead["variant_id"],
                "lead_p": lead["p"],
            }
        )
    return {
        "groups": merged,
        "n_groups": len(merged),
        "n_hits": len(hits),
        "distance": distance,
        "rule": (
            "single-linkage: consecutive significant variants closer than or equal to `distance` "
            "bp are merged"
        ),
    }, notes


def define_association_loci(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    threshold = opt_number(opts, "p_threshold", GENOME_WIDE, minimum=0, maximum=1)
    window = opt_number(opts, "window", 500_000, minimum=0, integer=True)
    flank = opt_number(opts, "flank", 500_000, minimum=0, integer=True)
    leads = leads_from_rows(rows, threshold, window)
    leads.sort(key=lambda r: (chrom_order(r["chrom"]), r["pos"]))
    loci: list[dict[str, Any]] = []
    for lead in leads:
        start, end = max(1, lead["pos"] - flank), lead["pos"] + flank
        if loci and loci[-1]["chrom"] == lead["chrom"] and start <= loci[-1]["end"]:
            locus = loci[-1]
            locus["end"] = max(locus["end"], end)
            locus["leads"].append(lead["variant_id"])
            if lead["p"] < locus["lead_p"]:
                locus["lead_variant"], locus["lead_p"] = lead["variant_id"], lead["p"]
        else:
            loci.append(
                {
                    "chrom": lead["chrom"],
                    "start": start,
                    "end": end,
                    "leads": [lead["variant_id"]],
                    "lead_variant": lead["variant_id"],
                    "lead_p": lead["p"],
                }
            )
    for locus in loci:
        inside = [
            r
            for r in rows
            if r["chrom"] == locus["chrom"]
            and r["pos"] is not None
            and locus["start"] <= r["pos"] <= locus["end"]
        ]
        locus["n_variants"] = len(inside)
        locus["n_significant"] = sum(
            1 for r in inside if r["p"] is not None and r["p"] <= threshold
        )
    return {
        "loci": loci,
        "n_loci": len(loci),
        "locus_definition": (
            "lead variant +/- flank bp (clamped at 1), overlapping loci on one chromosome merged"
        ),
    }, notes


def annotate_lead_snps(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "chrom", "pos", "p")
    genes_raw = require_option(
        opts,
        "genes",
        "annotation needs a gene table with name/chrom/start/end (no database is queried)",
    )
    genes = []
    for g in load_table(genes_raw, what="gene table"):
        name = g.get("name", g.get("gene", g.get("symbol")))
        chrom = normalize_chrom(g.get("chrom", g.get("chr")))
        try:
            start, end = to_number(g.get("start")), to_number(g.get("end"))
        except ValueError:
            start = end = None
        if name is None or chrom is None or start is None or end is None or start > end:
            raise InputFormatError(
                f"invalid gene record {dict(g)!r}.",
                code="INVALID_GENE_TABLE",
                hint="Each gene needs name, chrom, start <= end (and optionally strand).",
            )
        genes.append(
            {
                "name": str(name),
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": g.get("strand"),
            }
        )
    threshold = opt_number(opts, "p_threshold", GENOME_WIDE, minimum=0, maximum=1)
    window = opt_number(opts, "window", 500_000, minimum=0, integer=True)
    max_distance = opt_number(opts, "max_distance", 1_000_000, minimum=0, integer=True)
    flank = opt_number(opts, "flank", 100_000, minimum=0, integer=True)
    lead_ids = opts.get("leads")
    if lead_ids is not None:
        by_id = {r["variant_id"]: r for r in rows}
        missing = [i for i in lead_ids if i not in by_id]
        if missing:
            raise ConfigurationError(
                f"lead variants not in the table: {missing}",
                code="UNKNOWN_LEAD",
                hint="Use variant_id values present in the input.",
            )
        leads = [
            {
                k: by_id[i][k]
                for k in (
                    "variant_id",
                    "chrom",
                    "pos",
                    "effect_allele",
                    "other_allele",
                    "beta",
                    "se",
                    "p",
                    "eaf",
                )
            }
            for i in lead_ids
        ]
    else:
        leads = leads_from_rows(rows, threshold, window)
    out = []
    for lead in leads:

        def distance(g: Mapping[str, Any], pos: int = lead["pos"]) -> int:
            if g["start"] <= pos <= g["end"]:
                return 0
            return min(abs(pos - g["start"]), abs(pos - g["end"]))

        same = [g for g in genes if g["chrom"] == lead["chrom"]]
        ranked = sorted(same, key=lambda g: (distance(g), g["name"]))
        entry = dict(lead)
        entry["overlapping_genes"] = [g["name"] for g in same if distance(g) == 0]
        near = ranked[0] if ranked and distance(ranked[0]) <= max_distance else None
        entry["nearest_gene"] = None if near is None else near["name"]
        entry["distance_to_nearest_gene"] = None if near is None else distance(near)
        if near is None:
            entry["relation"] = "intergenic_no_gene_within_max_distance"
        elif distance(near) == 0:
            entry["relation"] = "within_gene"
        else:
            upstream_of_start = lead["pos"] < near["start"]
            if near["strand"] in ("+", "-", 1, -1, "1", "-1"):
                minus = near["strand"] in ("-", -1, "-1")
                entry["relation"] = "upstream" if upstream_of_start != minus else "downstream"
            else:
                entry["relation"] = "before_gene_start" if upstream_of_start else "after_gene_end"
        entry["genes_within_flank"] = [g["name"] for g in ranked if distance(g) <= flank]
        out.append(entry)
    return {
        "annotations": out,
        "n_leads": len(out),
        "method": "nearest gene by interval distance from the supplied gene table",
    }, notes


# --------------------------------------------------------------------------- fine mapping
def build_credible_set(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    rows, _, notes = load_sumstats(source, opts)
    _need(rows, "beta", "se")
    prior_sd = opt_number(opts, "prior_sd", 0.15, minimum=0, exclusive=True)
    coverage = opt_number(opts, "coverage", 0.95, minimum=0, maximum=1, exclusive=True)
    if opts.get("region") is not None:
        region = _parse_region(opts["region"])
        rows = [
            r
            for r in rows
            if r["chrom"] == region[0]
            and r["pos"] is not None
            and region[1] <= r["pos"] <= region[2]
        ]
    usable = [r for r in rows if r["beta"] is not None and r["se"] and r["se"] > 0]
    if len(usable) < 2:
        raise InputFormatError(
            "need at least two variants with beta and se > 0 to build a credible set.",
            code="TOO_FEW_VARIANTS",
            hint="Provide the summary statistics of all variants in one locus.",
        )
    w = prior_sd**2
    logs = []
    for r in usable:
        v = r["se"] ** 2
        z2 = (r["beta"] / r["se"]) ** 2
        logs.append(0.5 * (math.log(v / (v + w)) + (w / (v + w)) * z2))
    top = max(logs)
    weights = [math.exp(x - top) for x in logs]
    total = sum(weights)
    pips = [x / total for x in weights]
    order = sorted(range(len(usable)), key=lambda i: -pips[i])
    cumulative = 0.0
    chosen = []
    for i in order:
        chosen.append(i)
        cumulative += pips[i]
        if cumulative >= coverage:
            break
    variants = [
        {
            "variant_id": usable[i]["variant_id"],
            "chrom": usable[i]["chrom"],
            "pos": usable[i]["pos"],
            "log_abf": logs[i],
            "pip": pips[i],
        }
        for i in order
    ]
    return {
        "credible_set": [variants[k]["variant_id"] for k in range(len(chosen))],
        "set_size": len(chosen),
        "cumulative_pip": cumulative,
        "coverage_target": coverage,
        "prior_sd": prior_sd,
        "variants": variants[:1000],
        "method": (
            "Wakefield approximate Bayes factors, single causal variant per region, uniform "
            "prior; "
            "log ABF = 0.5*(ln(V/(V+W)) + W/(V+W)*z^2) with V = se^2, W = prior_sd^2"
        ),
    }, notes


# --------------------------------------------------------------------------- meta-analysis
def _group_studies(source: object, opts: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    studies: dict[str, list[dict[str, Any]]] = {}
    if (
        isinstance(source, Mapping)
        and source
        and all(
            isinstance(v, (Sequence, Mapping))
            and not isinstance(v, str)
            and (isinstance(v, Mapping) or (len(v) > 0 and isinstance(v[0], Mapping)))
            for v in source.values()
        )
        and not any(k in source for k in ("records", "rows", "table", "data"))
    ):
        for name, table in source.items():
            studies[str(name)] = standardize_rows(load_table(table, what=f"study {name}"))[0]
    elif (
        isinstance(source, Sequence)
        and not isinstance(source, (str, bytes))
        and source
        and all(
            isinstance(t, (Sequence, Mapping))
            and not isinstance(t, str)
            and not (isinstance(t, Mapping) and not any(k in t for k in ("records", "rows")))
            for t in source
        )
    ):
        for i, table in enumerate(source, start=1):
            studies[f"study{i}"] = standardize_rows(load_table(table, what=f"study {i}"))[0]
    else:
        records = load_table(source, what="meta-analysis table")
        column = opts.get("study_column", "study")
        if not records or column not in records[0]:
            raise InputFormatError(
                "provide several studies (list/dict of tables) or one table with a 'study' column.",
                code="NO_STUDIES",
                hint="Use {'study1': table1, 'study2': table2} or add a study column.",
            )
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for rec in records:
            grouped[str(rec[column])].append({k: v for k, v in rec.items() if k != column})
        for name, recs in grouped.items():
            studies[name] = standardize_rows(recs)[0]
    if len(studies) < 2:
        raise InputFormatError(
            "meta-analysis needs at least two studies.",
            code="TOO_FEW_STUDIES",
            hint="Provide two or more study tables.",
        )
    return studies


def _meta_inputs(
    source: object, opts: Mapping[str, Any]
) -> list[tuple[str, list[tuple[float, float, str]], dict[str, Any]]]:
    studies = _group_studies(source, opts)
    min_studies = opt_number(opts, "min_studies", 2, minimum=2, integer=True)
    keyed: dict[str, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
    for name, rows in studies.items():
        for r in rows:
            key = (
                r["variant_id"]
                if r["variant_id"]
                else (f"{r['chrom']}:{r['pos']}" if r["chrom"] and r["pos"] is not None else None)
            )
            if key is not None and r["beta"] is not None and r["se"] and r["se"] > 0:
                keyed[key].append((name, r))
    result = []
    skipped = 0
    for key, entries in keyed.items():
        ref = entries[0][1]
        items: list[tuple[float, float, str]] = []
        for name, r in entries:
            beta = r["beta"]
            if ref["effect_allele"] and r["effect_allele"]:
                if (r["effect_allele"], r["other_allele"]) == (
                    ref["effect_allele"],
                    ref["other_allele"],
                ):
                    pass
                elif (r["effect_allele"], r["other_allele"]) == (
                    ref["other_allele"],
                    ref["effect_allele"],
                ):
                    beta = -beta
                else:
                    skipped += 1
                    continue
            items.append((beta, r["se"], name))
        if len(items) >= min_studies:
            result.append((key, items, ref))
    return result


def _fixed(items: Sequence[tuple[float, float, str]]) -> dict[str, float]:
    w = [1.0 / se**2 for _, se, _ in items]
    sw = math.fsum(w)
    beta = math.fsum(wi * b for wi, (b, _, _) in zip(w, items, strict=False)) / sw
    q = math.fsum(wi * (b - beta) ** 2 for wi, (b, _, _) in zip(w, items, strict=False))
    df = len(items) - 1
    se = 1.0 / math.sqrt(sw)
    z = beta / se
    return {
        "beta": beta,
        "se": se,
        "z": z,
        "p": st.two_sided_p_from_z(z),
        "q": q,
        "q_df": df,
        "q_p": st.chi2_sf(q, df),
        "i2": max(0.0, (q - df) / q) if q > 0 else 0.0,
        "sum_w": sw,
        "sum_w2": math.fsum(x * x for x in w),
    }


def fixed_effect_meta_analysis(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    out = []
    for key, items, ref in _meta_inputs(source, opts):
        f = _fixed(items)
        out.append(
            {
                "variant_id": key,
                "chrom": ref["chrom"],
                "pos": ref["pos"],
                "effect_allele": ref["effect_allele"],
                "other_allele": ref["other_allele"],
                "n_studies": len(items),
                "beta": f["beta"],
                "se": f["se"],
                "z": f["z"],
                "p": f["p"],
                "q": f["q"],
                "q_p": f["q_p"],
                "i2": f["i2"],
            }
        )
    return {
        "records": out,
        "n_variants": len(out),
        "model": "fixed effect, inverse-variance weights w=1/se^2; Cochran's Q and I^2 reported",
    }, []


def random_effect_meta_analysis(
    source: object, opts: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    out = []
    for key, items, ref in _meta_inputs(source, opts):
        f = _fixed(items)
        denom = f["sum_w"] - f["sum_w2"] / f["sum_w"]
        tau2 = max(0.0, (f["q"] - f["q_df"]) / denom) if denom > 0 else 0.0
        w = [1.0 / (se**2 + tau2) for _, se, _ in items]
        sw = math.fsum(w)
        beta = math.fsum(wi * b for wi, (b, _, _) in zip(w, items, strict=False)) / sw
        se = 1.0 / math.sqrt(sw)
        z = beta / se
        out.append(
            {
                "variant_id": key,
                "chrom": ref["chrom"],
                "pos": ref["pos"],
                "effect_allele": ref["effect_allele"],
                "other_allele": ref["other_allele"],
                "n_studies": len(items),
                "beta": beta,
                "se": se,
                "z": z,
                "p": st.two_sided_p_from_z(z),
                "tau2": tau2,
                "q": f["q"],
                "q_p": f["q_p"],
                "i2": f["i2"],
            }
        )
    return {
        "records": out,
        "n_variants": len(out),
        "model": (
            "DerSimonian-Laird random effects: tau^2 = max(0, (Q-df)/(S1 - S2/S1)), weights "
            "1/(se^2+tau^2), normal-theory p"
        ),
    }, []


# --------------------------------------------------------------------------- SVG plot
def plot(source: object, opts: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    from ._impl_common import nice_ticks, svg_escape, svg_wrap

    kind = opt_choice(opts, "kind", "manhattan", ("manhattan", "qq"))
    width = int(opt_number(opts, "width", 900, minimum=300, maximum=4000, integer=True))
    height = int(opt_number(opts, "height", 450, minimum=200, maximum=3000, integer=True))
    max_points = opt_number(opts, "max_points", 20_000, minimum=10, integer=True)
    left, right, top, bottom = 60, 20, 30, 50
    body: list[str] = []
    px = left
    pw, ph = width - left - right, height - top - bottom
    if kind == "manhattan":
        data, notes = compute_manhattan_data(source, {**opts, "max_points": max_points})
        pts, total = data["points"], max(1, data["total_length"])
        ymax = (
            max(
                [pt["neglog10p"] for pt in pts if math.isfinite(pt["neglog10p"])]
                + [-math.log10(GENOME_WIDE)]
            )
            * 1.05
        )
        xs = lambda x: px + pw * x / total  # noqa: E731
        ys = lambda y: top + ph * (1 - min(y, ymax) / ymax)  # noqa: E731
        colors = ("#1f4e79", "#7f7f7f")
        order = {c["chrom"]: i for i, c in enumerate(data["chromosomes"])}
        for pt in pts:
            body.append(
                f'<circle cx="{xs(pt["x"]):.1f}" cy="{ys(pt["neglog10p"]):.1f}" r="2" '
                f'fill="{colors[order[pt["chrom"]] % 2]}"/>'
            )
        for name, thr, color in (
            ("genome-wide", -math.log10(GENOME_WIDE), "#c00000"),
            ("suggestive", -math.log10(SUGGESTIVE), "#e69f00"),
        ):
            body.append(
                f'<line x1="{left}" x2="{width - right}" y1="{ys(thr):.1f}" y2="{ys(thr):.1f}" '
                f'stroke="{color}" stroke-dasharray="4 3"><title>{name}</title></line>'
            )
        for c in data["chromosomes"]:
            body.append(
                f'<text x="{xs(c["center"]):.1f}" y="{height - bottom + 16}" font-size="10" '
                f'text-anchor="middle">{svg_escape(c["chrom"])}</text>'
            )
        xlabel, ylabel, n_points = "Chromosome", "-log10(p)", len(pts)
    else:
        data, notes = compute_qq_data(source, {**opts, "max_points": max_points})
        pts = data["points"]
        finite = [pt["observed"] for pt in pts if math.isfinite(pt["observed"])]
        lim = max(finite + [pt["expected"] for pt in pts]) * 1.05
        xs = lambda x: px + pw * x / lim  # noqa: E731
        ys = lambda y: top + ph * (1 - min(y, lim) / lim)  # noqa: E731
        body.append(
            f'<line x1="{xs(0):.1f}" y1="{ys(0):.1f}" x2="{xs(lim):.1f}" y2="{ys(lim):.1f}" '
            'stroke="#c00000"/>'
        )
        for pt in pts:
            body.append(
                f'<circle cx="{xs(pt["expected"]):.1f}" cy="{ys(pt["observed"]):.1f}" r="2" '
                'fill="#1f4e79"/>'
            )
        xlabel, ylabel, n_points = "Expected -log10(p)", "Observed -log10(p)", len(pts)
        ymax = lim
    body.append(
        f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="black"/>'
    )
    for t in nice_ticks(0, ymax, 5):
        y = top + ph * (1 - t / ymax)
        body.append(
            f'<text x="{left - 6}" y="{y + 3:.1f}" font-size="10" text-anchor="end">{t:g}</text>'
        )
    body.append(
        f'<text x="{left + pw / 2}" y="{height - 8}" font-size="12" '
        f'text-anchor="middle">{svg_escape(xlabel)}</text>'
    )
    body.append(
        f'<text x="14" y="{top + ph / 2}" font-size="12" text-anchor="middle" '
        f'transform="rotate(-90 14 {top + ph / 2})">{svg_escape(ylabel)}</text>'
    )
    svg = svg_wrap(width, height, body, f"GWAS {kind} plot")
    path = opts.get("path", opts.get("output"))
    written = (
        write_text_file(path, svg, overwrite=opt_bool(opts, "overwrite", False))
        if path is not None
        else None
    )
    return {"kind": kind, "svg": svg, "n_points_drawn": n_points, "path": written}, notes
