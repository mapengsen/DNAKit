"""Line-oriented format engines: PLINK ``.map``, PLINK ``.ped`` and PAF alignments."""

# ruff: noqa: E501, RUF012

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.io_formats._impl_common import (
    DEFAULT_MAX_RECORDS,
    Format,
    Opts,
    Payload,
    bad_input,
    bad_option,
    chrom_key,
    fmt_number,
    json_payload,
    line_table,
    parse_region,
    to_float,
    to_int,
)


def _public(row: Mapping[str, Any]) -> dict[str, object]:
    return {key: value for key, value in row.items() if not key.startswith("_")}


def _rows_of(records: object, where: str) -> list[Any]:
    if isinstance(records, Mapping):
        for key in ("records", "items", "rows"):
            if key in records:
                records = records[key]
                break
        else:
            records = [records]
    if isinstance(records, (str, bytes)) or not isinstance(records, Sequence):
        raise bad_input(f"{where}: expected a list of records.")
    return list(records)


def _text_payload(fmt: str, lines: list[str], count: int, **info: object) -> Payload:
    return Payload(
        fmt=fmt, text="".join(line + "\n" for line in lines), count=count, info=dict(info)
    )


def _field_problem(name: str, column: int | None, token: str, why: str) -> str:
    """One validation problem about a single field: its name, column and offending value."""

    where = f"field '{name}'" + (f" (column {column})" if column else "")
    return f"{where} has invalid value {token!r}: {why}"


def _int_problem(name: str, column: int | None, token: str) -> tuple[int | None, str | None]:
    try:
        return to_int(token), None
    except ValueError:
        return None, _field_problem(name, column, token, "not an integer")


def _float_problem(name: str, column: int | None, token: str) -> tuple[float | None, str | None]:
    try:
        return to_float(token), None
    except ValueError:
        return None, _field_problem(name, column, token, "not a finite number")


def _as_token(value: object) -> str:
    """Text form of a record value so that it goes through the same checks as a file column."""

    if isinstance(value, float) and value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return str(value)


def _describe(exc: Exception) -> str:
    """A readable message of an exception caught while reading record inputs (no Python reprs)."""

    if isinstance(exc, KeyError):
        return f"missing field {exc.args[0]!r}" if exc.args else "missing field"
    if isinstance(exc, IndexError):
        return "too few values"
    return str(exc)


def _cap(count: int, opts: Opts, noun: str) -> None:
    limit = opts.integer("max_records", DEFAULT_MAX_RECORDS, 1)
    if limit is not None and count > limit:
        raise bad_option(
            f"{count} {noun} exceed max_records={limit}.",
            hint="Raise max_records or page with stream_*.",
        )


# =============================================================================================
# PLINK .map
# =============================================================================================


class MapFormat(Format):
    name = "map"
    title = "PLINK .map variant file"
    noun = "variants"
    suffix = ".map"
    targets = ("map", "bed", "tsv", "json")
    options = {
        "filter": {
            "chromosomes": "list of chromosome names to keep",
            "exclude_chromosomes": "list of chromosome names to drop",
            "ids": "list of variant IDs to keep",
            "exclude_ids": "list of variant IDs to drop",
            "region": "'chrom:start-end' keep variants with start <= bp <= end (1-based, inclusive)",
            "min_position": "int: keep bp >= value",
            "max_position": "int: keep bp <= value",
            "min_cm": "float: keep genetic distance >= value",
            "max_cm": "float: keep genetic distance <= value",
        },
        "sort": {
            "key": "'position' (default: chromosome in natural order, then bp, then ID), 'id' or 'cm'",
            "reverse": "bool, default False",
        },
        "merge": {
            "duplicates": "'error' (default) or 'first': what to do when a variant ID occurs twice",
            "sort": "bool, default False: sort the merged variants by chromosome and bp",
        },
        "split": {
            "by": "'chromosome' (default; one part per chromosome) or 'size'",
            "size": "int: variants per part when by='size'",
        },
        "index": {},
        "summarize": {},
        "validate": {},
        "convert": {},
        "read": {},
        "write": {},
        "stream": {},
    }
    semantics = {
        "validate": "Each non-blank line needs 4 whitespace-separated columns: chromosome, variant ID, genetic "
        "distance in cM (finite number) and base-pair position (integer >= 0). Duplicate IDs are errors; "
        "duplicate (chromosome, bp) pairs are warnings. Every problem of a line is reported (line number, field name and column, "
        "offending value).",
        "index": "Maps each variant ID to its line number, byte offset and length in the text, and lists per-chromosome "
        "counts and bp ranges, so a variant can be located without re-parsing.",
        "convert": "Targets: 'map' (normalised, tab-separated), 'bed' (chrom, bp-1, bp, ID; BED is 0-based half-open), "
        "'tsv' (with header chrom, id, cm, pos) and 'json'.",
    }

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for number, (offset, length, line) in enumerate(line_table(data), 1):
            if not line.strip():
                continue
            fields = line.split()
            where = f"line {number}"
            if len(fields) != 4:
                self.bad(
                    errors, where, f"expected 4 columns (chrom, id, cM, bp), found {len(fields)}"
                )
                continue
            cm, cm_problem = _float_problem("cM", 3, fields[2])
            position, pos_problem = _int_problem("bp", 4, fields[3])
            if position is not None and position < 0:
                pos_problem = _field_problem("bp", 4, fields[3], "negative bp position")
            problems = [p for p in (cm_problem, pos_problem) if p]
            for problem in problems:
                self.bad(errors, where, problem)
            if problems:
                continue
            rows.append(
                {
                    "chrom": fields[0],
                    "id": fields[1],
                    "cm": cm,
                    "pos": position,
                    "_line": number,
                    "_offset": offset,
                    "_length": length,
                }
            )
        return rows

    def from_records(self, records: object, opts: Opts) -> list[dict[str, Any]]:
        rows = []
        for number, item in enumerate(_rows_of(records, "map"), 1):
            try:
                if isinstance(item, Mapping):
                    chrom, ident = str(item["chrom"]), str(item["id"])
                    cm_value, pos_value = item.get("cm", 0.0), item["pos"]
                else:
                    values = list(item)
                    if len(values) == 3:
                        values.insert(2, 0.0)
                    chrom, ident, cm_value, pos_value = (
                        str(values[0]),
                        str(values[1]),
                        values[2],
                        values[3],
                    )
                cm, cm_problem = _float_problem("cm", None, _as_token(cm_value))
                pos, pos_problem = _int_problem("pos", None, _as_token(pos_value))
                if cm_problem or pos_problem:
                    raise ValueError("; ".join(p for p in (cm_problem, pos_problem) if p))
                assert cm is not None and pos is not None
            except (KeyError, TypeError, ValueError, IndexError) as exc:
                raise bad_input(
                    f"map record {number} is malformed: {_describe(exc)}",
                    hint="Use {'chrom','id','cm','pos'} dicts or [chrom, id, cm, bp] lists.",
                ) from exc
            if pos < 0 or not chrom or not ident or any(c.isspace() for c in chrom + ident):
                raise bad_input(
                    f"map record {number} has an empty/whitespace name or negative position."
                )
            rows.append(
                {"chrom": chrom, "id": ident, "cm": cm, "pos": pos, "_line": number, "_offset": None, "_length": None}
            )
        return rows

    def records(self, doc: list[dict[str, Any]]) -> list[dict[str, object]]:
        return [_public(row) for row in doc]

    def dump(self, doc: list[dict[str, Any]], opts: Opts) -> Payload:
        return _text_payload(
            "map",
            [f"{r['chrom']}\t{r['id']}\t{fmt_number(r['cm'])}\t{r['pos']}" for r in doc],
            len(doc),
        )

    def validate_doc(
        self, doc: list[dict[str, Any]], opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        seen: dict[str, int] = {}
        positions: dict[tuple[str, int], int] = {}
        for row in doc:
            if row["id"] in seen:
                errors.append(
                    f"line {row['_line']}: duplicate variant ID {row['id']!r} (first on line {seen[row['id']]})"
                )
            seen.setdefault(row["id"], row["_line"])
            key = (row["chrom"], row["pos"])
            if key in positions:
                warnings.append(
                    f"line {row['_line']}: position {row['chrom']}:{row['pos']} repeats line {positions[key]}"
                )
            positions.setdefault(key, row["_line"])
        return {"variant_count": len(doc)}

    def filter(self, doc: list[dict[str, Any]], opts: Opts) -> list[dict[str, Any]]:
        keep_chrom, drop_chrom = opts.name_set("chromosomes"), opts.name_set("exclude_chromosomes")
        keep_ids, drop_ids = opts.name_set("ids"), opts.name_set("exclude_ids")
        region = opts.text("region")
        min_pos, max_pos = (
            opts.integer("min_position", None, 0),
            opts.integer("max_position", None, 0),
        )
        min_cm, max_cm = opts.number("min_cm"), opts.number("max_cm")
        reg = parse_region(region) if region else None
        out = []
        for row in doc:
            if keep_chrom is not None and row["chrom"] not in keep_chrom:
                continue
            if drop_chrom is not None and row["chrom"] in drop_chrom:
                continue
            if keep_ids is not None and row["id"] not in keep_ids:
                continue
            if drop_ids is not None and row["id"] in drop_ids:
                continue
            if reg and not (row["chrom"] == reg[0] and reg[1] <= row["pos"] <= reg[2]):
                continue
            if min_pos is not None and row["pos"] < min_pos:
                continue
            if max_pos is not None and row["pos"] > max_pos:
                continue
            if min_cm is not None and row["cm"] < min_cm:
                continue
            if max_cm is not None and row["cm"] > max_cm:
                continue
            out.append(row)
        return out

    def sort(self, doc: list[dict[str, Any]], opts: Opts) -> list[dict[str, Any]]:
        key = opts.text("key", "position", ["position", "id", "cm"])
        reverse = opts.boolean("reverse", False)
        keys = {
            "position": lambda r: (chrom_key(r["chrom"]), r["pos"], r["id"]),
            "id": lambda r: (r["id"], chrom_key(r["chrom"]), r["pos"]),
            "cm": lambda r: (r["cm"], chrom_key(r["chrom"]), r["pos"]),
        }
        return sorted(doc, key=keys[str(key)], reverse=reverse)

    def merge(self, docs: list[list[dict[str, Any]]], opts: Opts) -> list[dict[str, Any]]:
        policy = opts.text("duplicates", "error", ["error", "first"])
        order = opts.boolean("sort", False)
        out: list[dict[str, Any]] = []
        seen: set[str] = set()
        for doc in docs:
            for row in doc:
                if row["id"] in seen:
                    if policy == "first":
                        continue
                    raise bad_input(
                        f"Duplicate variant ID {row['id']!r} while merging maps.",
                        code="IO_FORMATS_DUPLICATE_ID",
                        hint="Pass duplicates='first' to keep the first occurrence.",
                    )
                seen.add(row["id"])
                out.append(row)
        return (
            sorted(out, key=lambda r: (chrom_key(r["chrom"]), r["pos"], r["id"])) if order else out
        )

    def split(
        self, doc: list[dict[str, Any]], opts: Opts
    ) -> list[tuple[str, list[dict[str, Any]]]]:
        by = opts.text("by", "chromosome", ["chromosome", "size"])
        size = opts.integer("size", None, 1)
        if by == "size":
            if size is None:
                raise bad_option("split_map with by='size' needs the size option.")
            return [(f"{i // size + 1:04d}", doc[i : i + size]) for i in range(0, len(doc), size)]
        groups: dict[str, list[dict[str, Any]]] = {}
        for row in doc:
            groups.setdefault(row["chrom"], []).append(row)
        return [(chrom, groups[chrom]) for chrom in sorted(groups, key=chrom_key)]

    def summarize(self, doc: list[dict[str, Any]], opts: Opts) -> dict[str, object]:
        chroms: dict[str, dict[str, object]] = {}
        for row in doc:
            entry = chroms.setdefault(
                row["chrom"], {"variant_count": 0, "min_pos": row["pos"], "max_pos": row["pos"]}
            )
            entry["variant_count"] += 1  # type: ignore[operator]
            entry["min_pos"] = min(entry["min_pos"], row["pos"])  # type: ignore[type-var]
            entry["max_pos"] = max(entry["max_pos"], row["pos"])  # type: ignore[type-var]
        ids = Counter(r["id"] for r in doc)
        ordered = sorted(doc, key=lambda r: (chrom_key(r["chrom"]), r["pos"]))
        return {
            "variant_count": len(doc),
            "chromosome_count": len(chroms),
            "chromosomes": {c: chroms[c] for c in sorted(chroms, key=chrom_key)},
            "duplicate_id_count": sum(1 for n in ids.values() if n > 1),
            "cm_min": min((r["cm"] for r in doc), default=None),
            "cm_max": max((r["cm"] for r in doc), default=None),
            "all_cm_zero": all(r["cm"] == 0 for r in doc),
            "sorted": [(r["chrom"], r["pos"]) for r in doc]
            == [(r["chrom"], r["pos"]) for r in ordered],
        }

    def index(self, doc: list[dict[str, Any]], opts: Opts) -> dict[str, object]:
        entries = {
            r["id"]: {
                "line": r["_line"],
                "offset": r["_offset"],
                "length": r["_length"],
                "chrom": r["chrom"],
                "pos": r["pos"],
            }
            for r in doc
        }
        chroms = self.summarize(doc, opts)["chromosomes"]
        first: dict[str, int] = {}
        for r in doc:
            first.setdefault(r["chrom"], r["_line"])
        for chrom, entry in chroms.items():  # type: ignore[union-attr]
            entry["first_line"] = first[chrom]
        counts = Counter(r["id"] for r in doc)
        duplicates = [ident for ident, n in counts.items() if n > 1]
        result: dict[str, object] = {
            "key": "id",
            "record_count": len(doc),
            "entries": entries,
            "chromosomes": chroms,
            "duplicate_ids": duplicates,
        }
        if duplicates:
            result["warnings"] = [
                f"{len(duplicates)} variant ID(s) occur on several lines ({', '.join(duplicates[:5])}"
                f"{', ...' if len(duplicates) > 5 else ''}): each entry keeps the last line of its ID, "
                f"while record_count ({len(doc)}) counts every line."
            ]
        return result

    def convert(self, doc: list[dict[str, Any]], target: str, opts: Opts) -> Payload:
        if target == "map":
            return self.dump(doc, opts)
        if target == "bed":
            return _text_payload(
                "bed",
                [
                    f"{r['chrom']}\t{r['pos'] - 1}\t{r['pos']}\t{r['id']}"
                    for r in doc
                    if r["pos"] >= 1
                ],
                len(doc),
            )
        if target == "tsv":
            lines = ["chrom\tid\tcm\tpos"] + [
                f"{r['chrom']}\t{r['id']}\t{fmt_number(r['cm'])}\t{r['pos']}" for r in doc
            ]
            return _text_payload("tsv", lines, len(doc))
        return json_payload("map", self.records(doc), len(doc))


# =============================================================================================
# PLINK .ped
# =============================================================================================

_SEX_NAMES = {"male": 1, "female": 2, "unknown": 0}
_PHENO_CLASSES = ("case", "control", "missing", "observed")


def _pheno_class(token: str) -> str:
    value = float(token)
    if value == 2:
        return "case"
    if value == 1:
        return "control"
    if value in (0, -9):
        return "missing"
    return "quantitative"


def _geno(item: object) -> tuple[str, str]:
    if isinstance(item, str):
        parts = item.split("/") if "/" in item else list(item)
    else:
        parts = [str(x) for x in item]  # type: ignore[attr-defined]
    if len(parts) != 2 or not all(parts):
        raise ValueError(f"genotype {item!r} must have exactly two alleles")
    return parts[0], parts[1]


def _missing(call: tuple[str, str]) -> bool:
    return call[0] == "0" or call[1] == "0"


class PedFormat(Format):
    name = "ped"
    title = "PLINK .ped pedigree/genotype file"
    noun = "individuals"
    suffix = ".ped"
    targets = ("ped", "tsv", "json")
    options = {
        "read": {
            "max_cells": "int, default 20000000: refuse files with more genotype cells (individuals x variants)"
        },
        "validate": {
            "variant_count": "int: expected number of variants (e.g. the line count of the matching .map)",
            "max_cells": "int, default 20000000",
        },
        "filter": {
            "families": "list of family IDs to keep",
            "exclude_families": "list of family IDs to drop",
            "individuals": "list of individual IDs to keep",
            "exclude_individuals": "list of individual IDs to drop",
            "sex": "keep one sex: 1/'male', 2/'female' or 0/'unknown'",
            "phenotype": "keep one class: 'case' (2), 'control' (1), 'missing' (0 or -9), "
            "'observed' (not missing) or 'quantitative' (any other value)",
            "founders_only": "bool, default False: keep individuals whose father and mother IDs are both 0",
            "max_missing_rate": "float in [0,1]: drop individuals whose fraction of missing genotypes exceeds it",
            "variant_indices": "list of 0-based variant columns to keep (apply the same to the .map file)",
            "max_variant_missing_rate": "float in [0,1]: drop variant columns whose missing-call fraction exceeds it "
            "(computed over the individuals that pass the other filters)",
            "max_cells": "int, default 20000000",
        },
        "sort": {
            "key": "'family' (default: family ID then individual ID), 'individual', 'sex' or 'phenotype'",
            "reverse": "bool, default False",
            "max_cells": "int, default 20000000",
        },
        "merge": {
            "duplicates": "'error' (default) or 'first': what to do when a (family, individual) pair occurs twice",
            "max_cells": "int, default 20000000",
        },
        "split": {
            "by": "'family' (default; one part per family) or 'size'",
            "size": "int: individuals per part when by='size'",
            "max_cells": "int, default 20000000",
        },
        "summarize": {"max_cells": "int, default 20000000"},
        "index": {"max_cells": "int, default 20000000"},
        "convert": {"max_cells": "int, default 20000000"},
        "stream": {"max_cells": "int, default 20000000"},
        "write": {"max_cells": "int, default 20000000"},
    }
    semantics = {
        "validate": "Rows need 6 + 2k whitespace-separated columns (family, individual, father, mother, sex, phenotype, "
        "then two alleles per variant) with the same k on every row, numeric phenotype and integer sex. Duplicate "
        "(family, individual) pairs are errors; sex outside {0,1,2} and parents absent from their family are warnings. "
        "Allele '0' marks a missing call. Every problem of a row is reported (line number, field name and column, offending value).",
        "index": "Maps 'family:individual' to its line number, byte offset and length and lists the individuals of each family.",
        "convert": "Targets: 'ped' (normalised), 'tsv' (header FID IID PAT MAT SEX PHENO V1..Vk, genotypes as 'a1/a2') and 'json'.",
        "filter": "Individual filters are combined with AND; variant filters drop genotype columns from every kept row.",
    }

    def _cells(self, opts: Opts) -> int:
        return opts.integer("max_cells", 20_000_000, 1) or 20_000_000

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> dict[str, Any]:
        rows: list[dict[str, Any]] = []
        k: int | None = None
        cells = 0
        limit = self._cells(opts)
        for number, (offset, length, line) in enumerate(line_table(data), 1):
            if not line.strip():
                continue
            fields = line.split()
            where = f"line {number}"
            problems: list[str] = []
            if len(fields) < 6 or (len(fields) - 6) % 2:
                problems.append(f"expected 6 + 2k columns, found {len(fields)}")
            else:
                width = (len(fields) - 6) // 2
                if k is None:
                    k = width
                elif width != k:
                    problems.append(f"has {width} variants but earlier rows have {k}")
            sex: int | None = None
            if len(fields) >= 6:
                sex, sex_problem = _int_problem("sex", 5, fields[4])
                _, pheno_problem = _float_problem("phenotype", 6, fields[5])
                problems.extend(p for p in (sex_problem, pheno_problem) if p)
            for problem in problems:
                self.bad(errors, where, problem)
            if problems:
                continue
            assert sex is not None
            cells += width
            if cells > limit:
                raise bad_option(
                    f"The .ped has more than max_cells={limit} genotype cells.",
                    hint="Raise max_cells.",
                )
            calls = [(fields[6 + 2 * i], fields[7 + 2 * i]) for i in range(width)]
            rows.append(
                {
                    "fid": fields[0],
                    "iid": fields[1],
                    "pat": fields[2],
                    "mat": fields[3],
                    "sex": sex,
                    "pheno": fields[5],
                    "genotypes": calls,
                    "_line": number,
                    "_offset": offset,
                    "_length": length,
                }
            )
        return {"rows": rows, "variants": k or 0}

    def from_records(self, records: object, opts: Opts) -> dict[str, Any]:
        rows, k = [], None
        for number, item in enumerate(_rows_of(records, "ped"), 1):
            try:
                if not isinstance(item, Mapping):
                    raise TypeError("record must be a mapping")
                calls = [_geno(g) for g in item.get("genotypes", [])]
                row = {
                    "fid": str(item["fid"]),
                    "iid": str(item["iid"]),
                    "pat": str(item.get("pat", "0")),
                    "mat": str(item.get("mat", "0")),
                    "sex": 0,
                    "pheno": str(item.get("pheno", "-9")),
                    "genotypes": calls,
                }
                sex, sex_problem = _int_problem("sex", None, _as_token(item.get("sex", 0)))
                _, pheno_problem = _float_problem("pheno", None, row["pheno"])
                if sex_problem or pheno_problem:
                    raise ValueError("; ".join(p for p in (sex_problem, pheno_problem) if p))
                assert sex is not None
                row["sex"] = sex
                for key in ("fid", "iid", "pat", "mat"):
                    if not row[key] or any(c.isspace() for c in row[key]):
                        raise ValueError(f"{key} must be a non-empty string without whitespace")
                if any(any(c.isspace() for c in a) for call in calls for a in call):
                    raise ValueError("alleles must not contain whitespace")
            except (KeyError, TypeError, ValueError) as exc:
                raise bad_input(
                    f"ped record {number} is malformed: {_describe(exc)}",
                    hint="Use {'fid','iid','pat','mat','sex','pheno','genotypes':[[a1,a2],...]} dicts.",
                ) from exc
            if k is None:
                k = len(calls)
            elif len(calls) != k:
                raise bad_input(f"ped record {number} has {len(calls)} genotypes, expected {k}.")
            row.update({"_line": number, "_offset": None, "_length": None})
            rows.append(row)
        return {"rows": rows, "variants": k or 0}

    def records(self, doc: dict[str, Any]) -> list[dict[str, object]]:
        return [{**_public(r), "genotypes": [list(c) for c in r["genotypes"]]} for r in doc["rows"]]

    def describe(self, doc: dict[str, Any]) -> dict[str, object]:
        return {"variant_count": doc["variants"]}

    def count(self, doc: dict[str, Any]) -> int:
        return len(doc["rows"])

    def _line(self, r: dict[str, Any]) -> str:
        alleles = "\t".join(f"{a}\t{b}" for a, b in r["genotypes"])
        head = f"{r['fid']}\t{r['iid']}\t{r['pat']}\t{r['mat']}\t{r['sex']}\t{r['pheno']}"
        return head + ("\t" + alleles if alleles else "")

    def dump(self, doc: dict[str, Any], opts: Opts) -> Payload:
        self._cells(opts)
        return _text_payload(
            "ped",
            [self._line(r) for r in doc["rows"]],
            len(doc["rows"]),
            variant_count=doc["variants"],
        )

    def validate_doc(
        self, doc: dict[str, Any], opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        expected = opts.integer("variant_count", None, 0)
        if expected is not None and doc["rows"] and doc["variants"] != expected:
            errors.append(f"file has {doc['variants']} variants but variant_count={expected}")
        seen: dict[tuple[str, str], int] = {}
        members = {(r["fid"], r["iid"]) for r in doc["rows"]}
        for r in doc["rows"]:
            key = (r["fid"], r["iid"])
            if key in seen:
                errors.append(
                    f"line {r['_line']}: duplicate individual {r['fid']}:{r['iid']} (first on line {seen[key]})"
                )
            seen.setdefault(key, r["_line"])
            if r["sex"] not in (0, 1, 2):
                warnings.append(
                    f"line {r['_line']}: sex code {r['sex']} is not 0, 1 or 2 (treated as unknown)"
                )
            for role in ("pat", "mat"):
                if r[role] != "0" and (r["fid"], r[role]) not in members:
                    warnings.append(
                        f"line {r['_line']}: {role} {r[role]!r} of {r['fid']}:{r['iid']} is not in family {r['fid']}"
                    )
        return {"individual_count": len(doc["rows"]), "variant_count": doc["variants"]}

    # -- helpers -------------------------------------------------------------------------
    @staticmethod
    def _miss_rate(row: dict[str, Any]) -> float:
        calls = row["genotypes"]
        return sum(1 for c in calls if _missing(c)) / len(calls) if calls else 0.0

    def filter(self, doc: dict[str, Any], opts: Opts) -> dict[str, Any]:
        families, drop_f = opts.name_set("families"), opts.name_set("exclude_families")
        people, drop_p = opts.name_set("individuals"), opts.name_set("exclude_individuals")
        sex = opts.get("sex")
        pheno = opts.text("phenotype", None, (*_PHENO_CLASSES, "quantitative"))
        founders = opts.boolean("founders_only", False)
        max_missing = opts.number("max_missing_rate", None, 0.0, 1.0)
        columns = opts.get("variant_indices")
        max_var_missing = opts.number("max_variant_missing_rate", None, 0.0, 1.0)
        self._cells(opts)
        wanted_sex: int | None = None
        if sex is not None:
            if isinstance(sex, str) and sex.lower() in _SEX_NAMES:
                wanted_sex = _SEX_NAMES[sex.lower()]
            elif isinstance(sex, int) and not isinstance(sex, bool) and sex in (0, 1, 2):
                wanted_sex = sex
            else:
                raise bad_option("option 'sex' must be 0/1/2 or 'unknown'/'male'/'female'.")
        rows = []
        for r in doc["rows"]:
            if (families is not None and r["fid"] not in families) or (
                drop_f is not None and r["fid"] in drop_f
            ):
                continue
            if (people is not None and r["iid"] not in people) or (
                drop_p is not None and r["iid"] in drop_p
            ):
                continue
            if wanted_sex is not None and (r["sex"] if r["sex"] in (1, 2) else 0) != wanted_sex:
                continue
            if pheno is not None:
                cls = _pheno_class(r["pheno"])
                if not (cls == pheno or (pheno == "observed" and cls != "missing")):
                    continue
            if founders and not (r["pat"] == "0" and r["mat"] == "0"):
                continue
            if max_missing is not None and self._miss_rate(r) > max_missing:
                continue
            rows.append(r)
        keep = list(range(doc["variants"]))
        if columns is not None:
            if (
                not isinstance(columns, Sequence)
                or isinstance(columns, str)
                or not all(
                    isinstance(c, int) and not isinstance(c, bool) and 0 <= c < doc["variants"]
                    for c in columns
                )
            ):
                raise bad_option(f"variant_indices must be integers in [0, {doc['variants']}).")
            keep = sorted(set(columns))
        if max_var_missing is not None and rows:
            keep = [
                j
                for j in keep
                if sum(1 for r in rows if _missing(r["genotypes"][j])) / len(rows)
                <= max_var_missing
            ]
        if len(keep) != doc["variants"]:
            rows = [{**r, "genotypes": [r["genotypes"][j] for j in keep]} for r in rows]
        return {"rows": rows, "variants": len(keep), "kept_variant_indices": keep}

    def sort(self, doc: dict[str, Any], opts: Opts) -> dict[str, Any]:
        key = opts.text("key", "family", ["family", "individual", "sex", "phenotype"])
        reverse = opts.boolean("reverse", False)
        keys = {
            "family": lambda r: (r["fid"], r["iid"]),
            "individual": lambda r: (r["iid"], r["fid"]),
            "sex": lambda r: (r["sex"], r["fid"], r["iid"]),
            "phenotype": lambda r: (float(r["pheno"]), r["fid"], r["iid"]),
        }
        return {
            "rows": sorted(doc["rows"], key=keys[str(key)], reverse=reverse),
            "variants": doc["variants"],
        }

    def merge(self, docs: list[dict[str, Any]], opts: Opts) -> dict[str, Any]:
        policy = opts.text("duplicates", "error", ["error", "first"])
        widths = {d["variants"] for d in docs if d["rows"]}
        if len(widths) > 1:
            raise bad_input(
                f"Cannot merge .ped files with different variant counts {sorted(widths)}.",
                hint="Subset every file to the same variants first (filter_ped variant_indices).",
            )
        seen, rows = set(), []
        for d in docs:
            for r in d["rows"]:
                key = (r["fid"], r["iid"])
                if key in seen:
                    if policy == "first":
                        continue
                    raise bad_input(
                        f"Duplicate individual {r['fid']}:{r['iid']} while merging.",
                        code="IO_FORMATS_DUPLICATE_ID",
                        hint="Pass duplicates='first' to keep the first row.",
                    )
                seen.add(key)
                rows.append(r)
        return {"rows": rows, "variants": widths.pop() if widths else 0}

    def split(self, doc: dict[str, Any], opts: Opts) -> list[tuple[str, dict[str, Any]]]:
        by = opts.text("by", "family", ["family", "size"])
        size = opts.integer("size", None, 1)
        rows = doc["rows"]
        if by == "size":
            if size is None:
                raise bad_option("split_ped with by='size' needs the size option.")
            return [
                (f"{i // size + 1:04d}", {"rows": rows[i : i + size], "variants": doc["variants"]})
                for i in range(0, len(rows), size)
            ]
        groups: dict[str, list[dict[str, Any]]] = {}
        for r in rows:
            groups.setdefault(r["fid"], []).append(r)
        return [(fid, {"rows": g, "variants": doc["variants"]}) for fid, g in groups.items()]

    def summarize(self, doc: dict[str, Any], opts: Opts) -> dict[str, object]:
        rows = doc["rows"]
        pheno = Counter(_pheno_class(r["pheno"]) for r in rows)
        sex = Counter({1: "male", 2: "female"}.get(r["sex"], "unknown") for r in rows)
        calls = sum(len(r["genotypes"]) for r in rows)
        missing = sum(1 for r in rows for c in r["genotypes"] if _missing(c))
        rates = [self._miss_rate(r) for r in rows]
        return {
            "individual_count": len(rows),
            "family_count": len({r["fid"] for r in rows}),
            "variant_count": doc["variants"],
            "founder_count": sum(1 for r in rows if r["pat"] == "0" and r["mat"] == "0"),
            "sex": {k: sex.get(k, 0) for k in ("male", "female", "unknown")},
            "phenotype": {
                k: pheno.get(k, 0) for k in ("case", "control", "missing", "quantitative")
            },
            "genotype_missing_rate": missing / calls if calls else 0.0,
            "individual_missing_rate": {
                "min": min(rates, default=None),
                "mean": sum(rates) / len(rates) if rates else None,
                "max": max(rates, default=None),
            },
        }

    def index(self, doc: dict[str, Any], opts: Opts) -> dict[str, object]:
        entries = {
            f"{r['fid']}:{r['iid']}": {
                "line": r["_line"],
                "offset": r["_offset"],
                "length": r["_length"],
            }
            for r in doc["rows"]
        }
        families: dict[str, list[str]] = {}
        for r in doc["rows"]:
            families.setdefault(r["fid"], []).append(r["iid"])
        return {
            "key": "family:individual",
            "record_count": len(entries),
            "entries": entries,
            "families": families,
        }

    def convert(self, doc: dict[str, Any], target: str, opts: Opts) -> Payload:
        if target == "ped":
            return self.dump(doc, opts)
        if target == "tsv":
            head = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO"] + [
                f"V{i + 1}" for i in range(doc["variants"])
            ]
            lines = ["\t".join(head)]
            for r in doc["rows"]:
                lines.append(
                    "\t".join(
                        [r["fid"], r["iid"], r["pat"], r["mat"], str(r["sex"]), r["pheno"]]
                        + [f"{a}/{b}" for a, b in r["genotypes"]]
                    )
                )
            return _text_payload("tsv", lines, len(doc["rows"]))
        return json_payload(
            "ped",
            {"variant_count": doc["variants"], "records": self.records(doc)},
            len(doc["rows"]),
        )


# =============================================================================================
# PAF
# =============================================================================================

_PAF_INT = (
    "query_length",
    "query_start",
    "query_end",
    "target_length",
    "target_start",
    "target_end",
    "matches",
    "block_length",
    "mapq",
)
_PAF_COLUMNS = (
    "query_name",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target_name",
    "target_length",
    "target_start",
    "target_end",
    "matches",
    "block_length",
    "mapq",
)
_TAG = re.compile(r"[A-Za-z][A-Za-z0-9]:[AifZHB]:")


def _paf_tag_value(raw: str) -> tuple[str, str, object]:
    tag, kind, value = raw.split(":", 2)
    if kind == "i":
        return tag, kind, to_int(value)
    if kind == "f":
        return tag, kind, to_float(value)
    if kind == "A" and len(value) != 1:
        raise ValueError(f"tag {tag} of type A needs exactly one character")
    return tag, kind, value


def _paf_row_problems(row: dict[str, Any]) -> list[str]:
    """Every consistency problem of a PAF row, checking whatever its (possibly partial) columns allow."""

    def has(*keys: str) -> bool:
        return all(key in row for key in keys)

    problems: list[str] = []
    for key in ("query_length", "target_length"):
        if key in row and row[key] <= 0:
            problems.append(_field_problem(key, None, str(row[key]), "must be positive"))
    if has("query_length", "query_start", "query_end") and not (
        0 <= row["query_start"] <= row["query_end"] <= row["query_length"]
    ):
        problems.append("query coordinates must satisfy 0 <= start <= end <= length")
    if has("target_length", "target_start", "target_end") and not (
        0 <= row["target_start"] <= row["target_end"] <= row["target_length"]
    ):
        problems.append("target coordinates must satisfy 0 <= start <= end <= length")
    if row["strand"] not in ("+", "-"):
        problems.append(_field_problem("strand", 5, row["strand"], "must be '+' or '-'"))
    if has("matches", "block_length") and not 0 <= row["matches"] <= row["block_length"]:
        problems.append("matches must satisfy 0 <= matches <= block_length")
    if "mapq" in row and not 0 <= row["mapq"] <= 255:
        problems.append(_field_problem("mapq", 12, str(row["mapq"]), "must be in [0, 255]"))
    return problems


def _paf_row_check(row: dict[str, Any]) -> str | None:
    problems = _paf_row_problems(row)
    return problems[0] if problems else None


def _identity(row: dict[str, Any]) -> float:
    return row["matches"] / row["block_length"] if row["block_length"] else 0.0


def _coverage(row: dict[str, Any]) -> float:
    return (row["query_end"] - row["query_start"]) / row["query_length"]


class PafFormat(Format):
    name = "paf"
    title = "PAF (Pairwise mApping Format) alignment file"
    noun = "alignments"
    suffix = ".paf"
    targets = ("paf", "bed", "tsv", "json")
    options = {
        "filter": {
            "min_mapq": "int: keep mapping quality >= value",
            "min_identity": "float: keep matches/block_length >= value",
            "min_block_length": "int: keep alignment block length >= value",
            "min_matches": "int: keep matching bases >= value",
            "min_query_coverage": "float: keep (query_end-query_start)/query_length >= value",
            "strand": "'+' or '-'",
            "query_names": "list of query names to keep",
            "exclude_query_names": "list to drop",
            "target_names": "list of target names to keep",
            "exclude_target_names": "list to drop",
            "target_region": "'name:start-end' keep alignments overlapping the 0-based half-open target interval",
            "query_region": "'name:start-end' keep alignments overlapping the 0-based half-open query interval",
            "tags": "mapping {tag: value}: keep alignments carrying every tag with that value (e.g. {'tp': 'P'})",
        },
        "sort": {
            "key": "'target' (default: target name, start, end), 'query' (query name, start, end), 'mapq', "
            "'identity', 'block_length' or 'matches'",
            "reverse": "bool, default False",
        },
        "merge": {
            "deduplicate": "bool, default False: drop alignment lines identical to an earlier one",
            "sort_by": "optional sort key as for sort_paf, applied to the merged records",
        },
        "split": {
            "by": "'query' (default), 'target', 'strand' or 'size'",
            "size": "int: alignments per part when by='size'",
        },
        "index": {"by": "'query' (default) or 'target': the name to index alignments by"},
        "summarize": {},
        "validate": {},
        "convert": {},
        "read": {},
        "write": {},
        "stream": {},
    }
    semantics = {
        "validate": "Each non-blank line needs the 12 mandatory tab-separated PAF columns (query name/length/start/end, strand, "
        "target name/length/start/end, matching bases, block length, mapping quality) with start <= end <= length, "
        "matches <= block length, strand +/- and mapq in [0,255]; optional columns must look like TAG:TYPE:VALUE. Every problem of a line is reported (line number, field name and column, "
        "offending value).",
        "summarize": "Counts alignments, distinct queries/targets, strands, aligned bases (sum of block lengths), identity = "
        "matches / block_length (mean and base-weighted overall), and a mapq histogram.",
        "index": "Maps each query (or target) name to its alignments' line number, byte offset and length, plus the interval "
        "on the indexed sequence, sorted by start.",
        "convert": "Targets: 'paf' (normalised), 'bed' (target interval: target, tstart, tend, query, mapq, strand), "
        "'tsv' (header plus the 12 columns and one column per tag) and 'json'.",
    }

    def parse(self, data: bytes, errors: list[str] | None, opts: Opts) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for number, (offset, length, line) in enumerate(line_table(data), 1):
            if not line.strip():
                continue
            fields = line.split("\t")
            where = f"line {number}"
            if len(fields) < 12:
                self.bad(
                    errors,
                    where,
                    f"expected at least 12 tab-separated columns, found {len(fields)}",
                )
                continue
            row: dict[str, Any] = {
                "query_name": fields[0],
                "strand": fields[4],
                "target_name": fields[5],
            }
            problems = []
            for column, key in enumerate(_PAF_COLUMNS, 1):
                if key in _PAF_INT:
                    value, problem = _int_problem(key, column, fields[column - 1])
                    if problem:
                        problems.append(problem)
                    else:
                        row[key] = value
            for position, raw in enumerate(fields[12:], 13):
                if not _TAG.match(raw):
                    problems.append(
                        _field_problem("optional", position, raw, "not TAG:TYPE:VALUE")
                    )
                    continue
                try:
                    _paf_tag_value(raw)
                except ValueError:
                    problems.append(
                        _field_problem(
                            f"tag {raw[:2]}", position, raw, "the value does not match its type"
                        )
                    )
            problems.extend(_paf_row_problems(row))  # checks whatever the numeric columns allow
            for problem in problems:
                self.bad(errors, where, problem)
            if problems:
                continue
            row.update(
                {"_tags": fields[12:], "_line": number, "_offset": offset, "_length": length}
            )
            rows.append(row)
        return rows

    def _tags(self, row: dict[str, Any]) -> dict[str, dict[str, object]]:
        out = {}
        for raw in row["_tags"]:
            tag, kind, value = _paf_tag_value(raw)
            out[tag] = {"type": kind, "value": value}
        return out

    def from_records(self, records: object, opts: Opts) -> list[dict[str, Any]]:
        rows = []
        for number, item in enumerate(_rows_of(records, "paf"), 1):
            try:
                if not isinstance(item, Mapping):
                    raise TypeError("record must be a mapping")
                row: dict[str, Any] = {
                    "query_name": str(item["query_name"]),
                    "strand": str(item["strand"]),
                    "target_name": str(item["target_name"]),
                }
                for key in _PAF_INT:
                    value = item[key]
                    if isinstance(value, bool) or not isinstance(value, int):
                        raise ValueError(f"{key} must be an integer")
                    row[key] = value
                if (
                    not row["query_name"]
                    or not row["target_name"]
                    or any("\t" in row[k] or "\n" in row[k] for k in ("query_name", "target_name"))
                ):
                    raise ValueError("names must be non-empty and free of tabs/newlines")
                problem = _paf_row_check(row)
                if problem:
                    raise ValueError(problem)
                raw_tags: list[str] = []
                tags = item.get("tags", {})
                if isinstance(tags, Mapping):
                    for tag, value in tags.items():
                        if isinstance(value, Mapping):
                            kind, val = value["type"], value["value"]
                        elif isinstance(value, bool):
                            raise ValueError("boolean tag values are not supported")
                        elif isinstance(value, int):
                            kind, val = "i", value
                        elif isinstance(value, float):
                            kind, val = "f", fmt_number(value)
                        else:
                            kind, val = "Z", value
                        raw_tags.append(f"{tag}:{kind}:{val}")
                else:
                    raw_tags = [str(t) for t in tags]
                for raw in raw_tags:
                    if not _TAG.match(raw) or "\t" in raw or "\n" in raw:
                        raise ValueError(f"tag {raw!r} is not TAG:TYPE:VALUE")
                    _paf_tag_value(raw)
                row["_tags"] = raw_tags
            except (KeyError, TypeError, ValueError) as exc:
                raise bad_input(
                    f"paf record {number} is malformed: {exc!r}",
                    hint="Use dicts with the 12 PAF fields (query_name, query_length, ..., mapq) and optional tags.",
                ) from exc
            row.update({"_line": number, "_offset": None, "_length": None})
            rows.append(row)
        return rows

    def records(self, doc: list[dict[str, Any]]) -> list[dict[str, object]]:
        return [{**{k: r[k] for k in _PAF_COLUMNS}, "tags": self._tags(r)} for r in doc]

    def _line(self, r: dict[str, Any]) -> str:
        return "\t".join([str(r[k]) for k in _PAF_COLUMNS] + list(r["_tags"]))

    def dump(self, doc: list[dict[str, Any]], opts: Opts) -> Payload:
        return _text_payload("paf", [self._line(r) for r in doc], len(doc))

    def validate_doc(
        self, doc: list[dict[str, Any]], opts: Opts, errors: list[str], warnings: list[str]
    ) -> dict[str, object]:
        return {"alignment_count": len(doc)}

    @staticmethod
    def _overlap(name: str, start: int, end: int, region: tuple[str, int, int]) -> bool:
        return name == region[0] and start < region[2] and end > region[1]

    def filter(self, doc: list[dict[str, Any]], opts: Opts) -> list[dict[str, Any]]:
        min_mapq, min_matches = (
            opts.integer("min_mapq", None, 0, 255),
            opts.integer("min_matches", None, 0),
        )
        min_block = opts.integer("min_block_length", None, 0)
        min_id, min_cov = (
            opts.number("min_identity", None, 0.0, 1.0),
            opts.number("min_query_coverage", None, 0.0, 1.0),
        )
        strand = opts.text("strand", None, ["+", "-"])
        qn, qx = opts.name_set("query_names"), opts.name_set("exclude_query_names")
        tn, tx = opts.name_set("target_names"), opts.name_set("exclude_target_names")
        treg = parse_region(opts.text("target_region")) if opts.has("target_region") else None
        qreg = parse_region(opts.text("query_region")) if opts.has("query_region") else None
        tags = opts.mapping("tags")
        wanted = {k: str(v) for k, v in (tags or {}).items()}
        out = []
        for r in doc:
            if min_mapq is not None and r["mapq"] < min_mapq:
                continue
            if min_matches is not None and r["matches"] < min_matches:
                continue
            if min_block is not None and r["block_length"] < min_block:
                continue
            if min_id is not None and _identity(r) < min_id:
                continue
            if min_cov is not None and _coverage(r) < min_cov:
                continue
            if strand is not None and r["strand"] != strand:
                continue
            if (qn is not None and r["query_name"] not in qn) or (
                qx is not None and r["query_name"] in qx
            ):
                continue
            if (tn is not None and r["target_name"] not in tn) or (
                tx is not None and r["target_name"] in tx
            ):
                continue
            if treg and not self._overlap(
                r["target_name"], r["target_start"], r["target_end"], treg
            ):
                continue
            if qreg and not self._overlap(r["query_name"], r["query_start"], r["query_end"], qreg):
                continue
            if wanted:
                have = {t: str(v["value"]) for t, v in self._tags(r).items()}
                if any(have.get(t) != v for t, v in wanted.items()):
                    continue
            out.append(r)
        return out

    _SORT_KEYS = ("target", "query", "mapq", "identity", "block_length", "matches")

    def _sorted(self, doc: list[dict[str, Any]], key: str, reverse: bool) -> list[dict[str, Any]]:
        keys = {
            "target": lambda r: (r["target_name"], r["target_start"], r["target_end"]),
            "query": lambda r: (r["query_name"], r["query_start"], r["query_end"]),
            "mapq": lambda r: r["mapq"],
            "identity": _identity,
            "block_length": lambda r: r["block_length"],
            "matches": lambda r: r["matches"],
        }
        return sorted(doc, key=keys[key], reverse=reverse)

    def sort(self, doc: list[dict[str, Any]], opts: Opts) -> list[dict[str, Any]]:
        return self._sorted(
            doc,
            str(opts.text("key", "target", list(self._SORT_KEYS))),
            opts.boolean("reverse", False),
        )

    def merge(self, docs: list[list[dict[str, Any]]], opts: Opts) -> list[dict[str, Any]]:
        dedupe = opts.boolean("deduplicate", False)
        sort_by = opts.text("sort_by", None, list(self._SORT_KEYS))
        out, seen = [], set()
        for doc in docs:
            for r in doc:
                if dedupe:
                    line = self._line(r)
                    if line in seen:
                        continue
                    seen.add(line)
                out.append(r)
        return self._sorted(out, sort_by, False) if sort_by else out

    def split(
        self, doc: list[dict[str, Any]], opts: Opts
    ) -> list[tuple[str, list[dict[str, Any]]]]:
        by = opts.text("by", "query", ["query", "target", "strand", "size"])
        size = opts.integer("size", None, 1)
        if by == "size":
            if size is None:
                raise bad_option("split_paf with by='size' needs the size option.")
            return [(f"{i // size + 1:04d}", doc[i : i + size]) for i in range(0, len(doc), size)]
        field = {"query": "query_name", "target": "target_name", "strand": "strand"}[str(by)]
        groups: dict[str, list[dict[str, Any]]] = {}
        for r in doc:
            groups.setdefault(r[field], []).append(r)
        return list(groups.items())

    def summarize(self, doc: list[dict[str, Any]], opts: Opts) -> dict[str, object]:
        block = sum(r["block_length"] for r in doc)
        matches = sum(r["matches"] for r in doc)
        ids = [_identity(r) for r in doc]
        return {
            "alignment_count": len(doc),
            "query_count": len({r["query_name"] for r in doc}),
            "target_count": len({r["target_name"] for r in doc}),
            "strand": {s: sum(1 for r in doc if r["strand"] == s) for s in ("+", "-")},
            "aligned_bases": block,
            "matching_bases": matches,
            "overall_identity": matches / block if block else None,
            "mean_identity": sum(ids) / len(ids) if ids else None,
            "min_identity": min(ids, default=None),
            "max_identity": max(ids, default=None),
            "mapq": {str(k): v for k, v in sorted(Counter(r["mapq"] for r in doc).items())},
            "mean_mapq": sum(r["mapq"] for r in doc) / len(doc) if doc else None,
        }

    def index(self, doc: list[dict[str, Any]], opts: Opts) -> dict[str, object]:
        by = opts.text("by", "query", ["query", "target"])
        name_f, start_f, end_f = (
            ("query_name", "query_start", "query_end")
            if by == "query"
            else ("target_name", "target_start", "target_end")
        )
        entries: dict[str, list[dict[str, int]]] = {}
        for r in doc:
            entries.setdefault(r[name_f], []).append(
                {
                    "line": r["_line"],
                    "offset": r["_offset"],
                    "length": r["_length"],
                    "start": r[start_f],
                    "end": r[end_f],
                }
            )
        for items in entries.values():
            items.sort(key=lambda e: (e["start"], e["end"], e["line"]))
        return {"key": f"{by}_name", "record_count": len(doc), "entries": entries}

    def convert(self, doc: list[dict[str, Any]], target: str, opts: Opts) -> Payload:
        if target == "paf":
            return self.dump(doc, opts)
        if target == "bed":
            return _text_payload(
                "bed",
                [
                    f"{r['target_name']}\t{r['target_start']}\t{r['target_end']}\t{r['query_name']}\t{r['mapq']}\t{r['strand']}"
                    for r in doc
                ],
                len(doc),
            )
        if target == "tsv":
            tags = sorted({t for r in doc for t in self._tags(r)})
            lines = ["\t".join(_PAF_COLUMNS + tuple(tags))]
            for r in doc:
                have = self._tags(r)
                lines.append(
                    "\t".join(
                        [str(r[k]) for k in _PAF_COLUMNS]
                        + [str(have[t]["value"]) if t in have else "" for t in tags]
                    )
                )
            return _text_payload("tsv", lines, len(doc))
        return json_payload("paf", self.records(doc), len(doc))
