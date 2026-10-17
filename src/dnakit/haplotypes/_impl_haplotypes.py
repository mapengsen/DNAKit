"""Native implementations behind the haplotype candidate functions.

Pure standard library, deterministic, bounded.  Genotypes are lists of calls; a call is
a VCF-style string (``"0/1"``, ``"1|0"``, ``"./."``), a pair ``(0, 1)``, an alt-allele
dosage ``0/1/2`` (biallelic) or ``None`` for missing.  Haplotypes are strings (or lists)
of per-site alleles, usually ``"0"``/``"1"`` characters.
"""

from __future__ import annotations

import itertools
import math
import random
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_SITES = 200_000
MAX_HAPLOTYPES = 100_000
MAX_LD_SITES = 2_000
MAX_EM_HET = 12
MAX_PEDIGREE_UNOBSERVED = 8



def options_dict(options: Mapping[str, object] | None) -> dict[str, object]:
    return dict(options or {})


def get_int(settings: Mapping[str, object], key: str, default: int | None, *, low: int, high: int) -> int | None:
    value = settings.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"options.{key} must be an integer in [{low}, {high}].", code="HAPLOTYPE_OPTION", hint=f"Set options['{key}'] to an integer in [{low}, {high}].")
    return value


def get_float(settings: Mapping[str, object], key: str, default: float, *, low: float, high: float) -> float:
    value = settings.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
        raise ConfigurationError(f"options.{key} must be a number in [{low}, {high}].", code="HAPLOTYPE_OPTION")
    return float(value)


def write_output(settings: Mapping[str, object], text: str) -> str | None:
    target = settings.get("output")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or "\x00" in str(target):
        raise ConfigurationError("options.output must be a file path.", code="HAPLOTYPE_OPTION")
    path = Path(target)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- genotypes


def parse_call(value: object, where: str = "genotype") -> tuple[tuple[int, ...], bool] | None:
    """Parse one diploid call into ``(alleles, phased)``; ``None`` means missing."""

    if value is None:
        return None
    if isinstance(value, bool):
        raise InputFormatError(f"{where}: a boolean is not a genotype.", code="GENOTYPE_FORMAT")
    if isinstance(value, int):
        if value not in (0, 1, 2):
            raise InputFormatError(f"{where}: integer genotypes are alt-allele dosages 0, 1 or 2; got {value}.", code="GENOTYPE_FORMAT")
        return ((0, 0), False) if value == 0 else ((0, 1), False) if value == 1 else ((1, 1), False)
    if isinstance(value, str):
        text = value.strip()
        if text in ("", ".", "./.", ".|.", "NA", "na", "NaN"):
            return None
        parts = re.split(r"([/|])", text)
        alleles = parts[::2]
        seps = set(parts[1::2])
        if len(alleles) != 2 or not all(a.isdigit() for a in alleles) or len(seps) != 1:
            if any(a == "." for a in alleles):
                return None
            raise InputFormatError(f"{where}: cannot parse genotype {value!r} (expected 'a/b' or 'a|b' with integer alleles).", code="GENOTYPE_FORMAT")
        return (int(alleles[0]), int(alleles[1])), seps == {"|"}
    if isinstance(value, Sequence) and len(value) == 2 and all(isinstance(a, int) and not isinstance(a, bool) and a >= 0 for a in value):
        return (int(value[0]), int(value[1])), False
    raise InputFormatError(f"{where}: unsupported genotype {value!r}.", code="GENOTYPE_FORMAT")


def parse_calls(value: object, name: str) -> list[tuple[tuple[int, ...], bool] | None]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise InputFormatError(f"{name} must be a list of genotype calls.", code="GENOTYPE_FORMAT")
    if len(value) > MAX_SITES:
        raise InputFormatError(f"{name} has more than {MAX_SITES} sites.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    return [parse_call(v, f"{name}[{i}]") for i, v in enumerate(value)]


def trio_calls(input: object, *, need_parents: bool = True) -> dict[str, Any]:
    if not isinstance(input, Mapping) or "child" not in input:
        raise InputFormatError("trio input must be a mapping with 'child', 'father' and 'mother' genotype lists.", code="TRIO_INPUT", hint="Pass {'child': [...], 'father': [...], 'mother': [...]}.")
    child = parse_calls(input["child"], "child")
    out: dict[str, Any] = {"child": child}
    for role in ("father", "mother"):
        if role in input:
            out[role] = parse_calls(input[role], role)
            if len(out[role]) != len(child):
                raise InputFormatError(f"{role} has {len(out[role])} sites but child has {len(child)}.", code="TRIO_INPUT")
        elif need_parents:
            raise InputFormatError(f"trio input is missing '{role}'.", code="TRIO_INPUT")
    if not child:
        raise InputFormatError("trio input has no sites.", code="TRIO_INPUT")
    labels = input.get("sites", input.get("positions", input.get("variants")))
    if labels is not None:
        if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence) or len(labels) != len(child):
            raise InputFormatError("'sites' must list one label per site.", code="TRIO_INPUT")
        out["labels"] = list(labels)
    return out


def assignments(child: tuple[int, ...], father: tuple[int, ...], mother: tuple[int, ...]) -> set[tuple[int, int]]:
    """All (paternal_allele, maternal_allele) pairs consistent with the child's genotype."""

    found: set[tuple[int, int]] = set()
    for a in set(father):
        for b in set(mother):
            if sorted((a, b)) == sorted(child):
                found.add((a, b))
    return found


def trio_inheritance(input: object) -> dict[str, Any]:
    trio = trio_calls(input)
    sites: list[dict[str, Any]] = []
    counts = Counter()
    bad: list[int] = []
    for i, (c, f, m) in enumerate(zip(trio["child"], trio["father"], trio["mother"])):
        entry: dict[str, Any] = {"index": i}
        if "labels" in trio:
            entry["site"] = trio["labels"][i]
        if c is None or f is None or m is None:
            entry["status"] = "missing"
        else:
            options = assignments(c[0], f[0], m[0])
            if not options:
                entry["status"] = "inconsistent"
                bad.append(i)
            else:
                entry["status"] = "consistent"
                if len(options) == 1:
                    (pat, mat) = next(iter(options))
                    entry["paternal_allele"], entry["maternal_allele"] = pat, mat
        counts[entry["status"]] += 1
        sites.append(entry)
    informative = counts["consistent"] + counts["inconsistent"]
    return {
        "site_count": len(sites),
        "informative_sites": informative,
        "consistent": counts["consistent"],
        "inconsistent": counts["inconsistent"],
        "missing": counts["missing"],
        "mendelian_error_rate": counts["inconsistent"] / informative if informative else None,
        "inconsistent_indices": bad,
        "sites": sites,
        "rule": "a site is consistent if the child's two alleles can be formed by one allele from each parent",
    }


def denovo_sites(input: object) -> dict[str, Any]:
    trio = trio_calls(input)
    denovo: list[dict[str, Any]] = []
    inconsistent: list[int] = []
    evaluated = 0
    for i, (c, f, m) in enumerate(zip(trio["child"], trio["father"], trio["mother"])):
        if c is None or f is None or m is None:
            continue
        evaluated += 1
        parental = set(f[0]) | set(m[0])
        novel = sorted(set(c[0]) - parental)
        if novel:
            entry = {"index": i, "child_genotype": "/".join(map(str, c[0])), "novel_alleles": novel}
            if "labels" in trio:
                entry["site"] = trio["labels"][i]
            denovo.append(entry)
        elif not assignments(c[0], f[0], m[0]):
            inconsistent.append(i)
    return {
        "evaluated_sites": evaluated,
        "denovo_count": len(denovo),
        "denovo": denovo,
        "other_mendelian_inconsistent_indices": inconsistent,
        "rule": "de novo = the child carries an allele that neither parent's genotype contains; other Mendelian violations (e.g. missing parental allele transmitted) are listed separately",
    }


def compound_hets(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    if not isinstance(input, Mapping) or "genes" not in input:
        raise InputFormatError("input must be a mapping with 'genes' (gene per site) and genotype lists.", code="COMPOUND_HET_INPUT", hint="Pass {'genes': [...], 'child': [...], 'father': [...], 'mother': [...]} or phased child calls only.")
    genes = input["genes"]
    if isinstance(genes, (str, bytes)) or not isinstance(genes, Sequence):
        raise InputFormatError("'genes' must be a list.", code="COMPOUND_HET_INPUT")
    has_parents = "father" in input and "mother" in input
    trio = trio_calls(input, need_parents=has_parents)
    if len(genes) != len(trio["child"]):
        raise InputFormatError("'genes' must list one gene per site.", code="COMPOUND_HET_INPUT")
    per_gene: dict[str, dict[str, list[int]]] = {}
    unresolved = 0
    for i, c in enumerate(trio["child"]):
        if c is None or c[0][0] == c[0][1]:
            continue  # missing or homozygous: not a heterozygous call
        entry = per_gene.setdefault(str(genes[i]), {"paternal": [], "maternal": [], "undetermined": []})
        if has_parents:
            f, m = trio["father"][i], trio["mother"][i]
            origin = assignments(c[0], f[0], m[0]) if f is not None and m is not None else set()
            if len(origin) != 1:
                entry["undetermined"].append(i)
                unresolved += 1
                continue
            pat, mat = next(iter(origin))
            first, second = pat, mat
        else:
            if not c[1]:
                raise InputFormatError(f"site {i}: child genotype is unphased and no parents were given.", code="COMPOUND_HET_INPUT", hint="Provide phased child calls ('0|1') or both parents.")
            first, second = c[0]
        if first != 0:
            entry["paternal"].append(i)
        if second != 0:
            entry["maternal"].append(i)
    min_variants = get_int(settings, "min_variants_per_haplotype", 1, low=1, high=1000) or 1
    hits = []
    for gene, entry in per_gene.items():
        if len(entry["paternal"]) >= min_variants and len(entry["maternal"]) >= min_variants:
            hits.append({"gene": gene, **entry})
    mode = "trio" if has_parents else "phased_child"
    labels = "paternal/maternal" if has_parents else "haplotype_1/haplotype_2"
    return {
        "mode": mode,
        "haplotype_labels": labels,
        "genes_with_het_variants": len(per_gene),
        "compound_heterozygous_genes": hits,
        "compound_heterozygous_count": len(hits),
        "undetermined_origin_sites": unresolved,
        "rule": "compound heterozygote = a gene with a heterozygous variant on each of the two parental haplotypes (in trans); trio mode assigns origin when only one parent carries the alt allele, phased mode uses the phase of the child's calls",
    }


def trio_phase(input: object) -> dict[str, Any]:
    trio = trio_calls(input)
    phased: list[str | None] = []
    status: list[str] = []
    for c, f, m in zip(trio["child"], trio["father"], trio["mother"]):
        if c is None:
            phased.append(None)
            status.append("missing")
            continue
        if f is None or m is None:
            phased.append("|".join(map(str, c[0])) if c[0][0] == c[0][1] else "/".join(map(str, sorted(c[0]))))
            status.append("homozygous" if c[0][0] == c[0][1] else "unresolved_parent_missing")
            continue
        options = assignments(c[0], f[0], m[0])
        if not options:
            phased.append("/".join(map(str, sorted(c[0]))))
            status.append("mendelian_inconsistent")
        elif len(options) == 1:
            pat, mat = next(iter(options))
            phased.append(f"{pat}|{mat}")
            status.append("homozygous" if pat == mat else "phased")
        else:
            phased.append("/".join(map(str, sorted(c[0]))))
            status.append("ambiguous")
    counts = Counter(status)
    het_total = counts["phased"] + counts["ambiguous"] + counts["unresolved_parent_missing"]
    return {
        "phased_child": phased,
        "status": status,
        "phased_het_sites": counts["phased"],
        "ambiguous_het_sites": counts["ambiguous"],
        "inconsistent_sites": counts["mendelian_inconsistent"],
        "fraction_het_sites_phased": counts["phased"] / het_total if het_total else None,
        "convention": "paternal|maternal",
        "rule": "Mendelian trio phasing: a heterozygous child site is phased when exactly one assignment of its two alleles to (father, mother) is compatible with the parents' genotypes",
    }


def infer_parentage(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    if not isinstance(input, Mapping) or "child" not in input or not isinstance(input.get("candidates"), Mapping) or not input["candidates"]:
        raise InputFormatError("input must be {'child': [...], 'candidates': {name: [...]}}.", code="PARENTAGE_INPUT")
    child = parse_calls(input["child"], "child")
    max_mismatches = get_int(settings, "max_mismatches", 0, low=0, high=10**9) or 0
    report: dict[str, dict[str, Any]] = {}
    for name, calls in input["candidates"].items():
        cand = parse_calls(calls, str(name))
        if len(cand) != len(child):
            raise InputFormatError(f"candidate {name!r} has {len(cand)} sites but the child has {len(child)}.", code="PARENTAGE_INPUT")
        compared = mismatches = 0
        for c, p in zip(child, cand):
            if c is None or p is None:
                continue
            compared += 1
            if not (set(c[0]) & set(p[0])):  # a parent must share at least one allele with its child
                mismatches += 1
        report[str(name)] = {
            "compared_sites": compared,
            "mismatches": mismatches,
            "mismatch_rate": mismatches / compared if compared else None,
            "excluded": mismatches > max_mismatches,
        }
    compatible = {n: r for n, r in report.items() if not r["excluded"] and r["compared_sites"] > 0}
    best_count = min((r["mismatches"] for r in compatible.values()), default=None)
    best = sorted(n for n, r in compatible.items() if r["mismatches"] == best_count) if compatible else []
    not_assessed = sorted(n for n, r in report.items() if r["compared_sites"] == 0)
    warnings = []
    if not_assessed:
        warnings.append(f"{len(not_assessed)} candidate(s) have no jointly called site with the child and could not be assessed (neither excluded nor "
                        f"not_excluded): {', '.join(not_assessed)}")
    return {
        "candidates": report,
        "not_excluded": sorted(compatible),
        "not_assessed": not_assessed,
        "warnings": warnings,
        "best_candidates": best,
        "max_mismatches": max_mismatches,
        "rule": "exclusion method: a candidate parent is excluded when the number of sites where parent and child share no allele (e.g. opposite homozygotes) exceeds max_mismatches; a candidate with no jointly called site cannot be assessed and is listed in not_assessed (never in not_excluded or best_candidates)",
    }


# --------------------------------------------------------------------------- pedigree likelihood


def _dosage(call: tuple[tuple[int, ...], bool] | None, where: str) -> int | None:
    if call is None:
        return None
    if any(a not in (0, 1) for a in call[0]):
        raise InputFormatError(f"{where}: only biallelic (0/1) genotypes are supported.", code="PEDIGREE_INPUT")
    return sum(call[0])


def pedigree_likelihood(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    if not isinstance(input, Mapping) or not isinstance(input.get("pedigree"), Sequence) or not isinstance(input.get("genotypes"), Mapping):
        raise InputFormatError("input must be {'pedigree': [{'id','father','mother'}, ...], 'genotypes': {id: [calls]}}.", code="PEDIGREE_INPUT")
    members: dict[str, tuple[str | None, str | None]] = {}
    for row in input["pedigree"]:  # type: ignore[union-attr]
        if not isinstance(row, Mapping) or "id" not in row:
            raise InputFormatError("each pedigree entry needs an 'id'.", code="PEDIGREE_INPUT")
        father, mother = row.get("father"), row.get("mother")
        father = None if father in (None, "", "0", 0) else str(father)
        mother = None if mother in (None, "", "0", 0) else str(mother)
        if (father is None) != (mother is None):
            raise InputFormatError(f"{row['id']}: give both parents or neither.", code="PEDIGREE_INPUT")
        members[str(row["id"])] = (father, mother)
    order: list[str] = []
    remaining = dict(members)
    while remaining:
        ready = [i for i, (f, m) in remaining.items() if (f is None or f not in remaining) and (m is None or m not in remaining)]
        if not ready:
            raise InputFormatError("the pedigree contains a cycle.", code="PEDIGREE_INPUT")
        for i in ready:
            order.append(i)
            del remaining[i]
    for i, (f, m) in members.items():
        for parent in (f, m):
            if parent is not None and parent not in members:
                raise InputFormatError(f"{i}: parent {parent!r} is not listed in the pedigree.", code="PEDIGREE_INPUT")
    genotypes = input["genotypes"]
    site_count = None
    calls: dict[str, list[int | None]] = {}
    for person in members:
        raw = genotypes.get(person, [])  # type: ignore[union-attr]
        parsed = [_dosage(c, f"{person}[{i}]") for i, c in enumerate(parse_calls(raw, person))] if raw else []
        calls[person] = parsed  # type: ignore[assignment]
        if parsed:
            if site_count is not None and len(parsed) != site_count:
                raise InputFormatError("all genotype lists must have the same length.", code="PEDIGREE_INPUT")
            site_count = len(parsed)
    if not site_count:
        raise InputFormatError("no genotypes given.", code="PEDIGREE_INPUT")
    freqs_raw = settings.get("allele_frequency", 0.5)
    if isinstance(freqs_raw, (int, float)) and not isinstance(freqs_raw, bool):
        freqs = [float(freqs_raw)] * site_count
    elif isinstance(freqs_raw, Sequence) and not isinstance(freqs_raw, (str, bytes)) and len(freqs_raw) == site_count:
        freqs = [float(x) for x in freqs_raw]
    else:
        raise ConfigurationError("options.allele_frequency must be a number or one number per site.", code="HAPLOTYPE_OPTION")
    if any(not 0.0 < p < 1.0 for p in freqs):
        raise ConfigurationError("allele frequencies must be in (0, 1).", code="HAPLOTYPE_OPTION")

    per_site: list[float] = []
    for s in range(site_count):
        p = freqs[s]
        observed = {person: (calls[person][s] if calls[person] else None) for person in members}
        unknown = [person for person in order if observed[person] is None]
        if len(unknown) > MAX_PEDIGREE_UNOBSERVED:
            raise InputFormatError(f"site {s}: more than {MAX_PEDIGREE_UNOBSERVED} individuals are ungenotyped; brute-force peeling is limited.", code="PEDIGREE_TOO_LARGE")
        total = 0.0
        for fill in itertools.product((0, 1, 2), repeat=len(unknown)):
            state = dict(observed)
            state.update(zip(unknown, fill))
            probability = 1.0
            for person in order:
                f, m = members[person]
                g = state[person]
                if f is None:
                    probability *= ((1 - p) ** 2, 2 * p * (1 - p), p * p)[g]  # type: ignore[index]
                else:
                    tf, tm = state[f] / 2.0, state[m] / 2.0  # type: ignore[index, operator]
                    probability *= ((1 - tf) * (1 - tm), tf * (1 - tm) + (1 - tf) * tm, tf * tm)[g]  # type: ignore[index]
                if probability == 0.0:
                    break
            total += probability
        per_site.append(total)
    zero = any(x == 0.0 for x in per_site)
    return {
        "site_likelihoods": per_site,
        "likelihood": math.prod(per_site),
        "log10_likelihood": None if zero else sum(math.log10(x) for x in per_site),
        "impossible": zero,
        "individuals": len(members),
        "sites": site_count,
        "model": "founders in Hardy-Weinberg proportions with the given alt-allele frequency, Mendelian transmission, no genotyping error, sites independent; ungenotyped members summed out exactly",
    }


# --------------------------------------------------------------------------- haplotype strings


def haplotype_list(input: object, *, min_count: int = 1) -> list[tuple[str, list[str]]]:
    """Return (name, [haplotype strings]) pairs; names are sample names or ``hap_i``."""

    if isinstance(input, Mapping):
        items = list(input.items())
        if items and all(isinstance(v, str) for _, v in items):
            return [(str(k), [v]) for k, v in items]
        if items and all(isinstance(v, Sequence) and not isinstance(v, (str, bytes)) and all(isinstance(x, str) for x in v) for _, v in items):
            return [(str(k), list(v)) for k, v in items]
        raise InputFormatError("haplotype mapping must map names to a haplotype or a list of haplotypes.", code="HAPLOTYPE_INPUT")
    if isinstance(input, Sequence) and not isinstance(input, (str, bytes)):
        if not all(isinstance(x, str) for x in input):
            raise InputFormatError("haplotypes must be strings of alleles.", code="HAPLOTYPE_INPUT")
        return [(f"hap_{i + 1}", [x]) for i, x in enumerate(input)]
    raise InputFormatError("haplotypes must be a list of strings or a {name: haplotype(s)} mapping.", code="HAPLOTYPE_INPUT")


def flat_haplotypes(input: object) -> list[str]:
    haps = [h for _, hs in haplotype_list(input) for h in hs]
    if not haps:
        raise InputFormatError("no haplotypes given.", code="HAPLOTYPE_INPUT")
    if len(haps) > MAX_HAPLOTYPES:
        raise InputFormatError(f"more than {MAX_HAPLOTYPES} haplotypes.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    return haps


def frequency_table(input: object) -> dict[str, Any]:
    if isinstance(input, Mapping) and input and all(isinstance(v, int) and not isinstance(v, bool) for v in input.values()):
        if any(v < 0 for v in input.values()):
            raise InputFormatError("haplotype counts must be non-negative.", code="HAPLOTYPE_INPUT")
        counts = Counter({str(k): v for k, v in input.items() if v > 0})
    else:
        counts = Counter(flat_haplotypes(input))
    total = sum(counts.values())
    if total == 0:
        raise InputFormatError("no haplotypes given.", code="HAPLOTYPE_INPUT")
    rows = [{"haplotype": h, "count": c, "frequency": c / total} for h, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    diversity = (total / (total - 1)) * (1 - sum((c / total) ** 2 for c in counts.values())) if total > 1 else None
    return {
        "total": total,
        "distinct": len(rows),
        "haplotypes": rows,
        "haplotype_diversity": diversity,
        "definition": "frequency = count / total chromosomes; haplotype diversity (Nei) = n/(n-1) * (1 - sum p_i^2)",
    }


def _check_equal_length(haps: Sequence[str]) -> int:
    lengths = {len(h) for h in haps}
    if len(lengths) != 1:
        raise InputFormatError("all haplotypes must have the same number of sites.", code="HAPLOTYPE_INPUT")
    return lengths.pop()


_MISSING = frozenset("?N.-")


def hamming(a: str, b: str, ignore_missing: bool = True) -> tuple[int, int]:
    """Return (differences, compared_sites)."""

    diff = compared = 0
    for x, y in zip(a, b):
        if ignore_missing and (x in _MISSING or y in _MISSING):
            continue
        compared += 1
        diff += x != y
    return diff, compared


def cluster_haplotypes(input: object, max_distance: int, ignore_missing: bool) -> dict[str, Any]:
    entries = [(name, h) for name, hs in haplotype_list(input) for h in hs]
    if not entries:
        raise InputFormatError("no haplotypes given.", code="HAPLOTYPE_INPUT")
    if len(entries) > 5_000:
        raise InputFormatError("clustering is limited to 5,000 haplotypes.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    _check_equal_length([h for _, h in entries])
    parent = list(range(len(entries)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(entries)):
        for j in range(i + 1, len(entries)):
            if hamming(entries[i][1], entries[j][1], ignore_missing)[0] <= max_distance:
                parent[find(j)] = find(i)
    groups: dict[int, list[int]] = {}
    for i in range(len(entries)):
        groups.setdefault(find(i), []).append(i)
    clusters = []
    for members in sorted(groups.values(), key=lambda m: m[0]):
        clusters.append({"size": len(members), "members": [entries[i][0] for i in members], "haplotypes": sorted({entries[i][1] for i in members}), "indices": members})
    return {"cluster_count": len(clusters), "clusters": clusters, "max_distance": max_distance, "linkage": "single (connected components of the graph joining haplotypes within max_distance mismatches)"}


def compare_haplotypes(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    mode = str(settings.get("mode", "hamming"))
    if mode == "switch_error":
        if not isinstance(input, Mapping) or "truth" not in input or "inferred" not in input:
            raise InputFormatError("switch_error mode needs {'truth': [h1, h2], 'inferred': [h1, h2]}.", code="HAPLOTYPE_INPUT")
        truth, inferred = input["truth"], input["inferred"]
        for pair in (truth, inferred):
            if not isinstance(pair, Sequence) or isinstance(pair, (str, bytes)) or len(pair) != 2 or not all(isinstance(x, str) for x in pair):
                raise InputFormatError("each of truth/inferred must be two haplotype strings.", code="HAPLOTYPE_INPUT")
        n = _check_equal_length([*truth, *inferred])
        orientation: list[bool] = []
        for i in range(n):
            t1, t2, i1, i2 = truth[0][i], truth[1][i], inferred[0][i], inferred[1][i]
            if t1 == t2 or i1 == i2 or {t1, t2} != {i1, i2}:
                continue
            orientation.append(t1 == i1)
        switches = sum(1 for a, b in zip(orientation, orientation[1:]) if a != b)
        return {
            "mode": "switch_error",
            "heterozygous_sites": len(orientation),
            "switches": switches,
            "switch_error_rate": switches / (len(orientation) - 1) if len(orientation) > 1 else None,
            "definition": "among heterozygous sites called identically in truth and inferred, a switch is a change of relative phase between consecutive sites; rate = switches / (sites - 1)",
        }
    if mode != "hamming":
        raise ConfigurationError("options.mode must be 'hamming' or 'switch_error'.", code="HAPLOTYPE_OPTION")
    ignore_missing = bool(settings.get("ignore_missing", True))
    entries = [(name if len(hs) == 1 else f"{name}_{k + 1}", h) for name, hs in haplotype_list(input) for k, h in enumerate(hs)]
    if len(entries) < 2:
        raise InputFormatError("compare needs at least two haplotypes.", code="HAPLOTYPE_INPUT")
    if len(entries) > 2_000:
        raise InputFormatError("comparison is limited to 2,000 haplotypes.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    _check_equal_length([h for _, h in entries])
    names = [n for n, _ in entries]
    matrix = [[hamming(a, b, ignore_missing)[0] for _, b in entries] for _, a in entries]
    result: dict[str, Any] = {"mode": "hamming", "names": names, "distance_matrix": matrix, "ignore_missing": ignore_missing}
    if len(entries) == 2:
        a, b = entries[0][1], entries[1][1]
        result["differing_positions"] = [i for i, (x, y) in enumerate(zip(a, b)) if x != y and not (ignore_missing and (x in _MISSING or y in _MISSING))]
        diff, compared = hamming(a, b, ignore_missing)
        result["distance"], result["compared_sites"] = diff, compared
    return result


# --------------------------------------------------------------------------- LD, blocks, tag SNPs


def allele_matrix(input: object) -> list[list[int]]:
    """Rows = haplotypes, columns = biallelic sites (0/1)."""

    if isinstance(input, Mapping):
        input = [v for _, hs in haplotype_list(input) for v in hs]
    if isinstance(input, (str, bytes)) or not isinstance(input, Sequence) or not input:
        raise InputFormatError("expected a list of haplotypes (strings of 0/1 or lists of 0/1).", code="HAPLOTYPE_INPUT")
    rows: list[list[int]] = []
    for r, row in enumerate(input):
        if isinstance(row, str):
            cells = list(row)
        elif isinstance(row, Sequence):
            cells = [str(x) for x in row]
        else:
            raise InputFormatError(f"haplotype {r} must be a string or list.", code="HAPLOTYPE_INPUT")
        if any(c not in ("0", "1") for c in cells):
            raise InputFormatError(f"haplotype {r} contains alleles other than 0/1 (missing data is not supported here).", code="HAPLOTYPE_INPUT")
        rows.append([int(c) for c in cells])
    if len({len(r) for r in rows}) != 1 or not rows[0]:
        raise InputFormatError("all haplotypes must have the same, non-zero number of sites.", code="HAPLOTYPE_INPUT")
    if len(rows) < 2:
        raise InputFormatError("at least two haplotypes are needed.", code="HAPLOTYPE_INPUT")
    if len(rows[0]) > MAX_LD_SITES:
        raise InputFormatError(f"LD analysis is limited to {MAX_LD_SITES} sites.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    return rows


def ld_pair(x: Sequence[int], y: Sequence[int]) -> tuple[float, float]:
    """Return ``(r2, |D'|)`` for two biallelic columns."""

    n = len(x)
    pa = sum(x) / n
    pb = sum(y) / n
    pab = sum(1 for a, b in zip(x, y) if a and b) / n
    d = pab - pa * pb
    denom = pa * (1 - pa) * pb * (1 - pb)
    r2 = d * d / denom if denom > 0 else 0.0
    dmax = min(pa * (1 - pb), (1 - pa) * pb) if d > 0 else min(pa * pb, (1 - pa) * (1 - pb))
    dprime = abs(d) / dmax if dmax > 0 else 0.0
    return r2, dprime


def ld_context(rows: Sequence[Sequence[int]], min_maf: float) -> tuple[list[int], list[list[int]]]:
    n = len(rows)
    cols = list(zip(*rows))
    keep = []
    for j, col in enumerate(cols):
        p = sum(col) / n
        if 0 < p < 1 and min(p, 1 - p) >= min_maf:
            keep.append(j)
    return keep, [list(cols[j]) for j in keep]


def ld_blocks(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    rows = allele_matrix(input)
    statistic = str(settings.get("statistic", "dprime"))
    if statistic not in ("dprime", "r2"):
        raise ConfigurationError("options.statistic must be 'dprime' or 'r2'.", code="HAPLOTYPE_OPTION")
    threshold = get_float(settings, "threshold", 0.8, low=0.0, high=1.0)
    min_snps = get_int(settings, "min_snps", 2, low=2, high=10_000) or 2
    min_maf = get_float(settings, "min_maf", 0.0, low=0.0, high=0.5)
    keep, cols = ld_context(rows, min_maf)
    positions = settings.get("positions")
    if positions is not None and (isinstance(positions, (str, bytes)) or not isinstance(positions, Sequence) or len(positions) != len(rows[0])):
        raise ConfigurationError("options.positions must list one coordinate per site.", code="HAPLOTYPE_OPTION")
    m = len(cols)
    cache: dict[tuple[int, int], float] = {}

    def ld(a: int, b: int) -> float:
        key = (a, b)
        if key not in cache:
            r2, dp = ld_pair(cols[a], cols[b])
            cache[key] = dp if statistic == "dprime" else r2
        return cache[key]

    blocks: list[dict[str, Any]] = []
    i = 0
    while i < m:
        j = i
        while j + 1 < m and all(ld(k, j + 1) >= threshold - 1e-12 for k in range(i, j + 1)):
            j += 1
        if j - i + 1 >= min_snps:
            sites = keep[i : j + 1]
            block: dict[str, Any] = {"first_site": sites[0], "last_site": sites[-1], "sites": sites, "snp_count": len(sites)}
            if positions is not None:
                block["start_position"], block["end_position"] = positions[sites[0]], positions[sites[-1]]
            blocks.append(block)
            i = j + 1
        else:
            i += 1
    return {
        "blocks": blocks,
        "block_count": len(blocks),
        "statistic": statistic,
        "threshold": threshold,
        "sites_total": len(rows[0]),
        "sites_analysed": m,
        "sites_excluded_monomorphic_or_low_maf": [j for j in range(len(rows[0])) if j not in set(keep)],
        "rule": "greedy left-to-right: a block is a maximal run of >= min_snps polymorphic sites in which every pair has LD (|D'| or r^2) >= threshold; pairwise LD from the haplotype frequencies of the given chromosomes",
    }


def tag_snps(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    rows = allele_matrix(input)
    threshold = get_float(settings, "r2_threshold", 0.8, low=0.0, high=1.0)
    min_maf = get_float(settings, "min_maf", 0.0, low=0.0, high=0.5)
    keep, cols = ld_context(rows, min_maf)
    m = len(cols)
    if m == 0:
        raise InputFormatError("no polymorphic sites to tag.", code="HAPLOTYPE_INPUT")
    neighbours: list[set[int]] = [{i} for i in range(m)]
    for a in range(m):
        for b in range(a + 1, m):
            if ld_pair(cols[a], cols[b])[0] >= threshold - 1e-12:
                neighbours[a].add(b)
                neighbours[b].add(a)
    uncovered = set(range(m))
    tags: list[int] = []
    bins: dict[int, list[int]] = {}
    while uncovered:
        best = max(sorted(uncovered), key=lambda i: (len(neighbours[i] & uncovered), -i))
        newly = sorted(neighbours[best] & uncovered)
        tags.append(best)
        bins[best] = newly
        uncovered -= set(newly)
    return {
        "tag_sites": [keep[t] for t in tags],
        "bins": {str(keep[t]): [keep[j] for j in bins[t]] for t in tags},
        "tag_count": len(tags),
        "sites_analysed": m,
        "r2_threshold": threshold,
        "excluded_monomorphic_or_low_maf": [j for j in range(len(rows[0])) if j not in set(keep)],
        "rule": "greedy tagging (Carlson et al. 2004): repeatedly choose the site with r^2 >= threshold to the most still-untagged sites (ties: lowest index), assign those sites to its bin",
    }


# --------------------------------------------------------------------------- EM reconstruction


def reconstruct_em(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    if not isinstance(input, Mapping) or not input:
        raise InputFormatError("input must map sample names to genotype lists (alt dosage 0/1/2 or '0/1' strings).", code="EM_INPUT")
    max_iter = get_int(settings, "max_iterations", 1000, low=1, high=100_000) or 1000
    tol = get_float(settings, "tolerance", 1e-8, low=1e-15, high=1.0)
    people: dict[str, list[tuple[int, int]]] = {}
    n_sites = None
    for name, calls in input.items():
        parsed = parse_calls(calls, str(name))
        if n_sites is None:
            n_sites = len(parsed)
        if len(parsed) != n_sites or not parsed:
            raise InputFormatError("all samples need the same, non-zero number of sites.", code="EM_INPUT")
        genos = []
        for i, c in enumerate(parsed):
            if c is None:
                raise InputFormatError(f"{name}[{i}]: missing genotypes are not supported by the EM reconstruction.", code="EM_INPUT", hint="Remove sites with missing data first.")
            if any(a not in (0, 1) for a in c[0]):
                raise InputFormatError(f"{name}[{i}]: only biallelic genotypes are supported.", code="EM_INPUT")
            genos.append((c[0][0], c[0][1]))
        people[str(name)] = genos
    candidates: dict[str, list[tuple[str, str]]] = {}
    for name, genos in people.items():
        hets = [i for i, (a, b) in enumerate(genos) if a != b]
        if len(hets) > MAX_EM_HET:
            raise InputFormatError(f"{name}: more than {MAX_EM_HET} heterozygous sites; reconstruction is limited to {MAX_EM_HET}.", code="EM_TOO_LARGE")
        base = [str(a) for a, _ in genos]
        pairs: list[tuple[str, str]] = []
        for flips in itertools.product((0, 1), repeat=max(len(hets) - 1, 0)):
            h1, h2 = list(base), list(base)
            # first het site fixed on haplotype 1 to enumerate each unordered pair once
            assignment = (0, *flips) if hets else ()
            for site, flip in zip(hets, assignment):
                a, b = genos[site]
                h1[site], h2[site] = (str(a), str(b)) if not flip else (str(b), str(a))
            x, y = "".join(h1), "".join(h2)
            pairs.append((min(x, y), max(x, y)))
        candidates[name] = sorted(set(pairs))
    universe = sorted({h for pairs in candidates.values() for pair in pairs for h in pair})
    if len(universe) > 4_096:
        raise InputFormatError("too many candidate haplotypes.", code="EM_TOO_LARGE")
    freq = {h: 1.0 / len(universe) for h in universe}
    n_chromosomes = 2 * len(people)
    converged = False
    iterations = 0
    loglik = 0.0
    for iterations in range(1, max_iter + 1):
        counts = {h: 0.0 for h in universe}
        loglik = 0.0
        for pairs in candidates.values():
            weights = [(2.0 if a != b else 1.0) * freq[a] * freq[b] for a, b in pairs]
            total = sum(weights)
            loglik += math.log(total)
            for (a, b), w in zip(pairs, weights):
                post = w / total
                counts[a] += post
                counts[b] += post
        new = {h: counts[h] / n_chromosomes for h in universe}
        delta = max(abs(new[h] - freq[h]) for h in universe)
        freq = new
        if delta < tol:
            converged = True
            break
    assignments_out = {}
    for name, pairs in candidates.items():
        weights = [(2.0 if a != b else 1.0) * freq[a] * freq[b] for a, b in pairs]
        total = sum(weights)
        best = max(range(len(pairs)), key=lambda k: (weights[k], [-ord(c) for c in pairs[k][0]]))
        assignments_out[name] = {"haplotypes": list(pairs[best]), "posterior": weights[best] / total if total else None, "alternatives": len(pairs)}
    # report the log-likelihood at the FINAL frequencies (the E step above evaluated the frequencies before the last M step)
    loglik = 0.0
    for pairs in candidates.values():
        loglik += math.log(sum((2.0 if a != b else 1.0) * freq[a] * freq[b] for a, b in pairs))
    ordered = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))
    return {
        "frequencies": [{"haplotype": h, "frequency": f} for h, f in ordered if f > 1e-12],
        "individuals": assignments_out,
        "iterations": iterations,
        "converged": converged,
        "log_likelihood": loglik,
        "model": "Excoffier-Slatkin (1995) EM for haplotype frequencies from unphased biallelic genotypes under Hardy-Weinberg, uniform initial frequencies; individual haplotypes = maximum posterior pair",
    }


# --------------------------------------------------------------------------- selection, simulation, plot


def select_haplotypes(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    table = frequency_table(input)
    min_frequency = get_float(settings, "min_frequency", 0.05, low=0.0, high=1.0)
    min_count = get_int(settings, "min_count", 1, low=1, high=10**9) or 1
    top_n = get_int(settings, "top_n", None, low=1, high=10**9)
    chosen = [r for r in table["haplotypes"] if r["frequency"] >= min_frequency - 1e-12 and r["count"] >= min_count]
    if top_n is not None:
        chosen = chosen[:top_n]
    return {
        "selected": chosen,
        "selected_count": len(chosen),
        "total_distinct": table["distinct"],
        "covered_frequency": sum(r["frequency"] for r in chosen),
        "min_frequency": min_frequency,
        "min_count": min_count,
        "top_n": top_n,
    }


def simulate_haplotypes(input: object, settings: Mapping[str, object]) -> dict[str, Any]:
    if not isinstance(input, Mapping) or not input:
        raise InputFormatError("input must map haplotypes to non-negative frequencies or counts.", code="HAPLOTYPE_INPUT")
    n = get_int(settings, "n_individuals", 10, low=1, high=100_000) or 10
    seed = get_int(settings, "seed", 0, low=0, high=2**32 - 1) or 0
    haps = sorted(str(k) for k in input)
    weights = []
    for h in haps:
        w = input[h]
        if isinstance(w, bool) or not isinstance(w, (int, float)) or w < 0 or not math.isfinite(w):
            raise InputFormatError("haplotype frequencies must be non-negative numbers.", code="HAPLOTYPE_INPUT")
        weights.append(float(w))
    total = sum(weights)
    if total <= 0:
        raise InputFormatError("haplotype frequencies must not all be zero.", code="HAPLOTYPE_INPUT")
    rng = random.Random(seed)
    people = [rng.choices(haps, weights=weights, k=2) for _ in range(n)]
    observed = Counter(h for pair in people for h in pair)
    return {
        "individuals": [{"id": f"sim_{i + 1}", "haplotypes": pair} for i, pair in enumerate(people)],
        "observed_haplotype_counts": dict(sorted(observed.items())),
        "input_frequencies": {h: w / total for h, w in zip(haps, weights)},
        "n_individuals": n,
        "seed": seed,
        "model": "random mating: each of the 2 chromosomes of an individual is drawn independently from the haplotype frequencies (python random.Random(seed)); no mutation or recombination",
    }


def haplotype_svg(input: object, cell: int = 14) -> tuple[str, int, int]:
    entries = [(name if len(hs) == 1 else f"{name}_{k + 1}", h) for name, hs in haplotype_list(input) for k, h in enumerate(hs)]
    if not entries:
        raise InputFormatError("no haplotypes given.", code="HAPLOTYPE_INPUT")
    n_sites = _check_equal_length([h for _, h in entries])
    if len(entries) * n_sites > 200_000:
        raise InputFormatError("visualisation is limited to 200,000 cells.", code="HAPLOTYPE_INPUT_TOO_LARGE")
    label_width = 8 * max(len(name) for name, _ in entries) + 12
    width, height = label_width + n_sites * cell + 8, len(entries) * cell + 8
    colors = {"0": "#dfe7ef", "1": "#c0392b"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    for r, (name, hap) in enumerate(entries):
        y = 4 + r * cell
        parts.append(f'<text x="4" y="{y + cell - 3}" font-size="{cell - 3}">{escape(name)}</text>')
        for c, allele in enumerate(hap):
            fill = colors.get(allele, "#999999")
            parts.append(f'<rect x="{label_width + c * cell}" y="{y}" width="{cell - 1}" height="{cell - 1}" fill="{fill}"><title>{escape(name)} site {c}: {escape(allele)}</title></rect>')
    parts.append("</svg>")
    return "\n".join(parts), len(entries), n_sites
