"""Native clinical-genomics computations behind ``dnakit.clinical``.

Pure standard library. Everything here works on evidence the caller supplies (genotypes, read counts,
ACMG criteria); nothing queries a database or a model, and nothing here is a validated diagnostic.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from itertools import product
from typing import Any

from dnakit.variants._impl_analysis import binom_sf, site_name
from dnakit.variants._impl_prioritize import clinvar_terms, consequence_terms, field_value, gene_of
from dnakit.variants._impl_vcf import (
    VariantSet,
    build_hgvs_genomic,
    call_ad,
    call_depth,
    call_gt,
    diploid_class,
    fail,
    is_called,
    load_variant_set,
    number_or_none,
    option_bool,
    option_choice,
    option_int,
    option_number,
    record_view,
    reference_from_options,
    sample_qc,
    site_qc,
    split_multiallelic_set,
    write_output,
)

CLASSES = ("Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign")
_ABBREV = {"P": CLASSES[0], "LP": CLASSES[1], "VUS": CLASSES[2], "LB": CLASSES[3], "B": CLASSES[4]}


# ---------------------------------------------------------------------------
# ACMG/AMP 2015 evidence combination
# ---------------------------------------------------------------------------
_CODE_RE = re.compile(
    r"^(PVS|PS|PM|PP|BA|BS|BP)(\d)(?:[_\-\s]*(very[_\s-]?strong|strong|moderate|supporting|stand[_\s-]?alone))?$",
    re.IGNORECASE,
)
_VALID_NUMBERS = {
    "PVS": range(1, 2),
    "PS": range(1, 5),
    "PM": range(1, 7),
    "PP": range(1, 6),
    "BA": range(1, 2),
    "BS": range(1, 5),
    "BP": range(1, 8),
}
_DEFAULT_STRENGTH = {
    "PVS": "very_strong",
    "PS": "strong",
    "PM": "moderate",
    "PP": "supporting",
    "BA": "stand_alone",
    "BS": "strong",
    "BP": "supporting",
}
_PATHOGENIC_STRENGTHS = {"very_strong", "strong", "moderate", "supporting"}
_BENIGN_STRENGTHS = {"stand_alone", "strong", "supporting"}


def parse_criterion(item: object) -> dict[str, str]:
    """Normalise one ACMG criterion given as ``'PM2'``, ``'PM2_Supporting'`` or ``{'code', 'strength'}``."""

    strength = None
    if isinstance(item, Mapping):
        code = item.get("code")
        strength = item.get("strength")
        if item.get("met") is False:
            return {}
    else:
        code = item
    if not isinstance(code, str):
        fail(
            f"criterion {item!r} must be a code string or a mapping with 'code'.",
            "CLINICAL_BAD_CRITERION",
            "Use codes such as 'PVS1', 'PM2_Supporting', 'BP4'.",
        )
    match = _CODE_RE.match(code.strip())
    if not match:
        fail(
            f"unknown ACMG/AMP criterion {code!r}.",
            "CLINICAL_BAD_CRITERION",
            "Valid codes: PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7 (optional _Strength suffix).",
        )
    family, number = match.group(1).upper(), int(match.group(2))
    if number not in _VALID_NUMBERS[family]:
        fail(
            f"criterion {family}{number} does not exist in ACMG/AMP 2015.",
            "CLINICAL_BAD_CRITERION",
            "Valid codes: PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7.",
        )
    suffix = match.group(3) or strength
    if suffix is not None:
        suffix = re.sub(r"[\s-]+", "_", str(suffix).strip().lower())
        suffix = {"verystrong": "very_strong", "standalone": "stand_alone"}.get(suffix, suffix)
    resolved = suffix or _DEFAULT_STRENGTH[family]
    pathogenic = family.startswith("P")
    if resolved not in (_PATHOGENIC_STRENGTHS if pathogenic else _BENIGN_STRENGTHS):
        fail(
            f"strength {resolved!r} is not defined for {family}{number} in the 2015 combining table.",
            "CLINICAL_BAD_CRITERION",
            "Pathogenic: very_strong/strong/moderate/supporting; benign: stand_alone/strong/supporting.",
        )
    return {
        "code": f"{family}{number}",
        "strength": resolved,
        "direction": "pathogenic" if pathogenic else "benign",
    }


_PATHOGENIC_RULES = (
    ("Pathogenic (i)(a): 1 Very strong + >=1 Strong", lambda c: c["vs"] >= 1 and c["s"] >= 1),
    ("Pathogenic (i)(b): 1 Very strong + >=2 Moderate", lambda c: c["vs"] >= 1 and c["m"] >= 2),
    (
        "Pathogenic (i)(c): 1 Very strong + 1 Moderate + 1 Supporting",
        lambda c: c["vs"] >= 1 and c["m"] >= 1 and c["p"] >= 1,
    ),
    ("Pathogenic (i)(d): 1 Very strong + >=2 Supporting", lambda c: c["vs"] >= 1 and c["p"] >= 2),
    ("Pathogenic (ii): >=2 Strong", lambda c: c["s"] >= 2),
    ("Pathogenic (iii)(a): 1 Strong + >=3 Moderate", lambda c: c["s"] >= 1 and c["m"] >= 3),
    (
        "Pathogenic (iii)(b): 1 Strong + 2 Moderate + >=2 Supporting",
        lambda c: c["s"] >= 1 and c["m"] >= 2 and c["p"] >= 2,
    ),
    (
        "Pathogenic (iii)(c): 1 Strong + 1 Moderate + >=4 Supporting",
        lambda c: c["s"] >= 1 and c["m"] >= 1 and c["p"] >= 4,
    ),
)
_LIKELY_PATHOGENIC_RULES = (
    ("Likely pathogenic (i): 1 Very strong + 1 Moderate", lambda c: c["vs"] >= 1 and c["m"] >= 1),
    ("Likely pathogenic (ii): 1 Strong + 1-2 Moderate", lambda c: c["s"] >= 1 and c["m"] >= 1),
    ("Likely pathogenic (iii): 1 Strong + >=2 Supporting", lambda c: c["s"] >= 1 and c["p"] >= 2),
    ("Likely pathogenic (iv): >=3 Moderate", lambda c: c["m"] >= 3),
    ("Likely pathogenic (v): 2 Moderate + >=2 Supporting", lambda c: c["m"] >= 2 and c["p"] >= 2),
    ("Likely pathogenic (vi): 1 Moderate + >=4 Supporting", lambda c: c["m"] >= 1 and c["p"] >= 4),
)
_BENIGN_RULES = (
    ("Benign (i): 1 Stand-alone (BA1)", lambda c: c["ba"] >= 1),
    ("Benign (ii): >=2 Strong", lambda c: c["bs"] >= 2),
)
_LIKELY_BENIGN_RULES = (
    ("Likely benign (i): 1 Strong + 1 Supporting", lambda c: c["bs"] >= 1 and c["bp"] >= 1),
    ("Likely benign (ii): >=2 Supporting", lambda c: c["bp"] >= 2),
)


def combine_acmg(criteria: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    """Apply Table 5 of Richards et al. 2015 (Genet Med 17:405) to normalised criteria.

    Pathogenic and likely pathogenic rules are tried first on the pathogenic evidence alone, then the benign
    / likely benign rules on the benign evidence alone. If both sides reach a classification the evidence is
    contradictory and the result is Uncertain significance; if neither does, the result is also
    Uncertain significance.
    """

    counts = Counter({"vs": 0, "s": 0, "m": 0, "p": 0, "ba": 0, "bs": 0, "bp": 0})
    key = {
        ("pathogenic", "very_strong"): "vs",
        ("pathogenic", "strong"): "s",
        ("pathogenic", "moderate"): "m",
        ("pathogenic", "supporting"): "p",
        ("benign", "stand_alone"): "ba",
        ("benign", "strong"): "bs",
        ("benign", "supporting"): "bp",
    }
    for crit in criteria:
        counts[key[(crit["direction"], crit["strength"])]] += 1
    path_side = next(
        (("Pathogenic", name) for name, test in _PATHOGENIC_RULES if test(counts)), None
    ) or next(
        (("Likely pathogenic", name) for name, test in _LIKELY_PATHOGENIC_RULES if test(counts)),
        None,
    )
    benign_side = next(
        (("Benign", name) for name, test in _BENIGN_RULES if test(counts)), None
    ) or next(
        (("Likely benign", name) for name, test in _LIKELY_BENIGN_RULES if test(counts)), None
    )
    if path_side and benign_side:
        label, rule = (
            "Uncertain significance",
            f"contradictory evidence: {path_side[1]} vs {benign_side[1]}",
        )
    elif path_side:
        label, rule = path_side
    elif benign_side:
        label, rule = benign_side
    else:
        label, rule = (
            "Uncertain significance",
            "criteria do not meet any pathogenic or benign combination",
        )
    return {
        "classification": label,
        "rule": rule,
        "evidence_counts": {
            "pathogenic_very_strong": counts["vs"],
            "pathogenic_strong": counts["s"],
            "pathogenic_moderate": counts["m"],
            "pathogenic_supporting": counts["p"],
            "benign_stand_alone": counts["ba"],
            "benign_strong": counts["bs"],
            "benign_supporting": counts["bp"],
        },
    }


def normalize_criteria(items: object) -> list[dict[str, str]]:
    """Parse and de-duplicate a list of criteria (a code may only appear once)."""

    if isinstance(items, str):
        items = [c for c in re.split(r"[,;\s]+", items) if c]
    if isinstance(items, Mapping) and not {"code"} <= set(items):
        items = [
            {"code": k, "strength": v} if isinstance(v, str) else k
            for k, v in items.items()
            if v not in (False, None)
        ]
    if not isinstance(items, Sequence) or isinstance(items, bytes):
        fail(
            "criteria must be a list of ACMG codes.",
            "CLINICAL_BAD_CRITERION",
            "Example: ['PVS1', 'PM2_Supporting'].",
        )
    parsed: dict[str, dict[str, str]] = {}
    for item in items:
        crit = parse_criterion(item)
        if not crit:
            continue
        if crit["code"] in parsed and parsed[crit["code"]] != crit:
            fail(
                f"criterion {crit['code']} is given twice with different strengths.",
                "CLINICAL_BAD_CRITERION",
                "List each criterion once.",
            )
        parsed[crit["code"]] = crit
    return list(parsed.values())


def classify_acmg(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    if (
        isinstance(value, Mapping)
        and isinstance(value.get("variants"), Sequence)
        and not isinstance(value.get("variants"), (str, bytes))
    ):
        results = []
        for index, variant in enumerate(value["variants"]):
            if not isinstance(variant, Mapping) or "criteria" not in variant:
                fail(
                    f"variant {index} needs a 'criteria' list.",
                    "CLINICAL_BAD_INPUT",
                    "Use {'variants': [{'id': 'v1', 'criteria': ['PVS1', 'PM2']}]}.",
                )
            crit = normalize_criteria(variant["criteria"])
            results.append({"id": variant.get("id", index), **combine_acmg(crit), "criteria": crit})
        return {"framework": _FRAMEWORK, "results": results}
    items = value["criteria"] if isinstance(value, Mapping) and "criteria" in value else value
    if items is None:
        fail(
            "input is required.",
            "CLINICAL_BAD_INPUT",
            "Pass the ACMG criteria met, e.g. ['PVS1', 'PM2_Supporting'].",
        )
    crit = normalize_criteria(items)
    result = {"framework": _FRAMEWORK, **combine_acmg(crit), "criteria": crit}
    if isinstance(value, Mapping) and "id" in value:
        result["id"] = value["id"]
    if not crit:
        result["warnings"] = [
            "no criteria were supplied; Uncertain significance is the default, not an assessment"
        ]
    return result


_FRAMEWORK = (
    "ACMG/AMP 2015 (Richards et al., Genet Med 17:405-424) Table 5 combining rules applied to caller-supplied criteria; "
    "criterion strengths may be adjusted with a suffix (e.g. PM2_Supporting). The criteria themselves are NOT evaluated here"
)


# ---------------------------------------------------------------------------
# reading ClinVar-style VCFs
# ---------------------------------------------------------------------------
def read_clinical(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    rows = []
    for rec in vs.records:
        row = {
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "id": rec["id"],
            "ref": rec["ref"],
            "alt": list(rec["alts"]),
            "gene": gene_of(rec),
            "clinvar_significance_raw": field_value(
                rec, ("clnsig", "clinvar", "clinvar_significance", "clin_sig")
            ),
            "clinvar_terms": clinvar_terms(rec),
            "review_status": field_value(rec, ("clnrevstat",)),
            "disease": [d for d in str(field_value(rec, ("clndn",)) or "").split("|") if d],
            "consequence": consequence_terms(rec),
            "filter": rec["filter"],
        }
        if vs.samples:
            row["genotypes"] = {
                n: c.get("GT") for n, c in zip(vs.samples, rec["calls"], strict=False)
            }
        rows.append(row)
    return {
        "samples": list(vs.samples),
        "n_variants": len(rows),
        "variants": rows,
        "method": "VCF/records parsed into clinical fields: ClinVar CLNSIG (normalised to lower-case terms), GENEINFO/GENE gene symbol, CLNREVSTAT, CLNDN, MC consequence",
    }


# ---------------------------------------------------------------------------
# trio analysis
# ---------------------------------------------------------------------------
def _trio_names(vs: VariantSet, options: Mapping[str, object]) -> tuple[int, int, int, list[str]]:
    warnings: list[str] = []
    names = {k: options.get(k) for k in ("child", "father", "mother")}
    if not any(names.values()):
        trio = options.get("trio")
        if isinstance(trio, Mapping):
            names = {k: trio.get(k) for k in ("child", "father", "mother")}
    if not any(names.values()):
        if len(vs.samples) != 3:
            fail(
                "trio sample names are required (options.child/father/mother).",
                "CLINICAL_BAD_OPTION",
                "Name the three samples, or pass exactly three samples in the order child, father, mother.",
                config=True,
            )
        names = dict(zip(("child", "father", "mother"), vs.samples, strict=False))
        warnings.append("sample roles were not given; assuming column order child, father, mother")
    for key, name in names.items():
        if not isinstance(name, str) or name not in vs.samples:
            fail(
                f"{key} sample {name!r} not found in {vs.samples[:10]}.",
                "CLINICAL_BAD_OPTION",
                "Use names from the VCF header.",
                config=True,
            )
    if len(set(names.values())) != 3:
        fail(
            "child, father and mother must be three different samples.",
            "CLINICAL_BAD_OPTION",
            "Use three distinct sample names.",
            config=True,
        )
    return (
        vs.samples.index(names["child"]),
        vs.samples.index(names["father"]),  # type: ignore[arg-type]
        vs.samples.index(names["mother"]),
        warnings,
    )  # type: ignore[arg-type]


def _alt_count(alleles: Sequence[int | None]) -> int | None:
    if not is_called(alleles):
        return None
    return sum(1 for a in alleles if a)


def trio_classify(
    vs: VariantSet, options: Mapping[str, object]
) -> tuple[list[dict[str, Any]], list[str]]:
    """Per biallelic-split site inheritance category for the child (see :func:`analyze_trio`)."""

    ci, fi, mi, warnings = _trio_names(vs, options)
    min_dp = option_number(options, "min_depth", 0, minimum=0) or 0
    min_gq = option_number(options, "min_genotype_quality", 0, minimum=0) or 0
    split = split_multiallelic_set(vs)
    rows = []
    for rec in split.records:
        calls = [rec["calls"][i] for i in (ci, fi, mi)]
        gts = [call_gt(c)[0] for c in calls]
        child, father, mother = gts
        if len(child) != 2 or not is_called(child):
            continue
        quality_ok = True
        for call in calls:
            dp, gq = call_depth(call), number_or_none(call.get("GQ"))
            if (min_dp and (dp is None or dp < min_dp)) or (min_gq and (gq is None or gq < min_gq)):
                quality_ok = False
        c_alt, f_alt, m_alt = _alt_count(child), _alt_count(father), _alt_count(mother)
        f_set = set(father) if is_called(father) else None
        m_set = set(mother) if is_called(mother) else None
        category: str | None
        if f_set is not None and m_set is not None:
            consistent = (child[0] in f_set and child[1] in m_set) or (
                child[1] in f_set and child[0] in m_set
            )
            if not consistent:
                category = "de_novo" if (c_alt and f_alt == 0 and m_alt == 0) else "mendelian_error"
            elif c_alt == 0:
                category = None
            elif c_alt == 2:
                category = "homozygous_both_parents_carriers"
            else:
                category = (
                    "inherited_biparental"
                    if f_alt and m_alt
                    else "inherited_paternal"
                    if f_alt
                    else "inherited_maternal"
                )
        else:
            category = None if c_alt == 0 else "undetermined_missing_parent"
        if category is None:
            continue
        rows.append(
            {
                "site": site_name(rec),
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "ref": rec["ref"],
                "alt": rec["alts"][0],
                "gene": gene_of(rec),
                "category": category,
                "quality_pass": quality_ok,
                "child_gt": calls[0].get("GT"),
                "father_gt": calls[1].get("GT"),
                "mother_gt": calls[2].get("GT"),
                "_calls": calls,
                "_rec": rec,
            }
        )
    return rows, warnings


def analyze_trio(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    rows, warnings = trio_classify(vs, options)
    require_quality = option_bool(options, "require_quality_pass", True)
    used = [r for r in rows if r["quality_pass"] or not require_quality]
    counts = Counter(r["category"] for r in used)
    by_gene: dict[str, dict[str, list[dict[str, Any]]]] = {}
    ungenned = 0
    for r in used:
        if r["category"] in {"inherited_paternal", "inherited_maternal"}:
            if r["gene"] is None:
                ungenned += 1
                continue
            by_gene.setdefault(r["gene"], {"paternal": [], "maternal": []})[
                "paternal" if r["category"] == "inherited_paternal" else "maternal"
            ].append(r)
    compound = []
    for gene in sorted(by_gene):
        for p, m in product(by_gene[gene]["paternal"], by_gene[gene]["maternal"]):
            compound.append(
                {"gene": gene, "paternal_variant": p["site"], "maternal_variant": m["site"]}
            )
    if ungenned:
        warnings.append(
            f"{ungenned} singly-inherited variants have no gene annotation and were not considered for compound heterozygosity"
        )
    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in used]
    return {
        "method": "per-variant inheritance category from child/parent genotypes (after splitting multi-allelic sites): de_novo (child carries ALT, both "
        "parents homozygous without it), mendelian_error (other genotype combinations impossible under Mendelian transmission), "
        "homozygous_both_parents_carriers, inherited_paternal / _maternal / _biparental, undetermined_missing_parent; compound-heterozygous "
        "candidates = within a gene, one child-heterozygous variant carried only by the father and one only by the mother. "
        "Autosomal diploid logic; X-linked and mitochondrial inheritance are not modelled",
        "category_counts": dict(sorted(counts.items())),
        "variants": clean,
        "compound_heterozygous_candidates": compound,
        "n_filtered_low_quality": len(rows) - len(used),
        "warnings": warnings,
    }


def de_novo(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    min_child_dp = option_number(options, "min_child_depth", 10, minimum=0) or 0
    min_parent_dp = option_number(options, "min_parent_depth", 10, minimum=0) or 0
    min_gq = option_number(options, "min_genotype_quality", 20, minimum=0) or 0
    min_vaf = option_number(options, "min_child_vaf", 0.25, minimum=0, maximum=1) or 0
    max_parent_alt = option_int(options, "max_parent_alt_reads", 0, minimum=0)
    rows, warnings = trio_classify(
        vs, {k: v for k, v in options.items() if k not in {"min_depth", "min_genotype_quality"}}
    )
    kept, rejected = [], Counter()
    for r in rows:
        if r["category"] != "de_novo":
            continue
        child, father, mother = r["_calls"]
        reason = None
        trio = (
            ("child", child, min_child_dp),
            ("father", father, min_parent_dp),
            ("mother", mother, min_parent_dp),
        )
        # fixed check order, the first failing check wins: depth (child, father, mother), then GQ (child, father, mother)
        for name, call, floor in trio:
            dp = call_depth(call)
            if reason is None and floor and (dp is None or dp < floor):
                reason = f"{name}_low_depth"
        for name, call, _floor in trio:
            gq = number_or_none(call.get("GQ"))
            if reason is None and min_gq and (gq is None or gq < min_gq):
                reason = f"{name}_low_gq"
        c_ad = call_ad(child)
        if reason is None and c_ad is not None and len(c_ad) == 2:
            total = sum(c_ad)
            if total and c_ad[1] / total < min_vaf:
                reason = "child_low_vaf"
        if reason is None:
            for name, call in (("father", father), ("mother", mother)):
                ad = call_ad(call)
                if ad is not None and len(ad) == 2 and ad[1] > max_parent_alt:  # type: ignore[operator]
                    reason = f"{name}_alt_reads"
        if reason:
            rejected[reason] += 1
            continue
        kept.append({k: v for k, v in r.items() if not k.startswith("_")} | {"child_ad": c_ad})
    return {
        "method": "de novo candidates = child carries the ALT allele while both parents are called homozygous for the other allele(s), then "
        "filtered by depth, GQ, child ALT fraction (from AD when present) and parental ALT reads; each rejected candidate is counted once in "
        "rejected_by under the FIRST failing check in the fixed order: child/father/mother depth, child/father/mother GQ, child ALT fraction, "
        "father/mother ALT reads (the genotype checks come first and define the candidate set)",
        "parameters": {
            "min_child_depth": min_child_dp,
            "min_parent_depth": min_parent_dp,
            "min_genotype_quality": min_gq,
            "min_child_vaf": min_vaf,
            "max_parent_alt_reads": max_parent_alt,
        },
        "n_candidates_before_filters": sum(1 for r in rows if r["category"] == "de_novo"),
        "rejected_by": dict(sorted(rejected.items())),
        "n_de_novo": len(kept),
        "variants": kept,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# ctDNA / MRD: tumour-informed detection from plasma allele counts
# ---------------------------------------------------------------------------
def _tracked_sites(value: object, options: Mapping[str, object]) -> list[dict[str, Any]]:
    items = value["sites"] if isinstance(value, Mapping) and "sites" in value else None
    if (
        items is None
        and isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and value
        and all(isinstance(x, Mapping) and ("alt_count" in x) for x in value)
    ):
        items = value
    sites = []
    if items is not None:
        for index, item in enumerate(items):
            if not isinstance(item, Mapping):
                fail(
                    f"site {index} must be a mapping.",
                    "CLINICAL_BAD_INPUT",
                    "Use {'alt_count': 3, 'depth': 5000}.",
                )
            alt, depth = item.get("alt_count"), item.get("depth")
            if (
                any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in (alt, depth))
                or depth == 0
                or alt > depth
            ):  # type: ignore[operator]
                fail(
                    f"site {index} needs integer 0 <= alt_count <= depth, depth > 0.",
                    "CLINICAL_BAD_INPUT",
                    "Provide alt_count and depth (plasma reads).",
                )
            sites.append(
                {
                    "id": item.get("id", f"site_{index + 1}"),
                    "alt_count": alt,
                    "depth": depth,
                    "error_rate": item.get("error_rate"),
                }
            )
        return sites
    vs = split_multiallelic_set(load_variant_set(value))
    if not vs.samples:
        fail(
            "plasma allele counts are required.",
            "CLINICAL_BAD_INPUT",
            "Provide sites [{'alt_count','depth'}] or a VCF with AD for the plasma sample.",
        )
    name = options.get("sample")
    idx = 0 if name is None else (vs.samples.index(name) if name in vs.samples else -1)
    if idx < 0:
        fail(
            f"sample {name!r} not found.",
            "CLINICAL_BAD_OPTION",
            "Use a sample name from the VCF.",
            config=True,
        )
    for rec in vs.records:
        ad = call_ad(rec["calls"][idx])
        if ad is None or len(ad) != 2 or sum(ad) == 0:
            continue
        sites.append(
            {"id": site_name(rec), "alt_count": ad[1], "depth": sum(ad), "error_rate": None}
        )
    if not sites:
        fail(
            "no AD found for the plasma sample.",
            "CLINICAL_BAD_INPUT",
            "Provide FORMAT/AD or explicit alt_count/depth.",
        )
    return sites


def tumor_informed_detection(
    value: object, options: Mapping[str, object], *, alpha: float, min_sites: int, label: str
) -> dict[str, Any]:
    error = option_number(options, "background_error_rate", 0.001, minimum=1e-9, maximum=0.1)
    alpha = option_number(options, "alpha", alpha, minimum=1e-300, maximum=1) or alpha
    site_alpha = option_number(options, "site_alpha", 0.01, minimum=0, maximum=1)
    min_alt = option_int(options, "min_alt_reads", 2, minimum=1)
    min_sites = option_int(options, "min_positive_sites", min_sites, minimum=1) or min_sites
    factor = option_number(options, "vaf_to_fraction_factor", 2.0, minimum=1.0, maximum=10.0)
    assert (
        error is not None and site_alpha is not None and min_alt is not None and factor is not None
    )
    sites = _tracked_sites(value, options)
    rows = []
    for s in sites:
        e = s["error_rate"] if s["error_rate"] is not None else error
        p = binom_sf(s["alt_count"], s["depth"], e)
        rows.append(
            {
                **s,
                "error_rate": e,
                "vaf": s["alt_count"] / s["depth"],
                "p_value": p,
                "positive": s["alt_count"] >= min_alt and p <= site_alpha,
            }
        )
    alt = sum(s["alt_count"] for s in sites)
    depth = sum(s["depth"] for s in sites)
    pooled_p = binom_sf(alt, depth, error)
    positives = sum(1 for r in rows if r["positive"])
    vaf = alt / depth
    detected = pooled_p <= alpha and positives >= min_sites
    return {
        "method": f"{label}: tumour-informed detection. The plasma ALT reads and depths at the tracked tumour variants are pooled; the pooled "
        "ALT count is tested against the background error rate with a one-sided binomial test (P(X >= alt | depth, error)). "
        "Detected = pooled p <= alpha AND at least min_positive_sites sites individually positive (p <= site_alpha, ALT >= min_alt_reads). "
        "Estimated tumour fraction = vaf_to_fraction_factor * max(0, pooled VAF - background error) (factor 2 = clonal heterozygous "
        "copy-neutral SNVs). The background error rate must be measured for the assay; the default is a placeholder",
        "parameters": {
            "background_error_rate": error,
            "alpha": alpha,
            "site_alpha": site_alpha,
            "min_alt_reads": min_alt,
            "min_positive_sites": min_sites,
        },
        "detected": detected,
        "n_sites": len(sites),
        "n_positive_sites": positives,
        "pooled_alt_reads": alt,
        "pooled_depth": depth,
        "pooled_vaf": vaf,
        "pooled_p_value": pooled_p,
        "estimated_tumor_fraction": min(1.0, factor * max(0.0, vaf - error)),
        "sites": rows,
    }


# ---------------------------------------------------------------------------
# cancer evolution tree (perfect phylogeny from a binary mutation matrix)
# ---------------------------------------------------------------------------
def evolution_tree(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    matrix = value["mutations"] if isinstance(value, Mapping) and "mutations" in value else value
    samples = (
        list(value["samples"])
        if isinstance(value, Mapping) and isinstance(value.get("samples"), Sequence)
        else None
    )
    if not isinstance(matrix, Mapping) or not matrix or len(matrix) > 5000:
        fail(
            "input must map mutation names to presence vectors.",
            "CLINICAL_BAD_INPUT",
            "Use {'mutations': {'TP53': [1,1,0], 'KRAS': [1,0,0]}, 'samples': ['S1','S2','S3']}.",
        )
    width = None
    sets: dict[str, frozenset[int]] = {}
    for name, vector in matrix.items():
        if isinstance(vector, Mapping):
            if samples is None:
                samples = list(vector)
            vector = [vector.get(s, 0) for s in samples]
        if (
            not isinstance(vector, Sequence)
            or isinstance(vector, (str, bytes))
            or any(v not in (0, 1, True, False) for v in vector)
        ):
            fail(
                f"mutation {name!r} must be a 0/1 vector.",
                "CLINICAL_BAD_INPUT",
                "Use presence/absence (1/0) per sample or cell.",
            )
        if width is None:
            width = len(vector)
        elif len(vector) != width:
            fail(
                "all mutation vectors must have the same length.",
                "CLINICAL_BAD_INPUT",
                "Give one entry per sample/cell for every mutation.",
            )
        sets[str(name)] = frozenset(i for i, v in enumerate(vector) if v)
    assert width is not None
    if samples is None:
        samples = [f"sample_{i + 1}" for i in range(width)]
    if len(samples) != width:
        fail(
            "samples does not match the vector length.",
            "CLINICAL_BAD_INPUT",
            "Provide one sample name per vector entry.",
        )
    absent = [m for m, s in sets.items() if not s]
    present = {m: s for m, s in sets.items() if s}
    names = sorted(present)
    violations = []
    for a, b in ((x, y) for i, x in enumerate(names) for y in names[i + 1 :]):
        sa, sb = present[a], present[b]
        if sa & sb and not (sa <= sb or sb <= sa):
            violations.append(
                {
                    "mutation_1": a,
                    "mutation_2": b,
                    "reason": "present sets partially overlap (four-gamete violation)",
                }
            )
    if violations:
        return {
            "method": _TREE_NOTE,
            "perfect_phylogeny": False,
            "violations": violations,
            "tree": None,
            "n_mutations": len(present),
            "warnings": [f"{len(absent)} mutations absent from all samples were ignored"]
            if absent
            else [],
        }
    groups: dict[frozenset[int], list[str]] = {}
    for m in names:
        groups.setdefault(present[m], []).append(m)
    ordered = sorted(groups, key=lambda s: (-len(s), sorted(s)))
    nodes = [
        {
            "id": "root",
            "label": "germline",
            "mutations": [],
            "samples": list(samples),
            "parent": None,
            "children": [],
            "depth": 0,
        }
    ]
    node_of: dict[frozenset[int], dict[str, Any]] = {}
    for s in ordered:
        parents = [t for t in node_of if s < t]
        parent = node_of[min(parents, key=len)] if parents else nodes[0]
        node = {
            "id": f"n{len(nodes)}",
            "label": "+".join(groups[s]),
            "mutations": groups[s],
            "samples": [samples[i] for i in sorted(s)],
            "parent": parent["id"],
            "children": [],
            "depth": parent["depth"] + 1,
        }
        parent["children"].append(node["id"])
        nodes.append(node)
        node_of[s] = node
    by_id = {n["id"]: n for n in nodes}

    def newick(node: Mapping[str, Any]) -> str:
        label = re.sub(r"[\s(),:;\[\]']", "_", node["label"])
        kids = ",".join(newick(by_id[c]) for c in node["children"])
        return f"({kids}){label}" if kids else label

    sample_node = {}
    for i, name in enumerate(samples):
        holders = [s for s in node_of if i in s]
        sample_node[name] = node_of[min(holders, key=len)]["id"] if holders else "root"
    return {
        "method": _TREE_NOTE,
        "perfect_phylogeny": True,
        "violations": [],
        "n_mutations": len(present),
        "nodes": nodes,
        "newick": newick(nodes[0]) + ";",
        "sample_to_node": sample_node,
        "warnings": [f"{len(absent)} mutations absent from all samples were ignored"]
        if absent
        else [],
    }


_TREE_NOTE = (
    "perfect-phylogeny (infinite-sites) tree from a binary mutation x sample matrix: mutations with identical sample sets share a node, "
    "a node's parent is the node with the smallest strictly larger sample set, and the sample sets must be pairwise nested or disjoint "
    "(otherwise a four-gamete violation is reported and no tree is built)"
)


# ---------------------------------------------------------------------------
# tumour clones: binomial mixture of ALT fractions
# ---------------------------------------------------------------------------
def _logsumexp(values: Sequence[float]) -> float:
    top = max(values)
    return top + math.log(sum(math.exp(v - top) for v in values))


def _fit_binomial_mixture(
    data: Sequence[tuple[int, int]], k: int
) -> tuple[float, list[float], list[float], list[list[float]]]:
    n = len(data)
    vafs = sorted(a / d for a, d in data)
    centers = [min(max(vafs[min(n - 1, int((j + 0.5) * n / k))], 1e-3), 1 - 1e-3) for j in range(k)]
    if k > 1:
        for j in range(1, k):  # break ties deterministically
            if abs(centers[j] - centers[j - 1]) < 1e-3:
                centers[j] = min(0.999, centers[j - 1] + 0.01 * j)
    weights = [1.0 / k] * k
    prev = -math.inf
    ll = -math.inf
    resp: list[list[float]] = []
    for _ in range(1000):
        resp = []
        ll = 0.0
        for a, d in data:
            logs = [
                math.log(max(w, 1e-300)) + a * math.log(p) + (d - a) * math.log1p(-p)
                for w, p in zip(weights, centers, strict=False)
            ]
            norm = _logsumexp(logs)
            ll += norm
            resp.append([math.exp(x - norm) for x in logs])
        for j in range(k):
            r_sum = sum(r[j] for r in resp)
            weights[j] = r_sum / n
            denom = sum(r[j] * d for r, (_, d) in zip(resp, data, strict=False))
            if r_sum > 1e-12 and denom > 0:
                centers[j] = min(
                    max(sum(r[j] * a for r, (a, _) in zip(resp, data, strict=False)) / denom, 1e-6),
                    1 - 1e-6,
                )
        if abs(ll - prev) < 1e-10:
            break
        prev = ll
    return ll, weights, centers, resp


def infer_clones(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    purity = option_number(options, "purity", 1.0, minimum=0.01, maximum=1.0)
    max_clones = option_int(options, "max_clones", 4, minimum=1, maximum=8)
    default_cn = option_number(options, "copy_number", 2, minimum=0, maximum=20)
    assert purity is not None and max_clones is not None and default_cn is not None
    sites = _tracked_sites(value, options)
    n = len(sites)
    if n > 20_000:
        fail("more than 20000 mutations.", "CLINICAL_BAD_INPUT", "Subsample the mutations.")
    raw_items = value["sites"] if isinstance(value, Mapping) and "sites" in value else value
    meta = (
        list(raw_items)
        if isinstance(raw_items, Sequence)
        and not isinstance(raw_items, (str, bytes))
        and raw_items
        and all(isinstance(x, Mapping) for x in raw_items)
        and len(raw_items) == n
        else [{}] * n
    )
    data = [(s["alt_count"], s["depth"]) for s in sites]
    if any(d == 0 for _, d in data):
        fail("depth must be positive.", "CLINICAL_BAD_INPUT", "Remove sites without coverage.")
    best = None
    scores = []
    for k in range(1, min(max_clones, n) + 1):
        ll, weights, centers, resp = _fit_binomial_mixture(data, k)
        bic = -2 * ll + (2 * k - 1) * math.log(n)
        scores.append({"k": k, "log_likelihood": ll, "bic": bic})
        if best is None or bic < best[0] - 1e-9:
            best = (bic, k, weights, centers, resp)
    assert best is not None
    _, k, weights, centers, resp = best
    assign = [max(range(k), key=lambda j: r[j]) for r in resp]
    clusters = []
    for j in range(k):
        members = [i for i in range(n) if assign[i] == j]
        if not members:
            continue
        ccfs = []
        for i in members:
            cn = number_or_none(meta[i].get("copy_number")) if meta[i] else None
            cn = default_cn if cn is None else cn
            mult = number_or_none(meta[i].get("multiplicity")) if meta[i] else None
            mult = 1.0 if mult is None else mult
            vaf = data[i][0] / data[i][1]
            ccfs.append(vaf * (purity * cn + 2 * (1 - purity)) / (purity * mult))
        pooled = sum(data[i][0] for i in members) / sum(data[i][1] for i in members)
        clusters.append(
            {
                "cluster_vaf": centers[j],
                "pooled_vaf": pooled,
                "weight": weights[j],
                "n_mutations": len(members),
                "mean_ccf": sum(ccfs) / len(ccfs),
                "mutations": [sites[i]["id"] for i in members],
                "_j": j,
            }
        )
    clusters.sort(key=lambda c: -c["cluster_vaf"])
    for rank, c in enumerate(clusters, 1):
        c["cluster"] = rank
        c["role"] = "clonal (highest ALT fraction)" if rank == 1 else "subclonal"
    mapping = {c.pop("_j"): c["cluster"] for c in clusters}
    return {
        "method": "binomial mixture model of ALT/depth counts fitted by EM (deterministic quantile start), number of clones chosen by "
        "the smallest BIC; CCF = VAF*(purity*CN + 2*(1-purity))/(purity*multiplicity) per mutation (CN default 2, multiplicity 1; "
        "copy_number/multiplicity may be given per mutation). It does not model copy-number changes or phylogenetic ordering",
        "purity": purity,
        "n_mutations": n,
        "n_clones": len(clusters),
        "model_selection": scores,
        "clones": clusters,
        "assignments": [{"id": sites[i]["id"], "cluster": mapping[assign[i]]} for i in range(n)],
    }


# ---------------------------------------------------------------------------
# clinical QC of a call set
# ---------------------------------------------------------------------------
def clinical_qc(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    vs = load_variant_set(value)
    min_scr = option_number(options, "min_sample_call_rate", 0.95, minimum=0, maximum=1)
    min_site = option_number(options, "min_site_call_rate", 0.9, minimum=0, maximum=1)
    min_hwe = option_number(options, "min_hwe_pvalue", 1e-6, minimum=0, maximum=1)
    min_dp = option_number(options, "min_mean_depth", None, minimum=0)
    max_hh = option_number(options, "max_het_hom_alt_ratio", None, minimum=0)
    titv = options.get("ti_tv_range")
    if titv is not None and (
        not isinstance(titv, Sequence)
        or len(titv) != 2
        or not all(isinstance(x, (int, float)) for x in titv)
    ):
        fail(
            "options.ti_tv_range must be [low, high].",
            "CLINICAL_BAD_OPTION",
            "Example: [2.0, 3.5].",
            config=True,
        )
    assert min_scr is not None and min_site is not None and min_hwe is not None
    samples = sample_qc(vs) if vs.samples else []
    for row in samples:
        flags = []
        if row["call_rate"] is not None and row["call_rate"] < min_scr:
            flags.append("low_call_rate")
        if min_dp is not None and (row["mean_depth"] is None or row["mean_depth"] < min_dp):
            flags.append("low_mean_depth")
        if (
            max_hh is not None
            and row["het_hom_alt_ratio"] is not None
            and row["het_hom_alt_ratio"] > max_hh
        ):
            flags.append("high_het_hom_ratio")
        if (
            titv is not None
            and row["ti_tv"] is not None
            and not (titv[0] <= row["ti_tv"] <= titv[1])
        ):
            flags.append("ti_tv_out_of_range")
        row["flags"] = flags
    sites = site_qc(vs) if vs.samples else []
    flagged_sites = []
    for row in sites:
        flags = []
        if row["call_rate"] is not None and row["call_rate"] < min_site:
            flags.append("low_call_rate")
        if row["hwe_pvalue"] is not None and row["hwe_pvalue"] < min_hwe:
            flags.append("hwe_violation")
        if flags:
            flagged_sites.append({**row, "flags": flags})
    return {
        "method": "call-set QC: per-sample call rate, het/hom-alt ratio, Ti/Tv and depth; per-site call rate and exact HWE p-value; flags only "
        "where a threshold is set (call-rate and HWE thresholds have defaults, the others are off unless given)",
        "thresholds": {
            "min_sample_call_rate": min_scr,
            "min_site_call_rate": min_site,
            "min_hwe_pvalue": min_hwe,
            "min_mean_depth": min_dp,
            "max_het_hom_alt_ratio": max_hh,
            "ti_tv_range": titv,
        },
        "n_sites": len(vs.records),
        "n_samples": len(vs.samples),
        "samples": samples,
        "samples_failing": [r["sample"] for r in samples if r["flags"]],
        "n_sites_flagged": len(flagged_sites),
        "flagged_sites": flagged_sites[:1000],
        "passed": not any(r["flags"] for r in samples),
    }


# ---------------------------------------------------------------------------
# clinical report
# ---------------------------------------------------------------------------
def _report_rows(
    value: object, options: Mapping[str, object]
) -> tuple[list[dict[str, Any]], dict[str, Any], list[str]]:
    meta: dict[str, Any] = {}
    warnings: list[str] = []
    items: Any = value
    if (
        isinstance(value, Mapping)
        and "variants" in value
        and isinstance(value["variants"], Sequence)
        and not isinstance(value["variants"], (str, bytes))
    ):
        meta = {k: v for k, v in value.items() if k != "variants"}
        items = value["variants"]
    ref = reference_from_options(options)
    accession = options.get("hgvs_accession")
    rows = []
    if (
        isinstance(items, Sequence)
        and not isinstance(items, (str, bytes))
        and items
        and all(isinstance(x, Mapping) for x in items)
        and not all(("chrom" in x or "CHROM" in x) and ("ref" in x or "REF" in x) for x in items)
    ):
        source = [dict(x) for x in items]
    else:
        vs = load_variant_set(items)
        source = []
        for rec in vs.records:
            view = record_view(rec, vs.samples)
            view["gene"] = view.get("gene") or gene_of(rec)
            if view.get("clinvar") is None:
                cv = clinvar_terms(rec)
                if cv:
                    view["clinvar"] = "/".join(cv)
            if view.get("consequence") is None:
                terms = consequence_terms(rec)
                if terms:
                    view["consequence"] = ", ".join(terms)
            if isinstance(accession, (Mapping, str)):
                acc = accession.get(rec["chrom"]) if isinstance(accession, Mapping) else accession
                if acc and len(rec["alts"]) == 1 and not rec["alts"][0].startswith("<"):
                    view["hgvs_g"] = build_hgvs_genomic(
                        acc, rec["pos"], rec["ref"], rec["alts"][0], ref, rec["chrom"]
                    )
            if vs.samples and "genotypes" in view:
                gt = next(iter(view["genotypes"].values())).get("GT")
                if gt:
                    alleles, _ = call_gt({"GT": gt})
                    cls = diploid_class(alleles)
                    view.setdefault(
                        "zygosity",
                        {
                            "het": "heterozygous",
                            "hom_alt": "homozygous",
                            "hom_ref": "reference",
                        }.get(cls),
                    )
            source.append(view)
    for index, row in enumerate(source):
        label = row.get("classification")
        if label is None and row.get("criteria"):
            label = combine_acmg(normalize_criteria(row["criteria"]))["classification"]
            row["classification_basis"] = "ACMG/AMP 2015 combination of supplied criteria"
        if isinstance(label, str):
            key = label.strip()
            label = _ABBREV.get(key.upper()) or next(
                (c for c in CLASSES if c.lower() == key.lower().replace("_", " ")), None
            )
            if label is None:
                fail(
                    f"variant {index}: unknown classification {key!r}.",
                    "CLINICAL_BAD_INPUT",
                    f"Use one of {list(CLASSES)} or P/LP/VUS/LB/B.",
                )
        if label is None:
            label = "Not classified"
        row["classification"] = label
        rows.append(row)
    return rows, meta, warnings


def _describe(row: Mapping[str, Any]) -> str:
    for key in ("hgvs", "hgvs_c", "hgvs_g", "hgvs_p"):
        if row.get(key):
            return str(row[key])
    if row.get("chrom") is not None and row.get("pos") is not None:
        alt = row.get("alt")
        alt = ",".join(alt) if isinstance(alt, (list, tuple)) else alt
        return f"{row['chrom']}:{row['pos']} {row.get('ref', '')}>{alt}"
    return str(row.get("id", "variant"))


def clinical_report(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    fmt = option_choice(options, "format", "markdown", ("markdown", "text"))
    include_vus = option_bool(options, "include_vus", True)
    title = options.get("title", "Clinical variant report")
    if not isinstance(title, str) or not title.strip():
        fail(
            "options.title must be a non-empty string.",
            "CLINICAL_BAD_OPTION",
            "Pass a report title.",
            config=True,
        )
    rows, meta, warnings = _report_rows(value, options)
    if not rows:
        fail("no variants to report.", "CLINICAL_BAD_INPUT", "Provide at least one variant.")
    by_class = Counter(r["classification"] for r in rows)
    order = [*CLASSES, "Not classified"]
    rows.sort(key=lambda r: (order.index(r["classification"]), str(r.get("gene") or "")))
    reportable = [r for r in rows if r["classification"] in CLASSES[:2]]
    vus = [r for r in rows if r["classification"] in {CLASSES[2], "Not classified"}]
    benign = [r for r in rows if r["classification"] in CLASSES[3:]]

    def line(r: Mapping[str, Any]) -> list[str]:
        evidence = ", ".join(
            c["code"] if isinstance(c, Mapping) else str(c)
            for c in (normalize_criteria(r["criteria"]) if r.get("criteria") else [])
        )
        return [
            str(r.get("gene") or "-"),
            _describe(r),
            str(r.get("zygosity") or "-"),
            r["classification"],
            evidence or str(r.get("clinvar") or "-"),
        ]

    heads = ["Gene", "Variant", "Zygosity", "Classification", "Evidence"]

    def table(items: Sequence[Mapping[str, Any]]) -> list[str]:
        if fmt == "markdown":
            return [
                "| " + " | ".join(heads) + " |",
                "|" + "---|" * len(heads),
                *("| " + " | ".join(c.replace("|", "/") for c in line(r)) + " |" for r in items),
                "",
            ]
        return [
            "; ".join(f"{h}={c}" for h, c in zip(heads, line(r), strict=False)) for r in items
        ] + [""]

    h1, h2 = ("# ", "## ") if fmt == "markdown" else ("", "")
    out = [f"{h1}{title}", ""]
    for key in ("sample", "patient", "patient_id", "indication", "date"):
        if key in meta:
            out.append(f"{key.replace('_', ' ').capitalize()}: {meta[key]}")
    out += [
        f"Variants evaluated: {len(rows)} ({', '.join(f'{v} {k}' for k, v in sorted(by_class.items()))})",
        "",
    ]
    out += [f"{h2}Reportable findings (Pathogenic / Likely pathogenic)", ""] + (
        table(reportable) if reportable else ["None.", ""]
    )
    if include_vus:
        out += [f"{h2}Variants of uncertain significance / not classified", ""] + (
            table(vus) if vus else ["None.", ""]
        )
    out += [
        f"{h2}Not reported",
        "",
        f"{len(benign)} benign / likely benign variants and {0 if include_vus else len(vus)} uncertain variants are not listed.",
        "",
    ]
    out += [
        f"{h2}Method and limitations",
        "",
        "Classifications are taken from the input or computed from the supplied ACMG/AMP criteria with the 2015 combining rules "
        "(Richards et al., Genet Med 2015). The criteria were not re-evaluated against databases. This automatically generated "
        "summary is not a validated diagnostic report and must be reviewed by a qualified clinical geneticist.",
        "",
    ]
    text = "\n".join(out)
    path = write_output(options, text)
    clean = []
    for r in rows:
        clean.append({k: v for k, v in r.items() if k not in {"genotypes", "info"}})
    return {
        "format": fmt,
        "report": text,
        "n_variants": len(rows),
        "classification_counts": dict(sorted(by_class.items())),
        "n_reportable": len(reportable),
        "variants": clean,
        "output_path": path,
        "warnings": warnings,
    }
