"""Transparent rule-based ranking of candidate variants from annotations already present on the records.

The score is a documented additive heuristic; it is not a validated clinical classifier. Only the fields
that exist on a record contribute, and every component is returned so the ranking can be audited.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any

from ._impl_vcf import (
    fail,
    load_variant_set,
    number_or_none,
    option_bool,
    option_int,
    option_number,
)

# Sequence Ontology consequence -> severity weight (0..1)
CONSEQUENCE_WEIGHT: dict[str, float] = {
    "transcript_ablation": 1.0,
    "splice_acceptor_variant": 1.0,
    "splice_donor_variant": 1.0,
    "stop_gained": 1.0,
    "frameshift_variant": 1.0,
    "start_lost": 1.0,
    "stop_lost": 0.8,
    "transcript_amplification": 0.6,
    "inframe_insertion": 0.6,
    "inframe_deletion": 0.6,
    "protein_altering_variant": 0.6,
    "missense_variant": 0.5,
    "splice_region_variant": 0.3,
    "incomplete_terminal_codon_variant": 0.2,
    "start_retained_variant": 0.1,
    "stop_retained_variant": 0.1,
    "synonymous_variant": 0.05,
    "5_prime_utr_variant": 0.1,
    "3_prime_utr_variant": 0.1,
}
_CONSEQ_KEYS = (
    "consequence",
    "most_severe_consequence",
    "effect",
    "variant_consequence",
    "ann_effect",
    "csq_consequence",
)
_CLINVAR_KEYS = ("clinvar", "clinvar_significance", "clnsig", "clin_sig", "clinical_significance")
_AF_KEYS = (
    "gnomad_af",
    "gnomad_popmax_af",
    "af_popmax",
    "popmax_af",
    "population_af",
    "max_af",
    "gnomad_exomes_af",
    "gnomad_genomes_af",
)
_CADD_KEYS = ("cadd", "cadd_phred", "cadd_score")
_REVEL_KEYS = ("revel", "revel_score")
_GENE_KEYS = ("gene", "gene_symbol", "symbol", "geneinfo", "gene_name")


def field_value(rec: Mapping[str, Any], keys: Sequence[str]) -> Any:
    """First present annotation among ``keys`` (case-insensitive) in the record's extras, then INFO."""

    for source in (rec.get("extra") or {}, rec.get("info") or {}):
        lowered = {str(k).lower(): v for k, v in source.items()}
        for key in keys:
            if key in lowered and lowered[key] not in (None, "", "."):
                return lowered[key]
    return None


def consequence_terms(rec: Mapping[str, Any]) -> list[str]:
    raw = field_value(rec, _CONSEQ_KEYS)
    if raw is None:
        mc = field_value(rec, ("mc",))
        if isinstance(mc, str):
            return [part.split("|")[-1].lower() for part in mc.split(",") if part]
        return []
    if isinstance(raw, (list, tuple)):
        return [str(x).lower() for x in raw]
    return [t.strip().lower() for t in re.split(r"[&|,;]", str(raw)) if t.strip()]


def clinvar_terms(rec: Mapping[str, Any]) -> list[str]:
    raw = field_value(rec, _CLINVAR_KEYS)
    if raw is None:
        return []
    if isinstance(raw, (list, tuple)):
        parts = [str(x) for x in raw]
    else:
        parts = re.split(r"[/|;,&]", str(raw))
    return [re.sub(r"[\s-]+", "_", p.strip().lower()) for p in parts if p.strip()]


def gene_of(rec: Mapping[str, Any]) -> str | None:
    raw = field_value(rec, _GENE_KEYS)
    if raw is None:
        return None
    return str(raw).split(":")[0].split("|")[0].split(",")[0] or None


def score_record(
    rec: Mapping[str, Any], *, max_af: float, panel: set[str] | None, exclude_benign: bool
) -> dict[str, Any]:
    """Score one record. Components: consequence (0-1), ClinVar (-1..1), rarity (0-0.5), predictors (0-0.5), panel (+0.5)."""

    parts: dict[str, float] = {}
    excluded: str | None = None
    frequency_unknown = False
    terms = consequence_terms(rec)
    if terms:
        parts["consequence"] = max(CONSEQUENCE_WEIGHT.get(t, 0.0) for t in terms)
    cv = clinvar_terms(rec)
    if cv:
        pathogenic = [
            t for t in cv if t in {"pathogenic", "likely_pathogenic", "pathogenic_low_penetrance"}
        ]
        benign = [t for t in cv if t in {"benign", "likely_benign"}]
        if pathogenic and not benign:
            parts["clinvar"] = 1.0 if "pathogenic" in pathogenic else 0.8
        elif benign and not pathogenic:
            parts["clinvar"] = -1.0
            if exclude_benign:
                excluded = "ClinVar benign / likely benign"
        elif any(t.startswith("uncertain") for t in cv):
            parts["clinvar"] = 0.1
    af_raw = field_value(rec, _AF_KEYS)
    af = number_or_none(af_raw) if not isinstance(af_raw, (list, tuple)) else None
    if af is None and af_raw is None:
        parts["rarity"] = (
            0.5  # absent from population data is treated as rare (flagged via frequency_unknown)
        )
        frequency_unknown = True
    elif af is not None:
        if af > max_af:
            excluded = excluded or f"population AF {af:g} > {max_af:g}"
        parts["rarity"] = 0.5 if af <= 1e-4 else 0.3 if af <= 1e-3 else 0.1 if af <= 1e-2 else 0.0
    predictors = []
    cadd = number_or_none(field_value(rec, _CADD_KEYS))
    if cadd is not None:
        predictors.append(min(max(cadd, 0.0), 40.0) / 40.0)
    revel = number_or_none(field_value(rec, _REVEL_KEYS))
    if revel is not None:
        predictors.append(min(max(revel, 0.0), 1.0))
    if predictors:
        parts["predictors"] = 0.5 * max(predictors)
    gene = gene_of(rec)
    if panel is not None and gene is not None and gene.upper() in panel:
        parts["panel_gene"] = 0.5
    return {
        "score": sum(parts.values()),
        "components": parts,
        "excluded_reason": excluded,
        "gene": gene,
        "frequency_unknown": frequency_unknown,
    }


def prioritize(value: object, options: Mapping[str, object]) -> dict[str, Any]:
    """Rank variants; returns the ranked list (excluded variants listed separately)."""

    max_af = option_number(options, "max_af", 0.05, minimum=0, maximum=1)
    top_n = option_int(options, "top_n", None, minimum=1)
    exclude_benign = option_bool(options, "exclude_benign", True)
    genes = options.get("genes")
    if isinstance(genes, str):
        genes = [genes]
    if genes is not None and (
        not isinstance(genes, Sequence) or any(not isinstance(g, str) for g in genes)
    ):
        fail(
            "options.genes must be a list of gene symbols.",
            "VARIANTS_BAD_OPTION",
            "Example: ['BRCA1', 'TP53'].",
            config=True,
        )
    panel = {g.upper() for g in genes} if genes else None
    assert max_af is not None
    vs = load_variant_set(value)
    if not vs.records:
        fail(
            "no variants to prioritise.",
            "VARIANTS_BAD_INPUT",
            "Provide at least one variant record.",
        )
    ranked, excluded = [], []
    warnings: list[str] = []
    for rec in vs.records:
        scored = score_record(rec, max_af=max_af, panel=panel, exclude_benign=exclude_benign)
        row = {
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": list(rec["alts"]),
            "id": rec["id"],
            **scored,
        }
        (excluded if scored["excluded_reason"] else ranked).append(row)
    n_unknown = sum(1 for r in (*ranked, *excluded) if r["frequency_unknown"])
    if n_unknown:
        warnings.append(
            f"{n_unknown} of {len(vs.records)} variants have no population-frequency annotation and received the full rarity score "
            "(0.5) as if absent from population databases; they are flagged frequency_unknown. Annotate allele frequencies to make the ranking reliable"
        )
    ranked.sort(key=lambda r: (-r["score"], r["chrom"], r["pos"]))
    for index, row in enumerate(ranked, 1):
        row["rank"] = index
    n_ranked = len(ranked)
    if top_n:
        ranked = ranked[:top_n]
    return {
        "method": "additive heuristic score: consequence severity (Sequence Ontology table, 0-1) + ClinVar (pathogenic 1.0, likely pathogenic 0.8, "
        "uncertain 0.1, benign excluded) + rarity (no AF annotation, flagged frequency_unknown, or AF<=1e-4: 0.5, <=1e-3: 0.3, <=1e-2: 0.1) + 0.5*max(CADD/40, REVEL) "
        "+ 0.5 for genes in options.genes; variants with population AF > max_af or ClinVar benign are excluded. Not a validated classifier.",
        "parameters": {
            "max_af": max_af,
            "exclude_benign": exclude_benign,
            "genes": sorted(panel) if panel else None,
        },
        "ranked": ranked,
        "excluded": excluded,
        "n_ranked": n_ranked,
        "n_excluded": len(excluded),
        "n_frequency_unknown": n_unknown,
        "warnings": warnings,
    }
