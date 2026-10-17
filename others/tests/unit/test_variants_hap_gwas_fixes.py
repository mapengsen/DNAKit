"""Regression tests for the defects fixed after documenting variants, clinical, haplotypes, assembly, gwas,
experimental_design and forensic_reproductive."""

from __future__ import annotations

import math

import numpy as np
import pytest

import dnakit.assembly as asm
import dnakit.clinical as clinical
import dnakit.experimental_design as design
import dnakit.forensic_reproductive as forensic
import dnakit.gwas as gwas
import dnakit.haplotypes as hap
import dnakit.variants as variants
from dnakit.exceptions import InputFormatError
from dnakit.gwas import _impl_models
from dnakit.variants import _impl_analysis


def full(envelope: dict) -> dict:
    assert envelope["status"] == "ok", envelope
    return envelope


def res(envelope: dict) -> dict:
    return full(envelope)["result"]


def vcf(samples, rows):
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
    return "\n".join(["##fileformat=VCFv4.2", "\t".join(cols), *("\t".join(r.split()) for r in rows)]) + "\n"


# --------------------------------------------------------------------------- 1. kinship docstring
def test_kinship_docstring_matches_code_thresholds() -> None:
    doc = variants.calculate_kinship_king.__doc__
    for cut in ("0.3535", "0.1768", "0.0884", "0.0442"):
        assert cut in doc
    assert "0.354" not in doc and "0.177 " not in doc
    assert "IBD-based" not in doc
    assert [_impl_analysis._degree(k) for k in (0.36, 0.3535, 0.2, 0.1768, 0.1, 0.05, 0.0)] == [
        "duplicate_or_monozygotic_twin", "first_degree", "first_degree", "second_degree",
        "second_degree", "third_degree", "unrelated"]


# --------------------------------------------------------------------------- 2. de novo rejected_by order
def test_de_novo_rejected_by_is_first_failing_check_in_fixed_order() -> None:
    text = vcf(["child", "father", "mother"], [
        # child depth AND father depth fail: the first check (child depth) wins (used to be overwritten by father_low_depth)
        "chr1 100 . A G . PASS . GT:DP:GQ:AD 0/1:5:99:3,2 0/0:4:99:4,0 0/0:30:99:30,0",
        # child GQ and mother depth fail: depth is checked before GQ
        "chr1 200 . A G . PASS . GT:DP:GQ:AD 0/1:30:5:15,15 0/0:30:99:30,0 0/0:3:99:3,0",
        # child GQ and father GQ fail: child first
        "chr1 300 . A G . PASS . GT:DP:GQ:AD 0/1:30:5:15,15 0/0:30:5:30,0 0/0:30:99:30,0",
        # depth fine, GQ fine, child VAF and father ALT reads fail: VAF comes first
        "chr1 400 . A G . PASS . GT:DP:GQ:AD 0/1:30:99:28,2 0/0:30:99:27,3 0/0:30:99:30,0",
        # father then mother ALT reads
        "chr1 500 . A G . PASS . GT:DP:GQ:AD 0/1:30:99:15,15 0/0:30:99:30,0 0/0:30:99:29,1",
        "chr1 600 . A G . PASS . GT:DP:GQ:AD 0/1:30:99:15,15 0/0:30:99:30,0 0/0:30:99:30,0",
    ])
    out = res(clinical.detect_de_novo_variants(text, options={"child": "child", "father": "father", "mother": "mother"}))
    assert out["rejected_by"] == {"child_low_depth": 1, "mother_low_depth": 1, "child_low_gq": 1, "child_low_vaf": 1, "mother_alt_reads": 1}
    assert [v["pos"] for v in out["variants"]] == [600]
    assert sum(out["rejected_by"].values()) + out["n_de_novo"] == out["n_candidates_before_filters"] == 6
    assert "FIRST failing check" in clinical.detect_de_novo_variants.__doc__


# --------------------------------------------------------------------------- 3. frequency_unknown
PRIORITY = [
    {"chrom": "1", "pos": 10, "ref": "A", "alt": "G", "gene": "BRCA1", "consequence": "missense_variant", "gnomad_af": 0.00005},
    {"chrom": "1", "pos": 20, "ref": "A", "alt": "G", "gene": "X", "consequence": "stop_gained"},
    {"chrom": "1", "pos": 30, "ref": "A", "alt": "G", "gene": "Y", "consequence": "stop_gained"},
]


@pytest.mark.parametrize("fn", [variants.prioritize_candidate_variants, clinical.prioritize])
def test_prioritize_flags_missing_frequency(fn) -> None:
    env = full(fn(PRIORITY))
    result = env["result"]
    flags = {r["pos"]: r["frequency_unknown"] for r in result["ranked"]}
    assert flags == {10: False, 20: True, 30: True}
    assert all(r["components"]["rarity"] == 0.5 for r in result["ranked"])  # score unchanged
    assert result["n_frequency_unknown"] == 2
    assert len(env["warnings"]) == 1 and "2 of 3 variants" in env["warnings"][0]
    # all annotated: no warning
    annotated = [dict(r, gnomad_af=0.0001) for r in PRIORITY]
    clean = full(fn(annotated))
    assert clean["warnings"] == [] and clean["result"]["n_frequency_unknown"] == 0
    assert not any(r["frequency_unknown"] for r in clean["result"]["ranked"])
    # an excluded variant lacking AF (benign ClinVar) is still counted
    benign = full(fn([dict(PRIORITY[1], clinvar="Benign")]))
    assert benign["result"]["excluded"][0]["frequency_unknown"] is True


# --------------------------------------------------------------------------- 4. haplotypes
def test_infer_parentage_reports_not_assessed() -> None:
    data = {"child": ["0/0", "1/1", None], "candidates": {"A": ["0/1", "1/1", "0/0"], "B": [None, None, "0/1"], "C": ["1/1", "0/0", None]}}
    env = full(hap.infer_parentage(data))
    out = env["result"]
    assert out["candidates"]["B"]["compared_sites"] == 0
    assert out["not_assessed"] == ["B"]
    assert out["not_excluded"] == ["A"] and out["best_candidates"] == ["A"]
    assert out["candidates"]["B"]["excluded"] is False and "B" not in out["best_candidates"]
    assert any("B" in w and "no jointly called site" in w for w in env["warnings"])
    clean = full(hap.infer_parentage({"child": ["0/0"], "candidates": {"A": ["0/1"]}}))
    assert clean["result"]["not_assessed"] == [] and clean["warnings"] == []


def _brute_loglik(genotypes: dict[str, list[int]], freq: dict[str, float]) -> float:
    """Independent brute force: sum over ALL ordered haplotype pairs consistent with each genotype."""
    total = 0.0
    for geno in genotypes.values():
        p = 0.0
        for h1 in freq:
            for h2 in freq:
                if all(int(a) + int(b) == g for a, b, g in zip(h1, h2, geno)):
                    p += freq[h1] * freq[h2]
        total += math.log(p)
    return total


def test_reconstruct_log_likelihood_is_at_final_frequencies() -> None:
    genotypes = {"a": [0, 0], "c": [2, 2], "e": [1, 1]}
    one = res(hap.reconstruct(genotypes, options={"max_iterations": 1}))
    freq = {r["haplotype"]: r["frequency"] for r in one["frequencies"]}
    # hand calculation: one EM step from the uniform start gives f(00)=f(11)=2.5/6, f(01)=f(10)=0.5/6
    f00, f01 = 2.5 / 6, 0.5 / 6
    assert freq == pytest.approx({"00": f00, "11": f00, "01": f01, "10": f01})
    hand = math.log(f00 ** 2) + math.log(f00 ** 2) + math.log(2 * f00 * f00 + 2 * f01 * f01)
    assert one["log_likelihood"] == pytest.approx(hand, rel=1e-12)
    uniform_start = 2 * math.log(1 / 16) + math.log(4 / 16)  # the OLD value (E step before the update)
    assert abs(one["log_likelihood"] - uniform_start) > 0.1
    assert one["log_likelihood"] == pytest.approx(_brute_loglik(genotypes, freq), rel=1e-12)
    # converged run: likelihood at the reported frequencies, not below the one-step value (EM never decreases it)
    done = res(hap.reconstruct(genotypes, options={"tolerance": 1e-12}))
    dfreq = {r["haplotype"]: r["frequency"] for r in done["frequencies"]}
    full_freq = {h: dfreq.get(h, 0.0) for h in ("00", "01", "10", "11")}
    assert done["log_likelihood"] == pytest.approx(_brute_loglik(genotypes, full_freq), rel=1e-9)
    assert done["log_likelihood"] >= one["log_likelihood"] - 1e-12
    assert done["log_likelihood"] == pytest.approx(2 * math.log(0.25) + math.log(0.5), rel=1e-4)  # limit f00=f11=1/2
    assert "FINAL" in hap.reconstruct.__doc__


# --------------------------------------------------------------------------- 5. pangenome openness wording
def test_pangenome_openness_model_string_and_docstring() -> None:
    sets = {"g1": list("abcdefgh"), "g2": list("abcd") + list("ijkl"), "g3": list("ab") + list("mn"), "g4": ["a", "o"], "g5": ["a", "p"]}
    out = res(asm.calculate_pangenome_openness(sets, options={"n_permutations": 0}))
    assert "open if alpha < 1, closed if alpha >= 1" in out["model"] and "alpha > 1" not in out["model"]
    assert out["open"] is (out["alpha"] < 1.0)
    assert "open if alpha < 1 and closed if alpha >= 1" in asm.calculate_pangenome_openness.__doc__


# --------------------------------------------------------------------------- 6. gwas warnings
def _he_inputs(seed: int) -> dict:
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, 0.4, size=(12, 40))
    gv = g.sum(1) * 0.1
    y1 = gv + rng.normal(size=12) * 0.5
    y2 = gv + rng.normal(size=12) * 0.5
    return {"genotypes": g.tolist(), "phenotype_1": y1.tolist(), "phenotype_2": y2.tolist()}


def test_genetic_correlation_warns_when_rg_outside_unit_interval() -> None:
    env = full(gwas.estimate_genetic_correlation(_he_inputs(2)))
    assert abs(env["result"]["rg"]) > 1  # kept, not clamped
    assert any("outside [-1, 1]" in w for w in env["warnings"])
    assert _impl_models._rg_warnings(0.99) == [] and _impl_models._rg_warnings(None) == [] and _impl_models._rg_warnings(-1.2)


def test_heritability_warns_outside_zero_one() -> None:
    env = full(gwas.estimate_genetic_correlation(_he_inputs(10)))
    assert any("h2_trait1" in w and "exceeds 1" in w for w in env["warnings"])
    inp = _he_inputs(10)
    single = full(gwas.estimate_snp_heritability({"genotypes": inp["genotypes"], "phenotype": inp["phenotype_1"]}))
    assert single["result"]["h2"] > 1 and any("exceeds 1" in w for w in single["warnings"])
    assert _impl_models._h2_range_warnings(0.5, "h2") == [] and "negative" in _impl_models._h2_range_warnings(-0.1, "h2")[0]


def _ldsc_rows(h2: float, n_samples: int = 10_000, m: int = 200) -> list[dict]:
    ell = [1.0 + 0.05 * i for i in range(m)]
    return [{"snp": f"s{i}", "a1": "A", "a2": "G", "z": math.sqrt(max(1e-6, 1 + n_samples * h2 / m * ld)), "n": n_samples, "ldscore": ld} for i, ld in enumerate(ell)]


def test_summary_path_always_warns_no_standard_errors() -> None:
    inrange = full(gwas.estimate_snp_heritability(_ldsc_rows(0.4)))
    assert 0 < inrange["result"]["h2"] < 1
    assert len(inrange["warnings"]) == 1 and "no standard errors" in inrange["warnings"][0]
    high = full(gwas.estimate_snp_heritability(_ldsc_rows(1.5)))
    assert any("exceeds 1" in w for w in high["warnings"]) and any("no standard errors" in w for w in high["warnings"])
    negative = full(gwas.estimate_snp_heritability(_ldsc_rows(-0.001)))
    assert negative["result"]["h2"] < 0 and any("negative" in w for w in negative["warnings"])
    rows1, rows2 = _ldsc_rows(0.4), _ldsc_rows(0.3)
    rg = full(gwas.estimate_genetic_correlation({"trait1": rows1, "trait2": rows2}, options={"ld_scores": {r["snp"]: r["ldscore"] for r in rows1}}))
    assert any("no standard errors" in w for w in rg["warnings"])


# --------------------------------------------------------------------------- 7. experimental design envelope
def test_check_export_validate_accept_envelope_and_inner_result() -> None:
    env = full(design.design({"design_type": "rcbd", "treatments": ["A", "B", "C"], "replicates": 3, "seed": 1}))
    inner = env["result"]
    assert res(design.check(env)) == res(design.check(inner)) == res(design.check(inner["layout"]))
    assert res(design.check(env))["complete_blocks"] is True
    a, b, c = (res(design.export(x, options={"format": "csv"})) for x in (env, inner, inner["layout"]))
    assert a["text"] == b["text"] == c["text"] and a["n_rows"] == 9
    assert res(design.validate(env))["valid"] is True
    bad = {"status": "error", "result": None, "error": {"code": "X"}}
    for fn in (design.check, design.export, design.validate):
        with pytest.raises(InputFormatError) as info:
            fn(bad)
        assert info.value.code == "INVALID_DESIGN_ENVELOPE"
    assert "full envelope" in design.check.__doc__ and "full envelope" in design.export.__doc__


# --------------------------------------------------------------------------- 8. indistinguishable relationships
def test_infer_warns_about_indistinguishable_relationships() -> None:
    f = {"L": {"a": 0.2, "b": 0.3, "c": 0.5}}
    inp = {"profile_1": {"L": ["a", "b"]}, "profile_2": {"L": ["a", "c"]}, "frequencies": f}
    env = full(forensic.infer(inp, options={"relationships": ["unrelated", "half_sibling", "grandparent_grandchild", "uncle_nephew", "full_sibling"]}))
    rows = {r["relationship"]: r for r in env["result"]["relationships"]}
    assert rows["half_sibling"]["posterior"] == pytest.approx(rows["grandparent_grandchild"]["posterior"]) == pytest.approx(rows["uncle_nephew"]["posterior"])
    assert len(env["warnings"]) == 1
    w = env["warnings"][0]
    assert "half_sibling, grandparent_grandchild, uncle_nephew" in w and "indistinguishable" in w
    two = full(forensic.infer(inp, options={"relationships": ["half_sibling", "uncle_nephew"]}))
    assert len(two["warnings"]) == 1
    default = full(forensic.infer(inp))  # only half_sibling among the three by default
    assert default["warnings"] == []
    assert "indistinguishable" in forensic.infer.__doc__
