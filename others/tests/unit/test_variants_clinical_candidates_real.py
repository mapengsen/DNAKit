"""Real (non-placeholder) behaviour of the dnakit.variants and dnakit.clinical candidate functions.

Every expected value is a hand calculation (derivations in the comments) or an independent brute-force formula.
"""

from __future__ import annotations

import decimal
import math
from fractions import Fraction

import pytest

from dnakit import clinical, variants
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.variants import _impl_analysis, _impl_vcf


def run(function, value, **options):
    out = function(value, options=options or None)
    assert out["status"] == "ok", out
    return out["result"]


def vcf(samples, rows, header=("##fileformat=VCFv4.2",)):
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] + (
        ["FORMAT", *samples] if samples else []
    )
    return "\n".join([*header, "\t".join(cols), *("\t".join(r.split()) for r in rows)]) + "\n"


# 4 sites, 3 samples
VCF_A = vcf(
    ["S1", "S2", "S3"],
    [
        "chr1 100 . A G 50 PASS DP=30 GT:DP:GQ 0/0:10:50 0/1:10:40 1/1:10:30",
        "chr1 200 . C T 20 q10 . GT 0/1 0/1 ./.",
        "chr1 300 . AT A 60 PASS . GT 0/0 0/1 0/0",
        "chr1 400 . G A,C 30 PASS . GT 1/2 0/1 0/0",
    ],
)

REF8 = "GCACACAT"  # G1 C2 A3 C4 A5 C6 A7 T8


# ---------------------------------------------------------------------------
# envelope
# ---------------------------------------------------------------------------
def test_envelope_names_and_status():
    out = variants.summarize(VCF_A)
    assert (out["module"], out["operation"], out["status"]) == ("variants", "summarize", "ok")
    out = clinical.classify(["PVS1", "PS1"])
    assert (out["module"], out["operation"], out["status"]) == ("clinical", "classify", "ok")


# ---------------------------------------------------------------------------
# classification / summary / Ti-Tv
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "ref,alt,kind,sub",
    [
        ("A", "G", "SNV", None),
        ("AT", "GT", "SNV", None),  # shared flank trimmed
        ("AC", "GT", "MNV", None),
        ("A", "AT", "indel", "insertion"),
        ("AT", "A", "indel", "deletion"),
        ("AT", "GC", "MNV", None),
        ("ATT", "G", "indel", "complex"),
        ("A", "A" + "T" * 60, "SV", "insertion"),
        ("A", "<DEL>", "SV", "DEL"),
        ("A", "A]chr2:5]", "SV", "BND"),
        ("A", "*", "other", "spanning_deletion"),
    ],
)
def test_allele_classification(ref, alt, kind, sub):
    info = _impl_vcf.classify_alleles(ref, alt)
    assert (info["type"], info["subtype"]) == (kind, sub)


def test_summarize_counts_and_titv():
    res = run(variants.summarize, VCF_A)
    # alleles: A>G, C>T, AT>A, G>A, G>C -> 4 SNV + 1 indel; Ti = A>G, C>T, G>A (3), Tv = G>C (1)
    assert res["type_counts"] == {"SNV": 4, "indel": 1}
    assert (res["transitions"], res["transversions"], res["ti_tv"]) == (3, 1, 3.0)
    assert res["n_records"] == 4 and res["n_alt_alleles"] == 5 and res["n_multiallelic_sites"] == 1
    assert res["filter_counts"] == {"PASS": 3, "q10": 1}
    s1, _s2, s3 = res["sample_qc"]
    assert (s1["n_het"], s1["n_hom_ref"], s1["call_rate"]) == (2, 2, 1.0)  # S1 = 0/0, 0/1, 0/0, 1/2
    assert (s3["n_missing"], s3["call_rate"], s3["n_hom_alt"]) == (1, 0.75, 1)


def test_summarize_mt_matches_summarize():
    assert (
        run(variants.summarize_mt, VCF_A)["type_counts"]
        == run(variants.summarize, VCF_A)["type_counts"]
    )


# ---------------------------------------------------------------------------
# HWE exact test
# ---------------------------------------------------------------------------
def _brute_hwe(aa, ab, bb):
    n = aa + ab + bb
    rare = 2 * min(aa, bb) + ab
    total = 2 * n
    probs = {}
    for het in range(rare % 2, rare + 1, 2):
        hom_r = (rare - het) // 2
        hom_c = n - het - hom_r
        if hom_c < 0:
            continue
        probs[het] = (
            Fraction(
                math.factorial(n),
                math.factorial(hom_r) * math.factorial(het) * math.factorial(hom_c),
            )
            * 2**het
            * Fraction(math.factorial(rare) * math.factorial(total - rare), math.factorial(total))
        )
    obs = probs[ab]
    return float(sum(p for p in probs.values() if p <= obs))


@pytest.mark.parametrize(
    "counts,expected",
    [
        ([1, 0, 1], 1 / 3),  # n=2: P(het=0)=1/3, P(het=2)=2/3 -> p = 1/3
        ([0, 2, 0], 1.0),  # observed het=2 is the most likely table -> p = 1/3+2/3
        ([0, 3, 0], 0.4),  # n=3, 3 copies each: P(het=1)=0.6, P(het=3)=0.4
        ([1, 1, 1], 1.0),  # P(het=1)=0.6 is the larger one
    ],
)
def test_hwe_hand_values(counts, expected):
    assert run(variants.calculate_hwe_pvalue, counts)["p_value"] == pytest.approx(expected)


@pytest.mark.parametrize(
    "counts", [[10, 20, 30], [50, 3, 47], [0, 40, 0], [5, 0, 95], [120, 160, 60]]
)
def test_hwe_matches_brute_force(counts):
    assert run(variants.calculate_hwe_pvalue, counts)["p_value"] == pytest.approx(
        _brute_hwe(*counts), rel=1e-9
    )


def test_hwe_from_vcf_and_mid_p():
    res = run(variants.calculate_hwe_pvalue, VCF_A)
    assert [round(s["p_value"], 6) for s in res["sites"][:1]] == [1.0]  # site 1: counts 1/1/1
    mid = run(variants.calculate_hwe_pvalue, [0, 3, 0], mid_p=True)["p_value"]
    assert mid == pytest.approx(0.4 - 0.5 * 0.4)


# ---------------------------------------------------------------------------
# normalisation / left alignment / HGVS
# ---------------------------------------------------------------------------
def recs(*items):
    return [{"chrom": "chr1", "pos": p, "ref": r, "alt": a} for p, r, a in items]


@pytest.mark.parametrize(
    "pos,ref,alt,exp",
    [
        (
            5,
            "ACA",
            "A",
            (1, "GCA", "G"),
        ),  # deleting CA (6-7) == deleting CA (2-3); leftmost anchor is G1
        (7, "A", "ACA", (1, "G", "GCA")),  # inserting CA after A7 == inserting CA after G1
        (7, "AT", "GT", (7, "A", "G")),  # MNV flank trimmed to a SNV
        (3, "A", "G", (3, "A", "G")),
    ],
)
def test_left_align_with_reference(pos, ref, alt, exp):
    res = run(variants.left_align, recs((pos, ref, alt)), reference=REF8)
    r = res["records"][0]
    assert (r["pos"], r["ref"], r["alt"]) == (exp[0], exp[1], [exp[2]])
    assert res["left_aligned"] is True


def test_left_align_at_sequence_start_pads_right():
    res = run(variants.left_align_variants, recs((1, "AA", "A")), reference="AAAT")
    r = res["records"][0]
    assert (r["pos"], r["ref"], r["alt"]) == (1, "AA", ["A"])


def test_normalize_without_reference_only_trims():
    res = run(variants.normalize, recs((100, "CAT", "CGT"), (5, "ACA", "A")))
    got = {(r["pos"], r["ref"], r["alt"][0]) for r in res["records"]}
    assert got == {(101, "A", "G"), (5, "ACA", "A")}
    assert res["left_aligned"] is False


def test_normalize_reference_mismatch_warns_and_keeps_record():
    out = variants.normalize_vcf(recs((5, "TTT", "T")), options={"reference": REF8})
    assert out["status"] == "ok" and out["warnings"]
    assert out["result"]["records"][0]["pos"] == 5


def test_normalize_multiallelic_joint_and_split_option():
    res = run(variants.normalize_mt, recs((5, "ACA", "A,ACACA")), reference=REF8)
    r = res["records"][0]
    # CA deletion and CA insertion both left-align to anchor G1
    assert (r["pos"], r["ref"], r["alt"]) == (1, "GCA", ["G", "GCACA"])
    split = run(
        variants.normalize, recs((5, "ACA", "A,ACACA")), reference=REF8, split_multiallelic=True
    )
    assert [r["alt"] for r in split["records"]] == [["G"], ["GCACA"]]


@pytest.mark.parametrize(
    "pos,ref,alt,hgvs",
    [
        (123, "A", "G", "NC_1.1:g.123A>G"),
        (10, "AC", "GT", "NC_1.1:g.10_11delinsGT"),
        (5, "AC", "A", "NC_1.1:g.6del"),
        (5, "A", "AGG", "NC_1.1:g.5_6insGG"),
    ],
)
def test_hgvs_without_reference(pos, ref, alt, hgvs):
    assert _impl_vcf.build_hgvs_genomic("NC_1.1", pos, ref, alt) == hgvs


def test_hgvs_with_reference_three_prime_shift_and_dup():
    ref = _impl_vcf.load_reference(REF8)
    assert (
        _impl_vcf.build_hgvs_genomic("NC_1.1", 1, "GCA", "G", ref, "chr1") == "NC_1.1:g.6_7del"
    )  # CA deletion 3' shifted
    assert (
        _impl_vcf.build_hgvs_genomic("NC_1.1", 1, "G", "GCA", ref, "chr1") == "NC_1.1:g.6_7dup"
    )  # insertion repeating CA
    assert _impl_vcf.build_hgvs_genomic("NC_1.1", 1, "G", "GT", ref, "chr1") == "NC_1.1:g.1_2insT"
    res = run(variants.normalize, recs((5, "ACA", "A")), reference=REF8, hgvs_accession="NC_1.1")
    assert res["records"][0]["hgvs_g"] == ["NC_1.1:g.6_7del"]


# ---------------------------------------------------------------------------
# split / merge / filter
# ---------------------------------------------------------------------------
def test_split_multiallelic_remaps_gt_and_ad():
    text = vcf(["S1", "S2"], ["chr1 400 . G A,C 30 PASS AC=2,1 GT:AD 1/2:5,3,2 0/1:6,4,0"])
    for fn in (variants.split_multiallelic_variants, variants.split_multiallelic_mt):
        res = run(fn, text)
        a, c = res["records"]
        assert (a["alt"], c["alt"]) == (["A"], ["C"])
        assert a["genotypes"]["S1"] == {"GT": "1/0", "AD": "5,3"} and c["genotypes"]["S1"] == {
            "GT": "0/1",
            "AD": "5,2",
        }
        assert a["genotypes"]["S2"]["GT"] == "0/1" and c["genotypes"]["S2"]["GT"] == "0/0"
        assert (a["info"]["AC"], c["info"]["AC"]) == ("2", "1")
        assert res["split_records"] == 1


MERGE_A = vcf(["S1"], ["chr1 100 . A G 50 PASS . GT 0/1", "chr1 300 . T C 40 PASS . GT 1/1"])
MERGE_B = vcf(["S2"], ["chr1 100 . A G 30 PASS . GT 1/1", "chr1 200 . C T 10 . . GT 0/1"])


@pytest.mark.parametrize("fn", [variants.merge, variants.merge_vcf_files, variants.merge_mt])
def test_merge_samples(fn):
    res = run(fn, [MERGE_A, MERGE_B])
    assert res["samples"] == ["S1", "S2"]
    rows = {r["pos"]: r for r in res["records"]}
    assert sorted(rows) == [100, 200, 300]
    assert (
        rows[100]["qual"] == 50
        and rows[100]["genotypes"]["S1"]["GT"] == "0/1"
        and rows[100]["genotypes"]["S2"]["GT"] == "1/1"
    )
    assert (
        rows[200]["genotypes"]["S1"]["GT"] == "./." and rows[300]["genotypes"]["S2"]["GT"] == "./."
    )


def test_merge_records_mode_dedupes_and_sorts():
    other = vcf(["S1"], ["chr1 100 . A G 50 PASS . GT 0/1", "chr1 50 . C A 5 PASS . GT 1/1"])
    res = run(variants.merge, [MERGE_A, other], mode="records")
    assert [r["pos"] for r in res["records"]] == [50, 100, 300]


def test_filter_criteria():
    res = run(variants.filter, VCF_A, pass_only=True)
    assert [r["pos"] for r in res["records"]] == [100, 300, 400] and res["reason_counts"] == {
        "not_pass": 1
    }
    assert [r["pos"] for r in run(variants.filter, VCF_A, min_qual=30)["records"]] == [
        100,
        300,
        400,
    ]  # QUAL 20 fails, 30 passes
    assert [r["pos"] for r in run(variants.filter, VCF_A, max_missing=0.2)["records"]] == [
        100,
        300,
        400,
    ]  # site 200: 1/3 missing
    assert [r["pos"] for r in run(variants.filter, VCF_A, variant_types=["indel"])["records"]] == [
        300
    ]
    assert [r["pos"] for r in run(variants.filter, VCF_A, regions=["chr1:150-350"])["records"]] == [
        200,
        300,
    ]
    # ALT frequencies: 3/6, 2/4, 1/6, 3/6
    assert [r["pos"] for r in run(variants.filter, VCF_A, max_af=0.2)["records"]] == [300]
    assert [r["pos"] for r in run(variants.filter_mt, VCF_A, min_af=0.4)["records"]] == [
        100,
        200,
        400,
    ]


def test_filter_flag_mode_and_genotype_masking():
    flagged = run(variants.filter, VCF_A, mode="flag", min_qual=40)
    assert {r["pos"]: r["filter"] for r in flagged["records"]} == {
        100: "PASS",
        200: "q10;low_qual",
        300: "PASS",
        400: "low_qual",
    }
    text = vcf(["S1", "S2"], ["chr1 1 . A G 50 PASS . GT:GQ 0/1:50 1/1:10"])
    res = run(variants.filter, text, min_genotype_quality=20)
    assert res["records"][0]["genotypes"]["S2"]["GT"] == "./." and res["genotypes_masked"] == 1


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def test_compare_vcfs_metrics_and_concordance():
    left = vcf(
        ["S1"],
        [
            "chr1 100 . A G . . . GT 0/1",
            "chr1 200 . C T . . . GT 0/1",
            "chr1 300 . AT A . . . GT 0/1",
        ],
    )
    right = vcf(
        ["S1"],
        [
            "chr1 100 . A G . . . GT 0/1",
            "chr1 300 . AT A . . . GT 1/1",
            "chr1 400 . G A . . . GT 0/1",
        ],
    )
    for fn in (variants.compare_vcfs, variants.compare_mt):
        res = run(fn, {"left": left, "right": right})
        # shared = {100, 300}; precision 2/3, recall 2/3, jaccard 2/4
        assert (res["shared_count"], res["only_left_count"], res["only_right_count"]) == (2, 1, 1)
        assert res["precision_right_vs_left"] == pytest.approx(2 / 3) and res[
            "recall_right_vs_left"
        ] == pytest.approx(2 / 3)
        assert res["f1"] == pytest.approx(2 / 3) and res["jaccard"] == pytest.approx(0.5)
        gc = res["genotype_concordance"]["S1"]
        assert (gc["concordant"], gc["discordant"], gc["genotype_concordance"]) == (1, 1, 0.5)
    assert (
        run(variants.compare_vcfs, {"left": left, "right": right}, match="position")["shared_count"]
        == 2
    )


def test_compare_vcfs_normalises_with_reference():
    left = recs((1, "GCA", "G"))
    right = recs((5, "ACA", "A"))
    assert run(variants.compare_vcfs, {"left": left, "right": right})["shared_count"] == 0
    assert (
        run(variants.compare_vcfs, {"left": left, "right": right}, reference=REF8)["shared_count"]
        == 1
    )


# ---------------------------------------------------------------------------
# genotype likelihood
# ---------------------------------------------------------------------------
def test_genotype_likelihood_hand_values():
    res = run(variants.calculate_genotype_likelihood, {"ref_count": 10, "alt_count": 0})["result"]
    # log10 L: 10*log10(0.99) = -0.04365 ; 10*log10(0.5) = -3.0103 ; 10*log10(0.01) = -20
    assert res["log10_likelihoods"] == pytest.approx([-0.043648, -3.010300, -20.0], abs=1e-5)
    assert (
        res["phred_scaled"] == [0, 30, 200]
        and res["best_genotype"] == "0/0"
        and res["genotype_quality"] == 30
    )
    res = run(variants.calculate_genotype_likelihood, {"ref_count": 5, "alt_count": 5})["result"]
    # 0/0 = 1/1 = 5 log10(.99) + 5 log10(.01) = -10.0218 ; 0/1 = -3.0103 -> PL 70
    assert (
        res["phred_scaled"] == [70, 0, 70]
        and res["best_genotype"] == "0/1"
        and res["genotype_quality"] == 70
    )


def test_genotype_likelihood_from_ad_and_base_qualities():
    text = vcf(["S1"], ["chr1 1 . A G 50 PASS . GT:AD 0/1:5,5"])
    assert run(variants.calculate_genotype_likelihood, text)["results"][0]["phred_scaled"] == [
        70,
        0,
        70,
    ]
    reads = {
        "reads": [{"allele": "alt", "base_quality": 20}, {"allele": "ref", "base_quality": 20}]
    }
    # e=0.01 each: 0/0 -> log10(.01)+log10(.99) ; 0/1 -> 2*log10(.5) ; 1/1 same as 0/0
    out = run(variants.calculate_genotype_likelihood, reads)["result"]
    assert out["best_genotype"] == "0/1"
    assert out["log10_likelihoods"][1] == pytest.approx(2 * math.log10(0.5))


# ---------------------------------------------------------------------------
# LD, ROH, IBD, haplotypes
# ---------------------------------------------------------------------------
def _dosage_vcf(rows):
    n = len(rows[0])
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}
    lines = [
        f"chr1 {100 + 10 * i} . A G . . . GT " + " ".join(gt[x] for x in row)
        for i, row in enumerate(rows)
    ]
    return vcf([f"S{k}" for k in range(n)], lines)


def test_ld_unphased_dosage_correlation():
    text = _dosage_vcf([[0, 1, 2, 1], [0, 1, 2, 1], [2, 1, 0, 1], [0, 0, 1, 1]])
    res = run(variants.calculate_ld_matrix, text)
    m = res["matrix"]
    assert m[0][1] == pytest.approx(1.0) and m[0][2] == pytest.approx(1.0)
    # v1 mean 1, v4 mean .5: cov sum = 1, sxx = 2, syy = 1 -> r = 1/sqrt(2) -> r^2 = 0.5
    assert m[0][3] == pytest.approx(0.5)
    assert run(variants.calculate_ld_matrix, text, statistic="r")["matrix"][0][2] == pytest.approx(
        -1.0
    )


def test_ld_phased_d_prime():
    text = vcf(
        ["S1", "S2"],
        [
            "chr1 1 . A G . . . GT 0|0 1|1",
            "chr1 2 . A G . . . GT 0|0 1|1",
            "chr1 3 . A G . . . GT 0|1 1|0",
        ],
    )
    # haplotypes (v1,v2): 00,00,11,11 -> D = 0.5 - 0.25 = 0.25, r2 = 1, D' = 1
    assert run(variants.calculate_ld_matrix, text)["matrix"][0][1] == pytest.approx(1.0)
    assert run(variants.calculate_ld_matrix, text, statistic="d_prime")["matrix"][0][
        1
    ] == pytest.approx(1.0)
    # v1 vs v3 haplotypes: (0,0)? s1: v1 0|0, v3 0|1 -> (0,0),(0,1); s2: 1|1, 1|0 -> (1,1),(1,0) -> D = 0
    assert run(variants.calculate_ld_matrix, text)["matrix"][0][2] == pytest.approx(0.0)


def test_roh_runs():
    gts = ["0/0", "1/1", "0/0", "1/1", "0/0", "1/1", "0/1", "0/0", "1/1", "0/0"]
    text = vcf(["S1"], [f"chr1 {1000 * (i + 1)} . A G . . . GT {g}" for i, g in enumerate(gts)])
    strict = run(
        variants.calculate_runs_of_homozygosity,
        text,
        min_snps=3,
        min_length_bp=1000,
        max_heterozygotes=0,
    )
    assert [(r["start"], r["end"], r["n_snps"]) for r in strict["runs"]] == [
        (1000, 6000, 6),
        (8000, 10000, 3),
    ]
    loose = run(
        variants.calculate_runs_of_homozygosity,
        text,
        min_snps=3,
        min_length_bp=1000,
        max_heterozygotes=1,
    )
    r = loose["runs"][0]
    assert (r["start"], r["end"], r["length_bp"], r["n_snps"], r["n_heterozygous"]) == (
        1000,
        10000,
        9001,
        10,
        1,
    )
    assert len(loose["runs"]) == 1 and loose["total_roh_bp"] == 9001


def test_kinship_king_robust():
    rows = [
        "0/1 0/1 0/1",
        "0/1 0/1 1/1",
        "0/0 0/0 1/1",
        "1/1 1/1 0/0",
        "0/1 0/1 0/0",
        "0/0 0/0 0/1",
    ]
    text = vcf(["S1", "S2", "S3"], [f"chr1 {i + 1} . A G . . . GT {r}" for i, r in enumerate(rows)])
    pairs = {
        (p["sample_1"], p["sample_2"]): p
        for p in run(variants.calculate_kinship_king, text)["pairs"]
    }
    # S1,S2: het-het 3, opposite 0, hets 3+3 -> 0.5
    p12 = pairs[("S1", "S2")]
    assert (
        p12["kinship"] == pytest.approx(0.5)
        and p12["ibs0_fraction"] == 0
        and p12["relationship"] == "duplicate_or_monozygotic_twin"
    )
    # S1,S3: het-het 1, opposite 2, hets 3 (S1) + 2 (S3) -> (1-4)/5 = -0.6
    p13 = pairs[("S1", "S3")]
    assert (
        p13["kinship"] == pytest.approx(-0.6)
        and p13["ibs0_fraction"] == pytest.approx(2 / 6)
        and p13["relationship"] == "unrelated"
    )


def test_haplotype_frequency_phased_counting():
    res = run(variants.calculate_haplotype_frequency, {"haplotypes": ["01", "01", "10", "00"]})
    freq = {h["haplotype"]: h["frequency"] for h in res["haplotypes"]}
    assert freq == {"0-1": 0.5, "1-0": 0.25, "0-0": 0.25}


def test_haplotype_frequency_em_fixed_point():
    # 10 x (00/00), 10 x (11/11), 10 double heterozygotes. With y=f(00)=f(11), x=f(01)=f(10) the EM fixed point obeys
    # q = y^2/(y^2+x^2), y = (20+10q)/60, x = 10(1-q)/60 ; q = 1 is the stable solution -> y = 0.5, x = 0.
    people = [[(0, 0), (0, 0)]] * 10 + [[(1, 1), (1, 1)]] * 10 + [[(0, 1), (0, 1)]] * 10
    genos = [["/".join(map(str, g)) for g in person] for person in people]
    res = run(variants.calculate_haplotype_frequency, {"genotypes": genos})
    freq = {h["haplotype"]: h["frequency"] for h in res["haplotypes"]}
    assert freq["0-0"] == pytest.approx(0.5, abs=1e-6) and freq["1-1"] == pytest.approx(
        0.5, abs=1e-6
    )
    assert freq["0-1"] == pytest.approx(0.0, abs=1e-6) and sum(freq.values()) == pytest.approx(1.0)
    # a single double heterozygote carries no phase information: all four haplotypes stay at 1/4
    one = run(variants.calculate_haplotype_frequency, {"genotypes": [["0/1", "0/1"]]})
    assert all(h["frequency"] == pytest.approx(0.25) for h in one["haplotypes"])


# ---------------------------------------------------------------------------
# copy number
# ---------------------------------------------------------------------------
def test_copy_number_profile():
    bins = [
        {"chrom": "chr1", "start": 1 + 100 * i, "end": 100 * (i + 1), "depth": d}
        for i, d in enumerate([30, 30, 30, 60, 15])
    ]
    res = run(variants.calculate_copy_number_profile, bins)
    assert res["median_depth"] == 30
    assert [b["rounded_copy_number"] for b in res["bins"]] == [2, 2, 2, 4, 1]
    assert [b["log2_ratio"] for b in res["bins"]] == [0, 0, 0, 1, -1]
    assert [(s["copy_number"], s["n_bins"], s["state"]) for s in res["segments"]] == [
        (2, 3, "neutral"),
        (4, 1, "gain"),
        (1, 1, "loss"),
    ]
    pure = run(
        variants.calculate_copy_number_profile, bins[:2] + [{**bins[2], "depth": 45}], purity=0.5
    )
    # purity 0.5: CN = 4*ratio - 2 -> ratio 1.5 -> 4
    assert pure["bins"][2]["copy_number"] == pytest.approx(4.0)


# ---------------------------------------------------------------------------
# low-frequency, mosaic, repeats, MEI, density
# ---------------------------------------------------------------------------
def test_binomial_tail_matches_exact_sum():
    n, a, p = 60, 7, 0.05
    exact = sum(math.comb(n, k) * p**k * (1 - p) ** (n - k) for k in range(a, n + 1))
    assert _impl_analysis.binom_sf(a, n, p) == pytest.approx(exact, rel=1e-9)
    exact_cdf = sum(math.comb(n, k) * 0.5**n for k in range(0, 20 + 1))
    assert _impl_analysis.binom_cdf(20, n, 0.5) == pytest.approx(exact_cdf, rel=1e-9)


AD_VCF = vcf(
    ["T"],
    [
        "chr1 10 . A G . . . GT:AD 0/0:950,50",  # VAF .05, depth 1000 -> detected
        "chr1 20 . A G . . . GT:AD 0/0:999,1",  # alt 1 < 3
        "chr1 30 . A G . . . GT:AD 0/0:45,5",  # depth 50 < 100
        "chr1 40 . A G . . . GT:AD 0/1:700,300",  # VAF .3 > .2
    ],
)


def test_low_frequency_detection():
    res = run(variants.detect_low_frequency_variants, AD_VCF)
    assert [v["pos"] for v in res["variants"]] == [10] and res["n_tested"] == 3
    assert res["variants"][0]["vaf"] == 0.05 and res["variants"][0]["p_value_vs_error"] < 1e-30


def test_mosaic_detection():
    text = vcf(
        ["T"],
        [
            "chr1 1 . A G . . . GT:AD 0/1:80,20",
            "chr1 2 . A G . . . GT:AD 0/1:50,50",
            "chr1 3 . A G . . . GT:AD 0/0:99,1",
        ],
    )
    res = run(variants.detect_mosaic_variants, text)
    assert [v["pos"] for v in res["variants"]] == [
        1
    ]  # VAF .2 in range, p vs 0.5 tiny; VAF .5 and alt=1 rejected


def test_repeat_expansion_longest_run():
    seq = "TTCAGCAGCAGCAGCAGAAACAGCAG"  # CAG x5 at 3..17, then AAA, then CAG x2
    res = run(variants.detect_repeat_expansions, seq, motif="CAG", expansion_threshold=4)[
        "results"
    ][0]
    assert (res["longest_pure_run_units"], res["start"], res["end"], res["expansion"]) == (
        5,
        3,
        17,
        True,
    )
    assert (
        run(variants.detect_repeat_expansions, seq, motif="CAG", expansion_threshold=6)["results"][
            0
        ]["expansion"]
        is False
    )
    rc = run(variants.detect_repeat_expansions, "GGCTGCTGCTGG", motif="CAG", both_strands=True)[
        "results"
    ]
    assert [r["longest_pure_run_units"] for r in rc] == [
        0,
        3,
    ]  # CTG x3 is the reverse complement of CAG


def test_mobile_element_insertions():
    text = vcf(
        [],
        [
            "chr1 1000 mei1 N <INS:ME:ALU> 40 PASS SVTYPE=INS;SVLEN=300",
            "chr2 50 del1 N <DEL> 40 PASS SVTYPE=DEL;SVLEN=-500",
            "chr3 70 mei2 N <INS:ME:L1> 40 PASS SVTYPE=INS",
        ],
    )
    res = run(variants.detect_mobile_element_insertions, text)
    assert [(i["id"], i["family"]) for i in res["insertions"]] == [
        ("mei1", "Alu"),
        ("mei2", "LINE-1"),
    ]
    assert [
        i["id"]
        for i in run(variants.detect_mobile_element_insertions, text, families=["Alu"])[
            "insertions"
        ]
    ] == ["mei1"]


@pytest.mark.parametrize("fn", [variants.plot_density, variants.plot_density_mt])
def test_density_bins(fn):
    res = run(
        fn,
        recs((100, "A", "G"), (500, "A", "G"), (1500, "A", "G"), (2500, "A", "G")),
        bin_size=1000,
    )
    assert [(b["start"], b["end"], b["count"]) for b in res["bins"]] == [
        (1, 1000, 2),
        (1001, 2000, 1),
        (2001, 3000, 1),
    ]
    assert res["bins"][0]["density_per_mb"] == 2000 and res["svg"].startswith("<svg")


# ---------------------------------------------------------------------------
# coding effect
# ---------------------------------------------------------------------------
CDS = "ATGGCCAAATAG"  # ATG GCC AAA TAG = M A K *


@pytest.mark.parametrize(
    "pos,ref,alt,conseq,hgvs",
    [
        (5, "C", "T", "missense_variant", "p.Ala2Val"),  # GCC -> GTC
        (6, "C", "T", "synonymous_variant", "p.Ala2="),  # GCC -> GCT
        (7, "A", "T", "stop_gained", "p.Lys3Ter"),  # AAA -> TAA
        (1, "A", "G", "start_lost", "p.Met1?"),  # ATG -> GTG
        (10, "T", "C", "stop_lost", "p.Ter4Glnext*?"),  # TAG -> CAG
        (3, "GG", "G", "frameshift_variant", "p.Ala2ProfsTer?"),  # ATG CCA AAT AG..
        (3, "GGCC", "G", "inframe_deletion", "p.Ala2del"),  # ATG AAA TAG
        (6, "C", "CGGC", "inframe_insertion", "p.Ala2_Lys3insGly"),  # ATG GCC GGC AAA TAG
    ],
)
def test_coding_effect(pos, ref, alt, conseq, hgvs):
    res = run(
        variants.annotate_coding_effect,
        {"cds": CDS, "variants": [{"cds_position": pos, "ref": ref, "alt": alt}]},
    )
    effect = res["effects"][0]
    assert (effect["consequence"], effect["hgvs_p"]) == (conseq, hgvs)


def test_coding_effect_string_variant_and_codons():
    res = run(variants.annotate_coding_effect, {"cds": CDS, "variants": ["c.5C>T"]})
    assert res["effects"][0]["codons"][0] == {
        "codon_index": 2,
        "ref_codon": "GCC",
        "alt_codon": "GTC",
        "ref_aa": "A",
        "alt_aa": "V",
    }


# ---------------------------------------------------------------------------
# haplotypes, imputation, phasing
# ---------------------------------------------------------------------------
def test_construct_haplotypes_with_reference():
    text = vcf(["S1"], ["chr1 2 . A G . . . GT 1|0", "chr1 5 . C T . . . GT 0|1"])
    res = run(variants.construct_haplotypes, text, reference="TAGCCAGT")["haplotypes"]["S1"]["chr1"]
    assert [a["allele"] for a in res["1"]["alleles"]] == ["G", "C"] and [
        a["allele"] for a in res["2"]["alleles"]
    ] == ["A", "T"]
    assert (res["1"]["sequence"], res["2"]["sequence"]) == (
        "GGCC",
        "AGCT",
    )  # window 2..5 of TAGCCAGT is AGCC


def test_construct_haplotypes_unphased_is_error_or_skipped():
    text = vcf(["S1"], ["chr1 2 . A G . . . GT 1/0"])
    with pytest.raises(ConfigurationError):
        variants.construct_haplotypes(text)
    assert run(variants.construct_haplotypes, text, on_unphased="skip")["haplotypes"]["S1"] == {}


def test_fill_missing_genotypes_with_mode_modal():
    text = vcf(["A", "B", "C", "D"], ["chr1 1 . A G . . . GT 0/1 0/1 ./. 0/0"])
    res = run(variants.fill_missing_genotypes_with_mode, text)
    assert res["imputed_genotypes"] == 1 and res["records"][0]["genotypes"]["C"]["GT"] == "0/1"


def test_phase_genotypes_by_transmission():
    text = vcf(
        ["c", "f", "m"],
        [
            "chr1 1 . A G . . . GT 0/1 0/0 0/1",  # ALT can only be maternal -> 0|1
            "chr1 2 . A G . . . GT 0/1 0/1 0/0",  # ALT only paternal -> 1|0
            "chr1 3 . A G . . . GT 0/1 0/1 0/1",  # ambiguous
            "chr1 4 . A G . . . GT 0/1 0/0 0/0",  # impossible (de novo) -> error
            "chr1 5 . A G . . . GT 1/1 0/1 0/1",  # homozygous
        ],
    )
    res = run(variants.phase_genotypes, text, trio={"child": "c", "father": "f", "mother": "m"})
    assert res["counts"] == {"phased": 2, "ambiguous": 1, "mendelian_error": 1, "homozygous": 1}
    gts = [line.split("\t")[9] for line in res["vcf"].splitlines() if not line.startswith("#")]
    assert gts == ["0|1", "1|0", "0/1", "0/1", "1|1"]


# ---------------------------------------------------------------------------
# prioritisation
# ---------------------------------------------------------------------------
PRIORITY = [
    {
        "chrom": "1",
        "pos": 10,
        "ref": "A",
        "alt": "G",
        "gene": "BRCA1",
        "consequence": "missense_variant",
        "gnomad_af": 0.00005,
        "clinvar": "Pathogenic",
        "cadd": 30,
    },
    {"chrom": "1", "pos": 20, "ref": "A", "alt": "G", "gene": "X", "consequence": "stop_gained"},
    {
        "chrom": "1",
        "pos": 30,
        "ref": "A",
        "alt": "G",
        "gene": "Y",
        "consequence": "synonymous_variant",
        "gnomad_af": 0.2,
    },
    {
        "chrom": "1",
        "pos": 40,
        "ref": "A",
        "alt": "G",
        "gene": "Z",
        "consequence": "missense_variant",
        "clinvar": "Benign",
    },
]


@pytest.mark.parametrize("fn", [variants.prioritize_candidate_variants, clinical.prioritize])
def test_prioritize_scores(fn):
    res = run(fn, PRIORITY)
    # variant 1: missense .5 + ClinVar pathogenic 1.0 + AF<=1e-4 .5 + .5*30/40 = 2.375 ; variant 2: stop_gained 1.0 + absent AF .5 = 1.5
    assert [(r["pos"], r["score"]) for r in res["ranked"]] == [(10, 2.375), (20, 1.5)]
    assert {r["pos"] for r in res["excluded"]} == {30, 40}
    boosted = run(fn, PRIORITY, genes=["brca1"])
    assert boosted["ranked"][0]["score"] == 2.875


# ---------------------------------------------------------------------------
# clinical
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "criteria,label",
    [
        (["PVS1", "PS1"], "Pathogenic"),
        (["PVS1", "PM1", "PM2"], "Pathogenic"),
        (["PVS1", "PM2", "PP1"], "Pathogenic"),
        (["PVS1", "PP1", "PP2"], "Pathogenic"),
        (["PS1", "PS3"], "Pathogenic"),
        (["PS1", "PM1", "PM2", "PM3"], "Pathogenic"),
        (["PS1", "PM1", "PM2", "PP1", "PP2"], "Pathogenic"),
        (["PS1", "PM1", "PP1", "PP2", "PP3", "PP4"], "Pathogenic"),
        (["PVS1", "PM2"], "Likely pathogenic"),
        (["PS1", "PM1"], "Likely pathogenic"),
        (["PS1", "PP1", "PP2"], "Likely pathogenic"),
        (["PM1", "PM2", "PM3"], "Likely pathogenic"),
        (["PM1", "PM2", "PP1", "PP2"], "Likely pathogenic"),
        (["PM1", "PP1", "PP2", "PP3", "PP4"], "Likely pathogenic"),
        (
            ["PVS1", "PM2_Supporting"],
            "Uncertain significance",
        ),  # PVS1 + one supporting is not in Table 5
        (["PM1", "PP1", "PP2", "PP3"], "Uncertain significance"),
        (["PM2"], "Uncertain significance"),
        (["PS1"], "Uncertain significance"),
        (["BA1"], "Benign"),
        (["BS1", "BS2"], "Benign"),
        (["BS1", "BP4"], "Likely benign"),
        (["BP4", "BP7"], "Likely benign"),
        (["BP4"], "Uncertain significance"),
        (["PVS1", "PS1", "BA1"], "Uncertain significance"),  # contradictory evidence
        (["PM2_Supporting", "PP3"], "Uncertain significance"),
        (["PS1_Moderate", "PS2_Moderate", "PM1"], "Likely pathogenic"),  # three moderates
    ],
)
def test_acmg_combination_table(criteria, label):
    assert run(clinical.classify, criteria)["classification"] == label


def test_acmg_batch_and_normalisation():
    res = run(
        clinical.classify,
        {
            "variants": [
                {"id": "a", "criteria": ["pvs1", "PM2"]},
                {"id": "b", "criteria": [{"code": "BP4"}, "BP7"]},
            ]
        },
    )
    assert [(r["id"], r["classification"]) for r in res["results"]] == [
        ("a", "Likely pathogenic"),
        ("b", "Likely benign"),
    ]
    one = run(clinical.classify, ["PVS1", "PM2"])
    assert (
        "Likely pathogenic (i)" in one["rule"]
        and one["evidence_counts"]["pathogenic_very_strong"] == 1
    )


TRIO = vcf(
    ["child", "father", "mother"],
    [
        "chr1 100 . A G . PASS GENE=G1 GT:DP:GQ:AD 0/1:30:99:15,15 0/0:30:99:30,0 0/0:30:99:30,0",
        "chr1 150 . A G . PASS GENE=G6 GT:DP:GQ:AD 0/1:30:99:28,2 0/0:30:99:30,0 0/0:30:99:30,0",
        "chr1 160 . A G . PASS GENE=G6 GT:DP:GQ:AD 0/1:30:99:15,15 0/0:28:99:25,3 0/0:30:99:30,0",
        "chr1 170 . A G . PASS GENE=G6 GT:DP:GQ:AD 0/1:5:99:3,2 0/0:30:99:30,0 0/0:30:99:30,0",
        "chr1 200 . C T . PASS GENE=G2 GT:DP:GQ 0/1:30:99 0/1:30:99 0/0:30:99",
        "chr1 300 . G A . PASS GENE=G2 GT:DP:GQ 0/1:30:99 0/0:30:99 0/1:30:99",
        "chr1 400 . T C . PASS GENE=G3 GT:DP:GQ 1/1:30:99 0/1:30:99 0/1:30:99",
        "chr1 500 . A C . PASS GENE=G4 GT:DP:GQ 0/0:30:99 1/1:30:99 0/0:30:99",
        "chr1 600 . A T . PASS GENE=G5 GT:DP:GQ 1/1:30:99 0/0:30:99 0/1:30:99",
    ],
)


def test_trio_inheritance_patterns_and_compound_het():
    res = run(clinical.analyze_trio_variants, TRIO)
    # de novo: 100, 150, 160, 170 (all parents 0/0) ; paternal 200 ; maternal 300 ; homozygous 400 ; errors 500 and 600
    assert res["category_counts"] == {
        "de_novo": 4,
        "homozygous_both_parents_carriers": 1,
        "inherited_maternal": 1,
        "inherited_paternal": 1,
        "mendelian_error": 2,
    }
    assert res["compound_heterozygous_candidates"] == [
        {"gene": "G2", "paternal_variant": "chr1:200:C:T", "maternal_variant": "chr1:300:G:A"}
    ]
    assert any(
        "assuming column order" in w for w in clinical.analyze_trio_variants(TRIO)["warnings"]
    )


def test_detect_de_novo_filters():
    res = run(
        clinical.detect_de_novo_variants, TRIO, child="child", father="father", mother="mother"
    )
    # 150: child VAF 2/30 < 0.25 ; 160: father ALT reads 3 > 0 ; 170: child depth 5 < 10 ; 100 passes
    assert [v["pos"] for v in res["variants"]] == [100]
    assert res["rejected_by"] == {"child_low_depth": 1, "child_low_vaf": 1, "father_alt_reads": 1}
    assert res["n_candidates_before_filters"] == 4


def test_trio_requires_names_when_not_three_samples():
    four = vcf(["a", "b", "c", "d"], ["chr1 1 . A G . . . GT 0/1 0/0 0/0 0/0"])
    with pytest.raises(ConfigurationError):
        clinical.analyze_trio_variants(four)
    with pytest.raises(ConfigurationError):
        clinical.detect_de_novo_variants(
            four, options={"child": "a", "father": "b", "mother": "zzz"}
        )


@pytest.mark.parametrize(
    "fn,sites,detected",
    [
        (clinical.detect_circulating_tumor_dna, [(30, 10000), (25, 8000), (40, 12000)], True),
        (clinical.detect_circulating_tumor_dna, [(8, 10000), (7, 8000), (12, 12000)], False),
        (
            clinical.detect_minimal_residual_disease,
            [(40, 10000)],
            False,
        ),  # one positive site < 2 required
        (clinical.detect_minimal_residual_disease, [(40, 10000), (35, 9000)], True),
    ],
)
def test_ctdna_and_mrd(fn, sites, detected):
    res = run(fn, [{"alt_count": a, "depth": d} for a, d in sites])
    assert res["detected"] is detected


def test_ctdna_numbers():
    res = run(
        clinical.detect_circulating_tumor_dna,
        [
            {"alt_count": 30, "depth": 10000},
            {"alt_count": 25, "depth": 8000},
            {"alt_count": 40, "depth": 12000},
        ],
    )
    assert (res["pooled_alt_reads"], res["pooled_depth"]) == (95, 30000)
    assert res["pooled_vaf"] == pytest.approx(95 / 30000) and res[
        "estimated_tumor_fraction"
    ] == pytest.approx(2 * (95 / 30000 - 0.001))
    assert res["n_positive_sites"] == 3
    # exact upper tail P(X >= 95 | n=30000, p=0.001) in 60-digit decimal arithmetic (terms beyond k=400 are < 1e-80)
    decimal.getcontext().prec = 60
    p, q = decimal.Decimal(1) / 1000, decimal.Decimal(999) / 1000
    term = decimal.Decimal(math.comb(30000, 95)) * p**95 * q ** (30000 - 95)
    exact = decimal.Decimal(0)
    for k in range(95, 400):
        exact += term
        term = term * (30000 - k) * p / ((k + 1) * q)
    assert res["pooled_p_value"] == pytest.approx(float(exact), rel=1e-6)
    ctdna_one = run(clinical.detect_circulating_tumor_dna, [{"alt_count": 40, "depth": 10000}])
    assert ctdna_one["detected"] is True and ctdna_one["pooled_p_value"] < 1e-9
    assert (
        run(
            clinical.detect_circulating_tumor_dna,
            [
                {"alt_count": 8, "depth": 10000},
                {"alt_count": 7, "depth": 8000},
                {"alt_count": 12, "depth": 12000},
            ],
        )["estimated_tumor_fraction"]
        == 0.0
    )


def test_ctdna_from_vcf_ad():
    text = vcf(
        ["plasma"], ["chr1 1 . A G . . . GT:AD 0/0:9960,40", "chr1 2 . A G . . . GT:AD 0/0:9965,35"]
    )
    res = run(clinical.detect_circulating_tumor_dna, text)
    assert (
        res["pooled_alt_reads"] == 75 and res["pooled_depth"] == 20000 and res["detected"] is True
    )


def test_cancer_evolution_tree_nesting():
    mutations = {
        "A": [1, 1, 1, 0],
        "B": [1, 1, 0, 0],
        "C": [1, 0, 0, 0],
        "D": [0, 0, 1, 0],
        "E": [0, 0, 0, 1],
    }
    res = run(
        clinical.build_cancer_evolution_tree,
        {"mutations": mutations, "samples": ["S1", "S2", "S3", "S4"]},
    )
    assert res["perfect_phylogeny"] is True
    parent = {n["label"]: n["parent"] for n in res["nodes"]}
    ids = {n["label"]: n["id"] for n in res["nodes"]}
    assert (
        parent["B"] == ids["A"]
        and parent["C"] == ids["B"]
        and parent["D"] == ids["A"]
        and parent["E"] == "root"
        and parent["A"] == "root"
    )
    assert res["newick"] == "(((C)B,D)A,E)germline;"
    assert res["sample_to_node"] == {"S1": ids["C"], "S2": ids["B"], "S3": ids["D"], "S4": ids["E"]}
    same = run(clinical.build_cancer_evolution_tree, {"mutations": {"M1": [1, 1], "M2": [1, 1]}})
    assert [n["label"] for n in same["nodes"][1:]] == ["M1+M2"]


def test_cancer_evolution_tree_reports_four_gamete_violation():
    res = run(clinical.build_cancer_evolution_tree, {"mutations": {"X": [1, 1, 0], "Y": [0, 1, 1]}})
    assert res["perfect_phylogeny"] is False and res["tree"] is None and len(res["violations"]) == 1


def test_infer_tumor_clones_two_populations():
    sites = [
        {"id": f"m{i}", "alt_count": a, "depth": 100}
        for i, a in enumerate([50, 52, 48, 25, 24, 26])
    ]
    res = run(clinical.infer_tumor_clones, sites)
    assert res["n_clones"] == 2
    top, low = res["clones"]
    # hard-separated clusters: cluster VAF = pooled ALT fraction (150/300 and 75/300); CCF = 2*VAF for purity 1, CN 2
    assert top["cluster_vaf"] == pytest.approx(0.5, abs=1e-6) and low[
        "cluster_vaf"
    ] == pytest.approx(0.25, abs=1e-6)
    assert top["mean_ccf"] == pytest.approx(1.0) and low["mean_ccf"] == pytest.approx(0.5)
    assert top["mutations"] == ["m0", "m1", "m2"] and top["role"].startswith("clonal")
    # purity 0.5, CN 2, multiplicity 1: CCF = VAF * (0.5*2 + 2*0.5) / 0.5 = 4*VAF -> 2.0 and 1.0
    half = run(clinical.infer_tumor_clones, sites, purity=0.5)["clones"]
    assert (half[0]["mean_ccf"], half[1]["mean_ccf"]) == (pytest.approx(2.0), pytest.approx(1.0))


CLINVAR_VCF = vcf(
    [],
    [
        "chr17 43094464 rs1 A G . . CLNSIG=Pathogenic/Likely_pathogenic;GENEINFO=BRCA1:672|BRCA1-AS1:100;"
        "CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Breast-ovarian_cancer|Other;MC=SO:0001583|missense_variant"
    ],
)


def test_clinical_read():
    row = run(clinical.read, CLINVAR_VCF)["variants"][0]
    assert row["gene"] == "BRCA1" and row["clinvar_terms"] == ["pathogenic", "likely_pathogenic"]
    assert row["disease"] == ["Breast-ovarian_cancer", "Other"] and row["consequence"] == [
        "missense_variant"
    ]
    assert row["review_status"] == "reviewed_by_expert_panel"


@pytest.mark.parametrize("fn", [clinical.generate_clinical_variant_report, clinical.report])
def test_clinical_report(fn):
    data = {
        "sample": "NA12878",
        "variants": [
            {
                "gene": "BRCA1",
                "hgvs": "NM_007294.4:c.5266dup",
                "zygosity": "heterozygous",
                "classification": "Pathogenic",
            },
            {"gene": "TTN", "hgvs": "c.1A>G", "criteria": ["PM2_Supporting", "PP3"]},
            {"gene": "APOE", "hgvs": "c.388T>C", "classification": "LB"},
        ],
    }
    res = run(fn, data, title="Test report")
    assert res["classification_counts"] == {
        "Pathogenic": 1,
        "Uncertain significance": 1,
        "Likely benign": 1,
    }
    text = res["report"]
    assert text.startswith("# Test report") and "Sample: NA12878" in text
    assert "| BRCA1 | NM_007294.4:c.5266dup | heterozygous | Pathogenic |" in text
    assert "Variants of uncertain significance" in text and "1 benign / likely benign" in text
    assert res["n_reportable"] == 1
    assert "Variants of uncertain significance" not in run(fn, data, include_vus=False)["report"]


def test_report_from_vcf_with_hgvs():
    res = run(clinical.report, CLINVAR_VCF, hgvs_accession="NC_000017.11")
    assert (
        res["variants"][0]["hgvs_g"] == "NC_000017.11:g.43094464A>G"
        and res["variants"][0]["gene"] == "BRCA1"
    )


def test_clinical_qc_flags():
    res = run(clinical.qc, VCF_A)
    flagged = {s["sample"]: s["flags"] for s in res["samples"]}
    assert flagged == {"S1": [], "S2": [], "S3": ["low_call_rate"]}  # S3 call rate 3/4 < 0.95
    assert [s["pos"] for s in res["flagged_sites"]] == [200] and res[
        "passed"
    ] is False  # site 200 call rate 2/3 < 0.9


# ---------------------------------------------------------------------------
# leftovers stay honest
# ---------------------------------------------------------------------------
NOT_IMPLEMENTED = [
    (variants, n)
    for n in [
        "annotate",
        "annotate_conservation_score",
        "annotate_mt",
        "annotate_pathogenicity_score",
        "annotate_population_frequency",
        "annotate_regulatory_effect",
        "annotate_splice_effect",
        "annotate_variant_consequence",
        "call",
        "call_mt",
        "detect_breakpoints",
        "genotype_structural_variants",
    ]
] + [(clinical, n) for n in ["annotate", "calculate", "call", "compare", "predict", "type_hla"]]


@pytest.mark.parametrize("module,name", NOT_IMPLEMENTED)
def test_unavailable_operations_report_not_implemented(module, name):
    out = getattr(module, name)("chr1\t1\t.\tA\tG")
    assert out["status"] == "not_implemented"
    assert out["required_dependency"] and out["message"] and "result" not in out


# ---------------------------------------------------------------------------
# invalid input
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fn,value,options,error",
    [
        (variants.normalize, "not a vcf", None, InputFormatError),
        (variants.summarize, 12345, None, InputFormatError),
        (
            variants.left_align,
            recs((5, "ACA", "A")),
            None,
            ConfigurationError,
        ),  # reference required
        (variants.calculate_hwe_pvalue, [-1, 0, 1], None, InputFormatError),
        (variants.merge, [MERGE_A], None, InputFormatError),
        (variants.merge, [MERGE_A, MERGE_A], None, InputFormatError),  # duplicate sample S1
        (variants.filter, VCF_A, {"min_qual": "high"}, ConfigurationError),
        (variants.filter, VCF_A, {"variant_types": ["weird"]}, ConfigurationError),
        (variants.detect_repeat_expansions, "CAGCAG", None, ConfigurationError),  # motif missing
        (
            variants.annotate_coding_effect,
            {"cds": "ATGGCC", "variants": [{"cds_position": 4, "ref": "T", "alt": "A"}]},
            None,
            InputFormatError,
        ),
        (
            variants.annotate_coding_effect,
            {"variants": [{"cds_position": 4, "ref": "G", "alt": "A"}]},
            None,
            InputFormatError,
        ),
        (
            variants.calculate_ld_matrix,
            recs((1, "A", "G"), (2, "A", "G")),
            None,
            InputFormatError,
        ),  # no genotypes
        (variants.phase_genotypes, VCF_A, None, ConfigurationError),  # trio missing
        (clinical.classify, ["PVS2"], None, InputFormatError),
        (clinical.classify, ["XX1"], None, InputFormatError),
        (clinical.classify, ["BA1_VeryStrong"], None, InputFormatError),
        (clinical.classify, None, None, InputFormatError),
        (
            clinical.build_cancer_evolution_tree,
            {"mutations": {"A": [1, 2]}},
            None,
            InputFormatError,
        ),
        (
            clinical.detect_circulating_tumor_dna,
            [{"alt_count": 5, "depth": 2}],
            None,
            InputFormatError,
        ),
        (
            clinical.report,
            {"variants": [{"gene": "G", "classification": "maybe"}]},
            None,
            InputFormatError,
        ),
    ],
)
def test_invalid_input_raises(fn, value, options, error):
    with pytest.raises(error) as info:
        fn(value, options=options)
    assert getattr(info.value, "code", None) and getattr(info.value, "hint", None)
