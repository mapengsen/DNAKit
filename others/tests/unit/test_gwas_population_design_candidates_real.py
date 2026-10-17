"""Real-implementation tests for the gwas, population, breeding, forensic and design candidates.

Every implemented function is checked against a hand-computed value, a textbook formula, or
scipy / numpy / statsmodels / scikit-learn computed independently of the DNAKit code.
"""

from __future__ import annotations

import gzip
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
scipy_stats = pytest.importorskip("scipy.stats")

import dnakit.breeding as breeding  # noqa: E402
import dnakit.experimental_design as design  # noqa: E402
import dnakit.forensic_reproductive as forensic  # noqa: E402
import dnakit.gwas as gwas  # noqa: E402
import dnakit.population as population  # noqa: E402
from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError  # noqa: E402
from dnakit.gwas import _impl_stats as st  # noqa: E402


def ok(envelope: dict) -> dict:
    assert envelope["status"] == "ok", envelope
    json.dumps(envelope, allow_nan=False)  # the Agent/MCP layer needs strict JSON
    return envelope["result"]


# =============================================================================== numerical helpers
def test_numerical_helpers_match_scipy() -> None:
    from scipy import special

    assert st.normal_ppf(0.975) == pytest.approx(scipy_stats.norm.ppf(0.975), abs=1e-12)
    assert st.normal_ppf(1e-10) == pytest.approx(scipy_stats.norm.ppf(1e-10), abs=1e-9)
    assert st.chi2_sf(7.3, 4) == pytest.approx(scipy_stats.chi2.sf(7.3, 4), rel=1e-10)
    assert st.chi2_ppf(0.95, 3) == pytest.approx(scipy_stats.chi2.ppf(0.95, 3), rel=1e-9)
    assert st.t_two_sided_p(2.3, 11) == pytest.approx(2 * scipy_stats.t.sf(2.3, 11), rel=1e-10)
    assert st.t_ppf(0.975, 8) == pytest.approx(scipy_stats.t.ppf(0.975, 8), rel=1e-9)
    assert st.betainc_regularized(2.5, 1.5, 0.3) == pytest.approx(
        special.betainc(2.5, 1.5, 0.3), rel=1e-10
    )
    assert st.ncx2_sf(10.0, 3, 6.5) == pytest.approx(scipy_stats.ncx2.sf(10.0, 3, 6.5), rel=1e-8)
    tc = scipy_stats.t.ppf(0.975, 38)
    expect = scipy_stats.nct.sf(tc, 38, 1.4) + scipy_stats.nct.cdf(-tc, 38, 1.4)
    assert st.nct_power(tc, 38, 1.4) == pytest.approx(expect, abs=1e-7)
    assert st.neglog10_p_from_z(40.0) == pytest.approx(
        -scipy_stats.norm.logsf(40.0) / math.log(10) - math.log10(2), abs=1e-6
    )


# ==================================================================================== gwas: tables
TSV = (
    "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tEAF\tN\n"
    "rs1\t1\t100\tA\tG\t0.5\t0.1\t5.7e-7\t0.3\t1000\n"
    "rs2\tchr2\t200\tC\tT\t-0.2\t0.1\tNA\t0.5\t1000\n"
)


def test_gwas_read_text_and_gzip(tmp_path: Path) -> None:
    result = ok(gwas.read(TSV))
    assert result["n_records"] == 2 and result["columns"][0] == "SNP"
    assert result["records"][0]["BETA"] == 0.5 and result["records"][1]["P"] is None
    path = tmp_path / "ss.tsv.gz"
    path.write_bytes(gzip.compress(TSV.encode()))
    assert ok(gwas.read(str(path), options={"standardize": True}))["records"][1]["chrom"] == "2"


def test_gwas_standardize_derives_columns() -> None:
    result = ok(
        gwas.standardize(
            [
                {
                    "rsid": "rsX",
                    "Chr": "chr23",
                    "pos": 5,
                    "ea": "a",
                    "nea": "g",
                    "OR": 2.0,
                    "se": 0.25,
                }
            ]
        )
    )
    rec = result["records"][0]
    assert rec["chrom"] == "X" and rec["effect_allele"] == "A" and rec["other_allele"] == "G"
    assert rec["beta"] == pytest.approx(math.log(2.0))
    assert rec["z"] == pytest.approx(math.log(2.0) / 0.25)
    assert rec["p"] == pytest.approx(2 * scipy_stats.norm.sf(abs(rec["z"])), rel=1e-9)
    only_p = ok(
        gwas.standardize(
            {"chr": [1], "bp": [10], "a1": ["A"], "a2": ["T"], "beta": [0.3], "p": [1e-4]}
        )
    )["records"][0]
    assert only_p["variant_id"] == "1:10:T:A"
    assert only_p["se"] == pytest.approx(0.3 / scipy_stats.norm.isf(0.5e-4), rel=1e-8)


def test_gwas_validate_reports_problems() -> None:
    good = ok(gwas.validate(TSV))
    assert good["valid"] is True
    bad = ok(
        gwas.validate(
            [
                {
                    "snp": "a",
                    "chr": 1,
                    "bp": 1,
                    "a1": "A",
                    "a2": "A",
                    "beta": 1,
                    "se": -1,
                    "p": 1.5,
                    "eaf": 2,
                },
                {"snp": "a", "chr": 1, "bp": 2, "a1": "A", "a2": "C", "beta": 1, "se": 1, "p": 0.5},
            ]
        )
    )
    assert bad["valid"] is False
    joined = " ".join(bad["errors"])
    assert (
        "outside [0, 1]" in joined
        and "standard errors" in joined
        and "identical" in joined
        and "frequencies" in joined
    )


def test_gwas_filter_counts() -> None:
    rows = [
        {
            "snp": f"s{i}",
            "chr": 1 + i // 3,
            "bp": 100 * (i + 1),
            "a1": "A",
            "a2": "AT" if i == 4 else "G",
            "beta": 0.1,
            "se": 0.1,
            "p": p,
            "eaf": e,
            "n": 100 + i,
        }
        for i, (p, e) in enumerate(
            [(0.001, 0.4), (0.2, 0.4), (0.04, 0.01), (0.03, 0.5), (0.5, 0.3), (1e-9, 0.2)]
        )
    ]
    result = ok(gwas.filter(rows, options={"max_p": 0.05, "min_maf": 0.05, "snp_only": True}))
    kept = {r["variant_id"] for r in result["records"]}
    assert kept == {"s0", "s3", "s5"}
    assert result["removed"] == {"p_above_threshold": 2, "maf_below_threshold": 1}
    region = ok(gwas.filter(rows, options={"region": "1:100-200"}))
    assert {r["variant_id"] for r in region["records"]} == {"s0", "s1"}
    with pytest.raises(ConfigurationError):
        gwas.filter(rows, options={"region": "nonsense"})


def test_gwas_harmonize_actions() -> None:
    ref = [
        {"snp": "a", "chr": 1, "bp": 1, "a1": "A", "a2": "G", "beta": 0, "se": 1, "p": 1},
        {"snp": "b", "chr": 1, "bp": 2, "a1": "A", "a2": "G", "beta": 0, "se": 1, "p": 1},
        {"snp": "c", "chr": 1, "bp": 3, "a1": "A", "a2": "G", "beta": 0, "se": 1, "p": 1},
        {
            "snp": "d",
            "chr": 1,
            "bp": 4,
            "a1": "A",
            "a2": "T",
            "beta": 0,
            "se": 1,
            "p": 1,
            "eaf": 0.8,
        },
        {
            "snp": "e",
            "chr": 1,
            "bp": 5,
            "a1": "A",
            "a2": "T",
            "beta": 0,
            "se": 1,
            "p": 1,
            "eaf": 0.5,
        },
        {"snp": "f", "chr": 1, "bp": 6, "a1": "A", "a2": "G", "beta": 0, "se": 1, "p": 1},
    ]
    tgt = [
        {
            "snp": "a",
            "chr": 1,
            "bp": 1,
            "a1": "G",
            "a2": "A",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.2,
        },  # swapped
        {
            "snp": "b",
            "chr": 1,
            "bp": 2,
            "a1": "T",
            "a2": "C",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.2,
        },  # strand flip
        {
            "snp": "c",
            "chr": 1,
            "bp": 3,
            "a1": "A",
            "a2": "C",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.2,
        },  # mismatch
        {
            "snp": "d",
            "chr": 1,
            "bp": 4,
            "a1": "T",
            "a2": "A",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.2,
        },  # palindromic, resolved by eaf
        {
            "snp": "e",
            "chr": 1,
            "bp": 5,
            "a1": "A",
            "a2": "T",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.5,
        },  # ambiguous palindrome
        {
            "snp": "f",
            "chr": 1,
            "bp": 6,
            "a1": "A",
            "a2": "G",
            "beta": 0.5,
            "se": 0.1,
            "p": 0.1,
            "eaf": 0.2,
        },
    ]  # identical
    result = ok(gwas.harmonize({"target": tgt, "reference": ref}))
    out = {r["variant_id"]: r for r in result["records"]}
    assert (
        out["a"]["beta"] == -0.5
        and out["a"]["eaf"] == pytest.approx(0.8)
        and out["a"]["effect_allele"] == "A"
    )
    assert out["a"]["harmonize_action"] == "swap"
    assert (
        out["b"]["beta"] == 0.5
        and out["b"]["effect_allele"] == "A"
        and out["b"]["harmonize_action"] == "strand_flip"
    )
    # palindromic T/A (eaf 0.2) vs reference A/T (eaf 0.8): same strand, labels swapped -> flip
    assert (
        out["d"]["beta"] == -0.5
        and out["d"]["harmonize_action"] == "swap"
        and out["d"]["eaf"] == pytest.approx(0.8)
    )
    assert "c" not in out and "e" not in out
    assert out["f"]["beta"] == 0.5 and out["f"]["harmonize_action"] == "none"
    assert (
        result["actions"]["allele_mismatch"] == 1 and result["actions"]["ambiguous_palindrome"] == 1
    )
    # reference given as an option
    assert ok(gwas.harmonize(tgt, options={"reference": ref}))["n_harmonized"] == 4
    with pytest.raises(ConfigurationError):
        gwas.harmonize(tgt)


# ============================================================================ gwas: QC / plot data
Z_VALUES = [0.3, -0.8, 1.1, 0.5, -1.9, 2.5, 0.1, -0.2, 3.2]


def _z_table() -> list[dict]:
    return [
        {
            "snp": f"s{i}",
            "chr": 1,
            "bp": 1000 * (i + 1),
            "a1": "A",
            "a2": "G",
            "z": z,
            "n": 500,
            "eaf": 0.3,
        }
        for i, z in enumerate(Z_VALUES)
    ]


def test_gwas_calculate_lambda_gc_and_adjusted_p() -> None:
    lam = float(np.median(np.square(Z_VALUES)) / scipy_stats.chi2.ppf(0.5, 1))
    summary = ok(gwas.calculate(_z_table()))
    assert summary["lambda_gc"] == pytest.approx(lam, rel=1e-9)
    assert summary["n_variants"] == 9 and summary["lead_variant"] == "s8"
    assert ok(gwas.calculate(_z_table(), options={"statistic": "lambda_gc"}))[
        "value"
    ] == pytest.approx(lam, rel=1e-9)
    assert (
        ok(
            gwas.calculate(_z_table(), options={"statistic": "n_significant", "significance": 0.05})
        )["value"]
        == 2
    )  # |z| >= 1.96: 2.5, 3.2
    sm = pytest.importorskip("statsmodels.stats.multitest")
    p = [2 * scipy_stats.norm.sf(abs(z)) for z in Z_VALUES]
    adj = ok(gwas.calculate(_z_table(), options={"statistic": "adjusted_p"}))["values"]
    assert [v["p_bonferroni"] for v in adj] == pytest.approx(
        sm.multipletests(p, method="bonferroni")[1], rel=1e-9
    )
    assert [v["q_bh"] for v in adj] == pytest.approx(
        sm.multipletests(p, method="fdr_bh")[1], rel=1e-9
    )


def test_gwas_apply_genomic_control() -> None:
    table = [{**r, "z": 2 * r["z"]} for r in _z_table()]
    lam = float(np.median(np.square([2 * z for z in Z_VALUES])) / scipy_stats.chi2.ppf(0.5, 1))
    assert lam > 1
    result = ok(gwas.apply_genomic_control(table))
    assert result["lambda_gc"] == pytest.approx(lam, rel=1e-9)
    for rec, z in zip(result["records"], [2 * z for z in Z_VALUES], strict=False):
        assert rec["p_gc"] == pytest.approx(scipy_stats.chi2.sf(z * z / lam, 1), rel=1e-8)
        assert rec["z_gc"] == pytest.approx(z / math.sqrt(lam))
    shrunk = [
        {**r, "z": r["z"] / 3} for r in _z_table()
    ]  # lambda = 1.4068 / 9 < 1 is not applied by default
    deflated = ok(gwas.apply_genomic_control(shrunk))
    assert deflated["lambda_gc"] < 1 and deflated["lambda_used"] == 1.0
    assert deflated["records"][0]["p_gc"] == pytest.approx(
        2 * scipy_stats.norm.sf(abs(Z_VALUES[0] / 3)), rel=1e-9
    )
    forced = ok(gwas.apply_genomic_control(shrunk, options={"only_inflated": False}))
    assert forced["lambda_used"] == pytest.approx(deflated["lambda_gc"])


def test_gwas_compute_qq_data() -> None:
    data = ok(
        gwas.compute_qq_data([{"snp": i, "p": p} for i, p in enumerate([0.5, 0.01, 0.2, 0.001])])
    )
    pts = sorted(data["points"], key=lambda r: r["rank"])
    assert [round(x["observed"], 6) for x in pts] == [
        3.0,
        2.0,
        round(-math.log10(0.2), 6),
        round(-math.log10(0.5), 6),
    ]
    assert [x["expected"] for x in pts] == pytest.approx(
        [-math.log10(v / 4) for v in (0.5, 1.5, 2.5, 3.5)]
    )
    thin = ok(
        gwas.compute_qq_data(
            [{"p": (i + 1) / 1001} for i in range(1000)], options={"max_points": 100}
        )
    )
    assert len(thin["points"]) == 100 and thin["downsampled"]


def test_gwas_compute_manhattan_data() -> None:
    table = [
        {"chr": 1, "bp": 100, "p": 1e-3},
        {"chr": 1, "bp": 500, "p": 0.5},
        {"chr": 2, "bp": 50, "p": 1e-9},
        {"chr": "X", "bp": 10, "p": 0.1},
    ]
    data = ok(gwas.compute_manhattan_data(table))
    assert [c["offset"] for c in data["chromosomes"]] == [0, 500, 550]
    xs = [p["x"] for p in data["points"]]
    assert xs == [100, 500, 550, 560]
    assert data["points"][2]["neglog10p"] == pytest.approx(9.0)
    assert data["thresholds"]["genome_wide"]["neglog10p"] == pytest.approx(-math.log10(5e-8))
    assert len(ok(gwas.compute_manhattan_data(table, options={"min_neglog10p": 2}))["points"]) == 2
    man = ok(gwas.compute_manhattan_data([{"chr": 1, "bp": 5, "p": 0.0}]))
    underflow = man["points"][0]["neglog10p"]
    assert underflow == pytest.approx(-math.log10(5e-324)) and math.isfinite(underflow)


def test_gwas_locus_plot_data_with_ld() -> None:
    table = [
        {"snp": f"v{i}", "chr": 3, "bp": 1000 * i, "p": p}
        for i, p in enumerate([0.1, 1e-8, 0.01, 0.5], start=1)
    ]
    table.append({"snp": "far", "chr": 3, "bp": 10_000_000, "p": 1e-3})
    ld = {
        "variant_ids": ["v1", "v2", "v3", "v4"],
        "matrix": [[1, 0.9, 0.5, 0.1], [0.9, 1, 0.7, 0.3], [0.5, 0.7, 1, 0], [0.1, 0.3, 0, 1]],
    }
    data = ok(gwas.compute_locus_plot_data(table, options={"flank": 5000, "ld": ld}))
    assert data["lead"]["variant_id"] == "v2" and len(data["points"]) == 4
    bands = {p["variant_id"]: (p["ld_band"], p["r2"]) for p in data["points"]}
    assert bands["v2"] == ("lead", 1.0)
    assert bands["v1"][0] == "0.8-1.0" and bands["v1"][1] == pytest.approx(0.81)
    assert bands["v3"][0] == "0.4-0.6" and bands["v4"][0] == "0.0-0.2"
    with pytest.raises(ConfigurationError):
        gwas.compute_locus_plot_data(table, options={"lead": "missing"})


# ====================================================================================== gwas: loci
LOCI = [
    {"snp": "a", "chr": 1, "bp": 1_000_000, "p": 1e-10, "beta": 1, "se": 0.1},
    {"snp": "b", "chr": 1, "bp": 1_200_000, "p": 1e-9, "beta": 1, "se": 0.1},
    {"snp": "c", "chr": 1, "bp": 3_000_000, "p": 1e-12, "beta": 1, "se": 0.1},
    {"snp": "d", "chr": 1, "bp": 3_100_000, "p": 1e-8, "beta": 1, "se": 0.1},
    {"snp": "e", "chr": 2, "bp": 100, "p": 1e-9, "beta": 1, "se": 0.1},
    {"snp": "f", "chr": 2, "bp": 200, "p": 0.2, "beta": 1, "se": 0.1},
]


def test_gwas_select_lead_snps_and_merge_and_loci() -> None:
    leads = ok(gwas.select_lead_snps(LOCI, options={"window": 500_000}))
    assert [lead["variant_id"] for lead in leads["leads"]] == ["a", "c", "e"]
    assert leads["leads"][0]["n_nearby_significant"] == 1
    merged = ok(gwas.merge_nearby_significant_hits(LOCI, options={"distance": 500_000}))
    assert [
        (g["chrom"], g["start"], g["end"], g["n_hits"], g["lead_variant"]) for g in merged["groups"]
    ] == [
        ("1", 1_000_000, 1_200_000, 2, "a"),
        ("1", 3_000_000, 3_100_000, 2, "c"),
        ("2", 100, 100, 1, "e"),
    ]
    loci = ok(gwas.define_association_loci(LOCI, options={"flank": 1_000_000}))
    # windows [1, 2e6] and [2e6, 4e6] touch at 2,000,000, so the leads a and c are merged
    assert [(x["chrom"], x["start"], x["end"], x["lead_variant"]) for x in loci["loci"]] == [
        ("1", 1, 4_000_000, "c"),
        ("2", 1, 1_000_100, "e"),
    ]
    assert loci["loci"][0]["n_significant"] == 4 and loci["loci"][0]["leads"] == ["a", "c"]
    narrow = ok(gwas.define_association_loci(LOCI, options={"flank": 100_000}))
    assert narrow["n_loci"] == 3


def test_gwas_annotate_lead_snps() -> None:
    genes = [
        {"name": "GENE1", "chrom": 1, "start": 900_000, "end": 1_050_000, "strand": "+"},
        {"name": "GENE2", "chrom": 1, "start": 3_200_000, "end": 3_300_000, "strand": "-"},
        {"name": "FAR", "chrom": 2, "start": 9_000_000, "end": 9_100_000},
    ]
    out = ok(gwas.annotate_lead_snps(LOCI, options={"genes": genes, "window": 500_000}))
    ann = {a["variant_id"]: a for a in out["annotations"]}
    assert (
        ann["a"]["nearest_gene"] == "GENE1"
        and ann["a"]["relation"] == "within_gene"
        and ann["a"]["distance_to_nearest_gene"] == 0
    )
    assert (
        ann["c"]["nearest_gene"] == "GENE2" and ann["c"]["distance_to_nearest_gene"] == 200_000
    )  # 3,200,000 - 3,000,000
    assert (
        ann["c"]["relation"] == "downstream"
    )  # minus-strand gene: positions before its start are downstream
    assert ann["e"]["nearest_gene"] is None
    assert ann["c"]["genes_within_flank"] == []  # default flank 100 kb < 200 kb
    wide = ok(
        gwas.annotate_lead_snps(
            LOCI, options={"genes": genes, "window": 500_000, "flank": 250_000, "leads": ["c"]}
        )
    )
    assert wide["annotations"][0]["genes_within_flank"] == ["GENE2"]
    with pytest.raises(ConfigurationError):
        gwas.annotate_lead_snps(LOCI)


def test_gwas_build_credible_set_matches_bayes_factors() -> None:
    betas, ses = [0.5, 0.45, 0.1, 0.02], [0.05, 0.05, 0.05, 0.05]
    table = [
        {"snp": f"v{i}", "beta": b, "se": s}
        for i, (b, s) in enumerate(zip(betas, ses, strict=False))
    ]
    w = 0.15**2
    bf = [
        scipy_stats.norm.pdf(b, 0, math.sqrt(s * s + w)) / scipy_stats.norm.pdf(b, 0, s)
        for b, s in zip(betas, ses, strict=False)
    ]
    pip = [x / sum(bf) for x in bf]
    result = ok(gwas.build_credible_set(table))
    got = {v["variant_id"]: v["pip"] for v in result["variants"]}
    assert [got[f"v{i}"] for i in range(4)] == pytest.approx(pip, rel=1e-9)
    order = sorted(range(4), key=lambda i: -pip[i])
    cum, expected = 0.0, []
    for i in order:
        expected.append(f"v{i}")
        cum += pip[i]
        if cum >= 0.95:
            break
    assert result["credible_set"] == expected
    assert ok(gwas.build_credible_set(table, options={"coverage": 0.5}))["set_size"] == 1
    with pytest.raises(InputFormatError):
        gwas.build_credible_set(table[:1])


# ============================================================================= gwas: meta-analysis
STUDIES = {
    "s1": [
        {"snp": "x", "a1": "A", "a2": "G", "beta": 0.30, "se": 0.10},
        {"snp": "y", "a1": "A", "a2": "G", "beta": 0.10, "se": 0.10},
    ],
    "s2": [
        {"snp": "x", "a1": "A", "a2": "G", "beta": 0.10, "se": 0.05},
        {"snp": "y", "a1": "G", "a2": "A", "beta": -0.12, "se": 0.08},
    ],
    "s3": [
        {"snp": "x", "a1": "A", "a2": "G", "beta": 0.45, "se": 0.12},
        {"snp": "y", "a1": "A", "a2": "G", "beta": 0.0, "se": 0.2},
    ],
}


def test_gwas_meta_analysis_matches_statsmodels() -> None:
    meta = pytest.importorskip("statsmodels.stats.meta_analysis")
    fixed = {r["variant_id"]: r for r in ok(gwas.fixed_effect_meta_analysis(STUDIES))["records"]}
    random_ = {r["variant_id"]: r for r in ok(gwas.random_effect_meta_analysis(STUDIES))["records"]}
    betas = {
        "x": ([0.30, 0.10, 0.45], [0.10, 0.05, 0.12]),
        "y": ([0.10, 0.12, 0.0], [0.10, 0.08, 0.2]),
    }
    for key, (b, s) in betas.items():
        res_f = meta.combine_effects(np.array(b), np.square(s), method_re="dl")
        frame = res_f.summary_frame()
        assert fixed[key]["beta"] == pytest.approx(frame.loc["fixed effect", "eff"], rel=1e-9)
        assert fixed[key]["se"] == pytest.approx(frame.loc["fixed effect", "sd_eff"], rel=1e-9)
        if res_f.tau2 > 0:
            assert random_[key]["beta"] == pytest.approx(
                frame.loc["random effect", "eff"], rel=1e-9
            )
            assert random_[key]["se"] == pytest.approx(
                frame.loc["random effect", "sd_eff"], rel=1e-9
            )
            assert random_[key]["tau2"] == pytest.approx(res_f.tau2, rel=1e-9)
        else:  # statsmodels keeps a negative tau^2; DerSimonian-Laird truncates it at 0
            assert random_[key]["tau2"] == 0.0
            assert random_[key]["beta"] == pytest.approx(fixed[key]["beta"], rel=1e-12)
        assert fixed[key]["q"] == pytest.approx(res_f.q, rel=1e-9)
        assert fixed[key]["i2"] == pytest.approx(max(0.0, (res_f.q - 2) / res_f.q), rel=1e-9)
        assert fixed[key]["p"] == pytest.approx(
            2 * scipy_stats.norm.sf(abs(fixed[key]["beta"] / fixed[key]["se"])), rel=1e-9
        )
    # the allele-swapped study-2 record for 'y' was sign-flipped before pooling (0.12 above)
    stacked = [{"study": name, **rec} for name, recs in STUDIES.items() for rec in recs]
    assert ok(gwas.fixed_effect_meta_analysis(stacked))["records"][0]["beta"] == pytest.approx(
        fixed["x"]["beta"]
    )
    with pytest.raises(InputFormatError):
        gwas.fixed_effect_meta_analysis({"only": STUDIES["s1"]})


# ========================================================================== gwas: individual-level
def _genotype_data(
    n: int = 120, m: int = 5, seed: int = 3
) -> tuple[list[list[int]], list[float], list[float]]:
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, [0.2, 0.35, 0.5, 0.15, 0.4][:m], size=(n, m))
    cov = rng.normal(size=n)
    y = 0.4 * g[:, 0] - 0.3 * g[:, 2] + 0.5 * cov + rng.normal(size=n)
    return g.tolist(), y.tolist(), cov.tolist()


def test_gwas_run_linear_matches_scipy_and_covariates_match_numpy() -> None:
    g, y, cov = _genotype_data()
    data = {"genotypes": g, "phenotype": y, "variant_ids": [f"v{j}" for j in range(5)]}
    res = ok(gwas.run(data))
    assert res["model"] == "linear"
    garr = np.array(g, dtype=float)
    for j, rec in enumerate(res["records"]):
        ref = scipy_stats.linregress(garr[:, j], y)
        assert rec["beta"] == pytest.approx(ref.slope, rel=1e-9) and rec["se"] == pytest.approx(
            ref.stderr, rel=1e-9
        )
        assert rec["p"] == pytest.approx(ref.pvalue, rel=1e-8) and rec["eaf"] == pytest.approx(
            garr[:, j].mean() / 2
        )
    adj = ok(gwas.run({**data, "covariates": [[c] for c in cov]}))
    x = np.column_stack([np.ones(len(y)), np.array(cov), garr[:, 0]])
    coef, *_ = np.linalg.lstsq(x, np.array(y), rcond=None)
    resid = np.array(y) - x @ coef
    sigma2 = resid @ resid / (len(y) - 3)
    se = math.sqrt(sigma2 * np.linalg.inv(x.T @ x)[2, 2])
    assert adj["records"][0]["beta"] == pytest.approx(coef[2], rel=1e-9)
    assert adj["records"][0]["se"] == pytest.approx(se, rel=1e-9)
    assert adj["records"][0]["p"] == pytest.approx(
        2 * scipy_stats.t.sf(abs(coef[2] / se), len(y) - 3), rel=1e-7
    )


def test_gwas_run_logistic_matches_statsmodels_and_handles_missing() -> None:
    sm = pytest.importorskip("statsmodels.api")
    g, y, _ = _genotype_data(n=300)
    status = [1 if v > np.median(y) else 0 for v in y]
    g[0][1] = None
    res = ok(gwas.run({"genotypes": g, "phenotype": status}))
    assert res["model"] == "logistic"
    garr = np.array([[np.nan if v is None else v for v in row] for row in g], dtype=float)
    for j in (0, 1, 2):
        keep = ~np.isnan(garr[:, j])
        fit = sm.Logit(np.array(status)[keep], sm.add_constant(garr[keep, j])).fit(disp=0)
        rec = res["records"][j]
        assert rec["beta"] == pytest.approx(fit.params[1], rel=1e-6) and rec["se"] == pytest.approx(
            fit.bse[1], rel=1e-6
        )
        assert rec["p"] == pytest.approx(fit.pvalues[1], rel=1e-5) and rec["n"] == int(keep.sum())
    mono = ok(gwas.run({"genotypes": [[1, 0]] * 6 + [[1, 2]] * 6, "phenotype": list(range(12))}))
    assert mono["records"][0]["status"] == "monomorphic_or_too_few_samples"


@pytest.mark.parametrize("backend", ["plink", "regenie", "gemma", "gcta"])
def test_gwas_run_external_backends_are_not_implemented(backend: str) -> None:
    g, y, _ = _genotype_data()
    envelope = gwas.run({"genotypes": g, "phenotype": y}, options={"backend": backend})
    assert envelope["status"] == "not_implemented" and envelope["required_dependency"]
    with pytest.raises(ConfigurationError):
        gwas.run({"genotypes": g, "phenotype": y}, options={"backend": "nonsense"})


def test_gwas_simulate_is_deterministic_and_consistent() -> None:
    a = ok(
        gwas.simulate_gwas_data(
            None, options={"n_samples": 60, "n_variants": 12, "n_causal": 3, "seed": 7}
        )
    )
    b = ok(gwas.simulate_gwas_data({"n_samples": 60, "n_variants": 12, "n_causal": 3, "seed": 7}))
    assert a["phenotype"] == b["phenotype"] and a["true_effects"] == b["true_effects"]
    assert (
        len(a["genotypes"]["genotypes"]) == 60
        and len(a["genotypes"]["genotypes"][0]) == 12
        and len(a["true_effects"]) == 3
    )
    assert all(v in (0, 1, 2) for row in a["genotypes"]["genotypes"] for v in row)
    g0 = np.array(a["genotypes"]["genotypes"])[:, 0]
    ref = scipy_stats.linregress(g0, a["phenotype"])
    assert a["summary_statistics"][0]["beta"] == pytest.approx(ref.slope, rel=1e-9)
    c = ok(gwas.simulate_gwas_data(None, options={"n_samples": 60, "n_variants": 12, "seed": 8}))
    assert c["phenotype"] != a["phenotype"]
    with pytest.raises(ConfigurationError):
        gwas.simulate_gwas_data(None, options={"heritability": 1.5})


# ======================================================================= gwas: joint / conditional
def _correlated_scan(n: int = 400) -> tuple[np.ndarray, np.ndarray, list[dict], dict]:
    rng = np.random.default_rng(11)
    base = rng.normal(size=(n, 4))
    x = np.column_stack(
        [
            base[:, 0],
            0.6 * base[:, 0] + 0.8 * base[:, 1],
            base[:, 2],
            0.5 * base[:, 1] + 0.5 * base[:, 3],
        ]
    )
    y = 0.3 * x[:, 0] + 0.2 * x[:, 3] + rng.normal(size=n)
    table = []
    for j in range(4):
        fit = scipy_stats.linregress(x[:, j], y)
        table.append({"snp": f"v{j}", "z": fit.slope / fit.stderr, "n": n})
    r = np.corrcoef(x, rowvar=False)
    ld = {"variant_ids": [f"v{j}" for j in range(4)], "matrix": r.tolist()}
    return x, y, table, ld


def test_gwas_joint_and_conditional_match_statsmodels_ols() -> None:
    sm = pytest.importorskip("statsmodels.api")
    x, y, table, ld = _correlated_scan()
    xs = (x - x.mean(0)) / x.std(0)
    ys = (y - y.mean()) / y.std()
    joint = ok(gwas.joint_gwas_analysis(table, options={"ld": ld, "snps": ["v0", "v1", "v3"]}))
    fit = sm.OLS(ys, sm.add_constant(xs[:, [0, 1, 3]])).fit()
    got = {r["variant_id"]: r for r in joint["records"]}
    for k, v in enumerate(["v0", "v1", "v3"], start=1):
        assert got[v]["joint_beta_std"] == pytest.approx(fit.params[k], rel=1e-7)
        assert got[v]["joint_se_std"] == pytest.approx(fit.bse[k], rel=1e-7)
        assert got[v]["joint_p"] == pytest.approx(fit.pvalues[k], rel=1e-5)
    assert joint["model_r2"] == pytest.approx(fit.rsquared, rel=1e-7)
    cond = ok(gwas.conditional_gwas(table, options={"ld": ld, "condition_on": ["v0"]}))
    by = {r["variant_id"]: r for r in cond["records"]}
    assert set(by) == {"v1", "v2", "v3"}
    for j in (1, 2, 3):
        f = sm.OLS(ys, sm.add_constant(xs[:, [0, j]])).fit()
        assert by[f"v{j}"]["conditional_beta_std"] == pytest.approx(f.params[2], rel=1e-7)
        assert by[f"v{j}"]["conditional_p"] == pytest.approx(f.pvalues[2], rel=1e-5)
    # LD can also come from reference genotypes
    ref = {
        "genotypes": [[0.0, 1.0], [1.0, 1.0], [2.0, 0.0], [0.0, 2.0], [1.0, 1.0]],
        "variant_ids": ["v0", "v1"],
    }
    with pytest.raises(ConfigurationError):
        gwas.conditional_gwas(table, options={"ld": ld})
    with pytest.raises(ConfigurationError):
        gwas.conditional_gwas(table, options={"genotypes": ref, "condition_on": ["nope"]})


# ======================================================================================= gwas: PRS
def _prs_data(n: int = 400, seed: int = 5) -> tuple[list[float], list[float], list[int]]:
    rng = np.random.default_rng(seed)
    score = rng.normal(size=n)
    quant = 0.5 * score + rng.normal(size=n)
    binary = (rng.random(n) < 1 / (1 + np.exp(-(-1.0 + 0.8 * score)))).astype(int)
    return score.tolist(), quant.tolist(), binary.tolist()


def test_gwas_evaluate_prs_quantitative_and_binary() -> None:
    sm = pytest.importorskip("statsmodels.api")
    metrics = pytest.importorskip("sklearn.metrics")
    score, quant, binary = _prs_data()
    q = ok(gwas.evaluate_prs({"score": score, "phenotype": quant}))
    ref = scipy_stats.pearsonr(score, quant)
    assert (
        q["trait"] == "quantitative"
        and q["r"] == pytest.approx(ref[0], rel=1e-9)
        and q["r2"] == pytest.approx(ref[0] ** 2, rel=1e-9)
    )
    assert q["p"] == pytest.approx(ref[1], rel=1e-6)
    assert q["slope_per_sd"] == pytest.approx(
        scipy_stats.linregress(np.array(score) / np.std(score, ddof=1), quant).slope, rel=1e-9
    )
    b = ok(gwas.evaluate_prs({"score": score, "phenotype": binary}))
    assert b["auc"] == pytest.approx(metrics.roc_auc_score(binary, score), rel=1e-12)
    z = (np.array(score) - np.mean(score)) / np.std(score, ddof=1)
    fit = sm.Logit(binary, sm.add_constant(z)).fit(disp=0)
    assert b["log_or"] == pytest.approx(fit.params[1], rel=1e-6) and b["p"] == pytest.approx(
        fit.pvalues[1], rel=1e-5
    )
    assert b["or_ci95"] == pytest.approx(np.exp(fit.conf_int()[1]).tolist(), rel=1e-6)
    ll0 = fit.llnull
    n = len(binary)
    assert b["nagelkerke_r2"] == pytest.approx(
        (1 - math.exp(2 / n * (ll0 - fit.llf))) / (1 - math.exp(2 * ll0 / n)), rel=1e-6
    )
    with pytest.raises(InputFormatError):
        gwas.evaluate_prs({"score": score, "phenotype": quant}, options={"trait": "binary"})


def test_gwas_calibrate_prs_risk() -> None:
    sm = pytest.importorskip("statsmodels.api")
    score, _, binary = _prs_data()
    res = ok(
        gwas.calibrate_prs_risk(
            {"score": score, "phenotype": binary}, options={"n_bins": 5, "new_scores": [0.0]}
        )
    )
    z = (np.array(score) - np.mean(score)) / np.std(score, ddof=1)
    fit = sm.Logit(binary, sm.add_constant(z)).fit(disp=0)
    assert res["intercept"] == pytest.approx(fit.params[0], rel=1e-6) and res[
        "slope_per_sd"
    ] == pytest.approx(fit.params[1], rel=1e-6)
    assert res["predicted_risk"] == pytest.approx(
        fit.predict(sm.add_constant(z)).tolist(), rel=1e-6
    )
    assert res["mean_predicted"] == pytest.approx(
        np.mean(binary), abs=1e-6
    )  # logistic MLE: calibration-in-the-large
    assert len(res["calibration_table"]) == 5 and sum(
        r["n"] for r in res["calibration_table"]
    ) == len(binary)
    assert res["brier_score"] == pytest.approx(
        np.mean((fit.predict(sm.add_constant(z)) - binary) ** 2), rel=1e-6
    )
    mu = np.mean(score)
    expected_new = 1 / (
        1 + math.exp(-(fit.params[0] + fit.params[1] * (0.0 - mu) / np.std(score, ddof=1)))
    )
    assert res["new_score_risk"] == pytest.approx([expected_new], rel=1e-6)
    prev = np.mean(binary)
    shifted = ok(
        gwas.calibrate_prs_risk({"score": score, "phenotype": binary}, options={"prevalence": 0.05})
    )
    assert shifted["intercept"] == pytest.approx(
        fit.params[0] + math.log(0.05 / 0.95) - math.log(prev / (1 - prev)), rel=1e-6
    )
    with pytest.raises(InputFormatError):
        gwas.calibrate_prs_risk({"score": score, "phenotype": score})


def test_gwas_clump_and_threshold_prs() -> None:
    # v1/v2 perfectly correlated (r^2 = 1), v3 weakly related to v1, v4 far away
    geno = [[0, 0, 0, 2], [1, 1, 1, 0], [2, 2, 1, 1], [0, 0, 2, 2], [1, 1, 0, 0], [2, 2, 2, 1]]
    ss = [
        {"snp": "v1", "chr": 1, "bp": 100, "beta": 0.5, "p": 1e-6},
        {"snp": "v2", "chr": 1, "bp": 150, "beta": 0.4, "p": 1e-4},
        {"snp": "v3", "chr": 1, "bp": 200, "beta": -0.2, "p": 0.03},
        {"snp": "v4", "chr": 1, "bp": 10_000_000, "beta": 0.1, "p": 0.3},
    ]
    data = {
        "summary_statistics": ss,
        "genotypes": geno,
        "variant_ids": ["v1", "v2", "v3", "v4"],
        "sample_ids": list("abcdef"),
    }
    res = ok(
        gwas.clump_and_threshold_prs(
            data, options={"p_thresholds": [1e-5, 0.05, 0.5], "clump_r2": 0.5}
        )
    )
    r_13 = np.corrcoef(np.array(geno)[:, 0], np.array(geno)[:, 2])[0, 1]
    assert r_13**2 < 0.5
    assert [c["index"] for c in res["clumps"]] == ["v1", "v3", "v4"] and res["clumps"][0][
        "members"
    ] == ["v2"]
    by_t = {t["p_threshold"]: t for t in res["thresholds"]}
    assert (
        by_t[1e-5]["variants"] == ["v1"]
        and by_t[0.05]["variants"] == ["v1", "v3"]
        and by_t[0.5]["n_variants"] == 3
    )
    expect = {s: 0.5 * row[0] - 0.2 * row[2] for s, row in zip("abcdef", geno, strict=False)}
    assert by_t[0.05]["scores"] == pytest.approx(expect)
    expect_all = {s: expect[s] + 0.1 * row[3] for s, row in zip("abcdef", geno, strict=False)}
    assert by_t[0.5]["scores"] == pytest.approx(expect_all)
    with pytest.raises(InputFormatError):
        gwas.clump_and_threshold_prs(
            {"summary_statistics": ss, "genotypes": geno, "variant_ids": ["q1", "q2", "q3", "q4"]}
        )


# ============================================================================== gwas: heritability
def test_gwas_heritability_and_genetic_correlation_he() -> None:
    from dnakit.genomics.breeding import genomic_relationship_matrix

    rng = np.random.default_rng(21)
    n, m = 60, 150
    g = rng.binomial(2, rng.uniform(0.2, 0.5, m), size=(n, m))
    beta = rng.normal(size=m)
    gv = ((g - g.mean(0)) / np.maximum(g.std(0), 1e-9)) @ beta / math.sqrt(m)
    y1 = gv + rng.normal(size=n) * 0.8
    y2 = 0.7 * gv + rng.normal(size=n)
    grm = np.array(genomic_relationship_matrix(g.tolist(), verbose=False)["relationship_matrix"])
    iu = np.triu_indices(n, 1)

    def slope(a, b):
        return scipy_stats.linregress(a, b).slope

    z1, z2 = (y1 - y1.mean()) / y1.std(ddof=1), (y2 - y2.mean()) / y2.std(ddof=1)
    h2 = ok(gwas.estimate_snp_heritability({"genotypes": g.tolist(), "phenotype": y1.tolist()}))
    assert h2["method"] == "haseman_elston_cross_product"
    assert h2["h2"] == pytest.approx(slope(grm[iu], np.outer(z1, z1)[iu]), rel=1e-8)
    rg = ok(
        gwas.estimate_genetic_correlation(
            {"genotypes": g.tolist(), "phenotype_1": y1.tolist(), "phenotype_2": y2.tolist()}
        )
    )
    cross = (np.outer(z1, z2) + np.outer(z2, z1)) / 2
    h1, h2b, cov = (
        slope(grm[iu], np.outer(z1, z1)[iu]),
        slope(grm[iu], np.outer(z2, z2)[iu]),
        slope(grm[iu], cross[iu]),
    )
    assert rg["genetic_covariance_std"] == pytest.approx(cov, rel=1e-8)
    assert rg["rg"] == pytest.approx(cov / math.sqrt(h1 * h2b), rel=1e-8)


def test_gwas_ldsc_recovers_exact_linear_model() -> None:
    n_samples, m = 10_000, 200
    ell = [1.0 + 0.05 * i for i in range(m)]
    h2, h2_b = 0.4, 0.3
    # chi2 = 1 + N*h2/M*ell exactly -> z = sqrt(chi2); the second trait has a different h2
    rows1 = [
        {
            "snp": f"s{i}",
            "a1": "A",
            "a2": "G",
            "z": math.sqrt(1 + n_samples * h2 / m * ld),
            "n": n_samples,
            "ldscore": ld,
        }
        for i, ld in enumerate(ell)
    ]
    est = ok(gwas.estimate_snp_heritability(rows1))
    assert est["method"] == "ldsc_unweighted"
    # the linear model is exact for chi2, so the OLS slope reproduces h2 (intercept 1)
    assert est["h2"] == pytest.approx(h2, rel=1e-9) and est["intercept"] == pytest.approx(
        1.0, abs=1e-9
    )
    fixed = ok(gwas.estimate_snp_heritability(rows1, options={"fix_intercept": True}))
    assert fixed["h2"] == pytest.approx(h2, rel=1e-9)
    rows2 = [
        {
            "snp": r["snp"],
            "a1": "G",
            "a2": "A",
            "z": -math.sqrt(1 + n_samples * h2_b / m * r["ldscore"]),
            "n": n_samples,
        }
        for r in rows1
    ]
    # trait 2 is allele-swapped: its z flips sign back after alignment
    out = ok(
        gwas.estimate_genetic_correlation(
            {"trait1": rows1, "trait2": rows2},
            options={"ld_scores": {r["snp"]: r["ldscore"] for r in rows1}},
        )
    )
    assert out["h2_trait1"] == pytest.approx(h2, rel=1e-9) and out["h2_trait2"] == pytest.approx(
        h2_b, rel=1e-9
    )
    zz = [
        r1["z"] * (-r2["z"]) for r1, r2 in zip(rows1, rows2, strict=False)
    ]  # trait 2 is allele-swapped, so its aligned z is -z
    ml, mz = np.mean(ell), np.mean(zz)
    slope_x = sum((x - ml) * (v - mz) for x, v in zip(ell, zz, strict=False)) / sum(
        (x - ml) ** 2 for x in ell
    )
    assert out["genetic_covariance"] == pytest.approx(slope_x * m / n_samples, rel=1e-9)
    assert out["rg"] == pytest.approx(out["genetic_covariance"] / math.sqrt(h2 * h2_b), rel=1e-9)
    with pytest.raises(ConfigurationError):
        gwas.estimate_snp_heritability(rows1[:5])


# ==================================================================================== gwas: output
def test_gwas_write_export_and_overwrite_protection(tmp_path: Path) -> None:
    out = tmp_path / "out.tsv.gz"
    res = ok(gwas.write(TSV, options={"path": str(out), "columns": ["variant_id", "beta", "p"]}))
    lines = res["text"].splitlines()
    assert lines[0] == "variant_id\tbeta\tp" and lines[1].startswith("rs1\t0.5\t5.7e-07")
    assert float(lines[2].split("\t")[2]) == pytest.approx(
        2 * scipy_stats.norm.sf(2.0)
    )  # missing p is derived from beta/se
    assert gzip.decompress(out.read_bytes()).decode() == res["text"]
    with pytest.raises(ConfigurationError):
        gwas.write(TSV, options={"path": str(out)})
    assert ok(gwas.write(TSV, options={"path": str(out), "overwrite": True}))["path"] == str(out)
    ldsc = ok(
        gwas.export(
            [{"snp": "r1", "a1": "A", "a2": "G", "z": 2.0, "n": 100}], options={"format": "ldsc"}
        )
    )
    assert ldsc["text"] == "SNP\tA1\tA2\tN\tZ\nr1\tA\tG\t100\t2.0\n"
    metal = ok(gwas.export(TSV, options={"format": "metal"}))["text"].splitlines()
    assert (
        metal[0].split("\t")[:3] == ["MarkerName", "Allele1", "Allele2"]
        and metal[1].split("\t")[4] == "0.5"
    )
    bed = ok(gwas.export(TSV, options={"format": "bed"}))["text"].splitlines()
    assert bed[0].split("\t")[:3] == ["chr1", "99", "100"]
    assert len(ok(gwas.export(TSV, options={"format": "jsonl"}))["text"].splitlines()) == 2


def test_gwas_plot_svg() -> None:
    table = [
        {"snp": f"s{i}", "chr": 1 + i % 3, "bp": 100 + i, "p": 10 ** -(1 + i % 9)}
        for i in range(30)
    ]
    man = ok(gwas.plot(table))
    root = ET.fromstring(man["svg"])
    circles = [e for e in root.iter() if e.tag.endswith("circle")]
    assert len(circles) == 30 == man["n_points_drawn"]
    qq = ok(gwas.plot(table, options={"kind": "qq"}))
    assert len([e for e in ET.fromstring(qq["svg"]).iter() if e.tag.endswith("circle")]) == 30
    with pytest.raises(ConfigurationError):
        gwas.plot(table, options={"kind": "pie"})


# ====================================================================================== population
GENO = [[0, 1, 2, 1], [1, 1, 2, 0], [2, 0, 2, 1], [1, 2, 2, 1], [0, 1, 2, 0], [None, 1, 2, 2]]


def test_population_calculate_allele_frequency_and_heterozygosity() -> None:
    res = ok(
        population.calculate({"genotypes": GENO, "populations": ["A", "A", "A", "B", "B", "B"]})
    )
    first = res["overall"][0]
    assert (
        first["n_called"] == 5
        and first["alt_allele_frequency"] == pytest.approx(4 / 10)
        and first["maf"] == pytest.approx(0.4)
    )
    assert res["overall"][2]["alt_allele_frequency"] == 1.0 and res["overall"][2]["maf"] == 0.0
    assert res["by_population"]["A"][0]["alt_allele_frequency"] == pytest.approx(3 / 6)
    het = ok(population.calculate(GENO, options={"statistic": "heterozygosity"}))["variants"]
    assert het[1]["observed_heterozygosity"] == pytest.approx(4 / 6)
    p = 6 / 12  # alt allele frequency of the second variant (dosages 1,1,0,2,1,1)
    assert het[1]["expected_heterozygosity"] == pytest.approx(2 * p * (1 - p))
    assert het[1]["expected_heterozygosity_unbiased"] == pytest.approx(2 * p * (1 - p) * 12 / 11)
    assert het[1]["f_statistic"] == pytest.approx(1 - (4 / 6) / (2 * p * (1 - p)))
    with pytest.raises(InputFormatError):
        population.calculate([[0.5, 1.0], [1.0, 2.0]], options={"statistic": "heterozygosity"})


def test_population_ld_from_haplotype_pool() -> None:
    pool = [(1, 1), (1, 1), (0, 0), (1, 0)]  # haplotypes (allele at locus 1, locus 2)
    genotypes = [
        [a[0] + b[0], a[1] + b[1]] for a in pool for b in pool
    ]  # random mating: all ordered pairs
    pa, pb, pab = 3 / 4, 2 / 4, 2 / 4
    d = pab - pa * pb
    r2 = d * d / (pa * (1 - pa) * pb * (1 - pb))
    pair = ok(population.calculate(genotypes, options={"statistic": "ld"}))["pairs"][0]
    assert pair["d_composite"] == pytest.approx(d) and pair["r2"] == pytest.approx(r2)
    assert (
        pair["d_prime"]
        == pytest.approx(d / min(pa * (1 - pb), (1 - pa) * pb))
        == pytest.approx(1.0)
    )
    assert pair["r2"] == pytest.approx(np.corrcoef(np.array(genotypes).T)[0, 1] ** 2)
    with pytest.raises(ConfigurationError):
        population.calculate(genotypes, options={"statistic": "ld", "max_variants": 1})


def _fst_data(
    seed: int = 4, sizes: tuple[int, ...] = (30, 20, 25), m: int = 40
) -> tuple[list[list[int]], list[str]]:
    rng = np.random.default_rng(seed)
    rows, pops = [], []
    for k, n in enumerate(sizes):
        freqs = np.clip(rng.beta(2 + k, 3, size=m), 0.05, 0.95)
        rows += rng.binomial(2, freqs, size=(n, m)).tolist()
        pops += [f"P{k}"] * n
    return rows, pops


def _anova_theta(rows: list[list[int]], pops: list[str]) -> tuple[float, list[float]]:
    """Weir & Cockerham theta from the 1984 ANOVA (independent of the DNAKit implementation)."""

    arr = np.array(rows, dtype=float)
    labels = sorted(set(pops))
    idx = [np.where(np.array(pops) == lab)[0] for lab in labels]
    r = len(labels)
    n_i = np.array([len(i) for i in idx], dtype=float)
    n_tot = n_i.sum()
    nc = (n_tot - (n_i**2).sum() / n_tot) / (r - 1)
    num = den = 0.0
    per = []
    for j in range(arr.shape[1]):
        x = [arr[i, j] for i in idx]
        p_i = np.array([xi.mean() / 2 for xi in x])
        p_bar = (n_i * p_i).sum() / n_tot
        ssp = (2 * n_i * (p_i - p_bar) ** 2).sum()
        ssi = sum((xi**2).sum() / 2 - 2 * len(xi) * pi**2 for xi, pi in zip(x, p_i, strict=False))
        ssg = sum(xi.sum() - (xi**2).sum() / 2 for xi in x)
        msp, msi, msg = ssp / (r - 1), ssi / (n_tot - r), ssg / n_tot
        a, b, c = (msp - msi) / (2 * nc), (msi - msg) / 2, msg
        per.append(a / (a + b + c))
        num += a
        den += a + b + c
    return num / den, per


def test_population_compare_weir_cockerham_matches_anova_and_hand_case() -> None:
    rows, pops = _fst_data()
    res = ok(population.compare({"genotypes": rows, "populations": pops}))
    theta, per = _anova_theta(rows, pops)
    assert res["fst"] == pytest.approx(theta, rel=1e-9)
    assert list(res["per_variant"].values()) == pytest.approx(per, rel=1e-9)
    assert set(res["pairwise"]) == {"P0|P1", "P0|P2", "P1|P2"}
    keep = [i for i, p in enumerate(pops) if p in ("P0", "P2")]
    assert res["pairwise"]["P0|P2"] == pytest.approx(
        _anova_theta([rows[i] for i in keep], [pops[i] for i in keep])[0], rel=1e-9
    )
    fixed = ok(
        population.compare({"genotypes": [[0], [0], [2], [2]], "populations": ["a", "a", "b", "b"]})
    )
    assert fixed["fst"] == pytest.approx(1.0)  # populations fixed for opposite alleles


def test_population_compare_hudson_matches_formula() -> None:
    rows, pops = _fst_data(sizes=(30, 20))
    res = ok(
        population.compare({"genotypes": rows, "populations": pops}, options={"method": "hudson"})
    )
    arr = np.array(rows, dtype=float)
    a, b = arr[:30], arr[30:]
    p1, p2 = a.mean(0) / 2, b.mean(0) / 2
    h1, h2 = 60, 40
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (h1 - 1) - p2 * (1 - p2) / (h2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    assert res["fst"] == pytest.approx(num.sum() / den.sum(), rel=1e-10)
    assert list(res["per_variant"].values()) == pytest.approx((num / den).tolist(), rel=1e-10)
    with pytest.raises(InputFormatError):
        population.compare({"genotypes": rows})
    with pytest.raises(InputFormatError):
        population.compare({"genotypes": rows, "populations": ["x"] * len(rows)})


def test_population_estimate_grm_kinship_inbreeding_pca() -> None:
    from dnakit.genomics.breeding import genomic_relationship_matrix
    from dnakit.genomics.population import genotype_pca

    rng = np.random.default_rng(2)
    g = rng.binomial(2, 0.4, size=(8, 30)).tolist()
    ref = genomic_relationship_matrix(g, verbose=False)["relationship_matrix"]
    assert np.allclose(ok(population.estimate(g))["relationship_matrix"], ref, rtol=1e-12)
    assert np.allclose(
        ok(population.estimate(g, options={"quantity": "kinship"}))["relationship_matrix"],
        np.array(ref) / 2,
        rtol=1e-12,
    )
    pca = ok(population.estimate(g, options={"quantity": "pca", "components": 2}))
    assert np.allclose(pca["scores"], genotype_pca(g, components=2, verbose=False)["scores"])
    # inbreeding by hand: 4 samples, 3 loci; allele frequencies from the sample
    small = [[0, 1, 2], [0, 1, 2], [2, 1, 0], [2, 0, 0]]
    inbreeding = ok(population.estimate(small, options={"quantity": "inbreeding"}))["samples"]
    exp_hom = []
    for j in range(3):
        p = sum(r[j] for r in small) / 8
        exp_hom.append(1 - 2 * p * (1 - p) * 8 / 7)
    for i, rec in enumerate(inbreeding):
        obs = sum(1 for v in small[i] if v in (0, 2))
        assert rec["f"] == pytest.approx((obs - sum(exp_hom)) / (3 - sum(exp_hom)))


def _brute_force_hwe_p(n_aa: int, n_ab: int, n_bb: int) -> float:
    n = n_aa + n_ab + n_bb
    na, nb = 2 * n_aa + n_ab, 2 * n_bb + n_ab

    def prob(h: int) -> float:
        aa, bb = (na - h) // 2, (nb - h) // 2
        return math.exp(
            h * math.log(2)
            + math.lgamma(n + 1)
            + math.lgamma(na + 1)
            + math.lgamma(nb + 1)
            - math.lgamma(aa + 1)
            - math.lgamma(h + 1)
            - math.lgamma(bb + 1)
            - math.lgamma(2 * n + 1)
        )

    feasible = [h for h in range(min(na, nb) % 2, min(na, nb) + 1, 2)]
    p_obs = prob(n_ab)
    return sum(prob(h) for h in feasible if prob(h) <= p_obs * (1 + 1e-9))


@pytest.mark.parametrize(
    "counts", [(2, 2, 2), (10, 3, 17), (25, 0, 5), (1, 8, 1), (40, 11, 9), (0, 6, 0)]
)
def test_population_hwe_exact_and_chisq(counts: tuple[int, int, int]) -> None:
    rows = [[0]] * counts[0] + [[1]] * counts[1] + [[2]] * counts[2]
    exact = ok(population.test(rows))["records"][0]
    assert exact["p"] == pytest.approx(_brute_force_hwe_p(*counts), rel=1e-8)
    chi = ok(population.test(rows, options={"method": "chisq"}))["records"][0]
    n = sum(counts)
    p = (2 * counts[2] + counts[1]) / (2 * n)
    expected = [n * (1 - p) ** 2, 2 * n * p * (1 - p), n * p * p]
    if min(expected) > 0:
        ref = scipy_stats.chisquare(counts, expected, ddof=0)
        assert chi["chi2"] == pytest.approx(ref.statistic, rel=1e-9) and chi["p"] == pytest.approx(
            scipy_stats.chi2.sf(ref.statistic, 1), rel=1e-9
        )


def test_population_differentiation_test_and_multiple_testing() -> None:
    rows, pops = _fst_data(sizes=(30, 20, 25), m=6)
    res = ok(
        population.test(
            {"genotypes": rows, "populations": pops}, options={"test": "differentiation"}
        )
    )
    arr = np.array(rows)
    for j, rec in enumerate(res["records"]):
        table = [
            [
                2 * (np.array(pops) == lab).sum() - arr[np.array(pops) == lab, j].sum(),
                arr[np.array(pops) == lab, j].sum(),
            ]
            for lab in ("P0", "P1", "P2")
        ]
        ref = scipy_stats.chi2_contingency(table, correction=False)
        assert (
            rec["chi2"] == pytest.approx(ref[0], rel=1e-9)
            and rec["p"] == pytest.approx(ref[1], rel=1e-8)
            and rec["df"] == 2
        )
    sm = pytest.importorskip("statsmodels.stats.multitest")
    ps = [r["p"] for r in res["records"]]
    assert [r["q_bh"] for r in res["records"]] == pytest.approx(
        sm.multipletests(ps, method="fdr_bh")[1]
    )


def test_population_export_formats(tmp_path: Path) -> None:
    data = {
        "genotypes": [[0, 1], [2, None]],
        "sample_ids": ["s1", "s2"],
        "variant_ids": ["v1", "v2"],
        "populations": ["p", "q"],
        "ref": ["A", "C"],
        "alt": ["G", "T"],
        "chrom": ["1", "1"],
        "pos": [10, 20],
    }
    assert ok(population.export(data))["text"] == "sample_id\tv1\tv2\ns1\t0\t1\ns2\t2\tNA\n"
    pm = ok(population.export(data, options={"format": "ped_map", "path": str(tmp_path / "x")}))
    assert pm["ped"].splitlines() == ["p s1 0 0 0 -9 A A C T", "q s2 0 0 0 -9 G G 0 0"]
    assert pm["map"].splitlines() == ["1\tv1\t0\t10", "1\tv2\t0\t20"]
    assert (tmp_path / "x.ped").read_text() == pm["ped"] and (tmp_path / "x.map").exists()
    frq = ok(population.export(data, options={"format": "frq"}))["text"].splitlines()
    assert frq[1] == "v1\tG\tA\t0.5\t4" and frq[2] == "v2\tT\tC\t0.5\t2"
    with pytest.raises(ConfigurationError):
        population.export(data, options={"format": "ped_map", "path": str(tmp_path / "x")})


def test_population_plot_svg() -> None:
    rows, pops = _fst_data(sizes=(12, 12), m=25)
    hist = ok(population.plot({"genotypes": rows}, options={"bins": 5}))
    assert (
        sum(hist["bin_counts"]) == 25
        and len([e for e in ET.fromstring(hist["svg"]).iter() if e.tag.endswith("rect")]) >= 6
    )
    pca = ok(population.plot({"genotypes": rows, "populations": pops}, options={"kind": "pca"}))
    assert len([e for e in ET.fromstring(pca["svg"]).iter() if e.tag.endswith("circle")]) == 24
    fst = ok(population.plot({"genotypes": rows, "populations": pops}, options={"kind": "fst"}))
    assert fst["n_points"] == 25


def test_population_simulate_balding_nichols_and_wright_fisher() -> None:
    sim = ok(
        population.simulate(
            None,
            options={
                "n_populations": 3,
                "n_samples_per_population": 150,
                "n_variants": 300,
                "fst": 0.1,
                "seed": 5,
            },
        )
    )
    again = ok(
        population.simulate(
            {
                "n_populations": 3,
                "n_samples_per_population": 150,
                "n_variants": 300,
                "fst": 0.1,
                "seed": 5,
            }
        )
    )
    assert sim["genotypes"]["genotypes"] == again["genotypes"]["genotypes"]
    assert len(sim["genotypes"]["populations"]) == 450 and set(sim["genotypes"]["populations"]) == {
        "pop1",
        "pop2",
        "pop3",
    }
    est = ok(population.compare(sim["genotypes"]))["fst"]
    assert est == pytest.approx(
        0.1, abs=0.02
    )  # the simulated differentiation is recovered by Weir & Cockerham
    wf = ok(
        population.simulate(
            None,
            options={
                "model": "wright_fisher",
                "population_size": 50,
                "generations": 20,
                "initial_frequency": 0.3,
                "replicates": 1500,
                "seed": 1,
            },
        )
    )
    assert wf["mean_final_frequency"] == pytest.approx(0.3, abs=0.01)
    assert wf["variance_final_frequency"] == pytest.approx(
        0.3 * 0.7 * (1 - (1 - 1 / 100) ** 20), rel=0.1
    )
    assert wf["expected_variance"] == pytest.approx(0.3 * 0.7 * (1 - 0.99**20))
    with pytest.raises(ConfigurationError):
        population.simulate(None, options={"fst": 0.0})


# ======================================================================================== breeding
def test_breeding_breeders_equation_and_intensity() -> None:
    z = scipy_stats.norm.isf(0.1)
    intensity = scipy_stats.norm.pdf(z) / 0.1  # 1.7550 for the top 10 %
    res = ok(
        breeding.calculate(
            {
                "heritability": 0.4,
                "proportion_selected": 0.1,
                "phenotypic_sd": 10,
                "generation_interval": 2,
            }
        )
    )
    assert (
        res["selection_intensity"]
        == pytest.approx(intensity, rel=1e-9)
        == pytest.approx(1.755, abs=5e-4)
    )
    assert res["selection_differential"] == pytest.approx(intensity * 10, rel=1e-9)
    assert res["response"] == pytest.approx(0.4 * intensity * 10, rel=1e-9) and res[
        "response_per_year"
    ] == pytest.approx(res["response"] / 2)
    assert ok(breeding.calculate({"heritability": 0.25, "selection_differential": 8.0}))[
        "response"
    ] == pytest.approx(2.0)
    half = ok(
        breeding.calculate({"heritability": 0.5, "proportion_selected": 0.5, "phenotypic_sd": 1})
    )
    assert half["selection_intensity"] == pytest.approx(0.7978845608, rel=1e-8)  # 2 * phi(0)
    with pytest.raises(InputFormatError):
        breeding.calculate({"selection_differential": 2.0})
    with pytest.raises(ConfigurationError):
        breeding.calculate({"heritability": 1.5, "selection_differential": 2.0})


def test_breeding_heritability_estimates() -> None:
    parents = [1.0, 2.0, 3.0, 4.0, 5.0]
    offspring = [2.2, 2.4, 2.6, 2.8, 3.0]  # slope exactly 0.2
    mid = ok(
        breeding.calculate(
            {"parent": parents, "offspring": offspring}, options={"quantity": "heritability"}
        )
    )
    assert mid["h2"] == pytest.approx(0.2) and mid["n_pairs"] == 5
    single = ok(
        breeding.calculate(
            {"parent": parents, "offspring": offspring},
            options={"quantity": "heritability", "relationship": "single_parent"},
        )
    )
    assert single["h2"] == pytest.approx(0.4)
    assert ok(
        breeding.calculate(
            {"response": 2.0, "selection_differential": 5.0}, options={"quantity": "heritability"}
        )
    )["h2"] == pytest.approx(0.4)
    assert ok(
        breeding.calculate(
            {"genetic_variance": 3.0, "phenotypic_variance": 10.0},
            options={"quantity": "heritability"},
        )
    )["h2"] == pytest.approx(0.3)
    with pytest.raises(InputFormatError):
        breeding.calculate(
            {"genetic_variance": 12.0, "phenotypic_variance": 10.0},
            options={"quantity": "heritability"},
        )


def test_breeding_pedigree_relationship_textbook_cases() -> None:
    ped = [
        {"id": "s", "sire": 0, "dam": 0},
        {"id": "d", "sire": 0, "dam": 0},
        {"id": "x1", "sire": "s", "dam": "d"},
        {"id": "x2", "sire": "s", "dam": "d"},  # full sibs
        {"id": "h", "sire": "s", "dam": 0},  # half sib of x1
        {"id": "fs", "sire": "x1", "dam": "x2"},  # full-sib mating -> F = 1/4
        {"id": "po", "sire": "s", "dam": "x1"},  # father x daughter -> F = 1/4
        {"id": "hs", "sire": "x1", "dam": "h"},
    ]  # half-sib mating -> F = 1/8
    res = ok(breeding.calculate(ped, options={"quantity": "pedigree_relationship"}))
    ids = res["ids"]
    a = np.array(res["matrix"])
    rel = lambda i, j: a[ids.index(i), ids.index(j)]  # noqa: E731
    assert (
        rel("x1", "x2") == pytest.approx(0.5)
        and rel("x1", "s") == pytest.approx(0.5)
        and rel("h", "x1") == pytest.approx(0.25)
    )
    assert res["inbreeding"]["fs"] == pytest.approx(0.25) and res["inbreeding"][
        "po"
    ] == pytest.approx(0.25)
    assert res["inbreeding"]["hs"] == pytest.approx(0.125) and res["inbreeding"]["x1"] == 0.0
    assert rel("fs", "fs") == pytest.approx(1.25)
    # the relationship matrix is symmetric and independent of the record order
    shuffled = ok(
        breeding.calculate(list(reversed(ped)), options={"quantity": "pedigree_relationship"})
    )
    assert shuffled["inbreeding"] == res["inbreeding"]
    with pytest.raises(InputFormatError):
        breeding.calculate(
            [{"id": "a", "sire": "b", "dam": 0}, {"id": "b", "sire": "a", "dam": 0}],
            options={"quantity": "pedigree_relationship"},
        )
    with pytest.raises(InputFormatError):
        breeding.calculate(
            [{"id": "a", "sire": "ghost", "dam": 0}], options={"quantity": "pedigree_relationship"}
        )


def _breeding_data(n: int = 30, m: int = 80, seed: int = 8) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, rng.uniform(0.2, 0.6, m), size=(n, m)).astype(float)
    y = (g - g.mean(0)) @ rng.normal(scale=0.3, size=m) + rng.normal(size=n)
    return g, y


def test_breeding_fit_equals_gblup_and_predict_variants() -> None:
    g, y = _breeding_data()
    n = len(y)
    h2 = 0.4
    p = g.mean(0) / 2
    z = g - g.mean(0)
    k = z @ z.T
    s = 2 * np.sum(p * (1 - p))
    lam_g = (1 - h2) / h2
    gblup = y.mean() + (k / s) @ np.linalg.solve(
        k / s + lam_g * np.eye(n), y - y.mean()
    )  # GBLUP with G = ZZ'/(2 sum pq)
    fitted = ok(
        breeding.fit({"genotypes": g.tolist(), "phenotypes": y.tolist()}, options={"h2": h2})
    )
    assert fitted["ridge_lambda"] == pytest.approx(s * lam_g, rel=1e-12)
    assert np.allclose(fitted["training_predictions"], gblup, rtol=1e-9, atol=1e-9)
    # predictions for new individuals: y_bar + G_new,train (G + lambda I)^-1 (y - y_bar)
    new = np.random.default_rng(1).binomial(2, 0.4, size=(5, g.shape[1])).astype(float)
    znew = new - g.mean(0)
    expected = y.mean() + (znew @ z.T / s) @ np.linalg.solve(
        k / s + lam_g * np.eye(n), y - y.mean()
    )
    via_model = ok(breeding.predict({"model": fitted, "genotypes": new.tolist()}))
    refit = ok(
        breeding.predict(
            {"train_genotypes": g.tolist(), "phenotypes": y.tolist(), "genotypes": new.tolist()},
            options={"h2": h2},
        )
    )
    assert np.allclose(via_model["predictions"], expected, rtol=1e-8, atol=1e-8)
    assert np.allclose(refit["predictions"], expected, rtol=1e-8, atol=1e-8)
    assert np.allclose(via_model["breeding_values"], np.array(expected) - y.mean())
    explicit = ok(
        breeding.fit(
            {"genotypes": g.tolist(), "phenotypes": y.tolist()}, options={"ridge_lambda": 5.0}
        )
    )
    assert explicit["ridge_lambda"] == 5.0 and explicit["ridge_lambda_source"] == "ridge_lambda"
    with pytest.raises(ConfigurationError):
        breeding.fit({"genotypes": g.tolist(), "phenotypes": y.tolist()})
    with pytest.raises(InputFormatError):
        breeding.predict({"model": {"marker_effects": []}, "genotypes": new.tolist()})


def test_breeding_select_truncation_and_smith_hazel_index() -> None:
    cands = [{"id": f"c{i}", "ebv": v} for i, v in enumerate([1.0, 5.0, 3.0, 4.0, 2.0])]
    res = ok(breeding.select(cands, options={"n_select": 2}))
    assert res["selected"] == ["c1", "c3"] and res["selection_differential"] == pytest.approx(
        4.5 - 3.0
    )
    assert ok(breeding.select(cands, options={"proportion": 0.4, "higher_is_better": False}))[
        "selected"
    ] == ["c0", "c4"]
    g = [[4.0, 1.0], [1.0, 2.0]]
    pmat = [[10.0, 3.0], [3.0, 6.0]]
    a = [2.0, 1.0]
    table = [
        {"id": f"i{k}", "t1": x, "t2": y}
        for k, (x, y) in enumerate([(5, 1), (2, 6), (4, 4), (0, 0)])
    ]
    index = ok(
        breeding.select(
            table,
            options={
                "method": "index",
                "traits": ["t1", "t2"],
                "economic_weights": a,
                "genetic_covariance": g,
                "phenotypic_covariance": pmat,
                "n_select": 1,
            },
        )
    )
    b = np.linalg.solve(pmat, np.array(g) @ a)
    assert index["index_weights"] == pytest.approx(b.tolist())
    scores = {r["id"]: r["score"] for r in index["ranking"]}
    assert scores["i2"] == pytest.approx(b[0] * 4 + b[1] * 4)
    assert index["accuracy"] == pytest.approx(
        math.sqrt(b @ np.array(pmat) @ b / (np.array(a) @ np.array(g) @ np.array(a)))
    )
    assert index["selected"] == [max(scores, key=scores.get)]
    with pytest.raises(ConfigurationError):
        breeding.select(cands)
    with pytest.raises(ConfigurationError):
        breeding.select(cands, options={"n_select": 9})


def test_breeding_optimize_matches_quadratic_programming() -> None:
    from scipy.optimize import minimize

    rng = np.random.default_rng(3)
    n = 8
    z = rng.normal(size=(n, 20))
    a = z @ z.T / 20 + 0.1 * np.eye(n)
    a = a / np.sqrt(np.outer(np.diag(a), np.diag(a)))
    u = rng.normal(size=n)
    cands = [{"id": f"c{i}", "ebv": u[i]} for i in range(n)]
    rel = {"ids": [f"c{i}" for i in range(n)], "matrix": a.tolist()}
    lam = 2.0
    res = ok(
        breeding.optimize({"candidates": cands, "relationship": rel}, options={"penalty": lam})
    )
    c = np.array([res["contributions"].get(f"c{i}", 0.0) for i in range(n)])
    ref = minimize(
        lambda x: -(x @ u - lam / 2 * x @ a @ x),
        np.full(n, 1 / n),
        method="SLSQP",
        bounds=[(0, 1)] * n,
        constraints=[{"type": "eq", "fun": lambda x: x.sum() - 1}],
        options={"ftol": 1e-14, "maxiter": 500},
    )
    assert c.sum() == pytest.approx(1.0) and (c >= -1e-12).all()
    assert c == pytest.approx(ref.x, abs=1e-5)
    assert res["objective"] == pytest.approx(-ref.fun, abs=1e-8) and res[
        "mean_ebv"
    ] == pytest.approx(c @ u)
    # with A = I only the best candidates are retained: c_i = max(0, u_i - nu) / lam
    ident = {"ids": rel["ids"], "matrix": np.eye(n).tolist()}
    flat = ok(
        breeding.optimize({"candidates": cands, "relationship": ident}, options={"penalty": 1.0})
    )
    kept = flat["contributions"]
    nu = (sum(u[int(k[1:])] for k in kept) - 1.0) / len(kept)
    assert all(kept[k] == pytest.approx(u[int(k[1:])] - nu) for k in kept)
    from_pedigree = ok(
        breeding.optimize(
            {
                "candidates": [
                    {"id": "a", "ebv": 1.0},
                    {"id": "b", "ebv": 1.0},
                    {"id": "c", "ebv": 0.9},
                ],
                "pedigree": [
                    {"id": "p", "sire": 0, "dam": 0},
                    {"id": "a", "sire": "p", "dam": 0},
                    {"id": "b", "sire": "p", "dam": 0},
                    {"id": "c", "sire": 0, "dam": 0},
                    {"id": "p2", "sire": 0, "dam": 0},
                ],
            },
            options={"penalty": 2.0},
        )
    )
    assert sum(from_pedigree["contributions"].values()) == pytest.approx(1.0)
    with pytest.raises(InputFormatError):
        breeding.optimize({"candidates": cands})


def test_breeding_simulate_is_deterministic_and_follows_the_breeders_equation() -> None:
    params = {
        "n_individuals": 3000,
        "n_loci": 120,
        "heritability": 0.5,
        "selection_proportion": 0.2,
        "generations": 3,
        "seed": 4,
    }
    a = ok(breeding.simulate(None, options=params))
    b = ok(breeding.simulate(params))
    assert a["generations"] == b["generations"] and len(a["generations"]) == 4
    first = a["generations"][0]
    assert first["realised_heritability_in_generation"] == pytest.approx(0.5, abs=0.06)
    # selection differential of the top 20 % of a normal trait: i * sigma_p, i = phi(z)/p
    i20 = scipy_stats.norm.pdf(scipy_stats.norm.isf(0.2)) / 0.2
    assert first["selection_differential"] == pytest.approx(
        i20 * math.sqrt(first["phenotypic_variance"]), rel=0.05
    )
    assert a["generations"][1]["observed_response"] == pytest.approx(
        first["predicted_response"], rel=0.2
    )
    assert a["generations"][3]["mean_genetic_value"] > a["generations"][0]["mean_genetic_value"]
    assert (
        a["generations"][3]["genetic_variance"] < a["generations"][0]["genetic_variance"]
    )  # selection reduces variance
    other = ok(breeding.simulate(None, options={**params, "seed": 5}))
    assert other["generations"][0]["mean_genetic_value"] != first["mean_genetic_value"]
    with pytest.raises(ConfigurationError):
        breeding.simulate(None, options={"heritability": 0.0})


def test_breeding_plot_svg() -> None:
    sim = ok(
        breeding.simulate(
            None, options={"n_individuals": 50, "n_loci": 10, "generations": 3, "seed": 1}
        )
    )
    svg = ok(breeding.plot(sim))["svg"]
    poly = [e for e in ET.fromstring(svg).iter() if e.tag.endswith("polyline")]
    assert len(poly) == 1 and len(poly[0].attrib["points"].split()) == 4
    two = ok(breeding.plot({"series": {"a": [1, 2, 3], "b": [3, 2, 1]}}))
    assert two["series"] == {"a": 3, "b": 3}
    with pytest.raises(InputFormatError):
        breeding.plot({"nothing": 1})


# =========================================================================== forensic_reproductive
FREQ = {
    "L1": {"12": 0.2, "13": 0.3, "14": 0.5},
    "L2": {"8": 0.1, "9": 0.9},
    "L3": {"5": 0.25, "6": 0.25, "7": 0.5},
}


def test_forensic_random_match_probability_product_rule_and_theta() -> None:
    profile = {"L1": [12, 13], "L2": [9, 9], "L3": ["7", "7"]}
    res = ok(forensic.calculate({"profile": profile, "frequencies": FREQ}))
    loci = {r["locus"]: r["genotype_probability"] for r in res["per_locus"]}
    assert loci == pytest.approx({"L1": 2 * 0.2 * 0.3, "L2": 0.9**2, "L3": 0.25})
    assert res["combined"] == pytest.approx(0.12 * 0.81 * 0.25) and res["one_in"] == pytest.approx(
        1 / (0.12 * 0.81 * 0.25)
    )
    th = 0.03
    corrected = ok(
        forensic.calculate({"profile": profile, "frequencies": FREQ}, options={"theta": th})
    )
    d = (1 + th) * (1 + 2 * th)
    expect = {
        "L1": 2 * (th + (1 - th) * 0.2) * (th + (1 - th) * 0.3) / d,
        "L2": (2 * th + (1 - th) * 0.9) * (3 * th + (1 - th) * 0.9) / d,
        "L3": (2 * th + (1 - th) * 0.5) * (3 * th + (1 - th) * 0.5) / d,
    }
    assert {r["locus"]: r["genotype_probability"] for r in corrected["per_locus"]} == pytest.approx(
        expect
    )
    assert (
        corrected["combined"] > res["combined"]
    )  # theta raises the match probability of homozygotes
    # frequency table as records, alleles as numbers or strings
    rows = [
        {"locus": "L1", "allele": a, "frequency": f}
        for a, f in {"12": 0.2, "13": 0.3, "14": 0.5}.items()
    ]
    assert ok(forensic.calculate({"profile": {"L1": [12.0, 13]}, "frequencies": rows}))[
        "combined"
    ] == pytest.approx(0.12)
    with pytest.raises(ConfigurationError):
        forensic.calculate({"profile": {"L1": [12, 99]}, "frequencies": FREQ})
    floor = ok(
        forensic.calculate(
            {"profile": {"L1": [12, 99]}, "frequencies": FREQ}, options={"min_frequency": 0.01}
        )
    )
    assert floor["combined"] == pytest.approx(2 * 0.2 * 0.01)
    with pytest.raises(InputFormatError):
        forensic.calculate({"profile": {"L1": [12, 13, 14]}, "frequencies": FREQ})


def test_forensic_likelihood_ratio_single_source() -> None:
    evidence = {"L1": [12, 13], "L2": [9, 9]}
    res = ok(
        forensic.calculate(
            {"evidence": evidence, "suspect": {"L1": [13, 12], "L2": [9, 9]}, "frequencies": FREQ},
            options={"quantity": "likelihood_ratio"},
        )
    )
    assert res["combined"] == pytest.approx(1 / (0.12 * 0.81)) and not res["excluded"]
    assert [r["lr"] for r in res["per_locus"]] == pytest.approx([1 / 0.12, 1 / 0.81])
    assert res["verbal_equivalent"].startswith("moderate support")  # LR ~ 10.3 on the ENFSI scale
    excl = ok(
        forensic.calculate(
            {"evidence": evidence, "suspect": {"L1": [12, 14], "L2": [9, 9]}, "frequencies": FREQ},
            options={"quantity": "likelihood_ratio"},
        )
    )
    assert excl["excluded"] and excl["combined"] == 0.0 and excl["per_locus"][0]["match"] is False
    with pytest.raises(InputFormatError):
        forensic.calculate(
            {"evidence": {"A": [1, 2]}, "suspect": {"B": [1, 2]}, "frequencies": FREQ},
            options={"quantity": "likelihood_ratio"},
        )


def test_forensic_paternity_index_textbook_cases() -> None:
    f = {"L": {"a": 0.2, "b": 0.3, "c": 0.5}}

    def pi(child: list[str], father: list[str], mother: list[str] | None = None) -> dict:
        data: dict = {"child": {"L": child}, "alleged_father": {"L": father}, "frequencies": f}
        if mother:
            data["mother"] = {"L": mother}
        return ok(forensic.calculate(data, options={"quantity": "paternity_index"}))

    assert pi(["a", "b"], ["b", "b"], ["a", "a"])["combined"] == pytest.approx(
        1 / 0.3
    )  # obligate allele b, father homozygous: 1/p
    assert pi(["a", "b"], ["b", "c"], ["a", "a"])["combined"] == pytest.approx(
        0.5 / 0.3
    )  # father heterozygous: 1/(2p)
    assert pi(["a", "b"], ["a", "b"], ["a", "b"])["combined"] == pytest.approx(
        1 / (0.2 + 0.3)
    )  # ambiguous: 1/(pa+pb)
    assert pi(["a", "a"], ["a", "c"], ["a", "a"])["combined"] == pytest.approx(0.5 / 0.2)
    exclusion = pi(["a", "b"], ["c", "c"], ["a", "a"])
    assert (
        exclusion["combined"] == 0.0
        and exclusion["exclusions"] == 1
        and exclusion["probability_of_paternity"] == 0.0
    )
    duo = pi(
        ["a", "b"], ["b", "c"]
    )  # no mother: P(child | father, random mother) / HW genotype probability
    assert duo["design"] == "duo" and duo["combined"] == pytest.approx(
        (0.5 * 0.2) / (2 * 0.2 * 0.3)
    )
    with pytest.raises(InputFormatError):  # child not consistent with the mother
        pi(["a", "b"], ["b", "b"], ["c", "c"])
    # combined index over loci and probability of paternity
    data = {
        "child": {"L1": [12, 13], "L2": [9, 9]},
        "alleged_father": {"L1": [12, 14], "L2": [8, 9]},
        "mother": {"L1": [13, 14], "L2": [9, 9]},
        "frequencies": FREQ,
    }
    res = ok(forensic.calculate(data, options={"quantity": "paternity_index", "prior": 0.5}))
    assert [r["pi"] for r in res["per_locus"]] == pytest.approx([0.5 / 0.2, 0.5 / 0.9])
    cpi = (0.5 / 0.2) * (0.5 / 0.9)
    assert res["combined"] == pytest.approx(cpi) and res[
        "probability_of_paternity"
    ] == pytest.approx(cpi / (cpi + 1))
    with pytest.raises(ConfigurationError):
        forensic.calculate(data, options={"quantity": "paternity_index", "theta": 0.01})


def test_forensic_kinship_index_matches_published_formulas() -> None:
    f = {"L": {"a": 0.2, "b": 0.3, "c": 0.5}}
    pa, pb = 0.2, 0.3

    def lr(g1: list[str], g2: list[str], rel: str) -> float:
        return ok(
            forensic.calculate(
                {"profile_1": {"L": g1}, "profile_2": {"L": g2}, "frequencies": f},
                options={"quantity": "kinship_index", "relationship": rel},
            )
        )["combined"]

    assert lr(["a", "a"], ["a", "a"], "full_sibling") == pytest.approx(
        0.25 + 0.5 / pa + 0.25 / pa**2
    )
    # full siblings, both ab: k0 + k1 (pa+pb)/(4 pa pb) + k2 / (2 pa pb)
    assert lr(["a", "b"], ["a", "b"], "full_sibling") == pytest.approx(
        0.25 + 0.5 * (pa + pb) / (4 * pa * pb) + 0.25 / (2 * pa * pb)
    )
    assert lr(["a", "a"], ["a", "a"], "parent_child") == pytest.approx(
        1 / pa
    )  # k1 only: p(a)/p(a)^2
    # parent ab, child ac: shared allele a (prob 1/2), other allele c from the population (0.5)
    assert lr(["a", "b"], ["a", "c"], "parent_child") == pytest.approx(0.5 * 0.5 / (2 * 0.2 * 0.5))
    # parent aa, child ab: shared allele a (prob 1), other allele b from the population
    assert lr(["a", "a"], ["a", "b"], "parent_child") == pytest.approx(pb / (2 * pa * pb))
    assert (
        lr(["a", "b"], ["c", "c"], "parent_child") == 0.0
    )  # no shared allele: impossible without mutation
    assert lr(["a", "b"], ["a", "b"], "identical_twin") == pytest.approx(1 / (2 * pa * pb))
    assert lr(["a", "a"], ["a", "a"], "unrelated") == pytest.approx(1.0)
    custom = ok(
        forensic.calculate(
            {"profile_1": {"L": ["a", "a"]}, "profile_2": {"L": ["a", "a"]}, "frequencies": f},
            options={"quantity": "kinship_index", "ibd": [0.5, 0.5, 0.0]},
        )
    )
    assert (
        custom["combined"] == pytest.approx(0.5 + 0.5 / pa) and custom["relationship"] == "custom"
    )
    with pytest.raises(ConfigurationError):
        forensic.calculate(
            {"profile_1": {"L": ["a", "a"]}, "profile_2": {"L": ["a", "a"]}, "frequencies": f},
            options={"quantity": "kinship_index", "ibd": [0.5, 0.2, 0.2]},
        )


def test_forensic_compare_infer_read_validate() -> None:
    a = {"L1": [12, 13], "L2": [9, 9], "L3": [5, 6]}
    b = {"L1": [13, 12], "L2": [8, 9], "L3": [7, 7], "L4": [1, 1]}
    cmp = ok(forensic.compare({"profile_1": a, "profile_2": b}))
    assert [r["status"] for r in cmp["per_locus"]] == [
        "identical",
        "one_allele_shared",
        "no_allele_shared",
    ]
    assert cmp["identical"] == 1 and cmp["one_allele_shared"] == 1 and cmp["no_allele_shared"] == 1
    assert (
        cmp["loci_only_in_profile_2"] == ["L4"]
        and cmp["consistent_with_same_source"] is False
        and cmp["mean_alleles_shared"] == pytest.approx(1.0)
    )
    withfreq = ok(forensic.compare({"profile_1": a, "profile_2": a, "frequencies": FREQ}))
    assert withfreq["consistent_with_same_source"] and withfreq[
        "random_match_probability_profile_1"
    ] == pytest.approx(0.12 * 0.81 * 2 * 0.25 * 0.25)

    f = {"L": {"a": 0.2, "b": 0.3, "c": 0.5}}
    inferred = ok(
        forensic.infer(
            {"profile_1": {"L": ["a", "b"]}, "profile_2": {"L": ["a", "b"]}, "frequencies": f},
            options={"relationships": ["unrelated", "full_sibling", "parent_child"]},
        )
    )
    lrs = {r["relationship"]: r["likelihood_ratio_vs_unrelated"] for r in inferred["relationships"]}
    assert lrs["unrelated"] == pytest.approx(1.0)
    norm = sum(lrs.values())
    assert {r["relationship"]: r["posterior"] for r in inferred["relationships"]} == pytest.approx(
        {k: v / norm for k, v in lrs.items()}
    )
    assert inferred["most_probable"] == max(lrs, key=lrs.get)
    biased = ok(
        forensic.infer(
            {"profile_1": {"L": ["a", "b"]}, "profile_2": {"L": ["c", "c"]}, "frequencies": f},
            options={
                "relationships": ["unrelated", "parent_child"],
                "priors": {"unrelated": 0.9, "parent_child": 0.1},
            },
        )
    )
    assert biased["most_probable"] == "unrelated" and biased["relationships"][0][
        "posterior"
    ] == pytest.approx(1.0)

    prof = ok(forensic.read("locus,allele1,allele2\nL1,12,13\nL2,9,\n"))
    assert prof["kind"] == "profile" and prof["profile"] == {"L1": ["12", "13"], "L2": ["9", "9"]}
    fr = ok(forensic.read("Locus\tAllele\tFrequency\nL1\t12\t0.4\nL1\t13\t0.6\n"))
    assert fr["kind"] == "frequencies" and fr["frequencies"] == {"L1": {"12": 0.4, "13": 0.6}}
    good = ok(forensic.validate({"frequencies": FREQ, "profile": {"L1": [12, 13]}}))
    assert good["valid"] and good["issues"] == []
    bad = ok(
        forensic.validate(
            {
                "frequencies": {"L1": {"12": 0.5, "13": 0.2}},
                "profile": {"L1": [12, 14], "L2": [1, 2, 3]},
                "child": {"L9": [1, 2]},
            }
        )
    )
    messages = " | ".join(i["message"] for i in bad["issues"])
    assert (
        not bad["valid"]
        and "sum to 0.7000" in messages
        and "not in the frequency table" in messages
        and "3 alleles" in messages
    )
    assert "locus missing" in messages
    one = ok(forensic.validate({"profile": {"L1": [12]}}))
    assert one["valid"] and one["n_warnings"] == 1


def test_forensic_report_and_plot(tmp_path: Path) -> None:
    data = {
        "evidence": {"L1": [12, 13], "L2": [9, 9]},
        "suspect": {"L1": [12, 13], "L2": [9, 9]},
        "frequencies": FREQ,
    }
    rep = ok(
        forensic.report(
            data,
            options={
                "quantity": "likelihood_ratio",
                "title": "Case 1",
                "path": str(tmp_path / "r.md"),
            },
        )
    )
    text = rep["text"]
    assert text.startswith("# Case 1") and "| L1 | 8.333 |" in text and "| L2 | 1.235 |" in text
    assert (
        "Combined: **10.29**" in text
        and "ENFSI 2015" in text
        and (tmp_path / "r.md").read_text() == text
    )
    assert rep["statistic"]["combined"] == pytest.approx(1 / (0.12 * 0.81))
    svg = ok(forensic.plot(data, options={"quantity": "likelihood_ratio"}))["svg"]
    assert (
        len([e for e in ET.fromstring(svg).iter() if e.tag.endswith("rect")]) == 4
    )  # background + plot frame + 2 bars
    prof = ok(forensic.plot({"profile": {"L1": [12, 13], "L2": [9, 9]}, "frequencies": FREQ}))
    assert prof["n_loci"] == 2


# ============================================================================= experimental_design
def _nct_power(
    d: float, n1: int, n2: int | None, alpha: float = 0.05, two_sided: bool = True
) -> float:
    if n2 is None:
        df, delta = n1 - 1, d * math.sqrt(n1)
    else:
        df, delta = n1 + n2 - 2, d * math.sqrt(n1 * n2 / (n1 + n2))
    if two_sided:
        tc = scipy_stats.t.ppf(1 - alpha / 2, df)
        return float(scipy_stats.nct.sf(tc, df, delta) + scipy_stats.nct.cdf(-tc, df, delta))
    return float(scipy_stats.nct.sf(scipy_stats.t.ppf(1 - alpha, df), df, delta))


@pytest.mark.parametrize(
    "test, d, alt, two_group",
    [
        ("two_sample_t", 0.5, "two-sided", True),
        ("two_sample_t", 0.8, "larger", True),
        ("one_sample_t", 0.4, "two-sided", False),
        ("paired_t", 0.3, "larger", False),
    ],
)
def test_design_t_test_sample_size_is_the_smallest_n_reaching_the_power(
    test: str, d: float, alt: str, two_group: bool
) -> None:
    res = ok(design.calculate({"effect_size": d, "test": test, "alternative": alt, "power": 0.8}))
    n = res["n"]
    two_sided = alt == "two-sided"
    assert _nct_power(d, n, n if two_group else None, two_sided=two_sided) >= 0.8
    assert _nct_power(d, n - 1, n - 1 if two_group else None, two_sided=two_sided) < 0.8
    assert res["achieved_power"] == pytest.approx(
        _nct_power(d, n, n if two_group else None, two_sided=two_sided), abs=1e-6
    )
    if test == "two_sample_t" and alt == "two-sided":
        assert (
            n == 64 and res["total"] == 128
        )  # the textbook value for d = 0.5, 80 % power, alpha = 0.05
    power = ok(
        design.calculate(
            {"effect_size": d, "test": test, "alternative": alt, "n": 25, "quantity": "power"}
        )
    )["power"]
    assert power == pytest.approx(
        _nct_power(d, 25, 25 if two_group else None, two_sided=two_sided), abs=1e-6
    )


def test_design_unequal_allocation_t_test() -> None:
    res = ok(design.calculate({"effect_size": 0.5, "ratio": 2.0}))
    n1, n2 = res["n"], res["n_group2"]
    assert n2 == 2 * n1 and _nct_power(0.5, n1, n2) >= 0.8 > _nct_power(0.5, n1 - 1, 2 * (n1 - 1))


def test_design_two_proportions_and_chi_square() -> None:
    p1, p2, alpha, power = 0.5, 0.75, 0.05, 0.9
    za, zb = scipy_stats.norm.isf(alpha / 2), scipy_stats.norm.ppf(power)
    exact = (
        za * math.sqrt((p1 + p2) * (1 - (p1 + p2) / 2))
        + zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
    ) ** 2 / (p1 - p2) ** 2
    res = ok(design.calculate({"test": "two_proportions", "p1": p1, "p2": p2, "power": power}))
    assert res["n_unrounded"] == pytest.approx(exact, rel=1e-12) == pytest.approx(76.7, abs=0.05)
    assert res["n_per_group"] == math.ceil(exact) == 77
    back = ok(
        design.calculate(
            {"test": "two_proportions", "p1": p1, "p2": p2, "n": 77, "quantity": "power"}
        )
    )["power"]
    manual = scipy_stats.norm.cdf(
        (abs(p1 - p2) * math.sqrt(77) - za * math.sqrt((p1 + p2) * (1 - (p1 + p2) / 2)))
        / math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
    )
    assert back == pytest.approx(manual, rel=1e-12) and back >= 0.9
    chi = ok(design.calculate({"test": "chi_square", "effect_size": 0.3, "df": 2, "power": 0.8}))
    crit = scipy_stats.chi2.ppf(0.95, 2)
    n = chi["n_total"]
    assert (
        scipy_stats.ncx2.sf(crit, 2, n * 0.09) >= 0.8 > scipy_stats.ncx2.sf(crit, 2, (n - 1) * 0.09)
    )
    assert n == 108  # Cohen (1988): w = 0.3, df = 2, power 0.80 -> N = 108
    assert ok(
        design.calculate(
            {"test": "chi_square", "effect_size": 0.3, "df": 2, "n": 100, "quantity": "power"}
        )
    )["power"] == pytest.approx(scipy_stats.ncx2.sf(crit, 2, 100 * 0.09), abs=1e-9)
    with pytest.raises(InputFormatError):
        design.calculate({"test": "two_proportions", "p1": 0.5, "p2": 0.5})
    with pytest.raises(ConfigurationError):
        design.calculate({"effect_size": 0.5, "alpha": 1.5})
    with pytest.raises(InputFormatError):
        design.calculate({"test": "two_sample_t"})


def test_design_lander_waterman_and_replicates() -> None:
    c = 150 * 1_000_000 / 5_000_000
    res = ok(
        design.calculate(
            {
                "quantity": "coverage",
                "genome_size": 5e6,
                "read_length": 150,
                "n_reads": 1e6,
                "target_fraction": 0.99,
            }
        )
    )
    assert res["coverage"] == pytest.approx(c) and res["fraction_covered"] == pytest.approx(
        1 - math.exp(-c)
    )
    assert res["expected_uncovered_bases"] == pytest.approx(5e6 * math.exp(-c)) and res[
        "expected_contigs"
    ] == pytest.approx(1e6 * math.exp(-c))
    assert res["coverage_needed_for_target_fraction"] == pytest.approx(math.log(100)) and res[
        "reads_needed_for_target_fraction"
    ] == math.ceil(math.log(100) * 5e6 / 150)
    five = ok(design.calculate({"quantity": "coverage", "coverage": 5.0, "min_depth": 3}))
    assert five["fraction_covered"] == pytest.approx(1 - math.exp(-5))
    assert five["fraction_covered_min_depth"] == pytest.approx(
        scipy_stats.poisson.sf(2, 5.0), rel=1e-12
    )
    assert ok(
        design.calculate(
            {"quantity": "coverage", "genome_size": 1000, "read_length": 100, "n_reads": 20}
        )
    )["expected_contig_length"] == pytest.approx(100 * ((math.exp(2.0) - 1) / 2.0 + 1 - 1))
    rep = ok(design.calculate({"quantity": "replicates", "sd": 10, "margin": 5}))
    n = rep["replicates"]
    assert (
        scipy_stats.t.ppf(0.975, n - 1) * 10 / math.sqrt(n)
        <= 5
        < scipy_stats.t.ppf(0.975, n - 2) * 10 / math.sqrt(n - 1)
    )
    assert (
        rep["replicates_normal_approximation"]
        == math.ceil((scipy_stats.norm.isf(0.025) * 2) ** 2)
        == 16
        and n == 18
    )
    with pytest.raises(InputFormatError):
        design.calculate({"quantity": "coverage"})


def test_design_layouts_crd_rcbd_latin_factorial() -> None:
    crd = ok(design.design({"treatments": ["A", "B", "C"], "replicates": 4, "seed": 3}))
    assert crd["n_plots"] == 12 and crd["replications"] == {"A": 4, "B": 4, "C": 4}
    assert (
        crd["layout"]
        == ok(design.design({"treatments": ["A", "B", "C"], "replicates": 4, "seed": 3}))["layout"]
    )
    assert (
        crd["layout"]
        != ok(design.design({"treatments": ["A", "B", "C"], "replicates": 4, "seed": 4}))["layout"]
    )
    rcbd = ok(
        design.design(
            {"treatments": ["A", "B", "C", "D"], "replicates": 3, "design_type": "rcbd", "seed": 1}
        )
    )
    for block in (1, 2, 3):
        assert sorted(r["treatment"] for r in rcbd["layout"] if r["block"] == block) == [
            "A",
            "B",
            "C",
            "D",
        ]
    latin = ok(
        design.design({"treatments": list("ABCD"), "design_type": "latin_square", "seed": 9})
    )
    grid = {(r["row"], r["column"]): r["treatment"] for r in latin["layout"]}
    for k in range(1, 5):
        assert sorted(grid[(k, c)] for c in range(1, 5)) == list("ABCD") and sorted(
            grid[(r, k)] for r in range(1, 5)
        ) == list("ABCD")
    fact = ok(
        design.design(
            {
                "design_type": "factorial",
                "factors": {"N": [0, 50], "P": ["lo", "hi", "mid"]},
                "replicates": 2,
                "seed": 2,
            }
        )
    )
    assert (
        fact["n_plots"] == 12
        and len(fact["replications"]) == 6
        and {r["N"] for r in fact["layout"]} == {"0", "50"}
    )
    blocked = ok(design.design({"design_type": "factorial", "factors": {"A": [0, 1]}, "blocks": 3}))
    assert blocked["n_plots"] == 6 and {r["block"] for r in blocked["layout"]} == {1, 2, 3}
    with pytest.raises(InputFormatError):
        design.design({"treatments": ["A", "A"]})
    with pytest.raises(InputFormatError):
        design.design({})


def test_design_generate_factorial_matrices() -> None:
    full = ok(design.generate({"factors": {"A": [-1, 1], "B": [-1, 1], "C": ["x", "y"]}}))
    assert [(r["A"], r["B"], r["C"]) for r in full["runs"]] == [
        ("-1", "-1", "x"),
        ("1", "-1", "x"),
        ("-1", "1", "x"),
        ("1", "1", "x"),
        ("-1", "-1", "y"),
        ("1", "-1", "y"),
        ("-1", "1", "y"),
        ("1", "1", "y"),
    ]
    frac = ok(
        design.generate({"kind": "fractional_factorial", "n_factors": 4, "generators": ["D=ABC"]})
    )
    assert frac["n_runs"] == 8 and all(r["D"] == r["A"] * r["B"] * r["C"] for r in frac["runs"])
    assert [r["A"] for r in frac["runs"]] == [-1, 1, -1, 1, -1, 1, -1, 1]
    checked = ok(
        design.check(
            frac["runs"],
            options={"treatment_column": "run", "factor_columns": ["A", "B", "C", "D"]},
        )
    )
    assert checked["orthogonal"] is True and all(
        v == 0 for v in checked["column_dot_products"].values()
    )
    cyc = ok(design.generate({"kind": "latin_square", "treatments": list("ABC")}))["square"]
    assert cyc == [["A", "B", "C"], ["B", "C", "A"], ["C", "A", "B"]]
    with pytest.raises(InputFormatError):
        design.generate({"kind": "fractional_factorial", "n_factors": 4, "generators": ["D=ABD"]})
    assert (
        ok(design.generate({"kind": "fractional_factorial", "n_factors": 4}))["n_runs"] == 16
    )  # no generator = full 2^4 design
    with pytest.raises(InputFormatError):
        design.generate({"kind": "fractional_factorial"})


def test_design_randomize_and_assign() -> None:
    items = list(range(20))
    res = ok(design.randomize(items, options={"seed": 11}))
    assert sorted(r["item"] for r in res["runs"]) == items and [
        r["run_order"] for r in res["runs"]
    ] == list(range(1, 21))
    assert [r["item"] for r in res["runs"]] == [
        r["item"] for r in ok(design.randomize(items, options={"seed": 11}))["runs"]
    ]
    assert [r["item"] for r in res["runs"]] != items
    rows = [{"id": i, "block": "b" + str(i // 5)} for i in range(15)]
    within = ok(design.randomize(rows, options={"seed": 1, "within": "block"}))["runs"]
    assert [r["block"] for r in within] == ["b0"] * 5 + ["b1"] * 5 + ["b2"] * 5
    units = [{"unit": f"u{i}", "sex": "M" if i % 2 else "F"} for i in range(20)]
    asg = ok(
        design.assign(units, options={"groups": ["ctl", "trt"], "stratify_by": "sex", "seed": 2})
    )
    for sex in ("M", "F"):
        counts = [
            sum(1 for a in asg["assignments"] if a["sex"] == sex and a["group"] == g)
            for g in ("ctl", "trt")
        ]
        assert counts == [5, 5]
    ratio = ok(design.assign(30, options={"groups": ["a", "b"], "ratio": [1, 2], "seed": 5}))
    assert ratio["group_sizes"] == {"a": 10, "b": 20} and ratio["n_units"] == 30
    assert ok(design.assign(["x", "y", "z"], options={"n_groups": 3}))["group_sizes"] == {
        "G1": 1,
        "G2": 1,
        "G3": 1,
    }
    with pytest.raises(ConfigurationError):
        design.assign(10, options={"groups": ["a"]})
    with pytest.raises(ConfigurationError):
        design.assign(10, options={"ratio": [1, 0]})
    with pytest.raises(ConfigurationError):
        design.randomize(rows, options={"within": "nonexistent"})


def test_design_check_validate_export(tmp_path: Path) -> None:
    rcbd = ok(
        design.design(
            {"treatments": list("ABC"), "design_type": "rcbd", "replicates": 2, "seed": 0}
        )
    )
    chk = ok(design.check(rcbd))
    assert chk["ok"] and chk["balanced"] and chk["complete_blocks"] and chk["n_blocks"] == 2
    broken = [r for r in rcbd["layout"] if not (r["block"] == 2 and r["treatment"] == "C")]
    bad = ok(design.check(broken))
    assert not bad["ok"] and bad["complete_blocks"] is False and not bad["balanced"]
    latin = ok(design.design({"treatments": list("ABC"), "design_type": "latin_square", "seed": 1}))
    assert ok(design.check(latin))["latin_square"] is True
    swapped = [dict(r) for r in latin["layout"]]
    swapped[0]["treatment"] = swapped[1]["treatment"]
    assert ok(design.check(swapped))["latin_square"] is False
    assert ok(design.validate(rcbd))["valid"]
    invalid = ok(design.validate([{"plot": 1, "treatment": "A"}, {"plot": 1, "treatment": ""}]))
    assert (
        not invalid["valid"]
        and any("empty" in e for e in invalid["errors"])
        and any("duplicated" in e for e in invalid["errors"])
    )
    assert not ok(design.validate([{"plot": 1}]))["valid"]
    exp = ok(design.export(rcbd, options={"path": str(tmp_path / "d.csv")}))
    assert exp["text"].splitlines()[0] == "plot,block,treatment" and exp["n_rows"] == 6
    assert (tmp_path / "d.csv").read_text() == exp["text"]
    assert (
        ok(design.export(rcbd, options={"format": "json", "columns": ["plot"]}))["text"].count(
            '"plot"'
        )
        == 6
    )
    with pytest.raises(ConfigurationError):
        design.export(rcbd, options={"path": str(tmp_path / "d.csv")})
    with pytest.raises(InputFormatError):
        design.check([{"x": 1}])


def test_design_simulated_power_agrees_with_exact_power() -> None:
    res = ok(design.simulate({"n": 20, "effect_size": 0.8, "n_simulations": 4000, "seed": 3}))
    assert res["analytic_power"] == pytest.approx(_nct_power(0.8, 20, 20), abs=1e-6)
    assert abs(res["empirical_power"] - res["analytic_power"]) < 4 * res["standard_error"] + 0.005
    again = ok(design.simulate({"n": 20, "effect_size": 0.8, "n_simulations": 4000, "seed": 3}))
    assert again["empirical_power"] == res["empirical_power"]
    one = ok(
        design.simulate(
            {"n": 15, "effect_size": 0.7, "test": "one_sample_t", "n_simulations": 3000, "seed": 1}
        )
    )
    assert (
        abs(one["empirical_power"] - _nct_power(0.7, 15, None)) < 4 * one["standard_error"] + 0.005
    )
    null = ok(design.simulate({"n": 20, "effect_size": 0.0, "n_simulations": 3000, "seed": 2}))
    assert null["empirical_power"] == pytest.approx(0.05, abs=0.02)  # type I error rate
    with pytest.raises(InputFormatError):
        design.simulate({"n": 20})


# ======================================================== leftovers, optional numpy, invalid input
@pytest.mark.parametrize(
    "module, name", [(population, "infer"), (forensic, "deconvolute"), (forensic, "predict")]
)
def test_leftover_functions_are_honestly_not_implemented(module: object, name: str) -> None:
    env = getattr(module, name)({"genotypes": GENO, "profile": {"L1": [1, 2]}})
    assert (
        env["status"] == "not_implemented"
        and env["required_dependency"]
        and env["message"]
        and "result" not in env
    )


def test_numpy_dependent_functions_report_not_implemented_without_numpy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import dnakit.genomics.breeding as genomics_breeding
    from dnakit.exceptions import BackendUnavailableError

    def missing(module: str, extra: str) -> None:
        raise BackendUnavailableError(f"{module} could not be imported", hint="install it")

    monkeypatch.setattr(genomics_breeding, "require", missing)
    g, y = _breeding_data()
    fit_env = breeding.fit({"genotypes": g.tolist(), "phenotypes": y.tolist()}, options={"h2": 0.5})
    assert fit_env["status"] == "not_implemented" and "numpy" in fit_env["required_dependency"]
    grm_env = population.estimate(g.tolist())
    assert grm_env["status"] == "not_implemented" and "numpy" in grm_env["required_dependency"]


def test_pure_python_grm_matches_numpy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    from dnakit.gwas import _impl_models

    g = np.random.default_rng(6).binomial(2, 0.35, size=(12, 40)).astype(float)
    rows = g.tolist()
    rows[2][3] = None
    with_numpy = _impl_models.vanraden_grm(rows)
    monkeypatch.setattr(_impl_models, "_numpy", lambda: None)
    pure = _impl_models.vanraden_grm(rows)
    assert (
        np.allclose(with_numpy["matrix"], pure["matrix"], rtol=1e-12, atol=1e-12)
        and with_numpy["n_markers"] == pure["n_markers"]
    )
    from dnakit.genomics.breeding import genomic_relationship_matrix

    assert np.allclose(
        pure["matrix"],
        genomic_relationship_matrix(rows, verbose=False)["relationship_matrix"],
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "call",
    [
        lambda: gwas.read(12345),
        lambda: gwas.standardize([]),
        lambda: gwas.standardize([{"foo": 1, "bar": 2}]),
        lambda: gwas.standardize([{"snp": "a", "p": "not-a-number"}]),
        lambda: gwas.calculate([{"snp": "a", "eaf": 0.3}], options={"statistic": "lambda_gc"}),
        lambda: gwas.run({"genotypes": [[0, 1], [1]], "phenotype": [1, 2]}),
        lambda: gwas.run({"genotypes": [[0, 3]] * 6, "phenotype": list(range(6))}),
        lambda: gwas.run("not a mapping"),
        lambda: gwas.evaluate_prs({"score": [1, 2], "phenotype": [1, 2]}),
        lambda: gwas.calculate([{"snp": "a", "p": 0.5}], options={"statistic": "bogus"}),
        lambda: gwas.calculate([{"snp": "a", "p": 0.5}], options="not a mapping"),
        lambda: population.calculate({"genotypes": [[0, 5], [1, 1]]}),
        lambda: population.calculate({"nothing": 1}),
        lambda: population.test([[0, 1], [1, 1]], options={"method": "bogus"}),
        lambda: breeding.select([], options={"n_select": 1}),
        lambda: forensic.calculate("text", options={}),
        lambda: forensic.calculate({"profile": {"L": [1, 2]}}),
        lambda: design.calculate({"effect_size": -1.0}),
        lambda: design.randomize([]),
    ],
)
def test_invalid_input_raises_dnakit_errors(call: object) -> None:
    with pytest.raises(DNAKitError) as info:
        call()  # type: ignore[operator]
    assert info.value.code and info.value.hint
