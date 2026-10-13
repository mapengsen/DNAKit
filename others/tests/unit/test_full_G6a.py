"""Hand-checkable and upstream-equivalence tests for the G6a variant wrappers.

Covers varcode (offline splice, SV, genotype and phasing helpers), scikit-allel
(population statistics on local arrays) and hgvs (indel normalization, CIGAR and repeat
arithmetic). Every test runs offline on inline or tmp_path inputs.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.variants import allel_functions as allel_api
from dnakit.variants import hgvs_functions as hgvs_api
from dnakit.variants import varcode_functions as varcode_api

np = pytest.importorskip("numpy")

_VARCODE_MODULE = "varcode"
_ALLEL_MODULE = "allel"
_HGVS_MODULE = "hgvs"


def _need(name: str):
    """Return the importable package or submodule, skipping the test when it is missing."""
    return pytest.importorskip(name)


# ---------------------------------------------------------------- module contracts


@pytest.mark.parametrize("module", [varcode_api, allel_api, hgvs_api])
def test_public_functions_take_keyword_only_verbose(module) -> None:
    for name in module.__all__:
        function = getattr(module, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_wrapped_counts_match_the_brief_cap() -> None:
    assert len(varcode_api.__all__) == 18
    assert len(allel_api.__all__) == 20
    assert len(hgvs_api.__all__) == 6
    for module in (varcode_api, allel_api, hgvs_api):
        assert len(set(module.__all__)) == len(module.__all__)


# ---------------------------------------------------------------- varcode


def test_enumerate_cryptic_exon_candidates_matches_hand_built_sequence() -> None:
    _need("varcode.cryptic_exons")
    # acceptor CAGG at local index 4, exon AAG ahead of donor GTAAGT at local index 41.
    sequence = "TTTT" + "CAGG" + "C" * 30 + "AAG" + "GTAAGT" + "TTTT"
    result = varcode_api.enumerate_cryptic_exon_candidates(
        "chr1", sequence, 1000, min_intron_length=10, verbose=False
    )
    assert result["count"] == 1
    candidate = result["candidates"][0]
    assert (candidate["interval_start"], candidate["interval_end"]) == (1007, 1041)
    assert candidate["donor_score"] == 1.0 and candidate["acceptor_score"] == 1.0
    with pytest.raises(ConfigurationError):
        varcode_api.enumerate_cryptic_exon_candidates("chr1", "ACGX", verbose=False)


def test_score_splice_site_donor_and_acceptor_anchors() -> None:
    varcode = _need("varcode.cryptic_exons") and _need("varcode")
    donor = varcode_api.score_splice_site("CAGGTAAGT", site="donor", verbose=False)
    assert donor["score"] == 1.0
    assert varcode_api.score_splice_site("CAGATAAGT", verbose=False)["score"] == 0.0
    assert donor["score"] == pytest.approx(varcode.cryptic_exons.score_donor("CAGGTAAGT"))
    acceptor = varcode_api.score_splice_site("CAGG", site="acceptor", verbose=False)
    assert acceptor["score"] == 1.0
    assert varcode_api.score_splice_site("CAAG", site="acceptor", verbose=False)["score"] == 0.0
    with pytest.raises(ConfigurationError):
        varcode_api.score_splice_site("CAGG", site="donor", verbose=False)


def test_parse_symbolic_sv_alt_reads_deletion_and_breakend() -> None:
    _need("varcode.sv_allele_parser")
    deletion = varcode_api.parse_symbolic_sv_alt(
        "chr1", 100, "N", "<DEL>", info={"END": 500}, verbose=False
    )
    assert deletion["symbolic"] is True
    assert (deletion["sv_type"], deletion["start"], deletion["end"]) == ("DEL", 100, 500)
    assert (deletion["affected_start"], deletion["affected_end"]) == (101, 500)
    breakend = varcode_api.parse_symbolic_sv_alt(
        "chr1", 100, "A", "A]chr1:500]", info={"SVTYPE": "DEL", "MATEID": "b2"}, verbose=False
    )
    assert breakend["sv_type"] == "BND"
    assert (breakend["mate_contig"], breakend["mate_start"]) == ("chr1", 500)
    plain = varcode_api.parse_symbolic_sv_alt("chr1", 100, "A", "G", verbose=False)
    assert plain["symbolic"] is False
    with pytest.raises(ConfigurationError):
        varcode_api.parse_symbolic_sv_alt(
            "chr1", 100, "N", "<DEL>", info={"END": 50}, verbose=False
        )


def test_breakend_helpers_follow_vcf_bracket_rules() -> None:
    _need(_VARCODE_MODULE)
    sides = varcode_api.breakend_alt_sides("A]chr1:500]", verbose=False)
    assert sides["valid"] is True
    assert (sides["this_side"], sides["mate_side"]) == ("left", "left")
    assert varcode_api.breakend_alt_sides("A", verbose=False)["valid"] is False
    mate = varcode_api.breakend_mate_position("A]chr1:500]", verbose=False)
    assert (mate["mate_contig"], mate["mate_position"], mate["mate_orientation"]) == (
        "chr1",
        500,
        "]]",
    )
    inserted = varcode_api.breakend_insert_bases("A", "AGT[chr1:500[", verbose=False)
    assert inserted["known"] is True and inserted["inserted_bases"] == "GT"
    unknown = varcode_api.breakend_insert_bases("AC", "AGT[chr1:500[", verbose=False)
    assert unknown["known"] is False


def test_classify_breakend_pair_types_a_deletion_and_rejects_plain_alleles() -> None:
    varcode = _need("varcode.structural_variant") and _need("varcode")
    first = {
        "contig": "chr1",
        "start": 100,
        "ref": "A",
        "alt": "A]chr1:500]",
        "info": {"SVTYPE": "DEL", "MATEID": "b2"},
    }
    second = {
        "contig": "chr1",
        "start": 500,
        "ref": "T",
        "alt": "[chr1:100[T",
        "info": {"SVTYPE": "DEL", "MATEID": "b1"},
    }
    typed = varcode_api.classify_breakend_pair(first, second, verbose=False)
    assert typed["typed"] is True
    assert (typed["sv_type"], typed["start"], typed["end"]) == ("DEL", 100, 499)
    first_event = varcode.sv_allele_parser.parse_symbolic_alt(
        "chr1", 100, "A", "A]chr1:500]", info=first["info"]
    )
    second_event = varcode.sv_allele_parser.parse_symbolic_alt(
        "chr1", 500, "T", "[chr1:100[T", info=second["info"]
    )
    assert varcode.structural_variant.typed_event_from_breakends(first_event, second_event) == (
        "DEL",
        100,
        499,
    )
    with pytest.raises(ConfigurationError):
        varcode_api.classify_breakend_pair(
            {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}, second, verbose=False
        )


def test_compare_sv_calls_groups_identical_calls_and_matches_upstream() -> None:
    varcode = _need("varcode.sv_comparison") and _need("varcode")
    rows = [
        {
            "call_id": "a",
            "caller": "x",
            "sample": "s",
            "build": "GRCh38",
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "<DEL>",
            "info": {"END": 500},
        },
        {
            "call_id": "b",
            "caller": "y",
            "sample": "s",
            "build": "GRCh38",
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "<DEL>",
            "info": "END=500",
        },
        {
            "call_id": "c",
            "caller": "x",
            "sample": "s",
            "build": "GRCh38",
            "chrom": "chr1",
            "pos": 5000,
            "ref": "A",
            "alt": "<DEL>",
            "info": {"END": 5400},
        },
    ]
    result = varcode_api.compare_sv_calls(rows, max_distance=100, verbose=False)
    groups = result["result"]["groups"]
    assert result["call_count"] == 3
    assert sorted(len(group["call_ids"]) for group in groups) == [1, 2]
    pair = next(group for group in groups if len(group["call_ids"]) == 2)
    assert pair["call_ids"] == ["a", "b"] and pair["callers"] == ["x", "y"]
    text_rows = [
        {**row, "info": "END=5400" if row["call_id"] == "c" else "END=500"} for row in rows
    ]
    direct = varcode.sv_comparison.compare_sv_calls(text_rows, max_distance=100)
    assert [g["call_ids"] for g in direct["groups"]] == [g["call_ids"] for g in groups]


def test_trim_shared_allele_flanks_hand_case_and_deletion() -> None:
    _need(_VARCODE_MODULE)
    result = varcode_api.trim_shared_allele_flanks("ACGT", "ACTT", verbose=False)
    assert (result["trimmed_ref"], result["trimmed_alt"], result["prefix"], result["suffix"]) == (
        "G",
        "T",
        "AC",
        "T",
    )
    deletion = varcode_api.trim_shared_allele_flanks("AAA", "AA", verbose=False)
    assert (deletion["trimmed_ref"], deletion["trimmed_alt"], deletion["prefix"]) == ("A", "", "AA")


def test_parse_vcf_genotype_handles_phase_and_missing_calls() -> None:
    _need(_VARCODE_MODULE)
    phased = varcode_api.parse_vcf_genotype("0|1", verbose=False)
    assert (phased["alleles"], phased["phased"], phased["ploidy"]) == ([0, 1], True, 2)
    missing = varcode_api.parse_vcf_genotype("./.", verbose=False)
    assert (missing["alleles"], missing["phased"]) == ([None, None], False)
    half = varcode_api.parse_vcf_genotype("1|.", verbose=False)
    assert half["alleles"] == [1, None]
    with pytest.raises(ConfigurationError):
        varcode_api.parse_vcf_genotype("x", verbose=False)


def test_normalize_nucleotide_allele_maps_gaps_and_uppercases() -> None:
    _need(_VARCODE_MODULE)
    assert varcode_api.normalize_nucleotide_allele(".", verbose=False)["is_empty"] is True
    assert varcode_api.normalize_nucleotide_allele("-", verbose=False)["normalized"] == ""
    assert varcode_api.normalize_nucleotide_allele("acgt", verbose=False)["normalized"] == "ACGT"


def test_variant_overlaps_interval_edges_and_insertions() -> None:
    _need(_VARCODE_MODULE)
    inside = varcode_api.variant_overlaps_interval(10, 3, 12, 20, verbose=False)
    assert inside["overlaps"] is True and inside["variant_end"] == 12
    # A 3-base variant at 10 covers 10..12, so an interval starting at 13 does not overlap.
    outside = varcode_api.variant_overlaps_interval(10, 3, 13, 20, verbose=False)
    assert outside["overlaps"] is False
    snv_neighbour = varcode_api.variant_overlaps_interval(10, 1, 11, 20, verbose=False)
    assert snv_neighbour["overlaps"] is False
    snv_hit = varcode_api.variant_overlaps_interval(10, 1, 10, 10, verbose=False)
    assert snv_hit["overlaps"] is True
    insertion = varcode_api.variant_overlaps_interval(10, 0, 10, 10, verbose=False)
    assert insertion["overlaps"] is True and insertion["variant_end"] is None


def test_four_gamete_phase_test_cis_trans_and_unknown() -> None:
    _need(_VARCODE_MODULE)
    cis = varcode_api.four_gamete_phase_test(10, 1, 2, 10, verbose=False)
    assert (cis["relation"], cis["in_cis"]) == ("cis", True)
    trans = varcode_api.four_gamete_phase_test(1, 10, 10, 1, verbose=False)
    assert (trans["relation"], trans["in_cis"]) == ("trans", False)
    unknown = varcode_api.four_gamete_phase_test(0, 0, 0, 0, verbose=False)
    assert (unknown["relation"], unknown["in_cis"]) == ("unknown", None)
    with pytest.raises(ConfigurationError):
        varcode_api.four_gamete_phase_test(-1, 1, 1, 1, verbose=False)


SOMATIC = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
GERMLINE = {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}


def test_detect_germline_allele_overlap_reports_matching_index() -> None:
    _need(_VARCODE_MODULE)
    result = varcode_api.detect_germline_allele_overlap(
        SOMATIC, [{"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}, GERMLINE], verbose=False
    )
    assert result["overlap"] is True and result["matched_indices"] == [0]
    assert result["germline_count"] == 2


def test_partition_germline_by_phase_uses_evidence_and_homozygous_calls() -> None:
    _need(_VARCODE_MODULE)
    unphased = varcode_api.partition_germline_by_phase(SOMATIC, [GERMLINE], verbose=False)
    assert unphased["cis"] == [] and unphased["trans"] == [] and unphased["unphased"] == [GERMLINE]
    homozygous = varcode_api.partition_germline_by_phase(
        SOMATIC, [GERMLINE], homozygous=[GERMLINE], verbose=False
    )
    assert homozygous["cis"] == [GERMLINE]
    cis = varcode_api.partition_germline_by_phase(
        SOMATIC,
        [GERMLINE],
        in_cis_evidence=[{"first": SOMATIC, "second": GERMLINE, "in_cis": True}],
        verbose=False,
    )
    assert cis["cis"] == [GERMLINE] and cis["unphased"] == []
    trans = varcode_api.partition_germline_by_phase(
        SOMATIC,
        [GERMLINE],
        in_cis_evidence=[{"first": SOMATIC, "second": GERMLINE, "in_cis": False}],
        verbose=False,
    )
    assert trans["trans"] == [GERMLINE] and trans["cis"] == []


def test_enumerate_phase_hypotheses_reports_phased_configuration() -> None:
    _need(_VARCODE_MODULE)
    result = varcode_api.enumerate_phase_hypotheses(
        SOMATIC,
        [GERMLINE],
        in_cis_evidence=[{"first": SOMATIC, "second": GERMLINE, "in_cis": True}],
        verbose=False,
    )
    assert result["count"] == 1
    hypothesis = result["hypotheses"][0]
    assert hypothesis["phase_state"] == "phased"
    assert hypothesis["cis"] == [GERMLINE] and hypothesis["trans"] == []
    with pytest.raises(ConfigurationError):
        varcode_api.enumerate_phase_hypotheses(
            SOMATIC,
            [GERMLINE, {"contig": "chr1", "start": 300, "ref": "G", "alt": "A"}],
            max_hypotheses=1,
            verbose=False,
        )


def test_load_maf_variants_reads_local_file(tmp_path: Path) -> None:
    _need(_VARCODE_MODULE)
    columns = [
        "Hugo_Symbol",
        "Entrez_Gene_Id",
        "Center",
        "NCBI_Build",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Strand",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "dbSNP_RS",
        "dbSNP_Val_Status",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "Match_Norm_Seq_Allele1",
        "Match_Norm_Seq_Allele2",
    ]
    rows = [
        {
            "Hugo_Symbol": "TP53",
            "NCBI_Build": "GRCh38",
            "Chromosome": "17",
            "Start_Position": "7676154",
            "End_Position": "7676154",
            "Strand": "+",
            "Variant_Classification": "Missense_Mutation",
            "Variant_Type": "SNP",
            "Reference_Allele": "C",
            "Tumor_Seq_Allele1": "C",
            "Tumor_Seq_Allele2": "T",
            "Tumor_Sample_Barcode": "TCGA-AA-0001-01A",
        },
        {
            "Hugo_Symbol": "KRAS",
            "NCBI_Build": "GRCh38",
            "Chromosome": "12",
            "Start_Position": "25245350",
            "End_Position": "25245350",
            "Strand": "+",
            "Variant_Classification": "Missense_Mutation",
            "Variant_Type": "SNP",
            "Reference_Allele": "C",
            "Tumor_Seq_Allele1": "C",
            "Tumor_Seq_Allele2": "A",
            "Tumor_Sample_Barcode": "TCGA-AA-0001-01A",
        },
    ]
    path = tmp_path / "sample.maf"
    path.write_text(
        "\t".join(columns)
        + "\n"
        + "".join("\t".join(row.get(column, "") for column in columns) + "\n" for row in rows),
        encoding="utf-8",
    )
    result = varcode_api.load_maf_variants(str(path), verbose=False)
    assert result["reference"] == "GRCh38"
    assert result["count"] == 2
    assert result["variants"][0] == {
        "contig": "12",
        "start": 25245350,
        "ref": "C",
        "alt": "A",
        "reference_name": "GRCh38",
    }
    with pytest.raises(ConfigurationError):
        varcode_api.load_maf_variants(str(tmp_path / "missing.maf"), verbose=False)


def test_check_vcf_sample_identity_matches_upstream_report(tmp_path: Path) -> None:
    _need("varcode.sample_identity")
    import varcode

    path = tmp_path / "pair.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/1:30\n"
        "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP\t1/1:30\t1/1:30\n",
        encoding="utf-8",
    )
    result = varcode_api.check_vcf_sample_identity(str(path), verbose=False)
    assert result["sample_count"] == 2
    pair = result["report"]["pairs"][0]
    assert pair["germline"]["shared_snps"] == 2
    assert pair["germline"]["exact_matches"] == 2
    samples = varcode.sample_identity.load_vcf_samples(path)
    direct = varcode.sample_identity.check_sample_identity(samples)
    assert direct["pairs"][0]["germline"] == pair["germline"]
    with pytest.raises(ConfigurationError):
        varcode_api.check_vcf_sample_identity(str(path), kind="tumour", verbose=False)
    with pytest.raises(ConfigurationError):
        varcode_api.check_vcf_sample_identity(str(path), config={"no_such_field": 1}, verbose=False)


# ---------------------------------------------------------------- scikit-allel

AC = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
POS = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
AC1 = [[4, 0], [3, 1], [2, 2], [0, 4]]
AC2 = [[0, 4], [1, 3], [2, 2], [4, 0]]


def test_windowed_tajima_d_matches_upstream_windows() -> None:
    allel = _need(_ALLEL_MODULE)
    result = allel_api.calculate_windowed_tajimas_d(POS, AC, window_size=500, verbose=False)
    values, windows, sites = allel.windowed_tajima_d(
        np.array(POS), allel.AlleleCountsArray(np.array(AC)), size=500
    )
    assert [row["start"] for row in result["windows"]] == [int(w[0]) for w in windows]
    assert [row["sites"] for row in result["windows"]] == [int(s) for s in sites]
    assert result["windows"][0]["tajima_d"] == pytest.approx(float(values[0]), abs=1e-9)


def test_windowed_nucleotide_diversity_matches_hand_calculation() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_windowed_nucleotide_diversity(
        POS, AC, window_size=500, verbose=False
    )
    first, second = result["windows"]
    # Per-site pairwise differences 0.5, 2/3, 0, 0.5, 0 summed over 500 callable bases.
    assert first["pi"] == pytest.approx((0.5 + 2 / 3 + 0 + 0.5 + 0) / 500)
    assert first["callable_bases"] == 500
    # Second window holds 401 bases (600..1000 inclusive).
    assert second["callable_bases"] == 401
    assert second["pi"] == pytest.approx((2 / 3 + 0.5 + 0.5 + 2 / 3 + 0) / 401)


def test_mean_pairwise_difference_hand_values() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_mean_pairwise_difference([[3, 1], [2, 2], [4, 0]], verbose=False)
    assert result["values"] == pytest.approx([0.5, 2 / 3, 0.0])


def test_hudson_fst_per_site_matches_formula() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_hudson_fst_per_site(AC1, AC2, verbose=False)
    # Second variant: p1 = 0.25, p2 = 0.75, n = 4 in both populations.
    p1, p2 = 0.25, 0.75
    numerator = (p1 - p2) ** 2 - p1 * (1 - p1) / 3 - p2 * (1 - p2) / 3
    denominator = p1 * (1 - p2) + p2 * (1 - p1)
    assert result["numerator"][1] == pytest.approx(numerator)
    assert result["denominator"][1] == pytest.approx(denominator)
    assert result["fst"] == pytest.approx((1 + 0.125 - 1 / 6 + 1) / 3.125)


def test_patterson_d_single_fixed_site_hand_value() -> None:
    _need(_ALLEL_MODULE)
    # Alt frequencies a=c=0 and b=d=1 give num = (a-b)(c-d) = 1 and den = 1 (allel's f4 form).
    result = allel_api.calculate_patterson_d([[4, 0]], [[0, 4]], [[4, 0]], [[0, 4]], verbose=False)
    assert result["numerator"] == [1.0]
    assert result["denominator"] == [1.0]
    assert result["d"] == pytest.approx(1.0)
    # Swapping C and D flips the sign of the single-site numerator.
    flipped = allel_api.calculate_patterson_d([[4, 0]], [[0, 4]], [[0, 4]], [[4, 0]], verbose=False)
    assert flipped["numerator"] == [-1.0]


def test_expected_heterozygosity_hand_values() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_expected_heterozygosity([[0.5, 0.5], [0.75, 0.25]], verbose=False)
    assert result["values"] == pytest.approx([0.5, 0.375])
    with pytest.raises(ConfigurationError):
        allel_api.calculate_expected_heterozygosity([[1.5, 0.5]], verbose=False)


def test_inbreeding_coefficient_hand_values() -> None:
    _need(_ALLEL_MODULE)
    genotypes = [
        [[0, 0], [0, 1]],
        [[1, 1], [0, 1]],
        [[0, 1], [1, 1]],
        [[0, 0], [0, 0]],
    ]
    result = allel_api.calculate_inbreeding_coefficient(genotypes, verbose=False)
    # Ho = 1/2 and He = 2 * (3/4) * (1/4) = 3/8 give F = 1 - (1/2) / (3/8) = -1/3.
    assert result["values"][:3] == pytest.approx([-1 / 3, -1 / 3, -1 / 3])
    assert result["values"][3] is None


def test_rogers_huff_r_hand_correlations() -> None:
    _need(_ALLEL_MODULE)
    dosage = [[0, 1], [2, 1], [1, 2], [0, 0]]
    result = allel_api.calculate_rogers_huff_r(dosage, verbose=False)
    pairs = {(row["variant_a"], row["variant_b"]): row["r"] for row in result["pairs"]}
    assert pairs[(0, 1)] == pytest.approx(-1.0)
    assert pairs[(0, 2)] == pytest.approx(1.0)
    assert pairs[(0, 3)] is None


def test_select_unlinked_variants_matches_upstream_mask() -> None:
    allel = _need(_ALLEL_MODULE)
    rng = np.random.default_rng(1)
    dosage = rng.integers(0, 3, size=(30, 8))
    result = allel_api.select_unlinked_variants(
        dosage.tolist(), size=10, step=5, threshold=0.1, verbose=False
    )
    mask = allel.locate_unlinked(dosage.astype("int8"), size=10, step=5, threshold=0.1)
    assert result["kept_indices"] == [int(i) for i in np.flatnonzero(mask)]


def test_pairwise_distance_hand_euclidean_matrix() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_pairwise_distance([[0, 1, 2], [1, 1, 2], [0, 0, 1]], verbose=False)
    matrix = result["matrix"]
    assert matrix[0][1] == pytest.approx(1.0)
    assert matrix[0][2] == pytest.approx(np.sqrt(6))
    assert matrix[1][2] == pytest.approx(np.sqrt(3))
    assert all(matrix[i][i] == 0 for i in range(3))


def _haplotypes(seed: int, variants: int = 40, samples: int = 12) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 2, size=(variants, samples))


def test_ihs_matches_upstream() -> None:
    allel = _need(_ALLEL_MODULE)
    haps = _haplotypes(3)
    positions = list(range(1000, 1000 + 100 * 40, 100))
    result = allel_api.calculate_ihs(haps.tolist(), positions, verbose=False)
    direct = allel.ihs(
        allel.HaplotypeArray(haps.astype("int8")), np.array(positions), use_threads=False
    )
    expected = [None if np.isnan(v) else round(float(v), 10) for v in direct]
    assert result["ihs"] == expected


def test_xpehh_matches_upstream() -> None:
    allel = _need(_ALLEL_MODULE)
    first, second = _haplotypes(4), _haplotypes(5)
    positions = np.array(range(1000, 1000 + 100 * 40, 100))
    result = allel_api.calculate_xpehh(
        first.tolist(), second.tolist(), positions.tolist(), verbose=False
    )
    direct = allel.xpehh(
        allel.HaplotypeArray(first.astype("int8")),
        allel.HaplotypeArray(second.astype("int8")),
        positions,
        use_threads=False,
    )
    expected = [None if np.isnan(v) else round(float(v), 10) for v in direct]
    assert result["xpehh"] == expected


def test_nsl_and_garud_h_match_upstream() -> None:
    allel = _need(_ALLEL_MODULE)
    haps = _haplotypes(6)
    nsl = allel_api.calculate_nsl(haps.tolist(), verbose=False)
    direct_nsl = allel.nsl(allel.HaplotypeArray(haps.astype("int8")), use_threads=False)
    assert nsl["nsl"] == [None if np.isnan(v) else round(float(v), 10) for v in direct_nsl]
    garud = allel_api.calculate_garud_h(haps.tolist(), verbose=False)
    h1, h12, h123, h2_h1 = allel.garud_h(allel.HaplotypeArray(haps.astype("int8")))
    assert garud["h1"] == pytest.approx(float(h1))
    assert garud["h12"] == pytest.approx(float(h12))
    assert garud["h123"] == pytest.approx(float(h123))
    assert garud["h2_h1"] == pytest.approx(float(h2_h1))


def test_ehh_decay_matches_upstream() -> None:
    allel = _need(_ALLEL_MODULE)
    haps = _haplotypes(7)
    result = allel_api.calculate_ehh_decay(haps.tolist(), verbose=False)
    direct = allel.ehh_decay(allel.HaplotypeArray(haps.astype("int8")))
    assert result["ehh"] == pytest.approx([round(float(v), 10) for v in direct])


def test_find_mendel_errors_counts_one_inconsistent_trio_site() -> None:
    _need(_ALLEL_MODULE)
    parents = [
        [[0, 0], [1, 1]],  # variant 0: mother 0/0, father 1/1
        [[0, 1], [0, 1]],  # variant 1: mother 0/1, father 0/1
        [[1, 1], [0, 0]],  # variant 2: mother 1/1, father 0/0
    ]
    progeny = [[[0, 1]], [[1, 1]], [[0, 0]]]
    result = allel_api.find_mendel_errors(parents, progeny, verbose=False)
    assert result["errors"] == [[0], [0], [1]]
    assert result["total_errors"] == 1


def test_phase_trio_by_transmission_matches_upstream_and_mendel() -> None:
    allel = _need(_ALLEL_MODULE)
    genotypes = [
        [[0, 1], [0, 0], [0, 1]],
        [[1, 1], [0, 1], [1, 1]],
        [[0, 0], [1, 0], [0, 0]],
        [[0, 1], [0, 1], [1, 1]],
    ]
    result = allel_api.phase_trio_by_transmission(genotypes, 2, verbose=False)
    direct = allel.phase_by_transmission(
        allel.GenotypeArray(np.array(genotypes, dtype="int8")), window_size=2
    )
    assert result["genotypes"] == direct.values.tolist()
    assert result["phased"] == direct.is_phased.tolist()
    # Phase-agnostic Mendelian check: the child carries one allele from each parent.
    for phased, (mother, father, _) in zip(result["genotypes"], genotypes, strict=False):
        child = sorted(phased[2])
        assert any(child == sorted([m, f]) for m in mother for f in father)


def test_call_roh_hmm_finds_long_homozygous_block() -> None:
    allel = _need(_ALLEL_MODULE)
    rng = np.random.default_rng(11)
    calls = rng.choice([[0, 0], [0, 1]], size=3000, p=[0.5, 0.5])
    calls[1000:2000] = 0  # long homozygous stretch
    positions = np.arange(1, 3001) * 10
    result = allel_api.call_roh_hmm(calls.tolist(), positions.tolist(), verbose=False)
    assert result["runs"], "expected at least one run of homozygosity"
    longest = max(result["runs"], key=lambda run: run["length"])
    assert longest["start"] <= 10_100 and longest["stop"] >= 19_900
    runs, froh = allel.roh_mhmm(
        allel.GenotypeVector(calls.astype("int8")), positions, contig_size=int(positions[-1])
    )
    assert [(int(r.start), int(r.stop)) for r in runs.itertuples()] == [
        (row["start"], row["stop"]) for row in result["runs"]
    ]
    assert result["froh"] == pytest.approx(float(froh))


def test_sequence_divergence_hand_formula() -> None:
    _need(_ALLEL_MODULE)
    result = allel_api.calculate_sequence_divergence(POS[:4], AC1, AC2, verbose=False)
    # Per-site between-population differences are 1, 0.625, 0.5, 1 over 400 - 100 + 1 bases.
    assert result["dxy"] == pytest.approx(3.125 / 301)
    assert (result["start"], result["stop"]) == (100, 400)


def test_population_branching_statistic_zero_for_identical_populations() -> None:
    allel = _need(_ALLEL_MODULE)
    result = allel_api.calculate_population_branching_statistic(
        AC, AC, AC, window_size=2, verbose=False
    )
    assert result["values"] == [0.0] * 5
    direct = allel.pbs(
        allel.AlleleCountsArray(np.array(AC)),
        allel.AlleleCountsArray(np.array(AC)),
        allel.AlleleCountsArray(np.array(AC)),
        window_size=2,
    )
    assert len(direct) == len(result["values"])


def test_invalid_allele_counts_raise_configuration_error() -> None:
    _need(_ALLEL_MODULE)
    with pytest.raises(ConfigurationError):
        allel_api.calculate_mean_pairwise_difference([[1, -2]], verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.calculate_windowed_nucleotide_diversity(
            [300, 200], AC[:2], window_size=50, verbose=False
        )


# ---------------------------------------------------------------- hgvs


def test_normalize_indel_alleles_left_and_right_shift() -> None:
    _need(_HGVS_MODULE)
    left = hgvs_api.normalize_indel_alleles("ACGTTTTA", 5, 6, ("T", ""), verbose=False)
    assert (left["start"], left["stop"], left["alleles"]) == (3, 4, ["T", ""])
    assert left["changed"] is True
    right = hgvs_api.normalize_indel_alleles(
        "ACGTTTTA", 3, 4, ("T", ""), left=False, bound=8, verbose=False
    )
    assert (right["start"], right["stop"]) == (6, 7)
    with pytest.raises(ConfigurationError):
        hgvs_api.normalize_indel_alleles("ACGT", 2, 1, ("A", ""), verbose=False)


def test_trim_common_prefix_strings_hand_case() -> None:
    _need(_HGVS_MODULE)
    result = hgvs_api.trim_common_prefix_strings(["ACGT", "ACTT"], verbose=False)
    assert result["trimmed_length"] == 2
    assert result["strings"] == ["GT", "TT"]


def test_build_transcript_cigar_joins_exons_with_introns() -> None:
    _need(_HGVS_MODULE)
    exons = [
        {"cigar": "10=", "alt_start_i": 0, "alt_end_i": 10},
        {"cigar": "5=", "alt_start_i": 20, "alt_end_i": 25},
    ]
    forward = hgvs_api.build_transcript_cigar(exons, 1, verbose=False)
    assert forward["cigar"] == "10=10N5="
    minus = hgvs_api.build_transcript_cigar(
        [{"cigar": "3M2=", "alt_start_i": 0, "alt_end_i": 5}, exons[1]], -1, verbose=False
    )
    assert minus["cigar"] == "2=3M15N5="
    with pytest.raises(ConfigurationError):
        hgvs_api.build_transcript_cigar(
            [{"cigar": "bad", "alt_start_i": 0, "alt_end_i": 1}], 1, verbose=False
        )


def test_validate_hgvs_type_accession_levels() -> None:
    _need(_HGVS_MODULE)
    import importlib

    importlib.import_module("hgvs.parser")
    validation = importlib.import_module("hgvs.utils.validation")
    valid = hgvs_api.validate_hgvs_type_accession("c", "NM_000088.3", verbose=False)
    assert valid["level"] == "VALID"
    assert (
        hgvs_api.validate_hgvs_type_accession("g", "NM_000088.3", verbose=False)["level"] == "ERROR"
    )
    assert (
        hgvs_api.validate_hgvs_type_accession("c", "NC_000001.11", verbose=False)["level"]
        == "ERROR"
    )
    level, _ = validation.validate_type_ac_pair("c", "NM_000088.3")
    assert level.name == valid["level"]
    with pytest.raises(ConfigurationError):
        hgvs_api.validate_hgvs_type_accession("x", "NM_000088.3", verbose=False)


def test_detect_and_assemble_repeat_blocks_round_trip() -> None:
    _need(_HGVS_MODULE)
    detected = hgvs_api.detect_repeat_blocks("ACACACGTTTT", verbose=False)
    assert detected["blocks"] == [
        {"repeat_unit": "AC", "repeat_count": 3, "block_size": 6, "block": "ACACAC"},
        {"repeat_unit": "T", "repeat_count": 4, "block_size": 4, "block": "TTTT"},
    ]
    assembled = hgvs_api.assemble_repeat_blocks("ACACACGTTTT", detected["blocks"], verbose=False)
    assert assembled["formatted"] == "AC[3]G[1]T[4]"
    with pytest.raises(ConfigurationError):
        hgvs_api.assemble_repeat_blocks("ACGT", [{"repeat_unit": "A"}], verbose=False)
