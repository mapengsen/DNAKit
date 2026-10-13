"""Hand-calculation and upstream-equivalence tests for the G18k review wrappers.

Covers varcode (nucleotide predicates, reference names, breakend junctions, layouts and local
VCF breakend pairing) and scikit-allel (population statistics, windows, SFS scaling and local
VCF archives). Every test runs offline on inline or tmp_path inputs.
"""

from __future__ import annotations

import inspect
import math
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.variants import allel_review_functions as allel_api
from dnakit.variants import varcode_review_functions as varcode_api

np = pytest.importorskip("numpy")

AC1 = [[4, 0], [2, 2], [0, 4], [3, 1], [4, 0], [1, 3]]
AC2 = [[0, 4], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2]]
AC3 = [[4, 0], [2, 2], [1, 3], [3, 1], [0, 4], [2, 2]]
POS = [10, 20, 30, 40, 50, 60]
HAPS = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 0]]
HAPS_B = [[0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 1, 0]]
GAC = [
    [[ref, 4 - ref] for ref in row]
    for row in [[4, 2, 0], [3, 1, 4], [2, 2, 1], [0, 4, 2], [4, 4, 0], [1, 3, 2]]
]

BND_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\tbnd_a\tA\tA]chr1:300]\t50\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT\t0/1
chr1\t300\tbnd_b\tT\t[chr1:100[T\t50\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT\t0/1
chr1\t500\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=600;SVLEN=-100\tGT\t0/1
"""

SMALL_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1
chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0
chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _need(name: str):
    return pytest.importorskip(name)


# ---------------------------------------------------------------- module contracts


@pytest.mark.parametrize("module", [varcode_api, allel_api])
def test_public_functions_take_keyword_only_verbose(module) -> None:
    for name in module.__all__:
        parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_export_counts_and_names_are_unique() -> None:
    assert len(varcode_api.__all__) == 10
    assert len(allel_api.__all__) == 31
    for module in (varcode_api, allel_api):
        assert len(set(module.__all__)) == len(module.__all__)


# ---------------------------------------------------------------- varcode: nucleotides


def test_is_purine_base_matches_hand_calls_and_rejects_extended_codes() -> None:
    _need("varcode")
    assert varcode_api.is_purine_base("a", verbose=False)["is_purine"] is True
    assert varcode_api.is_purine_base("G", verbose=False)["is_purine"] is True
    assert varcode_api.is_purine_base("C", verbose=False)["is_purine"] is False
    assert varcode_api.is_purine_base("N", allow_extended=True, verbose=False)["is_purine"] is False
    with pytest.raises(ConfigurationError):
        varcode_api.is_purine_base("N", verbose=False)
    with pytest.raises(ConfigurationError):
        varcode_api.is_purine_base("AC", verbose=False)


def test_has_only_standard_nucleotides_reports_one_based_positions() -> None:
    _need("varcode")
    mixed = varcode_api.has_only_standard_nucleotides("ACGTNacgt", verbose=False)
    assert mixed["standard"] is False
    assert mixed["non_standard_positions"] == [5]
    assert mixed["length"] == 9
    plain = varcode_api.has_only_standard_nucleotides("acgt", verbose=False)
    assert plain["standard"] is True
    assert plain["non_standard_positions"] == []
    with pytest.raises(ConfigurationError):
        varcode_api.has_only_standard_nucleotides("", verbose=False)
    with pytest.raises(ConfigurationError):
        varcode_api.has_only_standard_nucleotides("AC GT", verbose=False)


def test_is_mitochondrial_contig_recognises_common_names() -> None:
    _need("varcode")
    assert varcode_api.is_mitochondrial_contig("chrM", verbose=False)["mitochondrial"] is True
    assert varcode_api.is_mitochondrial_contig("MT", verbose=False)["mitochondrial"] is True
    assert varcode_api.is_mitochondrial_contig("chr1", verbose=False)["mitochondrial"] is False
    assert varcode_api.is_mitochondrial_contig(None, verbose=False)["mitochondrial"] is False


# ---------------------------------------------------------------- varcode: reference names


def test_reference_name_predicates_and_inference() -> None:
    _need("varcode")
    assert varcode_api.is_ucsc_reference_name("hg19", verbose=False)["ucsc"] is True
    assert varcode_api.is_ucsc_reference_name("GRCh37", verbose=False)["ucsc"] is False
    assert (
        varcode_api.infer_reference_assembly("/data/refs/GRCh37.fa", verbose=False)["assembly"]
        == "GRCh37"
    )
    assert varcode_api.infer_reference_assembly("mm10", verbose=False)["assembly"] == "mm10"
    with pytest.raises(ConfigurationError):
        varcode_api.infer_reference_assembly("no_such_genome", verbose=False)


# ---------------------------------------------------------------- varcode: breakends


def test_breakend_junction_span_matches_upstream_and_hand_reading() -> None:
    structural = _need("varcode.structural_variant")
    deletion = varcode_api.breakend_junction_span(
        {"contig": "chr1", "position": 100, "keeps": "left"},
        {"contig": "chr1", "position": 200, "keeps": "right"},
        verbose=False,
    )["span"]
    # Bases 101 to 199 are removed; bases 100 and 200 stay as the junction padding.
    assert deletion == {
        "sv_type": "DEL",
        "start": 100,
        "end": 199,
        "affected_start": 101,
        "affected_end": 199,
    }
    upstream = structural.junction_span_details(
        (structural.Breakend("chr1", 100, "right"), structural.Breakend("chr1", 200, "left"))
    )
    duplication = varcode_api.breakend_junction_span(
        {"contig": "chr1", "position": 100, "keeps": "right"},
        {"contig": "chr1", "position": 200, "keeps": "left"},
        verbose=False,
    )["span"]
    assert duplication["sv_type"] == upstream[0]
    assert (duplication["start"], duplication["end"]) == (upstream[1], upstream[2])
    unknown = varcode_api.breakend_junction_span(
        {"contig": "chr1", "position": 100},
        {"contig": "chr1", "position": 200, "keeps": "left"},
        verbose=False,
    )
    assert unknown["span"] is None
    cross = varcode_api.breakend_junction_span(
        {"contig": "chr1", "position": 100, "keeps": "left"},
        {"contig": "chr2", "position": 200, "keeps": "left"},
        verbose=False,
    )
    assert cross["span"] is None


def test_join_breakend_layouts_hand_calculation() -> None:
    _need("varcode")
    result = varcode_api.join_breakend_layouts(
        [{"contig": "chr1", "start": 100, "end": 110}],
        {"contig": "chr1", "position": 105, "keeps": "left"},
        [{"contig": "chr1", "start": 300, "end": 310}],
        {"contig": "chr1", "position": 302, "keeps": "right"},
        verbose=False,
    )
    # Left piece keeps 100..105 (6 bases); right piece keeps 302..310 (9 bases).
    assert [(s["start"], s["end"]) for s in result["segments"]] == [(100, 105), (302, 310)]
    assert result["length"] == 15
    with pytest.raises(ConfigurationError):
        varcode_api.join_breakend_layouts(
            [{"contig": "chr1", "start": 100, "end": 110}],
            {"contig": "chr1", "position": 105},
            [{"contig": "chr1", "start": 300, "end": 310}],
            {"contig": "chr1", "position": 302, "keeps": "right"},
            verbose=False,
        )


def test_project_intervals_onto_layout_hand_calculation() -> None:
    _need("varcode")
    result = varcode_api.project_intervals_onto_layout(
        [{"contig": "chr1", "start": 100, "end": 110}],
        [
            {"contig": "chr1", "start": 103, "end": 120},
            {"contig": "chr2", "start": 1, "end": 5},
        ],
        verbose=False,
    )
    # Overlap of 103..120 with 100..110 is 103..110, which is 8 bases.
    assert result["count"] == 1
    piece = result["pieces"][0]
    assert piece["segments"] == [
        {"contig": "chr1", "start": 103, "end": 110, "strand": "+", "length": 8}
    ]
    assert piece["length"] == 8
    with pytest.raises(ConfigurationError):
        varcode_api.project_intervals_onto_layout(
            [{"contig": "chr1", "start": 100}], [], verbose=False
        )


def test_pair_breakend_mates_merges_local_bnd_records(tmp_path: Path) -> None:
    _need("varcode")
    source = _write(tmp_path, "bnd.vcf", BND_VCF)
    result = varcode_api.pair_breakend_mates(str(source), verbose=False)
    assert result["input_count"] == 3
    assert result["output_count"] == 2
    assert result["paired_count"] == 1
    assert result["records"][0] == {
        "contig": "chr1",
        "start": 100,
        "ref": "A",
        "alt": "A]chr1:300]",
        "source_starts": [100, 300],
    }
    assert result["records"][1]["source_starts"] == []
    with pytest.raises(ConfigurationError):
        varcode_api.pair_breakend_mates(str(tmp_path / "missing.vcf"), verbose=False)


def test_rewrite_vcf_records_writes_a_fresh_file_once(tmp_path: Path) -> None:
    _need("varcode")
    source = _write(tmp_path, "small.vcf", SMALL_VCF)
    output = tmp_path / "rewritten.vcf"
    result = varcode_api.rewrite_vcf_records(str(source), str(output), verbose=False)
    assert result["variants"] == 3
    text = output.read_text(encoding="utf-8")
    assert text.startswith("##fileformat=VCFv4.2")
    assert "chr1\t10\trs1\tA\tG\t50.0\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1" in text
    with pytest.raises(ConfigurationError):
        varcode_api.rewrite_vcf_records(str(source), str(output), verbose=False)


# ---------------------------------------------------------------- allel: private and fixed


def test_locate_population_private_alleles_hand_calculation() -> None:
    _need("allel")
    result = allel_api.locate_population_private_alleles([AC1, AC2, AC3], verbose=False)
    # Variant 0: allele 1 is only in population 2. Variant 4: allele 0 is only in population 1.
    assert result["private_alleles"] == [
        {"variant_index": 0, "allele_index": 1},
        {"variant_index": 4, "allele_index": 0},
    ]
    assert result["private_variant_indices"] == [0, 4]
    with pytest.raises(ConfigurationError):
        allel_api.locate_population_private_alleles([AC1], verbose=False)


def test_locate_population_fixed_differences_hand_calculation() -> None:
    _need("allel")
    result = allel_api.locate_population_fixed_differences(AC1, AC2, verbose=False)
    # Variants 0, 2 and 4 have no shared allele; the others carry an allele in both.
    assert result["fixed_indices"] == [0, 2, 4]
    assert result["fixed_count"] == 3


# ---------------------------------------------------------------- allel: windows


def test_count_variants_in_windows_hand_calculation() -> None:
    _need("allel")
    result = allel_api.count_variants_in_windows(POS, window_size=20, verbose=False)
    # Windows of 20 bases from position 10: [10,29] holds 10 and 20; [30,49] holds 30 and 40.
    assert [(w["start"], w["stop"]) for w in result["windows"]][:2] == [(10, 29), (30, 49)]
    assert [w["variants"] for w in result["windows"]][:2] == [2, 2]
    assert result["step"] == 20
    with pytest.raises(ConfigurationError):
        allel_api.count_variants_in_windows([5, 3], window_size=10, verbose=False)


def test_windowed_values_per_base_hand_calculation() -> None:
    _need("allel")
    result = allel_api.windowed_values_per_base([1.0, 2.0], [[1, 3], [4, 6]], verbose=False)
    # Each window holds three bases, so the values are 1/3 and 2/3.
    assert [w["callable_bases"] for w in result["windows"]] == [3, 3]
    assert result["windows"][0]["value_per_base"] == pytest.approx(1 / 3)
    assert result["windows"][1]["value_per_base"] == pytest.approx(2 / 3)
    masked = allel_api.windowed_values_per_base(
        [1.0, 2.0],
        [[1, 3], [4, 6]],
        is_accessible=[True, False, True, True, True, True],
        verbose=False,
    )
    assert [w["callable_bases"] for w in masked["windows"]] == [2, 3]


def test_build_equally_accessible_windows_hand_calculation() -> None:
    _need("allel")
    result = allel_api.build_equally_accessible_windows(
        [True, True, False, True, True, True], window_size=2, verbose=False
    )
    # Accessible bases are 1, 2, 4, 5 and 6; two-base windows take {1, 2} and {4, 5}.
    assert [(w["start"], w["stop"]) for w in result["windows"]] == [(1, 2), (4, 5)]
    assert [w["accessible_bases"] for w in result["windows"]] == [2, 2]


# ---------------------------------------------------------------- allel: Fst and divergence


def test_between_population_pi_hand_calculation() -> None:
    _need("allel")
    values = allel_api.calculate_between_population_pi(AC1, AC2, verbose=False)["values"]
    # pi_between = (n1_ref * n2_alt + n1_alt * n2_ref) / (n1 * n2) with n1 = n2 = 4.
    assert values[0] == pytest.approx(1.0)
    assert values[1] == pytest.approx(0.5)
    assert values[3] == pytest.approx(0.625)


def test_windowed_fixed_difference_density_hand_calculation() -> None:
    _need("allel")
    result = allel_api.calculate_windowed_fixed_difference_density(
        POS, AC1, AC2, window_size=20, verbose=False
    )
    # The first window holds one fixed difference (position 10) over 20 callable bases.
    assert result["windows"][0]["sites"] == 2
    assert result["windows"][0]["callable_bases"] == 20
    assert result["windows"][0]["density"] == pytest.approx(1 / 20)


def test_patterson_fst_matches_upstream_and_fixed_variants() -> None:
    allel = _need("allel")
    result = allel_api.calculate_patterson_fst(AC1, AC2, verbose=False)
    numerator, denominator = allel.patterson_fst(
        allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2)
    )
    assert result["numerator"] == pytest.approx(numerator.tolist())
    assert result["denominator"] == pytest.approx(denominator.tolist())
    assert result["numerator"][0] == pytest.approx(1.0)
    assert result["fst"] == pytest.approx(float(numerator.sum() / denominator.sum()))


def test_windowed_patterson_fst_matches_upstream() -> None:
    allel = _need("allel")
    result = allel_api.calculate_windowed_patterson_fst(
        POS, AC1, AC2, window_size=20, verbose=False
    )
    fst, _windows, counts = allel.windowed_patterson_fst(
        allel.SortedIndex(POS), allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2), size=20
    )
    assert [w["fst"] for w in result["windows"]] == pytest.approx(fst.tolist())
    assert [w["sites"] for w in result["windows"]] == counts.tolist()


def test_blockwise_and_moving_patterson_fst_match_upstream() -> None:
    allel = _need("allel")
    blocked = allel_api.calculate_blockwise_patterson_fst(AC1, AC2, block_length=2, verbose=False)
    fst, error, blocks, jackknife = allel.blockwise_patterson_fst(
        allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2), 2
    )
    assert blocked["fst"] == pytest.approx(float(fst))
    assert blocked["standard_error"] == pytest.approx(float(error))
    assert blocked["block_values"] == pytest.approx(blocks.tolist())
    assert blocked["jackknife_values"] == pytest.approx(jackknife.tolist())
    moving = allel_api.calculate_moving_patterson_fst(AC1, AC2, window_size=2, verbose=False)
    upstream = allel.moving_patterson_fst(
        allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2), size=2
    )
    assert moving["values"] == pytest.approx(upstream.tolist())


def test_pairwise_dxy_is_symmetric_and_matches_upstream() -> None:
    allel = _need("allel")
    result = allel_api.calculate_pairwise_dxy(POS, GAC, verbose=False)
    matrix = result["matrix"]
    assert len(matrix) == 3
    for i in range(3):
        assert matrix[i][i] == 0.0
        for j in range(3):
            assert matrix[i][j] == matrix[j][i]
    condensed = allel.pairwise_dxy(allel.SortedIndex(POS), np.asarray(GAC))
    assert [matrix[0][1], matrix[0][2], matrix[1][2]] == pytest.approx(condensed.tolist())


def test_patterson_f2_and_f3_match_upstream() -> None:
    allel = _need("allel")
    f2 = allel_api.calculate_patterson_f2(AC1, AC2, verbose=False)
    upstream_f2 = allel.patterson_f2(allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2))
    assert f2["values"] == pytest.approx(upstream_f2.tolist())
    # Values are 1, -1/6, 1, 1/8, 1 and -1/12, so their mean is 2.875 / 6.
    assert f2["mean"] == pytest.approx(2.875 / 6)
    f3 = allel_api.calculate_patterson_f3(AC3, AC1, AC2, verbose=False)
    numerator, heterozygosity = allel.patterson_f3(
        allel.AlleleCountsArray(AC3), allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2)
    )
    assert f3["numerator"] == pytest.approx(numerator.tolist())
    assert f3["f3_star"] == pytest.approx(float(numerator.sum() / heterozygosity.sum()))


def test_blockwise_and_moving_patterson_f3_match_upstream() -> None:
    allel = _need("allel")
    blocked = allel_api.calculate_blockwise_patterson_f3(
        AC3, AC1, AC2, block_length=2, verbose=False
    )
    f3, error, z_score, blocks, _ = allel.blockwise_patterson_f3(
        allel.AlleleCountsArray(AC3),
        allel.AlleleCountsArray(AC1),
        allel.AlleleCountsArray(AC2),
        2,
    )
    assert blocked["f3"] == pytest.approx(float(f3))
    assert blocked["standard_error"] == pytest.approx(float(error))
    assert blocked["z_score"] == pytest.approx(float(z_score))
    assert blocked["block_values"] == pytest.approx(blocks.tolist())
    moving = allel_api.calculate_moving_patterson_f3(AC3, AC1, AC2, window_size=2, verbose=False)
    upstream = allel.moving_patterson_f3(
        allel.AlleleCountsArray(AC3), allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2), 2
    )
    assert moving["values"] == pytest.approx(upstream.tolist())


# ---------------------------------------------------------------- allel: haplotypes


def test_haplotype_diversity_hand_calculation() -> None:
    _need("allel")
    result = allel_api.calculate_haplotype_diversity(HAPS, verbose=False)
    # Four distinct haplotypes with frequency 1/4: (1 - 4/16) * 4/3 = 1.
    assert result["diversity"] == pytest.approx(1.0)
    assert result["haplotypes"] == 4


def test_paint_shared_haplotype_prefixes_matches_upstream() -> None:
    allel = _need("allel")
    result = allel_api.paint_shared_haplotype_prefixes(HAPS, verbose=False)
    painting, order = allel.voight_painting(np.asarray(HAPS, dtype="int8"))
    assert result["painting"] == painting.tolist()
    assert result["haplotype_order"] == order.tolist()


def test_moving_haplotype_diversity_and_xpnsl_match_upstream() -> None:
    allel = _need("allel")
    moving = allel_api.calculate_moving_haplotype_diversity(HAPS, window_size=3, verbose=False)
    upstream = allel.moving_haplotype_diversity(
        allel.HaplotypeArray(np.asarray(HAPS, dtype="int8")), 3
    )
    assert moving["values"] == pytest.approx(upstream.tolist())
    xpnsl = allel_api.calculate_xpnsl(HAPS, HAPS_B, verbose=False)
    expected = allel.xpnsl(
        allel.HaplotypeArray(np.asarray(HAPS, dtype="int8")),
        allel.HaplotypeArray(np.asarray(HAPS_B, dtype="int8")),
        use_threads=False,
    )
    assert xpnsl["xpnsl"] == pytest.approx(expected.tolist())


def test_paint_parental_transmission_hand_calculation() -> None:
    _need("allel")
    result = allel_api.paint_parental_transmission(
        [[0, 1], [1, 1], [0, 0], [0, -1]],
        [[0, 1, 0], [1, 0, 1], [0, 0, 0], [1, 1, 1]],
        verbose=False,
    )
    # Heterozygous parent: progeny alleles match haplotype 1 (code 1) or haplotype 2 (code 2).
    assert result["painting"][0] == [1, 2, 1]
    # Homozygous non-reference parent: non-reference allele on both haplotypes (code 4).
    assert result["painting"][1] == [4, 5, 4]
    # Homozygous reference parent: reference allele on both haplotypes (code 3).
    assert result["painting"][2] == [3, 3, 3]
    # Missing parental allele (code 6) when one parental haplotype is missing.
    assert result["painting"][3] == [6, 6, 6]


# ---------------------------------------------------------------- allel: selection and SFS


def test_standardize_score_by_allele_count_matches_upstream() -> None:
    allel = _need("allel")
    score = [0.1, -0.2, 0.3, 0.5, -0.4, 0.2, 0.0, 1.0]
    counts = [1, 2, 3, 4, 5, 6, 7, 8]
    result = allel_api.standardize_score_by_allele_count(score, counts, n_bins=2, verbose=False)
    values, bins = allel.standardize_by_allele_count(
        np.asarray(score), np.asarray(counts), n_bins=2, diagnostics=False
    )
    assert result["values"] == pytest.approx(values.tolist())
    assert result["bins"] == bins.tolist()
    with pytest.raises(ConfigurationError):
        allel_api.standardize_score_by_allele_count(score, counts, bins=[4, 2], verbose=False)


def test_moving_delta_tajimas_d_matches_upstream_and_reports_undefined_windows() -> None:
    allel = _need("allel")
    result = allel_api.calculate_moving_delta_tajimas_d(AC1, AC2, window_size=3, verbose=False)
    upstream = allel.moving_delta_tajima_d(
        allel.AlleleCountsArray(AC1), allel.AlleleCountsArray(AC2), 3
    )
    expected = [None if math.isnan(v) else pytest.approx(v) for v in upstream.tolist()]
    assert result["values"] == expected


def test_fold_joint_site_frequency_spectrum_hand_calculation() -> None:
    _need("allel")
    spectrum = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]
    result = allel_api.fold_joint_site_frequency_spectrum(
        spectrum, first_sample_size=2, second_sample_size=3, verbose=False
    )
    # Corner-to-corner sums: [0+4+3+7, 1+5+2+6] = [14, 14].
    assert result["folded"] == [[14, 14]]
    with pytest.raises(ConfigurationError):
        allel_api.fold_joint_site_frequency_spectrum(
            spectrum, first_sample_size=1, second_sample_size=3, verbose=False
        )


def test_scale_site_frequency_spectrum_hand_calculation() -> None:
    _need("allel")
    result = allel_api.scale_site_frequency_spectrum([5, 3, 2, 1], verbose=False)
    # Class i is multiplied by i: 0*5, 1*3, 2*2, 3*1.
    assert result["scaled"] == [0, 3, 4, 3]


def test_scale_folded_site_frequency_spectrum_hand_calculation() -> None:
    _need("allel")
    result = allel_api.scale_folded_site_frequency_spectrum([5, 3, 2], sample_size=6, verbose=False)
    # Class k is multiplied by k(n - k) / n with n = 6: 0, 1*5/6*3 = 2.5, 2*4/6*2 = 8/3.
    assert result["scaled"] == pytest.approx([0.0, 2.5, 8 / 3])


def test_scale_joint_site_frequency_spectrum_hand_calculation() -> None:
    _need("allel")
    result = allel_api.scale_joint_site_frequency_spectrum([[0, 1, 2], [3, 4, 5]], verbose=False)
    # Cell (i, j) is multiplied by i * j: row 1 becomes [0, 4*1*1, 5*1*2].
    assert result["scaled"] == [[0, 0, 0], [0, 4, 10]]


def test_scale_folded_joint_site_frequency_spectrum_hand_calculation() -> None:
    _need("allel")
    result = allel_api.scale_folded_joint_site_frequency_spectrum(
        [[0, 1, 2], [3, 4, 5]], first_sample_size=2, second_sample_size=3, verbose=False
    )
    # Cell (i, j) is multiplied by i j (n1 - i)(n2 - j): row 1 becomes [0, 4*1*1*1*2, 5*1*2*1*1].
    assert result["scaled"] == [[0, 0, 0], [0, 8, 10]]


# ---------------------------------------------------------------- allel: local VCF files


def test_read_vcf_header_summary_reads_samples_fields_and_contigs(tmp_path: Path) -> None:
    _need("allel")
    source = _write(tmp_path, "small.vcf", SMALL_VCF)
    result = allel_api.read_vcf_header_summary(str(source), verbose=False)
    assert result["samples"] == ["S1", "S2", "S3"]
    assert result["info_fields"] == ["DP"]
    assert result["format_fields"] == ["GT"]
    assert result["contigs"] == [{"id": "chr1", "length": 1000}]


def test_vcf_to_npz_and_hdf5_archives_report_shapes(tmp_path: Path) -> None:
    _need("allel")
    source = _write(tmp_path, "small.vcf", SMALL_VCF)
    fields = ["variants/POS", "calldata/GT"]
    npz = allel_api.vcf_to_npz_archive(
        str(source), str(tmp_path / "small.npz"), fields=fields, verbose=False
    )
    assert npz["arrays"] == {"calldata/GT": [3, 3, 2], "variants/POS": [3]}
    h5 = allel_api.vcf_to_hdf5_archive(
        str(source), str(tmp_path / "small.h5"), fields=fields, verbose=False
    )
    assert h5["datasets"] == {"calldata/GT": [3, 3, 2], "variants/POS": [3]}
    with pytest.raises(ConfigurationError):
        allel_api.vcf_to_npz_archive(str(source), str(tmp_path / "small.h5"), verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.vcf_to_npz_archive(
            str(source), str(tmp_path / "small.npz"), fields=fields, verbose=False
        )


# ---------------------------------------------------------------- validation


def test_allel_wrappers_reject_malformed_inputs() -> None:
    _need("allel")
    with pytest.raises(ConfigurationError):
        allel_api.calculate_patterson_fst([[1.5, 2.0]], [[1, 1]], verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.calculate_patterson_fst([[1, 1]], [[1, 1], [2, 2]], verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.count_variants_in_windows([9, 3], window_size=10, verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.calculate_haplotype_diversity([[0, 2], [1, 0]], verbose=False)
    with pytest.raises(ConfigurationError):
        allel_api.scale_folded_site_frequency_spectrum([1, 2, 3, 4], sample_size=2, verbose=False)
