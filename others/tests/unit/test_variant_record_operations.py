"""Tests for dependency-free VCF record operations."""

from __future__ import annotations

import pytest

from dnakit.variants import (
    calculate_minor_allele_frequency,
    compare,
    count,
    filter_variants_by_allele_frequency,
    filter_variants_by_depth,
    filter_variants_by_genotype_quality,
    filter_variants_by_missingness,
    filter_variants_by_quality,
    parse_vcf_records,
    split_multiallelic,
)


VCF = "\n".join(
    (
        "##fileformat=VCFv4.3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
        "chr1\t10\trs1\tA\tG,T\t50\tPASS\tDP=12\tGT\t0/1\t1/2",
    )
)


def _records():
    return [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50, "genotypes": ["0/1", "1/1"]},
        {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "qual": 10, "genotypes": ["0/0", "./."]},
    ]


def test_parse_vcf_records_and_split_remaps_genotypes():
    parsed = parse_vcf_records(VCF, verbose=False)
    assert parsed["result"]["record_count"] == 1
    record = parsed["result"]["records"][0]
    assert record["alt"] == ["G", "T"]
    assert record["genotypes"] == ["0/1", "1/2"]

    split = split_multiallelic([record], verbose=False)["result"]
    assert split["output_record_count"] == 2
    assert [row["alt"] for row in split["records"]] == ["G", "T"]
    assert [row["genotypes"] for row in split["records"]] == [["0/1", "1/0"], ["0/0", "0/1"]]


def test_count_and_exact_compare():
    records = _records()
    summary = count(records, verbose=False)["result"]
    assert summary["record_count"] == 2
    assert summary["alternate_allele_count"] == 2
    assert summary["allele_type_counts"] == {"snp": 2}

    result = compare(
        {"left": {"variants": records}, "right": {"variants": records[:1]}},
        verbose=False,
    )["result"]
    assert result["shared_allele_count"] == 1
    assert result["only_left"] == [{"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"}]
    assert result["jaccard"] == pytest.approx(0.5)


def test_minor_allele_frequency_includes_zero_count_alts():
    result = calculate_minor_allele_frequency(
        {"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "genotypes": ["0/0"]}]},
        verbose=False,
    )["result"]
    row = result["variants"][0]
    assert row["allele_frequencies"] == {"0": 1.0, "1": 0.0}
    assert row["minor_allele"] == "G"
    assert row["maf"] == 0.0

    tied = calculate_minor_allele_frequency(
        {"variants": [{"chrom": "chr1", "pos": 2, "ref": "A", "alt": "G", "genotypes": ["0/1"]}]},
        verbose=False,
    )["result"]["variants"][0]
    assert tied["maf"] == 0.5
    assert tied["minor_allele_index"] == 0  # deterministic tie break


def test_quality_missingness_and_frequency_filters():
    records = _records()
    quality = filter_variants_by_quality({"variants": records}, options={"min_quality": 20}, verbose=False)["result"]
    assert quality["output_record_count"] == 1
    missing = filter_variants_by_missingness({"variants": records}, options={"max_missingness": 0.4}, verbose=False)["result"]
    assert missing["output_record_count"] == 1
    frequency = filter_variants_by_allele_frequency(
        {"variants": records}, options={"min_frequency": 0.4}, verbose=False
    )["result"]
    assert frequency["output_record_count"] == 1


def test_depth_and_genotype_quality_filters():
    records = {
        "variants": [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "DP": 30, "GQ": 45},
            {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "DP": 5, "GQ": 10},
        ]
    }
    depth = filter_variants_by_depth(records, options={"min_depth": 10}, verbose=False)["result"]
    assert depth["output_record_count"] == 1
    genotype_quality = filter_variants_by_genotype_quality(
        records, options={"min_genotype_quality": 20}, verbose=False
    )["result"]
    assert genotype_quality["output_record_count"] == 1


def test_filters_require_explicit_thresholds_and_calls():
    with pytest.raises(Exception):
        filter_variants_by_quality({"variants": _records()}, verbose=False)
    with pytest.raises(Exception):
        calculate_minor_allele_frequency({"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G"}]}, verbose=False)
