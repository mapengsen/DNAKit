"""Tests for dnakit.sequencing_qc.picard_metrics_p08_cli_functions (Picard metrics tools).

The tests run the real Picard (3.5.0 tested) on a small synthetic dataset. They skip when Picard
(`picard` on PATH, or `PICARD_JAR` plus `java`) is not installed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

import _picard_p08_data as data

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import picard_metrics_p08_cli_functions as m

requires_picard = pytest.mark.skipif(
    not data.picard_available(),
    reason="Picard (picard on PATH, or PICARD_JAR with java) is not installed",
)


@pytest.fixture(scope="module")
def files(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return data.write_dataset(tmp_path_factory.mktemp("picard_metrics"))


def _ok(result: dict) -> dict:
    assert result["status"] == "completed", result["log_tail"][-1500:]
    assert result["backend"] == "picard"
    assert result["version"] is not None
    return result


# ---------------------------------------------------------------------------
# Behaviour that needs no Picard
# ---------------------------------------------------------------------------


def test_missing_picard_gives_install_hint(files: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError) as error:
        m.picard_duplicate_metrics(
            str(files["sam"]),
            str(tmp_path / "dm.txt"),
            picard_path=str(tmp_path / "no_picard"),
            verbose=False,
        )
    assert "picard" in str(error.value).lower()


def test_existing_output_is_never_overwritten(files: dict[str, Path], tmp_path: Path) -> None:
    target = tmp_path / "dm.txt"
    target.write_text("keep me")
    with pytest.raises(ConfigurationError):
        m.picard_duplicate_metrics(str(files["sam"]), str(target), verbose=False)
    assert target.read_text() == "keep me"


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (
            lambda f, t: m.picard_accumulate_metrics(
                [str(f["sam"])], str(t / "a"), metric_type="nope", verbose=False
            ),
            "metric_type",
        ),
        (lambda f, t: m.picard_accumulate_metrics([], str(t / "a"), verbose=False), "input_paths"),
        (
            lambda f, t: m.picard_rnaseq_metrics(
                str(f["sam"]), str(t / "r"), str(f["sam"]), "SOMETIMES", verbose=False
            ),
            "strand_specificity",
        ),
        (
            lambda f, t: m.picard_base_distribution_by_cycle(
                str(f["sam"]), str(t / "b"), stop_after=-1, verbose=False
            ),
            "stop_after",
        ),
        (lambda f, t: m.picard_check_illumina_directory(str(t), [0], "2T", verbose=False), "lanes"),
        (
            lambda f, t: m.picard_duplicate_metrics(
                str(t / "missing.sam"), str(t / "d"), verbose=False
            ),
            "input_path",
        ),
        (
            lambda f, t: m.picard_duplicate_metrics(
                str(f["sam"]), str(t / "d"), validation_stringency="LAX", verbose=False
            ),
            "validation_stringency",
        ),
    ],
)
def test_arguments_are_validated(
    files: dict[str, Path], tmp_path: Path, call, message: str
) -> None:
    with pytest.raises(ConfigurationError) as error:
        call(files, tmp_path)
    assert message in str(error.value)


# ---------------------------------------------------------------------------
# Alignment metrics
# ---------------------------------------------------------------------------


@requires_picard
def test_duplicate_metrics_counts_pairs_and_matches_the_cli(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_duplicate_metrics(str(files["sam"]), str(tmp_path / "dm.txt"), verbose=False)
    )
    row = result["metrics"][0]
    assert row["LIBRARY"] == "lib1"
    assert row["READ_PAIRS_EXAMINED"] == 22  # 22 templates were written to the SAM
    assert row["READ_PAIR_DUPLICATES"] == 0  # the SAM has no duplicate flags
    assert result["metrics_class"].endswith("DuplicationMetrics")


@requires_picard
def test_base_distribution_by_cycle_has_one_row_per_cycle_and_read(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_base_distribution_by_cycle(
            str(files["sam"]), str(tmp_path / "bd.txt"), verbose=False
        )
    )
    assert result["cycle_count"] == 2 * data.READ_LENGTH
    percentages = [
        result["cycles"][0][name] for name in ("PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N")
    ]
    assert sum(percentages) == pytest.approx(100.0, abs=0.01)
    aligned = _ok(
        m.picard_base_distribution_by_cycle(
            str(files["sam"]),
            str(tmp_path / "bd2.txt"),
            aligned_reads_only=True,
            pf_reads_only=False,
            verbose=False,
        )
    )
    assert aligned["cycle_count"] == 2 * data.READ_LENGTH


@requires_picard
def test_oxog_metrics_reports_contexts(files: dict[str, Path], tmp_path: Path) -> None:
    result = _ok(
        m.picard_oxog_metrics(
            str(files["sam"]), str(tmp_path / "ox.txt"), str(files["fasta"]), verbose=False
        )
    )
    assert result["metrics"], "at least one CpCG context is observed in the example reads"
    first = result["metrics"][0]
    assert first["SAMPLE_ALIAS"] == "sample1"
    assert len(first["CONTEXT"]) == 3 and first["CONTEXT"][1] in "CG"
    assert first["TOTAL_BASES"] >= first["TOTAL_SITES"] > 0
    restricted = _ok(
        m.picard_oxog_metrics(
            str(files["sam"]),
            str(tmp_path / "ox2.txt"),
            str(files["fasta"]),
            contexts=[first["CONTEXT"]],
            verbose=False,
        )
    )
    assert {row["CONTEXT"] for row in restricted["metrics"]} == {first["CONTEXT"]}


@requires_picard
def test_snvq_yield_metrics_needs_the_snv_quality_tags(
    files: dict[str, Path], tmp_path: Path
) -> None:
    plain = m.picard_snvq_yield_metrics(
        str(files["sam"]), str(tmp_path / "plain.txt"), verbose=False
    )
    assert plain["status"] == "failed"
    assert "NullPointerException" in plain["log_tail"] or "Exception" in plain["log_tail"]
    lines = []
    for line in files["sam"].read_text().splitlines():
        if line.startswith("@"):
            lines.append(line)
        else:
            quality = line.split("\t")[10]
            lines.append(line + "".join(f"\tq{base}:Z:{quality}" for base in "acgt"))
    tagged = tmp_path / "snvq.sam"
    tagged.write_text("\n".join(lines) + "\n")
    result = _ok(
        m.picard_snvq_yield_metrics(str(tagged), str(tmp_path / "snvq.txt"), verbose=False)
    )
    row = result["metrics"][0]
    assert row["TOTAL_READS"] == 44
    assert row["TOTAL_BASES"] == 44 * data.READ_LENGTH
    assert row["PF_READS"] == 44


@requires_picard
def test_rnaseq_metrics_partitions_all_aligned_bases(
    files: dict[str, Path], tmp_path: Path
) -> None:
    refflat = tmp_path / "genes.refflat"
    refflat.write_text(
        "geneA\ttxA\tchr1\t+\t100\t900\t100\t900\t2\t100,500,\t300,900,\n"
        "geneB\ttxB\tchr2\t-\t50\t700\t50\t700\t1\t50,\t700,\n"
    )
    result = _ok(
        m.picard_rnaseq_metrics(
            str(files["sam"]),
            str(tmp_path / "rna.txt"),
            str(refflat),
            "NONE",
            end_bias_bases=50,
            verbose=False,
        )
    )
    row = result["metrics"][0]
    assert row["PF_BASES"] == 44 * data.READ_LENGTH
    assert (
        row["CODING_BASES"] + row["UTR_BASES"] + row["INTRONIC_BASES"] + row["INTERGENIC_BASES"]
        == row["PF_ALIGNED_BASES"]
    )


@requires_picard
def test_sam_error_metrics_writes_one_file_per_stratifier(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_sam_error_metrics(
            str(files["sam"]), str(tmp_path / "err"), str(files["fasta"]), verbose=False
        )
    )
    names = {Path(path).name for path in result["files"]}
    assert {"err.error_by_all", "err.error_by_base_quality", "err.error_by_cycle"} <= names
    table = [
        line.split("\t")
        for line in (tmp_path / "err.error_by_all").read_text().splitlines()
        if "\t" in line
    ]
    header = next(row for row in table if row[0] == "ERROR_BASES")
    values = table[table.index(header) + 1]
    # Errors can never exceed the bases examined.
    assert int(values[0]) <= int(values[header.index("TOTAL_BASES")])


@requires_picard
def test_sequencing_artifact_metrics_and_conversion_to_oxog(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_sequencing_artifact_metrics(
            str(files["sam"]), str(tmp_path / "art"), str(files["fasta"]), verbose=False
        )
    )
    assert set(result["outputs"]) == {
        "pre_adapter_detail_metrics",
        "pre_adapter_summary_metrics",
        "bait_bias_detail_metrics",
        "bait_bias_summary_metrics",
        "error_summary_metrics",
    }
    assert result["pre_adapter_summary"]
    converted = _ok(
        m.picard_convert_artifact_to_oxog(
            str(tmp_path / "art"), str(tmp_path / "oxo"), verbose=False
        )
    )
    assert Path(converted["outputs"]["oxog_metrics"]).name == "oxo.oxog_metrics"
    assert converted["metrics"], "OxoG rows are derived from the CpCG contexts"
    with pytest.raises(ConfigurationError):
        m.picard_convert_artifact_to_oxog(
            str(tmp_path / "absent"), str(tmp_path / "oxo2"), verbose=False
        )


@requires_picard
def test_wgs_metrics_with_non_zero_coverage_has_two_categories(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_wgs_metrics_nonzero_coverage(
            str(files["sam"]),
            str(tmp_path / "wgs.txt"),
            str(files["fasta"]),
            theoretical_sensitivity_output=str(tmp_path / "ts.txt"),
            coverage_cap=100,
            verbose=False,
        )
    )
    categories = {row["CATEGORY"]: row for row in result["metrics"]}
    assert set(categories) == {"WHOLE_GENOME", "NON_ZERO_REGIONS"}
    assert categories["WHOLE_GENOME"]["GENOME_TERRITORY"] == 2000  # 1200 + 800 bp reference
    assert (
        categories["NON_ZERO_REGIONS"]["MEAN_COVERAGE"]
        > categories["WHOLE_GENOME"]["MEAN_COVERAGE"]
    )
    assert "theoretical_sensitivity_output" in result["outputs"]


@requires_picard
def test_chart_output_without_r_is_reported_as_a_failure(
    files: dict[str, Path], tmp_path: Path
) -> None:
    import shutil

    if shutil.which("Rscript") is not None:
        pytest.skip("Rscript is installed, so the chart would be drawn")
    result = m.picard_wgs_metrics_nonzero_coverage(
        str(files["sam"]),
        str(tmp_path / "wgs.txt"),
        str(files["fasta"]),
        chart_output=str(tmp_path / "chart.pdf"),
        verbose=False,
    )
    assert result["status"] == "failed"
    assert "chart" not in result["outputs"]


# ---------------------------------------------------------------------------
# Variant-calling metrics
# ---------------------------------------------------------------------------


@requires_picard
def test_variant_calling_metrics_and_accumulation(files: dict[str, Path], tmp_path: Path) -> None:
    one = _ok(
        m.picard_variant_calling_metrics(
            str(files["vcf"]), str(files["dbsnp"]), str(tmp_path / "vc1"), verbose=False
        )
    )
    summary = one["summary_metrics"][0]
    assert summary["TOTAL_SNPS"] == 8
    assert summary["NUM_IN_DB_SNP"] == 8 and summary["PCT_DBSNP"] == 1
    assert one["detail_metrics"][0]["SAMPLE_ALIAS"] == "sample1"
    _ok(
        m.picard_variant_calling_metrics(
            str(files["vcf"]), str(files["dbsnp"]), str(tmp_path / "vc2"), verbose=False
        )
    )
    merged = _ok(
        m.picard_accumulate_metrics(
            [str(tmp_path / "vc1"), str(tmp_path / "vc2")],
            str(tmp_path / "vc_all"),
            metric_type="variant_calling",
            verbose=False,
        )
    )
    assert merged["summary_metrics"][0]["TOTAL_SNPS"] == 16
    assert set(merged["outputs"]) == {
        "variant_calling_detail_metrics",
        "variant_calling_summary_metrics",
    }


@requires_picard
def test_arrays_variant_calling_metrics(files: dict[str, Path], tmp_path: Path) -> None:
    import random

    rng = random.Random(7)
    reference = data.build_reference(rng)
    sites = data.fingerprint_sites(reference, data.build_pairs(rng, reference))
    array_vcf = data.write_arrays_vcf(tmp_path / "arrays.vcf", reference, sites)
    indexed = data.index_vcf(array_vcf)
    result = _ok(
        m.picard_arrays_variant_calling_metrics(
            str(indexed), str(files["dbsnp"]), str(tmp_path / "arr"), verbose=False
        )
    )
    summary = result["summary_metrics"][0]
    assert (
        summary["NUM_ASSAYS"] == 8 and summary["NUM_CALLS"] == 8 and summary["NUM_IN_DB_SNP"] == 8
    )
    assert summary["CALL_RATE"] == 1
    assert len(result["detail_metrics"]) == 1
    assert "arrays_control_code_summary_metrics" in result["outputs"]
    broken = m.picard_arrays_variant_calling_metrics(
        str(files["vcf"]), str(files["dbsnp"]), str(tmp_path / "arr2"), verbose=False
    )
    assert broken["status"] == "failed"  # a plain VCF lacks the Infinium header lines


@requires_picard
def test_accumulate_quality_yield_metrics_doubles_the_totals(
    files: dict[str, Path], tmp_path: Path
) -> None:
    from dnakit.alignment import picard_sam_p08_cli_functions as sam

    for name in ("a", "b"):
        _ok(
            sam.picard_downsample_sam(
                str(files["sam"]),
                str(tmp_path / f"{name}.sam"),
                metrics_file=str(tmp_path / f"{name}.txt"),
                verbose=False,
            )
        )
    merged = _ok(
        m.picard_accumulate_metrics(
            [str(tmp_path / "a.txt"), str(tmp_path / "b.txt")],
            str(tmp_path / "sum.txt"),
            verbose=False,
        )
    )
    assert merged["metrics"][0]["TOTAL_READS"] == 88
    assert merged["metrics"][0]["TOTAL_BASES"] == 88 * data.READ_LENGTH


@requires_picard
def test_format_contamination_metrics(tmp_path: Path) -> None:
    baf = tmp_path / "bafregress.txt"
    baf.write_text(
        "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\nS1\t0.012\t0.004\t3.0\t0.01\t0.99\t900\n"
    )
    result = _ok(
        m.picard_format_contamination_metrics(str(baf), str(tmp_path / "baf"), verbose=False)
    )
    row = result["metrics"][0]
    assert row["SAMPLE"] == "S1" and row["ESTIMATE"] == 0.012 and row["NHOM"] == 900
    assert row["LOG10_PVAL"] == -2
    assert Path(result["outputs"]["bafregress_metrics"]).name == "baf.bafregress_metrics"
    verify = tmp_path / "verify.txt"
    verify.write_text(
        "ID\t%Mix\tLLK\tLLK0\n--\n1\t0.05\t-120.5\t-125.5\n".replace("--\n", "----\n")
    )
    other = _ok(
        m.picard_format_contamination_metrics(
            str(verify), str(tmp_path / "vid"), source="verifyidintensity", verbose=False
        )
    )
    assert other["metrics"][0]["PCT_MIX"] == 0.05 and other["metrics"][0]["LLK"] == -120.5
    wrong = m.picard_format_contamination_metrics(str(verify), str(tmp_path / "bad"), verbose=False)
    assert wrong["status"] == "failed"


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


@requires_picard
def test_read_group_checksum_changes_with_the_read_groups(
    files: dict[str, Path], tmp_path: Path
) -> None:
    first = _ok(
        m.picard_read_group_checksum(str(files["sam"]), str(tmp_path / "one.txt"), verbose=False)
    )
    again = _ok(
        m.picard_read_group_checksum(str(files["sam"]), str(tmp_path / "two.txt"), verbose=False)
    )
    assert len(first["checksum"]) == 32 and first["checksum"] == again["checksum"]
    changed = tmp_path / "changed.sam"
    changed.write_text(files["sam"].read_text().replace("SM:sample1", "SM:sample2"))
    other = _ok(
        m.picard_read_group_checksum(str(changed), str(tmp_path / "three.txt"), verbose=False)
    )
    assert other["checksum"] != first["checksum"]


@requires_picard
def test_check_terminator_block_tells_complete_from_truncated_files(
    files: dict[str, Path], tmp_path: Path
) -> None:
    bam = data.sam_to_bam(files["sam"], tmp_path / "ok.bam")
    good = _ok(m.picard_check_terminator_block(str(bam), verbose=False))
    assert good["valid"] is True and good["return_code"] == 0
    truncated = tmp_path / "cut.bam"
    truncated.write_bytes(bam.read_bytes()[:1000])
    bad = _ok(m.picard_check_terminator_block(str(truncated), verbose=False))
    assert bad["valid"] is False and bad["return_code"] == 100


@requires_picard
def test_check_illumina_directory(tmp_path: Path) -> None:
    run = data.write_illumina_run(tmp_path / "run")
    good = _ok(
        m.picard_check_illumina_directory(
            str(run), [1], "2T", data_types=["BaseCalls"], verbose=False
        )
    )
    assert good["valid"] is True
    (run / "L001" / "C2.1" / "s_1_1101.bcl").unlink()
    bad = m.picard_check_illumina_directory(
        str(run), [1], "2T", data_types=["BaseCalls"], verbose=False
    )
    assert bad["valid"] is False and bad["status"] == "failed"
    with pytest.raises(ConfigurationError):
        m.picard_check_illumina_directory(
            str(run), [1], "2T", data_types=["Nonsense"], verbose=False
        )


# ---------------------------------------------------------------------------
# Fingerprinting
# ---------------------------------------------------------------------------


@requires_picard
def test_fingerprint_metrics_cover_every_haplotype(files: dict[str, Path], tmp_path: Path) -> None:
    result = _ok(
        m.picard_fingerprint_metrics(
            [str(files["sam"])],
            str(files["haplotype_map"]),
            str(tmp_path / "fp.txt"),
            calculate_by="SAMPLE",
            verbose=False,
        )
    )
    row = result["metrics"][0]
    assert row["SAMPLE_ALIAS"] == "sample1"
    assert row["HAPLOTYPES"] == 8
    assert 0 < row["HAPLOTYPES_WITH_EVIDENCE"] <= 8
    with pytest.raises(ConfigurationError):
        m.picard_fingerprint_metrics(
            [str(files["sam"])],
            str(files["haplotype_map"]),
            str(tmp_path / "fp2.txt"),
            calculate_by="CHROMOSOME",
            verbose=False,
        )


@requires_picard
def test_extract_fingerprint_writes_a_vcf_with_likelihoods(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_extract_fingerprint(
            str(files["sam"]),
            str(files["haplotype_map"]),
            str(tmp_path / "fp.vcf"),
            str(files["fasta"]),
            sample_alias="S9",
            verbose=False,
        )
    )
    assert result["sites"] == 8
    text = (tmp_path / "fp.vcf").read_text()
    assert "\tS9" in text and "PL" in text


@requires_picard
def test_check_fingerprint_compares_with_known_genotypes(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_check_fingerprint(
            str(files["sam"]),
            str(files["genotypes"]),
            str(files["haplotype_map"]),
            str(tmp_path / "chk"),
            verbose=False,
        )
    )
    summary = result["summary_metrics"][0]
    assert summary["SAMPLE"] == "sample1" and summary["HAPLOTYPES_WITH_GENOTYPES"] == 8
    assert isinstance(summary["LOD_EXPECTED_SAMPLE"], float)
    assert result["detail_metrics"], "one detail row per examined haplotype"
    missing = m.picard_check_fingerprint(
        str(files["sam"]),
        str(files["genotypes"]),
        str(files["haplotype_map"]),
        str(tmp_path / "chk2"),
        expected_sample_alias="nobody",
        verbose=False,
    )
    assert missing["status"] == "failed" and missing["return_code"] == 1


@requires_picard
def test_crosscheck_fingerprints_and_clustering(files: dict[str, Path], tmp_path: Path) -> None:
    same = _ok(
        m.picard_crosscheck_fingerprints(
            [str(files["sam"])],
            str(files["haplotype_map"]),
            str(tmp_path / "x1.txt"),
            crosscheck_by="SAMPLE",
            verbose=False,
        )
    )
    assert same["mismatch_detected"] is False and same["return_code"] == 0
    assert {row["RESULT"] for row in same["metrics"]} == {"EXPECTED_MATCH"}
    both = _ok(
        m.picard_crosscheck_fingerprints(
            [str(files["sam"]), str(files["genotypes"])],
            str(files["haplotype_map"]),
            str(tmp_path / "x2.txt"),
            crosscheck_by="FILE",
            matrix_output=str(tmp_path / "matrix.txt"),
            verbose=False,
        )
    )
    assert both["mismatch_detected"] == (both["return_code"] == 1)
    assert len(both["metrics"]) == 4  # every ordered pair of two groups
    assert "matrix" in both["outputs"]
    clustered = _ok(
        m.picard_cluster_crosscheck_metrics(
            str(tmp_path / "x2.txt"), str(tmp_path / "cl.txt"), lod_threshold=0.0, verbose=False
        )
    )
    assert {row["CLUSTER"] for row in clustered["metrics"]} >= {0}
    assert all("CLUSTER_SIZE" in row for row in clustered["metrics"])
