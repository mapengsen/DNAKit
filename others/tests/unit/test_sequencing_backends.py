"""Contract tests for cutadapt, umi-tools and pyslow5 sequencing QC backends."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc.backends import (
    calculate_expected_errors,
    cluster_umi_counts,
    extract_read_barcodes,
    find_adapter_match,
    find_quality_trim_bounds,
    map_barcodes_to_whitelist,
    read_slow5_signal,
    summarize_slow5_signals,
    trim_fastq_adapters,
    write_slow5_signals,
)

ADAPTER = "AGATCGGAAGAGC"
FASTQ_TEXT = (
    "@r1\nACGTTTACGTAGATCGGAAGAGCTTT\n+\n" + "I" * 26 + "\n"
    "@r2\nGGGGGGGGGGGGGGGG\n+\n" + "I" * 16 + "\n"
    "@r3\nAGATCGGAAGAGCACGT\n+\n" + "I" * 17 + "\n"
)


def _write(tmp_path: Path, text: str, name: str = "reads.fastq") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_find_adapter_match_matches_hand_calculation() -> None:
    # The 13-base adapter starts at index 10 of a 26-base read and has no errors.
    result = find_adapter_match("ACGTTTACGTAGATCGGAAGAGCTTT", ADAPTER, verbose=False)
    assert result["matched"] is True
    assert (result["start"], result["stop"]) == (10, 23)
    assert result["errors"] == 0
    assert result["score"] == 13  # cutadapt scores +1 per matched base
    assert result["trimmed_sequence"] == "ACGTTTACGT"


def test_find_adapter_match_agrees_with_cutadapt_adapter_object() -> None:
    adapters = pytest.importorskip("cutadapt.adapters")
    read = "TTTTACGTAGATCGGAAGAGCAAA"
    expected = adapters.BackAdapter(ADAPTER, max_errors=0.1, min_overlap=3).match_to(read)

    result = find_adapter_match(read, ADAPTER, verbose=False)

    assert (result["start"], result["stop"]) == (expected.rstart, expected.rstop)
    assert result["errors"] == expected.errors
    assert result["trimmed_sequence"] == expected.trimmed(read)


def test_find_adapter_match_five_prime_and_anywhere_positions() -> None:
    pytest.importorskip("cutadapt")
    five = find_adapter_match("GGGGACGT", "GGGG", position="5prime", max_errors=0.0, verbose=False)
    assert (five["start"], five["stop"]) == (0, 4)
    assert five["trimmed_sequence"] == "ACGT"

    anywhere = find_adapter_match(
        "TTTTACGTAGATCGGAAGAGCAAA", ADAPTER, position="anywhere", verbose=False
    )
    assert (anywhere["start"], anywhere["stop"]) == (8, 21)
    assert anywhere["trimmed_sequence"] == "TTTTACGT"


def test_find_adapter_match_error_rate_controls_mismatches() -> None:
    # One mismatch (T for C at adapter position 8) is allowed at rate 0.1 but not at 0.
    read = "ACGTTTACGTAGATCGGTAGAGCTTT"
    tolerant = find_adapter_match(read, ADAPTER, max_errors=0.1, verbose=False)
    strict = find_adapter_match(read, ADAPTER, max_errors=0.0, verbose=False)
    assert tolerant["matched"] is True and tolerant["errors"] == 1
    assert strict["matched"] is False
    assert strict["trimmed_sequence"] == read


def test_find_adapter_match_rejects_bad_input() -> None:
    pytest.importorskip("cutadapt")
    with pytest.raises(ConfigurationError, match="non-empty DNA sequence"):
        find_adapter_match("", ADAPTER, verbose=False)
    with pytest.raises(ConfigurationError, match="position must be one of"):
        find_adapter_match("ACGT", ADAPTER, position="middle", verbose=False)
    with pytest.raises(ConfigurationError, match="max_errors"):
        find_adapter_match("ACGT", ADAPTER, max_errors=1.0, verbose=False)
    with pytest.raises(ConfigurationError, match="max_errors"):
        find_adapter_match("ACGT", ADAPTER, max_errors=True, verbose=False)
    with pytest.raises(ConfigurationError, match="min_overlap"):
        find_adapter_match("ACGT", ADAPTER, min_overlap=True, verbose=False)
    with pytest.raises(ConfigurationError, match="unsupported DNA symbols"):
        find_adapter_match("ACGU", ADAPTER, verbose=False)


def test_trim_fastq_adapters_counts_reads_and_writes_retained_records(tmp_path: Path) -> None:
    pytest.importorskip("cutadapt")
    source = _write(tmp_path, FASTQ_TEXT)
    output = tmp_path / "trimmed.fastq"

    result = trim_fastq_adapters(source, output, ADAPTER, verbose=False)

    # r1 loses 16 bases (13 adapter + 3 tail), r3 loses all 17 bases, r2 is untouched.
    assert result["reads"] == 3
    assert result["reads_with_adapter"] == 2
    assert result["reads_too_short"] == 1  # r3 becomes empty and is discarded
    assert result["reads_written"] == 2
    assert result["bases_removed"] == 33
    assert output.read_text(encoding="utf-8") == (
        "@r1\nACGTTTACGT\n+\n" + "I" * 10 + "\n" + "@r2\n" + "G" * 16 + "\n+\n" + "I" * 16 + "\n"
    )


def test_trim_fastq_adapters_enforces_limits_and_keeps_existing_files(tmp_path: Path) -> None:
    pytest.importorskip("cutadapt")
    source = _write(tmp_path, FASTQ_TEXT)

    output = tmp_path / "limited.fastq"
    with pytest.raises(ConfigurationError, match="max_reads=2"):
        trim_fastq_adapters(source, output, ADAPTER, max_reads=2, verbose=False)
    assert not output.exists()  # the partial output is removed

    existing = _write(tmp_path, "keep me\n", name="existing.fastq")
    with pytest.raises(ConfigurationError, match="already exists"):
        trim_fastq_adapters(source, existing, ADAPTER, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep me\n"


def test_find_quality_trim_bounds_matches_hand_calculation() -> None:
    pytest.importorskip("cutadapt")
    # 'I' is Q40 and '!' is Q0. With cutoff 20 the four Q0 bases at the 3' end are removed.
    back = find_quality_trim_bounds("IIII!!!!", cutoff_back=20, verbose=False)
    assert (back["start"], back["stop"], back["retained_length"]) == (0, 4, 4)

    front = find_quality_trim_bounds("!!!!IIII", cutoff_front=20, cutoff_back=0, verbose=False)
    assert (front["start"], front["stop"]) == (4, 8)

    untouched = find_quality_trim_bounds("IIII!!!!", cutoff_back=0, verbose=False)
    assert (untouched["start"], untouched["stop"]) == (0, 8)


def test_find_quality_trim_bounds_rejects_bad_input() -> None:
    pytest.importorskip("cutadapt")
    with pytest.raises(ConfigurationError, match="Phred"):
        find_quality_trim_bounds("", verbose=False)
    with pytest.raises(ConfigurationError, match="outside the Phred"):
        find_quality_trim_bounds("II\x01", verbose=False)
    with pytest.raises(ConfigurationError, match="phred_offset"):
        find_quality_trim_bounds("IIII", phred_offset=20, verbose=False)
    with pytest.raises(ConfigurationError, match="cutoff_back"):
        find_quality_trim_bounds("IIII", cutoff_back=-1, verbose=False)


def test_calculate_expected_errors_matches_hand_calculation() -> None:
    pytest.importorskip("cutadapt")
    # '+' is Q10 (0.1 expected errors each) and 'I' is Q40 (0.0001 each).
    result = calculate_expected_errors("++II", verbose=False)
    assert result["expected_errors"] == pytest.approx(2 * 0.1 + 2 * 1e-4)
    assert result["length"] == 4


def test_missing_cutadapt_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Submodules are cached after earlier imports, so block them as well as the package.
    for name in ("cutadapt", "cutadapt.adapters", "cutadapt.modifiers", "cutadapt.qualtrim"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError, match="bio-cutadapt"):
        find_adapter_match("ACGT", ADAPTER, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-cutadapt"):
        calculate_expected_errors("IIII", verbose=False)
    source = _write(tmp_path, FASTQ_TEXT)
    with pytest.raises(BackendUnavailableError, match="bio-cutadapt"):
        trim_fastq_adapters(source, tmp_path / "x.fastq", ADAPTER, verbose=False)


def _slow5_records() -> list[dict]:
    # pA = (raw + offset) * range / digitisation
    common = {"digitisation": 8192.0, "offset": 10.0, "range": 1443.0, "sampling_rate": 4000.0}
    return [
        {"read_id": "r1", "signal": [500, 510, 520, 530], **common},
        {"read_id": "r2", "signal": [100, 200, 300], **common},
    ]


def test_write_and_summarize_slow5_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("pyslow5")
    path = tmp_path / "signals.slow5"
    written = write_slow5_signals(_slow5_records(), path, verbose=False)
    assert written["reads_written"] == 2
    assert written["signal_samples_written"] == 7

    summary = summarize_slow5_signals(path, verbose=False)
    assert summary["reads"] == 2
    assert summary["read_groups"] == 1
    assert summary["signal_samples"] == 7
    assert (summary["min_signal_samples"], summary["max_signal_samples"]) == (3, 4)
    assert summary["mean_signal_samples"] == pytest.approx(3.5)
    assert summary["sampling_rates"] == [4000.0]
    assert summary["header_sample_frequency"] == "4000"


def test_summarize_slow5_agrees_with_direct_pyslow5_scan(tmp_path: Path) -> None:
    pyslow5 = pytest.importorskip("pyslow5")
    path = tmp_path / "signals.slow5"
    write_slow5_signals(_slow5_records(), path, verbose=False)
    handle = pyslow5.Open(str(path), "r")
    direct = [len(read["signal"]) for read in handle.seq_reads()]
    handle.close()

    summary = summarize_slow5_signals(path, verbose=False)

    assert summary["reads"] == len(direct)
    assert summary["signal_samples"] == sum(direct)


def test_read_slow5_signal_returns_raw_and_picoamp_values(tmp_path: Path) -> None:
    pytest.importorskip("pyslow5")
    path = tmp_path / "signals.slow5"
    write_slow5_signals(_slow5_records(), path, verbose=False)

    raw = read_slow5_signal(path, "r2", verbose=False)
    assert raw["signal"] == [100, 200, 300]
    assert raw["units"] == "raw DAC"

    picoamps = read_slow5_signal(path, "r1", pA=True, verbose=False)
    expected = [(value + 10.0) * 1443.0 / 8192.0 for value in [500, 510, 520, 530]]
    assert picoamps["units"] == "pA"
    assert picoamps["signal"] == pytest.approx(expected, rel=1e-6)


def test_blow5_output_round_trips_and_no_index_is_written(tmp_path: Path) -> None:
    pytest.importorskip("pyslow5")
    path = tmp_path / "signals.blow5"
    write_slow5_signals(_slow5_records(), path, verbose=False)
    assert read_slow5_signal(path, "r1", verbose=False)["signal"] == [500, 510, 520, 530]
    assert sorted(item.name for item in tmp_path.iterdir()) == ["signals.blow5"]


def test_slow5_reads_reject_missing_reads_and_limits(tmp_path: Path) -> None:
    pytest.importorskip("pyslow5")
    path = tmp_path / "signals.slow5"
    write_slow5_signals(_slow5_records(), path, verbose=False)
    with pytest.raises(ConfigurationError, match="was not found"):
        read_slow5_signal(path, "missing", verbose=False)
    with pytest.raises(ConfigurationError, match="max_samples=3"):
        read_slow5_signal(path, "r1", max_samples=3, verbose=False)
    with pytest.raises(ConfigurationError, match="max_reads=1"):
        summarize_slow5_signals(path, max_reads=1, verbose=False)


def test_write_slow5_rejects_bad_records_and_keeps_existing_files(tmp_path: Path) -> None:
    pytest.importorskip("pyslow5")
    good = _slow5_records()
    with pytest.raises(ConfigurationError, match="integers in"):
        write_slow5_signals([{**good[0], "signal": [40000]}], tmp_path / "a.slow5", verbose=False)
    with pytest.raises(ConfigurationError, match="duplicate read_id"):
        write_slow5_signals([good[0], good[0]], tmp_path / "b.slow5", verbose=False)
    with pytest.raises(ConfigurationError, match="share one sampling_rate"):
        write_slow5_signals(
            [good[0], {**good[1], "sampling_rate": 5000.0}], tmp_path / "c.slow5", verbose=False
        )
    with pytest.raises(ConfigurationError, match="read_id"):
        write_slow5_signals([{**good[0], "read_id": "a b"}], tmp_path / "d.slow5", verbose=False)
    with pytest.raises(ConfigurationError, match=r"\.slow5 or \.blow5"):
        write_slow5_signals(good, tmp_path / "e.txt", verbose=False)
    assert sorted(item.name for item in tmp_path.iterdir()) == []

    existing = tmp_path / "keep.slow5"
    existing.write_text("keep\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        write_slow5_signals(good, existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep\n"


def test_malformed_slow5_is_rejected_before_pyslow5_runs(tmp_path: Path) -> None:
    # pyslow5 1.5.0 terminates the interpreter on malformed headers, so the wrapper checks first.
    pytest.importorskip("pyslow5")
    text = tmp_path / "text.slow5"
    text.write_text("hello\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="does not start with a SLOW5 or BLOW5"):
        summarize_slow5_signals(text, verbose=False)
    empty = tmp_path / "empty.slow5"
    empty.write_bytes(b"")
    with pytest.raises(ConfigurationError, match="does not start with a SLOW5 or BLOW5"):
        read_slow5_signal(empty, "r1", verbose=False)
    wrong_suffix = tmp_path / "signals.txt"
    wrong_suffix.write_text("#slow5_version\t0.2.0\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"\.slow5 or \.blow5 file"):
        summarize_slow5_signals(wrong_suffix, verbose=False)


def test_missing_pyslow5_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "signals.slow5"
    monkeypatch.setitem(sys.modules, "pyslow5", None)
    with pytest.raises(BackendUnavailableError, match="bio-pyslow5"):
        summarize_slow5_signals(path, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-pyslow5"):
        write_slow5_signals(_slow5_records(), path, verbose=False)


UMI_COUNTS = {"AAAA": 10, "AAAT": 3, "CCCC": 5, "AAAC": 1}
READ = "ACGTTTAGCCATTGGGAAA"  # cell ACGT, UMI TTA, discarded GC, insert CATTGGGAAA
READ_PATTERN = r"^(?P<cell_1>.{4})(?P<umi_1>.{3})(?P<discard_1>.{2})"
WHITELIST = ["AAAA", "CCCC", "GGGG"]
OBSERVED = ["AAAT", "AAAG", "CCCA", "GGGG", "TTTT"]


def test_cluster_umi_counts_matches_hand_calculation() -> None:
    pytest.importorskip("umi_tools")
    # AAAA (10) joins AAAT (3) and AAAC (1): 10 >= 2*3-1 and 10 >= 2*1-1. CCCC is distant.
    result = cluster_umi_counts(UMI_COUNTS, verbose=False)
    assert result["clusters"] == [["AAAA", "AAAT", "AAAC"], ["CCCC"]]
    assert result["representatives"] == ["AAAA", "CCCC"]
    assert result["cluster_counts"] == [14, 5]


def test_directional_rule_separates_similar_counts_that_adjacency_merges() -> None:
    pytest.importorskip("umi_tools")
    # 3 >= 2*3-1 is false, so the directional method does not join the two UMIs.
    counts = {"AAAA": 3, "AAAT": 3}
    directional = cluster_umi_counts(counts, method="directional", verbose=False)
    adjacency = cluster_umi_counts(counts, method="adjacency", verbose=False)
    assert sorted(directional["clusters"]) == [["AAAA"], ["AAAT"]]
    assert sorted(adjacency["clusters"][0]) == ["AAAA", "AAAT"]


def test_cluster_umi_counts_matches_direct_umi_tools_call() -> None:
    umi_tools = pytest.importorskip("umi_tools")
    counts = {b"AAAA": 10, b"AAAT": 3, b"CCCC": 5, b"AAAC": 1}
    expected = umi_tools.UMIClusterer(cluster_method="directional")(counts, 1)
    result = cluster_umi_counts(UMI_COUNTS, verbose=False)
    assert result["clusters"] == [[item.decode("ascii") for item in group] for group in expected]


def test_unique_method_keeps_every_umi_separate() -> None:
    pytest.importorskip("umi_tools")
    result = cluster_umi_counts(UMI_COUNTS, method="unique", verbose=False)
    assert sorted(result["representatives"]) == ["AAAA", "AAAC", "AAAT", "CCCC"]
    assert all(len(group) == 1 for group in result["clusters"])


def test_cluster_umi_counts_rejects_bad_input() -> None:
    pytest.importorskip("umi_tools")
    with pytest.raises(ConfigurationError, match="same length"):
        cluster_umi_counts({"AAAA": 1, "AAA": 2}, verbose=False)
    with pytest.raises(ConfigurationError, match="positive integer"):
        cluster_umi_counts({"AAAA": True}, verbose=False)
    with pytest.raises(ConfigurationError, match="only ACGTN"):
        cluster_umi_counts({"AAXA": 1}, verbose=False)
    with pytest.raises(ConfigurationError, match="method must be one of"):
        cluster_umi_counts(UMI_COUNTS, method="percentile", verbose=False)
    with pytest.raises(ConfigurationError, match="max_umis"):
        cluster_umi_counts(UMI_COUNTS, max_umis=3, verbose=False)
    with pytest.raises(ConfigurationError, match="max_distance"):
        cluster_umi_counts(UMI_COUNTS, max_distance=True, verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty mapping"):
        cluster_umi_counts({}, verbose=False)


def test_extract_read_barcodes_matches_hand_calculation() -> None:
    pytest.importorskip("umi_tools")
    quality = "I" * len(READ)
    result = extract_read_barcodes(READ, quality, READ_PATTERN, verbose=False)
    assert result["matched"] is True
    assert result["cell_barcode"] == "ACGT"
    assert result["umi"] == "TTA"
    assert result["umi_qualities"] == "III"
    assert result["read_sequence"] == "CATTGGGAAA"
    assert result["read_qualities"] == "I" * 10

    kept = extract_read_barcodes(READ, quality, READ_PATTERN, discard=False, verbose=False)
    assert kept["read_sequence"] == "GCCATTGGGAAA"  # discard bases stay when discard=False


def test_extract_read_barcodes_matches_direct_umi_tools_call() -> None:
    extract_methods = pytest.importorskip("umi_tools.extract_methods")
    umi_methods = pytest.importorskip("umi_tools.umi_methods")
    regex = pytest.importorskip("regex")
    match = regex.compile(READ_PATTERN).match(READ)
    record = umi_methods.Record("read", READ, "I" * len(READ))
    cell, _cell_quals, umi, _umi_quals, new_seq, _new_quals = extract_methods.ExtractBarcodes(
        record, match, extract_cell=True, extract_umi=True, discard=True
    )

    result = extract_read_barcodes(READ, "I" * len(READ), READ_PATTERN, verbose=False)

    assert (result["cell_barcode"], result["umi"], result["read_sequence"]) == (cell, umi, new_seq)


def test_extract_read_barcodes_reports_misses_and_rejects_bad_patterns() -> None:
    pytest.importorskip("umi_tools")
    missed = extract_read_barcodes("GG", "II", READ_PATTERN, verbose=False)
    assert missed["matched"] is False
    assert (missed["cell_barcode"], missed["umi"]) == ("", "")
    assert missed["read_sequence"] == "GG"

    with pytest.raises(ConfigurationError, match="cell_<n> or umi_<n>"):
        extract_read_barcodes(READ, "I" * len(READ), r"^(?P<discard_1>.{2})", verbose=False)
    with pytest.raises(ConfigurationError, match="not a valid regular expression"):
        extract_read_barcodes(READ, "I" * len(READ), r"^(?P<cell_1>.{4}", verbose=False)
    # The UMI group is optional, so it does not take part in a four-base read.
    with pytest.raises(ConfigurationError, match="participate in every match"):
        extract_read_barcodes("ACGT", "IIII", r"^(?P<cell_1>.{4})(?P<umi_1>.{3})?", verbose=False)


def test_map_barcodes_to_whitelist_matches_hand_calculation() -> None:
    pytest.importorskip("umi_tools")
    result = map_barcodes_to_whitelist(OBSERVED, WHITELIST, verbose=False)
    # AAAT and AAAG are one mismatch from AAAA, CCCA is one from CCCC, GGGG is already
    # whitelisted, and TTTT has no candidate.
    assert result["corrected"] == {"AAAA": ["AAAG", "AAAT"], "CCCC": ["CCCA"]}
    assert result["corrected_barcodes"] == 3
    assert result["exact_matches"] == 1
    assert result["unassigned_barcodes"] == ["TTTT"]


def test_map_barcodes_leaves_ambiguous_barcodes_unassigned() -> None:
    pytest.importorskip("umi_tools")
    # AAAT is one mismatch from both AAAA and AAAC, so it is not assigned.
    result = map_barcodes_to_whitelist(["AAAT"], ["AAAA", "AAAC"], verbose=False)
    assert result["corrected"] == {}
    assert result["unassigned_barcodes"] == ["AAAT"]


def test_map_barcodes_matches_direct_umi_tools_mapping() -> None:
    whitelist_methods = pytest.importorskip("umi_tools.whitelist_methods")
    expected = whitelist_methods.getErrorCorrectMapping(
        sorted(set(OBSERVED)), sorted(set(WHITELIST)), threshold=1
    )
    result = map_barcodes_to_whitelist(OBSERVED, WHITELIST, verbose=False)
    assert result["corrected"] == {str(key): sorted(value) for key, value in expected.items()}


def test_map_barcodes_rejects_bad_input() -> None:
    pytest.importorskip("umi_tools")
    with pytest.raises(ConfigurationError, match="threshold must be an integer"):
        map_barcodes_to_whitelist(OBSERVED, WHITELIST, threshold=0, verbose=False)
    with pytest.raises(ConfigurationError, match="same length"):
        map_barcodes_to_whitelist(["AAAA"], ["AAA"], verbose=False)
    with pytest.raises(ConfigurationError, match="only ACGT"):
        map_barcodes_to_whitelist(["AANA"], WHITELIST, verbose=False)


def test_missing_umi_tools_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in (
        "umi_tools",
        "umi_tools.extract_methods",
        "umi_tools.umi_methods",
        "umi_tools.whitelist_methods",
    ):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError, match="bio-umi-tools"):
        cluster_umi_counts(UMI_COUNTS, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-umi-tools"):
        extract_read_barcodes(READ, "I" * len(READ), READ_PATTERN, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-umi-tools"):
        map_barcodes_to_whitelist(OBSERVED, WHITELIST, verbose=False)
