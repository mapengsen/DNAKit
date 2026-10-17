"""Regression tests for the fixes in sequencing_qc, singlecell_spatial, transcriptomics, ml_models, visualization and external_adapters."""

from __future__ import annotations

import pytest

import dnakit.external_adapters as ea
import dnakit.ml_models as ml
import dnakit.sequencing_qc as qc
import dnakit.singlecell_spatial as sc
import dnakit.visualization._candidates as viz
from dnakit.exceptions import InputFormatError


def ok(envelope: dict) -> dict:
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# ---------------------------------------------------------------- sequencing_qc


def test_calculate_quality_percentiles_have_no_float_noise() -> None:
    reads = [("r1", "ACGT", "5I?+"), ("r2", "ACGT", "I5+?"), ("r3", "ACGT", "?+I5"), ("r4", "ACGT", "+?5I")]
    result = ok(qc.calculate_quality(reads))
    values = [result["per_read_mean_quality"]["median"]]
    for row in result["per_position"]:
        values += [row[key] for key in ("mean", "median", "q1", "q3", "p10", "p90")]
    assert all(v == round(v, 6) for v in values)


def test_percentile_rounding_removes_noise_in_sequencing_qc_output() -> None:
    # 14 reads with Phred 6..19: the 90th percentile position is 11.7 -> 17.700000000000003 without rounding
    reads = [(f"r{i}", "A", chr(33 + 6 + i)) for i in range(14)]
    row = ok(qc.calculate_quality(reads))["per_position"][0]
    assert row["p90"] == 17.7
    assert ok(qc.calculate_yield(reads))["reads"] == 14


@pytest.mark.parametrize(
    "call",
    [
        lambda pair: qc.assign_sample_barcodes(pair, options={"samples": {"s": "ACG"}}),
        lambda pair: qc.calculate_barcode_quality(pair),
        lambda pair: qc.estimate_read_error_model(pair),
        lambda pair: qc.extract_umi_sequences(pair),
        lambda pair: qc.calculate_index_hopping(pair, options={"samples": {"s": "A+C"}}),
        lambda pair: qc.collapse_umi_reads(pair),
        lambda pair: qc.demultiplex_barcoded_reads(pair, options={"samples": {"s": "ACG"}}),
    ],
)
def test_functions_with_single_read_set_input_document_and_reject_pairs(call) -> None:
    pair = {"read1": [("a", "ACGTACGT", "IIIIIIII")], "read2": [("a", "ACGTACGT", "IIIIIIII")]}
    with pytest.raises(InputFormatError) as caught:
        call(pair)
    assert caught.value.code == "INVALID_INPUT"


@pytest.mark.parametrize(
    "function",
    [qc.assign_sample_barcodes, qc.calculate_barcode_quality, qc.estimate_read_error_model, qc.extract_umi_sequences,
     qc.calculate_index_hopping, qc.calculate_index_hopping_tenx, qc.collapse_umi_reads, qc.demultiplex_barcoded_reads],
)
def test_docstrings_state_that_a_pair_is_not_accepted(function) -> None:
    assert "not accepted" in function.__doc__ or "not a {'read1'" in function.__doc__


def test_detect_platform_tenx_docstring_states_length_only() -> None:
    assert "Only the read-1 length is used" in qc.detect_platform_tenx.__doc__


# ---------------------------------------------------------------- singlecell_spatial


def test_build_cell_peak_matrix_reports_n_peaks_not_n_genes() -> None:
    fragments = "chr1\t110\t120\tAAA\nchr1\t160\t170\tBBB\n"
    result = ok(sc.build_cell_peak_matrix(fragments, options={"peaks": ["chr1:100-200", "chr1:150-300", "chr2:0-50"]}))
    assert result["n_peaks"] == 3 and result["n_cells"] == 2
    assert "n_genes" not in result


def test_extract_barcodes_separates_extracted_unmatched_and_too_short() -> None:
    reads = ["AAAACCCC", "AAAACCCA", "GGGGGGGG", "TTTTTTTT", "AA"]
    result = ok(sc.extract_barcodes(reads, options={"barcode_length": 4, "whitelist": ["AAAA", "TTTT"]}))
    assert result["n_reads"] == 5
    assert result["n_extracted"] == 3  # AAAA, AAAA, TTTT
    assert result["n_unmatched"] == 1  # GGGG
    assert result["n_too_short"] == 1
    assert result["n_reads"] == result["n_extracted"] + result["n_unmatched"] + result["n_too_short"]
    assert result["whitelist"]["unmatched"] == result["n_unmatched"]


def test_extract_barcodes_without_whitelist_has_zero_unmatched() -> None:
    result = ok(sc.extract_barcodes(["AAAACCCC", "GGGGTTTT"], options={"barcode_length": 4}))
    assert (result["n_reads"], result["n_extracted"], result["n_unmatched"], result["n_too_short"]) == (2, 2, 0, 0)
    spatial = ok(sc.extract_barcodes_spatial(["AAAACCCC", "GGGGTTTT"], options={"barcode_length": 4, "whitelist": ["AAAA"]}))
    assert (spatial["n_extracted"], spatial["n_unmatched"]) == (1, 1)


@pytest.mark.parametrize("function", [sc.remove_doublets, sc.remove_doublets_spatial])
def test_remove_doublets_does_not_suggest_installing_scrublet(function) -> None:
    envelope = function({"g": [1, 2]})
    assert envelope["status"] == "not_implemented"
    assert "pip install" not in envelope["required_dependency"]
    assert "DNAKit does not call" in envelope["required_dependency"]
    assert "deconvolution" not in envelope["message"]
    assert "DNAKit does not call" in envelope["message"]


def test_differential_docstrings_describe_group1_versus_group2() -> None:
    for function in (sc.find_differential_features, sc.find_differential_features_spatial):
        assert "log2(mean(value, group1)" in function.__doc__ and "GROUP_TOO_SMALL" in function.__doc__
    assert "versus control" in sc.calculate_perturbation_effect.__doc__ or "control) + pseudocount" in sc.calculate_perturbation_effect.__doc__
    assert "skipped with a warning" in sc.calculate_markers.__doc__


# ---------------------------------------------------------------- ml_models


def test_validate_message_uses_acgtn_order() -> None:
    out = ok(ml.validate([{"sequence": "ACGX", "label": 0}]))
    message = out["errors"][0]["message"]
    assert message.endswith("outside the alphabet ACGTN"), message
    custom = ok(ml.validate([{"sequence": "ACGX", "label": 0}], options={"alphabet": "tgca"}))
    assert custom["errors"][0]["message"].endswith("alphabet TGCA")


def test_calibrate_ece_before_is_null_outside_unit_interval() -> None:
    out = ok(ml.calibrate({"y_true": [0, 1, 0, 1], "y_score": [-2.0, 1.0, -1.0, 3.0]}))
    assert out["ece_before"] is None and out["ece_after_on_fit_data"] is not None
    inside = ok(ml.calibrate({"y_true": [0, 1, 0, 1], "y_score": [0.1, 0.9, 0.4, 0.6]}))
    assert inside["ece_before"] is not None


def test_explain_accepts_scorer_as_alias_of_model() -> None:
    data = {"X": [[1.0, 0.0], [2.0, 1.0], [3.0, 0.0], [4.0, 1.0]], "y": [1.0, 2.0, 3.0, 4.0], "feature_names": ["a", "b"]}
    via_model = ok(ml.explain(data, options={"model": lambda rows: [r[0] for r in rows], "seed": 1}))
    via_scorer = ok(ml.explain(data, options={"scorer": lambda rows: [r[0] for r in rows], "seed": 1}))
    assert via_model == via_scorer
    assert "scorer" in ml.explain.__doc__


# ---------------------------------------------------------------- visualization


def test_plot_circular_map_forces_circular_whatever_the_topology() -> None:
    from dnakit.visualization import plot

    features = [{"type": "gene", "start": 1, "end": 5}]
    linear = ok(plot("ACGTACGTAC", options={"kind": "circular_map", "features": features, "topology": "linear"}))
    circular = ok(plot("ACGTACGTAC", options={"kind": "circular_map", "features": features, "topology": "circular"}))
    assert linear["svg"] == circular["svg"]
    assert "forced to topology" in viz.plot.__doc__


# ---------------------------------------------------------------- external_adapters


def argv(command, **options) -> list[str]:
    return ok(ea.build_command(command, options=options))["argv"]


def test_build_command_subcommand_precedes_flags() -> None:
    mapping = {"program": "samtools", "args": ["x.bam"], "flags": {"-b": True}}
    assert argv(mapping, subcommand="view") == ["samtools", "view", "-b", "x.bam"]
    assert argv(mapping, subcommand=["index", "build"]) == ["samtools", "index", "build", "-b", "x.bam"]


def test_build_command_without_subcommand_keeps_the_old_order() -> None:
    assert argv({"program": "samtools", "args": ["view", "x.bam"], "flags": {"-b": True}}) == ["samtools", "-b", "view", "x.bam"]


def test_build_command_subcommand_with_option_flags_and_other_forms() -> None:
    mapping = {"program": "tool", "flags": {"-a": True}, "args": ["one"]}
    assert argv(mapping, subcommand="sub", flags={"-z": 1}, args=["two"]) == ["tool", "sub", "-a", "one", "-z", "1", "two"]
    assert argv("samtools x.bam", subcommand="view") == ["samtools", "view", "x.bam"]
    assert argv(["samtools", "x.bam"], subcommand=["view"]) == ["samtools", "view", "x.bam"]


@pytest.mark.parametrize("bad", [[], "", [""], 5.5 + 0j, [None]])
def test_build_command_rejects_invalid_subcommand(bad) -> None:
    with pytest.raises(Exception) as caught:
        ea.build_command("samtools", options={"subcommand": bad})
    assert getattr(caught.value, "code", "") in {"INVALID_OPTION", "INVALID_COMMAND"}


def test_run_uses_the_same_argv_as_build_command() -> None:
    command = {"program": "echo", "args": ["tail"], "flags": {"-n": True}}
    built = argv(command, subcommand="head")
    ran = ok(ea.run(command, options={"subcommand": "head"}))
    assert ran["argv"][1:] == built[1:] == ["head", "-n", "tail"]
    assert ran["stdout"] == "head -n tail\n"


# ---------------------------------------------------------------- transcriptomics


def test_translation_efficiency_gene_lengths_only_relabel_the_unit() -> None:
    import dnakit.transcriptomics as tx

    ribo = "gene,s1,s2\ng1,100,200\ng2,300,400\n"
    rna = "gene,s1,s2\ng1,50,100\ng2,600,800\n"
    plain = ok(tx.calculate_translation_efficiency(ribo, options={"rna": rna}))
    sized = ok(tx.calculate_translation_efficiency(ribo, options={"rna": rna, "gene_lengths": {"g1": 1000, "g2": 4000}}))
    assert plain["genes"] == sized["genes"]
    assert "CPM" in plain["definition"] and "RPKM" in sized["definition"]
    assert "no effect" in tx.calculate_translation_efficiency.__doc__
