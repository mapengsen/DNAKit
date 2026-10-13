"""Upstream checks behind the G18d review for kipoi, genomic-benchmarks and gReLU.

Every wrapped gReLU function is tested against the real upstream call, or against a hand
calculation. Models are tiny randomly initialized LightningModel objects built on CPU; no
test downloads weights, genomes or model-zoo sources. The kipoi listing functions are
replaced with mocks that fail if called, because the first call clones a 202 MB repository.
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import numpy as np
import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.predictions import (
    genomic_benchmarks_review_functions,
    kipoi_review_functions,
)
from dnakit.predictions import (
    grelu_review_functions as grelu,
)


def _random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _conv_model(seed: int = 0):
    lightning = pytest.importorskip("grelu.lightning")
    torch = pytest.importorskip("torch")
    torch.manual_seed(seed)
    model = lightning.LightningModel(
        model_params={
            "model_type": "ConvModel",
            "n_tasks": 2,
            "stem_channels": 8,
            "stem_kernel_size": 5,
            "n_conv": 1,
            "channel_init": 8,
            "kernel_size": 3,
            "final_pool_func": "avg",
        },
        train_params={"task": "regression", "loss": "mse"},
    )
    model.eval()
    model.data_params["tasks"] = {"name": ["t0", "t1"]}
    model.data_params["train"] = {"seq_len": 60}
    return model


def _transformer_model(seed: int = 0):
    lightning = pytest.importorskip("grelu.lightning")
    torch = pytest.importorskip("torch")
    torch.manual_seed(seed)
    model = lightning.LightningModel(
        model_params={
            "model_type": "ConvTransformerModel",
            "n_tasks": 2,
            "stem_channels": 8,
            "stem_kernel_size": 5,
            "n_conv": 1,
            "channel_init": 8,
            "kernel_size": 3,
            "n_transformers": 1,
            "n_heads": 2,
            "key_len": 4,
        },
        train_params={"task": "regression", "loss": "mse"},
    )
    model.eval()
    return model


def _write_genome(folder: Path, sequence: str) -> Path:
    fasta = folder / "genome.fa"
    fasta.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
    # Index line: name, length, byte offset of the first base, bases per line, bytes per line.
    (folder / "genome.fa.fai").write_text(
        f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n", encoding="utf-8"
    )
    (folder / "genome.fa.sizes").write_text(f"chr1\t{len(sequence)}\n", encoding="utf-8")
    return fasta


# ---------------------------------------------------------------------------
# Module contracts
# ---------------------------------------------------------------------------


def test_kipoi_and_genomic_benchmarks_review_modules_export_no_wrappers() -> None:
    assert list(kipoi_review_functions.__all__) == []
    assert list(genomic_benchmarks_review_functions.__all__) == []


def test_grelu_review_module_exports_twenty_wrappers_with_keyword_verbose() -> None:
    import inspect

    assert len(grelu.__all__) == 20
    assert len(set(grelu.__all__)) == len(grelu.__all__)
    for name in grelu.__all__:
        function = getattr(grelu, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_kipoi_listing_functions_are_mocked_and_never_called(monkeypatch) -> None:
    kipoi = pytest.importorskip("kipoi")
    calls: list[str] = []

    def refuse(name: str):
        def _fail(*args, **kwargs):
            calls.append(name)
            raise AssertionError(f"{name} must not run: it clones the model-zoo repository")

        return _fail

    monkeypatch.setattr(kipoi, "list_sources", refuse("kipoi.list_sources"))
    monkeypatch.setattr(kipoi.config, "list_sources", refuse("kipoi.config.list_sources"))
    # The review module must not touch the mocked listing functions.
    assert list(kipoi_review_functions.__all__) == []
    assert calls == []


# ---------------------------------------------------------------------------
# Intervals, chromosomes and coverage tables (no model)
# ---------------------------------------------------------------------------


def test_read_narrowpeak_intervals_reads_ten_columns(tmp_path: Path) -> None:
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text(
        "chr1\t100\t200\tp1\t500\t.\t5.5\t3.2\t1.1\t50\n"
        "chr2\t300\t420\tp2\t80\t+\t2.0\t-1\t-1\t60\n",
        encoding="utf-8",
    )
    result = grelu.read_narrowpeak_intervals(peaks, verbose=False)
    assert result["count"] == 2
    first = result["peaks"][0]
    # Hand check: summit 50 bases after start 100 gives the 0-based summit position.
    assert (first["chrom"], first["start"], first["end"]) == ("chr1", 100, 200)
    assert first["summit"] == 50 and first["signal"] == 5.5 and first["strand"] == "."
    assert result["peaks"][1]["pvalue"] == -1


def test_read_narrowpeak_intervals_rejects_short_rows(tmp_path: Path) -> None:
    short = tmp_path / "short.narrowPeak"
    short.write_text("chr1\t100\t200\tp1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        grelu.read_narrowpeak_intervals(short, verbose=False)


def test_expand_chromosome_names_counts_match_hand_count() -> None:
    assert grelu.expand_chromosome_names("autosomes", verbose=False)["count"] == 22
    assert grelu.expand_chromosome_names("autosomesX", verbose=False)["count"] == 23
    xy = grelu.expand_chromosome_names("autosomesXY", verbose=False)
    assert xy["count"] == 24 and xy["chromosomes"][-2:] == ["chrX", "chrY"]
    explicit = grelu.expand_chromosome_names(["chrM", "chr1"], verbose=False)
    assert explicit["chromosomes"] == ["chrM", "chr1"] and explicit["shortcut"] is None
    with pytest.raises(ConfigurationError):
        grelu.expand_chromosome_names("autosomesW", verbose=False)


def test_check_genomic_interval_table_accepts_layout_and_rejects_floats() -> None:
    good = [{"chrom": "chr1", "start": 1, "end": 5}]
    assert grelu.check_genomic_interval_table(good, verbose=False)["valid"] is True
    float_start = [{"chrom": "chr1", "start": 1.5, "end": 5}]
    assert grelu.check_genomic_interval_table(float_start, verbose=False)["valid"] is False
    wrong_order = [{"start": 1, "chrom": "chr1", "end": 5}]
    assert grelu.check_genomic_interval_table(wrong_order, verbose=False)["valid"] is False


def test_filter_interval_coverage_cutoff_and_top_follow_hand_calculation() -> None:
    intervals = [
        {"chrom": "chr1", "start": 0, "end": 5},
        {"chrom": "chr1", "start": 5, "end": 9},
        {"chrom": "chr2", "start": 1, "end": 4},
    ]
    coverage = [[0.5, 3.0, 1.0]]
    cutoff = grelu.filter_interval_coverage(coverage, intervals, cutoff=1.0, verbose=False)
    assert [item["index"] for item in cutoff["intervals"]] == [1, 2]
    top = grelu.filter_interval_coverage(coverage, intervals, method="top", cutoff=1, verbose=False)
    assert [item["index"] for item in top["intervals"]] == [1]
    # Mean over two tasks: [[0, 2], [2, 0]] averages to [1, 1]; both pass a cutoff of 1.
    averaged = grelu.filter_interval_coverage(
        [[0.0, 2.0], [2.0, 0.0]], intervals[:2], cutoff=1.0, verbose=False
    )
    assert [item["index"] for item in averaged["intervals"]] == [0, 1]
    none = grelu.filter_interval_coverage([[0.5, 0.2, 0.1]], intervals, cutoff=5.0, verbose=False)
    assert none["count"] == 0


def test_add_negative_intervals_labels_classes_and_zeroes_negative_columns() -> None:
    positives = [{"chrom": "chr1", "start": 0, "end": 5}, {"chrom": "chr1", "start": 5, "end": 9}]
    negatives = [{"chrom": "chr3", "start": 0, "end": 7}]
    result = grelu.add_negative_intervals([[1.0, 2.0]], positives, negatives, verbose=False)
    assert result["count"] == 3 and result["negative_count"] == 1
    assert [item["class"] for item in result["intervals"]] == ["positive", "positive", "negative"]
    assert result["coverage"][0][2] == 0.0
    assert result["coverage"][0][:2] == [1.0, 2.0]


# ---------------------------------------------------------------------------
# Motif sites and PWM scores (offline numeric checks)
# ---------------------------------------------------------------------------


def test_score_motif_sites_averages_attribution_over_hit() -> None:
    sequences = {"s0": "ACGTAC"}
    attributions = {
        "s0": [[1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    }
    hits = [{"motif": "m1", "sequence": "s0", "start": 0, "end": 4, "strand": "+"}]
    result = grelu.score_motif_sites(hits, sequences, attributions, verbose=False)
    # Hand check: the hit covers four bases with attribution 1, 1, 1, 1 in the A row only;
    # the mean over all four bases and four rows is (4 * 1) / 16 = 0.25.
    assert result["hits"][0]["site_attr_score"] == pytest.approx(0.25)


def test_score_motif_pwm_products_matches_hand_sum() -> None:
    sequences = ["ACGT"]
    # PWM with A probability 1 at both positions; attributions give A = [2, 3] and 0 elsewhere.
    pwm = [[1, 1], [0, 0], [0, 0], [0, 0]]
    attributions = [[[2, 3, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]
    hits = [{"motif": "m1", "sequence": "seq0", "start": 0, "end": 2, "strand": "+"}]
    result = grelu.score_motif_pwm_products(
        hits, sequences, attributions, {"m1": pwm}, verbose=False
    )
    # Hand check: sum over bases of attribution * PWM is [2, 3]; its mean is 2.5.
    assert result["hits"][0]["motif_attr_score"] == pytest.approx(2.5)


def test_read_tomtom_matches_applies_q_value_threshold(tmp_path: Path) -> None:
    (tmp_path / "tomtom.tsv").write_text(
        "Target_ID\tOptimal_offset\tp-value\tE-value\tq-value\tOverlap\t"
        "Query_consensus\tTarget_consensus\tOrientation\n"
        "MotifA\t0\t0.001\t0.01\t0.02\t8\tACGTACGT\tACGTACGA\t+\n"
        "MotifB\t1\t0.2\t0.5\t0.4\t7\tACGTACGT\tTTTTACGT\t-\n",
        encoding="utf-8",
    )
    result = grelu.read_tomtom_matches(tmp_path, verbose=False)
    assert result["count"] == 1
    assert result["matches"][0]["Target_ID"] == "MotifA"


def test_read_tomtom_matches_rejects_folder_without_tomtom_file(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        grelu.read_tomtom_matches(tmp_path, verbose=False)


def test_read_modisco_motifs_reads_pattern_from_h5(tmp_path: Path) -> None:
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "modisco_report.h5"
    with h5py.File(path, "w") as handle:
        group = handle.create_group("pos_patterns").create_group("pattern_0")
        # Counts are all A, so the normalized PPM is a row of ones in the A row.
        group.create_dataset("sequence", data=np.tile([4.0, 0.0, 0.0, 0.0], (6, 1)))
        group.create_dataset("contrib_scores", data=np.ones((6, 4)))
    result = grelu.read_modisco_motifs(path, verbose=False)
    assert result["count"] == 1
    motif = result["motifs"][0]
    assert motif["name"] == "pos_pattern_0"
    ppm = np.asarray(motif["ppm"])
    assert ppm.shape[0] == 4
    assert np.allclose(ppm[0], 1.0) and np.allclose(ppm[1:], 0.0)


# ---------------------------------------------------------------------------
# Model-based experiments (tiny random models, CPU)
# ---------------------------------------------------------------------------


def test_get_sequence_attributions_shape_and_global_rng_restored() -> None:
    model = _conv_model()
    sequences = [_random_dna(60, 1), _random_dna(60, 2)]
    np.random.seed(123)
    before = np.random.get_state()[1].copy()
    result = grelu.get_sequence_attributions(
        model, sequences, method="deepshap", n_shuffles=2, seed=0, verbose=False
    )
    after = np.random.get_state()[1]
    assert np.array_equal(before, after)
    assert result["count"] == 2 and result["length"] == 60
    assert np.asarray(result["attributions"]).shape == (2, 4, 60)


def test_get_sequence_attributions_rejects_unequal_lengths_and_bad_task() -> None:
    model = _conv_model()
    with pytest.raises(ConfigurationError):
        grelu.get_sequence_attributions(model, ["ACGT", "ACGTAC"], verbose=False)
    with pytest.raises(ConfigurationError):
        grelu.get_sequence_attributions(model, _random_dna(60, 3), task_index=5, verbose=False)


def test_in_silico_mutagenesis_reference_base_scores_zero_under_subtract() -> None:
    model = _conv_model()
    sequence = _random_dna(60, 4)
    result = grelu.in_silico_mutagenesis_scores(
        model, sequence, compare_func="subtract", start=0, end=5, verbose=False
    )
    assert result["positions"] == 5
    for position, row in enumerate(result["scores"]):
        reference = "ACGT".index(sequence[position])
        # Subtract compares each mutant with the reference, so the reference base is zero.
        assert row[reference] == pytest.approx(0.0, abs=1e-6)


def test_get_attention_scores_rows_are_probabilities() -> None:
    model = _transformer_model()
    sequence = _random_dna(60, 5)
    result = grelu.get_attention_scores(model, sequence, block_idx=0, verbose=False)
    attention = np.asarray(result["attention"])
    assert result["shape"] == [2, 60, 60]
    # Hand check: each attention row is a softmax, so it sums to one.
    assert np.allclose(attention.sum(axis=-1), 1.0, atol=1e-4)


def test_marginalize_pattern_insertions_returns_effect_per_shuffle() -> None:
    model = _conv_model()
    result = grelu.marginalize_pattern_insertions(
        model,
        [_random_dna(60, 6)],
        ["GATTACA"],
        compare_func="subtract",
        n_shuffles=2,
        verbose=False,
    )
    assert result["shape"] == [1, 2, 1, 1, 1]
    assert len(result["effects"]) == 1


def test_marginalize_pattern_insertions_rejects_zero_shuffles() -> None:
    model = _conv_model()
    with pytest.raises(ConfigurationError):
        grelu.marginalize_pattern_insertions(
            model, [_random_dna(60, 6)], ["GATTACA"], n_shuffles=0, verbose=False
        )


def test_marginalize_pattern_spacing_returns_one_row_per_spacing() -> None:
    model = _conv_model()
    result = grelu.marginalize_pattern_spacing(
        model,
        _random_dna(60, 7),
        "GATTACA",
        "CCC",
        stride=20,
        n_shuffles=1,
        compare_func="subtract",
        verbose=False,
    )
    assert len(result["spacings"]) == 3
    assert result["shape"][2] == 3


def test_shuffle_tile_effects_tile_table_matches_hand_count() -> None:
    model = _conv_model()
    result = grelu.shuffle_tile_effects(
        model, [_random_dna(60, 8)], 10, stride=10, n_shuffles=1, verbose=False
    )
    # Hand check: 60 bases in tiles of 10 with stride 10 give six tiles starting at 0..50.
    assert result["tile_count"] == 6
    assert [tile["start"] for tile in result["tiles"]] == [0, 10, 20, 30, 40, 50]


def test_grelu_predict_variant_effects_checks_reference_and_identity_ratio(tmp_path: Path) -> None:
    model = _conv_model()
    sequence = _random_dna(300, 9)
    fasta = _write_genome(tmp_path, sequence)
    base = sequence[99]
    result = grelu.grelu_predict_variant_effects(
        model,
        [{"chrom": "chr1", "pos": 100, "ref": base, "alt": base}],
        fasta,
        seq_len=60,
        compare_func="divide",
        verbose=False,
    )
    # Hand check: an alternate allele equal to the reference gives identical predictions,
    # so the divide ratio is exactly one.
    assert result["effects"][0][0] == pytest.approx(1.0)
    wrong = "A" if base != "A" else "C"
    with pytest.raises(ConfigurationError):
        grelu.grelu_predict_variant_effects(
            model,
            [{"chrom": "chr1", "pos": 100, "ref": wrong, "alt": base}],
            fasta,
            seq_len=60,
            verbose=False,
        )


def test_grelu_marginalize_variants_returns_test_table(tmp_path: Path) -> None:
    model = _conv_model()
    sequence = _random_dna(300, 10)
    fasta = _write_genome(tmp_path, sequence)
    base = sequence[99]
    result = grelu.grelu_marginalize_variants(
        model,
        [{"chrom": "chr1", "pos": 100, "ref": base, "alt": "A" if base != "A" else "C"}],
        fasta,
        seq_len=60,
        n_shuffles=2,
        seed=0,
        verbose=False,
    )
    assert result["count"] == 1
    assert result["result"]


def test_design_sequences_evolve_changes_only_allowed_positions() -> None:
    model = _conv_model()
    sequence = _random_dna(30, 11)
    result = grelu.design_sequences_evolve(
        [sequence], model, max_iter=2, positions=[0, 1, 2], verbose=False
    )
    assert result["count"] > 0
    for record in result["records"]:
        designed = record.get("seq")
        if designed:
            changed = [
                index
                for index, (a, b) in enumerate(zip(sequence, designed, strict=False))
                if a != b
            ]
            assert set(changed) <= {0, 1, 2}


def test_design_sequences_evolve_requires_task_names() -> None:
    model = _conv_model()
    model.data_params.pop("tasks")
    with pytest.raises(ConfigurationError):
        grelu.design_sequences_evolve([_random_dna(30, 12)], model, max_iter=1, verbose=False)


def test_design_sequence_ledidi_edits_only_mutable_positions() -> None:
    model = _conv_model()
    sequence = _random_dna(60, 13)
    result = grelu.design_sequence_ledidi(
        sequence,
        model,
        max_iter=30,
        positions=[0, 1, 2, 3, 4],
        edit_penalty=0.0,
        verbose=False,
    )
    assert result["count"] == 16
    for design in result["designs"]:
        assert len(design) == len(sequence)
        changed = [
            index for index, (a, b) in enumerate(zip(sequence, design, strict=True)) if a != b
        ]
        assert set(changed) <= {0, 1, 2, 3, 4}


def test_run_modisco_motifs_reports_no_patterns_for_random_sequences(tmp_path: Path) -> None:
    model = _conv_model()
    out_dir = tmp_path / "modisco_out"
    sequences = [_random_dna(60, 100 + index) for index in range(40)]
    result = grelu.run_modisco_motifs(
        model, sequences, out_dir=out_dir, window=60, n_shuffles=2, seed=0, verbose=False
    )
    assert result["patterns_found"] is False
    assert result["files"] == []
    assert not out_dir.exists()


def test_run_modisco_motifs_refuses_non_empty_output_folder(tmp_path: Path) -> None:
    model = _conv_model()
    (tmp_path / "keep.txt").write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        grelu.run_modisco_motifs(
            model,
            [_random_dna(60, 200)],
            out_dir=tmp_path,
            verbose=False,
        )


def test_run_modisco_motifs_leaves_working_directory_unchanged(tmp_path: Path) -> None:
    cwd = os.getcwd()
    model = _conv_model()
    out_dir = tmp_path / "modisco_cwd"
    grelu.run_modisco_motifs(
        model,
        [_random_dna(60, 300 + index) for index in range(40)],
        out_dir=out_dir,
        n_shuffles=2,
        verbose=False,
    )
    assert os.getcwd() == cwd
