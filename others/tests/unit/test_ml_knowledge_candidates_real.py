"""Hand-verified tests for the real ml_models and knowledge_graph candidate implementations."""

from __future__ import annotations

import csv
import hashlib
import inspect
import json
import math
from pathlib import Path

import pytest

from dnakit import knowledge_graph as kg
from dnakit import ml_models as ml
from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError

V = {"verbose": False}


def run(function, value=None, **options):
    """Call a candidate function and return its envelope."""
    return function(value, options=options, verbose=False)


def result(function, value=None, **options):
    envelope = run(function, value, **options)
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# ---------------------------------------------------------------- ml_models: metrics

BINARY = {
    "y_true": [0, 0, 1, 1, 1, 0],
    "y_pred": [0, 1, 1, 1, 0, 0],
    "y_score": [0.1, 0.6, 0.8, 0.7, 0.4, 0.2],
}


def test_evaluate_binary_classification_hand_values():
    out = result(ml.evaluate, BINARY)
    m = out["metrics"]
    # confusion: TN=2 FP=1 FN=1 TP=2
    assert out["confusion_matrix"] == {"classes": ["0", "1"], "matrix": [[2, 1], [1, 2]]}
    assert m["accuracy"] == pytest.approx(4 / 6)
    assert m["precision"] == pytest.approx(2 / 3)
    assert m["recall"] == pytest.approx(2 / 3)
    assert m["f1"] == pytest.approx(2 / 3)
    assert m["specificity"] == pytest.approx(2 / 3)
    assert m["mcc"] == pytest.approx((2 * 2 - 1 * 1) / math.sqrt(3 * 3 * 3 * 3))  # = 1/3
    # AUROC: 8 of the 9 positive/negative pairs are ordered correctly
    assert m["auroc"] == pytest.approx(8 / 9)
    # AP: 1/3*1 + 1/3*1 + 1/3*(3/4)
    assert m["auprc"] == pytest.approx(1 / 3 + 1 / 3 + 0.25)
    assert m["log_loss"] == pytest.approx(
        -(math.log(0.9) + math.log(0.4) + math.log(0.8) + math.log(0.7) + math.log(0.4) + math.log(0.8)) / 6
    )


def test_evaluate_multiclass_hand_values():
    out = result(ml.evaluate, {"y_true": [0, 1, 2, 2], "y_pred": [0, 2, 2, 2]})
    m = out["metrics"]
    assert out["confusion_matrix"]["matrix"] == [[1, 0, 0], [0, 0, 1], [0, 0, 2]]
    assert m["accuracy"] == pytest.approx(0.75)
    assert m["balanced_accuracy"] == pytest.approx((1 + 0 + 1) / 3)
    assert m["f1_macro"] == pytest.approx((1.0 + 0.0 + 0.8) / 3)
    # Gorodkin: (3*4 - (1*1 + 1*0 + 2*3)) / (sqrt(16-10) * sqrt(16-6))
    assert m["mcc"] == pytest.approx(5 / math.sqrt(60))
    assert out["per_class"]["2"]["precision"] == pytest.approx(2 / 3)


def test_evaluate_scores_only_uses_threshold_and_string_labels():
    out = result(
        ml.evaluate,
        {"y_true": ["neg", "neg", "pos", "pos"], "y_score": [0.2, 0.6, 0.4, 0.9]},
        threshold=0.5,
    )
    assert out["metrics"]["positive_label"] == "pos"
    assert out["confusion_matrix"]["matrix"] == [[1, 1], [1, 1]]
    assert out["metrics"]["auroc"] == pytest.approx(0.75)


def test_evaluate_regression_matches_textbook_values():
    # the worked example of scikit-learn's mean_squared_error / r2_score documentation
    m = result(ml.evaluate, {"y_true": [3, -0.5, 2, 7], "y_pred": [2.5, 0.0, 2, 8]})["metrics"]
    assert m["mae"] == pytest.approx(0.5)
    assert m["mse"] == pytest.approx(0.375)
    assert m["rmse"] == pytest.approx(math.sqrt(0.375))
    assert m["r2"] == pytest.approx(0.9486081370449679)
    assert m["spearman"] == pytest.approx(1.0)
    assert m["pearson"] == pytest.approx(0.9848696184482703)


def test_evaluate_correlations_and_constant_target():
    out = result(ml.evaluate, {"y_true": [1.5, 2.5, 3.5, 9.5], "y_pred": [4.0, 3.0, 2.0, 1.0]}, task="regression")
    # Sxy = -12.5, Sxx = 38.75, Syy = 5 -> r = -12.5 / sqrt(38.75 * 5)
    assert out["metrics"]["pearson"] == pytest.approx(-12.5 / math.sqrt(38.75 * 5))
    assert out["metrics"]["spearman"] == pytest.approx(-1.0)
    constant = result(ml.evaluate, {"y_true": [2.0, 2.0], "y_pred": [1.0, 3.0]}, task="regression")["metrics"]
    assert constant["pearson"] is None and constant["r2"] == 0.0


def test_evaluate_matches_scikit_learn_with_ties():
    sk = pytest.importorskip("sklearn.metrics")
    y = [1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1]
    s = [0.9, 0.9, 0.8, 0.5, 0.5, 0.3, 0.3, 0.3, 0.7, 0.1, 0.5, 0.95]
    pred = [int(v >= 0.5) for v in s]
    m = result(ml.evaluate, {"y_true": y, "y_pred": pred, "y_score": s})["metrics"]
    assert m["auroc"] == pytest.approx(sk.roc_auc_score(y, s))
    assert m["auprc"] == pytest.approx(sk.average_precision_score(y, s))
    assert m["mcc"] == pytest.approx(sk.matthews_corrcoef(y, pred))
    assert m["f1"] == pytest.approx(sk.f1_score(y, pred))
    assert m["cohen_kappa"] == pytest.approx(sk.cohen_kappa_score(y, pred))
    assert m["log_loss"] == pytest.approx(sk.log_loss(y, s))
    mc_true, mc_pred = [0, 1, 2, 1, 0, 2, 2, 1, 0, 0], [0, 2, 2, 1, 1, 2, 0, 1, 0, 0]
    mc = result(ml.evaluate, {"y_true": mc_true, "y_pred": mc_pred})["metrics"]
    assert mc["mcc"] == pytest.approx(sk.matthews_corrcoef(mc_true, mc_pred))
    assert mc["f1_macro"] == pytest.approx(sk.f1_score(mc_true, mc_pred, average="macro"))
    assert mc["f1_weighted"] == pytest.approx(sk.f1_score(mc_true, mc_pred, average="weighted"))
    assert mc["balanced_accuracy"] == pytest.approx(sk.balanced_accuracy_score(mc_true, mc_pred))
    rt = [0.3, 1.2, -0.7, 5.5, 2.2, 2.2]
    rp = [0.5, 1.0, -0.2, 4.5, 2.0, 2.9]
    r = result(ml.evaluate, {"y_true": rt, "y_pred": rp})["metrics"]
    assert r["r2"] == pytest.approx(sk.r2_score(rt, rp))
    assert r["rmse"] == pytest.approx(math.sqrt(sk.mean_squared_error(rt, rp)))
    assert r["explained_variance"] == pytest.approx(sk.explained_variance_score(rt, rp))


def test_evaluate_reads_table_files_and_filters_metrics(tmp_path: Path):
    path = tmp_path / "pred.csv"
    path.write_text("label,prediction\n1,1\n0,1\n1,1\n0,0\n", encoding="utf-8")
    out = result(ml.evaluate, str(path), metrics=["accuracy", "mcc"])
    assert out["metrics"]["accuracy"] == 0.75 and set(out["metrics"]) == {"accuracy", "mcc"}
    # TP=2 FP=1 FN=0 TN=1
    expected = (2 * 1 - 1 * 0) / math.sqrt((2 + 1) * (2 + 0) * (1 + 1) * (1 + 0))
    assert out["metrics"]["mcc"] == pytest.approx(expected) == pytest.approx(1 / math.sqrt(3))


def test_check_reports_problems_without_raising():
    ok = result(ml.check, BINARY)
    assert ok["valid"] and ok["n"] == 6 and ok["warnings"] == []
    mismatch = result(ml.check, {"y_true": [0, 1, 1], "y_pred": [0, 1]})
    assert not mismatch["valid"] and mismatch["errors"][0]["code"] == "ML_LENGTH_MISMATCH"
    odd = result(ml.check, {"y_true": [0, 1], "y_pred": [0, 2], "y_score": [0.2, 1.7]})
    assert odd["valid"] and len(odd["warnings"]) == 2
    bad = result(ml.check, {"y_true": [1.5, 2.5], "y_pred": ["a", "b"]}, task="regression")
    assert not bad["valid"] and bad["errors"][0]["code"] == "NON_NUMERIC"


def test_validate_dataset_finds_errors_and_leakage():
    rows = [
        {"id": "a", "sequence": "ACGT", "label": 1, "split": "train"},
        {"id": "b", "sequence": "ACGT", "label": 1, "split": "test"},  # leakage
        {"id": "c", "sequence": "AXGT", "label": 0, "split": "train"},  # bad character
        {"id": "c", "sequence": "TTTT", "label": "", "split": "train"},  # duplicate id, missing label
        {"id": "e", "sequence": "GGGG", "label": 0, "split": "test"},
    ]
    out = result(ml.validate, rows, split_column="split")
    codes = {e["code"] for e in out["errors"]}
    assert codes == {"SPLIT_LEAKAGE", "INVALID_CHARACTERS", "DUPLICATE_ID", "MISSING_LABEL"}
    assert not out["valid"]
    assert out["summary"]["label_counts"] == {"0": 2, "1": 2}
    assert out["summary"]["duplicate_sequence_rows"] == 1
    clean = result(ml.validate, rows[:1] + rows[4:], split_column="split")
    assert clean["valid"] and clean["errors"] == []
    conflict = result(ml.validate, [{"sequence": "AC", "label": 0}, {"sequence": "ac", "label": 1}])
    assert [e["code"] for e in conflict["errors"]] == ["CONFLICTING_LABELS"]


# ---------------------------------------------------------------- ml_models: calibration, uncertainty

def test_isotonic_calibration_hand_values():
    out = result(
        ml.calibrate,
        {"y_true": [0, 0, 1, 0, 1, 1], "y_score": [0.1, 0.2, 0.35, 0.4, 0.8, 0.9]},
        method="isotonic",
        apply_to=[0.15, 0.375, 0.6, 0.0, 1.0],
    )
    # pool-adjacent-violators: [0, 0 | 1, 0 -> 0.5 | 1, 1]; linear interpolation between blocks, clipped outside
    assert out["model"]["parameters"]["blocks"] == [
        {"x_first": 0.1, "x_last": 0.2, "value": 0.0},
        {"x_first": 0.35, "x_last": 0.4, "value": 0.5},
        {"x_first": 0.8, "x_last": 0.9, "value": 1.0},
    ]
    assert out["calibrated"] == pytest.approx([0.0, 0.5, 0.5 + 0.5 * (0.6 - 0.4) / 0.4, 0.0, 1.0])
    # ECE before: |0-.1| + |0-.2| + |1-.35| + |0-.4| + |1-.8| + |1-.9| over 6 samples (one sample per bin)
    assert out["ece_before"] == pytest.approx(1.65 / 6)
    assert out["ece_after_on_fit_data"] == pytest.approx(0.0)


def test_platt_calibration_symmetric_case_has_closed_form():
    # scores -1 and +1 with one positive: smoothed targets 1/3 and 2/3 give A = -ln 2, B = 0
    out = result(ml.calibrate, {"y_true": [0, 1], "y_score": [-1.0, 1.0]}, method="platt")
    params = out["model"]["parameters"]
    assert params["A"] == pytest.approx(-math.log(2), abs=1e-6)
    assert params["B"] == pytest.approx(0.0, abs=1e-6)
    assert out["calibrated"] == pytest.approx([1 / 3, 2 / 3], abs=1e-6)


def test_platt_calibration_matches_scikit_learn():
    cal = pytest.importorskip("sklearn.calibration")
    np = pytest.importorskip("numpy")
    scores = [-2.0, -1.5, -0.4, -0.1, 0.2, 0.3, 0.9, 1.4, 2.0, -0.8]
    labels = [0, 0, 0, 1, 0, 1, 1, 1, 1, 1]
    out = result(ml.calibrate, {"y_true": labels, "y_score": scores}, method="platt")
    a, b = cal._sigmoid_calibration(np.array(scores), np.array(labels))
    assert out["model"]["parameters"]["A"] == pytest.approx(a, abs=1e-3)  # p = 1 / (1 + exp(a * f + b))
    assert out["model"]["parameters"]["B"] == pytest.approx(b, abs=1e-3)
    expected = 1 / (1 + np.exp(a * np.array(scores) + b))
    assert out["calibrated"] == pytest.approx(list(expected), abs=1e-3)


def test_export_calibration_model_writes_json(tmp_path: Path):
    target = tmp_path / "cal.json"
    out = result(
        ml.export_calibration_model, {"y_true": [0, 0, 1, 1], "y_score": [0.1, 0.4, 0.6, 0.9]},
        method="isotonic", output=str(target),
    )
    saved = json.loads(target.read_text(encoding="utf-8"))
    assert saved["type"] == "dnakit.calibration" and saved["method"] == "isotonic"
    assert out["output_path"] == str(target)
    assert out["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest()
    with pytest.raises(ConfigurationError, match="ML_OUTPUT_EXISTS"):
        run(ml.export_calibration_model, {"y_true": [0, 1], "y_score": [0.1, 0.9]}, output=str(target))


def test_prediction_uncertainty_classification_and_ensembles():
    single = result(ml.calculate_prediction_uncertainty, {"probabilities": [[0.5, 0.5], [1.0, 0.0]]})
    assert single["samples"][0]["entropy"] == pytest.approx(math.log(2))
    assert single["samples"][0]["margin"] == pytest.approx(0.0)
    assert single["samples"][1]["entropy"] == pytest.approx(0.0)
    assert single["samples"][1]["margin"] == pytest.approx(1.0)
    base2 = result(ml.calculate_prediction_uncertainty, {"probabilities": [[0.5, 0.25, 0.25]]}, log_base=2)
    assert base2["samples"][0]["entropy"] == pytest.approx(1.5)
    # two members that fully disagree: predictive entropy ln 2, expected entropy 0, mutual information ln 2
    ensemble = result(
        ml.calculate_prediction_uncertainty, {"probabilities": [[[1.0, 0.0]], [[0.0, 1.0]]]}
    )["samples"][0]
    assert ensemble["mean_probabilities"] == [0.5, 0.5]
    assert ensemble["predictive_entropy"] == pytest.approx(math.log(2))
    assert ensemble["expected_entropy"] == pytest.approx(0.0)
    assert ensemble["mutual_information"] == pytest.approx(math.log(2))
    agree = result(ml.calculate_prediction_uncertainty, {"probabilities": [[[0.5, 0.5]], [[0.5, 0.5]]]})["samples"][0]
    assert agree["mutual_information"] == pytest.approx(0.0, abs=1e-12)


def test_prediction_uncertainty_regression_ensemble():
    out = result(ml.calculate_prediction_uncertainty, {"predictions": [[1.0, 2.0], [3.0, 6.0]]})
    assert [s["mean"] for s in out["samples"]] == [2.0, 4.0]
    assert [s["variance"] for s in out["samples"]] == [1.0, 4.0]
    assert [s["std"] for s in out["samples"]] == [1.0, 2.0]
    sample = result(ml.calculate_prediction_uncertainty, {"predictions": [[1.0], [3.0]]}, ddof=1)
    assert sample["samples"][0]["variance"] == pytest.approx(2.0)


def test_conformal_uncertainty_model_quantile_rank(tmp_path: Path):
    t = [float(i) for i in range(1, 10)]
    p = [0.0] * 9  # residuals 1..9
    model = result(ml.export_conformal_uncertainty_model, {"y_true": t, "y_pred": p}, alpha=0.2)["model"]
    assert model["rank"] == 8 and model["quantile"] == 8.0  # ceil(10 * 0.8) = 8th smallest
    small = result(ml.export_conformal_uncertainty_model, {"y_true": t, "y_pred": p}, alpha=0.05)["model"]
    assert small["rank"] is None and small["quantile"] is None
    cls = result(
        ml.export_conformal_uncertainty_model,
        {"y_true": ["a", "b", "a", "b"], "probabilities": [[0.9, 0.1], [0.2, 0.8], [0.6, 0.4], [0.3, 0.7]]},
        alpha=0.5, output=str(tmp_path / "u.json"),
    )
    # nonconformity 1 - p(true) = [0.1, 0.2, 0.4, 0.3]; k = ceil(5 * 0.5) = 3 -> 0.3
    assert cls["model"]["quantile"] == pytest.approx(0.3)
    assert json.loads((tmp_path / "u.json").read_text(encoding="utf-8"))["task"] == "classification"


# ---------------------------------------------------------------- ml_models: tables and splits

def test_leaderboard_uses_competition_ranking_and_direction():
    rows = [{"model": "a", "score": 0.9}, {"model": "b", "score": 0.8}, {"model": "c", "score": 0.8}, {"model": "d", "score": 0.5}]
    out = result(ml.export_leaderboard, rows, format="json")
    assert [r["rank"] for r in out["ranking"]] == [1, 2, 2, 4]
    assert out["higher_is_better"] is True
    rmse = result(ml.export_leaderboard, [{"model": "x", "rmse": 3.0}, {"model": "y", "rmse": 1.0}], metric="rmse")
    assert [r["model"] for r in rmse["ranking"]] == ["y", "x"]
    assert rmse["content"].splitlines()[0] == "rank\tmodel\trmse"


def test_benchmark_ranks_per_dataset_and_summarises():
    rows = [
        {"model": "A", "dataset": "d1", "value": 0.9}, {"model": "B", "dataset": "d1", "value": 0.8},
        {"model": "A", "dataset": "d2", "value": 0.7}, {"model": "B", "dataset": "d2", "value": 0.8},
        {"model": "A", "dataset": "d3", "value": 0.6}, {"model": "B", "dataset": "d3", "value": 0.1},
    ]
    out = result(ml.export_benchmark, rows)
    summary = {r["model"]: r for r in out["summary"]}
    assert summary["A"]["mean_rank"] == pytest.approx((1 + 2 + 1) / 3)
    assert summary["B"]["mean_rank"] == pytest.approx((2 + 1 + 2) / 3)
    assert summary["A"]["wins"] == 2 and summary["B"]["wins"] == 1
    assert summary["A"]["mean_value"] == pytest.approx((0.9 + 0.7 + 0.6) / 3)
    assert out["summary"][0]["model"] == "A"
    with pytest.raises(InputFormatError, match="ML_DUPLICATE_RESULT"):
        run(ml.export_benchmark, rows + [rows[0]])


def test_export_dataset_formats(tmp_path: Path):
    rows = [{"id": "s1", "sequence": "ACGT", "label": 1}, {"id": "s2", "sequence": "GG", "label": 0}]
    csv_text = result(ml.export_dataset, rows)["content"]
    assert csv_text == "id,sequence,label\ns1,ACGT,1\ns2,GG,0\n"
    assert result(ml.export_dataset, rows, format="fasta")["content"] == ">s1\nACGT\n>s2\nGG\n"
    assert result(ml.export_dataset, rows, format="tsv", columns=["label", "id"])["content"] == "label\tid\n1\ts1\n0\ts2\n"
    path = tmp_path / "out.jsonl"
    info = result(ml.export_dataset, {"id": ["x", "y"], "sequence": ["A", "C"]}, output=str(path))
    assert info["n_records"] == 2
    assert [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()] == [
        {"id": "x", "sequence": "A"}, {"id": "y", "sequence": "C"},
    ]


def _records(n):
    return [{"id": f"r{i}", "sequence": "ACGT", "group": f"g{i // 3}"} for i in range(n)]


@pytest.mark.parametrize(
    ("n", "expected"),
    [(10, {"train": 8, "validation": 1, "test": 1}), (15, {"train": 12, "validation": 2, "test": 1}), (7, {"train": 5, "validation": 1, "test": 1})],
)
def test_splits_use_largest_remainder_counts_and_partition_the_rows(n, expected):
    rows = _records(n)
    parts = {name: result(getattr(ml, f"export_{name}_split"), rows, seed=3, format="json") for name in ("train", "validation", "test")}
    assert parts["train"]["counts"] == expected
    all_indices = sorted(i for p in parts.values() for i in p["row_indices"])
    assert all_indices == list(range(n))  # every row exactly once
    assert {name: p["n_records"] for name, p in parts.items()} == expected
    again = result(ml.export_train_split, rows, seed=3, format="json")
    assert again["row_indices"] == parts["train"]["row_indices"]  # deterministic
    other = result(ml.export_train_split, _records(40), seed=4)["row_indices"]
    assert other != result(ml.export_train_split, _records(40), seed=5)["row_indices"]


def test_group_split_never_splits_a_group_and_existing_split_column_wins():
    rows = _records(30)
    seen = {}
    for name in ("train", "validation", "test"):
        for index in result(getattr(ml, f"export_{name}_split"), rows, group_column="group", seed=1)["row_indices"]:
            seen.setdefault(rows[index]["group"], set()).add(name)
    assert all(len(v) == 1 for v in seen.values()) and len(seen) == 10
    fixed = [{"id": "a", "split": "train"}, {"id": "b", "split": "val"}, {"id": "c", "split": "test"}, {"id": "d", "split": "TRAIN"}]
    assert result(ml.export_validation_split, fixed)["row_indices"] == [1]
    assert result(ml.export_train_split, fixed)["row_indices"] == [0, 3]
    with pytest.raises(InputFormatError, match="ML_UNKNOWN_SPLIT"):
        run(ml.export_test_split, [{"split": "holdout"}])
    with pytest.raises(ConfigurationError, match="INVALID_SPLIT_RATIOS"):
        run(ml.export_train_split, _records(5), train_ratio=0.5)


# ---------------------------------------------------------------- ml_models: sequence operations

def test_encode_one_hot_integer_and_kmer():
    onehot = result(ml.encode, "ACGT")
    assert onehot["encodings"][0]["values"] == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    assert onehot["feature_names"] == ["A", "C", "G", "T"]
    assert result(ml.encode, ["ACGT", "TTA"], method="integer")["encodings"][1]["values"] == [3, 3, 0]
    kmers = result(ml.encode, "ACGTAC", method="kmer", k=2)
    names = kmers["feature_names"]
    counts = dict(zip(names, kmers["encodings"][0]["values"], strict=True))
    assert counts["AC"] == 2 and counts["CG"] == 1 and counts["GT"] == 1 and counts["TA"] == 1
    assert sum(counts.values()) == 5
    with pytest.raises(DNAKitError):
        run(ml.encode, "ACGN")  # ambiguity_policy 'error' is the one-hot default
    assert result(ml.encode, ">s1\nACG\n>s2\nTT\n", method="integer")["encodings"][1]["id"] == "s2"


def test_export_onehot_and_kmer_matrices(tmp_path: Path):
    tsv = result(ml.export_onehot_embedding, [{"id": "s1", "sequence": "AC"}], format="tsv")["content"]
    assert tsv.splitlines() == ["id\tposition\tA\tC\tG\tT", "s1\t0\t1\t0\t0\t0", "s1\t1\t0\t1\t0\t0"]
    path = tmp_path / "oh.json"
    info = result(ml.export_onehot_embedding, ["GT"], output=str(path))
    assert json.loads(path.read_text(encoding="utf-8"))["sequences"]["seq1"] == [[0, 0, 1, 0], [0, 0, 0, 1]]
    assert info["total_positions"] == 2
    kmer = result(ml.export_kmer_embedding, ["AAAA", "ACGT"], k=1)["content"].splitlines()
    assert kmer[0] == "id\tA\tC\tG\tT"
    assert kmer[1] == "seq1\t1\t0\t0\t0"
    assert kmer[2] == "seq2\t0.25\t0.25\t0.25\t0.25"
    counts = result(ml.export_kmer_embedding, ["AAAA"], k=2, mode="count", format="json")
    vector = json.loads(counts["content"])["vectors"]["seq1"]
    assert vector[0] == 3 and sum(vector) == 3  # AA x 3


def test_export_tokenizer_vocabulary(tmp_path: Path):
    out = result(ml.export_tokenizer, None, k=2)
    vocab = out["tokenizer"]["model"]["vocab"]
    assert out["vocab_size"] == 5 + 16 == len(vocab)
    assert vocab["[PAD]"] == 0 and vocab["[MASK]"] == 4
    assert vocab["AA"] == 5 and vocab["AC"] == 6 and vocab["TT"] == 20
    path = tmp_path / "tok.json"
    result(ml.export_tokenizer, None, k=1, special_tokens=["<pad>"], output=str(path))
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["model"]["vocab"] == {"<pad>": 0, "A": 1, "C": 2, "G": 3, "T": 4}
    with pytest.raises(ConfigurationError):
        run(ml.export_tokenizer, None, k=9)


def test_mask_dna_tokens_counts_positions_and_determinism():
    seq = "ACGT" * 5
    out = result(ml.mask_dna_tokens, seq, mask_rate=0.25)["sequences"][0]
    assert len(out["masked_positions"]) == 5 and out["tokens"] == list(seq)
    for position, (orig, masked) in enumerate(zip(out["tokens"], out["masked_tokens"], strict=True)):
        assert masked == ("[MASK]" if position in out["masked_positions"] else orig)
    assert out["labels"] == {str(p): seq[p] for p in out["masked_positions"]}
    assert result(ml.mask_dna_tokens, seq, mask_rate=0.25)["sequences"][0] == out
    kmer = result(ml.mask_dna_tokens, "ACGTACGTAC", token_size=3, mask_rate=0.5, mask_token="N")["sequences"][0]
    assert kmer["tokens"] == ["ACG", "TAC", "GTA", "C"] and len(kmer["masked_positions"]) == 2
    none = result(ml.mask_dna_tokens, "ACGT", mask_rate=0.0)["sequences"][0]
    assert none["masked_positions"] == [] and none["masked_tokens"] == list("ACGT")


def test_mutagenesis_and_attribution_with_scorer():
    def gc(seq):
        return float(sum(base in "GC" for base in seq))

    ism = result(ml.mutagenesis, "AAA", scorer=gc)
    assert ism["reference_score"] == 0.0 and ism["n_variants"] == 9
    assert {(e["position"], e["alt"]): e["delta"] for e in ism["effects"]}[(1, "G")] == 1.0
    assert all(p["mean_delta"] == pytest.approx(2 / 3) for p in ism["per_position"])

    def a_count(seq):
        return float(seq.count("A"))

    # ACG: position 0 (A) -> every substitution loses the A: attribution 1 - 0 = 1
    # position 1 (C) -> T,G keep 1 A, A gives 2 A's: attribution 1 - (2 + 1 + 1) / 3 = -1/3
    attribution = result(ml.attribute, "ACG", scorer=a_count)
    assert attribution["attributions"] == pytest.approx([1.0, -1 / 3, -1 / 3])
    only = result(ml.attribute, "ACG", scorer=a_count, positions=[1])
    assert only["positions"] == [1] and only["attributions"] == pytest.approx([-1 / 3])
    with pytest.raises(InputFormatError, match="ML_AMBIGUOUS_SEQUENCE"):
        run(ml.mutagenesis, "ACN", scorer=gc)
    with pytest.raises(InputFormatError, match="ML_BAD_SCORER_OUTPUT"):
        run(ml.mutagenesis, "AC", scorer=lambda s: float("nan"))


def test_permutation_importance_of_a_known_model():
    x = [[float(i), float((i * 7) % 5)] for i in range(12)]
    y = [row[0] for row in x]
    out = result(ml.explain, {"X": x, "y": y, "feature_names": ["signal", "noise"]}, model=lambda rows: [r[0] for r in rows], seed=1)
    by_name = {i["feature"]: i for i in out["importances"]}
    assert out["baseline_score"] == pytest.approx(1.0)
    assert by_name["noise"]["importance_mean"] == pytest.approx(0.0, abs=1e-12)  # model ignores x1
    assert by_name["signal"]["importance_mean"] > 0.5
    again = result(ml.explain, {"X": x, "y": y}, model=lambda rows: [r[0] for r in rows], seed=1)
    assert again["importances"][0]["importance_mean"] == pytest.approx(by_name["signal"]["importance_mean"])
    acc = result(ml.explain, {"X": [[0, 1], [1, 1], [0, 0], [1, 0]], "y": [0, 1, 0, 1]}, model=lambda r: [row[0] for row in r], task="classification")
    assert acc["baseline_score"] == 1.0 and acc["metric"] == "accuracy"


def test_export_interpretability_tsv_and_json():
    tsv = result(ml.export_interpretability, {"sequence": "ACG", "attributions": [0.5, -1.0, 0.0]})["content"]
    assert tsv == "position\tbase\tattribution\n0\tA\t0.5\n1\tC\t-1\n2\tG\t0\n"
    computed = result(ml.export_interpretability, "ACG", scorer=lambda s: float(s.count("A")), format="json")
    assert json.loads(computed["content"])["attributions"] == pytest.approx([1.0, -1 / 3, -1 / 3])
    with pytest.raises(InputFormatError, match="ML_LENGTH_MISMATCH"):
        run(ml.export_interpretability, {"sequence": "AC", "attributions": [1.0]})


# ---------------------------------------------------------------- ml_models: delegations

@pytest.fixture
def trained(tmp_path: Path):
    pytest.importorskip("sklearn")
    import random

    rng = random.Random(0)
    rows = []
    for i in range(40):
        letters, label = ("AC", "a") if i % 2 == 0 else ("GT", "b")
        rows.append({"id": f"r{i}", "sequence": "".join(rng.choice(letters) for _ in range(30)), "label": label})
    data = tmp_path / "data.csv"
    with data.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "sequence", "label"])
        writer.writeheader()
        writer.writerows(rows)
    options = {
        "task_type": "classification", "model": "logistic_regression", "output_dir": str(tmp_path / "model"),
        "mapping": {"sequence_field": "sequence", "label_fields": ["label"], "id_field": "id"},
    }
    return data, options, rows


def test_train_predict_load_roundtrip_on_separable_sequences(trained):
    data, options, rows = trained
    trained_out = result(ml.train, str(data), **options)
    assert trained_out["delegate"].endswith("train_custom_property_model")
    assert trained_out["training"]["finetune_method"] == "machine_learning"
    assert trained_out["training"]["metrics"]["test"]["accuracy"] == 1.0  # A/C-only vs G/T-only is separable
    predicted = result(ml.predict, str(data), model_dir=options["output_dir"], trust_model=True)["prediction"]
    assert [p["label"] for p in predicted["predictions"]] == [r["label"] for r in rows]
    manifest = result(ml.load, options["output_dir"])
    assert manifest["manifest"]["task_type"] == "classification" and manifest["artifact_present"]
    assert manifest["weights_loaded"] is False
    with pytest.raises(ConfigurationError, match="UNTRUSTED_MACHINE_LEARNING_ARTIFACT"):
        run(ml.predict, str(data), model_dir=options["output_dir"])


def test_train_validation_errors(trained, tmp_path: Path):
    data, options, _ = trained
    with pytest.raises(ConfigurationError, match="output_dir"):
        run(ml.train, str(data), task_type="classification")
    with pytest.raises(ConfigurationError, match="task_type"):
        run(ml.train, str(data), output_dir=str(tmp_path / "x"))
    with pytest.raises(ConfigurationError, match="UNKNOWN_ML_OPTION"):
        run(ml.train, str(data), bogus=1, **options)
    with pytest.raises(InputFormatError, match="ML_INPUT_NOT_FOUND"):
        run(ml.train, str(tmp_path / "missing.csv"), **options)
    deep = run(ml.train, str(data), finetune_method="full", **options)
    assert deep["status"] == "not_implemented"


def test_build_describes_untrained_estimator():
    pytest.importorskip("sklearn")
    out = result(ml.build, None, task_type="classification", model="svm", model_params={"kernel": "linear", "C": 2.0})
    assert out["estimator_class"].endswith("SVC") and out["trained"] is False
    assert out["params"]["kernel"] == "linear" and out["params"]["C"] == 2.0
    reg = result(ml.build, {"task_type": "regression", "model": "ridge"})
    assert reg["estimator_class"].endswith("Ridge")


def test_search_filters_the_model_catalogue():
    names = lambda out: {m["name"] for m in out["matches"]}  # noqa: E731
    assert "random_forest" in names(result(ml.search, "random forest"))
    assert "random_forest" in names(result(ml.search, "random_forest"))
    regression = result(ml.search, "", kind="machine_learning_estimator", task="regression", limit=1000)
    assert "ridge" in names(regression) and "logistic_regression" not in names(regression)
    assert "dnabert2" in names(result(ml.search, "dnabert2", kind="dna_embedding_model"))
    assert result(ml.search, "zzzz-no-such-model")["n_matches"] == 0
    with pytest.raises(ConfigurationError):
        run(ml.search, "x", kind="nonsense")


class StubBackend:
    """Deterministic embedding backend: [fraction A, fraction G, 1]."""

    def extract(self, sequences, *, show_progress=False):
        return [[s.count("A") / len(s), s.count("G") / len(s), 1.0] for s in sequences]


def test_sequence_embeddings_delegate_to_representation_extraction():
    out = result(ml.generate_sequence_embeddings, ["AACC", "GGTT"], representation_backend=StubBackend(), model="dnabert2")
    assert out["embeddings"] == [[0.5, 0.0, 1.0], [0.0, 0.5, 1.0]]
    assert out["record_ids"] == ["seq1", "seq2"] and out["dimension"] == 3
    assert out["delegate"].endswith("extract_representations")


def test_finetune_linear_probe_with_frozen_embeddings(trained, tmp_path: Path):
    pytest.importorskip("torch")
    data, options, _ = trained
    probe = {k: v for k, v in options.items() if k not in {"model", "output_dir"}}
    out = result(
        ml.finetune, str(data), representation_backend=StubBackend(), model="dnabert2", output_dir=str(tmp_path / "probe"),
        epochs=30, workers=0, learning_rate=0.05, **probe,
    )
    assert out["training"]["finetune_method"] == "linear_probe"
    assert out["training"]["metrics"]["test"]["accuracy"] == 1.0  # fraction of A / G separates the classes
    assert (tmp_path / "probe" / "head.npz").is_file()


def test_export_torchscript_roundtrip(tmp_path: Path):
    torch = pytest.importorskip("torch")
    net = torch.nn.Sequential(torch.nn.Linear(3, 2), torch.nn.ReLU())
    target = tmp_path / "net.pt"
    out = result(ml.export_torchscript_model, net, example_input=[[1.0, 2.0, 3.0]], output=str(target))
    assert out["roundtrip_max_abs_difference"] == 0.0 and out["input_shape"] == [1, 3]
    reloaded = torch.jit.load(str(target))
    sample = torch.tensor([[0.5, -1.0, 2.0]])
    assert torch.allclose(reloaded(sample), net(sample))
    with pytest.raises(ConfigurationError, match="ML_OUTPUT_EXISTS"):
        run(ml.export_torchscript_model, net, example_input=[[1.0, 2.0, 3.0]], output=str(target))
    with pytest.raises(InputFormatError, match="ML_NOT_A_TORCH_MODULE"):
        run(ml.export_torchscript_model, "not a model", example_input=[[1.0]], output=str(tmp_path / "x.pt"))


# ---------------------------------------------------------------- ml_models: leftovers and invalid input

ML_NOT_IMPLEMENTED = [
    "pretrain", "serve", "quantize_model", "run_cross_species_transfer", "export_model_onnx", "export_onnx_model",
    "export_chromatin_model", "export_cnn_model", "export_expression_model", "export_graph_model",
    "export_long_context_model", "export_pangenome_model", "export_phenotype_model", "export_protein_model",
    "export_rnn_model", "export_splice_model", "export_transformer_model", "export_variant_effect_model",
]


@pytest.mark.parametrize("name", ML_NOT_IMPLEMENTED)
def test_ml_models_leftovers_report_not_implemented(name):
    envelope = getattr(ml, name)({"anything": 1}, verbose=False)
    assert envelope["status"] == "not_implemented"
    assert envelope["required_dependency"] and envelope["message"]
    assert "result" not in envelope


@pytest.mark.parametrize(
    ("name", "value"),
    [("mutagenesis", "ACGT"), ("attribute", "ACGT"), ("explain", {"X": [[1.0]], "y": [1.0]}),
     ("export_interpretability", "ACGT"), ("generate_sequence_embeddings", ["ACGT"]), ("finetune", "missing.csv")],
)
def test_model_dependent_functions_without_a_model_report_not_implemented(name, value):
    assert getattr(ml, name)(value, verbose=False)["status"] == "not_implemented"


def test_every_ml_function_is_accounted_for():
    implemented = set(ml.__all__) - set(ML_NOT_IMPLEMENTED)
    assert len(ml.__all__) == 48 and len(implemented) == 30
    for name in ml.__all__:  # all keep the original signature
        parameters = list(inspect.signature(getattr(ml, name)).parameters)
        assert parameters[:2] == ["input", "options"]


@pytest.mark.parametrize(
    ("function", "value", "options", "error"),
    [
        (ml.evaluate, {"y_true": [1, 2], "y_pred": [1]}, {}, "ML_LENGTH_MISMATCH"),
        (ml.evaluate, {"y_true": [1, 2]}, {}, "ML_COLUMN_NOT_FOUND"),
        (ml.evaluate, "not a path", {}, "INVALID_ML_INPUT|ML_INPUT"),
        (ml.evaluate, {"y_true": [0, 1], "y_pred": [0, 1]}, {"metrics": ["nope"]}, "UNKNOWN_ML_METRIC"),
        (ml.calibrate, {"y_true": [1, 1], "y_score": [0.2, 0.9]}, {}, "ML_AMBIGUOUS_POSITIVE|CALIBRATION_ONE_CLASS"),
        (ml.calculate_prediction_uncertainty, {"probabilities": [[0.9, 0.3]]}, {}, "INVALID_ML_INPUT"),
        (ml.calculate_prediction_uncertainty, {"x": 1}, {}, "ML_COLUMN_NOT_FOUND"),
        (ml.export_leaderboard, [{"model": "a"}], {}, "ML_COLUMN_NOT_FOUND"),
        (ml.encode, [], {}, "INVALID_ML_INPUT"),
        (ml.mask_dna_tokens, "ACGT", {"mask_rate": 2}, "INVALID_ML_OPTION"),
        (ml.export_dataset, [], {}, "INVALID_ML_INPUT"),
        (ml.export_dataset, 42, {}, "INVALID_ML_INPUT"),
    ],
)
def test_ml_invalid_input_raises_dnakit_errors(function, value, options, error):
    with pytest.raises(DNAKitError) as caught:
        function(value, options=options, verbose=False)
    assert any(code in caught.value.code for code in error.split("|"))


def test_options_must_be_a_mapping():
    with pytest.raises(ConfigurationError, match="INVALID_ML_OPTIONS"):
        ml.evaluate(BINARY, options=["bad"], verbose=False)
    with pytest.raises(ConfigurationError, match="INVALID_GRAPH_OPTIONS"):
        kg.build_knowledge_graph([["a", "r", "b"]], options=["bad"], verbose=False)


# ---------------------------------------------------------------- knowledge_graph

TRIPLES = [
    ["BRCA1", "associated_with", "breast cancer"],
    ["TP53", "regulates", "MDM2"],
    ["MDM2", "inhibits", "TP53"],
    ["BRCA1", "interacts_with", "TP53"],
]
CHAIN = [["A", "next", "B"], ["B", "next", "C"], ["C", "next", "D"]]


def edge_set(graph):
    return {(e["source"], e["relation"], e["target"]) for e in graph["edges"]}


def test_build_merges_duplicates_and_accepts_every_input_form(tmp_path: Path):
    out = result(kg.build_knowledge_graph, TRIPLES + [TRIPLES[0]])
    assert out["n_nodes"] == 4 and out["n_edges"] == 4 and out["duplicate_edges"] == 1
    assert out["relation_counts"] == {"associated_with": 1, "inhibits": 1, "interacts_with": 1, "regulates": 1}
    merged = next(e for e in out["graph"]["edges"] if e["relation"] == "associated_with")
    assert merged["count"] == 2
    mappings = result(kg.build_knowledge_graph, [
        {"subject": "a", "predicate": "p", "object": "b", "weight": 0.5},
        {"head": "b", "relation": "q", "tail": "c"},
        {"from": "c", "to": "a"},
    ])
    assert edge_set(mappings["graph"]) == {("a", "p", "b"), ("b", "q", "c"), ("c", "related_to", "a")}
    assert next(e for e in mappings["graph"]["edges"] if e["relation"] == "p")["weight"] == 0.5
    text = "subject\tpredicate\tobject\nx\tr\ty\ny\tr\tz\n"
    assert result(kg.build_knowledge_graph, text)["n_edges"] == 2
    path = tmp_path / "edges.tsv"
    path.write_text("a\tr\tb\n# comment\nb\tr\tc\n", encoding="utf-8")
    assert result(kg.build_knowledge_graph, str(path))["n_nodes"] == 3
    graph = result(kg.build_knowledge_graph, {"nodes": [{"id": "g1", "type": "gene"}, "lonely"], "edges": [["g1", "r", "g2"]]})
    assert graph["n_nodes"] == 3 and graph["graph"]["nodes"][0] == {"id": "g1", "type": "gene"}


def test_undirected_graphs_treat_edge_endpoints_as_unordered():
    out = result(kg.build_knowledge_graph, [["a", "r", "b"], ["b", "r", "a"]], directed=False)
    assert out["n_edges"] == 1 and out["duplicate_edges"] == 1 and out["directed"] is False
    loops = result(kg.build_knowledge_graph, [["a", "r", "a"], ["a", "r", "b"]], drop_self_loops=True)
    assert loops["n_edges"] == 1 and loops["dropped_self_loops"] == 1


def test_create_with_seed_nodes_and_edges():
    empty = result(kg.create, None, name="demo")
    assert empty["n_nodes"] == 0 and empty["graph"]["name"] == "demo"
    seeded = result(kg.create, None, nodes=["a", {"id": "b", "type": "x"}], edges=[["a", "r", "b"]], directed=False)
    assert seeded["n_nodes"] == 2 and seeded["directed"] is False and seeded["n_edges"] == 1


def test_add_graph_relationship_does_not_mutate_and_handles_duplicates():
    original = [["a", "r", "b"]]
    out = result(kg.add_graph_relationship, original, subject="b", predicate="r", object="c")
    assert out["n_edges"] == 2 and out["added"] == [{"source": "b", "relation": "r", "target": "c"}]
    assert original == [["a", "r", "b"]]
    again = result(kg.add_graph_relationship, original, relationships=[["a", "r", "b"], ["x", "r", "y"]])
    assert again["already_present"] == [{"source": "a", "relation": "r", "target": "b"}] and again["n_edges"] == 2
    with pytest.raises(InputFormatError, match="DUPLICATE_EDGE"):
        run(kg.add_graph_relationship, original, subject="a", predicate="r", object="b", on_duplicate="error")
    with pytest.raises(InputFormatError, match="UNKNOWN_NODE"):
        run(kg.add_graph_relationship, original, subject="a", predicate="r", object="zz", create_missing_nodes=False)
    with pytest.raises(ConfigurationError):
        run(kg.add_graph_relationship, original)


def test_update_removes_nodes_with_their_edges_and_sets_properties():
    out = result(
        kg.update, TRIPLES, remove_nodes=["MDM2"], set_node_properties={"TP53": {"type": "gene"}},
        set_edge_properties=[{"source": "BRCA1", "relation": "interacts_with", "target": "TP53", "properties": {"score": 0.9}}],
    )
    assert out["nodes_removed"] == 1 and out["edges_removed"] == 2 and out["n_edges"] == 2
    assert {n["id"]: n for n in out["graph"]["nodes"]}["TP53"]["type"] == "gene"
    assert next(e for e in out["graph"]["edges"] if e["relation"] == "interacts_with")["score"] == 0.9
    removed = result(kg.update, TRIPLES, remove_edges=[["TP53", "regulates", "MDM2"]])
    assert removed["n_edges"] == 3 and removed["n_nodes"] == 4
    with pytest.raises(InputFormatError, match="UNKNOWN_NODE"):
        run(kg.update, TRIPLES, remove_nodes=["nope"])
    assert result(kg.update, TRIPLES, remove_nodes=["nope"], ignore_missing=True)["n_nodes"] == 4
    with pytest.raises(ConfigurationError, match="no change"):
        run(kg.update, TRIPLES)


def test_link_pairs_and_shared_property():
    out = result(kg.link, [["a", "r", "b"]], pairs=[["a", "c"], ["a", "b"]], relation="same_as")
    assert out["n_created"] == 2 and ("a", "same_as", "c") in edge_set(out["graph"])
    graph = {"nodes": [{"id": "n1", "symbol": "TP53"}, {"id": "n2", "symbol": "tp53"}, {"id": "n3", "symbol": " TP53 "}, {"id": "n4", "symbol": "MDM2"}]}
    shared = result(kg.link, graph, by="symbol")
    assert shared["n_created"] == 3  # C(3, 2) pairs among n1, n2, n3
    assert {(c["source"], c["target"]) for c in shared["created"]} == {("n1", "n2"), ("n1", "n3"), ("n2", "n3")}
    with pytest.raises(InputFormatError, match="GRAPH_GROUP_TOO_LARGE"):
        run(kg.link, graph, by="symbol", max_group_size=2)
    with pytest.raises(ConfigurationError):
        run(kg.link, graph)


def test_map_renames_and_merges_nodes():
    triples = [["TP53", "regulates", "MDM2"], ["P53", "binds", "DNA"], ["P53", "regulates", "MDM2"]]
    out = result(kg.map, triples, mapping={"TP53": "ENSG1", "P53": "ENSG1", "MDM2": "ENSG2"})
    assert out["n_nodes_merged"] == 1 and out["n_ids_renamed"] == 3
    assert edge_set(out["graph"]) == {("ENSG1", "regulates", "ENSG2"), ("ENSG1", "binds", "DNA")}
    assert next(e for e in out["graph"]["edges"] if e["relation"] == "regulates")["count"] == 2
    dropped = result(kg.map, triples, mapping={"TP53": "X", "MDM2": "Y"}, unmapped="drop")
    assert edge_set(dropped["graph"]) == {("X", "regulates", "Y")} and dropped["n_nodes_dropped"] == 2
    with pytest.raises(InputFormatError, match="UNMAPPED_NODE"):
        run(kg.map, triples, mapping={"TP53": "X"}, unmapped="error")
    renamed = result(kg.map, triples[:1], mapping=[{"from": "TP53", "to": "T"}], relation_mapping={"regulates": "controls"})
    assert edge_set(renamed["graph"]) == {("T", "controls", "MDM2")}
    with pytest.raises(ConfigurationError):
        run(kg.map, triples)


def test_aggregate_degree_statistics_hand_values():
    graph = {"nodes": ["A", "B", "C", "D", "E", "F"], "edges": [["A", "r", "B"], ["B", "s", "C"]]}
    out = result(kg.aggregate, graph, top_k=2)
    # total degrees A=1 B=2 C=1 D=E=F=0
    assert out["degree"]["total"] == {"min": 0, "max": 2, "mean": 4 / 6, "median": 0.5}
    assert out["degree"]["in"]["max"] == 1 and out["degree"]["out"]["max"] == 1
    assert out["density"] == pytest.approx(2 / (6 * 5))
    assert out["n_isolated_nodes"] == 3 and out["n_weakly_connected_components"] == 4
    assert out["largest_component_size"] == 3
    assert out["top_hubs"][0] == {"id": "B", "degree": 2, "in_degree": 1, "out_degree": 1} and len(out["top_hubs"]) == 2
    assert out["relation_counts"] == {"r": 1, "s": 1}
    typed = result(kg.aggregate, {"nodes": [{"id": "a", "type": "gene"}, {"id": "b", "type": "gene"}, {"id": "c"}], "edges": []})
    assert typed["node_type_counts"] == {"gene": 2, "unknown": 1} and typed["density"] == 0.0
    undirected = result(kg.aggregate, [["a", "r", "b"]], directed=False)
    assert undirected["density"] == 1.0


def test_query_neighbors_depth_direction_and_shortest_path():
    out = result(kg.query_graph_neighbors, CHAIN, node="A", depth=2, direction="out")
    assert [(n["id"], n["distance"]) for n in out["neighbors"]] == [("B", 1), ("C", 2)]
    assert out["neighbors"][1]["path"] == ["A", "B", "C"]
    assert result(kg.query_graph_neighbors, CHAIN, node="C", direction="in", depth=5)["neighbors"][0]["id"] == "B"
    both = result(kg.query_graph_neighbors, CHAIN, node="B", direction="both")
    assert {n["id"] for n in both["neighbors"]} == {"A", "C"}
    assert result(kg.query_graph_neighbors, CHAIN, node="D", direction="out")["n_neighbors"] == 0
    path = result(kg.query_graph_neighbors, CHAIN, node="A", direction="out", target="D")["shortest_path"]
    assert path["found"] and path["length"] == 3 and path["nodes"] == ["A", "B", "C", "D"]
    assert [e["target"] for e in path["edges"]] == ["B", "C", "D"]
    assert not result(kg.query_graph_neighbors, CHAIN, node="D", direction="out", target="A")["shortest_path"]["found"]
    assert result(kg.query_graph_neighbors, CHAIN, node="D", direction="both", target="A")["shortest_path"]["length"] == 3
    shortcut = CHAIN + [["A", "jump", "D"]]
    assert result(kg.query_graph_neighbors, shortcut, node="A", direction="out", target="D")["shortest_path"]["length"] == 1
    filtered = result(kg.query_graph_neighbors, shortcut, node="A", direction="out", relations=["next"], depth=3)
    assert [n["id"] for n in filtered["neighbors"]] == ["B", "C", "D"] and filtered["neighbors"][2]["distance"] == 3
    sub = result(kg.query_graph_neighbors, CHAIN, node="A", direction="out", depth=1, include_subgraph=True)["subgraph"]
    assert edge_set(sub) == {("A", "next", "B")}
    undirected = result(kg.query_graph_neighbors, [["a", "r", "b"], ["c", "r", "b"]], node="a", depth=2, directed=False)
    assert {n["id"]: n["distance"] for n in undirected["neighbors"]} == {"b": 1, "c": 2}
    with pytest.raises(InputFormatError, match="UNKNOWN_NODE"):
        run(kg.query_graph_neighbors, CHAIN, node="Z")
    with pytest.raises(ConfigurationError):
        run(kg.query_graph_neighbors, CHAIN)


def test_search_nodes_by_text_and_type():
    graph = {"nodes": [
        {"id": "BRCA1", "type": "gene", "name": "breast cancer 1"},
        {"id": "BRCA2", "type": "gene", "name": "breast cancer 2"},
        {"id": "breast cancer", "type": "disease"},
    ]}
    assert [m["id"] for m in result(kg.search, graph, query="brca")["matches"]] == ["BRCA1", "BRCA2"]
    assert [m["id"] for m in result(kg.search, graph, query="BRCA1", match="exact")["matches"]] == ["BRCA1"]
    assert result(kg.search, graph, query="brca", case_sensitive=True)["n_matches"] == 0
    assert [m["id"] for m in result(kg.search, graph, query="breast", node_type="disease")["matches"]] == ["breast cancer"]
    assert [m["id"] for m in result(kg.search, graph, query="breast", fields=["name"])["matches"]] == ["BRCA1", "BRCA2"]
    prefix = result(kg.search, graph, query="bre", match="prefix", limit=1)
    assert prefix["n_matches"] == 3 and prefix["truncated"] and len(prefix["matches"]) == 1  # ids and names


@pytest.mark.parametrize("fmt", ["json", "jsonl", "tsv", "graphml"])
def test_export_then_read_roundtrip(tmp_path: Path, fmt):
    graph = {"nodes": [{"id": "g1", "type": "gene"}, {"id": "d<1>&"}], "edges": [{"source": "g1", "relation": "assoc", "target": "d<1>&"}]}
    target = tmp_path / f"graph.{fmt}"
    out = result(kg.export, graph, output=str(target))
    assert out["format"] == fmt and out["output_path"] == str(target) and target.is_file()
    back = result(kg.read, str(target))
    assert edge_set(back["graph"]) == {("g1", "assoc", "d<1>&")}
    assert back["n_nodes"] == 2
    if fmt != "tsv":
        assert {n["id"]: n for n in back["graph"]["nodes"]}["g1"]["type"] == "gene"
    if fmt == "tsv":
        assert target.read_text(encoding="utf-8") == "source\trelation\ttarget\ng1\tassoc\td<1>&\n"
    with pytest.raises(ConfigurationError, match="GRAPH_OUTPUT_EXISTS"):
        run(kg.export, graph, output=str(target))
    assert result(kg.export, graph, output=str(target), overwrite=True)["bytes"] > 0


def test_export_inline_content_and_notes():
    tsv = result(kg.export, {"nodes": ["lonely"], "edges": [["a", "r", "b"]]}, format="tsv")
    assert "isolated" in tsv["note"] and tsv["content"].splitlines()[1] == "a\tr\tb"
    graphml = result(kg.export, [["a", "r", "b"]], format="graphml")["content"]
    assert 'edgedefault="directed"' in graphml and 'source="a" target="b"' in graphml
    undirected = result(kg.export, [["a", "r", "b"]], format="graphml", directed=False)["content"]
    assert 'edgedefault="undirected"' in undirected
    assert result(kg.read, undirected)["directed"] is False
    with pytest.raises(ConfigurationError):
        run(kg.export, TRIPLES, format="pdf")


def test_read_rejects_unsafe_or_unsupported_input(tmp_path: Path):
    bomb = '<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "aaaa">]><graphml><graph/></graphml>'
    with pytest.raises(InputFormatError, match="UNSAFE_XML"):
        run(kg.read, bomb)
    other = tmp_path / "graph.xlsx"
    other.write_text("x", encoding="utf-8")
    with pytest.raises(InputFormatError, match="UNSUPPORTED_GRAPH_FORMAT"):
        run(kg.read, str(other))
    with pytest.raises(InputFormatError, match="INVALID_EDGE"):
        run(kg.read, [["only-one"]])
    with pytest.raises(InputFormatError, match="GRAPH_INPUT_TOO_LARGE"):
        run(kg.read, [["a", "r", "b"], ["b", "r", "c"]], max_edges=1)
    with pytest.raises(InputFormatError):
        run(kg.read, 42)
    with pytest.raises(InputFormatError, match="INVALID_NODE_ID"):
        run(kg.read, [["a", "r", ""]])


def test_validate_reports_each_kind_of_problem():
    graph = {
        "nodes": [{"id": "a"}, {"id": "a"}, {"id": "b"}, {"id": "iso"}],
        "edges": [["a", "r", "b"], ["a", "r", "b"], ["b", "r", "zzz"], ["b", "bad", "a"], ["a", "r", "a"]],
    }
    out = result(kg.validate, graph, allowed_relations=["r"], allow_self_loops=False)
    errors = {e["code"] for e in out["errors"]}
    assert errors == {"DUPLICATE_NODE", "DANGLING_EDGE", "UNKNOWN_RELATION", "SELF_LOOP"}
    assert {w["code"] for w in out["warnings"]} == {"DUPLICATE_EDGE", "ISOLATED_NODE"}
    assert not out["valid"]
    clean = result(kg.validate, [["a", "r", "b"]])
    assert clean["valid"] and clean["errors"] == [] and clean["n_nodes"] == 2
    loops_ok = result(kg.validate, [["a", "r", "a"]])
    assert loops_ok["valid"] and loops_ok["warnings"][0]["code"] == "SELF_LOOP"
    empty_relation = result(kg.validate, [{"source": "a", "relation": " ", "target": "b"}])
    assert [e["code"] for e in empty_relation["errors"]] == ["MISSING_RELATION"]


def test_version_is_order_independent_and_matches_its_definition():
    a = result(kg.version, TRIPLES)
    b = result(kg.version, list(reversed(TRIPLES)))
    assert a["version_id"] == b["version_id"] and a["short_id"] == a["version_id"][7:19]
    changed = result(kg.version, TRIPLES + [["x", "r", "y"]], previous_version=a["version_id"], label="v2")
    assert changed["version_id"] != a["version_id"] and changed["changed_since_previous"] is True and changed["label"] == "v2"
    assert result(kg.version, TRIPLES, previous_version=a["version_id"])["changed_since_previous"] is False
    # recompute the documented definition independently for a one-edge graph
    canonical = {"directed": True, "nodes": [["a", {}], ["b", {}]], "edges": [["a", "r", "b", {}]]}
    expected = "sha256:" + hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert result(kg.version, [["a", "r", "b"]])["version_id"] == expected
    assert result(kg.version, TRIPLES, directed=False)["version_id"] != a["version_id"]


def test_snapshot_writes_a_reproducible_file(tmp_path: Path):
    target = tmp_path / "snap.json"
    out = result(kg.snapshot_database_version, TRIPLES, label="release-1", timestamp="2026-01-01T00:00:00Z", output=str(target))
    saved = json.loads(target.read_text(encoding="utf-8"))
    assert saved["version_id"] == out["version_id"] == result(kg.version, TRIPLES)["version_id"]
    assert saved["label"] == "release-1" and saved["timestamp"] == "2026-01-01T00:00:00Z"
    assert len(saved["graph"]["edges"]) == 4
    inline = result(kg.snapshot_database_version, TRIPLES)
    assert "timestamp" not in inline["snapshot"] and inline["output_path"] is None
    with pytest.raises(ConfigurationError, match="GRAPH_OUTPUT_EXISTS"):
        run(kg.snapshot_database_version, TRIPLES, output=str(target))
    # a snapshot can itself be compared with the live graph
    assert result(kg.compare_database_versions, saved["graph"], other=TRIPLES)["identical"] is True


def test_compare_versions_added_removed_changed():
    old = {"nodes": [{"id": "a", "type": "x"}, {"id": "b"}], "edges": [["a", "r", "b"], ["a", "s", "b"]]}
    new = {"nodes": [{"id": "a", "type": "y"}, {"id": "b"}, {"id": "c"}], "edges": [["a", "r", "b"], ["b", "r", "c"], {"source": "a", "relation": "s", "target": "b", "w": 2}]}
    out = result(kg.compare_database_versions, {"old": old, "new": new})
    assert out["counts"] == {"nodes_added": 1, "nodes_removed": 0, "nodes_changed": 1, "edges_added": 1, "edges_removed": 0, "edges_changed": 1}
    assert out["added_nodes"] == ["c"] and out["changed_nodes"][0]["new"] == {"type": "y"}
    assert (out["added_edges"][0]["source"], out["added_edges"][0]["target"]) == ("b", "c")
    assert out["changed_edges"][0]["new"] == {"w": 2} and not out["identical"]
    same = result(kg.compare_database_versions, [old, old])
    assert same["identical"] and same["old_version_id"] == same["new_version_id"]
    via_other = result(kg.compare_database_versions, new, other=old)
    assert via_other["counts"]["nodes_added"] == 0 and via_other["counts"]["nodes_removed"] == 1
    with pytest.raises(InputFormatError):
        run(kg.compare_database_versions, old)


def test_data_lineage_order_roots_leaves_and_trace():
    records = [
        {"output": "aligned.bam", "inputs": ["reads.fq", "ref.fa"], "step": "bwa"},
        {"output": "variants.vcf", "inputs": "aligned.bam", "step": "caller"},
        {"output": "report.pdf", "inputs": ["variants.vcf", "ref.fa"], "tool": "report"},
    ]
    out = result(kg.build_data_lineage, records, trace="report.pdf")
    assert sorted(out["roots"]) == ["reads.fq", "ref.fa"]
    assert out["leaves"] == ["report.pdf"]
    order = out["topological_order"]
    assert order.index("reads.fq") < order.index("aligned.bam") < order.index("variants.vcf") < order.index("report.pdf")
    assert ("variants.vcf", "derived_from", "aligned.bam") in edge_set(out["graph"])
    assert next(e for e in out["graph"]["edges"] if e["source"] == "aligned.bam" and e["target"] == "reads.fq")["step"] == "bwa"
    ancestors = {a["id"]: a["distance"] for a in out["upstream_of"]["ancestors"]}
    assert ancestors == {"variants.vcf": 1, "ref.fa": 1, "aligned.bam": 2, "reads.fq": 3}
    with pytest.raises(InputFormatError, match="LINEAGE_CYCLE"):
        run(kg.build_data_lineage, [{"output": "a", "inputs": ["b"]}, {"output": "b", "inputs": ["a"]}])
    with pytest.raises(InputFormatError, match="LINEAGE_CYCLE"):
        run(kg.build_data_lineage, [{"output": "a", "inputs": ["a"]}])
    with pytest.raises(InputFormatError, match="INVALID_LINEAGE"):
        run(kg.build_data_lineage, [{"inputs": ["a"]}])


def test_access_policy_filters_by_clearance():
    graph = {
        "nodes": [
            {"id": "gene", "classification": "public", "note": "n1"},
            {"id": "cohort", "classification": "internal", "note": "n2"},
            {"id": "patients", "classification": "secret"},
        ],
        "edges": [
            {"source": "gene", "relation": "in", "target": "cohort", "classification": "internal"},
            {"source": "cohort", "relation": "has", "target": "patients", "classification": "public"},
            {"source": "gene", "relation": "x", "target": "cohort", "classification": "secret"},
        ],
    }
    internal = result(kg.apply_access_policy, graph, clearance="internal", redact_properties=["note"])
    assert {n["id"] for n in internal["graph"]["nodes"]} == {"gene", "cohort"}
    assert internal["n_nodes_removed"] == 1 and internal["removed_nodes"] == ["patients"]
    assert edge_set(internal["graph"]) == {("gene", "in", "cohort")}
    assert internal["n_edges_removed"] == 2
    assert all("note" not in n for n in internal["graph"]["nodes"])
    public = result(kg.apply_access_policy, graph)
    assert [n["id"] for n in public["graph"]["nodes"]] == ["gene"] and public["n_edges"] == 0
    top = result(kg.apply_access_policy, graph, clearance="secret", deny_nodes=["cohort"])
    assert {n["id"] for n in top["graph"]["nodes"]} == {"gene", "patients"} and top["n_edges"] == 0
    assert result(kg.apply_access_policy, TRIPLES)["n_nodes"] == 0  # unclassified items fail closed
    assert result(kg.apply_access_policy, TRIPLES, unclassified="allow")["n_edges"] == 4
    with pytest.raises(ConfigurationError, match="UNKNOWN_CLEARANCE"):
        run(kg.apply_access_policy, graph, clearance="top")
    with pytest.raises(InputFormatError, match="UNKNOWN_CLASSIFICATION"):
        run(kg.apply_access_policy, {"nodes": [{"id": "a", "classification": "weird"}]})


def test_knowledge_graph_function_inventory():
    assert len(kg.__all__) == 18
    # pseudonymize_genomic_metadata_fields replaced the misleadingly named earlier de-identification function
    redacted = kg.pseudonymize_genomic_metadata_fields({"name": "Jane", "kept": 1}, verbose=False)
    assert redacted["status"] == "ok" and redacted["result"]["metadata"] == {"name": "[REDACTED]", "kept": 1}
    for name in kg.__all__:
        parameters = list(inspect.signature(getattr(kg, name)).parameters)
        assert parameters[:2] == ["input", "options"]
