"""Hand-verified tests for the real singlecell_spatial and transcriptomics functions.

Expected values are either computed by hand (comments show the arithmetic) or taken from an independent
oracle (scipy / numpy / scikit-learn) applied to the textbook definition.
"""

from __future__ import annotations

import json
import math
import stat
import sys
from pathlib import Path

import pytest

import dnakit.singlecell_spatial as sc
import dnakit.transcriptomics as tx
from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError
from dnakit.singlecell_spatial import _impl_embed, _impl_matrix, _impl_stats

np = pytest.importorskip("numpy")
scipy_stats = pytest.importorskip("scipy.stats")
sklearn_metrics = pytest.importorskip("sklearn.metrics")


def ok(envelope: dict) -> dict:
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


def column(result_matrix: dict, cell: str) -> list[float]:
    j = result_matrix["cells"].index(cell)
    return [row[j] for row in result_matrix["values"]]


# ---------------------------------------------------------------------------------------------
# input forms (shared parser)
# ---------------------------------------------------------------------------------------------

MAPPING = {"g1": [1, 0, 3], "g2": [0, 0, 2]}


def test_all_matrix_forms_give_the_same_matrix(tmp_path: Path) -> None:
    cells = ["c1", "c2", "c3"]
    csv_text = "gene,c1,c2,c3\ng1,1,0,3\ng2,0,0,2\n"
    tsv_path = tmp_path / "m.tsv"
    tsv_path.write_text(csv_text.replace(",", "\t"))
    mtx = tmp_path / "tenx"
    mtx.mkdir()
    # entries are (row, col, value), 1-based: g1/c1=1, g1/c3=3, g2/c3=2
    (mtx / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 3 3\n1 1 1\n1 3 3\n2 3 2\n")
    (mtx / "features.tsv").write_text("id1\tg1\tGene Expression\nid2\tg2\tGene Expression\n")
    (mtx / "barcodes.tsv").write_text("c1\nc2\nc3\n")
    forms = [
        (MAPPING, {"cells": cells}),
        ({"g1": {"c1": 1, "c3": 3}, "g2": {"c3": 2}}, {"cells": cells}),
        ({"genes": ["g1", "g2"], "cells": cells, "counts": [[1, 0, 3], [0, 0, 2]]}, {}),
        ([[1, 0, 3], [0, 0, 2]], {"genes": ["g1", "g2"], "cells": cells}),
        (np.array([[1, 0, 3], [0, 0, 2]]), {"genes": ["g1", "g2"], "cells": cells}),
        (csv_text, {}),
        (tsv_path, {}),
        (str(tsv_path), {}),
        (mtx, {}),
        ([[1, 0], [0, 0], [3, 2]], {"genes": ["g1", "g2"], "cells": cells, "orientation": "cells_by_genes"}),
        ([{"g1": 1, "g2": 0}, {"g1": 0, "g2": 0}, {"g1": 3, "g2": 2}], {"cells": cells, "orientation": "cells_by_genes"}),
    ]
    for value, options in forms:
        parsed = _impl_matrix.parse_matrix(value, options)
        assert parsed.genes == ["g1", "g2"], value
        assert parsed.cells == cells, value
        assert parsed.values == [[1, 0, 3], [0, 0, 2]], value


def test_pandas_dataframe_is_accepted() -> None:
    pd = pytest.importorskip("pandas")
    frame = pd.DataFrame({"c1": [1, 0], "c2": [3, 2]}, index=["g1", "g2"])
    parsed = _impl_matrix.parse_matrix(frame, {})
    assert (parsed.genes, parsed.cells, parsed.values) == (["g1", "g2"], ["c1", "c2"], [[1, 3], [0, 2]])


def test_read_summary() -> None:
    result = ok(sc.read(MAPPING, options={"cells": ["c1", "c2", "c3"]}))
    # 5 non-zero... entries: g1 has 1,3 and g2 has 2 -> 3 non-zeros of 6; total counts 1+3+2 = 6
    assert (result["n_genes"], result["n_cells"], result["total_counts"], result["nonzero_entries"]) == (2, 3, 6.0, 3)
    assert result["sparsity"] == 0.5
    assert result["matrix"]["values"] == [[1.0, 0.0, 3.0], [0.0, 0.0, 2.0]]
    assert ok(tx.read(MAPPING))["n_cells"] == 3


def test_read_spatial_with_coordinates_mapping_and_visium_file(tmp_path: Path) -> None:
    cells = ["c1", "c2", "c3"]
    coords = {"c1": [0, 0], "c2": [10, 0], "c3": [10, 5]}
    result = ok(sc.read_spatial(MAPPING, options={"cells": cells, "coordinates": coords}))
    assert result["x_range"] == [0.0, 10.0] and result["y_range"] == [0.0, 5.0]
    assert result["coordinates"]["c3"] == [10.0, 5.0]
    # headerless Visium v1 file: barcode,in_tissue,array_row,array_col,pxl_row,pxl_col -> x = pxl_col, y = pxl_row
    positions = tmp_path / "tissue_positions_list.csv"
    positions.write_text("c1,1,0,0,100,200\nc2,1,0,1,110,220\nc3,0,0,2,120,240\n")
    result = ok(sc.read_spatial(MAPPING, options={"cells": cells, "coordinates": str(positions)}))
    assert result["coordinates"]["c2"] == [220.0, 110.0]


# ---------------------------------------------------------------------------------------------
# validate / qc / normalize
# ---------------------------------------------------------------------------------------------


def test_validate_reports_problems() -> None:
    good = ok(sc.validate(MAPPING))
    assert good["valid"] is True and good["errors"] == []
    # one all-zero cell and one all-zero gene are warnings, not errors
    zero = ok(sc.validate({"a": [1, 0], "b": [0, 0]}))
    assert zero["valid"] is True and zero["empty_cells"] == 1 and zero["empty_genes"] == 1
    fractional = ok(sc.validate({"a": [0.5, 1]}))
    assert fractional["valid"] is False and "non-integer" in fractional["errors"][0]
    assert ok(sc.validate({"a": [0.5, 1]}, options={"require_integer": False}))["valid"] is True
    negative = ok(sc.validate({"a": [-1, 1]}))
    assert negative["valid"] is False and "1 negative" in negative["errors"][0]
    duplicate = ok(sc.validate("gene,c1,c2\na,1,2\na,3,4\n"))
    assert duplicate["valid"] is False and "duplicate gene" in duplicate["errors"][0]
    ragged = ok(sc.validate("gene,c1,c2\na,1\n"))
    assert ragged["valid"] is False


def test_validate_spatial_needs_coordinates_for_every_spot() -> None:
    cells = ["c1", "c2", "c3"]
    assert ok(sc.validate_spatial(MAPPING, options={"cells": cells}))["valid"] is False
    complete = ok(sc.validate_spatial(MAPPING, options={"cells": cells, "coordinates": {"c1": [0, 0], "c2": [1, 1], "c3": [2, 2]}}))
    assert complete["valid"] is True
    partial = ok(sc.validate_spatial(MAPPING, options={"cells": cells, "coordinates": {"c1": [0, 0]}}))
    assert partial["valid"] is False and "no coordinates" in partial["errors"][0]


def test_qc_metrics_and_thresholds() -> None:
    # cell x: MT-A=2 B=4 C=4 (total 10, 3 genes, 20 % mito); cell y: B=5 C=5 (total 10, 2 genes, 0 % mito)
    matrix = {"MT-A": [2, 0], "B": [4, 5], "C": [4, 5]}
    result = ok(sc.qc(matrix, options={"cells": ["x", "y"], "max_pct_mito": 10, "min_genes": 3}))
    per_cell = {row["cell"]: row for row in result["per_cell"]}
    assert per_cell["x"]["total_counts"] == 10 and per_cell["x"]["n_genes_detected"] == 3
    assert per_cell["x"]["pct_mito"] == pytest.approx(20.0)
    assert per_cell["y"]["pct_mito"] == 0.0 and per_cell["y"]["n_genes_detected"] == 2
    assert per_cell["x"]["fails"] == ["max_pct_mito"] and per_cell["y"]["fails"] == ["min_genes"]
    assert result["n_cells_failing"] == 2 and result["n_mito_genes"] == 1
    assert result["median_genes_detected"] == 2.5
    spatial = ok(sc.qc_spatial(matrix, options={"cells": ["x", "y"], "coordinates": {"x": [0, 1], "y": [4, 3]}}))
    assert spatial["x_range"] == [0.0, 4.0] and spatial["n_cells_failing"] == 0


@pytest.mark.parametrize("function", [sc.normalize, sc.normalize_spatial])
def test_normalize_log1p_total_and_cpm(function) -> None:
    matrix = {"a": [1, 3], "b": [3, 1]}  # both cells total 4
    log1p = ok(function(matrix, options={"target_sum": 100}))["matrix"]
    # value = count / 4 * 100 -> 25 and 75; ln(1 + x)
    assert log1p["values"][0] == pytest.approx([math.log(26), math.log(76)])
    assert log1p["values"][1] == pytest.approx([math.log(76), math.log(26)])
    assert ok(function(matrix, options={"method": "total", "target_sum": 100}))["matrix"]["values"] == [[25.0, 75.0], [75.0, 25.0]]
    assert ok(function(matrix, options={"method": "cpm"}))["matrix"]["values"] == [[250000.0, 750000.0], [750000.0, 250000.0]]


def test_normalize_zero_cell_warns_and_writes_output(tmp_path: Path) -> None:
    out = tmp_path / "norm.tsv"
    envelope = sc.normalize({"a": [0, 2], "b": [0, 2]}, options={"method": "total", "target_sum": 10, "output_path": str(out)})
    assert envelope["warnings"] and ok(envelope)["matrix"]["values"] == [[0.0, 5.0], [0.0, 5.0]]
    assert out.read_text().splitlines()[0] == "gene\tcell_1\tcell_2"
    with pytest.raises(ConfigurationError):  # never overwrites
        sc.normalize({"a": [1, 2]}, options={"output_path": str(out)})


# ---------------------------------------------------------------------------------------------
# statistics (against scipy)
# ---------------------------------------------------------------------------------------------


def test_stats_match_scipy() -> None:
    x, y = [1.0, 2.5, 3.1, 4.8, 2.2], [3.3, 5.1, 6.2, 5.8, 7.9, 6.1]
    t, df, p = _impl_stats.welch_t_test(x, y)
    ref = scipy_stats.ttest_ind(x, y, equal_var=False)
    assert t == pytest.approx(ref.statistic, rel=1e-12) and p == pytest.approx(ref.pvalue, rel=1e-9)
    ints_a, ints_b = [0, 0, 1, 3, 3, 5], [1, 2, 2, 4, 6, 6, 7]
    u, p_u = _impl_stats.mann_whitney_u(ints_a, ints_b)
    ref_u = scipy_stats.mannwhitneyu(ints_a, ints_b, method="asymptotic", alternative="two-sided")
    assert u == ref_u.statistic and p_u == pytest.approx(ref_u.pvalue, rel=1e-9)
    assert _impl_stats.fisher_exact_two_sided(30, 30, 70, 20) == pytest.approx(scipy_stats.fisher_exact([[30, 30], [70, 20]])[1], rel=1e-9)
    pvals = [0.01, 0.04, 0.03, 0.2, 0.5]
    assert _impl_stats.benjamini_hochberg(pvals) == pytest.approx(scipy_stats.false_discovery_control(pvals))
    assert _impl_stats.spearman([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]) == pytest.approx(scipy_stats.spearmanr([1, 2, 3, 4, 5], [5, 6, 7, 8, 7])[0])
    assert _impl_stats.quantile([1, 2, 3, 4, 10], 0.75) == pytest.approx(float(np.quantile([1, 2, 3, 4, 10], 0.75)))


# ---------------------------------------------------------------------------------------------
# differential features, markers, perturbations
# ---------------------------------------------------------------------------------------------

DE_MATRIX = {
    # cells: a1 a2 a3 b1 b2 b3 (library sizes are all 100 so CP10k = count * 100)
    "UP": [20, 30, 25, 5, 0, 10],
    "DOWN": [5, 5, 10, 40, 35, 45],
    "FLAT": [75, 65, 65, 55, 65, 45],
}
DE_CELLS = ["a1", "a2", "a3", "b1", "b2", "b3"]
DE_GROUPS = ["A", "A", "A", "B", "B", "B"]


def _de_oracle(gene: str, method: str, target: float):
    row = DE_MATRIX[gene]
    scaled = [v / 100.0 * target for v in row]
    first = [math.log2(v + 1) for v in scaled[:3]]
    second = [math.log2(v + 1) for v in scaled[3:]]
    if method == "welch":
        return scipy_stats.ttest_ind(first, second, equal_var=False).pvalue
    return scipy_stats.mannwhitneyu(first, second, method="asymptotic").pvalue


@pytest.mark.parametrize("function", [sc.find_differential_features, sc.find_differential_features_spatial])
@pytest.mark.parametrize("method", ["welch", "wilcoxon"])
def test_find_differential_features(function, method) -> None:
    result = ok(function(DE_MATRIX, options={"cells": DE_CELLS, "groups": DE_GROUPS, "method": method}))
    by_gene = {row["gene"]: row for row in result["genes"]}
    assert result["group1"] == "A" and result["group2"] == "B" and result["n_genes_tested"] == 3
    for gene in DE_MATRIX:
        assert by_gene[gene]["p_value"] == pytest.approx(_de_oracle(gene, method, 10000.0), rel=1e-9)
    # UP: group A mean CP10k = (2000 + 3000 + 2500)/3 = 2500; B = (500 + 0 + 1000)/3 = 500; log2fc = log2(2501) - log2(501)
    assert by_gene["UP"]["log2fc"] == pytest.approx(math.log2(2501) - math.log2(501))
    assert by_gene["UP"]["mean_group1"] == pytest.approx(2500.0) and by_gene["UP"]["pct_group2"] == pytest.approx(2 / 3)
    assert by_gene["DOWN"]["log2fc"] < 0
    qs = _impl_stats.benjamini_hochberg([by_gene[g]["p_value"] for g in by_gene])
    assert [by_gene[g]["q_value"] for g in by_gene] == pytest.approx(qs)
    assert [row["gene"] for row in result["genes"]] == sorted(by_gene, key=lambda g: by_gene[g]["p_value"])


def test_differential_features_explicit_groups_and_errors() -> None:
    labels = ["A", "A", "A", "B", "B", "B"]
    three = ["A", "A", "B", "B", "C", "C"]
    options = {"cells": DE_CELLS}
    result = ok(sc.find_differential_features(DE_MATRIX, options={**options, "groups": three, "group1": "C", "group2": "A", "method": "wilcoxon"}))
    assert (result["n_group1"], result["n_group2"]) == (2, 2)
    with pytest.raises(ConfigurationError):
        sc.find_differential_features(DE_MATRIX, options={**options, "groups": three})  # needs group1/group2
    with pytest.raises(InputFormatError):
        sc.find_differential_features(DE_MATRIX, options={**options, "groups": labels, "group1": "A", "group2": "Z"})
    with pytest.raises(InputFormatError):  # Welch needs two cells per group
        sc.find_differential_features(DE_MATRIX, options={**options, "groups": ["A", "B", "B", "B", "B", "B"]})


@pytest.mark.parametrize("function", [sc.calculate_markers, sc.calculate_markers_spatial])
def test_calculate_markers(function) -> None:
    result = ok(function(DE_MATRIX, options={"cells": DE_CELLS, "groups": DE_GROUPS, "method": "welch"}))
    markers = result["markers"]
    # only genes up in the group are reported: UP and (slightly) FLAT (raw means 68.3 vs 55) are up in A, DOWN is up in B
    assert {row["gene"] for row in markers["A"]["genes"]} == {"UP", "FLAT"}
    assert [row["gene"] for row in markers["A"]["genes"]] == sorted(["UP", "FLAT"], key=lambda g: _de_oracle(g, "welch", 10000.0))
    assert [row["gene"] for row in markers["B"]["genes"]] == ["DOWN"]
    up = next(row for row in markers["A"]["genes"] if row["gene"] == "UP")
    assert up["p_value"] == pytest.approx(_de_oracle("UP", "welch", 10000.0), rel=1e-9)
    # B versus rest is the mirror image of A versus rest: same p-value, opposite log2fc
    assert markers["B"]["genes"][0]["p_value"] == pytest.approx(_de_oracle("DOWN", "welch", 10000.0), rel=1e-9)
    both = ok(function(DE_MATRIX, options={"cells": DE_CELLS, "groups": DE_GROUPS, "only_positive": False}))["markers"]["A"]["genes"]
    assert {row["gene"] for row in both} == {"UP", "DOWN", "FLAT"}


def test_markers_skip_tiny_group_with_warning() -> None:
    envelope = sc.calculate_markers(DE_MATRIX, options={"cells": DE_CELLS, "groups": ["A", "A", "A", "B", "B", "C"]})
    assert "C" not in envelope["result"]["markers"] and any("'C' skipped" in w for w in envelope["warnings"])


PERT_MATRIX = {
    # every cell has total 100; T is the target gene
    "T": [20, 40, 10, 5],
    "OTHER": [80, 60, 90, 95],
}
PERT_CELLS = ["ctl1", "ctl2", "g1", "g2"]
PERT_LABELS = ["control", "control", "sgT", "sgT"]


def test_analyze_crispr_perturbations_knockdown() -> None:
    result = ok(sc.analyze_crispr_perturbations(PERT_MATRIX, options={"cells": PERT_CELLS, "perturbations": PERT_LABELS, "targets": {"sgT": "T"}}))
    guide = result["guides"]["sgT"]
    # control CP10k of T: 2000, 4000 -> mean 3000; guide: 1000, 500 -> mean 750; knockdown = 100 * (1 - 750/3000) = 75 %
    assert guide["mean_normalized_control"] == pytest.approx(3000.0) and guide["mean_normalized_guide"] == pytest.approx(750.0)
    assert guide["percent_knockdown"] == pytest.approx(75.0)
    assert guide["log2fc"] == pytest.approx(math.log2(751) - math.log2(3001))
    assert result["cells_per_label"] == {"control": 2, "sgT": 2}
    missing = sc.analyze_crispr_perturbations(PERT_MATRIX, options={"cells": PERT_CELLS, "perturbations": PERT_LABELS, "targets": {"sgT": "NOPE"}})
    assert ok(missing)["guides"]["sgT"]["target_in_matrix"] is False and missing["warnings"]


def test_calculate_perturbation_effect() -> None:
    result = ok(sc.calculate_perturbation_effect(PERT_MATRIX, options={"cells": PERT_CELLS, "perturbations": PERT_LABELS, "method": "wilcoxon"}))
    genes = {row["gene"]: row for row in result["effects"]["sgT"]["genes"]}
    assert result["n_control_cells"] == 2 and result["effects"]["sgT"]["n_cells"] == 2
    assert genes["T"]["log2fc"] == pytest.approx(math.log2(751) - math.log2(3001))
    first = [math.log2(v + 1) for v in (1000.0, 500.0)]
    second = [math.log2(v + 1) for v in (2000.0, 4000.0)]
    assert genes["T"]["p_value"] == pytest.approx(scipy_stats.mannwhitneyu(first, second, method="asymptotic").pvalue)
    with pytest.raises(InputFormatError):
        sc.calculate_perturbation_effect(PERT_MATRIX, options={"cells": PERT_CELLS, "perturbations": PERT_LABELS, "control": "nope"})


# ---------------------------------------------------------------------------------------------
# marker-based annotation, cell cycle, interactions
# ---------------------------------------------------------------------------------------------

CELL_TYPE_MATRIX = {
    # each cell total 100 -> with target_sum=100 the values are ln(1 + count)
    "CD3D": [40, 50, 0, 2],
    "CD3E": [30, 40, 1, 0],
    "CD79A": [1, 0, 45, 50],
    "MS4A1": [0, 2, 35, 40],
    "ACTB": [29, 8, 19, 8],
}


@pytest.mark.parametrize("function", [sc.annotate_cells, sc.annotate_cells_spatial, sc.infer_cell_type_from_markers])
def test_marker_annotation(function) -> None:
    options = {"cells": ["t1", "t2", "b1", "b2"], "target_sum": 100,
               "markers": {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1", "NOT_IN_MATRIX"]}}
    envelope = function(CELL_TYPE_MATRIX, options=options)
    result = ok(envelope)
    assert result["assignments"] == {"t1": "T cell", "t2": "T cell", "b1": "B cell", "b2": "B cell"}
    # hand score of T cell in t1: z-scores of ln(1+x) of CD3D = [ln41, ln51, ln1, ln3] ...
    values = [math.log(1 + x) for x in CELL_TYPE_MATRIX["CD3D"]]
    mean = sum(values) / 4
    sd = math.sqrt(sum((v - mean) ** 2 for v in values) / 4)
    z_d = (values[0] - mean) / sd
    values_e = [math.log(1 + x) for x in CELL_TYPE_MATRIX["CD3E"]]
    mean_e = sum(values_e) / 4
    sd_e = math.sqrt(sum((v - mean_e) ** 2 for v in values_e) / 4)
    z_e = (values_e[0] - mean_e) / sd_e
    assert result["scores"]["T cell"]["t1"] == pytest.approx((z_d + z_e) / 2)
    assert any("NOT_IN_MATRIX" in w for w in envelope["warnings"])
    # an impossible threshold leaves everything unassigned
    strict = ok(function(CELL_TYPE_MATRIX, options={**options, "min_score": 100}))
    assert set(strict["assignments"].values()) == {"Unassigned"}


def test_cell_cycle_score_hand_computed() -> None:
    # every cell has total 10 and target_sum=10, so values are ln(1 + count); n_bins=1 makes the control set
    # "all genes that are not in the scored set".
    matrix = {
        "S1": [6, 0, 1],
        "S2": [4, 0, 1],
        "M1": [0, 6, 1],
        "M2": [0, 4, 1],
        "X": [0, 0, 6],
    }
    options = {"cells": ["s_cell", "m_cell", "g1_cell"], "s_genes": ["S1", "S2"], "g2m_genes": ["M1", "M2"],
               "target_sum": 10, "n_bins": 1, "ctrl_size": 100}
    result = ok(sc.calculate_cell_cycle_score(matrix, options=options))
    ln = math.log1p
    # (ln below is log1p) S score of s_cell = mean(ln(1+6), ln(1+4)) - mean over controls (M1, M2, X) = (ln1 + ln1 + ln1)/3 = 0
    assert result["s_score"]["s_cell"] == pytest.approx((ln(6) + ln(4)) / 2)
    # S score of m_cell = 0 - mean(ln7, ln5, 0) = -(ln7 + ln5)/3
    assert result["s_score"]["m_cell"] == pytest.approx(-(ln(6) + ln(4) + 0.0) / 3)
    assert result["g2m_score"]["m_cell"] == pytest.approx((ln(6) + ln(4)) / 2)
    # g1_cell has counts 1,1,1,1,6: S score = ln2 - mean(ln2, ln2, ln7) = (ln2 - ln7)/3 < 0, same for G2M -> G1
    assert result["s_score"]["g1_cell"] == pytest.approx(ln(1) - (2 * ln(1) + ln(6)) / 3)
    assert result["phase"] == {"s_cell": "S", "m_cell": "G2M", "g1_cell": "G1"}
    assert result["phase_counts"] == {"G1": 1, "S": 1, "G2M": 1}
    with pytest.raises(ConfigurationError):
        sc.calculate_cell_cycle_score(matrix, options={k: v for k, v in options.items() if k != "s_genes"})


INTERACTION_MATRIX = {
    # cells a1 a2 (type A) and b1 b2 (type B), totals 20 and target_sum 20 -> values ln(1 + count)
    "L": [9, 3, 0, 0],
    "R": [0, 0, 5, 7],
    "F": [11, 17, 15, 13],
}


def test_cell_interactions_mean_product() -> None:
    options = {"cells": ["a1", "a2", "b1", "b2"], "cell_types": ["A", "A", "B", "B"], "pairs": [["L", "R"]], "target_sum": 20}
    result = ok(sc.calculate_cell_interactions(INTERACTION_MATRIX, options=options))
    top = result["interactions"][0]
    # sender A: mean(ln10, ln4); receiver B: mean(ln6, ln8)
    ligand = (math.log(10) + math.log(4)) / 2
    receptor = (math.log(6) + math.log(8)) / 2
    assert (top["sender"], top["receiver"]) == ("A", "B")
    assert top["score"] == pytest.approx(ligand * receptor)
    assert result["n_interactions"] == 4 and result["interactions"][-1]["score"] == pytest.approx(0.0)
    assert ok(sc.calculate_cell_interactions(INTERACTION_MATRIX, options={**options, "min_pct": 0.5}))["n_interactions"] == 1
    skipped = sc.calculate_cell_interactions(INTERACTION_MATRIX, options={**options, "pairs": [["L", "ZZZ"]]})
    assert skipped["result"]["n_interactions"] == 0 and skipped["warnings"]


def test_cell_interactions_spatial_only_counts_cells_within_radius() -> None:
    coordinates = {"a1": [0, 0], "a2": [10, 0], "b1": [1, 0], "b2": [50, 0]}
    options = {"cells": ["a1", "a2", "b1", "b2"], "cell_types": ["A", "A", "B", "B"], "pairs": [["L", "R"]], "target_sum": 20,
               "coordinates": coordinates, "radius": 2}
    result = ok(sc.calculate_cell_interactions_spatial(INTERACTION_MATRIX, options=options))
    # only (a1, b1) and (b1, a1) are within radius 2; A -> B uses ligand(a1) * receptor(b1)
    a_to_b = next(r for r in result["interactions"] if r["sender"] == "A" and r["receiver"] == "B")
    assert a_to_b["n_cell_pairs"] == 1
    assert a_to_b["score"] == pytest.approx(math.log(10) * math.log(6))
    assert result["n_interactions"] == 2  # A->B and B->A pairs exist; no A->A or B->B pair is within radius
    with pytest.raises(ConfigurationError):
        sc.calculate_cell_interactions_spatial(INTERACTION_MATRIX, options={k: v for k, v in options.items() if k != "radius"})


# ---------------------------------------------------------------------------------------------
# ambient RNA, demultiplexing, deconvolution, integration, velocity, CNV
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("function", [sc.remove_ambient, sc.remove_ambient_spatial])
def test_remove_ambient_subtracts_expected_soup(function) -> None:
    matrix = {"g1": [6, 10], "g2": [3, 5], "g3": [1, 85]}  # e1 = empty droplet (profile .6/.3/.1), c1 = real cell (total 100)
    options = {"cells": ["e1", "c1"], "empty_cells": ["e1"], "contamination": 0.1}
    result = ok(function(matrix, options=options))
    # soup = 0.1 * 100 * [0.6, 0.3, 0.1] = [6, 3, 1] -> corrected [4, 2, 84]
    assert result["matrix"]["cells"] == ["c1"]
    assert [row[0] for row in result["matrix"]["values"]] == pytest.approx([4.0, 2.0, 84.0])
    assert result["counts_removed"] == pytest.approx(10.0)
    explicit = ok(function(matrix, options={"cells": ["e1", "c1"], "ambient_profile": {"g1": 6, "g2": 3, "g3": 1}, "contamination": 0.1}))
    assert column(explicit["matrix"], "c1") == pytest.approx([4.0, 2.0, 84.0])
    assert column(explicit["matrix"], "e1") == pytest.approx([6 - 0.1 * 10 * 0.6, 3 - 0.1 * 10 * 0.3, 1 - 0.1 * 10 * 0.1])
    by_size = ok(function(matrix, options={"cells": ["e1", "c1"], "empty_max_counts": 10, "contamination": 0.1}))
    assert by_size["matrix"]["cells"] == ["c1"]
    with pytest.raises(ConfigurationError):
        function(matrix, options={"cells": ["e1", "c1"], "empty_cells": ["e1"]})  # contamination is required


@pytest.mark.parametrize("function", [sc.demultiplex, sc.demultiplex_spatial])
def test_demultiplex_rule(function) -> None:
    hashtags = {
        "H1": [20, 1, 10, 30, 0, 3],
        "H2": [1, 2, 9, 0, 50, 3],
        "H3": [0, 1, 0, 2, 10, 3],
    }
    cells = ["single", "neg", "dbl", "single2", "single3", "tie"]
    result = ok(function(hashtags, options={"cells": cells}))
    # single: top 20, second 1 -> H1; neg: top 2 < 5; dbl: 10 vs 9 (>= 5, 10 < 3 * 9); single2: second 2 < 5 -> H1;
    # single3: 50 vs 10 -> 50 >= 30 -> H2; tie: top 3 < 5 -> Negative
    assert result["assignments"] == {"single": "H1", "neg": "Negative", "dbl": "Doublet", "single2": "H1", "single3": "H2", "tie": "Negative"}
    assert result["counts"] == {"Doublet": 1, "H1": 2, "H2": 1, "Negative": 2}


@pytest.mark.parametrize("function", [sc.deconvolute_spots, sc.deconvolute_spots_spatial, tx.deconvolute])
def test_deconvolution_recovers_known_mixture(function) -> None:
    signature = {"g1": [5, 0], "g2": [5, 2], "g3": [0, 4], "g4": [0, 4]}  # types T1, T2; both columns sum to 10
    # s1 = 0.3 * T1 + 0.7 * T2 = [1.5, 1.5 + 1.4, 2.8, 2.8]; s2 = pure T1; s3 = pure T2
    mixture = {"g1": [1.5, 5.0, 0.0], "g2": [2.9, 5.0, 2.0], "g3": [2.8, 0.0, 4.0], "g4": [2.8, 0.0, 4.0]}
    options = {"cells": ["s1", "s2", "s3"], "signature": {"genes": list(signature), "cells": ["T1", "T2"], "counts": list(signature.values())}}
    result = ok(function(mixture, options=options))
    assert result["cell_types"] == ["T1", "T2"]
    assert result["fractions"]["s1"]["T1"] == pytest.approx(0.3, abs=1e-9) and result["fractions"]["s1"]["T2"] == pytest.approx(0.7, abs=1e-9)
    assert result["fractions"]["s2"] == pytest.approx({"T1": 1.0, "T2": 0.0}, abs=1e-9)
    assert result["fractions"]["s3"] == pytest.approx({"T1": 0.0, "T2": 1.0}, abs=1e-9)
    assert max(result["rmse"].values()) < 1e-9


def test_deconvolution_from_labelled_reference() -> None:
    reference = {"g1": [5, 5, 0, 0], "g2": [5, 5, 2, 2], "g3": [0, 0, 4, 4], "g4": [0, 0, 4, 4]}
    options = {"cells": ["m1"], "reference": {"genes": list(reference), "cells": ["r1", "r2", "r3", "r4"], "counts": list(reference.values())},
               "reference_labels": ["T1", "T1", "T2", "T2"]}
    mixture = {"g1": [1.5], "g2": [2.9], "g3": [2.8], "g4": [2.8]}
    result = ok(sc.deconvolute_spots(mixture, options=options))
    assert result["fractions"]["m1"]["T1"] == pytest.approx(0.3, abs=1e-9)
    assert tx.deconvolute(mixture, options={**options})["status"] == "ok"
    with pytest.raises(ConfigurationError):
        sc.deconvolute_spots(mixture, options={"cells": ["m1"]})


@pytest.mark.parametrize("function", [sc.integrate, sc.integrate_single_cell_batches, sc.integrate_spatial])
def test_integrate_batch_mean_centring(function) -> None:
    # values [1, 3 | 5, 9]: grand mean 4.5, batch means 2 and 7 -> shifts +2.5 and -2.5
    result = ok(function({"g": [1, 3, 5, 9]}, options={"batches": ["a", "a", "b", "b"], "log_normalize": False}))
    assert result["matrix"]["values"] == [[3.5, 5.5, 2.5, 6.5]]
    assert result["batches"] == {"a": 2, "b": 2}


def test_integrate_equals_ols_residuals_plus_grand_mean() -> None:
    matrix = {"g1": [10, 20, 30, 5, 8, 2], "g2": [1, 0, 4, 9, 7, 12], "g3": [50, 40, 45, 41, 39, 60]}
    batches = ["x", "x", "x", "y", "y", "y"]
    result = ok(sc.integrate(matrix, options={"batches": batches}))
    log = sc.normalize(matrix)["result"]["matrix"]["values"]
    design = np.array([[1, 0], [1, 0], [1, 0], [0, 1], [0, 1], [0, 1]], dtype=float)
    for row, corrected in zip(log, result["matrix"]["values"]):
        y = np.array(row)
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        expected = y - design @ beta + y.mean()
        assert corrected == pytest.approx(expected)
    with pytest.raises(InputFormatError):
        sc.integrate(matrix, options={"batches": ["x"] * 6})


@pytest.mark.parametrize("function", [sc.infer_velocity, sc.infer_velocity_spatial])
def test_velocity_steady_state(function) -> None:
    spliced = {"g": [1, 2, 3, 4]}
    unspliced = {"genes": ["g"], "cells": ["c1", "c2", "c3", "c4"], "counts": [[2, 4, 6, 9]]}
    options = {"cells": ["c1", "c2", "c3", "c4"], "unspliced": unspliced}
    # quantile=1 fits all cells: gamma = (2 + 8 + 18 + 36) / (1 + 4 + 9 + 16) = 64/30
    everything = ok(function(spliced, options={**options, "quantile": 1.0}))
    gamma = 64 / 30
    assert everything["gamma"]["g"] == pytest.approx(gamma)
    assert everything["matrix"]["values"][0] == pytest.approx([2 - gamma, 4 - 2 * gamma, 6 - 3 * gamma, 9 - 4 * gamma])
    # quantile=0.25: only the cell with the top s (4) or top u (9) is used -> gamma = 4*9 / 16 = 2.25
    extreme = ok(function(spliced, options={**options, "quantile": 0.25}))
    assert extreme["gamma"]["g"] == pytest.approx(2.25)
    assert extreme["matrix"]["values"][0] == pytest.approx([2 - 2.25, 4 - 4.5, 6 - 6.75, 0.0])
    with pytest.raises(InputFormatError):
        function(spliced, options={"cells": options["cells"], "unspliced": {"h": [1, 2, 3, 4]}})


@pytest.mark.parametrize("function", [sc.infer_cnv, sc.infer_cnv_spatial])
def test_infer_cnv_hand_computed(function) -> None:
    # columns all sum to 20 and target_sum=20, so CPM = count; log2(count + 1) is 2 for 3, 3 for 7 and 1 for 1
    matrix = {"a": [3, 3, 7], "b": [3, 3, 7], "c": [3, 3, 1], "d": [3, 3, 3], "z": [8, 8, 2]}
    positions = {"a": ["chr1", 100], "b": ["chr1", 200], "c": ["chr1", 300], "d": ["chr2", 50]}
    options = {"cells": ["r1", "r2", "t1"], "groups": ["ref", "ref", "test"], "reference_group": "ref", "gene_positions": positions,
               "window": 3, "target_sum": 20}
    result = ok(function(matrix, options=options))
    assert result["matrix"]["genes"] == ["a", "b", "c", "d"]  # z has no position
    # test cell relative to the reference (all 2): [1, 1, -1, 0]; moving average (window 3, shrinking at ends) on chr1 = [1, 1/3, 0]; chr2 = [0]
    # median of [1, 1/3, 0, 0] = 1/6 -> [5/6, 1/6, -1/6, -1/6]; reference cells stay 0
    assert column(result["matrix"], "t1") == pytest.approx([5 / 6, 1 / 6, -1 / 6, -1 / 6])
    assert column(result["matrix"], "r1") == pytest.approx([0, 0, 0, 0])
    assert result["cnv_score"]["t1"] == pytest.approx(1 / 3) and result["cnv_score"]["r1"] == 0.0
    with pytest.raises(ConfigurationError):
        function(matrix, options={k: v for k, v in options.items() if k != "reference_group"})
    with pytest.raises(InputFormatError):
        function(matrix, options={**options, "reference_group": "nope"})


# ---------------------------------------------------------------------------------------------
# numpy-based: PCA, clustering, trajectory, label transfer
# ---------------------------------------------------------------------------------------------

BLOBS = {
    "g1": [50, 48, 52, 2, 3, 1],
    "g2": [3, 2, 4, 40, 45, 41],
    "g3": [10, 12, 11, 9, 11, 10],
}
BLOB_CELLS = ["a1", "a2", "a3", "b1", "b2", "b3"]


def test_pca_matches_sklearn() -> None:
    decomposition = pytest.importorskip("sklearn.decomposition")
    parsed = _impl_matrix.parse_matrix(BLOBS, {"cells": BLOB_CELLS})
    points, ratio = _impl_embed.pca_embedding(parsed, {"n_pcs": 2})
    log = np.array(sc.normalize(BLOBS)["result"]["matrix"]["values"]).T
    reference = decomposition.PCA(n_components=2, svd_solver="full").fit(log)
    assert np.allclose(np.abs(points), np.abs(reference.transform(log)))
    assert np.allclose(ratio, reference.explained_variance_ratio_[:2])


def test_silhouette_and_kmeans_match_oracles() -> None:
    points = np.array([[0.0, 0.0], [0.2, 0.1], [0.1, 0.3], [5.0, 5.0], [5.2, 4.9], [4.9, 5.3]])
    labels, inertia = _impl_embed.kmeans(points, 2)
    assert labels.tolist() == [0, 0, 0, 1, 1, 1]
    assert _impl_embed.silhouette(points, labels) == pytest.approx(sklearn_metrics.silhouette_score(points, labels))
    manual = sum(((points[labels == c] - points[labels == c].mean(axis=0)) ** 2).sum() for c in (0, 1))
    assert inertia == pytest.approx(manual)
    three = np.array([[0, 0], [0, 1], [10, 0], [10, 1], [20, 0], [20, 1.0]])
    labels3, _ = _impl_embed.kmeans(three, 3)
    assert labels3.tolist() == [0, 0, 1, 1, 2, 2]


def test_cluster_separates_two_populations() -> None:
    result = ok(sc.cluster(BLOBS, options={"cells": BLOB_CELLS, "n_clusters": 2, "n_pcs": 2, "return_embedding": True}))
    labels = result["labels"]
    assert len({labels["a1"], labels["a2"], labels["a3"]}) == 1 and len({labels["b1"], labels["b2"], labels["b3"]}) == 1
    assert labels["a1"] != labels["b1"] and result["cluster_sizes"] == {"0": 3, "1": 3}
    assert len(result["embedding"]["a1"]) == 2
    auto = ok(sc.cluster(BLOBS, options={"cells": BLOB_CELLS, "n_pcs": 2}))
    assert auto["n_clusters"] == 2 and max(auto["silhouette_by_k"], key=auto["silhouette_by_k"].get) == "2"
    with pytest.raises(ConfigurationError):
        sc.cluster(BLOBS, options={"cells": BLOB_CELLS, "n_clusters": 99})


def test_spatial_smoothing_formula_and_cluster_spatial() -> None:
    points = np.array([[0.0], [10.0], [20.0], [30.0]])
    coords = [(0, 0), (1, 0), (2, 0), (3, 0)]
    # k=1: nearest neighbours (ties to the lower index) are 1, 0, 1, 2 -> smoothed = 0.5 * own + 0.5 * neighbour
    smoothed = _impl_embed.spatial_smooth(points, coords, 1, 0.5)
    assert smoothed.ravel().tolist() == [5.0, 5.0, 15.0, 25.0]
    # two spatial blocks with different expression
    matrix = {"g1": [50, 48, 52, 2, 3, 1], "g2": [3, 2, 4, 40, 45, 41], "g3": [10, 12, 11, 9, 11, 10]}
    coordinates = {"a1": [0, 0], "a2": [1, 0], "a3": [2, 0], "b1": [10, 0], "b2": [11, 0], "b3": [12, 0]}
    result = ok(sc.cluster_spatial(matrix, options={"cells": BLOB_CELLS, "coordinates": coordinates, "n_clusters": 2, "n_pcs": 2, "n_neighbors": 2}))
    labels = result["labels"]
    assert labels["a1"] == labels["a2"] == labels["a3"] != labels["b1"] and labels["b1"] == labels["b2"] == labels["b3"]
    with pytest.raises(ConfigurationError):
        sc.cluster_spatial(matrix, options={"cells": BLOB_CELLS, "n_clusters": 2})  # coordinates are required


def test_missing_numpy_gives_not_implemented(monkeypatch: pytest.MonkeyPatch) -> None:
    def no_numpy():
        raise _impl_matrix.MissingDependency("numpy")

    monkeypatch.setattr(_impl_embed, "numpy_module", no_numpy)
    envelope = sc.cluster(BLOBS, options={"cells": BLOB_CELLS, "n_clusters": 2})
    assert envelope["status"] == "not_implemented" and "numpy" in envelope["required_dependency"]


def test_dijkstra_prefers_cheaper_indirect_path() -> None:
    graph = {0: {1: 1.0, 2: 5.0}, 1: {0: 1.0, 2: 1.5}, 2: {0: 5.0, 1: 1.5}, 3: {}}
    assert _impl_embed.dijkstra(graph, 0, 4) == [0.0, 1.0, 2.5, math.inf]


def test_infer_trajectory_complete_graph_equals_normalised_euclidean_distance() -> None:
    # all columns sum to 10 and target_sum=10 -> values are ln(1 + count); with k = n - 1 the kNN graph is complete,
    # so shortest paths equal direct Euclidean distances (triangle inequality).
    matrix = {"g1": [0, 2, 5, 9], "g2": [10, 8, 5, 1]}
    cells = ["c1", "c2", "c3", "c4"]
    vectors = [(math.log1p(a), math.log1p(b)) for a, b in zip(matrix["g1"], matrix["g2"])]
    distances = [math.dist(vectors[0], v) for v in vectors]
    result = ok(sc.infer_trajectory(matrix, options={"cells": cells, "root": "c1", "n_neighbors": 3, "target_sum": 10}))
    assert [result["pseudotime"][c] for c in cells] == pytest.approx([d / max(distances) for d in distances])


def test_infer_trajectory_disconnected_cells_are_null_and_spatial_chain() -> None:
    matrix = {"g1": [0, 1, 20, 21], "g2": [20, 19, 0, 1]}
    cells = ["c1", "c2", "c3", "c4"]
    result = sc.infer_trajectory(matrix, options={"cells": cells, "root": "c1", "n_neighbors": 1, "target_sum": 21})
    assert result["result"]["pseudotime"]["c3"] is None and result["warnings"]
    # spatial: spots on a line, k=1 gives the chain c1-c2-c3; weights are expression distances
    chain = {"g1": [0, 4, 8], "g2": [8, 4, 0]}
    coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2.5, 0]}
    options = {"cells": ["c1", "c2", "c3"], "root": "c1", "coordinates": coordinates, "n_neighbors": 1, "target_sum": 8}
    vectors = [(math.log1p(a), math.log1p(b)) for a, b in zip(chain["g1"], chain["g2"])]
    d12, d23 = math.dist(vectors[0], vectors[1]), math.dist(vectors[1], vectors[2])
    result = ok(sc.infer_trajectory_spatial(chain, options=options))
    assert result["pseudotime"] == pytest.approx({"c1": 0.0, "c2": d12 / (d12 + d23), "c3": 1.0})
    with pytest.raises(InputFormatError):
        sc.infer_trajectory(matrix, options={"cells": cells, "root": "zzz"})


@pytest.mark.parametrize("function", [sc.map_cells, sc.map_cells_spatial, sc.map_single_cell_to_reference])
def test_label_transfer_majority_vote(function) -> None:
    reference = {"gA": [40, 38, 3, 2], "gB": [2, 4, 45, 41], "gC": [8, 8, 2, 7]}
    options = {
        "cells": ["q1", "q2"],
        "reference": {"genes": list(reference), "cells": ["r1", "r2", "r3", "r4"], "counts": list(reference.values())},
        "reference_labels": ["X", "X", "Y", "Y"], "n_neighbors": 3,
    }
    query = {"gA": [41, 1], "gB": [3, 44], "gC": [6, 5]}
    result = ok(function(query, options=options))
    # q1 looks like X cells: its 3 nearest reference cells are r1, r2 (X) and one Y -> X with confidence 2/3
    assert result["assignments"]["q1"]["label"] == "X" and result["assignments"]["q1"]["confidence"] == pytest.approx(2 / 3)
    assert result["assignments"]["q2"]["label"] == "Y" and result["assignments"]["q2"]["nearest_reference_cell"] in {"r3", "r4"}
    assert result["n_shared_genes"] == 3
    with pytest.raises(InputFormatError):
        function({"zz": [1, 2], "yy": [3, 4]}, options=options)  # no shared genes


# ---------------------------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------------------------


def test_plot_embedding_draws_one_circle_per_cell(tmp_path: Path) -> None:
    embedding = {"c1": [0, 0], "c2": [1, 1], "c3": [2, 0]}
    result = ok(sc.plot_embedding(embedding, options={"labels": ["a", "b", "a"], "title": "T"}))
    assert result["svg"].count("<circle") == 3 and result["n_points"] == 3 and ">T<" in result["svg"]
    assert result["svg"].count("#0072B2") >= 2 and "#D55E00" in result["svg"]  # two label colours + legend entries
    out = tmp_path / "e.svg"
    written = ok(sc.plot_embedding("cell,x,y\nc1,0,0\nc2,1,1\n", options={"output_path": str(out)}))
    assert out.read_text() == written["svg"] and written["n_points"] == 2
    with pytest.raises(ConfigurationError):
        sc.plot_embedding(embedding, options={"output_path": str(out)})
    with pytest.raises(InputFormatError):
        sc.plot_embedding({"c1": [1], "c2": [2]})


def test_plot_embedding_spatial_colours_by_value() -> None:
    result = ok(sc.plot_embedding_spatial({"s1": [0, 0], "s2": [5, 5]}, options={"values": {"s1": 0.0, "s2": 9.0}}))
    # the lowest value gets the lightest colour, the highest the darkest
    assert "#f7fbff" in result["svg"] and "#08306b" in result["svg"]


def test_plot_heatmap_and_spatial_gene_map() -> None:
    result = ok(sc.plot_heatmap(DE_MATRIX, options={"cells": DE_CELLS, "max_genes": 2}))
    # raw-count variances: DOWN 346.7, UP 140, FLAT 106.7
    assert result["genes_shown"] == ["DOWN", "UP"]
    assert result["svg"].count("<rect") == 1 + 2 * 6  # background + 2 genes x 6 cells
    spatial = ok(sc.plot_heatmap_spatial({"G": [0, 10, 5]}, options={"cells": ["s1", "s2", "s3"], "gene": "G", "coordinates": {"s1": [0, 0], "s2": [1, 0], "s3": [2, 0]}, "target_sum": 10}))
    assert spatial["svg"].count("<circle") == 3 and spatial["max_value"] == pytest.approx(math.log1p(10))
    with pytest.raises(InputFormatError):
        sc.plot_heatmap_spatial({"G": [1, 2]}, options={"gene": "H", "coordinates": [[0, 0], [1, 1]]})


# ---------------------------------------------------------------------------------------------
# reads, UMIs, record -> matrix
# ---------------------------------------------------------------------------------------------

READ1 = "ACGTACGTACGTACGT"
UMI1 = "AAACCCGGGTTT"
FASTQ = f"@r1\n{READ1}{UMI1}\n+\n{'I' * 28}\n@r2\n{READ1}{UMI1[::-1]}\n+\n{'I' * 28}\n@r3\nTTTT\n+\nIIII\n"


@pytest.mark.parametrize("function", [sc.extract_barcodes, sc.extract_barcodes_spatial])
def test_extract_barcodes(function) -> None:
    result = ok(function(FASTQ))
    assert result["barcode_counts"] == {READ1: 2} and result["n_reads"] == 3 and result["n_too_short"] == 1
    assert result["reads"][0] == {"read": "r1", "barcode": READ1}
    whitelist = ["ACGTACGTACGTACGA", "TTTTTTTTTTTTTTTT"]  # one mismatch from the observed barcode
    corrected = ok(function(FASTQ, options={"whitelist": whitelist}))
    assert corrected["barcode_counts"] == {"ACGTACGTACGTACGA": 2}
    assert corrected["whitelist"] == {"size": 2, "exact": 0, "corrected": 2, "unmatched": 0}
    unmatched = ok(function(FASTQ, options={"whitelist": ["GGGGGGGGGGGGGGGG"]}))
    assert unmatched["whitelist"]["unmatched"] == 2 and unmatched["barcode_counts"] == {}
    plain = ok(function([f"{READ1}{UMI1}"], options={"barcode_length": 4}))
    assert plain["barcode_counts"] == {"ACGT": 1}


@pytest.mark.parametrize("function", [sc.extract_umis, sc.extract_umis_spatial])
def test_extract_umis(function) -> None:
    result = ok(function(FASTQ))
    assert result["umi_counts"] == {UMI1: 1, UMI1[::-1]: 1}
    assert result["reads"][1]["umi"] == "TTTGGGCCCAAA"
    fasta = ok(function(f">x\n{READ1}\n{UMI1}\n", options={"umi_start": 16, "umi_length": 4}))
    assert fasta["umi_counts"] == {"AAAC": 1}
    assert ok(function({"r": f"{READ1}{UMI1}"}, options={"umi_length": 6, "umi_start": 20}))["umi_counts"] == {"CCGGGT": 1}


@pytest.mark.parametrize("function", [sc.correct_umis, sc.correct_umis_spatial])
def test_correct_umis_algorithms(function) -> None:
    counts = {"AAAA": 10, "AAAT": 3, "AAAC": 6, "CCCC": 2}
    # directional: AAAT (3 <= (10 + 1) / 2) is absorbed by AAAA; AAAC (6 > 5.5) is not; CCCC is far away
    directional = ok(function(counts))
    assert directional["n_unique_before"] == 4 and directional["n_unique_after"] == 3
    assert directional["corrected"] == {"AAAT": "AAAA"}
    # adjacency: AAAA absorbs both neighbours (AAAT, AAAC) in one step; cluster: same component
    for method in ("adjacency", "cluster"):
        result = ok(function(counts, options={"method": method}))
        assert result["corrected"] == {"AAAT": "AAAA", "AAAC": "AAAA"} and result["n_unique_after"] == 2
    assert ok(function(counts, options={"method": "unique"}))["n_unique_after"] == 4
    # list input counts occurrences
    assert ok(function(["AAAA"] * 4 + ["AAAT"]))["n_unique_after"] == 1
    # chained neighbours: B is absorbed by A, and C (adjacent only to B) is absorbed transitively by directional/cluster but not by adjacency
    chain = {"AAAA": 100, "AAAT": 40, "AATT": 15}
    assert ok(function(chain, options={"method": "directional"}))["n_unique_after"] == 1
    assert ok(function(chain, options={"method": "adjacency"}))["n_unique_after"] == 2
    assert ok(function(chain, options={"method": "cluster", "threshold": 1}))["n_unique_after"] == 1
    assert ok(function(chain, options={"method": "adjacency", "threshold": 2}))["n_unique_after"] == 1
    assert ok(function("AAAA\t5\nAAAT\t1\n"))["n_unique_after"] == 1


def test_build_cell_gene_matrix_counts_distinct_umis() -> None:
    records = [["c1", "g1", "AAA"], ["c1", "g1", "AAA"], ["c1", "g1", "CCC"], ["c2", "g1", "AAA"], ["c2", "g2", "GGG"]]
    result = ok(sc.build_cell_gene_matrix(records))
    assert result["matrix"]["genes"] == ["g1", "g2"] and result["matrix"]["cells"] == ["c1", "c2"]
    assert result["matrix"]["values"] == [[2.0, 1.0], [0.0, 1.0]]
    text = "barcode\tgene\tcount\nc1\tg1\t3\nc1\tg2\t1\nc2\tg1\t2\n"
    assert ok(sc.build_cell_gene_matrix(text))["matrix"]["values"] == [[3.0, 2.0], [1.0, 0.0]]
    mapped = [{"cell": "c1", "gene": "g1"}, {"cell": "c1", "gene": "g1"}]
    assert ok(sc.build_cell_gene_matrix(mapped))["matrix"]["values"] == [[2.0]]
    with pytest.raises(InputFormatError):
        sc.build_cell_gene_matrix([["c1"]])


def test_build_cell_peak_matrix_fragments_and_cut_sites() -> None:
    peaks = ["chr1:100-200", "chr1:150-300", "chr2:0-50"]
    fragments = (
        "# comment\n"
        "chr1\t90\t110\tAAA\n"      # overlaps peak 1 only (110 > 100, 90 < 150... fragment end 110 < 150)
        "chr1\t160\t170\tAAA\n"     # inside both peak 1 and peak 2
        "chr1\t250\t400\tBBB\t3\n"  # overlaps peak 2 only, weight 3
        "chr2\t60\t80\tBBB\n"       # outside the chr2 peak
        "chr1\t500\t600\tCCC\n"     # outside every peak
    )
    result = ok(sc.build_cell_peak_matrix(fragments, options={"peaks": peaks}))
    assert result["matrix"]["genes"] == peaks and result["matrix"]["cells"] == ["AAA", "BBB", "CCC"]
    assert result["matrix"]["values"] == [[2.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 0.0]]
    # cut sites: AAA 90-110 -> bases 90 and 109: 109 in peak 1; 160-170 -> 160, 169 both in peaks 1 and 2 (each counted)
    cuts = ok(sc.build_cell_peak_matrix(fragments, options={"peaks": peaks, "count": "cut_sites"}))
    assert cuts["matrix"]["values"] == [[3.0, 0.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 0.0]]
    with pytest.raises(ConfigurationError):
        sc.build_cell_peak_matrix(fragments)
    with pytest.raises(InputFormatError):
        sc.build_cell_peak_matrix("chr1\t10\n", options={"peaks": peaks})


# ---------------------------------------------------------------------------------------------
# transcriptomics
# ---------------------------------------------------------------------------------------------

BULK = {"g1": [10, 20], "g2": [30, 60], "g3": [100, 100]}


def test_bulk_normalisation_methods() -> None:
    opts = {"cells": ["s1", "s2"]}
    cpm = ok(tx.normalize(BULK, options=opts))["matrix"]["values"]
    assert cpm[0] == pytest.approx([10 / 140 * 1e6, 20 / 180 * 1e6])
    lengths = {"g1": 1000, "g2": 2000, "g3": 500}
    # TPM: rates (per kb) s1 = [10, 15, 200], s2 = [20, 30, 200]; s1 sum 225, s2 sum 250
    tpm = ok(tx.normalize(BULK, options={**opts, "method": "tpm", "gene_lengths": lengths}))["matrix"]["values"]
    assert [row[0] for row in tpm] == pytest.approx([10 / 225e-6, 15 / 225e-6, 200 / 225e-6])
    assert [row[1] for row in tpm] == pytest.approx([20 / 250e-6, 30 / 250e-6, 200 / 250e-6])
    assert sum(row[0] for row in tpm) == pytest.approx(1e6)
    # FPKM: count / (kb * millions): g1 s1 = 10 / (1 * 140e-6)
    fpkm = ok(tx.normalize(BULK, options={**opts, "method": "fpkm", "gene_lengths": lengths}))["matrix"]["values"]
    assert fpkm[0][0] == pytest.approx(10 / (1.0 * 140e-6)) and fpkm[1][1] == pytest.approx(60 / (2.0 * 180e-6))
    # DESeq2 median of ratios: geometric means 14.142, 42.426, 100; ratios s1 = [.7071, .7071, 1] -> median .7071; s2 -> 1.4142
    mor = ok(tx.normalize(BULK, options={**opts, "method": "median_of_ratios"}))
    assert mor["size_factors"]["s1"] == pytest.approx(math.sqrt(0.5)) and mor["size_factors"]["s2"] == pytest.approx(math.sqrt(2.0))
    assert mor["matrix"]["values"][0] == pytest.approx([10 / math.sqrt(0.5), 20 / math.sqrt(2.0)])
    # upper quartile: s1 q75 of [10, 30, 100] = 65, s2 of [20, 60, 100] = 80, mean 72.5
    uq = ok(tx.normalize(BULK, options={**opts, "method": "upper_quartile"}))
    assert uq["size_factors"] == pytest.approx({"s1": 65 / 72.5, "s2": 80 / 72.5})
    logged = ok(tx.normalize(BULK, options={**opts, "log2": True}))["matrix"]["values"]
    assert logged[0][0] == pytest.approx(math.log2(cpm[0][0] + 1))
    with pytest.raises(ConfigurationError):
        tx.normalize(BULK, options={"method": "tpm"})  # needs gene_lengths
    with pytest.raises(InputFormatError):
        tx.normalize({"a": [1, 0], "b": [0, 1]}, options={"method": "median_of_ratios"})  # no gene positive everywhere


def test_bulk_qc_and_detect() -> None:
    matrix = {"MT-X": [5, 0, 1], "A": [50, 3, 0], "B": [45, 0, 0], "C": [0, 7, 0]}
    cells = ["s1", "s2", "s3"]
    qc = ok(tx.qc(matrix, options={"cells": cells, "min_library_size": 5, "max_pct_mito": 4, "min_count": 1}))
    s1 = qc["samples"]["s1"]
    assert s1["library_size"] == 100 and s1["genes_detected"] == 3 and s1["pct_mito"] == pytest.approx(5.0)
    assert s1["top_gene_fraction"] == pytest.approx(0.5) and s1["top10_genes_fraction"] == pytest.approx(1.0)
    assert qc["samples"]["s2"]["library_size"] == 10 and qc["samples"]["s2"]["fraction_zero_genes"] == pytest.approx(0.5)
    assert qc["failing_samples"] == ["s1", "s3"]  # s1: 5 % mito > 4; s3: library 1 < 5 (and mito 100 %)
    detect = ok(tx.detect(matrix, options={"cells": cells, "min_count": 3, "min_samples": 2}))
    assert detect["detected_genes"] == ["A"] and detect["per_sample_detected"] == {"s1": 3, "s2": 2, "s3": 0}  # only A reaches 3 reads in two samples
    detect = ok(tx.detect(matrix, options={"cells": cells, "min_count": 3}))
    assert detect["detected_genes"] == ["MT-X", "A", "B", "C"] and detect["n_detected"] == 4
    cpm = ok(tx.detect(matrix, options={"cells": cells, "min_cpm": 400000, "min_samples": 1}))
    # CPM >= 4e5: s1 A (0.5e6) and B (0.45e6); s2 C 0.7e6; s3 MT-X 1e6
    assert cpm["detected_genes"] == ["MT-X", "A", "B", "C"]


def test_bulk_differential_expression_matches_scipy() -> None:
    cells = ["a1", "a2", "a3", "b1", "b2", "b3"]
    result = ok(tx.calculate_differential_expression(DE_MATRIX, options={"cells": cells, "groups": DE_GROUPS, "target_sum": 100}))
    by_gene = {row["gene"]: row for row in result["genes"]}
    for gene in DE_MATRIX:
        row = DE_MATRIX[gene]
        first = [math.log2(v + 1) for v in row[:3]]
        second = [math.log2(v + 1) for v in row[3:]]
        assert by_gene[gene]["p_value"] == pytest.approx(scipy_stats.ttest_ind(first, second, equal_var=False).pvalue, rel=1e-9)
    # default scaling is CPM (target 1e6) on libraries of 100
    default = ok(tx.calculate_differential_expression(DE_MATRIX, options={"cells": cells, "groups": DE_GROUPS, "method": "wilcoxon"}))
    assert default["target_sum"] == 1e6 and default["genes"][0]["gene"] in DE_MATRIX
    raw = ok(tx.calculate_differential_expression(DE_MATRIX, options={"cells": cells, "groups": DE_GROUPS, "normalization": "none"}))
    up = next(r for r in raw["genes"] if r["gene"] == "UP")
    assert up["log2fc"] == pytest.approx(math.log2(25 + 1) - math.log2(5 + 1))  # means 25 and 5 of the raw counts


def test_differential_splicing_psi_and_fisher() -> None:
    inclusion = {"e1": [10, 20, 30, 40], "e2": [5, 5, 5, 5]}
    exclusion = {"genes": ["e1", "e2"], "cells": ["s1", "s2", "s3", "s4"], "counts": [[10, 20, 10, 10], [5, 5, 5, 5]]}
    options = {"cells": ["s1", "s2", "s3", "s4"], "exclusion": exclusion, "groups": ["A", "A", "B", "B"], "min_reads": 10}
    result = ok(tx.calculate_differential_splicing(inclusion, options=options))
    events = {e["event"]: e for e in result["events"]}
    # PSI: group A samples 10/20 = .5 and 20/40 = .5; group B 30/40 = .75 and 40/50 = .8 -> mean .775
    assert events["e1"]["psi_group1"] == pytest.approx(0.5) and events["e1"]["psi_group2"] == pytest.approx(0.775)
    assert events["e1"]["delta_psi"] == pytest.approx(-0.275)
    assert (events["e1"]["inclusion_group1"], events["e1"]["exclusion_group1"], events["e1"]["inclusion_group2"], events["e1"]["exclusion_group2"]) == (30, 30, 70, 20)
    assert events["e1"]["p_value"] == pytest.approx(scipy_stats.fisher_exact([[30, 30], [70, 20]])[1], rel=1e-9)
    assert events["e2"]["p_value"] == pytest.approx(1.0) and events["e2"]["delta_psi"] == pytest.approx(0.0)
    assert events["e1"]["q_value"] == pytest.approx(min(1.0, events["e1"]["p_value"] * 2))
    lengths = ok(tx.calculate_differential_splicing(inclusion, options={**options, "inclusion_length": 2.0}))
    first = next(e for e in lengths["events"] if e["event"] == "e1")
    assert first["psi_group1"] == pytest.approx(((10 / 2) / (10 / 2 + 10) + (20 / 2) / (20 / 2 + 20)) / 2)
    sparse = ok(tx.calculate_differential_splicing(inclusion, options={**options, "min_reads": 10_000}))
    assert all(e["p_value"] is None for e in sparse["events"])
    with pytest.raises(ConfigurationError):
        tx.calculate_differential_splicing(inclusion, options={k: v for k, v in options.items() if k != "exclusion"})
    with pytest.raises(InputFormatError):
        tx.calculate_differential_splicing({"e1": [1.5, 2, 3, 4]}, options={**options, "exclusion": {"e1": [1, 2, 3, 4]}})


def test_translation_efficiency_hand_values() -> None:
    ribo = {"g1": [10], "g2": [30]}
    rna = {"g1": [20], "g2": [20]}
    result = ok(tx.calculate_translation_efficiency(ribo, options={"cells": ["s"], "rna": rna}))
    # both libraries total 40: TE = (10/40) / (20/40) = 0.5 and (30/40) / (20/40) = 1.5
    assert result["genes"]["g1"]["te"] == {"s": pytest.approx(0.5)} and result["genes"]["g2"]["mean_te"] == pytest.approx(1.5)
    logged = ok(tx.calculate_translation_efficiency(ribo, options={"cells": ["s"], "rna": rna, "log2": True}))
    assert logged["genes"]["g1"]["te"]["s"] == pytest.approx(-1.0) and logged["genes"]["g2"]["te"]["s"] == pytest.approx(math.log2(1.5))
    # a gene without RNA reads has no TE; libraries are normalised separately (ribo total 40, RNA total 20 here)
    zero = ok(tx.calculate_translation_efficiency(ribo, options={"cells": ["s"], "rna": {"g1": [0], "g2": [20]}}))
    assert zero["genes"]["g1"]["te"]["s"] is None and zero["genes"]["g2"]["te"]["s"] == pytest.approx((30 / 40) / (20 / 20))
    named = ok(tx.calculate_translation_efficiency({"g1": [10, 5], "g2": [30, 5]}, options={"cells": ["x", "y"], "rna": {"genes": ["g1", "g2"], "cells": ["y", "x"], "counts": [[5, 20], [5, 20]]}}))
    assert named["genes"]["g1"]["te"]["x"] == pytest.approx(0.5) and named["genes"]["g1"]["te"]["y"] == pytest.approx(1.0)
    with pytest.raises(ConfigurationError):
        tx.calculate_translation_efficiency(ribo)


def test_compare_correlations_match_scipy() -> None:
    matrix = {"g1": [1, 2], "g2": [5, 3], "g3": [9, 8], "g4": [0, 1], "g5": [20, 25]}
    other = {"g1": [2], "g2": [4], "g3": [7], "g4": [1], "g5": [30]}
    x = [math.log2(v + 1) for v in (1, 5, 9, 0, 20)]
    y = [math.log2(v + 1) for v in (2, 4, 7, 1, 30)]
    pearson = ok(tx.compare(matrix, options={"cells": ["a", "b"], "other": other}))
    assert pearson["correlation"]["a"]["cell_1"] == pytest.approx(scipy_stats.pearsonr(x, y)[0])
    spearman = ok(tx.compare(matrix, options={"cells": ["a", "b"], "other": other, "method": "spearman"}))
    assert spearman["correlation"]["a"]["cell_1"] == pytest.approx(scipy_stats.spearmanr(x, y)[0])
    within = ok(tx.compare(matrix, options={"cells": ["a", "b"], "log2": False}))
    assert within["correlation"]["a"]["a"] == pytest.approx(1.0)
    assert within["correlation"]["a"]["b"] == pytest.approx(scipy_stats.pearsonr([1, 5, 9, 0, 20], [2, 3, 8, 1, 25])[0])
    with pytest.raises(InputFormatError):
        tx.compare(matrix, options={"cells": ["a", "b"], "other": {"zz": [1]}})


def test_bulk_plot_heatmap_and_scatter() -> None:
    heat = ok(tx.plot(BULK, options={"cells": ["s1", "s2"]}))
    assert heat["svg"].count("<rect") == 1 + 3 * 2
    scatter = ok(tx.plot(BULK, options={"cells": ["s1", "s2"], "kind": "scatter", "x": "s1", "y": "s2"}))
    assert scatter["svg"].count("<circle") == 3
    assert scatter["pearson_r"] == pytest.approx(scipy_stats.pearsonr([math.log2(v + 1) for v in (10, 30, 100)], [math.log2(v + 1) for v in (20, 60, 100)])[0])
    with pytest.raises(ConfigurationError):
        tx.plot(BULK, options={"kind": "scatter"})
    with pytest.raises(InputFormatError):
        tx.plot(BULK, options={"kind": "scatter", "x": "nope", "y": "cell_2"})


GTF = (
    "#!genome-build test\n"
    'chr17\tHAVANA\tgene\t7661779\t7687550\t.\t-\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53"; gene_biotype "protein_coding";\n'
    'chr17\tHAVANA\texon\t7661779\t7661957\t.\t-\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53";\n'
    'chr1\tENSEMBL\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "ABC1"; gene_biotype "lncRNA";\n'
)


def test_annotate_from_gtf(tmp_path: Path) -> None:
    gtf = tmp_path / "a.gtf"
    gtf.write_text(GTF)
    result = ok(tx.annotate(["TP53", "ENSG00000000001.7", "ENSG00000141510", "NOPE"], options={"gtf": str(gtf)}))
    tp53 = result["annotations"]["TP53"]
    assert (tp53["chromosome"], tp53["start"], tp53["end"], tp53["strand"], tp53["biotype"]) == ("chr17", 7661779, 7687550, "-", "protein_coding")
    assert result["annotations"]["ENSG00000141510"]["gene_name"] == "TP53"
    assert result["annotations"]["ENSG00000000001.7"]["gene_name"] == "ABC1"  # version suffix of the query is ignored
    assert result["annotations"]["NOPE"] is None and result["not_found"] == ["NOPE"] and result["n_found"] == 3
    import gzip
    zipped = tmp_path / "a.gtf.gz"
    zipped.write_bytes(gzip.compress(GTF.encode()))
    assert ok(tx.annotate({"ABC1": [1]}, options={"gtf": str(zipped)}))["annotations"]["ABC1"]["start"] == 100
    with pytest.raises(ConfigurationError):
        tx.annotate(["TP53"])
    with pytest.raises(InputFormatError):
        tx.annotate(["TP53"], options={"gtf": str(tmp_path / "missing.gtf")})
    empty = tmp_path / "e.gtf"
    empty.write_text("# nothing\n")
    with pytest.raises(InputFormatError):
        tx.annotate(["TP53"], options={"gtf": str(empty)})


# ---------------------------------------------------------------------------------------------
# external programs (fake executables on PATH)
# ---------------------------------------------------------------------------------------------


def _fake_tool(directory: Path, name: str, body: str) -> Path:
    path = directory / name
    path.write_text(f"#!{sys.executable}\nimport json, os, sys\nargs = sys.argv[1:]\nopen(os.path.join(os.path.dirname(__file__), '{name}.argv'), 'w').write(json.dumps(args))\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _argv(directory: Path, name: str) -> list[str]:
    return json.loads((directory / f"{name}.argv").read_text())


@pytest.fixture
def tool_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    directory = tmp_path / "bin"
    directory.mkdir()
    monkeypatch.setenv("PATH", str(directory))
    return directory


def test_run_featurecounts_builds_command_and_parses_table(tool_dir: Path, tmp_path: Path) -> None:
    _fake_tool(tool_dir, "featureCounts", r"""
out = args[args.index('-o') + 1]
open(out, 'w').write('# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\nG1\tc\t1\t9\t+\t9\t5\nG2\tc\t1\t9\t+\t9\t7\n')
open(out + '.summary', 'w').write('x')
""")
    bam, gtf = tmp_path / "a.bam", tmp_path / "a.gtf"
    bam.write_bytes(b"x")
    gtf.write_text(GTF)
    out = tmp_path / "counts.txt"
    result = ok(tx.run_featurecounts(str(bam), options={"annotation": str(gtf), "output": str(out), "strand": 2, "paired": True, "threads": 3}))
    assert _argv(tool_dir, "featureCounts") == ["-a", str(gtf), "-o", str(out), "-T", "3", "-s", "2", "-t", "exon", "-g", "gene_id", "-p", "--countReadPairs", str(bam)]
    assert result["n_genes"] == 2 and result["samples"] == ["a.bam"] and result["assigned_counts"] == {"a.bam": 12.0}
    # quantify dispatches BAM input to featureCounts
    out2 = tmp_path / "counts2.txt"
    assert ok(tx.quantify([str(bam)], options={"annotation": str(gtf), "output": str(out2)}))["method"] == "featurecounts"
    with pytest.raises(ConfigurationError):  # output already exists
        tx.run_featurecounts(str(bam), options={"annotation": str(gtf), "output": str(out)})
    with pytest.raises(InputFormatError):
        tx.run_featurecounts(str(tmp_path / "nope.bam"), options={"annotation": str(gtf), "output": str(tmp_path / "x.txt")})


def test_run_salmon_and_kallisto(tool_dir: Path, tmp_path: Path) -> None:
    _fake_tool(tool_dir, "salmon", r"""
out = args[args.index('-o') + 1]
os.makedirs(out)
open(os.path.join(out, 'quant.sf'), 'w').write('Name\tLength\tEffectiveLength\tTPM\tNumReads\nT1\t100\t80\t900000\t90\nT2\t200\t180\t100000\t10\n')
""")
    _fake_tool(tool_dir, "kallisto", r"""
out = args[args.index('-o') + 1]
os.makedirs(out)
open(os.path.join(out, 'abundance.tsv'), 'w').write('target_id\tlength\teff_length\test_counts\ttpm\nT1\t100\t80\t50\t250000\nT2\t200\t180\t150\t750000\n')
""")
    r1, r2, index = tmp_path / "r1.fq", tmp_path / "r2.fq", tmp_path / "idx"
    r1.write_text("@a\nACGT\n+\nIIII\n")
    r2.write_text("@a\nACGT\n+\nIIII\n")
    index.mkdir()
    out = tmp_path / "salmon_out"
    salmon = ok(tx.run_salmon_quantification([str(r1), str(r2)], options={"index": str(index), "output_dir": str(out), "threads": 2}))
    assert _argv(tool_dir, "salmon") == ["quant", "-i", str(index), "-l", "A", "-o", str(out), "-p", "2", "--validateMappings", "-1", str(r1), "-2", str(r2)]
    assert salmon["n_transcripts"] == 2 and salmon["total_tpm"] == pytest.approx(1e6) and salmon["total_reads_assigned"] == 100.0
    assert salmon["top_transcripts"][0]["Name"] == "T1"
    kindex = tmp_path / "k.idx"
    kindex.write_text("x")
    kout = tmp_path / "kallisto_out"
    kallisto = ok(tx.run_kallisto_quantification(str(r1), options={"index": str(kindex), "output_dir": str(kout), "fragment_length": 200, "fragment_sd": 20}))
    assert _argv(tool_dir, "kallisto") == ["quant", "-i", str(kindex), "-o", str(kout), "-t", "1", "--single", "-l", "200.0", "-s", "20.0", str(r1)]
    assert kallisto["total_est_counts"] == 200.0 and kallisto["top_transcripts"][0]["target_id"] == "T2"
    with pytest.raises(ConfigurationError):  # single-end kallisto needs fragment statistics
        tx.run_kallisto_quantification(str(r1), options={"index": str(kindex), "output_dir": str(tmp_path / "k2")})
    # quantify prefers salmon for FASTQ input
    assert ok(tx.quantify(str(r1), options={"index": str(index), "output_dir": str(tmp_path / "q1")}))["method"] == "salmon"
    assert ok(tx.quantify(str(r1), options={"method": "kallisto", "index": str(kindex), "output_dir": str(tmp_path / "q2"), "fragment_length": 180, "fragment_sd": 30}))["method"] == "kallisto"
    with pytest.raises(ConfigurationError):
        tx.quantify(str(r1), options={"method": "featurecounts"})


def test_tool_failure_and_align_commands(tool_dir: Path, tmp_path: Path) -> None:
    _fake_tool(tool_dir, "STAR", r"""
prefix = args[args.index('--outFileNamePrefix') + 1]
open(prefix + 'Aligned.sortedByCoord.out.bam', 'wb').write(b'BAM')
""")
    _fake_tool(tool_dir, "hisat2", r"""
open(args[args.index('-S') + 1], 'w').write('@HD\n')
""")
    r1 = tmp_path / "r1.fq.gz"
    r1.write_bytes(b"x")
    genome = tmp_path / "genome"
    genome.mkdir()
    star = ok(tx.align(str(r1), options={"index": str(genome), "output_prefix": str(tmp_path / "s_"), "threads": 4}))
    assert star["aligner"] == "star" and star["alignment_bytes"] == 3
    assert _argv(tool_dir, "STAR") == ["--runThreadN", "4", "--genomeDir", str(genome), "--readFilesIn", str(r1), "--outFileNamePrefix", str(tmp_path / "s_"),
                                        "--outSAMtype", "BAM", "SortedByCoordinate", "--readFilesCommand", "zcat"]
    prefix = tmp_path / "hs"
    (tmp_path / "hs.1.ht2").write_text("x")
    hisat = ok(tx.align(str(r1), options={"aligner": "hisat2", "index": str(prefix), "output_prefix": str(tmp_path / "h")}))
    assert hisat["aligner"] == "hisat2" and _argv(tool_dir, "hisat2") == ["-p", "1", "-x", str(prefix), "-U", str(r1), "-S", str(tmp_path / "h.sam")]
    _fake_tool(tool_dir, "salmon", "sys.stderr.write('boom')\nsys.exit(3)")
    index = tmp_path / "idx"
    index.mkdir()
    from dnakit.exceptions import BackendExecutionError

    with pytest.raises(BackendExecutionError, match="boom"):
        tx.run_salmon_quantification(str(r1), options={"index": str(index), "output_dir": str(tmp_path / "fail")})


# ---------------------------------------------------------------------------------------------
# honest not_implemented
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "call",
    [
        lambda: sc.remove_doublets(MAPPING),
        lambda: sc.remove_doublets_spatial(MAPPING),
        lambda: tx.predict_protein_function("MKT"),
        lambda: tx.predict_protein_localization("MKT"),
    ],
)
def test_models_and_simulation_stay_not_implemented(call) -> None:
    envelope = call()
    assert envelope["status"] == "not_implemented" and envelope["required_dependency"] and "result" not in envelope


def test_external_programs_missing_from_path_are_not_implemented(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    reads, bam, gtf, index = tmp_path / "r.fq", tmp_path / "a.bam", tmp_path / "a.gtf", tmp_path / "idx"
    for path in (reads, bam, gtf):
        path.write_text("x")
    index.mkdir()
    calls = {
        "align": lambda: tx.align(str(reads), options={"index": str(index), "output_prefix": str(tmp_path / "o_")}),
        "run_featurecounts": lambda: tx.run_featurecounts(str(bam), options={"annotation": str(gtf), "output": str(tmp_path / "c.txt")}),
        "run_salmon_quantification": lambda: tx.run_salmon_quantification(str(reads), options={"index": str(index), "output_dir": str(tmp_path / "s")}),
        "run_kallisto_quantification": lambda: tx.run_kallisto_quantification(str(reads), options={"index": str(reads), "output_dir": str(tmp_path / "k"), "fragment_length": 1, "fragment_sd": 1}),
        "quantify": lambda: tx.quantify(str(reads), options={}),
    }
    for name, call in calls.items():
        envelope = call()
        assert envelope["status"] == "not_implemented", name
        assert "PATH" in envelope["required_dependency"] and "result" not in envelope, name


# ---------------------------------------------------------------------------------------------
# invalid input
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "function, value, options",
    [
        (sc.normalize, None, {}),
        (sc.normalize, {}, {}),
        (sc.normalize, {"g": [1, 2], "h": [1]}, {}),  # ragged
        (sc.normalize, {"g": [1, -2]}, {}),  # negative counts
        (sc.normalize, {"g": [1, "x"]}, {}),  # non-numeric
        (sc.normalize, {"g": [1, float("nan")]}, {}),
        (sc.normalize, "gene,c1,c1\ng,1,2\n", {}),  # duplicate cell names
        (sc.normalize, [[1, 2]], {"genes": ["a", "b"]}),  # wrong number of names
        (sc.normalize, {"g": [1, 2]}, {"method": "nope"}),
        (sc.normalize, {"g": [1, 2]}, {"target_sum": -1}),
        (sc.normalize, {"g": [1, 2]}, {"orientation": "sideways"}),
        (sc.normalize, {"g": [1, 2]}, {"max_matrix_elements": 1}),
        (sc.qc, {"g": [1, 2]}, {"max_pct_mito": 500}),
        (sc.calculate_markers, DE_MATRIX, {}),  # groups are required
        (sc.calculate_markers, DE_MATRIX, {"groups": ["A", "B"]}),  # wrong number of labels
        (sc.annotate_cells, DE_MATRIX, {"markers": {"T": ["nope"]}}),
        (sc.annotate_cells, DE_MATRIX, {}),
        (sc.correct_umis, ["AAAA", "AAA"], {}),
        (sc.correct_umis, ["AXAA"], {}),
        (sc.correct_umis, [], {}),
        (sc.extract_barcodes, ["ACGZ"], {}),
        (sc.extract_barcodes, "@r\nACGT\n+\n", {}),
        (sc.build_cell_gene_matrix, [], {}),
        (sc.remove_ambient, DE_MATRIX, {"contamination": 1.5}),
        (sc.demultiplex, {"H1": [1, 2]}, {}),
        (sc.infer_cnv, DE_MATRIX, {"gene_positions": [1, 2]}),
        (tx.normalize, BULK, {"method": "cpm", "log2": "yes"}),
        (tx.detect, BULK, {"min_samples": 5}),
        (tx.compare, {"g": [1, 2]}, {"method": "kendall"}),
    ],
)
def test_invalid_inputs_raise_dnakit_errors(function, value, options) -> None:
    with pytest.raises(DNAKitError) as caught:
        function(value, options=options)
    assert caught.value.code and caught.value.hint


def test_options_must_be_a_mapping() -> None:
    with pytest.raises(ConfigurationError):
        sc.normalize(MAPPING, options=["method"])  # type: ignore[arg-type]


def test_public_signatures_are_unchanged() -> None:
    import inspect

    for module in (sc, tx):
        for name in module.__all__:
            function = getattr(module, name)
            if not hasattr(function, "__wrapped__"):
                continue
            params = list(inspect.signature(function.__wrapped__).parameters)
            if name in {"aggregate_pseudobulk", "calculate_cell_qc_metrics", "calculate_gene_qc_metrics", "filter_genes_by_min_cells", "rank_marker_genes",
                        "read_anndata_summary", "select_highly_variable_genes", "summarize_mudata", "train_scvi_model"}:
                continue
            assert params == ["input", "options"], name
            assert function.__doc__ and function.__doc__.strip().splitlines()[0].endswith("."), name
