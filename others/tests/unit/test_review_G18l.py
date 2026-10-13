"""Unit tests for the G18l review wrappers around cooler, cooltools and pairtools.

Every wrapped function is tested against a hand calculation or a direct upstream call on a
small synthetic cooler. The tests run offline and write their inputs into ``tmp_path``.
"""

from __future__ import annotations

import ctypes
import warnings

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("pandas")
np = pytest.importorskip("numpy")


def _cooler_fixture(tmp_path, name="synthetic.cool"):
    """Two chromosomes (2000 and 1500 bp) at 100 bp, with a planted enriched pixel, balanced."""

    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    pd = pytest.importorskip("pandas")
    rng = np.random.default_rng(0)
    chromsizes = pd.Series({"chr1": 2000, "chr2": 1500})
    bins = cooler.binnify(chromsizes, 100)
    rows = []
    for i in range(len(bins)):
        for j in range(i, len(bins)):
            distance = abs(i - j)
            same = bins.chrom[i] == bins.chrom[j]
            count = rng.poisson(6.0 / (1 + distance) if same else 0.5)
            if distance == 0:
                count += 3
            if count > 0:
                rows.append((i, j, int(count)))
    rows.append((2, 9, 40))
    pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"]).sort_values(
        ["bin1_id", "bin2_id"]
    )
    path = tmp_path / name
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
    cooler.balance_cooler(cooler.Cooler(str(path)), store=True)
    return path


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


# --- cooler -----------------------------------------------------------------------------------


def test_infer_cooler_binsize_uniform_and_non_uniform(tmp_path):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    from dnakit.epigenome.cooler_review_functions import infer_cooler_binsize

    path = _cooler_fixture(tmp_path)
    result = infer_cooler_binsize(path, verbose=False)
    assert result["binsize"] == 100
    assert result["uniform"] is True
    assert result["bin_count"] == 35

    bins = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr1"], "start": [0, 100, 150], "end": [100, 150, 250]}
    )
    uneven = tmp_path / "uneven.cool"
    cooler.create_cooler(
        str(uneven),
        bins,
        pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [1]}),
        ensure_sorted=True,
    )
    result = infer_cooler_binsize(uneven, verbose=False)
    assert result["binsize"] is None
    assert result["uniform"] is False


def test_create_single_cell_cooler_hand_counts(tmp_path):
    cooler = pytest.importorskip("cooler")
    from dnakit.epigenome.cooler_review_functions import (
        create_single_cell_cooler,
        list_cooler_groups,
        list_scool_cells,
    )

    pairs = _write(
        tmp_path / "cell.pairs",
        "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
        "r1\tchr1\t100\tchr1\t300\t+\t-\tUU\n"
        "r2\tchr1\t100\tchr1\t300\t+\t-\tUU\n"
        "r3\tchr2\t10\tchr2\t20\t+\t-\tUU\n"
        "r4\tchr1\t500\tchr2\t100\t+\t-\tUU\n",
    )
    output = tmp_path / "cells.scool"
    result = create_single_cell_cooler(
        output,
        {"cellA": pairs, "cellB": pairs},
        {"chr1": 1000, "chr2": 800},
        100,
        verbose=False,
    )
    assert result["bin_count"] == 18
    assert result["cell_count"] == 2
    cell_a = next(c for c in result["cells"] if c["name"] == "cellA")
    # chr1 bins 1 and 3 hold two pairs; chr2 bin 10 twice; chr1 bin 5 with chr2 bin 11
    assert cell_a["pair_count"] == 4
    assert cell_a["nonzero_pixels"] == 3

    table = cooler.Cooler(f"{output}::/cells/cellA").pixels()[:]
    observed = sorted(zip(table["bin1_id"], table["bin2_id"], table["count"], strict=True))
    assert observed == [(1, 3, 2), (5, 11, 1), (10, 10, 1)]

    groups = list_cooler_groups(output, verbose=False)["groups"]
    assert groups == ["/cells/cellA", "/cells/cellB"]
    cells = list_scool_cells(output, verbose=False)
    assert [c["name"] for c in cells["cells"]] == ["cellA", "cellB"]
    assert cooler.api.list_coolers(str(output)) == groups


def test_list_cooler_groups_matches_upstream_on_plain_cooler(tmp_path):
    cooler = pytest.importorskip("cooler")
    from dnakit.epigenome.cooler_review_functions import list_cooler_groups

    path = _cooler_fixture(tmp_path)
    result = list_cooler_groups(path, verbose=False)
    assert result["groups"] == sorted(cooler.api.list_coolers(str(path)))
    assert result["groups"] == ["/"]


def test_list_scool_cells_rejects_plain_cooler(tmp_path):
    from dnakit.epigenome.cooler_review_functions import list_scool_cells

    path = _cooler_fixture(tmp_path)
    with pytest.raises(ConfigurationError):
        list_scool_cells(path, verbose=False)


def test_create_single_cell_cooler_rejects_out_of_range_pair(tmp_path):
    from dnakit.epigenome.cooler_review_functions import create_single_cell_cooler

    pairs = _write(
        tmp_path / "bad.pairs",
        "r1\tchr1\t5000\tchr1\t10\t+\t-\tUU\n",
    )
    with pytest.raises(ConfigurationError):
        create_single_cell_cooler(
            tmp_path / "bad.scool", {"cellA": pairs}, {"chr1": 1000}, 100, verbose=False
        )


# --- cooltools: track, view and bin utilities -------------------------------------------------


def test_align_track_to_cooler_bins_counts_finite_weights(tmp_path):
    pytest.importorskip("cooltools")
    cooler = pytest.importorskip("cooler")
    from dnakit.epigenome.cooltools_review_functions import align_track_to_cooler_bins

    path = _cooler_fixture(tmp_path)
    track = _write(
        tmp_path / "track.bedgraph",
        "chr1\t0\t100\t1.0\nchr1\t100\t200\t2.0\nchr1\t200\t300\t3.0\n"
        "chr2\t0\t100\t0.5\nchr2\t100\t200\t0.6\n",
    )
    result = align_track_to_cooler_bins(track, path, verbose=False)
    bins = cooler.Cooler(str(path)).bins()[:]
    weights = {(str(row.chrom), int(row.start)): row.weight for row in bins.itertuples(index=False)}
    track_bins = [("chr1", 0), ("chr1", 100), ("chr1", 200), ("chr2", 0), ("chr2", 100)]
    expected_valued = sum(1 for key in track_bins if np.isfinite(weights[key]))
    assert result["bin_count"] == 35
    assert result["valued_bin_count"] == expected_valued
    value_by_bin = {(row["chrom"], row["start"]): row["value"] for row in result["bins"]}
    if np.isfinite(weights[("chr1", 100)]):
        assert value_by_bin[("chr1", 100)] == pytest.approx(2.0)


def test_make_track_view_hand_extent(tmp_path):
    from dnakit.epigenome.cooltools_review_functions import make_track_view

    track = _write(
        tmp_path / "track.bedgraph",
        "chr1\t0\t100\t1.0\nchr1\t200\t300\t3.0\nchr2\t50\t120\t0.5\n",
    )
    result = make_track_view(track, verbose=False)
    regions = {row["chrom"]: (row["start"], row["end"]) for row in result["regions"]}
    assert regions == {"chr1": (0, 300), "chr2": (50, 120)}
    assert result["region_count"] == 2


def test_mask_track_bad_bins_matches_weight_table(tmp_path):
    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import mask_track_bad_bins

    path = _cooler_fixture(tmp_path)
    track = _write(
        tmp_path / "track.bedgraph",
        "chr1\t0\t100\t1.0\nchr1\t100\t200\t2.0\nchr2\t0\t100\t3.0\n",
    )
    result = mask_track_bad_bins(track, path, verbose=False)
    bins = cooler.Cooler(str(path)).bins()[:]
    masked_bins = 0
    for chrom, start in [("chr1", 0), ("chr1", 100), ("chr2", 0)]:
        row = bins[(bins["chrom"].astype(str) == chrom) & (bins["start"] == start)]
        if not np.isfinite(row["weight"].to_numpy()[0]):
            masked_bins += 1
    assert result["row_count"] == 3
    assert result["masked_count"] == masked_bins


def test_assign_paired_features_to_view_names_each_side(tmp_path):
    from dnakit.epigenome.cooltools_review_functions import assign_paired_features_to_view

    features = _write(
        tmp_path / "features.bedpe",
        "chr1\t100\t150\tchr1\t700\t750\nchr1\t500\t550\tchr2\t50\t80\n",
    )
    view = _write(
        tmp_path / "view.bed",
        "chr1\t0\t400\tA\nchr1\t600\t1000\tB\nchr2\t0\t800\tX\n",
    )
    result = assign_paired_features_to_view(features, view, verbose=False)
    first, second = result["features"]
    assert (first["region1"], first["region2"]) == ("A", "B")
    assert second["region1"] is None
    assert second["region2"] == "X"
    assert result["assigned_count"] == 1


def test_assign_bins_to_regions_hand_spans():
    from dnakit.epigenome.cooltools_review_functions import assign_bins_to_regions

    result = assign_bins_to_regions(
        [0, 4, 5, 11, 12],
        [
            {"name": "a", "bin_start": 0, "bin_end": 5},
            {"name": "b", "bin_start": 5, "bin_end": 12},
        ],
        verbose=False,
    )
    assert [row["region"] for row in result["assignments"]] == ["a", "a", "b", "b", None]
    assert result["assigned_count"] == 4


def test_assign_bins_to_regions_rejects_overlap():
    from dnakit.epigenome.cooltools_review_functions import assign_bins_to_regions

    with pytest.raises(ConfigurationError):
        assign_bins_to_regions(
            [1],
            [
                {"name": "a", "bin_start": 0, "bin_end": 5},
                {"name": "b", "bin_start": 4, "bin_end": 9},
            ],
            verbose=False,
        )


def test_expand_features_to_bins_hand_calculation(tmp_path):
    from dnakit.epigenome.cooltools_review_functions import expand_features_to_bins

    features = _write(tmp_path / "features.bed", "chr1\t300\t500\nchr2\t100\t200\n")
    result = expand_features_to_bins(features, flank=100, resolution=50, verbose=False)
    spans = [(row["chrom"], row["start"], row["end"]) for row in result["features"]]
    # centre 400 -> floor(400/50)=8, flank 100//50=2 -> bins 6..10, edges 300 and 550
    # centre 150 -> floor(150/50)=3 -> bins 1..5, edges 50 and 300
    assert spans == [("chr1", 300, 550), ("chr2", 50, 300)]


# --- cooltools: pairwise sums, expected and observed/expected ----------------------------------


def test_pairwise_diagonal_sums_counts_diagonals(tmp_path):
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import pairwise_diagonal_sums

    path = _cooler_fixture(tmp_path)
    result = pairwise_diagonal_sums(path, verbose=False)
    # chr1 has 20 bins (diagonals 0..19), chr2 has 15 bins (diagonals 0..14)
    assert result["row_count"] == 35
    chr1_rows = [row for row in result["rows"] if row["region1"] == "chr1"]
    assert len(chr1_rows) == 20


def test_pairwise_block_sums_one_trans_pair(tmp_path):
    cooltools_api = pytest.importorskip("cooltools.api.expected")
    cooler = pytest.importorskip("cooler")
    from cooltools.lib.common import make_cooler_view

    from dnakit.epigenome.cooltools_review_functions import pairwise_block_sums

    path = _cooler_fixture(tmp_path)
    result = pairwise_block_sums(path, verbose=False)
    assert result["row_count"] == 1
    clr = cooler.Cooler(str(path))
    direct = cooltools_api.blocksum_pairwise(clr, make_cooler_view(clr), clr_weight_name="weight")
    row = result["rows"][0]
    assert row["region1"] == "chr1" and row["region2"] == "chr2"
    assert row["count.sum"] == pytest.approx(float(direct["count.sum"].iloc[0]))
    assert row["n_valid"] == int(direct["n_valid"].iloc[0])


def test_expected_full_table_row_count_and_fast_variant(tmp_path):
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import expected_full_table

    path = _cooler_fixture(tmp_path)
    slow = expected_full_table(path, verbose=False)
    fast = expected_full_table(path, fast=True, verbose=False)
    # 20 + 15 cis diagonals and one trans block
    assert slow["row_count"] == 36
    assert fast["row_count"] == 36
    assert [row["expected"] for row in slow["rows"][2:4]] == pytest.approx(
        [row["expected"] for row in fast["rows"][2:4]]
    )


def test_observed_over_expected_pixels_is_balanced_over_expected(tmp_path):
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import observed_over_expected_pixels

    path = _cooler_fixture(tmp_path)
    result = observed_over_expected_pixels(path, verbose=False)
    checked = 0
    for row in result["pixels"]:
        if row["expected"] is None or row["balanced"] is None or row["expected"] == 0:
            continue
        assert row["oe"] == pytest.approx(row["balanced"] / row["expected"], rel=1e-9)
        checked += 1
    assert checked > 0


def test_smooth_contact_decay_matches_upstream(tmp_path):
    pytest.importorskip("cooltools")
    cooler = pytest.importorskip("cooler")
    import cooltools.api.expected as expected_api
    from cooltools.lib.common import make_cooler_view

    from dnakit.epigenome.cooltools_review_functions import smooth_contact_decay

    path = _cooler_fixture(tmp_path)
    clr = cooler.Cooler(str(path))
    view = make_cooler_view(clr)
    cvd = expected_api.expected_cis(
        clr, view_df=view, intra_only=True, smooth=False, clr_weight_name="weight"
    )
    direct = expected_api.genomewide_smooth_cvd(
        cvd, cols=expected_api.DEFAULT_BALANCED_CVD_COLS, suffix=".smoothed"
    )
    result = smooth_contact_decay(path, mode="genomewide", verbose=False)
    by_key = {(row["region1"], row["region2"], row["dist"]): row for row in result["rows"]}
    for _, line in direct.iterrows():
        key = (line["region1"], line["region2"], int(line["dist"]))
        value = line["balanced.avg.smoothed"]
        if np.isfinite(value):
            assert by_key[key]["balanced.avg.smoothed"] == pytest.approx(float(value))
    per_region = smooth_contact_decay(path, mode="per_region", verbose=False)
    assert per_region["row_count"] == result["row_count"]


def test_call_chromatin_dots_matches_direct_upstream(tmp_path):
    pytest.importorskip("cooltools")
    cooler = pytest.importorskip("cooler")
    import cooltools.api.dotfinder as dotfinder
    import cooltools.api.expected as expected_api
    from cooltools.lib.common import make_cooler_view
    from cooltools.lib.numutils import get_kernel

    from dnakit.epigenome.cooltools_review_functions import call_chromatin_dots

    path = _cooler_fixture(tmp_path)
    clr = cooler.Cooler(str(path))
    view = make_cooler_view(clr)
    expected = expected_api.expected_cis(
        clr, view_df=view, intra_only=True, smooth=False, clr_weight_name="weight"
    )
    kernels = {
        kind: get_kernel(3, 1, kind) for kind in ("donut", "vertical", "horizontal", "lowleft")
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = dotfinder.dots(
            clr,
            expected,
            expected_value_col="balanced.avg",
            clr_weight_name="weight",
            view_df=view,
            kernels=kernels,
            n_lambda_bins=40,
            nproc=1,
        )
    result = call_chromatin_dots(path, kernel_half_width=3, donut_inner_half_width=1, verbose=False)
    assert result["dot_count"] == len(direct)
    with pytest.raises(ConfigurationError):
        call_chromatin_dots(path, kernel_half_width=1, verbose=False)


# --- cooltools: matrices, saddle strength, decoys and eigenvectors ----------------------------


def test_saddle_compartment_strength_hand_value():
    from dnakit.epigenome.cooltools_review_functions import saddle_compartment_strength

    sums = [[1, 2, 3, 4], [2, 5, 6, 7], [3, 6, 8, 9], [4, 7, 9, 10]]
    counts = [[2, 2, 2, 2]] * 4
    result = saddle_compartment_strength(sums, counts, verbose=False)
    first = result["strengths"][0]
    # k = 1: intra = (S00 + S33) / (C00 + C33) = 11 / 4; inter = (S03 + S30) / (C03 + C30) = 8 / 4
    assert first == {"k": 1, "strength": pytest.approx(1.375)}
    assert len(result["strengths"]) == 3


def test_adaptive_coarsegrain_matrix_matches_upstream():
    pytest.importorskip("cooltools")
    from cooltools.lib import numutils

    from dnakit.epigenome.cooltools_review_functions import adaptive_coarsegrain_matrix

    rng = np.random.default_rng(3)
    matrix = rng.random((9, 9))
    counts = np.full((9, 9), 2.0)
    result = adaptive_coarsegrain_matrix(
        matrix.tolist(), counts.tolist(), cutoff=3, max_levels=2, min_shape=4, verbose=False
    )
    direct = numutils.adaptive_coarsegrain(matrix, counts, cutoff=3, max_levels=2, min_shape=4)
    assert np.allclose(np.array(result["matrix"], dtype=float), direct)


def test_replace_cis_with_decoys_changes_only_masked_block():
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import replace_cis_with_decoys

    rng = np.random.default_rng(4)
    base = rng.random((12, 12))
    matrix = (base + base.T).tolist()
    mask = [[1 if i < 6 and j < 6 else 0 for j in range(12)] for i in range(12)]
    result = replace_cis_with_decoys(matrix, mask, seed=0, verbose=False)
    changed = np.array(result["matrix"], dtype=float)
    original = np.array(matrix)
    assert result["replaced_pixels"] == 36
    assert not np.allclose(changed[:6, :6], original[:6, :6])
    assert np.array_equal(changed[6:, 6:], original[6:, 6:])
    again = replace_cis_with_decoys(matrix, mask, seed=0, verbose=False)
    assert again["matrix"] == result["matrix"]
    with pytest.raises(ConfigurationError):
        replace_cis_with_decoys(matrix, [[2] * 12] * 12, verbose=False)


def _trans_matrix():
    rng = np.random.default_rng(5)
    base = rng.random((12, 12))
    return (base + base.T).tolist()


def test_dense_trans_eigenvectors_is_repeatable_and_matches_upstream():
    pytest.importorskip("cooltools")
    from cooltools.api import eigdecomp

    from dnakit.epigenome.cooltools_review_functions import (
        _isolated_random,
        dense_trans_eigenvectors,
    )

    matrix = _trans_matrix()
    first = dense_trans_eigenvectors(matrix, [0, 6, 12], n_eigs=2, verbose=False)
    second = dense_trans_eigenvectors(matrix, [0, 6, 12], n_eigs=2, verbose=False)
    assert first["eigenvalues"] == second["eigenvalues"]
    assert first["eigenvectors"] == second["eigenvectors"]
    with _isolated_random(0, eigen=True):
        direct_values, _ = eigdecomp.trans_eig(
            np.array(matrix), [0, 6, 12], n_eigs=2, perc_top=99.95, perc_bottom=1
        )
    assert np.allclose(np.abs(first["eigenvalues"]), np.abs(direct_values), rtol=1e-6)
    with pytest.raises(ConfigurationError):
        dense_trans_eigenvectors(matrix, [0, 5, 11], verbose=False)


def test_dense_trans_eigenvectors_restores_numpy_state():
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import dense_trans_eigenvectors

    np.random.seed(7)
    reference = np.random.get_state()[1].copy()
    np.random.seed(7)
    dense_trans_eigenvectors(_trans_matrix(), [0, 6, 12], n_eigs=2, verbose=False)
    after = np.random.get_state()[1]
    assert np.array_equal(reference, after)


def test_dense_trans_eigenvectors_restores_c_random_stream():
    pytest.importorskip("cooltools")
    from dnakit.epigenome.cooltools_review_functions import dense_trans_eigenvectors

    libc = ctypes.CDLL(None)
    libc.srandom.argtypes = [ctypes.c_uint]
    libc.random.restype = ctypes.c_long
    libc.srandom(5)
    stream = [libc.random() for _ in range(4)]
    libc.srandom(5)
    libc.random()
    dense_trans_eigenvectors(_trans_matrix(), [0, 6, 12], n_eigs=2, verbose=False)
    tail = [libc.random() for _ in range(3)]
    assert tail == stream[1:]


# --- pairtools --------------------------------------------------------------------------------

_SAM = (
    "@HD\tVN:1.6\tSO:queryname\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "@SQ\tSN:chr2\tLN:800\n"
    "r1\t65\tchr1\t101\t60\t10M\t=\t101\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
    "r1\t129\tchr2\t301\t60\t10M\t=\t301\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
    "r2\t65\tchr2\t51\t60\t10M\t=\t51\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
    "r2\t129\tchr2\t201\t60\t10M\t=\t201\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
)


def test_sam_to_pairs_hand_pairs(tmp_path):
    pytest.importorskip("pairtools")
    pytest.importorskip("pysam")
    from dnakit.epigenome.pairtools_review_functions import sam_to_pairs

    sam = _write(tmp_path / "reads.sam", _SAM)
    output = tmp_path / "reads.pairs"
    result = sam_to_pairs(sam, output, verbose=False)
    assert result["pair_count"] == 2
    assert result["chromosome_count"] == 2
    body = [line for line in output.read_text().splitlines() if not line.startswith("#")]
    assert body == ["r1\tchr1\t101\tchr2\t301\t+\t+\tUU", "r2\tchr2\t51\tchr2\t201\t+\t+\tUU"]
    header = [line for line in output.read_text().splitlines() if line.startswith("#chromsize")]
    assert header == ["#chromsize: chr1 1000", "#chromsize: chr2 800"]


def test_filter_pairs_by_coverage_keeps_and_removes(tmp_path):
    pytest.importorskip("pairtools")
    from dnakit.epigenome.pairtools_review_functions import filter_pairs_by_coverage

    pairs = _write(
        tmp_path / "cov.pairs",
        "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
        "r1\tchr1\t100\tchr1\t300\t+\t-\tUU\n"
        "r2\tchr1\t100\tchr1\t300\t+\t-\tUU\n"
        "r3\tchr1\t100\tchr1\t300\t+\t-\tUU\n"
        "r4\tchr2\t10\tchr2\t20\t+\t-\tUU\n"
        "r5\tchr1\t500\tchr2\t100\t+\t-\tUU\n",
    )
    everything = filter_pairs_by_coverage(pairs, tmp_path / "all.pairs", max_cov=100, verbose=False)
    assert everything["kept_count"] == 5
    assert everything["removed_count"] is None
    strict = filter_pairs_by_coverage(
        pairs,
        tmp_path / "strict.pairs",
        method="max",
        max_dist=500,
        max_cov=2,
        removed_path=tmp_path / "removed.pairs",
        verbose=False,
    )
    assert strict["kept_count"] + strict["removed_count"] == 5
    assert strict["removed_count"] >= 1


def test_mark_sam_duplicate_sets_flag_1024_and_dd_tag():
    pytest.importorskip("pairtools")
    from dnakit.epigenome.pairtools_review_functions import mark_sam_duplicate

    line = "r1\t99\tchr1\t101\t60\t10M\t=\t201\t110\tACGTACGTAC\tIIIIIIIIII\tYt:Z:UU"
    result = mark_sam_duplicate(line, verbose=False)
    assert result["flag_before"] == 99
    assert result["flag_after"] == 1123
    assert result["duplicate_flag_set"] is True
    assert result["sam"].endswith("Yt:Z:DD")
    again = mark_sam_duplicate(result["sam"], verbose=False)
    assert again["flag_after"] == 1123
    with pytest.raises(ConfigurationError):
        mark_sam_duplicate("r1\t99\tchr1", verbose=False)


def test_chromsizes_from_sam_header_hand_values(tmp_path):
    pytest.importorskip("pairtools")
    pytest.importorskip("pysam")
    from dnakit.epigenome.pairtools_review_functions import chromsizes_from_sam_header

    sam = _write(tmp_path / "reads.sam", _SAM)
    result = chromsizes_from_sam_header(sam, verbose=False)
    assert result["chromosomes"] == [
        {"name": "chr1", "length": 1000},
        {"name": "chr2", "length": 800},
    ]
    assert result["chromosome_count"] == 2
