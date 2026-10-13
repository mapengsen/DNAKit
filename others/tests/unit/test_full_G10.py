"""Tests for the G10 wrappers: cooler, cooltools and pairtools file and matrix functions.

Each wrapper is compared with a hand calculation or with a direct call to the upstream
package on the same input. The UCSC download path is tested with a mocked HTTP client.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest

from dnakit.epigenome import cooler_functions as CF
from dnakit.epigenome import cooltools_functions as TF
from dnakit.epigenome import pairtools_functions as PF
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

# Two chromosomes, binsize 100: chr1 has bins 0-9, chr2 has bins 10-15.
# Pixels: (0,1)x2, (0,2)x1, (5,5)x1, (1,10)x1, (10,10)x1, (12,14)x1 -> 6 nonzero, 7 contacts.
TINY_SIZES = {"chr1": 1000, "chr2": 600}
TINY_PAIRS = """## pairs format v1.0
#chromsize: chr1 1000
#chromsize: chr2 600
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1\tchr1\t10\tchr1\t150\t+\t-\tUU
r2\tchr1\t20\tchr1\t160\t+\t-\tUU
r3\tchr1\t15\tchr1\t250\t+\t-\tUU
r4\tchr1\t510\tchr1\t520\t+\t-\tUU
r5\tchr1\t120\tchr2\t30\t+\t-\tUR
r6\tchr2\t40\tchr2\t90\t+\t-\tUU
r7\tchr2\t200\tchr2\t420\t-\t+\tUU
"""

SAM_TEXT = (
    "@HD\tVN:1.6\tSO:queryname\n"
    "@SQ\tSN:chrA\tLN:1000\n"
    "@SQ\tSN:chrB\tLN:800\n"
    "r1\t65\tchrA\t101\t60\t10M\t=\t101\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
    "r1\t129\tchrB\t301\t60\t10M\t=\t301\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10\tXS:i:0\n"
)


def _random_pairs(seed: int = 3, contacts: int = 60) -> str:
    """Deterministic pairs for two chromosomes with local contacts and a few trans contacts."""

    import numpy as np

    rng = np.random.default_rng(seed)
    sizes = {"chr1": 1000, "chr2": 800}
    lines = [
        "## pairs format v1.0",
        "#chromsize: chr1 1000",
        "#chromsize: chr2 800",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    for index in range(contacts):
        chrom = "chr1" if index % 3 else "chr2"
        length = sizes[chrom]
        pos1 = int(rng.integers(1, length))
        pos2 = min(length - 1, pos1 + int(rng.exponential(120)) + 1)
        other = chrom
        if index % 11 == 0:
            other = "chr1" if chrom == "chr2" else "chr2"
            pos2 = int(rng.integers(1, sizes[other]))
        lines.append(f"r{index}\t{chrom}\t{pos1}\t{other}\t{pos2}\t+\t-\tUU")
    return "\n".join(lines) + "\n"


@pytest.fixture(scope="module")
def workdir(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("g10")


def _cool(workdir: Path, name: str, text: str = TINY_PAIRS, sizes=TINY_SIZES) -> Path:
    pairs = workdir / f"{name}.pairs"
    pairs.write_text(text, encoding="utf-8")
    output = workdir / f"{name}.cool"
    build_contact_matrix(pairs, output, chromsizes=sizes, binsize=100, verbose=False)
    return output


def _balanced(workdir: Path, name: str, text: str = TINY_PAIRS, sizes=TINY_SIZES) -> Path:
    output = _cool(workdir, name, text, sizes)
    balance_contact_matrix(output, min_nnz=1, verbose=False)
    return output


@pytest.fixture(scope="module")
def tiny_cool(workdir: Path) -> Path:
    pytest.importorskip("cooler")
    return _cool(workdir, "tiny")


@pytest.fixture(scope="module")
def tiny_balanced(workdir: Path) -> Path:
    pytest.importorskip("cooltools")
    return _balanced(workdir, "tiny_bal")


@pytest.fixture(scope="module")
def random_balanced(workdir: Path) -> Path:
    pytest.importorskip("cooltools")
    pairs_text = _random_pairs()
    return _balanced(workdir, "random_bal", pairs_text, {"chr1": 1000, "chr2": 800})


# --- cooler_functions -------------------------------------------------------------------


def test_read_chrom_sizes_matches_upstream_and_keeps_file_order(workdir: Path) -> None:
    cooler_util = pytest.importorskip("cooler.util")
    path = workdir / "sizes.tsv"
    path.write_text("chr1\t1000\nchr2\t600\n", encoding="utf-8")
    result = CF.read_chrom_sizes(path, all_names=True, verbose=False)
    direct = cooler_util.read_chromsizes(str(path), all_names=True)
    assert result["chromosomes"] == [
        {"name": "chr1", "length": 1000},
        {"name": "chr2", "length": 600},
    ]
    assert [row["name"] for row in result["chromosomes"]] == list(direct.index)
    assert result["status"] == "completed" and result["operation"] == "read_chrom_sizes"


def test_read_chrom_sizes_name_filter_and_empty_filter(workdir: Path) -> None:
    pytest.importorskip("cooler.util")
    path = workdir / "sizes_filter.tsv"
    path.write_text("chr1\t1000\nchr2\t600\nchrUn\t50\n", encoding="utf-8")
    only_one = CF.read_chrom_sizes(path, name_patterns=["^chr1$"], verbose=False)
    assert only_one["chromosomes"] == [{"name": "chr1", "length": 1000}]
    with pytest.raises(ConfigurationError, match="name filter"):
        CF.read_chrom_sizes(path, name_patterns=["^chrZ$"], verbose=False)
    with pytest.raises(ConfigurationError, match="Invalid name pattern"):
        CF.read_chrom_sizes(path, name_patterns=["("], verbose=False)


def test_fetch_chrom_sizes_uses_ucsc_url_and_mocked_client(monkeypatch) -> None:
    pytest.importorskip("cooler.util")
    import pandas.io.common as pandas_common

    class _Response(io.BytesIO):
        def __init__(self, data: bytes) -> None:
            super().__init__(data)
            self.headers: dict[str, str] = {}

    seen: list[str] = []
    payload = gzip.compress(b"chr1\t248956422\nchr2\t242193529\nchrUn_gl000220\t161802\n")

    def fake_urlopen(request, *args, **kwargs):
        seen.append(getattr(request, "full_url", request))
        return _Response(payload)

    monkeypatch.setattr(pandas_common, "urlopen", fake_urlopen)
    result = CF.fetch_chrom_sizes("hg38", verbose=False)
    assert seen == ["http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz"]
    assert result["database"] == "hg38"
    assert result["chromosomes"] == [
        {"name": "chr1", "length": 248956422},
        {"name": "chr2", "length": 242193529},
    ]
    assert result["chromosome_count"] == 2


def test_fetch_chrom_sizes_network_error_is_reported(monkeypatch) -> None:
    pytest.importorskip("cooler.util")
    import pandas.io.common as pandas_common

    def failing_urlopen(request, *args, **kwargs):
        raise OSError("network unreachable")

    monkeypatch.setattr(pandas_common, "urlopen", failing_urlopen)
    with pytest.raises(BackendUnavailableError, match="could not be downloaded"):
        CF.fetch_chrom_sizes("hg38", verbose=False)
    with pytest.raises(ConfigurationError, match="UCSC assembly name"):
        CF.fetch_chrom_sizes("hg 38", verbose=False)


def test_parse_genomic_region_hand_values() -> None:
    pytest.importorskip("cooler.util")
    result = CF.parse_genomic_region("chr1:100-200", verbose=False)
    assert (result["chrom"], result["start"], result["end"], result["length"]) == (
        "chr1",
        100,
        200,
        100,
    )
    whole = CF.parse_genomic_region("chr2", chromsizes=TINY_SIZES, verbose=False)
    assert (whole["start"], whole["end"]) == (0, 600)
    with pytest.raises(ConfigurationError, match="not valid"):
        CF.parse_genomic_region("chr1:100-2000", chromsizes=TINY_SIZES, verbose=False)


def test_query_bed_intervals_by_hand_and_upstream_order(workdir: Path) -> None:
    pytest.importorskip("cooler.util")
    bed = workdir / "intervals.bed"
    bed.write_text("chr1\t0\t150\tx\nchr1\t200\t300\ty\nchr2\t10\t20\tz\n", encoding="utf-8")
    overlap = CF.query_bed_intervals(bed, "chr1:100-250", chromsizes=TINY_SIZES, verbose=False)
    assert overlap["intervals"] == [
        {"chrom": "chr1", "start": 0, "end": 150},
        {"chrom": "chr1", "start": 200, "end": 300},
    ]
    empty = CF.query_bed_intervals(bed, "chr1:500-600", chromsizes=TINY_SIZES, verbose=False)
    assert empty["interval_count"] == 0
    unsorted = workdir / "unsorted.bed"
    unsorted.write_text("chr1\t200\t300\nchr1\t0\t150\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="sorted by start"):
        CF.query_bed_intervals(unsorted, "chr1", chromsizes=TINY_SIZES, verbose=False)


def test_rename_cooler_chromosomes_changes_only_the_copy(tiny_cool: Path, workdir: Path) -> None:
    cooler = pytest.importorskip("cooler")
    output = workdir / "renamed.cool"
    result = CF.rename_cooler_chromosomes(tiny_cool, output, {"chr1": "1"}, verbose=False)
    assert result["chromosomes_after"] == ["1", "chr2"]
    assert cooler.Cooler(str(output)).chromnames == ["1", "chr2"]
    assert cooler.Cooler(str(tiny_cool)).chromnames == ["chr1", "chr2"]
    with pytest.raises(ConfigurationError, match="not chromosomes"):
        CF.rename_cooler_chromosomes(tiny_cool, workdir / "bad.cool", {"chrZ": "1"}, verbose=False)


def test_coarsen_contact_matrix_doubles_binsize_and_keeps_contacts(
    tiny_cool: Path, workdir: Path
) -> None:
    cooler = pytest.importorskip("cooler")
    output = workdir / "coarse.cool"
    result = CF.coarsen_contact_matrix(tiny_cool, output, factor=2, verbose=False)
    assert (result["input_binsize"], result["output_binsize"]) == (100, 200)
    assert result["output_total_contacts"] == 7  # hand: 7 contacts before and after
    direct = cooler.Cooler(str(output))
    assert int(direct.info["sum"]) == 7
    assert int(direct.info["nnz"]) == result["output_nonzero_pixels"]


def test_merge_contact_matrices_doubles_contacts_of_identical_bins(
    workdir: Path, tiny_cool: Path
) -> None:
    pytest.importorskip("cooler")
    twin = _cool(workdir, "twin")
    output = workdir / "merged.cool"
    result = CF.merge_contact_matrices([tiny_cool, twin], output, verbose=False)
    assert result["input_total_contacts"] == 14
    assert result["output_total_contacts"] == 14
    assert result["output_nonzero_pixels"] == 6


def test_zoomify_contact_matrix_writes_each_resolution(tiny_cool: Path, workdir: Path) -> None:
    cooler = pytest.importorskip("cooler")
    output = workdir / "zoom.mcool"
    result = CF.zoomify_contact_matrix(
        tiny_cool, output, resolutions=[100, 200, 400], verbose=False
    )
    assert [level["binsize"] for level in result["resolutions"]] == [100, 200, 400]
    assert all(level["total_contacts"] == 7 for level in result["resolutions"])
    direct = cooler.Cooler(f"{output}::/resolutions/200")
    assert direct.binsize == 200


def test_load_fasta_lengths_matches_sequence_lengths(workdir: Path) -> None:
    pytest.importorskip("cooler.util")
    pytest.importorskip("pyfaidx")
    fasta = workdir / "genome.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 600 + "\n", encoding="utf-8")
    result = CF.load_fasta_lengths(fasta, ["chr2", "chr1"], verbose=False)
    assert result["sequences"] == [
        {"name": "chr2", "length": 600},
        {"name": "chr1", "length": 1000},
    ]
    with pytest.raises(ConfigurationError, match="not in the FASTA"):
        CF.load_fasta_lengths(fasta, ["chrQ"], verbose=False)


def test_read_cooler_pixels_matches_upstream_rows(tiny_cool: Path) -> None:
    cooler = pytest.importorskip("cooler")
    result = CF.read_cooler_pixels(tiny_cool, limit=3, verbose=False)
    direct = cooler.Cooler(str(tiny_cool)).pixels()[0:3]
    assert result["nonzero_pixels"] == 6  # hand: six nonzero pixels in the tiny matrix
    assert result["columns"] == ["bin1_id", "bin2_id", "count"]
    assert [row["bin2_id"] for row in result["pixels"]] == list(direct["bin2_id"])
    joined = CF.read_cooler_pixels(tiny_cool, limit=1, join=True, verbose=False)
    assert joined["pixels"][0]["chrom1"] == "chr1" and joined["pixels"][0]["start1"] == 0


def test_read_cooler_chroms_lists_names_and_lengths(tiny_cool: Path) -> None:
    pytest.importorskip("cooler")
    result = CF.read_cooler_chroms(tiny_cool, verbose=False)
    assert result["chromosomes"] == [
        {"name": "chr1", "length": 1000},
        {"name": "chr2", "length": 600},
    ]


def test_cooler_wrappers_refuse_existing_output(tiny_cool: Path) -> None:
    pytest.importorskip("cooler")
    with pytest.raises(ConfigurationError, match="will not be overwritten"):
        CF.coarsen_contact_matrix(tiny_cool, tiny_cool, factor=2, verbose=False)


# --- cooltools_functions ----------------------------------------------------------------


def test_calculate_expected_trans_matches_direct_and_hand(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    result = TF.calculate_expected_trans(random_balanced, verbose=False)
    direct = cooltools.expected_trans(
        cooler.Cooler(str(random_balanced)), clr_weight_name="weight", nproc=1
    )
    assert result["pair_count"] == len(direct) == 1
    row = result["pairs"][0]
    assert (row["region1"], row["region2"]) == ("chr1", "chr2")
    assert row["count.sum"] == pytest.approx(float(direct["count.sum"].iloc[0]))
    assert row["balanced.sum"] == pytest.approx(float(direct["balanced.sum"].iloc[0]))


def test_calculate_expected_trans_needs_balanced_weights(tiny_cool: Path) -> None:
    pytest.importorskip("cooltools")
    with pytest.raises(ConfigurationError, match="balance_contact_matrix"):
        TF.calculate_expected_trans(tiny_cool, verbose=False)


def test_calculate_compartments_matches_direct_eigs(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    from cooltools.lib import common as cooltools_common

    result = TF.calculate_compartments(random_balanced, n_eigs=2, verbose=False)
    clr = cooler.Cooler(str(random_balanced))
    with TF._reproducible_eigen():
        values, vectors = cooltools.eigs_cis(
            clr, view_df=cooltools_common.make_cooler_view(clr), n_eigs=2, clr_weight_name="weight"
        )
    assert len(result["eigenvalues"]) == 2  # hand: one row per chromosome
    assert len(result["eigenvectors"]) == int(clr.info["nbins"])  # hand: one row per bin
    first = result["eigenvalues"][0]
    assert first["eigval1"] == pytest.approx(float(values.loc[0, "eigval1"]))
    # eigenvector signs follow the wrapper convention, so compare magnitudes
    assert abs(result["eigenvectors"][0]["E1"]) == pytest.approx(abs(float(vectors.iloc[0]["E1"])))


def test_calculate_trans_compartments_matches_direct_eigs(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    result = TF.calculate_trans_compartments(random_balanced, n_eigs=1, verbose=False)
    with TF._reproducible_eigen():
        values, vectors = cooltools.eigs_trans(
            cooler.Cooler(str(random_balanced)), n_eigs=1, clr_weight_name="weight"
        )
    assert result["eigenvalues"][0]["eigval1"] == pytest.approx(float(values["eigval1"].iloc[0]))
    assert len(result["eigenvectors"]) == len(vectors)
    assert abs(result["eigenvectors"][0]["E1"]) == pytest.approx(abs(float(vectors.iloc[0]["E1"])))


def test_calculate_saddle_has_padded_matrix_and_direct_sum(
    random_balanced: Path, workdir: Path
) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    from cooltools.lib import common as cooltools_common

    lines = [f"chr1\t{i * 100}\t{(i + 1) * 100}\t{-1 + 2 * i / 9:.6f}" for i in range(10)]
    lines += [f"chr2\t{i * 100}\t{(i + 1) * 100}\t{1 - 2 * i / 7:.6f}" for i in range(8)]
    track = workdir / "track_saddle.bed"
    track.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = TF.calculate_saddle(
        random_balanced, track, n_bins=3, qrange=(0.0, 1.0), min_diag=1, verbose=False
    )
    assert result["matrix_size"] == 5  # hand: three bins plus one padding row/column per side
    assert result["interaction_sum"][0] == [0.0] * 5  # padding row is empty for qrange (0, 1)
    import pandas as pd

    table = pd.read_csv(track, sep="\t", header=None, names=["chrom", "start", "end", "value"])
    clr = cooler.Cooler(str(random_balanced))
    view = cooltools_common.make_cooler_view(clr)
    expected = cooltools.expected_cis(clr, view_df=view, clr_weight_name="weight", nproc=1)
    direct_sum, _ = cooltools.saddle(
        clr,
        expected,
        table,
        "cis",
        3,
        qrange=(0.0, 1.0),
        view_df=view,
        clr_weight_name="weight",
        expected_value_col="balanced.avg",
        min_diag=1,
    )
    assert result["interaction_sum"][2][2] == pytest.approx(float(direct_sum[2, 2]))


def test_sample_contacts_exact_count_and_restores_numpy_state(
    random_balanced: Path, workdir: Path
) -> None:
    cooler = pytest.importorskip("cooler")
    np = pytest.importorskip("numpy")
    before = np.random.get_state()
    result = TF.sample_contacts(random_balanced, workdir / "exact40.cool", count=40, verbose=False)
    after = np.random.get_state()
    assert result["total_contacts_after"] == 40  # exact sampling keeps the requested count
    assert result["method"] == "exact"
    assert int(cooler.Cooler(str(workdir / "exact40.cool")).info["sum"]) == 40
    assert before[0] == after[0] and (before[1] == after[1]).all() and before[2:] == after[2:]


def test_sample_contacts_binomial_is_reproducible_with_seed(
    random_balanced: Path, workdir: Path
) -> None:
    pytest.importorskip("cooler")
    first = TF.sample_contacts(
        random_balanced, workdir / "bin_a.cool", fraction=0.5, seed=5, verbose=False
    )
    second = TF.sample_contacts(
        random_balanced, workdir / "bin_b.cool", fraction=0.5, seed=5, verbose=False
    )
    assert first["total_contacts_after"] == second["total_contacts_after"]
    assert first["method"] == "binomial"
    with pytest.raises(ConfigurationError, match="exactly one"):
        TF.sample_contacts(random_balanced, workdir / "neither.cool", verbose=False)


def test_calculate_virtual_4c_matches_direct_profile(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    result = TF.calculate_virtual_4c(random_balanced, "chr1:0-100", verbose=False)
    direct = cooltools.virtual4c(
        cooler.Cooler(str(random_balanced)), ("chr1", 0, 100), clr_weight_name="weight", nproc=1
    )
    assert result["bin_count"] == len(direct)
    assert result["viewpoint"] == {"chrom": "chr1", "start": 0, "end": 100}
    first = next(row for row in result["profile"] if row["start"] == 100)
    assert first["balanced"] == pytest.approx(float(direct["balanced"].iloc[1]))


def test_calculate_directionality_index_matches_direct_call(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    result = TF.calculate_directionality_index(random_balanced, window_bp=300, verbose=False)
    direct = cooltools.directionality(
        cooler.Cooler(str(random_balanced)), window_bp=300, balance="weight"
    )
    assert result["bin_count"] == len(direct)
    value_column = "directionality_index_300"
    index = next(i for i, value in enumerate(direct[value_column]) if value == value)
    assert result["bins"][index][value_column] == pytest.approx(
        float(direct[value_column].iloc[index])
    )
    assert result["bins"][0][value_column] is None  # NaN is returned as None


def test_pileup_features_mean_matches_direct_pileup(random_balanced: Path, workdir: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    np = pytest.importorskip("numpy")
    import pandas as pd
    from cooltools.lib import common as cooltools_common

    features = workdir / "features_pileup.bed"
    features.write_text("chr1\t300\t400\nchr1\t500\t600\nchr2\t200\t300\n", encoding="utf-8")
    result = TF.pileup_features(random_balanced, features, flank_bp=200, verbose=False)
    clr = cooler.Cooler(str(random_balanced))
    table = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr2"], "start": [300, 500, 200], "end": [400, 600, 300]}
    )
    stack = cooltools.pileup(
        clr,
        table,
        view_df=cooltools_common.make_cooler_view(clr),
        expected_df=None,
        flank=200,
        clr_weight_name="weight",
        nproc=1,
    )
    with np.errstate(all="ignore"):
        mean = np.nanmean(stack, axis=0)
    assert result["feature_count"] == 3
    assert result["snippet_size"] == 5
    finite = [(r, c) for r in range(5) for c in range(5) if mean[r, c] == mean[r, c]]
    row, col = finite[0]
    assert result["mean_snippet"][row][col] == pytest.approx(float(mean[row, col]))
    assert result["mean_snippet"][0][0] is None  # NaN is returned as None


def test_calculate_logbin_expected_matches_direct_chain(random_balanced: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    from cooltools.api import expected as expected_api
    from cooltools.lib import common as cooltools_common

    result = TF.calculate_logbin_expected(
        random_balanced, bins_per_order_magnitude=4, min_nvalid=1, min_count=1, verbose=False
    )
    clr = cooler.Cooler(str(random_balanced))
    expected = cooltools.expected_cis(
        clr, view_df=cooltools_common.make_cooler_view(clr), clr_weight_name="weight", nproc=1
    )
    binned, slopes, edges = expected_api.logbin_expected(
        expected,
        summary_name="balanced.sum",
        bins_per_order_magnitude=4,
        min_nvalid=1,
        min_count=1,
    )
    assert result["binned_rows"] == len(binned)
    assert result["bin_edges"] == [int(edge) for edge in edges]
    assert len(result["slopes"]) == len(slopes)


def test_rearrange_contact_matrix_reorders_and_keeps_contacts(
    random_balanced: Path, workdir: Path
) -> None:
    cooler = pytest.importorskip("cooler")
    view = workdir / "rearrange_view.tsv"
    view.write_text(
        "chrom\tstart\tend\tname\tnew_chrom\tstrand\n"
        "chr1\t0\t1000\tchr1\tchrX\t+\n"
        "chr2\t0\t800\tchr2\tchrX\t-\n",
        encoding="utf-8",
    )
    result = TF.rearrange_contact_matrix(
        random_balanced, view, workdir / "rearranged.cool", verbose=False
    )
    direct = cooler.Cooler(str(workdir / "rearranged.cool"))
    assert result["chromosomes_after"] == ["chrX"]
    assert int(direct.info["nbins"]) == int(cooler.Cooler(str(random_balanced)).info["nbins"])
    assert result["output_total_contacts"] == int(cooler.Cooler(str(random_balanced)).info["sum"])


def test_align_windows_to_bins_hand_value_and_direct_call() -> None:
    cooltools = pytest.importorskip("cooltools")
    result = TF.align_windows_to_bins(
        100, ["chr1", "chr1"], [450, 800], flank_bp=200, verbose=False
    )
    first = result["windows"][0]
    # hand: centre 450 +- 200 is [250, 650]; the bin edges give [200, 700), bins 2 to 6
    assert (first["start"], first["end"], first["lo"], first["hi"]) == (200, 700, 2, 7)
    direct = cooltools.api.snipping.make_bin_aligned_windows(
        100, ["chr1", "chr1"], [450, 800], flank_bp=200
    )
    assert [row["end"] for row in result["windows"]] == list(direct["end"])


def test_assign_features_to_view_by_hand(workdir: Path) -> None:
    pytest.importorskip("cooltools")
    features = workdir / "assign_features.bed"
    features.write_text("chr1\t300\t400\nchr2\t10\t20\n", encoding="utf-8")
    view = workdir / "assign_view.bed"
    view.write_text("chr1\t0\t1000\tchr1\n", encoding="utf-8")
    result = TF.assign_features_to_view(features, view, verbose=False)
    assert result["features"][0]["region"] == "chr1"
    assert result["features"][1]["region"] is None
    assert result["unassigned_count"] == 1


def test_make_chromosome_view_matches_cooler_sizes(tiny_cool: Path) -> None:
    pytest.importorskip("cooltools")
    result = TF.make_chromosome_view(tiny_cool, verbose=False)
    assert result["regions"] == [
        {"chrom": "chr1", "start": 0, "end": 1000, "name": "chr1"},
        {"chrom": "chr2", "start": 0, "end": 600, "name": "chr2"},
    ]


def test_read_viewframe_bed_accepts_valid_and_rejects_duplicate_names(workdir: Path) -> None:
    pytest.importorskip("cooltools")
    good = workdir / "view_good.bed"
    good.write_text("chr1\t0\t1000\tchr1\nchr2\t0\t600\tchr2\n", encoding="utf-8")
    result = TF.read_viewframe_bed(good, verbose=False)
    assert result["region_count"] == 2
    bad = workdir / "view_bad.bed"
    bad.write_text("chr1\t0\t1000\tdup\nchr2\t0\t600\tdup\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="not valid"):
        TF.read_viewframe_bed(bad, verbose=False)


def test_read_expected_table_reads_cis_expected_file(tiny_balanced: Path, workdir: Path) -> None:
    cooler = pytest.importorskip("cooler")
    cooltools = pytest.importorskip("cooltools")
    from cooltools.lib import common as cooltools_common

    clr = cooler.Cooler(str(tiny_balanced))
    expected = cooltools.expected_cis(
        clr, view_df=cooltools_common.make_cooler_view(clr), clr_weight_name="weight", nproc=1
    )
    columns = ["region1", "region2", "dist", "n_valid", "count.avg", "balanced.avg"]
    path = workdir / "expected.tsv"
    expected[columns].to_csv(path, sep="\t", index=False)
    result = TF.read_expected_table(path, limit=3, verbose=False)
    assert result["row_count"] == len(expected)
    assert result["rows"][0]["dist"] == 0
    assert result["columns"] == columns


def test_observed_over_expected_all_ones_is_one(random_balanced: Path) -> None:
    pytest.importorskip("cooltools")
    matrix = [[1.0] * 4 for _ in range(4)]
    result = TF.observed_over_expected_matrix(matrix, verbose=False)
    assert all(
        value == pytest.approx(1.0) for row in result["observed_over_expected"] for value in row
    )
    with pytest.raises(ConfigurationError, match="symmetric"):
        TF.observed_over_expected_matrix([[1.0, 2.0], [3.0, 1.0]], verbose=False)


def test_iterative_correction_all_ones_hand_and_asymmetric_direct() -> None:
    numutils = pytest.importorskip("cooltools.lib.numutils")
    ones = TF.iterative_correction_dense([[1.0] * 4 for _ in range(4)], verbose=False)
    assert ones["converged"] is True and ones["iterations"] == 0  # hand: already balanced
    assert ones["bias"] == [1.0, 1.0, 1.0, 1.0]
    rows = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    asym = TF.iterative_correction_dense(rows, symmetric=False, max_iter=200, verbose=False)
    direct, bias_row, _, _ = numutils.iterative_correction_asymmetric(
        __import__("numpy").asarray(rows), max_iter=200, tol=1e-5
    )
    assert asym["corrected"][0][0] == pytest.approx(float(direct[0, 0]))
    assert asym["bias"]["row"][0] == pytest.approx(float(bias_row[0]))
    with pytest.raises(ConfigurationError, match="symmetric matrices only"):
        TF.iterative_correction_dense(rows, symmetric=False, ignore_diags=1, verbose=False)


def test_dense_cis_eigenvectors_match_direct_decomposition() -> None:
    np = pytest.importorskip("numpy")
    eigdecomp = pytest.importorskip("cooltools.api.eigdecomp")
    rng = np.random.default_rng(11)
    raw = rng.random((8, 8))
    matrix = (raw + raw.T) / 2 + 0.1
    result = TF.dense_cis_eigenvectors(matrix.tolist(), n_eigs=2, ignore_diags=0, verbose=False)
    with TF._reproducible_eigen():
        values, vectors = eigdecomp.cis_eig(matrix, n_eigs=2, ignore_diags=0)
    assert len(result["eigenvalues"]) == 2
    assert len(result["eigenvectors"][0]) == 8
    assert result["eigenvalues"][0] == pytest.approx(float(values[0]))
    assert abs(result["eigenvectors"][1][3]) == pytest.approx(abs(float(vectors[1][3])))
    assert max(result["eigenvectors"][0], key=abs) > 0  # convention: largest entry positive


def test_expected_from_dense_matrix_hand_diagonal_sums() -> None:
    pytest.importorskip("cooltools")
    matrix = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
    result = TF.expected_from_dense_matrix(matrix, ignore_diags=0, verbose=False)
    by_dist = {row["dist"]: row for row in result["diagonals"]}
    # hand: diagonal 0 sums to 1+4+6 = 11 over 3 pixels; diagonal 1: 2+5 = 7 over 2; diagonal 2: 3
    assert by_dist[0]["balanced.sum"] == pytest.approx(11.0) and by_dist[0]["n_valid"] == 3
    assert by_dist[1]["balanced.sum"] == pytest.approx(7.0) and by_dist[1]["n_valid"] == 2
    assert by_dist[2]["balanced.sum"] == pytest.approx(3.0) and by_dist[2]["n_valid"] == 1


def test_digitize_track_matches_direct_and_hand_classes(workdir: Path) -> None:
    saddle = pytest.importorskip("cooltools.api.saddle")
    import pandas as pd

    path = workdir / "digitize.bed"
    path.write_text(
        "chr1\t0\t100\t0.0\nchr1\t100\t200\t1.0\nchr1\t200\t300\t2.0\nchr1\t300\t400\t3.0\n"
    )
    result = TF.digitize_track(path, n_bins=2, qrange=(0.0, 1.0), verbose=False)
    table = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "value"])
    digitized, edges = saddle.digitize(table, 2, qrange=(0.0, 1.0))
    assert [row["class"] for row in result["intervals"]] == list(digitized["value.d"].astype(int))
    assert result["intervals"][0]["class"] == 1 and result["intervals"][1]["class"] == 1  # hand
    assert result["bin_edges"] == pytest.approx(list(edges))


# --- pairtools_functions ----------------------------------------------------------------


def test_find_restriction_fragment_hand_values() -> None:
    pytest.importorskip("pairtools")
    sites = {"chrA": [0, 300, 700, 1000], "chrB": [0, 800]}
    middle = PF.find_restriction_fragment(sites, "chrA", 450, verbose=False)
    assert (middle["fragment_index"], middle["start"], middle["end"]) == (1, 300, 700)
    assert middle["length"] == 400
    first = PF.find_restriction_fragment(sites, "chrA", 0, verbose=False)
    last = PF.find_restriction_fragment(sites, "chrA", 1000, verbose=False)
    assert (first["fragment_index"], last["fragment_index"]) == (0, 2)
    with pytest.raises(ConfigurationError, match="position must be"):
        PF.find_restriction_fragment(sites, "chrA", 1200, verbose=False)


def test_compute_contact_scaling_counts_match_hand_and_direct(workdir: Path) -> None:
    scaling = pytest.importorskip("pairtools.lib.scaling")
    pairs = workdir / "scaling.pairs"
    pairs.write_text(TINY_PAIRS, encoding="utf-8")
    result = PF.compute_contact_scaling(
        pairs, chromsizes=TINY_SIZES, max_distance=1000, bins_per_decade=2, verbose=False
    )
    # hand: six cis pairs (distances 140, 140, 235, 10, 50, 220) and one trans pair
    assert result["cis_pairs"] == 6
    assert sum(row["n_pairs"] for row in result["bins"]) == 6
    assert result["trans_pairs"] == 1
    cis, _ = scaling.compute_scaling(
        str(pairs), chromsizes=TINY_SIZES, dist_range=(1, 1000), n_dist_bins_decade=2, nproc_in=1
    )
    assert int(cis["n_pairs"].sum()) == 6


def test_normalize_scaling_curve_hand_mean_window() -> None:
    pytest.importorskip("pairtools")
    result = PF.normalize_scaling_curve(
        [1, 10, 100], [0.5, 0.25, 0.1], norm_window=(1, 100), verbose=False
    )
    # hand: mean of 0.5, 0.25, 0.1 in the window is 0.85 / 3; each value is divided by it
    assert result["factor"] == pytest.approx(0.85 / 3)
    values = [point["frequency"] for point in result["points"]]
    assert values == pytest.approx([0.5 / (0.85 / 3), 0.25 / (0.85 / 3), 0.1 / (0.85 / 3)])


def test_assign_genomic_regions_by_hand() -> None:
    pytest.importorskip("pairtools")
    positions = [
        {"chrom": "chrA", "pos": 50},
        {"chrom": "chrA", "pos": 1500},
        {"chrom": "chrB", "pos": 10},
        {"chrom": "chrC", "pos": 5},
    ]
    regions = [
        {"chrom": "chrA", "start": 0, "end": 1000},
        {"chrom": "chrA", "start": 1000, "end": 2000},
        {"chrom": "chrB", "start": 0, "end": 800},
    ]
    result = PF.assign_genomic_regions(positions, regions, verbose=False)
    indices = [row["region_index"] for row in result["positions"]]
    assert indices == [0, 1, 0, None]
    assert result["assigned_count"] == 3
    assert result["positions"][1]["region_start"] == 1000


def test_find_read_mismatches_by_hand() -> None:
    pytest.importorskip("pairtools")
    result = PF.find_read_mismatches(
        "ACTT",
        [30, 31, 32, 33],
        [(0, 100, "A"), (1, 101, "C"), (2, 102, "G"), (3, 103, "T")],
        verbose=False,
    )
    # hand: read base 2 is T but the reference base is G; the other three bases agree
    assert result["mismatch_count"] == 1
    assert result["mismatches"][0] == {
        "ref_base": "G",
        "read_base": "T",
        "quality": 32,
        "ref_pos": 102,
        "read_pos": 2,
    }
    with pytest.raises(ConfigurationError, match="one value per base"):
        PF.find_read_mismatches("AC", [30], [(0, 1, "A")], verbose=False)


def test_parse_sam_alignments_matches_hand_coordinates(workdir: Path) -> None:
    pytest.importorskip("pysam")
    pytest.importorskip("pairtools.lib.parse_pysam")
    sam = workdir / "sample.sam"
    sam.write_text(SAM_TEXT, encoding="utf-8")
    result = PF.parse_sam_alignments(sam, min_mapq=30, max_records=2, verbose=False)
    first, second = result["alignments"]
    # hand: a 10M read starting at 1-based 101 spans 101-110 on the + strand
    assert (first["chrom"], first["pos5"], first["pos3"], first["strand"]) == (
        "chrA",
        101,
        110,
        "+",
    )
    assert first["mapq"] == 60 and first["cigar"] == "10M" and first["read_id"] == "r1"
    assert (second["chrom"], second["pos5"], second["pos3"]) == ("chrB", 301, 310)
    no_header = workdir / "no_header.sam"
    no_header.write_text("r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="readable SAM or BAM"):
        PF.parse_sam_alignments(no_header, verbose=False)


def test_read_pairs_table_counts_rows_and_header(workdir: Path) -> None:
    pairsio = pytest.importorskip("pairtools.lib.pairsio")
    pairs = workdir / "table.pairs"
    pairs.write_text(TINY_PAIRS, encoding="utf-8")
    result = PF.read_pairs_table(pairs, limit=2, verbose=False)
    assert result["row_count"] == 7  # hand: seven pair lines
    assert result["chromsizes"] == TINY_SIZES
    assert len(result["rows"]) == 2
    direct_frame, _, _ = pairsio.read_pairs(str(pairs), nproc=1)
    assert len(direct_frame) == result["row_count"]


def test_make_pairs_header_matches_direct_header() -> None:
    headerops = pytest.importorskip("pairtools.lib.headerops")
    result = PF.make_pairs_header(chromsizes=TINY_SIZES, assembly="hg38", verbose=False)
    direct = headerops.make_standard_pairsheader(assembly="hg38", chromsizes=dict(TINY_SIZES))
    assert result["lines"] == direct
    assert "#chromsize: chr1 1000" in result["lines"]
    assert "#genome_assembly: hg38" in result["lines"]


def test_subset_pairs_header_keeps_one_chromosome() -> None:
    pytest.importorskip("pairtools")
    header = PF.make_pairs_header(chromsizes=TINY_SIZES, verbose=False)["lines"]
    result = PF.subset_pairs_header(header, ["chr2"], verbose=False)
    chromsize_lines = [line for line in result["lines"] if line.startswith("#chromsize")]
    assert chromsize_lines == ["#chromsize: chr2 600"]  # hand: chr1 line is dropped
    assert result["output_line_count"] == len(header) - 1


def test_read_chrom_order_numbers_file_order_then_extra(workdir: Path) -> None:
    pytest.importorskip("pairtools")
    path = workdir / "order.tsv"
    path.write_text("chr1\t1000\nchr2\t600\n", encoding="utf-8")
    result = PF.read_chrom_order(path, sam_chroms=["chr2", "chr1", "chrC"], verbose=False)
    assert result["chromosomes"] == [
        {"chrom": "chr1", "order": 1},
        {"chrom": "chr2", "order": 2},
        {"chrom": "chrC", "order": 3},
    ]


def test_assign_haplotype_phase_xb_and_xa_by_hand() -> None:
    phase_module = pytest.importorskip("pairtools.lib.phase")
    xb = PF.assign_haplotype_phase(
        "chr1_1", tag_value="", primary_score=50, secondary_score=40, variant="XB", verbose=False
    )
    assert (xb["phase"], xb["chromosome_base"], xb["scores"]) == ("0", "chr1", [50, 40, -1])
    direct = phase_module.phase_side_XB("chr1_1", "", 50, 40, ("_1", "_2"))
    assert xb["phase"] == direct[0]
    xa = PF.assign_haplotype_phase(
        "chr1_2",
        tag_value="chr1_1,100,50M,2;",
        primary_score=50,
        secondary_score=50,
        variant="XA",
        mismatches=1,
        verbose=False,
    )
    # hand: the alternative hit lies on the other homologue, so the side is homologous
    assert xa["phase"] == "." and xa["chromosome_base"] == "chr1"
    with pytest.raises(ConfigurationError, match="mismatches"):
        PF.assign_haplotype_phase(
            "chr1_2",
            tag_value="",
            primary_score=50,
            secondary_score=40,
            variant="XA",
            verbose=False,
        )
