"""Offline parity tests for P38; optional packages are checked per test family."""

import pytest


@pytest.mark.parametrize("region", ["chr1", "chr1:10-20", "chr1:0-1"])
def test_cooler_region_and_bins(tmp_path, region):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    from dnakit.epigenome import cooler_cooltools_p38_functions as m

    path = tmp_path / "tiny.cool"
    bins = cooler.binnify(pd.Series({"chr1": 40}), 10)
    pixels = pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [0, 2], "count": [2, 1]})
    cooler.create_cooler(str(path), bins, pixels)
    direct = cooler.Cooler(str(path))
    got = m.cooler_region_bin_ids(str(path), region, verbose=False)
    assert got["offset"] == direct.offset(region)
    assert got["extent"] == list(direct.extent(region))
    table = m.query_cooler_bins(str(path), region=region, verbose=False)
    assert table["bin_count"] == len(direct.bins().fetch(region))
    assert m.detect_cooler_file_type(str(path), verbose=False)["file_type"] == "cooler"


@pytest.mark.parametrize("binsize", [5000, 10000, 25000])
def test_cooltools_kernels_match_upstream(binsize):
    np = pytest.importorskip("numpy")
    pytest.importorskip("cooltools")
    from cooltools.api.dotfinder import recommend_kernels

    from dnakit.epigenome import cooler_cooltools_p38_functions as m

    got = m.build_dot_calling_kernels(binsize=binsize, verbose=False)
    for key, value in recommend_kernels(binsize).items():
        np.testing.assert_array_equal(got["kernels"][key], value)
    assert m.validate_dot_calling_kernels(got["kernels"], binsize, verbose=False)["valid"]


@pytest.mark.parametrize("kappa", [1.0, 2.0, 0.0])
def test_dendropy_hky85_matches_upstream(kappa):
    np = pytest.importorskip("numpy")
    pytest.importorskip("dendropy")
    from dendropy.model.discrete import Hky85

    from dnakit.phylogeny import dendropy_p38_functions as m

    expected = Hky85(kappa=kappa, base_freqs=[0.25] * 4)
    got = m.calculate_hky85_rate_matrix(kappa=kappa, verbose=False)
    np.testing.assert_allclose(got["qmatrix"], expected.qmatrix(1.0))
    assert got["corrected_substitution_rate"] == expected.corrected_substitution_rate(1.0)


@pytest.mark.parametrize(
    "sequences",
    [
        {"a": "ACGT", "b": "ACGA", "c": "ACGG"},
        {"a": "A", "b": "A"},
        {"a": "AA", "b": "AT", "c": "TT"},
    ],
)
def test_dendropy_folded_spectrum_matches_upstream(sequences):
    dendropy = pytest.importorskip("dendropy")
    from dnakit.phylogeny import dendropy_p38_functions as m

    expected = dendropy.DnaCharacterMatrix.from_dict(sequences).folded_site_frequency_spectrum(
        False
    )
    assert (
        m.calculate_folded_site_frequency_spectrum(sequences, verbose=False)["spectrum"] == expected
    )


@pytest.mark.parametrize("sequence", ["ACGT", "A", "GGGGAAAA"])
def test_dnachisel_rule_matches_upstream(sequence):
    dc = pytest.importorskip("dnachisel")
    from dnakit.molecular_lab import dnachisel_p38_functions as m

    direct = dc.DnaOptimizationProblem(
        sequence=sequence,
        objectives=[dc.SequenceLengthBounds(min_length=2, max_length=6)],
        logger=None,
    )
    expected = next(iter(direct.objectives_evaluations()))
    got = m.evaluate_dnachisel_design_rule(
        sequence,
        "SequenceLengthBounds",
        parameters={"min_length": 2, "max_length": 6},
        verbose=False,
    )
    assert got["evaluation"]["score"] == expected.score


@pytest.mark.parametrize("start,end,strand", [(0, 4, "+"), (2, 6, "-"), (0, 1, "+")])
def test_selene_sequence_and_encoding_match_upstream(tmp_path, start, end, strand):
    np = pytest.importorskip("numpy")
    pytest.importorskip("selene_sdk")
    from selene_sdk.sequences import Genome

    from dnakit.predictions import selene_sdk_p38_functions as m

    path = tmp_path / "tiny.fa"
    path.write_text(">chr1\nACGTACGT\n")
    direct = Genome(str(path))
    seq = m.get_selene_genome_window(str(path), "chr1", start, end, strand=strand, verbose=False)
    assert seq["sequence"] == direct.get_sequence_from_coords(
        "chr1", start, end, strand=strand, pad=False
    )
    got = m.get_selene_genome_window(
        str(path), "chr1", start, end, strand=strand, output="encoding", verbose=False
    )
    np.testing.assert_array_equal(
        got["encoding"],
        direct.get_encoding_from_coords("chr1", start, end, strand=strand, pad=False),
    )


@pytest.mark.parametrize("n_cells", [1, 2, 3])
def test_snapatac_import_matches_upstream(tmp_path, n_cells):
    snap = pytest.importorskip("snapatac2")
    from dnakit.regulation import snapatac2_p38_functions as m

    path = tmp_path / "fragments.tsv"
    path.write_text("".join(f"chr1\t{10 * i}\t{10 * i + 5}\tcell{i}\t1\n" for i in range(n_cells)))
    direct = snap.pp.import_fragments(
        path, {"chr1": 1000}, min_num_fragments=1, chrM=["chrM", "M"], n_jobs=1, tempdir=tmp_path
    )
    got = m.import_atac_fragments(
        str(path), {"chr1": 1000}, str(tmp_path / "out.h5ad"), min_num_fragments=1, verbose=False
    )
    assert got["cell_count"] == direct.n_obs
    assert sorted(row["n_fragment"] for row in got["qc"]) == sorted(
        direct.obs["n_fragment"].to_list()
    )
