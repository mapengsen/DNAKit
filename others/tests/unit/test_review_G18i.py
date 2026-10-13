"""Tests for the G18i review wrappers in dnakit.regulation (gimmemotifs, tangermeme, snapATAC2)."""

from __future__ import annotations

import importlib
import json
import math
import random
import tempfile
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import numpy
import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.regulation import gimmemotifs_review_functions as gm
from dnakit.regulation import snapatac2_review_functions as sa
from dnakit.regulation import tangermeme_review_functions as tm

# Contigs must be much longer than gimmemotifs' 128-base background windows (about 4 kb here).
_GM_GENOME = {
    "chr1": "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70,
    "chr2": "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60,
}
_GM_PFM = ">GATA\n0 0 10 0\n0 0 0 10\n10 0 0 0\n0 10 0 0\n"


def _gm_write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _gm_fasta(path: Path, records: dict[str, str]) -> Path:
    lines = []
    for name, text in records.items():
        lines.append(f">{name}")
        lines.extend(text[start : start + 60] for start in range(0, len(text), 60))
    return _gm_write(path, "\n".join(lines) + "\n")


def test_gm_fasta_record_count_matches_upstream_and_hand_count(tmp_path: Path) -> None:
    utils = pytest.importorskip("gimmemotifs.utils")
    fasta = _gm_fasta(tmp_path / "three.fa", {"a": "ACGT", "b": "GGCC", "c": "TTAA"})
    result = gm.fasta_record_count(fasta)
    assert result["record_count"] == 3 == utils.number_of_seqs_in_file(str(fasta))
    assert result["operation"] == "fasta_record_count"
    assert result["provenance"]["package"] == "gimmemotifs"


def test_gm_equal_length_motif_pair_pads_with_background_and_matches_upstream() -> None:
    utils = pytest.importorskip("gimmemotifs.utils")
    np = pytest.importorskip("numpy")
    pfm_a = [[10, 0, 0, 0], [0, 10, 0, 0]]
    pfm_b = [[0, 0, 10, 0]]
    result = gm.equal_length_motif_pair(pfm_a, pfm_b, pos=0)
    assert result["length"] == 2
    assert result["pfm_a"] == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]
    # Hand value: the missing second row of pfm_b is the uniform background.
    assert result["pfm_b"] == [[0.0, 0.0, 1.0, 0.0], [0.25, 0.25, 0.25, 0.25]]
    direct_a, direct_b = utils.make_equal_length(
        np.asarray(pfm_a, dtype=float) / 10.0, np.asarray([[0.0, 0.0, 1.0, 0.0]]), 0
    )
    np.testing.assert_allclose(result["pfm_a"], direct_a)
    np.testing.assert_allclose(result["pfm_b"], direct_b)


def test_gm_equal_length_motif_pair_rejects_bad_truncate() -> None:
    with pytest.raises(ConfigurationError):
        gm.equal_length_motif_pair([[1, 0, 0, 0]], [[0, 1, 0, 0]], truncate="middle")


def test_gm_combine_summit_files_keeps_strongest_summit_and_matches_upstream(
    tmp_path: Path,
) -> None:
    preprocessing = pytest.importorskip("gimmemotifs.preprocessing")
    genome = _gm_fasta(tmp_path / "genome.fa", _GM_GENOME)
    first = _gm_write(
        tmp_path / "summits1.bed", "chr1\t100\t101\tp1\t9.0\nchr2\t500\t501\tp2\t7.0\n"
    )
    second = _gm_write(
        tmp_path / "summits2.bed", "chr1\t110\t111\tp3\t8.0\nchr1\t2000\t2001\tp4\t6.0\n"
    )
    result = gm.combine_summit_files([first, second], genome, window=200)
    # Hand value: chr1:100 (value 9) and chr1:110 (value 8) fall in one 200-base window.
    assert [(s["chrom"], s["start"], s["end"]) for s in result["summits"]] == [
        ("chr1", 100, 101),
        ("chr1", 2000, 2001),
        ("chr2", 500, 501),
    ]
    direct = preprocessing.combine_peaks([str(first), str(second)], str(genome), window=200)
    assert [int(v) for v in direct["start"]] == [s["start"] for s in result["summits"]]


def test_gm_peak_coverage_table_matches_hand_mean_and_upstream(tmp_path: Path) -> None:
    pyBigWig = pytest.importorskip("pyBigWig")
    preprocessing = pytest.importorskip("gimmemotifs.preprocessing")
    peaks = _gm_write(tmp_path / "peaks.bed", "chr1\t100\t101\tp1\t9.0\nchr2\t500\t501\tp2\t7.0\n")
    track = str(tmp_path / "cov.bw")
    bw = pyBigWig.open(track, "w")
    bw.addHeader([("chr1", 4000), ("chr2", 3000)])
    bw.addEntries(["chr1"], [0], ends=[4000], values=[2.0])
    bw.close()
    result = gm.peak_coverage_table(peaks, [track], window=200)
    assert result["tracks"] == ["cov"]
    # Hand value: coverage 2.0 with log(1 + x) gives log(3) at the chr1 peak; chr2 has no data.
    assert result["values"][0][0] == pytest.approx(math.log(3.0))
    assert result["values"][1][0] == pytest.approx(0.0)
    direct = preprocessing.coverage_table(str(peaks), [track], window=200, ncpus=1)
    assert result["values"][0][0] == pytest.approx(float(direct.iloc[0, 0]))


def test_gm_scan_motifs_to_file_writes_hits_and_refuses_overwrite(tmp_path: Path) -> None:
    scanner = pytest.importorskip("gimmemotifs.scanner")
    genome = _gm_fasta(tmp_path / "genome.fa", _GM_GENOME)
    pfm = _gm_write(tmp_path / "motifs.pfm", _GM_PFM)
    output = tmp_path / "hits.bed"
    result = gm.scan_motifs_to_file(genome, pfm, output, cutoff=0.0)
    assert result["hit_count"] == len(result["hits"]) > 0
    assert all(hit["motif"] == "GATA" for hit in result["hits"])
    assert all(hit["strand"] in ("+", "-") for hit in result["hits"])
    starts = [(hit["chrom"], hit["start"]) for hit in result["hits"]]
    assert starts == sorted(starts)
    direct_file = tmp_path / "direct.bed"
    scanner.scan_to_file(
        str(genome),
        pfmfile=str(pfm),
        filepath_or_buffer=str(direct_file),
        cutoff=0.0,
        bed=True,
        ncpus=1,
        zscore=False,
        gcnorm=False,
        progress=False,
    )
    direct_rows = [
        line.split("\t")
        for line in direct_file.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(direct_rows) == result["hit_count"]
    with pytest.raises(ConfigurationError):
        gm.scan_motifs_to_file(genome, pfm, output, cutoff=0.0)


def test_gm_scan_motifs_to_file_rejects_bad_cutoff(tmp_path: Path) -> None:
    genome = _gm_fasta(tmp_path / "genome.fa", _GM_GENOME)
    pfm = _gm_write(tmp_path / "motifs.pfm", _GM_PFM)
    with pytest.raises(ConfigurationError):
        gm.scan_motifs_to_file(genome, pfm, tmp_path / "x.bed", cutoff=1.5)


def test_gm_scan_region_scores_keeps_every_region_and_matches_upstream(tmp_path: Path) -> None:
    scanner = pytest.importorskip("gimmemotifs.scanner")
    np = pytest.importorskip("numpy")
    genome = _gm_fasta(tmp_path / "genome.fa", _GM_GENOME)
    pfm = _gm_write(tmp_path / "motifs.pfm", _GM_PFM)
    regions = _gm_write(tmp_path / "regions.txt", "chr1:100-300\nchr2:400-600\n")
    result = gm.scan_region_scores(regions, genome, pfm)
    assert result["regions"] == ["chr1:100-300", "chr2:400-600"]
    assert result["motifs"] == ["GATA"]
    assert result["region_count"] == 2
    table = _gm_write(tmp_path / "table.tsv", "region\nchr1:100-300\nchr2:400-600\n")
    frame = scanner.scan_regionfile_to_table(
        str(table),
        str(genome),
        "score",
        pfmfile=str(pfm),
        ncpus=1,
        zscore=False,
        gc=False,
        progress=False,
    )
    direct = [float(frame.loc[region, "GATA"]) for region in result["regions"]]
    np.testing.assert_allclose(result["scores"], [[value] for value in direct])


def test_gm_scan_region_scores_rejects_malformed_region(tmp_path: Path) -> None:
    genome = _gm_fasta(tmp_path / "genome.fa", _GM_GENOME)
    pfm = _gm_write(tmp_path / "motifs.pfm", _GM_PFM)
    regions = _gm_write(tmp_path / "regions.txt", "chr1 100 300\n")
    with pytest.raises(ConfigurationError):
        gm.scan_region_scores(regions, genome, pfm)


def test_gm_motif_activity_scores_matches_direct_moap(tmp_path: Path) -> None:
    moap = pytest.importorskip("gimmemotifs.maelstrom.moap")
    labels = _gm_write(
        tmp_path / "labels.txt",
        "label\tgroup\nr0\t1\nr1\t1\nr2\t0\nr3\t0\nr4\t1\nr5\t0\n",
    )
    counts = _gm_write(
        tmp_path / "counts.txt",
        "region\tGATA\tTTGA\nr0\t3\t0\nr1\t2\t1\nr2\t0\t0\nr3\t0\t1\nr4\t4\t0\nr5\t0\t2\n",
    )
    result = gm.motif_activity_scores(labels, counts)
    assert result["classes"] == ["0", "1"]
    assert result["motifs"] == ["GATA", "TTGA"]
    direct = moap.moap(
        str(labels),
        method="hypergeom",
        scoring="count",
        outfile=None,
        motiffile=str(counts),
        ncpus=1,
        zscore=False,
        gc=False,
        progress=False,
    )
    assert result["values"][0][1] == pytest.approx(float(direct.loc["GATA", 1]))
    assert result["values"][1][0] == pytest.approx(float(direct.loc["TTGA", 0]))


def test_gm_motif_score_metric_roc_auc_hand_value_and_upstream() -> None:
    rocmetrics = pytest.importorskip("gimmemotifs.rocmetrics")
    np = pytest.importorskip("numpy")
    foreground = [3, 4, 5, 6]
    background = [1, 2, 3, 4]
    result = gm.motif_score_metric(foreground, background, metric="roc_auc")
    # Hand value: 14 of the 16 foreground-background pairs count, with the two ties at half weight.
    assert result["value"] == pytest.approx(0.875)
    direct = rocmetrics.roc_auc(np.asarray(foreground, float), np.asarray(background, float))
    assert result["value"] == pytest.approx(float(direct))


def test_gm_motif_score_metric_curve_and_family_members() -> None:
    rocmetrics = pytest.importorskip("gimmemotifs.rocmetrics")
    np = pytest.importorskip("numpy")
    foreground = [3.0, 4.0, 5.0, 6.0]
    background = [1.0, 2.0, 3.0, 4.0]
    curve = gm.motif_score_metric(foreground, background, metric="roc_values")
    x_values, y_values = rocmetrics.roc_values(np.asarray(foreground), np.asarray(background))
    assert curve["curve"]["fpr"] == pytest.approx([float(v) for v in x_values])
    assert curve["curve"]["tpr"] == pytest.approx([float(v) for v in y_values])
    enrichment = gm.motif_score_metric(foreground, background, metric="enr_at_fpr", fpr=0.25)
    direct = rocmetrics.enr_at_fpr(np.asarray(foreground), np.asarray(background), fpr=0.25)
    assert enrichment["value"] == pytest.approx(float(direct))


def test_gm_motif_score_metric_rejects_bad_arguments() -> None:
    with pytest.raises(ConfigurationError):
        gm.motif_score_metric([1.0], [2.0], metric="not_a_metric")
    with pytest.raises(ConfigurationError):
        gm.motif_score_metric([1.0, math.nan], [2.0])
    with pytest.raises(ConfigurationError):
        gm.motif_score_metric([1.0], [2.0], fpr=0.0)


def test_gm_motif_position_ks_matches_upstream_and_significance_is_minus_log10() -> None:
    rocmetrics = pytest.importorskip("gimmemotifs.rocmetrics")
    np = pytest.importorskip("numpy")
    foreground = [10, 20, 30, 40, 50]
    background = [100, 110, 120, 130, 140]
    pvalue = gm.motif_position_ks(foreground, background, output="pvalue")
    direct = rocmetrics.ks_pvalue(np.asarray(foreground, float), np.asarray(background, float))
    assert pvalue["value"] == pytest.approx(float(direct))
    significance = gm.motif_position_ks(foreground, background, output="significance")
    assert significance["value"] == pytest.approx(-math.log10(pvalue["value"]))


def test_gm_gimmemotifs_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = importlib.import_module

    def fake_import(name: str, package: str | None = None):  # type: ignore[no-untyped-def]
        if name.startswith("gimmemotifs"):
            raise ImportError("simulated missing gimmemotifs")
        return real_import(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", fake_import)
    with pytest.raises(BackendUnavailableError, match="bio-gimmemotifs"):
        gm.motif_score_metric([1.0], [2.0])


_TM_SEQUENCE = "ACGTACGTAA"


def test_tm_shuffle_region_keeps_outside_bases_and_region_composition() -> None:
    result = tm.shuffle_region(_TM_SEQUENCE, 2, 6, n=3, random_state=0)
    assert result["n"] == 3 and len(result["sequences"]) == 3
    for sequence in result["sequences"]:
        # Hand values: positions 0-1 ("AC") and 6-9 ("GTAA") are outside the region.
        assert sequence[:2] == "AC"
        assert sequence[6:] == "GTAA"
        assert Counter(sequence[2:6]) == Counter("GTAC")


def test_tm_shuffle_region_matches_direct_upstream_with_same_seed() -> None:
    ersatz = pytest.importorskip("tangermeme.ersatz")
    utils = pytest.importorskip("tangermeme.utils")
    x = utils.one_hot_encode(_TM_SEQUENCE)[None].float()
    direct = ersatz.shuffle(x, 2, 6, n=2, random_state=5)
    expected = ["".join("ACGT"[i] for i in direct[0, k].argmax(dim=0).tolist()) for k in range(2)]
    result = tm.shuffle_region(_TM_SEQUENCE, 2, 6, n=2, random_state=5)
    assert result["sequences"] == expected


def test_tm_randomize_region_draws_from_probabilities() -> None:
    # Hand value: probs (0, 0, 1, 0) puts G in every drawn position (2 to 5).
    result = tm.randomize_region(_TM_SEQUENCE, 2, 6, probs=[0, 0, 1, 0], n=1, random_state=0)
    assert result["sequences"] == ["ACGGGGGTAA"]
    assert result["probs"] == [0.0, 0.0, 1.0, 0.0]


def test_tm_delete_region_removes_bases_by_hand() -> None:
    # Hand value: removing positions 2 to 5 ("GTAC") leaves "AC" + "GTAA".
    result = tm.delete_region(_TM_SEQUENCE, 2, 6)
    assert result["sequence"] == "ACGTAA"
    assert result["deleted"] == 4 and result["length"] == 6


def test_tm_random_one_hot_sequences_hand_and_upstream() -> None:
    utils = pytest.importorskip("tangermeme.utils")
    fixed = tm.random_one_hot_sequences(2, 5, probs=[1, 0, 0, 0], random_state=0)
    assert fixed["sequences"] == ["AAAAA", "AAAAA"]
    result = tm.random_one_hot_sequences(2, 6, random_state=7)
    direct = utils.random_one_hot((2, 4, 6), probs=[[0.25, 0.25, 0.25, 0.25]], random_state=7)
    expected = ["".join("ACGT"[i] for i in direct[k].argmax(dim=0).tolist()) for k in range(2)]
    assert result["sequences"] == expected


def test_tm_place_motif_set_places_motifs_at_spacing() -> None:
    # Hand value: GG at 1-2, then TT starts 2 bases after GG ends (position 5).
    result = tm.place_motif_set("AAAAAAAAAA", ["GG", "TT"], [2], start=1)
    assert result["sequence"] == "AGGAATTAAA"
    assert result["spacing"] == [2]


def test_tm_place_motif_set_rejects_wrong_spacing_count() -> None:
    with pytest.raises(ConfigurationError):
        tm.place_motif_set("AAAAAAAAAA", ["GG", "TT"], [], start=1)


def test_tm_local_dinucleotide_shuffle_keeps_composition_and_matches_upstream() -> None:
    ersatz = pytest.importorskip("tangermeme.ersatz")
    utils = pytest.importorskip("tangermeme.utils")
    torch = pytest.importorskip("torch")
    text = "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 25
    result = tm.local_dinucleotide_shuffle(
        [text, text[::-1]], n=2, bin_size=512, min_bin_size=256, random_state=0
    )
    assert len(result["sequences"]) == 2 and all(len(row) == 2 for row in result["sequences"])
    for row, source in zip(result["sequences"], [text, text[::-1]], strict=True):
        for shuffled in row:
            assert Counter(shuffled) == Counter(source)
    x = torch.stack([utils.one_hot_encode(text).float(), utils.one_hot_encode(text[::-1]).float()])
    direct = ersatz.local_dinucleotide_shuffle(
        x, n=2, bin_size=512, min_bin_size=256, random_state=0, verbose=False
    )
    expected = [
        ["".join("ACGT"[i] for i in direct[r, k].argmax(dim=0).tolist()) for k in range(2)]
        for r in range(2)
    ]
    assert result["sequences"] == expected


def test_tm_pwm_consensus_sequence_takes_largest_column() -> None:
    # Hand value: the largest column of each row is A, G then T.
    pwm = [[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.1, 0.1, 0.1, 0.7]]
    assert tm.pwm_consensus_sequence(pwm)["consensus"] == "AGT"


def test_tm_pwm_consensus_rejects_all_zero_row() -> None:
    with pytest.raises(ConfigurationError):
        tm.pwm_consensus_sequence([[0, 0, 0, 0]])


def test_tm_example_to_sequence_coords_adds_locus_start() -> None:
    # Hand value: example positions 5 to 10 in a locus starting at 100 are chr1:105-110.
    result = tm.example_to_sequence_coords(
        [{"example_idx": 0, "start": 5, "end": 10}],
        [{"chrom": "chr1", "start": 100, "end": 200}],
    )
    assert result["intervals"] == [{"chrom": "chr1", "start": 105, "end": 110}]


def test_tm_one_hot_to_fasta_file_writes_records_and_refuses_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    result = tm.one_hot_to_fasta_file(["ACGT", "TTGA"], output, headers=["a", "b"])
    assert result["record_count"] == 2 and result["headers"] == ["a", "b"]
    lines = [
        line.strip() for line in output.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    records = {}
    current = None
    for line in lines:
        if line.startswith(">"):
            current = line[1:].strip()
            records[current] = ""
        else:
            records[current] += line
    assert records == {"a": "ACGT", "b": "TTGA"}
    with pytest.raises(ConfigurationError):
        tm.one_hot_to_fasta_file(["ACGT"], output)


def test_tm_gapped_kmer_features_matches_direct_upstream() -> None:
    kmers = pytest.importorskip("tangermeme.kmers")
    utils = pytest.importorskip("tangermeme.utils")
    torch = pytest.importorskip("torch")
    sequences = ["ACGTTGATAAGCCGTAGGAT", "TTGATAAGGGCCATTAGAAG"]
    result = tm.gapped_kmer_features(
        sequences, min_k=2, max_k=2, max_gap=1, max_len=3, max_gkmers=5
    )
    x = torch.stack([utils.one_hot_encode(seq).float() for seq in sequences])
    direct = kmers.gapped_kmers(x, min_k=2, max_k=2, max_gap=1, max_len=3, max_gkmers=5).tocoo()
    expected = sorted(
        [int(r), int(c), int(v)]
        for r, c, v in zip(direct.row, direct.col, direct.data, strict=True)
    )
    assert result["counts"] == expected
    assert result["feature_count"] == direct.shape[1]


def test_tm_tangermeme_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = importlib.import_module

    def fake_import(name: str, package: str | None = None):  # type: ignore[no-untyped-def]
        if name.startswith("tangermeme"):
            raise ImportError("simulated missing tangermeme")
        return real_import(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", fake_import)
    with pytest.raises(BackendUnavailableError) as info:
        tm.delete_region(_TM_SEQUENCE, 2, 6)
    assert "bio-tangermeme" in str(info.value) + str(getattr(info.value, "hint", ""))


_TM2_SEQUENCES = ["ACGTACGTAA", "TTCCGGAATT"]
_TM2_GATA_CONTEXT = "TTCCGGAATTGATACCGGTTAACCGG"


def _tm2_torch():
    return pytest.importorskip("torch")


def _tm2_count_model(base: int):
    """Model that returns how many times one base (0=A, 1=C, 2=G, 3=T) occurs in a sequence."""

    torch = _tm2_torch()

    class CountBase(torch.nn.Module):
        def forward(self, x):
            return x[:, base, :].sum(dim=-1, keepdim=True)

    return CountBase().eval()


def _tm2_linear_model(length: int):
    """Linear model with fixed weights (no random draws, so no global seed is touched)."""

    torch = _tm2_torch()
    layer = torch.nn.Linear(4 * length, 2)
    with torch.no_grad():
        layer.weight.copy_(torch.linspace(-1.0, 1.0, 2 * 4 * length).reshape(2, 4 * length))
        layer.bias.zero_()

    class Linear(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.flatten = torch.nn.Flatten()
            self.layer = layer

        def forward(self, x):
            return self.layer(self.flatten(x))

    return Linear().eval()


def _tm2_random_sequences(count: int, length: int, seed: int) -> list[str]:
    rng = numpy.random.RandomState(seed)
    return ["".join(rng.choice(list("ACGT"), size=length)) for _ in range(count)]


def _tm2_one_hot(sequences: list[str]):
    torch = _tm2_torch()
    utils = pytest.importorskip("tangermeme.utils")
    return torch.stack([utils.one_hot_encode(text).float() for text in sequences])


def _tm2_flat(value):
    if isinstance(value, list):
        flat = []
        for item in value:
            flat.extend(_tm2_flat(item))
        return flat
    return [value]


def _tm2_assert_close(actual, expected) -> None:
    left = _tm2_flat(actual)
    right = _tm2_flat(expected)
    assert len(left) == len(right)
    assert left == pytest.approx(right, rel=1e-5, abs=1e-6)


def _tm2_motif_db() -> dict[str, list[list[float]]]:
    """GATA plus 30 random PWMs (one row per position) from a fixed seed; 25 or more targets."""

    gata = [[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]
    rng = numpy.random.RandomState(7)
    database = {"GATA": gata}
    for index in range(30):
        database[f"rand{index}"] = rng.dirichlet(numpy.ones(4), size=6).tolist()
    return database


def test_tm2_ablate_annotations_composition_model_is_unchanged_by_shuffle() -> None:
    # Hand values: the model counts C. Each sequence has two Cs, and shuffling a span keeps its
    # bases, so the count is 2 before and after every ablation.
    annotations = [
        {"example_idx": 0, "start": 0, "end": 4},
        {"example_idx": 1, "start": 2, "end": 6},
    ]
    result = tm.ablate_annotations(
        _tm2_count_model(1), _TM2_SEQUENCES, annotations, n=3, random_state=0
    )
    assert result["operation"] == "ablate_annotations"
    assert result["before"] == [[2.0], [2.0]]
    assert result["after"] == [[[2.0], [2.0], [2.0]], [[2.0], [2.0], [2.0]]]
    json.dumps(result)


def test_tm2_ablate_annotations_matches_direct_upstream() -> None:
    ablate = pytest.importorskip("tangermeme.ablate")
    torch = _tm2_torch()
    sequences = _tm2_random_sequences(2, 30, seed=3)
    model = _tm2_linear_model(30)
    annotations = [
        {"example_idx": 0, "start": 0, "end": 5},
        {"example_idx": 1, "start": 10, "end": 20},
    ]
    result = tm.ablate_annotations(model, sequences, annotations, n=2, random_state=5)
    direct = ablate.ablate_annotations(
        model,
        _tm2_one_hot(sequences),
        torch.tensor([[0, 0, 5], [1, 10, 20]]),
        n=2,
        random_state=5,
    )
    _tm2_assert_close(result["before"], direct.y_befores[:, 0].tolist())
    _tm2_assert_close(result["after"], direct.y_afters[:, 0].tolist())


def test_tm2_ablate_annotations_rejects_span_past_end_and_training_mode() -> None:
    with pytest.raises(ConfigurationError):
        tm.ablate_annotations(
            _tm2_count_model(1),
            _TM2_SEQUENCES,
            [{"example_idx": 0, "start": 2, "end": 11}],
        )
    training = _tm2_count_model(1)
    training.train()
    with pytest.raises(ConfigurationError):
        tm.ablate_annotations(training, _TM2_SEQUENCES, [{"example_idx": 0, "start": 0, "end": 4}])


def test_tm2_annotate_seqlets_names_gata_and_matches_direct_upstream() -> None:
    annotate = pytest.importorskip("tangermeme.annotate")
    pandas = pytest.importorskip("pandas")
    torch = _tm2_torch()
    database = _tm2_motif_db()
    seqlets = [
        {"example_idx": 0, "start": 10, "end": 14},
        {"example_idx": 0, "start": 0, "end": 6},
    ]
    result = tm.annotate_seqlets([_TM2_GATA_CONTEXT], seqlets, database, n_nearest=2)
    # Hand value: positions 10 to 13 read GATA, the first motif in the database, so it is the
    # best hit for the first seqlet with a very small p-value.
    assert result["matches"][0][0]["motif"] == "GATA"
    assert result["matches"][0][0]["p_value"] < 1e-6
    frame = pandas.DataFrame([[0, 10, 14], [0, 0, 6]], columns=["example_idx", "start", "end"])
    pwms = {name: torch.tensor(rows, dtype=torch.float32).T for name, rows in database.items()}
    idxs, pvals = annotate.annotate_seqlets(
        _tm2_one_hot([_TM2_GATA_CONTEXT]), frame, pwms, n_nearest=2, n_jobs=1
    )
    names = list(database)
    for row, hits in enumerate(result["matches"]):
        assert [hit["motif"] for hit in hits] == [names[index] for index in idxs[row].tolist()]
        _tm2_assert_close([hit["p_value"] for hit in hits], pvals[row].tolist())


def test_tm2_annotate_seqlets_rejects_all_zero_pwm_row() -> None:
    with pytest.raises(ConfigurationError):
        tm.annotate_seqlets(
            [_TM2_GATA_CONTEXT],
            [{"example_idx": 0, "start": 0, "end": 4}],
            {"bad": [[0, 0, 0, 0]]},
        )


def test_tm2_count_seqlet_annotations_hand_counts_and_direct_upstream() -> None:
    annotate = pytest.importorskip("tangermeme.annotate")
    torch = _tm2_torch()
    pairs = [[0, 0], [0, 1], [0, 1], [1, 0]]
    # Hand values: example 0 has annotation 0 once and annotation 1 twice; example 1 has
    # annotation 0 once.
    assert tm.count_seqlet_annotations(pairs, verbose=False)["counts"] == [[1, 2], [1, 0]]
    assert tm.count_seqlet_annotations(pairs, dim=0)["counts"] == [2, 2]
    assert tm.count_seqlet_annotations(pairs, dim=1)["counts"] == [3, 1]
    direct = annotate.count_annotations(torch.tensor(pairs), dtype=torch.int64)
    assert tm.count_seqlet_annotations(pairs)["counts"] == direct.tolist()
    with pytest.raises(ConfigurationError):
        tm.count_seqlet_annotations(pairs, dim=2)


def test_tm2_pairwise_annotation_counts_pairs_mode_counts_unique_pairs() -> None:
    annotate = pytest.importorskip("tangermeme.annotate")
    torch = _tm2_torch()
    # Hand values: example 0 holds annotations {0, 1} once each,
    # so the pair 0-1 is counted once and mirrored.
    result = tm.pairwise_annotation_counts([[0, 0], [0, 1], [0, 1], [1, 0]])
    assert result["mode"] == "pairs" and result["max_distance"] is None
    assert result["counts"] == [[0, 1], [1, 0]]
    # Hand values: example 0 holds {0, 1, 2}, giving the three unique pairs 0-1, 0-2 and 1-2.
    triple = tm.pairwise_annotation_counts([[0, 0], [0, 1], [0, 2], [0, 2]])
    assert triple["counts"] == [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
    direct = annotate.pairwise_annotations(torch.tensor([[0, 0], [0, 1], [0, 2], [0, 2]]))
    assert triple["counts"] == direct.tolist()


def test_tm2_pairwise_annotation_counts_spacing_mode_hand_gap_and_direct_upstream() -> None:
    annotate = pytest.importorskip("tangermeme.annotate")
    torch = _tm2_torch()
    rows = [
        {"example_idx": 0, "annotation": 0, "start": 0, "end": 4},
        {"example_idx": 0, "annotation": 1, "start": 6, "end": 10},
    ]
    result = tm.pairwise_annotation_counts(rows, max_distance=5)
    # Hand value: the gap between [0, 4) and [6, 10) is 2 bases, so bin 2 counts the pair once in
    # each direction when symmetric.
    counts = result["counts"]
    assert counts[0][1][2] == 1 and counts[1][0][2] == 1
    assert sum(sum(sum(row) for row in plane) for plane in counts) == 2
    table = torch.tensor(
        [[0, 0, 0, 4], [0, 1, 6, 10]],
        dtype=torch.int64,
    )
    direct = annotate.pairwise_annotations_spacing(table, max_distance=5, dtype=torch.int64)
    assert counts == direct.tolist()


def test_tm2_pairwise_annotation_counts_rejects_overlap_and_bad_symmetric() -> None:
    overlapping = [
        {"example_idx": 0, "annotation": 0, "start": 0, "end": 6},
        {"example_idx": 0, "annotation": 1, "start": 4, "end": 8},
    ]
    with pytest.raises(BackendExecutionError):
        tm.pairwise_annotation_counts(overlapping, max_distance=5)
    with pytest.raises(ConfigurationError):
        tm.pairwise_annotation_counts([[0, 0]], symmetric="yes")


def test_tm2_screen_sequence_design_reaches_zero_count_and_matches_upstream() -> None:
    design = pytest.importorskip("tangermeme.design")
    torch = _tm2_torch()
    model = _tm2_count_model(0)
    result = tm.screen_sequence_design(model, 6, [0.0], n_best=1, max_iter=100, random_state=0)
    # Hand value: the output is the count of A and the target is 0, so a zero squared error means
    # the screened sequence has no A at all.
    assert "A" not in result["sequences"][0]
    assert result["squared_error"] == [0.0]
    assert result["predictions"] == [[0.0]]
    direct = design.screen(
        model,
        (4, 6),
        y=torch.tensor([[0.0]]),
        n_best=1,
        max_iter=100,
        random_state=0,
        device="cpu",
    )
    assert result["sequences"] == ["".join("ACGT"[i] for i in direct[0].argmax(dim=0).tolist())]


def test_tm2_screen_sequence_design_keeps_best_first_for_two_candidates() -> None:
    design = pytest.importorskip("tangermeme.design")
    torch = _tm2_torch()

    class Weighted(torch.nn.Module):
        def __init__(self):
            super().__init__()
            weights = numpy.random.RandomState(0).normal(size=(4, 12))
            self.weights = torch.tensor(weights, dtype=torch.float32)

        def forward(self, x):
            return (x * self.weights).sum(dim=(1, 2)).unsqueeze(1)

    model = Weighted().eval()
    result = tm.screen_sequence_design(model, 12, [0.0], n_best=2, max_iter=3, random_state=0)
    assert len(result["sequences"]) == 2
    assert result["squared_error"] == sorted(result["squared_error"])
    direct = design.screen(
        model,
        (4, 12),
        y=torch.tensor([[0.0]]),
        n_best=2,
        max_iter=3,
        random_state=0,
        device="cpu",
    )
    expected = ["".join("ACGT"[i] for i in row.argmax(dim=0).tolist()) for row in direct]
    assert result["sequences"] == expected


def test_tm2_greedy_marginalize_design_places_motifs_at_max_spacing() -> None:
    design = pytest.importorskip("tangermeme.design")
    utils = pytest.importorskip("tangermeme.utils")
    torch = _tm2_torch()
    model = _tm2_count_model(3)
    background = "A" * 60
    # Hand value: the model counts T and the target is 2, so one TT is needed. The search starts at
    # 60 // 2 - 1 = 29.
    single = tm.greedy_marginalize_design(
        model, background, [2.0], ["TT"], reverse_complement=False
    )
    assert single["sequence"] == "TT"
    # Hand value: target 4 needs a second TT. Its search window starts at 29 - 3 - 2 = 24, the
    # first position that gains, so the two motifs sit 3 bases apart.
    double = tm.greedy_marginalize_design(
        model, background, [4.0], ["TT"], reverse_complement=False, max_spacing=3
    )
    assert double["sequence"] == "TTNNNTT"
    direct = design.greedy_marginalize(
        model,
        utils.one_hot_encode(background).float()[None],
        torch.tensor([[4.0]]),
        ["TT"],
        reverse_complement=False,
        max_spacing=3,
        device="cpu",
    )
    columns = direct.sum(dim=0).tolist()
    indices = direct.argmax(dim=0).tolist()
    text = "".join(
        "ACGT"[i] if value > 0.5 else "N" for i, value in zip(indices, columns, strict=True)
    )
    assert text == double["sequence"]


def test_tm2_greedy_marginalize_design_rejects_short_background() -> None:
    with pytest.raises(ConfigurationError):
        tm.greedy_marginalize_design(_tm2_count_model(3), "A" * 10, [2.0], ["TT"], max_spacing=12)


def test_tm2_marginalize_annotations_inserts_span_at_centre() -> None:
    marginalize = pytest.importorskip("tangermeme.marginalize")
    torch = _tm2_torch()
    model = _tm2_count_model(1)
    sources = ["TTCCCCTTTT"]
    backgrounds = ["AAAAAAAAAA"]
    annotations = [{"example_idx": 0, "start": 2, "end": 6}]
    result = tm.marginalize_annotations(model, sources, backgrounds, annotations)
    # Hand value: CCCC (4 bases) goes to position 10 // 2 - 4 // 2 = 3 of AAAAAAAAAA, so the model
    # that counts C goes from 0 to 4.
    assert result["insert_start"] == [3]
    assert result["before"] == [[[0.0]]]
    assert result["after"] == [[[4.0]]]
    direct = marginalize.marginalize_annotations(
        model,
        _tm2_one_hot(sources),
        _tm2_one_hot(backgrounds),
        torch.tensor([[0, 2, 6]]),
    )
    _tm2_assert_close(result["before"], direct.y_befores.tolist())
    _tm2_assert_close(result["after"], direct.y_afters.tolist())


def test_tm2_pisa_attribution_matches_upstream_and_keeps_global_rng() -> None:
    pisa = pytest.importorskip("tangermeme.pisa")
    torch = _tm2_torch()
    sequences = _tm2_random_sequences(2, 30, seed=3)
    model = _tm2_linear_model(30)
    numpy_before = numpy.random.get_state()
    torch_before = torch.get_rng_state()
    result = tm.pisa_attribution(model, sequences, n_shuffles=3, random_state=0)
    numpy_after = numpy.random.get_state()
    assert numpy.array_equal(numpy_before[1], numpy_after[1]) and numpy_before[2] == numpy_after[2]
    assert torch.equal(torch_before, torch.get_rng_state())
    assert result["shape"] == [2, 2, 4, 30]
    attributions = numpy.array(result["attributions"])
    observed = _tm2_one_hot(sequences).numpy()
    # Hand value: attributions are multiplied by the observed one-hot bases, so every base that is
    # absent at a position has an attribution of exactly zero.
    absent = numpy.broadcast_to(observed[:, None, :, :] == 0, attributions.shape)
    assert numpy.all(attributions[absent] == 0)
    direct = pisa.pisa(model, _tm2_one_hot(sequences), n_shuffles=3, random_state=0, device="cpu")
    assert numpy.allclose(attributions, direct.numpy(), atol=1e-6)


def test_tm2_tfmodisco_seqlets_attribution_is_sum_of_trimmed_span() -> None:
    seqlet = pytest.importorskip("tangermeme.seqlet")
    torch = _tm2_torch()
    rng = numpy.random.RandomState(0)
    attr = rng.normal(0, 0.1, size=(3, 60))
    attr[0, 20:41] += 1.0
    attr[1, 25:46] += 1.0
    attr[2, 10:31] += 1.0
    result = tm.tfmodisco_seqlets(attr.tolist())
    assert [row["example_idx"] for row in result["seqlets"]] == [0, 1, 2]
    for row in result["seqlets"]:
        # Hand value: seqlets are 41 bases long for a window of 21, so 10 bases are trimmed from
        # each end and the sum covers the planted block at 20 to 41 of its example.
        trim = (row["end"] - row["start"] - 21) // 2
        expected = attr[row["example_idx"], row["start"] + trim : row["end"] - trim].sum()
        assert row["attribution"] == pytest.approx(expected, rel=1e-5)
    direct = seqlet.tfmodisco_seqlets(torch.tensor(attr, dtype=torch.float32))
    assert [[row["example_idx"], row["start"], row["end"]] for row in result["seqlets"]] == direct[
        ["example_idx", "start", "end"]
    ].values.tolist()


def test_tm2_tfmodisco_seqlets_rejects_window_longer_than_rows() -> None:
    with pytest.raises(ConfigurationError):
        tm.tfmodisco_seqlets([[0.1, 0.2, 0.3]], window_size=5)


def test_tm2_new_functions_are_exported_in_sorted_all() -> None:
    names = [
        "ablate_annotations",
        "annotate_seqlets",
        "count_seqlet_annotations",
        "greedy_marginalize_design",
        "marginalize_annotations",
        "pairwise_annotation_counts",
        "pisa_attribution",
        "screen_sequence_design",
        "tfmodisco_seqlets",
    ]
    assert set(names) <= set(tm.__all__)
    assert list(tm.__all__) == sorted(tm.__all__)
    for name in names:
        assert callable(getattr(tm, name))


_SA_SIZES = {"chr1": 3000, "chr2": 2000}
_SA_GROUPS = {"cellA": "G1", "cellB": "G2", "cellC": "G1"}
# Hand counts: cellA (3) and cellC (2) form G1; cellB (2) forms G2. Coordinates are distinct.
_SA_FRAGMENT_LINES = (
    "chr1\t100\t300\tcellA\t1",
    "chr1\t400\t600\tcellA\t1",
    "chr1\t700\t900\tcellA\t1",
    "chr2\t50\t250\tcellB\t1",
    "chr1\t1000\t1200\tcellB\t1",
    "chr2\t1200\t1400\tcellC\t1",
    "chr1\t2000\t2200\tcellC\t1",
)
_SA_CONTACTS = (
    "cellA\tchr1\t10\tchr1\t40\t1\n"
    "cellA\tchr1\t120\tchr1\t950\t2\n"
    "cellB\tchr2\t300\tchr2\t700\t1\n"
    "cellB\tchr1\t50\tchr2\t800\t1\n"
)
_SA_VALUES_A = "chrom\tpos\tmethyl\tunmethyl\nchr1\t10\t3\t7\nchr1\t20\t0\t4\nchr2\t5\t2\t0\n"
_SA_VALUES_B = "chrom\tpos\tmethyl\tunmethyl\nchr1\t10\t1\t1\n"
_SA_GTF = (
    'chr1\tt\tgene\t1000\t1500\t.\t+\t.\tgene_id "G1"; gene_name "GATA1";\n'
    "chr1\tt\ttranscript\t1000\t1500\t.\t+\t.\t"
    'gene_id "G1"; transcript_id "T1"; gene_name "GATA1";\n'
    'chr1\tt\texon\t1000\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GATA1";\n'
)
_SA_NETWORK_GTF = (
    'chr1\tt\tgene\t900\t1500\t.\t+\t.\tgene_id "G1"; gene_name "GATA1";\n'
    "chr1\tt\ttranscript\t900\t1500\t.\t+\t.\t"
    'gene_id "G1"; transcript_id "T1"; gene_name "GATA1";\n'
    'chr1\tt\texon\t900\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GATA1";\n'
    'chr2\tt\tgene\t400\t800\t.\t+\t.\tgene_id "G2"; gene_name "OTHER";\n'
    'chr2\tt\ttranscript\t400\t800\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "OTHER";\n'
    'chr2\tt\texon\t400\t800\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "OTHER";\n'
)
_SA_MOTIF = "GATAAGGATC"
_SA_MEME = (
    "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n"
    "A 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF GATA1\n"
    "letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0\n"
    + "\n".join(
        " " + " ".join("1.0" if "ACGT"[k] == base else "0.0" for k in range(4))
        for base in _SA_MOTIF
    )
    + "\n"
)


def _sa_write_fragments(path: Path) -> Path:
    path.write_text("\n".join(_SA_FRAGMENT_LINES) + "\n", encoding="utf-8")
    return path


def _sa_write_bam(path: Path, reads: Sequence[tuple[str, str, int, int, int]]) -> Path:
    """Write read pairs as (name, barcode, reference id, 0-based start, template length)."""

    pysam = pytest.importorskip("pysam")
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 3000}, {"SN": "chr2", "LN": 2000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, barcode, tid, start, tlen in reads:
            for read1 in (True, False):
                record = pysam.AlignedSegment()
                record.query_name = name
                record.query_sequence = "A" * 50
                record.flag = 0x1 | 0x2 | (0x40 if read1 else 0x80) | (0 if read1 else 0x10)
                record.reference_id = tid
                record.reference_start = start if read1 else start + tlen - 50
                record.mapping_quality = 60
                record.cigarstring = "50M"
                record.next_reference_id = tid
                record.next_reference_start = start + tlen - 50 if read1 else start
                record.template_length = tlen if read1 else -tlen
                record.query_qualities = pysam.qualitystring_to_array("I" * 50)
                record.set_tag("CB", barcode, value_type="Z")
                out.write(record)
    return path


def _sa_write_tenx_bam(path: Path) -> Path:
    """30 cells with 300 pairs inside a peak and 120 background barcodes with 3 pairs each."""

    rng = random.Random(7)
    reads = []
    for index in range(30):
        cell = f"C{index:03d}"
        for k in range(300):
            reads.append((f"{cell}_{k}", cell, 0, rng.randint(900, 1400), rng.randint(160, 240)))
    for index in range(120):
        cell = f"B{index:03d}"
        for k in range(3):
            tid = rng.choice([0, 1])
            reads.append((f"{cell}_{k}", cell, tid, rng.randint(0, 1700), rng.randint(160, 240)))
    rng.shuffle(reads)
    return _sa_write_bam(path, reads)


def _sa_write_clustered_fragments(path: Path) -> Path:
    """60 cellA fragments on chr1 near 1000, 60 cellB on chr2 near 500, 20 random cellC."""

    rng = random.Random(3)
    rows = []
    for _ in range(60):
        start = 1000 + rng.randint(0, 200)
        rows.append(f"chr1\t{start}\t{start + 150}\tcellA\t1")
    for _ in range(60):
        start = 500 + rng.randint(0, 100)
        rows.append(f"chr2\t{start}\t{start + 150}\tcellB\t1")
    for _ in range(20):
        start = rng.randint(0, 1800)
        rows.append(f"chr1\t{start}\t{start + 150}\tcellC\t1")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _sa_motif_inputs(directory: Path) -> tuple[Path, Path, Path]:
    """A genome with the motif planted at chr1:1100-1110, a MEME file and a GTF."""

    rng = random.Random(11)
    seq1 = [rng.choice("ACGT") for _ in range(3000)]
    seq1[1100:1110] = list(_SA_MOTIF)
    seq2 = [rng.choice("ACGT") for _ in range(2000)]
    fasta = directory / "genome.fa"
    with fasta.open("w", encoding="utf-8") as handle:
        for name, sequence in (("chr1", "".join(seq1)), ("chr2", "".join(seq2))):
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")
    meme = directory / "motifs.meme"
    meme.write_text(_SA_MEME, encoding="utf-8")
    gtf = directory / "network.gtf"
    gtf.write_text(_SA_NETWORK_GTF, encoding="utf-8")
    return gtf, meme, fasta


def test_sa_export_fragment_files_writes_hand_counts_and_matches_upstream(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    fragments = _sa_write_fragments(tmp_path / "fragments.tsv")
    out = tmp_path / "exported"
    out.mkdir()
    result = sa.export_fragment_files(fragments, _SA_SIZES, _SA_GROUPS, out, prefix="sa_")
    assert [(item["file"], item["fragment_count"]) for item in result["groups"]] == [
        ("sa_G1.bed", 5),
        ("sa_G2.bed", 2),
    ]
    direct_dir = tmp_path / "direct"
    direct_dir.mkdir()
    adata = snap.pp.import_fragments(
        str(fragments),
        _SA_SIZES,
        is_paired=True,
        min_num_fragments=1,
        file=str(tmp_path / "direct.h5ad"),
        sorted_by_barcode=False,
        n_jobs=1,
    )
    try:
        labels = [_SA_GROUPS[barcode] for barcode in adata.obs_names]
        snap.ex.export_fragments(
            adata,
            groupby=labels,
            out_dir=str(direct_dir),
            prefix="sa_",
            suffix=".bed",
            compression=None,
        )
    finally:
        adata.close()
    for name in ("sa_G1.bed", "sa_G2.bed"):
        assert (out / name).read_text(encoding="utf-8") == (direct_dir / name).read_text(
            encoding="utf-8"
        )


def test_sa_export_fragment_files_refuses_missing_barcode_and_existing_output(
    tmp_path: Path,
) -> None:
    fragments = _sa_write_fragments(tmp_path / "fragments.tsv")
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(ConfigurationError, match="missing"):
        sa.export_fragment_files(fragments, _SA_SIZES, {"cellA": "G1"}, out)
    assert list(out.iterdir()) == []
    sa.export_fragment_files(fragments, _SA_SIZES, _SA_GROUPS, out, verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        sa.export_fragment_files(fragments, _SA_SIZES, _SA_GROUPS, out)


def test_sa_export_coverage_tracks_matches_upstream_bedgraph(tmp_path: Path) -> None:
    snap = pytest.importorskip("snapatac2")
    fragments = _sa_write_fragments(tmp_path / "fragments.tsv")
    out = tmp_path / "tracks"
    out.mkdir()
    result = sa.export_coverage_tracks(
        fragments, _SA_SIZES, _SA_GROUPS, out, bin_size=100, output_format="bedgraph"
    )
    assert [item["file"] for item in result["tracks"]] == ["G1.bedgraph", "G2.bedgraph"]
    direct_dir = tmp_path / "direct"
    direct_dir.mkdir()
    adata = snap.pp.import_fragments(
        str(fragments),
        _SA_SIZES,
        is_paired=True,
        min_num_fragments=1,
        file=str(tmp_path / "direct.h5ad"),
        sorted_by_barcode=False,
        n_jobs=1,
    )
    try:
        labels = [_SA_GROUPS[barcode] for barcode in adata.obs_names]
        snap.ex.export_coverage(
            adata,
            groupby=labels,
            bin_size=100,
            normalization=None,
            max_frag_length=None,
            counting_strategy="fragment",
            out_dir=str(direct_dir),
            prefix="",
            suffix=".bedgraph",
            output_format="bedgraph",
            compression=None,
            n_jobs=1,
        )
    finally:
        adata.close()
    for name in ("G1.bedgraph", "G2.bedgraph"):
        assert (out / name).read_text(encoding="utf-8") == (direct_dir / name).read_text(
            encoding="utf-8"
        )


def test_sa_export_coverage_tracks_bigwig_lists_chromosomes_with_fragments(
    tmp_path: Path,
) -> None:
    pyBigWig = pytest.importorskip("pyBigWig")
    fragments = _sa_write_fragments(tmp_path / "fragments.tsv")
    out = tmp_path / "tracks"
    out.mkdir()
    result = sa.export_coverage_tracks(
        fragments,
        _SA_SIZES,
        _SA_GROUPS,
        out,
        bin_size=100,
        output_format="bigwig",
        prefix="sa_",
    )
    assert result["track_count"] == 2
    bigwig = pyBigWig.open(str(out / "sa_G1.bw"))
    try:
        # G1 holds fragments on chr1 and chr2 (cellC), so both chromosomes are in the header.
        assert bigwig.chroms() == {"chr1": 3000, "chr2": 2000}
    finally:
        bigwig.close()


def test_sa_fragments_from_bam_matches_hand_fragments_and_upstream(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    bam = _sa_write_bam(
        tmp_path / "reads.bam",
        [
            ("r1", "AAAA", 0, 100, 300),
            ("r2", "AAAA", 0, 600, 200),
            ("r3", "CCCC", 1, 50, 200),
        ],
    )
    out = tmp_path / "fragments.tsv"
    result = sa.fragments_from_bam(bam, out)
    # Hand values: a fragment starts 4 bases after the 0-based read start and ends 5 bases before
    # the template end (start + template length), so 100 + 300 - 5 = 395.
    assert out.read_text(encoding="utf-8").splitlines() == [
        "chr1\t104\t395\tAAAA\t1",
        "chr1\t604\t795\tAAAA\t1",
        "chr2\t54\t245\tCCCC\t1",
    ]
    assert result["fragment_count"] == 3
    assert result["barcode_count"] == 2
    assert result["qc"]["sequenced_read_pairs"] == 3.0
    direct = tmp_path / "direct.tsv"
    snap.pp.make_fragment_file(
        str(bam),
        str(direct),
        is_paired=True,
        barcode_tag="CB",
        min_mapq=30,
        chrM=["chrM", "M"],
        source=None,
        compression=None,
    )
    assert direct.read_text(encoding="utf-8") == out.read_text(encoding="utf-8")


def test_sa_import_contact_pairs_matches_upstream_and_hand_counts(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    contacts = tmp_path / "contacts.tsv"
    contacts.write_text(_SA_CONTACTS, encoding="utf-8")
    sizes = {"chr1": 1000, "chr2": 1000}
    result = sa.import_contact_pairs(contacts, sizes, tmp_path / "contacts.h5ad", bin_size=500)
    # Hand values: four contact pairs with counts 1 + 2 + 1 + 1; two cells; two 500-base bins per
    # chromosome give four bins and 4 x 4 = 16 bin pairs.
    assert result["cell_count"] == 2
    assert result["contact_pair_count"] == 4
    assert result["contact_count"] == 5
    assert result["feature_count"] == 16
    assert (tmp_path / "contacts.h5ad").is_file()
    direct = snap.pp.import_contacts(
        str(contacts), sizes, sorted_by_barcode=False, bin_size=500, file=None
    )
    assert direct.n_obs == result["cell_count"]
    assert direct.n_vars == result["feature_count"]
    assert int(direct.X.sum()) == result["contact_count"]


def test_sa_import_contact_pairs_rejects_position_outside_chromosome(
    tmp_path: Path,
) -> None:
    contacts = tmp_path / "contacts.tsv"
    contacts.write_text("cellA\tchr1\t10\tchr2\t1000\t1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="pos"):
        sa.import_contact_pairs(contacts, {"chr1": 1000, "chr2": 1000}, tmp_path / "bad.h5ad")
    assert not (tmp_path / "bad.h5ad").exists()


def test_sa_import_base_pair_values_counts_match_files_and_upstream(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    directory = tmp_path / "values"
    directory.mkdir()
    (directory / "cellA.tsv").write_text(_SA_VALUES_A, encoding="utf-8")
    (directory / "cellB.tsv").write_text(_SA_VALUES_B, encoding="utf-8")
    sizes = {"chr1": 100, "chr2": 50}
    result = sa.import_base_pair_values(directory, sizes, tmp_path / "values.h5ad")
    # Hand values: cellA has three base rows and cellB has one.
    assert result["cells"] == [
        {"cell": "cellA", "value_count": 3},
        {"cell": "cellB", "value_count": 1},
    ]
    assert result["total_value_count"] == 4
    assert result["reference_length"] == 150
    direct = snap.pp.import_values(directory, sizes, file=None)
    counts = dict(zip(direct.obs_names, direct.obs["num_values"].to_list(), strict=True))
    assert counts == {"cellA": 3, "cellB": 1}


def test_sa_import_base_pair_values_rejects_bad_header_and_extra_files(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "values"
    directory.mkdir()
    (directory / "cellA.tsv").write_text("pos\tchrom\tmethyl\tunmethyl\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="header"):
        sa.import_base_pair_values(directory, {"chr1": 100}, tmp_path / "bad.h5ad")
    (directory / "cellA.tsv").write_text(_SA_VALUES_B, encoding="utf-8")
    (directory / "notes.txt").write_text("not a cell\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"\.tsv"):
        sa.import_base_pair_values(directory, {"chr1": 100}, tmp_path / "bad.h5ad")


def test_sa_tenx_atac_metrics_matches_upstream_recipe_and_hand_counts(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    bam = _sa_write_tenx_bam(tmp_path / "reads.bam")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_SA_GTF, encoding="utf-8")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t900\t1400\tp1\t10\t.\n", encoding="utf-8")
    result = sa.tenx_atac_metrics(
        bam, tmp_path / "fragments.tsv", tmp_path / "qc.h5ad", _SA_SIZES, gtf, peaks
    )
    metrics = result["metrics"]
    # Hand values: 30 x 300 + 120 x 3 = 9360 read pairs; one peak line in the BED file.
    assert metrics["Sequencing"]["Sequenced_read_pairs"] == 9360
    assert metrics["Targeting"]["Number_of_peaks"] == 1
    direct = snap.pp.recipe_10x_metrics(
        str(bam),
        str(tmp_path / "direct.tsv"),
        str(tmp_path / "direct.h5ad"),
        peaks=str(peaks),
        barcode_tag="CB",
        source=None,
        chrom_sizes=_SA_SIZES,
        gene_anno=str(gtf),
        sorted_by_barcode=False,
        n_jobs=1,
    )
    assert metrics == json.loads(json.dumps(direct, default=float))


def test_sa_call_peaks_macs3_matches_upstream_and_cluster_bounds(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    fragments = _sa_write_clustered_fragments(tmp_path / "clustered.tsv")
    before = tempfile.tempdir
    result = sa.call_peaks_macs3(fragments, _SA_SIZES, _SA_GROUPS)
    assert tempfile.tempdir == before
    assert [(item["group"], item["peak_count"]) for item in result["groups"]] == [
        ("G1", 1),
        ("G2", 1),
    ]
    # Hand bounds: cellA fragments cover chr1:1000-1350 and cellB fragments cover chr2:500-650.
    g1_peak, g2_peak = result["groups"][0]["peaks"][0], result["groups"][1]["peaks"][0]
    assert g1_peak["chrom"] == "chr1" and g1_peak["start"] <= 1000 and g1_peak["end"] >= 1350
    assert g2_peak["chrom"] == "chr2" and g2_peak["start"] <= 500 and g2_peak["end"] >= 650
    try:
        adata = snap.pp.import_fragments(
            str(fragments),
            _SA_SIZES,
            is_paired=True,
            min_num_fragments=1,
            file=str(tmp_path / "macs.h5ad"),
            sorted_by_barcode=False,
            n_jobs=1,
        )
        try:
            labels = [_SA_GROUPS[barcode] for barcode in adata.obs_names]
            frames = snap.tl.macs3(
                adata,
                groupby=labels,
                qvalue=0.05,
                shift=-100,
                extsize=200,
                inplace=False,
                n_jobs=1,
            )
            direct = {
                name: frame.sort(["chrom", "start", "end"]).to_dicts()
                for name, frame in frames.items()
            }
        finally:
            adata.close()
    finally:
        # Upstream MACS3 leaves tempfile.tempdir pointing at its removed scratch directory.
        tempfile.tempdir = before
    assert [item["peaks"] for item in result["groups"]] == [direct["G1"], direct["G2"]]


def test_sa_add_tf_binding_edges_matches_upstream_and_hand_motif_site(
    tmp_path: Path,
) -> None:
    snap = pytest.importorskip("snapatac2")
    gtf, meme, fasta = _sa_motif_inputs(tmp_path)
    regions = ["chr1:1000-1200", "chr1:2000-2200", "chr2:500-700"]
    result = sa.add_tf_binding_edges(regions, gtf, meme, fasta, pvalue=1e-5)
    # Hand values: the motif sits at chr1:1100-1110, inside chr1:1000-1200 only; there are three
    # region nodes and two gene nodes.
    assert result["motif_count"] == 1
    assert result["region_count"] == 3
    assert result["gene_count"] == 2
    assert result["binding_edges"] == [{"motif": "GATA1", "region": "chr1:1000-1200"}]
    network = snap.tl.init_network_from_annotation(
        regions,
        str(gtf),
        upstream=250_000,
        downstream=250_000,
        id_type="gene_name",
        coding_gene_only=True,
    )
    snap.tl.add_tf_binding(
        network,
        motifs=snap.read_motifs(str(meme)),
        genome_fasta=str(fasta),
        pvalue=1e-5,
    )
    direct = sorted(
        (network[src].id, network[dst].id)
        for src, dst, _data in network.edge_index_map().values()
        if network[src].type == "motif" and network[dst].type == "region"
    )
    assert direct == [(item["motif"], item["region"]) for item in result["binding_edges"]]


def test_sa_add_tf_binding_edges_rejects_bad_region_and_pvalue(tmp_path: Path) -> None:
    gtf, meme, fasta = _sa_motif_inputs(tmp_path)
    with pytest.raises(ConfigurationError, match="chrom:start-end"):
        sa.add_tf_binding_edges(["chr1:200-100"], gtf, meme, fasta)
    with pytest.raises(ConfigurationError, match="pvalue"):
        sa.add_tf_binding_edges(["chr1:1000-1200"], gtf, meme, fasta, pvalue=0.0)


def test_sa_snapatac2_missing_package_names_the_extra(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fragments = _sa_write_fragments(tmp_path / "fragments.tsv")
    real_import = importlib.import_module

    def fake_import(name: str, package: str | None = None):  # type: ignore[no-untyped-def]
        if name.startswith("snapatac2"):
            raise ImportError("simulated missing snapatac2")
        return real_import(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", fake_import)
    with pytest.raises(BackendUnavailableError) as info:
        sa.export_fragment_files(fragments, _SA_SIZES, _SA_GROUPS, tmp_path)
    assert "bio-snapatac2" in str(info.value) + str(getattr(info.value, "hint", ""))


def test_review_g18i_public_names_are_sorted_and_take_verbose() -> None:
    import importlib
    import inspect

    for module_name in (
        "dnakit.regulation.gimmemotifs_review_functions",
        "dnakit.regulation.tangermeme_review_functions",
        "dnakit.regulation.snapatac2_review_functions",
    ):
        module = importlib.import_module(module_name)
        names = list(module.__all__)
        assert names == sorted(names), module_name
        for name in names:
            parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, (module_name, name)
            assert parameter.default is True, (module_name, name)
