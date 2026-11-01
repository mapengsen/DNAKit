from __future__ import annotations

import random

import pytest

pytest.importorskip("gimmemotifs")
pytest.importorskip("genomepy")
np = pytest.importorskip("numpy")

from gimmemotifs.fasta import Fasta  # noqa: E402
from gimmemotifs.motif import read_motifs  # noqa: E402
from gimmemotifs.scanner import Scanner  # noqa: E402

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.regulation import gimmemotifs_p28_data_functions as module  # noqa: E402

PFM_FILE = (
    ">acgt\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
    ">weak\n6 2 1 1\n1 6 2 1\n1 1 6 2\n2 1 1 6\n3 3 2 2\n"
)


def _random_fasta(path, count=40, length=120, seed=0):
    rng = random.Random(seed)
    path.write_text(
        "".join(
            f">r{i}\n{''.join(rng.choice('ACGT') for _ in range(length))}\n" for i in range(count)
        )
    )
    return path


def _genome(path):
    rng = random.Random(5)
    path.write_text(
        "".join(f">chr{c}\n{''.join(rng.choice('ACGT') for _ in range(4000))}\n" for c in (1, 2))
    )
    return path


def test_hardmask_and_random_subset(tmp_path):
    source = tmp_path / "soft.fa"
    source.write_text(">a\nACGTacgtnNAC\n>b\nggggCCCC\n")
    result = module.gimmemotifs_hardmask_fasta(source, tmp_path / "hard.fa", verbose=False)
    assert (tmp_path / "hard.fa").read_text() == ">a\nACGTNNNNNNAC\n>b\nNNNNCCCC\n"
    assert result["masked_bases"] == 9
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_hardmask_fasta(source, tmp_path / "hard.fa", verbose=False)

    fasta = _random_fasta(tmp_path / "r.fa")
    first = module.gimmemotifs_random_fasta_subset(fasta, 5, seed=2, verbose=False)
    again = module.gimmemotifs_random_fasta_subset(fasta, 5, seed=2, verbose=False)
    assert first["sequences"] == again["sequences"] and len(first["sequences"]) == 5
    random.seed(2)
    expected = Fasta(str(fasta)).get_random(5)
    assert first["sequences"] == dict(expected.items())
    windows = module.gimmemotifs_random_fasta_subset(
        fasta, 3, length=30, seed=1, output_fasta=tmp_path / "w.fa", verbose=False
    )
    assert list(windows["sequences"]) == ["random1", "random2", "random3"]
    assert all(len(seq) == 30 for seq in windows["sequences"].values())
    assert (tmp_path / "w.fa").read_text().count(">") == 3
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_random_fasta_subset(fasta, 41, verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_random_fasta_subset(fasta, 3, length=500, verbose=False)


def test_sequence_dict(tmp_path):
    fasta = _random_fasta(tmp_path / "r.fa", count=3)
    result = module.gimmemotifs_sequence_dict(fasta, minsize=100, verbose=False)
    assert list(result["sequences"]) == ["r0", "r1", "r2"]
    assert module.gimmemotifs_sequence_dict(["ACGT", "GGAA"], verbose=False)["sequences"] == {
        "1": "ACGT",
        "2": "GGAA",
    }
    with pytest.raises(Exception, match="shorter"):
        module.gimmemotifs_sequence_dict(fasta, minsize=500, verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_sequence_dict([], verbose=False)


STATS = {
    "m_1": {"bg": {"roc_auc": 0.9, "recall_at_fdr": 0.5}},
    "m_2": {"bg": {"roc_auc": 0.6, "recall_at_fdr": 0.7}},
    "m_3": {"bg": {"roc_auc": 0.7, "recall_at_fdr": 0.1}},
}


def test_rank_and_write_stats(tmp_path):
    ranked = module.gimmemotifs_rank_motif_stats(STATS, verbose=False)
    assert ranked["mean_rank"] == {"m_1": 2.5, "m_2": 2.0, "m_3": 1.5}
    assert ranked["order"] == ["m_1", "m_2", "m_3"]
    single = module.gimmemotifs_rank_motif_stats(STATS, metrics=["recall_at_fdr"], verbose=False)
    assert single["order"] == ["m_2", "m_1", "m_3"]
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_rank_motif_stats(STATS, metrics=["nope"], verbose=False)

    written = module.gimmemotifs_write_motif_stats(
        STATS, str(tmp_path / "stats.{}.txt"), header="# run\n", verbose=False
    )
    lines = (tmp_path / "stats.bg.txt").read_text().splitlines()
    assert lines[:2] == ["# run", "Motif\trecall_at_fdr\troc_auc"]
    assert lines[2] == "m\t0.5\t0.9"
    assert written["outputs"] == {"bg": str(tmp_path / "stats.bg.txt")}
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_write_motif_stats(STATS, str(tmp_path / "stats.{}.txt"), verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_write_motif_stats(STATS, str(tmp_path / "plain.txt"), verbose=False)


def test_thresholds_and_counts(tmp_path):
    pfm = tmp_path / "m.pfm"
    pfm.write_text(PFM_FILE)
    background = _random_fasta(tmp_path / "bg.fa", count=60, length=100, seed=1)
    sequences = _random_fasta(tmp_path / "seq.fa", count=10, length=100, seed=2)

    result = module.gimmemotifs_motif_thresholds(pfm, background, fpr=0.5, seed=3, verbose=False)
    again = module.gimmemotifs_motif_thresholds(pfm, background, fpr=0.5, seed=3, verbose=False)
    assert result["thresholds"] == again["thresholds"]
    weak = result["thresholds"]["weak"]
    scanner = Scanner(ncpus=1)
    scanner.set_motifs(str(pfm))
    best = scanner._threshold_from_seqs(read_motifs(str(pfm)), Fasta(str(background)))["weak"]
    assert best.min() <= weak["threshold"] <= best.max()
    assert weak["threshold"] < weak["max_score"] and 0 < weak["fraction"] < 1
    # the threshold is the 50% quantile of the best background scores
    scanner.set_threshold(threshold=weak["fraction"])
    reaching = sum(
        1 for _, counts in scanner.count(Fasta(str(background)), nreport=1) if counts[1] > 0
    )
    assert 0 < reaching <= int((best >= weak["threshold"] - 1e-9).sum())
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_motif_thresholds(pfm, background, fpr=1.5, verbose=False)

    counts = module.gimmemotifs_count_motif_matches(sequences, pfm, threshold=0.8, verbose=False)
    scanner = Scanner(ncpus=1)
    scanner.set_motifs(str(pfm))
    scanner.set_threshold(threshold=0.8)
    expected = [list(c) for _, c in scanner.count(Fasta(str(sequences)))]
    assert [list(row.values()) for row in counts["counts"].values()] == expected
    assert list(counts["total"].values()) == list(np.sum(np.array(expected), 0))
    assert counts["motifs"] == ["acgt", "weak"]
    by_fpr = module.gimmemotifs_count_motif_matches(
        sequences, pfm, background_fasta=background, fpr=0.5, seed=3, verbose=False
    )
    assert by_fpr["cutoff_mode"] == "fpr" and by_fpr["total"]["weak"] >= 0
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_count_motif_matches(sequences, pfm, fpr=0.1, verbose=False)


def test_prepare_denovo_input(tmp_path):
    fasta = _random_fasta(tmp_path / "in.fa", count=20, length=60)
    result = module.gimmemotifs_prepare_denovo_input(
        fasta, tmp_path / "fa", fraction=0.5, seed=1, verbose=False
    )
    assert result["records_per_file"] == {
        "localization.fa": 10,
        "prediction.fa": 10,
        "validation.fa": 10,
    }
    same = module.gimmemotifs_prepare_denovo_input(
        fasta, tmp_path / "fa2", fraction=0.5, seed=1, verbose=False
    )
    assert (tmp_path / "fa" / "prediction.fa").read_text() == (
        tmp_path / "fa2" / "prediction.fa"
    ).read_text()
    assert same["parameters"]["size"] == 200

    genome = _genome(tmp_path / "genome.fa")
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "".join(f"chr1\t{100 + 150 * i}\t{300 + 150 * i}\tp{i}\t{100 - i}\n" for i in range(10))
    )
    result = module.gimmemotifs_prepare_denovo_input(
        bed,
        tmp_path / "bed",
        input_type="bed",
        genome_fasta=genome,
        size=100,
        fraction=0.5,
        lsize=140,
        seed=1,
        verbose=False,
    )
    assert result["records_per_file"]["input.bed"] == 10
    assert result["records_per_file"]["prediction.fa"] == 5
    first = (tmp_path / "bed" / "prediction.fa").read_text().splitlines()[1]
    assert len(first) == 100
    assert len((tmp_path / "bed" / "localization.fa").read_text().splitlines()[1]) == 140

    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text(
        "".join(
            f"chr2\t{100 + 150 * i}\t{300 + 150 * i}\tp{i}\t0\t.\t5\t-1\t-1\t100\n"
            for i in range(8)
        )
    )
    result = module.gimmemotifs_prepare_denovo_input(
        peaks,
        tmp_path / "np",
        input_type="narrowpeak",
        genome_fasta=genome,
        size=60,
        fraction=0.5,
        seed=1,
        verbose=False,
    )
    assert result["records_per_file"]["prediction.fa"] == 4
    bed_text = (tmp_path / "np" / "input.bed").read_text().splitlines()[0].split("\t")
    assert (bed_text[1], bed_text[2]) == ("170", "230")
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_prepare_denovo_input(
            bed, tmp_path / "x", input_type="bed", verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_prepare_denovo_input(fasta, tmp_path / "fa", verbose=False)


class _Response:
    def __init__(self, text):
        self._lines = text.encode().splitlines(keepends=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        return iter(self._lines)


HOMER_TEXT = (
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\t6.04\t-1782\t0\t9805\n"
    "0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n0.010\t0.002\t0.965\t0.023\n"
    ">GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer\t6.1\t-12\t0\t9\n"
    "0.1\t0.1\t0.7\t0.1\n0.1\t0.7\t0.1\t0.1\n0.7\t0.1\t0.1\t0.1\n"
)


def test_download_and_annotate_homer_with_mock(tmp_path, monkeypatch):
    import gimmemotifs.db as db

    monkeypatch.setattr(db, "urlopen", lambda url: _Response(HOMER_TEXT))
    result = module.gimmemotifs_download_motif_database("homer", tmp_path / "db", verbose=False)
    assert sorted(result["outputs"]) == ["HOMER.motif2factors.txt", "HOMER.pfm"]
    pfm = tmp_path / "db" / "HOMER.pfm"
    assert pfm.read_text().count(">") == 2
    table = (tmp_path / "db" / "HOMER.motif2factors.txt").read_text().splitlines()
    assert table[0] == "Motif\tFactor\tEvidence\tCurated"
    assert table[1].endswith("\tPU.1\tChIP-seq\tY")
    annotated = module.gimmemotifs_annotate_motif_factors(
        "homer", pfm, output_dir=tmp_path / "out", verbose=False
    )
    assert annotated["n_motifs"] == 2 and annotated["n_annotated"] == 2
    assert [row[0][0] for row in annotated["annotation"].values()] == ["PU.1", "Gfi1b"]
    assert (tmp_path / "out" / "HOMER.motif2factors.txt").read_text() == "\n".join(table) + "\n"
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_download_motif_database("factorbook", tmp_path / "db", verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_annotate_motif_factors("rsat", pfm, verbose=False)


def test_annotate_hocomoco_and_jaspar_with_mock(tmp_path, monkeypatch):
    import gimmemotifs.db as db

    pfm = tmp_path / "h.pfm"
    pfm.write_text(">AHR_HUMAN.H11MO.0.B\n1 2 3 4\n4 3 2 1\n>MA0004.1_Arnt\n1 2 3 4\n")
    urls: list[str] = []

    def fake_open(url):
        urls.append(url)
        return _Response("AHR_HUMAN.H11MO.0.B\tAHR\nMA0004.1_Arnt\tARNT\n")

    monkeypatch.setattr(db, "urlopen", fake_open)
    result = module.gimmemotifs_annotate_motif_factors(
        "hocomoco", pfm, hocomoco_group="MOUSE", verbose=False
    )
    assert "MOUSE" in urls[0]
    assert result["annotation"]["AHR_HUMAN.H11MO.0.B"] == [["AHR", "ChIP-seq", "Y"]]

    monkeypatch.setattr(db, "get_jaspar_motif_info", lambda motif_id: {"type": "SELEX"})
    jaspar = module.gimmemotifs_annotate_motif_factors("jaspar", pfm, verbose=False)
    assert jaspar["annotation"]["MA0004.1_Arnt"] == [["Arnt", "SELEX", "Y"]]


def test_narrowpeak_to_bed(tmp_path):
    peaks = tmp_path / "p.narrowPeak"
    peaks.write_text(
        "# comment\n"
        "chr1\t100\t300\tp0\t0\t.\t12.5\t-1\t-1\t80\n"
        "chr2\t500\t900\tp1\t0\t.\t7.0\t-1\t-1\t-1\n"
    )
    plain = module.gimmemotifs_narrowpeak_to_bed(peaks, tmp_path / "plain.bed", verbose=False)
    assert plain["intervals"] == [["chr1", 100, 300, "12.5"], ["chr2", 500, 900, "7.0"]]
    summit = module.gimmemotifs_narrowpeak_to_bed(
        peaks, tmp_path / "summit.bed", size=60, verbose=False
    )
    assert summit["intervals"] == [["chr1", 150, 210, "12.5"], ["chr2", 670, 730, "7.0"]]
    assert (tmp_path / "summit.bed").read_text().startswith("chr1\t150\t210\t12.5")
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_narrowpeak_to_bed(peaks, tmp_path / "summit.bed", verbose=False)
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\t1\t2\n")
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_narrowpeak_to_bed(bad, tmp_path / "x.bed", verbose=False)
