from __future__ import annotations

import pytest

pytest.importorskip("gimmemotifs")
np = pytest.importorskip("numpy")

from gimmemotifs.comparison import MotifComparer  # noqa: E402
from gimmemotifs.fasta import Fasta  # noqa: E402
from gimmemotifs.motif import Motif  # noqa: E402

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.regulation import gimmemotifs_p28_motif_functions as module  # noqa: E402

PFM_A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
PFM_B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]
FASTA = ">s1\nTTACGTTTACGTAA\n>s2\nGGGGACGTCCCC\n>s3\nCCCCCCCCCC\n"


def _fasta(tmp_path):
    path = tmp_path / "s.fa"
    path.write_text(FASTA)
    return path


@pytest.mark.parametrize(
    ("kwargs", "match", "metric"),
    [
        ({}, "total", "wic"),
        ({"metric": "pcc"}, "total", "pcc"),
        ({"metric": "akl", "match": "partial"}, "partial", "akl"),
        ({"metric": "ssd", "match": "subtotal"}, "subtotal", "ssd"),
        ({"metric": "ed", "combine": "sum"}, "total", "ed"),
    ],
)
def test_compare_motif_pair_matches_upstream(kwargs, match, metric):
    result = module.gimmemotifs_compare_motif_pair(PFM_A, PFM_B, verbose=False, **kwargs)
    expected = MotifComparer().compare_motifs(
        Motif(PFM_A), Motif(PFM_B), match, metric, kwargs.get("combine", "mean")
    )
    assert result["value"] == pytest.approx(expected[0])
    assert (result["position"], result["strand"]) == (expected[1], expected[2])


def test_compare_motif_pair_pvalue_and_errors(tmp_path):
    result = module.gimmemotifs_compare_motif_pair(PFM_A, PFM_B, pvalue=True, verbose=False)
    assert result["value_kind"] == "pvalue"
    assert 0 <= result["value"] <= 1
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_compare_motif_pair(
            PFM_A, PFM_B, metric="akl", pvalue=True, verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_compare_motif_pair(
            PFM_A, PFM_B, metric="pcc", pvalue=True, verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_compare_motif_pair(PFM_A, PFM_B, metric="bogus", verbose=False)


def test_score_distribution_file_feeds_pvalue(tmp_path):
    pfm = tmp_path / "set.pfm"
    pfm.write_text(
        ">a\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
        ">b\n0 10 0 0\n0 0 10 0\n0 0 0 10\n10 0 0 0\n"
        ">c\n5 5 0 0\n0 5 5 0\n0 0 5 5\n5 0 0 5\n"
        ">d\n10 0 0 0\n10 0 0 0\n0 0 10 0\n0 0 10 0\n0 10 0 0\n"
    )
    out = tmp_path / "dist.txt"
    result = module.gimmemotifs_motif_score_distribution(pfm, out, metric="pcc", verbose=False)
    assert result["n_motifs"] == 4
    lengths = {(row["length_a"], row["length_b"]) for row in result["distribution"]}
    assert lengths == {(4, 4), (4, 5), (5, 4), (5, 5)}
    assert out.read_text().count("\n") == 4
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_motif_score_distribution(pfm, out, verbose=False)
    scored = module.gimmemotifs_compare_motif_pair(
        PFM_A, PFM_B, metric="pcc", pvalue=True, score_distribution=out, verbose=False
    )
    assert 0 <= scored["value"] <= 1


def test_pcc_and_ic_match_upstream():
    first, second = Motif(PFM_A), Motif(PFM_B)
    best = module.gimmemotifs_motif_pcc(PFM_A, PFM_B, verbose=False)
    expected = first.max_pcc(second)
    assert best["score"] == pytest.approx(expected[0])
    assert (best["position"], best["strand"]) == (expected[1], expected[2])
    fixed = module.gimmemotifs_motif_pcc(PFM_A, PFM_B, position=-1, verbose=False)
    assert fixed["score"] == pytest.approx(first.pcc(first.ppm, second.ppm, -1))
    no_rc = module.gimmemotifs_motif_pcc(PFM_A, PFM_B, revcomp=False, verbose=False)
    assert no_rc["strand"] == 1
    matrix = module.gimmemotifs_motif_ic_similarity(PFM_A, PFM_B, verbose=False)
    assert matrix["score"] == pytest.approx(first.matrix_ic(first.ppm, second.ppm)[0])
    assert matrix["information_content_a"] == pytest.approx(first.information_content)
    aligned = module.gimmemotifs_motif_ic_similarity(
        PFM_A, PFM_B, mode="aligned", position=1, verbose=False
    )
    assert aligned["score"] == pytest.approx(first.ic(first.ppm.tolist(), second.ppm.tolist(), 1))


def test_motif_properties_and_transforms():
    motif = Motif(PFM_A)
    span = module.gimmemotifs_motif_score_range(PFM_A, verbose=False)
    assert (span["min_score"], span["max_score"]) == (
        pytest.approx(motif.min_score),
        pytest.approx(motif.max_score),
    )
    assert (
        module.gimmemotifs_motif_consensus(PFM_A, verbose=False)["consensus"]
        == motif.to_consensusv2()
        == "ACGT"
    )
    reverse = module.gimmemotifs_motif_reverse_complement(PFM_A, verbose=False)
    assert reverse["pfm"] == [
        [7.0, 1.0, 1.0, 1.0],
        [2.0, 8.0, 0.0, 0.0],
        [0.0, 0.0, 10.0, 0.0],
        [0.0, 0.0, 0.0, 10.0],
    ]
    first = module.gimmemotifs_motif_shuffle(PFM_A, seed=3, verbose=False)
    again = module.gimmemotifs_motif_shuffle(PFM_A, seed=3, verbose=False)
    assert first["pfm"] == again["pfm"]
    assert sorted(map(tuple, first["pfm"])) == sorted(
        map(tuple, [[float(v) for v in r] for r in PFM_A])
    )
    ppm = module.gimmemotifs_pfm_to_ppm(PFM_A, verbose=False)["ppm"]
    assert ppm == [list(map(float, r)) for r in Motif().pfm_to_ppm(PFM_A)]
    assert all(sum(row) == pytest.approx(1.0) for row in ppm)
    pseudo = module.gimmemotifs_ppm_pseudocount([[1, 0, 0, 0]], verbose=False)["ppm"][0]
    assert pseudo[0] == pytest.approx((1 + 1e-6) / (1 + 4e-6))
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_pfm_to_ppm([[1, 2, 3]], verbose=False)


def test_format_motif_matches_upstream():
    motif = Motif(PFM_A)
    motif.id = "m1"
    for fmt, expected in [
        ("ppm", motif.to_ppm()),
        ("pfm", motif.to_pfm()),
        ("meme", motif.to_meme()),
        ("motevo", motif.to_motevo()),
        ("transfac", motif.to_transfac()),
    ]:
        text = module.gimmemotifs_format_motif(
            PFM_A, output_format=fmt, motif_id="m1", verbose=False
        )["text"]
        assert text == expected
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_format_motif(PFM_A, output_format="jaspar", verbose=False)


def test_average_motifs_and_sampling():
    result = module.gimmemotifs_average_motifs(
        PFM_A, PFM_B, position=-1, orientation=-1, verbose=False
    )
    expected = Motif(PFM_A).average_motifs(Motif(PFM_B), -1, -1)
    assert result["pfm"] == [[float(v) for v in row] for row in expected.pfm]
    assert result["consensus"] == expected.consensus
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_average_motifs(PFM_A, PFM_B, orientation=2, verbose=False)
    sampled = module.gimmemotifs_sample_motif_sequences(PFM_A, 5, seed=4, verbose=False)
    assert sampled["sequences"] == Motif(PFM_A).sample(5, rng=np.random.default_rng(4))
    assert all(len(s) == 4 and s[:2] == "AC" for s in sampled["sequences"])


def test_pwm_scan_outputs_match_upstream(tmp_path):
    fasta = _fasta(tmp_path)
    motif = Motif([[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]])
    fa = Fasta(str(fasta))
    scores = module.gimmemotifs_pwm_scan(fasta, motif.pfm.tolist(), verbose=False)
    assert scores["matches"] == {
        k: list(map(float, v)) for k, v in motif.pwm_scan_score(fa).items()
    }
    hits = module.gimmemotifs_pwm_scan(
        fasta, motif.pfm.tolist(), output_format="hits", cutoff=0.95, verbose=False
    )
    upstream = motif.scan_all(fa, cutoff=0.95)
    assert {
        k: [(h["position"], h["score"], h["strand"]) for h in v] for k, v in hits["matches"].items()
    } == {k: [(p, float(s), t) for p, s, t in v] for k, v in upstream.items()}
    positions = module.gimmemotifs_pwm_scan(
        fasta, motif.pfm.tolist(), output_format="positions", cutoff=0.95, verbose=False
    )
    assert positions["matches"]["s1"] == [2, 8, 2, 8]
    gff = tmp_path / "hits.gff"
    written = module.gimmemotifs_pwm_scan(
        fasta,
        motif.pfm.tolist(),
        output_format="gff",
        gff_path=gff,
        cutoff=0.95,
        motif_id="acgt",
        verbose=False,
    )
    lines = gff.read_text().splitlines()
    assert written["outputs"]["gff"] == str(gff)
    assert lines[0].split("\t")[:5] == ["s1", "pfmscan", "misc_feature", "2", "6"]
    assert 'motif_name "acgt"' in lines[0]
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_pwm_scan(fasta, motif.pfm.tolist(), output_format="gff", verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_pwm_scan(
            fasta, motif.pfm.tolist(), output_format="gff", gff_path=gff, verbose=False
        )


def test_consensus_scan_matches_upstream(tmp_path):
    fasta = _fasta(tmp_path)
    result = module.gimmemotifs_consensus_scan(fasta, PFM_A, verbose=False)
    expected = Motif(PFM_A).consensus_scan(Fasta(str(fasta)))
    assert result["matches"] == {k: list(map(float, v)) for k, v in expected.items()}
    assert result["consensus"] == "ACGT"
    assert result["matches"]["s1"] == [4.0, 10.0]
