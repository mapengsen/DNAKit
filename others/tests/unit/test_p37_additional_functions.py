"""Offline Chamaeleo, seqme and contig-reader checks with real optional backends."""

import itertools

import pytest


@pytest.mark.parametrize("bits", [[[0, 0]], [[0, 0, 0, 1, 1, 0, 1, 1]], [[1, 1], [0, 1]]])
def test_chamaeleo_direct_matches_upstream(bits):
    pytest.importorskip("Chamaeleo")
    from Chamaeleo.methods.default import BaseCodingAlgorithm

    from dnakit.molecular_lab import chamaeleo_p37_functions as m

    expected = BaseCodingAlgorithm().silicon_to_carbon(bits, bit_size=sum(map(len, bits)))
    got = m.chamaeleo_encode_direct(bits, verbose=False)
    assert got["dna"] == ["".join(row) for row in expected["dna"]]
    assert m.chamaeleo_decode_direct(got["dna"], verbose=False)["bit_segments"] == bits


@pytest.mark.parametrize("seq", ["ACGT", "AAAA", "A"])
def test_chamaeleo_screen_matches_upstream(seq):
    pytest.importorskip("Chamaeleo")
    from Chamaeleo.utils import screen

    from dnakit.molecular_lab import chamaeleo_p37_functions as m

    got = m.chamaeleo_screen_dna(seq, max_homopolymer=3, max_gc_deviation=0.1, verbose=False)
    assert got["passes"] == bool(screen.check(seq, 3, 0.6))


@pytest.mark.parametrize("sequences", [["ACGT", "GGGG", "AAAA"], ["A", "G"], ["ACGT", "ACGT"]])
def test_seqme_property_mean_matches_upstream(sequences):
    np = pytest.importorskip("numpy")
    pytest.importorskip("seqme")
    from seqme.metrics import ID

    from dnakit.predictions import seqme_p37_functions as m

    def predictor(seqs):
        return np.asarray([(s.count("G") + s.count("C")) / len(s) for s in seqs])

    direct = ID(predictor, "gc_content", "maximize", deviation="se", estimate="unbiased")(sequences)
    got = m.seqme_property_mean(sequences, "gc_content", verbose=False)
    assert got["value"] == pytest.approx(direct.value)
    assert got["deviation"] == pytest.approx(direct.deviation)


@pytest.mark.parametrize("sequences", [["ACGT", "GGGG", "AAAA"], ["A", "G"], ["ACGT", "ACGT"]])
def test_seqme_mmd_matches_upstream(sequences):
    pytest.importorskip("seqme")
    from seqme.metrics import MMD
    from seqme.models import KmerFrequencyEmbedding

    from dnakit.predictions import seqme_p37_functions as m

    embedder = KmerFrequencyEmbedding(["".join(k) for k in itertools.product("ACGT", repeat=1)])
    reference = ["AT", "GC", "GT"]
    expected = MMD(reference, embedder, estimate="biased", sigma=10, scale=1000)(sequences)
    assert m.seqme_mmd(sequences, reference, kmer_size=1, verbose=False)["value"] == pytest.approx(
        expected.value
    )


@pytest.mark.parametrize("sequences", [["ACGT"], ["A", "TTT"], ["N"]])
def test_contig_reader_matches_biopython(tmp_path, sequences):
    pytest.importorskip("Bio")
    from Bio import SeqIO

    from dnakit.visualization import pygenomeviz_p37_functions as m

    path = tmp_path / "seq.fa"
    path.write_text("".join(f">seq{i}\n{s}\n" for i, s in enumerate(sequences)))
    expected = list(SeqIO.parse(path, "fasta"))
    got = m.read_contig_sequences(path, verbose=False)
    assert [r["sequence"] for r in got["contigs"]] == [str(r.seq) for r in expected]
