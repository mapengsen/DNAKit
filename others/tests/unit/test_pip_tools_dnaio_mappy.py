"""Contract tests for the optional dnaio and mappy adapters."""

from __future__ import annotations

from types import SimpleNamespace


def test_dnaio_summary_is_bounded_and_json_serializable(monkeypatch, tmp_path) -> None:
    import dnakit.integrations.pip_tools as pip_tools

    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")

    class FakeReader:
        def __iter__(self):
            return iter(
                [
                    SimpleNamespace(name="r1", sequence="ACGT", qualities="IIII"),
                    SimpleNamespace(name="r2", sequence="GGCA", qualities=None),
                ]
            )

        def close(self):
            return None

    class FakeDnaio:
        @staticmethod
        def open(path):
            assert path == str(source.resolve())
            return FakeReader()

    def fake_require(module, extra):
        assert (module, extra) == ("dnaio", "bio-dnaio")
        return FakeDnaio

    monkeypatch.setattr(pip_tools, "require", fake_require)
    result = pip_tools.summarize_dnaio_fastx(source, preview_records=2, verbose=False)

    assert result["tool"] == "dnaio"
    assert result["records"] == 2
    assert result["fastq_records"] == 1
    assert result["total_bases"] == 8
    assert result["gc_bases"] == 5
    assert result["preview"] == [{"name": "r1", "length": 4}, {"name": "r2", "length": 4}]
    assert result["provenance"]["license"] == "MIT"


def test_mappy_mapping_serializes_allowlisted_hits(monkeypatch, tmp_path) -> None:
    import dnakit.integrations.pip_tools as pip_tools

    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fa"
    reference.write_text(">chr1\nAACCGGTT\n", encoding="ascii")
    reads.write_text(">r1\nACCG\n", encoding="ascii")

    class FakeHit:
        ctg = "chr1"
        r_st = 1
        r_en = 5
        q_st = 0
        q_en = 4
        strand = 1
        mapq = 60
        mlen = 4
        blen = 4
        cigar_str = "4M"
        cs = ":4"
        MD = "4"

    class FakeAligner:
        def __init__(self, path, *, preset, n_threads):
            assert path == str(reference.resolve())
            assert preset == "map-ont"
            assert n_threads == 1

        def map(self, sequence, *, cs, MD):
            assert sequence == "ACCG"
            assert cs is True
            assert MD is True
            return iter([FakeHit()])

    class FakeMappy:
        Aligner = FakeAligner

        @staticmethod
        def fastx_read(path):
            assert path == str(reads.resolve())
            return iter([("r1", "ACCG", "")])

    def fake_require(module, extra):
        assert (module, extra) == ("mappy", "bio-mappy")
        return FakeMappy

    monkeypatch.setattr(pip_tools, "require", fake_require)
    result = pip_tools.map_mappy_sequences(
        reference,
        reads,
        preset="map-ont",
        max_reads=1,
        max_hits_per_read=1,
        cs=True,
        md=True,
        verbose=False,
    )

    assert result["tool"] == "mappy"
    assert result["reads_mapped"] == 1
    assert result["hit_count"] == 1
    assert result["hits"][0]["contig"] == "chr1"
    assert result["hits"][0]["reference_start"] == 1
    assert result["hits"][0]["cs"] == ":4"
    assert result["hits"][0]["md"] == "4"
    assert result["provenance"]["citation"] == "10.1093/bioinformatics/bty191"
