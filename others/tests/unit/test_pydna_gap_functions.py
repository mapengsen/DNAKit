"""Tests for dnakit.molecular_lab.pydna_gap_functions (codon index, SnapGene, gel, BLAST, Tm).

Offline: every input is written into tmp_path or built in memory. The SnapGene fixture is a .dna
file written with sgffp (pydna's own SnapGene dependency). The BLAST and NEB wrappers are tested
with fakes that replace the pydna attribute they call, so nothing reaches the network.
"""

from __future__ import annotations

import importlib.util
import inspect
import socket
import sys
from pathlib import Path

import pytest

pytest.importorskip("pydna")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import pydna_gap_functions as gap

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
CODING = "ATGGCTAAAGCTTAA"  # ATG, GCT, AAA, GCT, TAA: ATG and the stop are not scored


class _NoConnectSocket(socket.socket):
    def connect(self, *args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    def connect_ex(self, *args: object, **kwargs: object) -> int:
        raise AssertionError("network access is not allowed in this test")


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", _NoConnectSocket)


def write_snapgene_file(path: Path, sequence: str, *, circular: bool) -> Path:
    """Write a SnapGene .dna file with one sequence block and no cloning history."""

    sgffp = pytest.importorskip("sgffp")
    topology = "circular" if circular else "linear"
    obj = sgffp.SgffObject.new(sequence=sequence, topology=topology)
    path.write_bytes(sgffp.SgffWriter.to_bytes(obj))
    return path


def test_public_functions_take_keyword_only_verbose() -> None:
    assert gap.__all__
    for name in gap.__all__:
        parameter = inspect.signature(getattr(gap, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


# Codon adaptation index.


def test_cai_with_explicit_weights_is_the_geometric_mean_of_scored_codons() -> None:
    weights = {"GCT": 0.5, "AAA": 0.25}
    result = gap.pydna_codon_adaptation_index(CODING, weights=weights, verbose=False)
    # Scored codons: GCT, AAA, GCT. (0.5 * 0.25 * 0.5) ** (1/3) = 0.397
    assert result["cai"] == 0.397
    assert result["codon_count"] == 5
    assert result["scored_codon_count"] == 3


def test_cai_with_explicit_weights_rejects_missing_codons_and_bad_sequences() -> None:
    with pytest.raises(ConfigurationError, match=r"lack"):
        gap.pydna_codon_adaptation_index(CODING, weights={"GCT": 0.5}, verbose=False)
    with pytest.raises(ConfigurationError, match=r"positive"):
        gap.pydna_codon_adaptation_index(CODING, weights={"GCT": 0.0, "AAA": 1.0}, verbose=False)
    with pytest.raises(ConfigurationError, match=r"whole number of codons"):
        gap.pydna_codon_adaptation_index("ATGGC", weights={"GCT": 1.0}, verbose=False)
    with pytest.raises(ConfigurationError, match=r"only A, C, G and T"):
        gap.pydna_codon_adaptation_index("ATGNNN", weights={"NNN": 1.0}, verbose=False)
    with pytest.raises(ConfigurationError, match=r"no codon left"):
        gap.pydna_codon_adaptation_index("ATGTGGTAA", weights={"GCT": 1.0}, verbose=False)


def test_cai_organism_route_needs_cai2(monkeypatch: pytest.MonkeyPatch) -> None:
    real_find_spec = importlib.util.find_spec

    def without_cai2(name: str, *args: object, **kwargs: object) -> object:
        return None if name == "cai2" else real_find_spec(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", without_cai2)
    with pytest.raises(BackendUnavailableError, match=r"cai2"):
        gap.pydna_codon_adaptation_index(CODING, verbose=False)


# SnapGene history files.


def test_parse_snapgene_history_reads_the_root_sequence(tmp_path: Path) -> None:
    sequence = "ATGCGTACGTTAGCATGCGTACGTTAGCATGCGTACGTTAGC"
    path = write_snapgene_file(tmp_path / "construct.dna", sequence, circular=True)
    result = gap.pydna_parse_snapgene_history(str(path), verbose=False)
    assert result["name"] == "construct"
    assert result["length"] == len(sequence)
    assert result["sequence"] == sequence
    assert result["circular"] is True
    assert result["has_history"] is False
    assert result["feature_count"] == 0


def test_parse_snapgene_history_guards(tmp_path: Path) -> None:
    other = tmp_path / "x.gb"
    other.write_text("LOCUS\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"ending in \.dna"):
        gap.pydna_parse_snapgene_history(str(other), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.pydna_parse_snapgene_history(str(tmp_path / "missing.dna"), verbose=False)


# Gel images and the size model.


def test_predict_relative_front_passes_through_the_ladder_knots() -> None:
    result = gap.pydna_predict_gel_relative_front([10000, 5000, 30000], verbose=False)
    assert result["ladder"] == "GeneRuler_1kb_plus"
    assert result["relative_front"][0] == pytest.approx(0.04, abs=1e-9)
    assert result["relative_front"][1] == pytest.approx(0.154, abs=1e-9)
    assert result["relative_front"][2] is None  # above the largest ladder band
    assert result["ladder_sizes"][-1] == 20000


def test_predict_relative_front_rejects_unknown_ladder() -> None:
    with pytest.raises(ConfigurationError, match=r"ladder"):
        gap.pydna_predict_gel_relative_front([1000], ladder="NoSuchLadder", verbose=False)


def test_draw_gel_writes_a_png_with_one_column_per_lane(tmp_path: Path) -> None:
    target = tmp_path / "gel" / "lanes.png"
    result = gap.pydna_draw_gel(
        str(target), lanes=[[1000, 500], [3000]], gel_length=300, verbose=False
    )
    assert result["path"] == str(target.resolve())
    assert target.read_bytes()[:8] == PNG_SIGNATURE
    from matplotlib.image import imread

    pixels = imread(str(target))[..., :3]
    assert (pixels.max(axis=-1) > 0.5).mean() > 0.01  # bright bands on the black gel
    data = result["data"]
    assert data["lane_count"] == 2 and data["band_count"] == 3
    assert data["image_width"] == 60 + 60 * 2
    assert data["image_height"] == 300


def test_draw_gel_default_draws_the_ladder(tmp_path: Path) -> None:
    result = gap.pydna_draw_gel(str(tmp_path / "ladder.png"), verbose=False)
    assert result["data"]["lane_count"] == 1
    assert result["data"]["lanes"][0][0] == 20000


def test_draw_gel_guards(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"outside"):
        gap.pydna_draw_gel(str(tmp_path / "a.png"), lanes=[[5]], verbose=False)
    with pytest.raises(ConfigurationError, match=r"lanes"):
        gap.pydna_draw_gel(str(tmp_path / "b.png"), lanes=[[1000]] * 21, verbose=False)
    existing = tmp_path / "c.png"
    existing.write_bytes(PNG_SIGNATURE)
    with pytest.raises(ConfigurationError, match=r"already exists"):
        gap.pydna_draw_gel(str(existing), verbose=False)
    with pytest.raises(ConfigurationError, match=r"\.png"):
        gap.pydna_draw_gel(str(tmp_path / "d.pdf"), verbose=False)


def test_missing_pydna_is_a_backend_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pydna.ladders", None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        gap.pydna_predict_gel_relative_front([1000], verbose=False)
    assert "bio-assembly" in excinfo.value.hint


# Network wrappers, tested with fakes.


def test_blast_best_nt_hit_calls_the_pydna_attribute(monkeypatch: pytest.MonkeyPatch) -> None:
    import pydna.sequence_picker as picker
    from pydna.dseqrecord import Dseqrecord

    calls: list[str] = []

    def fake_genbank_accession(sequence: str) -> Dseqrecord:
        calls.append(sequence)
        return Dseqrecord("ATGCCC", circular=False, description="FAKE000001.1 REGION: 1..6")

    monkeypatch.setattr(picker, "genbank_accession", fake_genbank_accession)
    result = gap.pydna_blast_best_nt_hit("atgcccaaa", name="probe", verbose=False)
    assert calls == ["ATGCCCAAA"]
    assert result["service"].startswith("NCBI BLAST")
    assert result["result"]["sequence"] == "ATGCCC"
    assert result["result"]["query_name"] == "probe"
    assert result["result"]["description"] == "FAKE000001.1 REGION: 1..6"


def test_blast_best_nt_hit_guards() -> None:
    with pytest.raises(ConfigurationError, match=r"only A, C, G and T"):
        gap.pydna_blast_best_nt_hit("ACGTX", verbose=False)
    with pytest.raises(ConfigurationError, match=r"limit"):
        gap.pydna_blast_best_nt_hit("A" * 20_001, verbose=False)


def test_neb_melting_temperature_calls_the_pydna_attribute(monkeypatch: pytest.MonkeyPatch) -> None:
    import pydna.tm as tm

    calls: list[tuple[str, float, str]] = []

    def fake_tm_neb(primer: str, conc: float = 0.5, prodcode: str = "q5-0") -> float:
        calls.append((primer, conc, prodcode))
        return 58.4

    monkeypatch.setattr(tm, "tm_neb", fake_tm_neb)
    result = gap.pydna_neb_melting_temperature("gattacagg", conc=0.25, verbose=False)
    assert calls == [("GATTACAGG", 0.25, "q5-0")]
    assert result["source"] == "https://tmapi.neb.com/tm"
    assert result["result"]["tm_celsius"] == 58.4
    assert result["result"]["conc_uM"] == 0.25


def test_neb_melting_temperature_guards() -> None:
    with pytest.raises(ConfigurationError, match=r"primer"):
        gap.pydna_neb_melting_temperature("GATTACAXX", verbose=False)
    with pytest.raises(ConfigurationError, match=r"conc"):
        gap.pydna_neb_melting_temperature("GATTACA", conc=0.0, verbose=False)
    with pytest.raises(ConfigurationError, match=r"prodcode"):
        gap.pydna_neb_melting_temperature("GATTACA", prodcode="q5 0", verbose=False)
