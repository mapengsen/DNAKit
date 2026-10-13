"""Contract tests for optional gimmemotifs, tangermeme and snapatac2 regulation adapters."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.regulation.backends import (
    count_peak_fragments,
    count_tile_matrix,
    deep_lift_shap_attribution,
    marginalize_motif_effect,
    predict_sequences,
    read_fragment_qc,
    read_motif_file,
    scan_motif_hits,
    score_motif_kmers,
)

# PFM counts for the consensus GTTT (columns A, C, G, T); rows are positions.
GTTT_PFM = [[0, 0, 10, 0], [0, 0, 0, 10], [0, 0, 0, 10], [0, 0, 0, 10]]
PFM_FILE = ">MOT1\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
TRANSFAC_FILE = (
    "ID  MOT2\nP0 A C G T\n01 10 0 0 0 A\n02 0 10 0 0 C\n03 0 0 10 0 G\n04 0 0 0 10 T\nXX\n//\n"
)
JASPAR_FILE = ">MOT3\nA [ 10 0 0 0 ]\nC [ 0 10 0 0 ]\nG [ 0 0 10 0 ]\nT [ 0 0 0 10 ]\n"

CHROM_SIZES = {"chr1": 2000, "chr2": 2000, "chrM": 1000}
# Five columns: chrom, start, end, barcode, read count. cellA has a 2-read fragment and one
# chrM fragment; cellB has a 3-read fragment on chr2.
FRAGMENTS = (
    "chr1\t100\t300\tcellA\t1\n"
    "chr1\t120\t330\tcellA\t2\n"
    "chr1\t600\t800\tcellA\t1\n"
    "chrM\t10\t200\tcellA\t1\n"
    "chr1\t650\t900\tcellB\t1\n"
    "chr2\t50\t250\tcellB\t1\n"
    "chr2\t1200\t1400\tcellB\t3\n"
)
PEAKS = "chr1\t50\t400\nchr1\t590\t950\nchr2\t1100\t1500\n"


def _flat(value: object) -> list[float]:
    """Flatten nested lists so pytest.approx can compare tensors element by element."""

    if isinstance(value, list):
        return [number for item in value for number in _flat(item)]
    return [float(value)]  # type: ignore[arg-type]


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _hand_log_odds(ppm: list[list[float]], kmer: str) -> float:
    """Sum of ln(ppm / 0.25 + 0.01) over positions, the gimmemotifs log-odds rule."""

    bases = "ACGT"
    return sum(
        math.log(ppm[index][bases.index(base)] / 0.25 + 0.01) for index, base in enumerate(kmer)
    )


def test_read_motif_file_reads_each_format_with_the_same_counts(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    from gimmemotifs.motif import read_motifs

    for name, fmt, text in (
        ("motifs.pfm", "pfm", PFM_FILE),
        ("motifs.tf", "transfac", TRANSFAC_FILE),
        ("motifs.jaspar", "jaspar", JASPAR_FILE),
    ):
        path = _write(tmp_path, name, text)
        result = read_motif_file(path, fmt=fmt)
        direct = read_motifs(str(path), fmt=fmt, as_dict=False)[0]
        assert result["motif_count"] == 1
        motif = result["motifs"][0]
        assert motif["consensus"] == "ACGT"
        assert motif["length"] == 4
        assert motif["pfm"] == [
            [10.0, 0.0, 0.0, 0.0],
            [0.0, 10.0, 0.0, 0.0],
            [0.0, 0.0, 10.0, 0.0],
            [0.0, 0.0, 0.0, 10.0],
        ]
        assert motif["pfm"] == [[float(v) for v in row] for row in direct.pfm]
        assert motif["max_score"] == pytest.approx(float(direct.max_score))
    assert result["operation"] == "read_motif_file"
    assert result["provenance"]["package"] == "gimmemotifs"


def test_read_motif_file_rejects_bad_format_and_large_files(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    path = _write(tmp_path, "motifs.pfm", PFM_FILE)
    with pytest.raises(ConfigurationError, match="fmt must be one of"):
        read_motif_file(path, fmt="meme")
    with pytest.raises(ConfigurationError, match="max_bytes"):
        read_motif_file(path, max_bytes=10)
    with pytest.raises(ConfigurationError, match="max_motifs"):
        read_motif_file(path, max_motifs=0)


def test_scan_motif_hits_matches_hand_positions_and_direct_upstream(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    from gimmemotifs.motif import Motif

    sequences = {"s": "GGAAACGG", "t": "GTTTAAAA"}
    result = scan_motif_hits(sequences, {"GTTT": GTTT_PFM}, threshold=0.9)
    # GTTT occurs on the reverse strand of "GGAAACGG" (AAAC at 0-based 2..5) and forward in "t".
    assert [
        (hit["sequence"], hit["start"], hit["end"], hit["strand"]) for hit in result["hits"]
    ] == [
        ("s", 2, 6, "-"),
        ("t", 0, 4, "+"),
    ]
    motif = Motif(GTTT_PFM)
    expected = _hand_log_odds(motif.ppm.tolist(), "GTTT")
    assert result["hits"][0]["score"] == pytest.approx(expected)
    assert result["hits"][0]["score"] == pytest.approx(float(motif.score_kmer("GTTT")))

    forward_only = scan_motif_hits(sequences, {"GTTT": GTTT_PFM}, scan_rc=False)
    assert [hit["sequence"] for hit in forward_only["hits"]] == ["t"]
    assert forward_only["hit_count"] == 1


def test_scan_motif_hits_threshold_and_hit_limit(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    sequence = {"repeat": "GTTT" * 5}
    result = scan_motif_hits(sequence, {"GTTT": GTTT_PFM}, threshold=1.0, scan_rc=False)
    assert [hit["start"] for hit in result["hits"]] == [0, 4, 8, 12, 16]
    with pytest.raises(ConfigurationError, match="max_hits_per_sequence=2"):
        scan_motif_hits(sequence, {"GTTT": GTTT_PFM}, max_hits_per_sequence=2, scan_rc=False)


def test_scan_motif_hits_rejects_bad_input() -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    motifs = {"GTTT": GTTT_PFM}
    with pytest.raises(ConfigurationError, match="non-ACGT"):
        scan_motif_hits({"s": "GTNTT"}, motifs)
    with pytest.raises(ConfigurationError, match="At least one sequence"):
        scan_motif_hits({}, motifs)
    with pytest.raises(ConfigurationError, match="max_sequences"):
        scan_motif_hits({"a": "ACGT", "b": "ACGT"}, motifs, max_sequences=1)
    with pytest.raises(ConfigurationError, match="threshold"):
        scan_motif_hits({"s": "ACGT"}, motifs, threshold=True)
    with pytest.raises(ConfigurationError, match="threshold"):
        scan_motif_hits({"s": "ACGT"}, motifs, threshold=1.5)
    with pytest.raises(ConfigurationError, match="scan_rc"):
        scan_motif_hits({"s": "ACGT"}, motifs, scan_rc="yes")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="motifs must be"):
        scan_motif_hits({"s": "ACGT"}, {})
    with pytest.raises(ConfigurationError, match="non-negative"):
        scan_motif_hits({"s": "ACGT"}, {"bad": [[-1, 0, 0, 1]]})


def test_score_motif_kmers_matches_formula_and_upstream() -> None:
    pytest.importorskip("gimmemotifs.motif", exc_type=ImportError)
    from gimmemotifs.motif import Motif

    result = score_motif_kmers(GTTT_PFM, ["GTTT", "AAAA"])
    ppm = Motif(GTTT_PFM).ppm.tolist()
    scores = {item["kmer"]: item["score"] for item in result["scores"]}
    assert scores["GTTT"] == pytest.approx(_hand_log_odds(ppm, "GTTT"))
    assert scores["AAAA"] == pytest.approx(_hand_log_odds(ppm, "AAAA"))
    assert scores["AAAA"] == pytest.approx(float(Motif(GTTT_PFM).score_kmer("AAAA")))
    best = sum(math.log(max(row) / 0.25 + 0.01) for row in ppm)
    assert result["max_score"] == pytest.approx(best)
    assert result["consensus"] == "GTTT"
    with pytest.raises(ConfigurationError, match="length 4"):
        score_motif_kmers(GTTT_PFM, ["GTT"])
    with pytest.raises(ConfigurationError, match="non-ACGT"):
        score_motif_kmers(GTTT_PFM, ["GTTN"])
    with pytest.raises(ConfigurationError, match="max_kmers"):
        score_motif_kmers(GTTT_PFM, ["GTTT", "GTTT"], max_kmers=1)


def test_missing_gimmemotifs_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, "motifs.pfm", PFM_FILE)
    monkeypatch.setitem(sys.modules, "gimmemotifs", None)
    monkeypatch.setitem(sys.modules, "gimmemotifs.motif", None)
    with pytest.raises(BackendUnavailableError, match="bio-gimmemotifs"):
        read_motif_file(path)
    with pytest.raises(BackendUnavailableError, match="bio-gimmemotifs"):
        scan_motif_hits({"s": "ACGT"}, {"GTTT": GTTT_PFM})


def test_read_fragment_qc_matches_hand_counts_and_upstream(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    import snapatac2 as snap

    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    result = read_fragment_qc(fragments, CHROM_SIZES, min_num_fragments=1)
    cells = {cell["barcode"]: cell for cell in result["cells"]}
    # cellA: 3 non-chrM fragments; reads 1+2+1+1 with one duplicate read; chrM is 1 of 4 fragments.
    assert cells["cellA"]["n_fragment"] == 3
    assert cells["cellA"]["frac_dup"] == pytest.approx(0.2)
    assert cells["cellA"]["frac_mito"] == pytest.approx(0.25)
    # cellB: 3 fragments; reads 1+1+3 with two duplicate reads; no chrM fragment.
    assert cells["cellB"]["n_fragment"] == 3
    assert cells["cellB"]["frac_dup"] == pytest.approx(0.4)
    assert cells["cellB"]["frac_mito"] == pytest.approx(0.0)

    direct = snap.pp.import_fragments(
        fragments,
        CHROM_SIZES,
        is_paired=True,
        min_num_fragments=1,
        file=tmp_path / "direct.h5ad",
        sorted_by_barcode=False,
        tempdir=tmp_path,
        n_jobs=1,
    )
    try:
        direct_counts = dict(
            zip(direct.obs_names, direct.obs["n_fragment"][:].to_list(), strict=True)
        )
    finally:
        direct.close()
    assert {name: cells[name]["n_fragment"] for name in cells} == direct_counts


def test_count_tile_matrix_matches_hand_counts_and_upstream(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    import snapatac2 as snap

    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    result = count_tile_matrix(fragments, CHROM_SIZES, bin_size=500, min_num_fragments=1)
    assert result["shape"] == [2, 8]
    assert result["features"][:2] == ["chr1:0-500", "chr1:500-1000"]
    assert result["features"][4] == "chr2:0-500"
    # cellA: two fragments in chr1:0-500 and one in chr1:500-1000; chrM is excluded.
    # cellB: chr1:500-1000, chr2:0-500 and chr2:1000-1500, one fragment each.
    assert result["entries"] == [[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 4, 1], [1, 6, 1]]
    assert result["total_count"] == 6

    direct = snap.pp.import_fragments(
        fragments,
        CHROM_SIZES,
        is_paired=True,
        min_num_fragments=1,
        file=tmp_path / "tile.h5ad",
        sorted_by_barcode=False,
        tempdir=tmp_path,
        n_jobs=1,
    )
    try:
        snap.pp.add_tile_matrix(direct, bin_size=500, exclude_chroms=["chrM", "chrY", "M", "Y"])
        dense = direct.X[:].toarray()
    finally:
        direct.close()
    rebuilt = [[0] * 8 for _ in range(2)]
    for row, column, value in result["entries"]:
        rebuilt[row][column] = value
    assert rebuilt == dense.tolist()


def test_count_tile_matrix_keeps_chrM_when_not_excluded(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    result = count_tile_matrix(
        fragments, CHROM_SIZES, bin_size=500, exclude_chroms=(), min_num_fragments=1
    )
    assert result["shape"] == [2, 10]
    assert [0, 8, 1] in result["entries"]  # chrM:0-500 receives the cellA mitochondrial fragment


def test_count_tile_matrix_rejects_bad_options(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    with pytest.raises(ConfigurationError, match="bin_size"):
        count_tile_matrix(fragments, CHROM_SIZES, bin_size=True)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="bin_size"):
        count_tile_matrix(fragments, CHROM_SIZES, bin_size=0)
    with pytest.raises(ConfigurationError, match="counting_strategy"):
        count_tile_matrix(fragments, CHROM_SIZES, counting_strategy="reads")
    with pytest.raises(ConfigurationError, match="max_entries"):
        count_tile_matrix(fragments, CHROM_SIZES, min_num_fragments=1, max_entries=2)


def test_count_peak_fragments_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    peaks = _write(tmp_path, "peaks.bed", PEAKS)
    result = count_peak_fragments(fragments, CHROM_SIZES, peaks, min_num_fragments=1)
    assert result["features"] == ["chr1:50-400", "chr1:590-950", "chr2:1100-1500"]
    assert result["barcodes"] == ["cellA", "cellB"]
    # cellA: two fragments in the first peak and one in the second; cellB: one in each of the
    # second and third peaks (the 3-read fragment still counts once under the fragment strategy).
    assert result["entries"] == [[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 2, 1]]
    assert result["total_count"] == 5


def test_fragment_and_peak_inputs_are_validated(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2", exc_type=ImportError)
    cases = {
        "columns": "chr1\t100\t300\tcellA\n",
        "unknown_chrom": "chr9\t100\t300\tcellA\t1\n",
        "span": "chr1\t300\t100\tcellA\t1\n",
        "count": "chr1\t100\t300\tcellA\t0\n",
        "barcode": "chr1\t100\t300\t\t1\n",
        "comment": "# header\nchr1\t100\t300\tcellA\t1\n",
    }
    for name, text in cases.items():
        path = _write(tmp_path, f"{name}.tsv", text)
        with pytest.raises(ConfigurationError):
            read_fragment_qc(path, CHROM_SIZES, min_num_fragments=1)
    with pytest.raises(ConfigurationError, match="does not exist"):
        read_fragment_qc(tmp_path / "missing.tsv", CHROM_SIZES)
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    with pytest.raises(ConfigurationError, match="max_fragments"):
        read_fragment_qc(fragments, CHROM_SIZES, max_fragments=2)
    with pytest.raises(ConfigurationError, match="max_cells"):
        read_fragment_qc(fragments, CHROM_SIZES, min_num_fragments=1, max_cells=1)
    bed_two_columns = _write(tmp_path, "peaks2.bed", "chr1\t50\n")
    with pytest.raises(ConfigurationError, match="needs chrom, start and end"):
        count_peak_fragments(fragments, CHROM_SIZES, bed_two_columns)
    peaks = _write(tmp_path, "peaks.bed", PEAKS)
    with pytest.raises(ConfigurationError, match="max_peaks"):
        count_peak_fragments(fragments, CHROM_SIZES, peaks, max_peaks=2)


def test_upstream_base_exception_becomes_backend_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    snap = pytest.importorskip("snapatac2", exc_type=ImportError)

    class FakePanic(BaseException):
        """Stands in for a pyo3 panic, which derives from BaseException."""

    def boom(*_args: object, **_kwargs: object) -> None:
        raise FakePanic("simulated upstream panic")

    monkeypatch.setattr(snap.pp, "import_fragments", boom)
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    with pytest.raises(BackendExecutionError, match="simulated upstream panic"):
        read_fragment_qc(fragments, CHROM_SIZES)


def test_missing_snapatac2_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fragments = _write(tmp_path, "fragments.tsv", FRAGMENTS)
    monkeypatch.setitem(sys.modules, "snapatac2", None)
    with pytest.raises(BackendUnavailableError, match="bio-snapatac2"):
        read_fragment_qc(fragments, CHROM_SIZES)


def _count_a_model(torch: object) -> object:
    """Tiny model whose output is the number of A bases (one-hot channel 0) per sequence."""

    class CountA(torch.nn.Module):  # type: ignore[attr-defined, misc]
        def forward(self, inputs: object) -> object:
            return inputs[:, 0, :].sum(dim=1, keepdim=True)  # type: ignore[index]

    return CountA().eval()


def test_predict_sequences_matches_direct_model_call() -> None:
    torch = pytest.importorskip("torch", exc_type=ImportError)
    pytest.importorskip("tangermeme.predict", exc_type=ImportError)
    torch.manual_seed(0)
    model = torch.nn.Sequential(torch.nn.Flatten(), torch.nn.Linear(4 * 6, 2)).eval()
    sequences = {"a": "ACGTAC", "b": "TTTTGG"}
    result = predict_sequences(model, sequences)
    from tangermeme.utils import one_hot_encode

    inputs = torch.stack([one_hot_encode(text) for text in sequences.values()]).float()
    with torch.no_grad():
        direct = model(inputs).tolist()
    assert result["shape"] == [2, 2]
    assert _flat(result["predictions"]) == pytest.approx(_flat(direct))


def test_marginalize_motif_effect_matches_hand_counts() -> None:
    torch = pytest.importorskip("torch", exc_type=ImportError)
    pytest.importorskip("tangermeme.marginalize", exc_type=ImportError)
    model = _count_a_model(torch)
    # Ten C bases; AAAA is written at the centre, positions 3..6 (10 // 2 - 4 // 2 = 3).
    result = marginalize_motif_effect(model, {"c": "CCCCCCCCCC", "ac": "ACCCCCCCCC"}, "AAAA")
    assert result["start"] == 3
    assert result["y_before"] == [[0.0], [1.0]]
    assert result["y_after"] == [[4.0], [5.0]]
    assert result["delta"] == [[4.0], [4.0]]
    explicit = marginalize_motif_effect(model, {"c": "CCCCCCCCCC"}, "AAAA", start=6)
    assert explicit["start"] == 6
    assert explicit["delta"] == [[4.0]]


def test_deep_lift_shap_attribution_matches_direct_call() -> None:
    torch = pytest.importorskip("torch", exc_type=ImportError)
    pytest.importorskip("tangermeme.deep_lift_shap", exc_type=ImportError)
    from tangermeme.deep_lift_shap import deep_lift_shap
    from tangermeme.utils import one_hot_encode

    model = _count_a_model(torch)
    sequences = {"a": "ACGTTGCAAC", "b": "GGGATCCATG"}
    result = deep_lift_shap_attribution(model, sequences, n_shuffles=4, random_state=3)
    inputs = torch.stack([one_hot_encode(text) for text in sequences.values()]).float()
    direct = deep_lift_shap(model, inputs, target=0, n_shuffles=4, random_state=3)
    assert result["shape"] == [2, 10, 4]
    assert _flat(result["attributions"]) == pytest.approx(_flat(direct.permute(0, 2, 1).tolist()))
    again = deep_lift_shap_attribution(model, sequences, n_shuffles=4, random_state=3)
    assert again["attributions"] == result["attributions"]


def test_tangermeme_helpers_validate_models_and_inputs() -> None:
    torch = pytest.importorskip("torch", exc_type=ImportError)
    pytest.importorskip("tangermeme.predict", exc_type=ImportError)
    model = _count_a_model(torch)
    with pytest.raises(ConfigurationError, match=r"torch\.nn\.Module"):
        predict_sequences(object(), {"s": "ACGT"})
    model.train()
    with pytest.raises(ConfigurationError, match="eval"):
        predict_sequences(model, {"s": "ACGT"})
    model.eval()
    with pytest.raises(ConfigurationError, match="same length"):
        predict_sequences(model, {"a": "ACGT", "b": "ACG"})
    with pytest.raises(ConfigurationError, match="max_positions"):
        deep_lift_shap_attribution(model, {"s": "ACGTACGT"}, max_positions=4)
    with pytest.raises(ConfigurationError, match="longer"):
        marginalize_motif_effect(model, {"s": "ACGT"}, "ACGTA")


def test_missing_tangermeme_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    torch = pytest.importorskip("torch", exc_type=ImportError)
    model = _count_a_model(torch)
    monkeypatch.setitem(sys.modules, "tangermeme.predict", None)
    with pytest.raises(BackendUnavailableError, match="bio-tangermeme"):
        predict_sequences(model, {"s": "ACGT"})
