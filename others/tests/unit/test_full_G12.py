"""Tests for the gimmemotifs, tangermeme and snapATAC2 wrappers in dnakit.regulation (group G12)."""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.regulation import gimmemotifs_functions as gm
from dnakit.regulation import snapatac2_functions as sa
from dnakit.regulation import tangermeme_functions as tm

GENOME = {"chr1": 3000, "chr2": 2000}
SIZES = {"chr1": 3000, "chr2": 2000, "chrM": 1000}

FRAGMENTS = """\
chr1\t100\t300\tcellA\t1
chr1\t120\t330\tcellA\t2
chr1\t600\t800\tcellA\t1
chr1\t1050\t1200\tcellA\t1
chr2\t50\t250\tcellA\t1
chrM\t10\t200\tcellA\t1
chr1\t650\t900\tcellB\t1
chr2\t50\t250\tcellB\t1
chr2\t1200\t1400\tcellB\t3
chr1\t2000\t2200\tcellB\t1
chr1\t150\t310\tcellB\t1
chr2\t100\t400\tcellC\t1
chr1\t2400\t2600\tcellC\t1
chr2\t1500\t1700\tcellC\t1
chr1\t120\t200\tcellC\t1
chr1\t1100\t1280\tcellC\t1
"""

GTF_ROWS = (
    ("chr1", "gene", 90, 400, "+", 'gene_id "G1"; gene_name "GENEA";'),
    ("chr1", "transcript", 90, 400, "+", 'gene_id "G1"; transcript_id "T1"; gene_name "GENEA";'),
    ("chr1", "exon", 90, 400, "+", 'gene_id "G1"; transcript_id "T1"; gene_name "GENEA";'),
    ("chr1", "gene", 1000, 1300, "-", 'gene_id "G2"; gene_name "GENEB";'),
    ("chr1", "transcript", 1000, 1300, "-", 'gene_id "G2"; transcript_id "T2"; gene_name "GENEB";'),
    ("chr1", "exon", 1000, 1300, "-", 'gene_id "G2"; transcript_id "T2"; gene_name "GENEB";'),
    ("chr2", "gene", 400, 1600, "+", 'gene_id "G3"; gene_name "GENEC";'),
    ("chr2", "transcript", 400, 1600, "+", 'gene_id "G3"; transcript_id "T3"; gene_name "GENEC";'),
    ("chr2", "exon", 400, 1600, "+", 'gene_id "G3"; transcript_id "T3"; gene_name "GENEC";'),
)
GTF = "".join(
    f"{chrom}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
    for chrom, kind, start, end, strand, attributes in GTF_ROWS
)


def _flat(values: object) -> list[float]:
    """Flatten nested lists so pytest.approx can compare them."""

    if isinstance(values, (list, tuple)):
        return [item for value in values for item in _flat(value)]
    return [float(values)]  # type: ignore[arg-type]


def _random_genome(seed: int = 7) -> dict[str, str]:
    rng = random.Random(seed)
    return {
        name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in GENOME.items()
    }


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        for name, text in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(text), 60):
                handle.write(text[start : start + 60] + "\n")
    return path


def _genome_file(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    records = _random_genome()
    return _write_fasta(tmp_path / "genome.fa", records), records


def _fragment_file(tmp_path: Path) -> Path:
    path = tmp_path / "fragments.tsv"
    path.write_text(FRAGMENTS, encoding="utf-8")
    return path


def _direct_import(snap, tmp_path: Path, fragments: Path):
    """Import fragments with snapATAC2 directly, for comparison with the wrapper."""

    return snap.pp.import_fragments(
        fragments,
        SIZES,
        is_paired=True,
        min_num_fragments=1,
        file=tmp_path / "direct.h5ad",
        sorted_by_barcode=False,
        n_jobs=1,
    )


def _gtf_file(tmp_path: Path) -> Path:
    path = tmp_path / "genes.gtf"
    path.write_text(GTF, encoding="utf-8")
    return path


def _pfm_file(tmp_path: Path) -> Path:
    path = tmp_path / "motifs.pfm"
    path.write_text(
        ">GTTT_demo\n0 0 10 0\n0 0 0 10\n0 0 0 10\n0 0 0 10\n"
        ">ACGT_demo\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n",
        encoding="utf-8",
    )
    return path


def _consensus_pfm_file(tmp_path: Path, consensus: str = "GATAAG", counts: int = 1000) -> Path:
    """A PFM where the consensus base gets all counts, written in gimmemotifs' PFM layout."""

    path = tmp_path / "gata.pfm"
    rows = []
    for base in consensus:
        values = [0, 0, 0, 0]
        values["ACGT".index(base)] = counts
        rows.append(" ".join(str(v) for v in values))
    path.write_text(">GATA_consensus\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def _summit_rows(tmp_path: Path) -> Path:
    path = tmp_path / "summits.bed"
    path.write_text(
        "chr1\t250\t251\tsummit1\t9.5\nchr2\t600\t601\tsummit2\t7.0\n", encoding="utf-8"
    )
    return path


def _narrowpeak(tmp_path: Path) -> Path:
    path = tmp_path / "peaks.narrowPeak"
    path.write_text(
        "chr1\t100\t400\tp1\t10\t+\t200\t5.0\t3.0\t150\n"
        "chr1\t1000\t1300\tp2\t8\t+\t120\t4.0\t2.0\t100\n"
        "chr2\t500\t900\tp3\t9\t+\t300\t4.5\t2.5\t200\n",
        encoding="utf-8",
    )
    return path


def _peak_bed(tmp_path: Path) -> Path:
    path = tmp_path / "peaks.bed"
    path.write_text(
        "chr1\t100\t400\tp1\t10\t+\nchr1\t1000\t1300\tp2\t8\t+\nchr2\t500\t900\tp3\t9\t+\n",
        encoding="utf-8",
    )
    return path


def _gene_bed6(tmp_path: Path) -> Path:
    path = tmp_path / "genes6.bed"
    path.write_text(
        "chr1\t500\t900\tG1\t0\t+\nchr1\t2000\t2500\tG2\t0\t-\nchr2\t400\t800\tG3\t0\t+\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# gimmemotifs
# ---------------------------------------------------------------------------


def test_motif_from_consensus_matches_hand_counts() -> None:
    pytest.importorskip("gimmemotifs")
    result = gm.motif_from_consensus("GATAAG", n_sites=100, verbose=False)
    assert result["status"] == "completed"
    assert result["consensus"] == "GATAAG"
    assert result["pfm"] == [
        [0, 0, 100, 0],
        [100, 0, 0, 0],
        [0, 0, 0, 100],
        [100, 0, 0, 0],
        [100, 0, 0, 0],
        [0, 0, 100, 0],
    ]
    assert result["provenance"]["package"] == "gimmemotifs"
    with pytest.raises(ConfigurationError):
        gm.motif_from_consensus("GATXQ", verbose=False)


def test_motif_from_alignment_counts_and_direct_upstream_probabilities() -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.motif import motif_from_align

    sequences = ["ACGT", "ACGA", "ACGT"]
    result = gm.motif_from_alignment(sequences, verbose=False)
    # Column 4 holds T, A, T: two T and one A.
    assert result["pfm"] == [[3, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [1, 0, 0, 2]]
    direct = motif_from_align(sequences)
    assert result["consensus"] == str(direct.consensus)
    assert result["ppm"][3][0] == pytest.approx(float(direct.ppm[3][0]), abs=1e-5)
    with pytest.raises(ConfigurationError):
        gm.motif_from_alignment(["ACGT", "ACG"], verbose=False)


def test_compare_motif_positions_matches_direct_akl_and_family_metrics() -> None:
    pytest.importorskip("gimmemotifs")
    import numpy as np
    from gimmemotifs import comparison
    from gimmemotifs.utils import ppm_pseudocount

    pfm_a = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0]]
    pfm_b = [[10, 0, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
    ppm_a = ppm_pseudocount(np.asarray(pfm_a, dtype=float) / 10.0)
    ppm_b = ppm_pseudocount(np.asarray(pfm_b, dtype=float) / 10.0)
    direct = [comparison.akl(a, b) for a, b in zip(ppm_a, ppm_b, strict=True)]
    result = gm.compare_motif_positions(pfm_a, pfm_b, metric="akl", verbose=False)
    assert result["position_scores"] == pytest.approx([float(v) for v in direct])
    assert result["position_scores"][0] == pytest.approx(10.0, abs=1e-3)
    chisq = gm.compare_motif_positions(pfm_a, pfm_b, metric="chisq", verbose=False)
    ssd = gm.compare_motif_positions(pfm_a, pfm_b, metric="ssd", verbose=False)
    assert chisq["metric"] == "chisq" and ssd["metric"] == "ssd"
    assert ssd["position_scores"][0] == pytest.approx(2.0, abs=1e-3)
    assert ssd["position_scores"][1] == pytest.approx(2.0 - 2.0 * 1.0, abs=1e-3)
    with pytest.raises(ConfigurationError):
        gm.compare_motif_positions(pfm_a, pfm_b, metric="euclid", verbose=False)


def test_motif_score_correlation_matches_direct_seqcor() -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.comparison import seqcor
    from gimmemotifs.motif import Motif

    matrix_a = [[0, 0, 10, 0], [0, 0, 0, 10], [10, 0, 0, 0], [0, 0, 0, 10]]
    matrix_b = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
    score, position, strand = seqcor(Motif(matrix_a), Motif(matrix_b))
    result = gm.motif_score_correlation(matrix_a, matrix_b, verbose=False)
    assert result["score"] == pytest.approx(float(score))
    assert result["position"] == int(position)
    assert result["strand"] == int(strand)
    assert result["length_a"] == 4 and result["length_b"] == 4


def test_scan_best_motif_hits_places_forward_and_reverse_sites(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    motifs = _consensus_pfm_file(tmp_path)
    forward = "A" * 10 + "GATAAG" + "C" * 14
    reverse = "T" * 7 + "CTTATC" + "T" * 15  # reverse complement of GATAAG at index 7
    fasta = tmp_path / "scan.fa"
    fasta.write_text(f">fwd\n{forward}\n>rev\n{reverse}\n", encoding="utf-8")
    result = gm.scan_best_motif_hits(fasta, motifs, verbose=False)
    assert result["sequence_count"] == 2 and result["motif_count"] == 1
    by_sequence = {hit["sequence"]: hit for hit in result["hits"]}
    assert by_sequence["fwd"]["start"] == 10
    assert by_sequence["fwd"]["end"] == 16
    assert by_sequence["fwd"]["strand"] == "+"
    assert by_sequence["rev"]["strand"] == "-"
    assert by_sequence["rev"]["start"] == 7
    assert by_sequence["rev"]["end"] == 13
    # Hits follow FASTA order, then motif ID.
    assert [hit["sequence"] for hit in result["hits"]] == ["fwd", "rev"]


def test_scan_best_motif_hits_agrees_with_direct_upstream(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.scanner import scan_to_best_match

    motifs = _pfm_file(tmp_path)
    fasta = tmp_path / "seqs.fa"
    fasta.write_text(
        ">f1\nACGTTGATAAGCCGTAGGATAAGTTCAG\n>f2\nTTGATAAGGGCCATTAGAAGTCCTAGGA\n", encoding="utf-8"
    )
    direct = scan_to_best_match(
        str(fasta),
        pfmfile=str(motifs),
        ncpus=1,
        score=False,
        zscore=False,
        gc=False,
        progress=False,
    )
    result = gm.scan_best_motif_hits(fasta, motifs, verbose=False)
    for hit in result["hits"]:
        index = ["f1", "f2"].index(hit["sequence"])
        score, start, strand = direct[hit["motif"]][index]
        assert hit["score"] == pytest.approx(float(score))
        assert hit["start"] == int(start)
        assert hit["strand"] == ("+" if int(strand) > 0 else "-")


def test_motif_enrichment_stats_matches_direct_calc_stats(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.stats import calc_stats

    motifs = _pfm_file(tmp_path)
    fg = tmp_path / "fg.fa"
    bg = tmp_path / "bg.fa"
    fg.write_text(
        ">f1\nACGTTGATAAGCCGTAGGATAAGTTCAG\n>f2\nTTGATAAGGGCCATTAGAAGTCCTAGGA\n", encoding="utf-8"
    )
    bg.write_text(
        ">b1\nGGGCCCAAATTTGGGCCCAAATTTGGGC\n>b2\nATATATCGCGCGTATATAGCGCGTATAT\n", encoding="utf-8"
    )
    direct = calc_stats(
        fg_file=str(fg), bg_file=str(bg), motifs=str(motifs), zscore=False, gc=False, ncpus=1
    )
    result = gm.motif_enrichment_stats(fg, bg, motifs, verbose=False)
    assert [row["motif"] for row in result["motifs"]] == sorted(direct)
    for row in result["motifs"]:
        assert row["stats"]["roc_auc"] == pytest.approx(float(direct[row["motif"]]["roc_auc"]))
        assert 0.0 <= row["stats"]["roc_auc"] <= 1.0


def test_cluster_motif_set_matches_direct_tree(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.motif.cluster import cluster_motifs

    motifs = _pfm_file(tmp_path)
    root = cluster_motifs(
        str(motifs),
        match="total",
        metric="wic",
        combine="mean",
        pval=True,
        threshold=0.95,
        progress=False,
        ncpus=1,
    )
    expected = sorted(sorted(str(m.id) for m in members) for _rep, members in root.getResult())
    result = gm.cluster_motif_set(motifs, threshold=0.95, verbose=False)
    assert sorted(item["members"] for item in result["clusters"]) == expected
    assert result["motif_count"] == 2
    single = tmp_path / "single.pfm"
    single.write_text(">only\n1 0 0 0\n0 1 0 0\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        gm.cluster_motif_set(single, verbose=False)


def test_read_peak_summits_reads_summit_and_narrowpeak_rows(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    summits = _summit_rows(tmp_path)
    result = gm.read_peak_summits(summits, verbose=False)
    assert result["peak_count"] == 2
    first = result["peaks"][0]
    assert (first["chrom"], first["start"], first["end"], first["name"]) == (
        "chr1",
        250,
        251,
        "summit1",
    )
    assert first["log_value"] == pytest.approx(math.log1p(9.5))
    narrow = gm.read_peak_summits(_narrowpeak(tmp_path), verbose=False)
    assert narrow["peak_count"] == 3
    assert narrow["peaks"][0]["start"] == 250  # 100 + summit offset 150
    with pytest.raises(ConfigurationError):
        gm.read_peak_summits(_peak_bed(tmp_path), verbose=False)


def test_summit_intervals_from_narrowpeak_centres_on_summit(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    result = gm.summit_intervals_from_narrowpeak(_narrowpeak(tmp_path), size=50, verbose=False)
    # Summit 100 + 150 = 250; width 50 starts at 250 - 25 = 225.
    assert result["intervals"][0]["chrom"] == "chr1"
    assert (result["intervals"][0]["start"], result["intervals"][0]["end"]) == (225, 275)
    assert result["interval_count"] == 3
    with pytest.raises(ConfigurationError):
        gm.summit_intervals_from_narrowpeak(_peak_bed(tmp_path), size=50, verbose=False)


def test_equal_size_intervals_centres_on_midpoint(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    result = gm.equal_size_intervals(_peak_bed(tmp_path), size=100, verbose=False)
    spans = [(item["chrom"], item["start"], item["end"]) for item in result["intervals"]]
    # Midpoints 250, 1150 and 700 give intervals of 100 bases starting 50 before the midpoint.
    assert spans == [("chr1", 200, 300), ("chr1", 1100, 1200), ("chr2", 650, 750)]
    full = gm.equal_size_intervals(_narrowpeak(tmp_path), size=100, verbose=False)
    trimmed = gm.equal_size_intervals(
        _narrowpeak(tmp_path), size=100, keep_columns=False, verbose=False
    )
    assert len(full["intervals"][0]["columns"]) == 7
    assert len(trimmed["intervals"][0]["columns"]) == 3


def test_random_genomic_regions_matches_direct_upstream(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from gimmemotifs.background import create_random_genomic_bedfile

    genome, _records = _genome_file(tmp_path)
    direct_path = tmp_path / "direct.bed"
    create_random_genomic_bedfile(str(direct_path), str(genome), 100, 5, seed=1)
    direct = sorted(
        (fields[0], int(fields[1]), int(fields[2]))
        for fields in (line.split("\t") for line in direct_path.read_text().splitlines())
    )
    result = gm.random_genomic_regions(genome, 100, 5, seed=1, verbose=False)
    observed = [(item["chrom"], item["start"], item["end"]) for item in result["intervals"]]
    assert observed == direct
    assert all(end - start == 100 for _chrom, start, end in observed)
    assert all(end <= GENOME[chrom] for chrom, _start, end in observed)


def test_promoter_regions_hand_windows_and_restores_random_state(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    before = random.getstate()
    result = gm.promoter_regions(_gene_bed6(tmp_path), 100, 3, seed=1, verbose=False)
    assert random.getstate() == before
    windows = [
        (item["chrom"], item["start"], item["end"], item["strand"]) for item in result["intervals"]
    ]
    # Plus-strand genes: 100 bases before the start. Minus-strand gene: 100 bases after the end.
    assert windows == [
        ("chr1", 400, 500, "+"),
        ("chr1", 2500, 2600, "-"),
        ("chr2", 300, 400, "+"),
    ]
    with pytest.raises(ConfigurationError):
        gm.promoter_regions(_summit_rows(tmp_path), 100, 3, verbose=False)


def test_gimmemotifs_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("gimmemotifs")
    for name in ("gimmemotifs", "gimmemotifs.motif"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError) as info:
        gm.motif_from_consensus("GATAAG", verbose=False)
    assert "bio-gimmemotifs" in (info.value.hint or "")


# ---------------------------------------------------------------------------
# tangermeme
# ---------------------------------------------------------------------------


def _linear_model(length: int = 20):
    """Output is the sum over positions of W[base, position] for the observed one-hot bases."""

    torch = pytest.importorskip("torch")
    weights = torch.tensor(
        [[((base * 7 + pos * 3) % 11) / 4.0 - 1.0 for pos in range(length)] for base in range(4)],
        dtype=torch.float32,
    )

    class Linear(torch.nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.register_buffer("weights", weights)

        def forward(self, x, *args):
            return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

    model = Linear().eval()
    return model, weights


SEQS = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}


def test_substitution_effect_matches_weight_difference() -> None:
    pytest.importorskip("tangermeme")
    model, weights = _linear_model()
    result = tm.substitution_effect(
        model, SEQS, [{"sequence": "s1", "position": 3, "base": "G"}], verbose=False
    )
    # Original base at s1[3] is T; the change is W[G, 3] - W[T, 3].
    expected = float(weights[2, 3] - weights[3, 3])
    assert result["delta"][0][0] == pytest.approx(expected, abs=1e-5)
    assert result["delta"][1][0] == pytest.approx(0.0, abs=1e-6)
    assert result["variant_count"] == 1


def test_saturation_mutagenesis_matches_hand_centred_weights() -> None:
    pytest.importorskip("tangermeme")
    model, weights = _linear_model()
    result = tm.saturation_mutagenesis(model, SEQS, start=0, end=6, verbose=False)
    assert result["shape"] == [2, 4, 6]
    observed = {"s1": SEQS["s1"], "s2": SEQS["s2"]}
    for index, name in enumerate(result["sequences"]):
        for position in range(6):
            base = "ACGT".index(observed[name][position])
            column = weights[:, position]
            expected = float(column[base] - column.mean())
            assert result["attributions"][index][base][position] == pytest.approx(
                expected, abs=1e-5
            )
            for other in range(4):
                if other != base:
                    assert result["attributions"][index][other][position] == pytest.approx(0.0)


def test_saturation_mutagenesis_matches_direct_upstream() -> None:
    pytest.importorskip("tangermeme")
    import torch
    from tangermeme.saturation_mutagenesis import saturation_mutagenesis
    from tangermeme.utils import one_hot_encode

    model, _weights = _linear_model()
    inputs = torch.stack([one_hot_encode(SEQS["s1"]), one_hot_encode(SEQS["s2"])]).float()
    direct = saturation_mutagenesis(model, inputs, start=2, end=9, batch_size=8)
    result = tm.saturation_mutagenesis(model, SEQS, start=2, end=9, batch_size=8, verbose=False)
    assert _flat(result["attributions"]) == pytest.approx(_flat(direct.tolist()), abs=1e-6)


def test_insertion_and_deletion_effects_match_direct_upstream() -> None:
    pytest.importorskip("tangermeme")
    import torch
    from tangermeme.utils import one_hot_encode
    from tangermeme.variant_effect import deletion_effect, insertion_effect

    model, _weights = _linear_model(length=20)
    inserted = tm.insertion_effect(
        model, SEQS, [{"sequence": "s2", "position": 4, "base": "C"}], verbose=False
    )
    inputs = torch.stack([one_hot_encode(SEQS["s1"]), one_hot_encode(SEQS["s2"])]).float()
    direct = insertion_effect(model, inputs, torch.tensor([[1, 4, 1]]))
    assert _flat(inserted["y_after"]) == pytest.approx(_flat(direct.y_after.tolist()), abs=1e-6)
    # Deletions need a window longer than the model input: 21 bases for a 20-base model.
    model21, _ = _linear_model(length=20)
    long_seqs = {"s1": SEQS["s1"] + "A", "s2": SEQS["s2"] + "A"}
    deleted = tm.deletion_effect(
        model21, long_seqs, [{"sequence": "s1", "position": 5}], verbose=False
    )
    inputs21 = torch.stack(
        [one_hot_encode(long_seqs["s1"]), one_hot_encode(long_seqs["s2"])]
    ).float()
    direct_del = deletion_effect(model21, inputs21, torch.tensor([[0, 5]]))
    assert _flat(deleted["y_after"]) == pytest.approx(_flat(direct_del.y_after.tolist()), abs=1e-6)


def test_ablate_region_and_spacing_match_direct_upstream() -> None:
    pytest.importorskip("tangermeme")
    import torch
    from tangermeme.ablate import ablate
    from tangermeme.ersatz import dinucleotide_shuffle
    from tangermeme.space import space
    from tangermeme.utils import one_hot_encode

    model, _weights = _linear_model()
    inputs = torch.stack([one_hot_encode(SEQS["s1"]), one_hot_encode(SEQS["s2"])]).float()
    ablated = tm.ablate_region(model, SEQS, start=2, end=14, n=3, random_state=5, verbose=False)
    direct = ablate(
        model, inputs, start=2, end=14, n=3, shuffle_fn=dinucleotide_shuffle, random_state=5
    )
    assert _flat(ablated["y_after"]) == pytest.approx(_flat(direct.y_after.tolist()), abs=1e-6)
    assert len(ablated["y_after"]) == 2 and len(ablated["y_after"][0]) == 3

    spaced = tm.motif_pair_spacing(model, SEQS, ["GG", "TT"], [[0], [2]], start=2, verbose=False)
    direct_space = space(model, inputs, ["GG", "TT"], torch.tensor([[0], [2]]), start=2)
    for row, expected in zip(spaced["y_after"], direct_space.y_afters, strict=True):
        assert _flat(row) == pytest.approx(_flat(expected.tolist()), abs=1e-6)


def test_design_greedy_substitution_matches_direct_upstream() -> None:
    pytest.importorskip("tangermeme")
    import torch
    from tangermeme.design import greedy_substitution
    from tangermeme.utils import characters, one_hot_encode

    model, _weights = _linear_model()
    sequence = SEQS["s1"]
    result = tm.design_greedy_substitution(
        model, sequence, [5.0], motifs=["GGG"], max_iter=2, verbose=False
    )
    designed = greedy_substitution(
        model,
        one_hot_encode(sequence)[None].float(),
        y=torch.tensor([[5.0]]),
        motifs=["GGG"],
        max_iter=2,
        tol=0.001,
    )
    assert result["designed_sequence"] == characters(designed[0])
    assert len(result["designed_sequence"]) == len(sequence)


def test_call_seqlets_finds_block_and_matches_upstream() -> None:
    pytest.importorskip("tangermeme")
    import numpy as np
    import torch
    from tangermeme.seqlet import recursive_seqlets

    rng = np.random.RandomState(0)
    values = rng.randn(20, 300) * 0.5
    values[0, 100:112] = 8.0
    result = tm.call_seqlets(values.tolist(), threshold=0.01, verbose=False)
    direct = recursive_seqlets(torch.tensor(values).float(), threshold=0.01)
    assert result["seqlet_count"] == len(direct)
    first = result["seqlets"][0]
    assert (first["example_index"], first["start"], first["end"]) == (0, 97, 115)
    assert first["start"] <= 100 and first["end"] >= 112
    assert first["attribution"] == pytest.approx(float(direct.iloc[0]["attribution"]), rel=1e-5)
    assert first["p_value"] == pytest.approx(float(direct.iloc[0]["p-value"]), rel=1e-5)


def test_kmer_counts_hand_values_and_weighted_sums() -> None:
    pytest.importorskip("tangermeme")
    counts = tm.kmer_counts({"s": "ACGT"}, 2, verbose=False)
    assert [(item["kmer"], item["value"]) for item in counts["entries"]] == [
        ("AC", 1.0),
        ("CG", 1.0),
        ("GT", 1.0),
    ]
    weighted = tm.kmer_counts({"s": "ACGT"}, 2, scores={"s": [1.0, 2.0, 3.0, 4.0]}, verbose=False)
    # Each 2-mer adds the sum of the two scores it spans: AC = 1 + 2, CG = 2 + 3, GT = 3 + 4.
    assert [item["value"] for item in weighted["entries"]] == [3.0, 5.0, 7.0]
    assert weighted["weighted"] is True


def test_extract_loci_sequences_cuts_the_midpoint_window(tmp_path: Path) -> None:
    pytest.importorskip("tangermeme")
    genome, records = _genome_file(tmp_path)
    loci = tmp_path / "loci.bed"
    loci.write_text("chr1\t100\t120\nchr2\t500\t520\n", encoding="utf-8")
    result = tm.extract_loci_sequences(loci, genome, in_window=10, out_window=8, verbose=False)
    assert result["shape"] == [2, 4, 10]
    # Midpoints are 110 and 510; a 10-base window starts 5 bases before each midpoint.
    assert result["loci"][0]["sequence"] == records["chr1"][105:115]
    assert result["loci"][1]["sequence"] == records["chr2"][505:515]


def test_matched_background_loci_is_deterministic(tmp_path: Path) -> None:
    pytest.importorskip("tangermeme")
    genome, _records = _genome_file(tmp_path)
    loci = tmp_path / "loci.bed"
    loci.write_text("chr1\t100\t120\nchr2\t500\t520\n", encoding="utf-8")
    first = tm.matched_background_loci(
        loci, genome, in_window=40, out_window=20, random_state=0, verbose=False
    )
    second = tm.matched_background_loci(
        loci, genome, in_window=40, out_window=20, random_state=0, verbose=False
    )
    assert first["loci"] == second["loci"]
    assert first["background_count"] == first["input_locus_count"] == 2
    assert all(item["end"] - item["start"] == 40 for item in first["loci"])


def test_pwm_information_content_hand_bits() -> None:
    pytest.importorskip("tangermeme")
    result = tm.pwm_information_content(
        [[0.25, 0.25, 0.25, 0.25], [1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.0, 0.0]], verbose=False
    )
    assert result["information_content"] == pytest.approx([0.0, 2.0, 1.0], abs=1e-6)
    assert result["total"] == pytest.approx(3.0, abs=1e-6)


def test_place_motif_substitutes_and_inserts_by_hand() -> None:
    pytest.importorskip("tangermeme")
    substituted = tm.place_motif(SEQS, "GGG", start=5, verbose=False)
    assert substituted["sequences"][0]["sequence"] == "ACGTTGGGAAGCCGTAGGAT"
    inserted = tm.place_motif(SEQS, "GGG", start=5, mode="insert", verbose=False)
    # Insert GGG before base 5, then keep the first 20 bases.
    assert inserted["sequences"][0]["sequence"] == "ACGTTGGGGATAAGCCGTAG"
    with pytest.raises(ConfigurationError):
        tm.place_motif(SEQS, "GGG", mode="delete", verbose=False)


def test_tangermeme_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("tangermeme")
    for name in ("tangermeme", "tangermeme.kmers", "tangermeme.utils"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError) as info:
        tm.kmer_counts({"s": "ACGT"}, 2, verbose=False)
    assert "bio-tangermeme" in (info.value.hint or "")


# ---------------------------------------------------------------------------
# snapATAC2
# ---------------------------------------------------------------------------


def test_tss_enrichment_matches_direct_upstream(tmp_path: Path) -> None:
    snap = pytest.importorskip("snapatac2")
    fragments = _fragment_file(tmp_path)
    gtf = _gtf_file(tmp_path)
    result = sa.tss_enrichment(fragments, SIZES, gtf, verbose=False)
    assert [cell["barcode"] for cell in result["cells"]] == ["cellA", "cellB", "cellC"]
    adata = _direct_import(snap, tmp_path, fragments)
    try:
        direct = snap.metrics.tsse(
            adata, gtf, exclude_chroms=["chrM", "M"], inplace=False, n_jobs=1
        )
    finally:
        adata.close()
    assert [cell["tsse"] for cell in result["cells"]] == pytest.approx(
        list(direct["tsse"]), abs=1e-9
    )
    assert result["library_tsse"] == pytest.approx(float(direct["library_tsse"]), abs=1e-9)
    assert result["frac_overlap_TSS"] == pytest.approx(float(direct["frac_overlap_TSS"]), abs=1e-9)


def test_fragment_in_peaks_fraction_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2")
    fragments = _fragment_file(tmp_path)
    peaks = _peak_bed(tmp_path)
    result = sa.fragment_in_peaks_fraction(fragments, SIZES, peaks, verbose=False)
    values = {cell["barcode"]: cell["frip"] for cell in result["cells"]}
    # cellA: 3 of 6 fragments in peaks; cellB: 1 of 5; cellC: 2 of 5.
    assert values["cellA"] == pytest.approx(0.5)
    assert values["cellB"] == pytest.approx(0.2)
    assert values["cellC"] == pytest.approx(0.4)


def test_fragment_size_distribution_matches_hand_lengths(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2")
    fragments = _fragment_file(tmp_path)
    result = sa.fragment_size_distribution(fragments, SIZES, max_recorded_size=400, verbose=False)
    counts = result["counts"]
    assert result["fragment_count"] == 16
    assert (
        counts[200] == 8
    )  # 100-300, 600-800, 50-250 x2, 1200-1400, 2000-2200, 2400-2600, 1500-1700
    assert counts[210] == 1  # 120-330
    assert counts[190] == 1  # chrM 10-200
    assert counts[80] == 1  # 120-200
    assert counts[250] == 1  # 650-900


def test_chromosome_fragment_counts_match_hand_tallies(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2")
    fragments = _fragment_file(tmp_path)
    result = sa.chromosome_fragment_counts(fragments, SIZES, verbose=False)
    assert result["chromosomes"] == ["chr1", "chr2", "chrM"]
    rows = {cell["barcode"]: cell["counts"] for cell in result["cells"]}
    assert rows == {"cellA": [4.0, 1.0, 1.0], "cellB": [3.0, 2.0, 0.0], "cellC": [3.0, 2.0, 0.0]}


def test_gene_activity_counts_match_direct_upstream(tmp_path: Path) -> None:
    snap = pytest.importorskip("snapatac2")
    fragments = _fragment_file(tmp_path)
    gtf = _gtf_file(tmp_path)
    result = sa.gene_activity_counts(fragments, SIZES, gtf, upstream=500, verbose=False)
    assert result["genes"] == ["GENEA", "GENEB", "GENEC"]
    adata = _direct_import(snap, tmp_path, fragments)
    try:
        matrix = snap.pp.make_gene_matrix(
            adata,
            gtf,
            inplace=False,
            file=tmp_path / "direct_genes.h5ad",
            upstream=500,
            downstream=0,
        )
        try:
            values = matrix.X[:]
            dense = values.toarray() if hasattr(values, "toarray") else values
            direct = [
                [row, column, int(dense[row, column])]
                for row in range(dense.shape[0])
                for column in range(dense.shape[1])
                if dense[row, column] != 0
            ]
            genes = [str(item) for item in matrix.var_names]
        finally:
            matrix.close()
    finally:
        adata.close()
    assert genes == result["genes"]
    assert result["entries"] == direct


def test_regulatory_regions_to_genes_links_overlapping_genes(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2")
    gtf = _gtf_file(tmp_path)
    result = sa.regulatory_regions_to_genes(
        ["chr1:100-300", "chr1:1000-1200", "chr2:500-700"],
        gtf,
        upstream=250,
        downstream=250,
        verbose=False,
    )
    pairs = [(edge["region"], edge["gene"]) for edge in result["edges"]]
    assert pairs == [
        ("chr1:100-300", "GENEA"),
        ("chr1:1000-1200", "GENEB"),
        ("chr2:500-700", "GENEC"),
    ]
    assert all(edge["distance"] >= 0 for edge in result["edges"])
    with pytest.raises(ConfigurationError):
        sa.regulatory_regions_to_genes(["chr1:300-100"], gtf, verbose=False)


def test_region_motif_enrichment_matches_direct_upstream(tmp_path: Path) -> None:
    snap = pytest.importorskip("snapatac2")
    genome, _records = _genome_file(tmp_path)
    motifs = tmp_path / "gata.meme"
    motifs.write_text(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF GATA_demo\nletter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0\n"
        " 0.0 0.0 1.0 0.0\n 0.0 0.0 0.0 1.0\n 1.0 0.0 0.0 0.0\n"
        " 0.0 0.0 0.0 1.0\n 0.0 0.0 0.0 1.0\n 0.0 0.0 0.0 1.0\n",
        encoding="utf-8",
    )
    regions = {"g1": ["chr1:100-300", "chr1:1000-1200"], "g2": ["chr2:500-700"]}
    result = sa.region_motif_enrichment(regions, motifs, genome, verbose=False)
    direct = snap.tl.motif_enrichment(
        snap.read_motifs(str(motifs)), regions, genome_fasta=genome, background=None, method=None
    )
    assert [group["group"] for group in result["groups"]] == ["g1", "g2"]
    for group in result["groups"]:
        row = group["rows"][0]
        expected = direct[group["group"]].to_dicts()[0]
        assert row["p-value"] == pytest.approx(float(expected["p-value"]))
        assert row["id"] == expected["id"]


def test_merge_peak_sets_centres_on_summits(tmp_path: Path) -> None:
    pytest.importorskip("snapatac2")
    group_a = [
        {
            "chrom": "chr1",
            "start": 100,
            "end": 400,
            "name": "p",
            "score": 10,
            "strand": ".",
            "signal_value": 5.0,
            "p_value": 3.0,
            "q_value": 2.0,
            "summit": 200,
        }
    ]
    group_b = [dict(group_a[0], chrom="chr2", start=60, end=300, summit=100)]
    result = sa.merge_peak_sets({"A": group_a, "B": group_b}, GENOME, half_width=50, verbose=False)
    spans = [
        (item["chrom"], item["start"], item["end"], item["groups"]) for item in result["peaks"]
    ]
    # Summit 300 with half_width 50 gives [250, 351); summit 160 gives [110, 211).
    assert spans == [("chr1", 250, 351, ["A"]), ("chr2", 110, 211, ["B"])]


def test_snapatac2_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "snapatac2", None)
    with pytest.raises(BackendUnavailableError) as info:
        sa.merge_peak_sets({"A": []}, GENOME, verbose=False)
    assert "bio-snapatac2" in (info.value.hint or "")


def test_public_names_are_sorted_and_take_verbose() -> None:
    import inspect

    for module in (gm, tm, sa):
        assert module.__all__ == sorted(module.__all__)
        for name in module.__all__:
            parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is True
