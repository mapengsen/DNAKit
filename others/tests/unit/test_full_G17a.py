"""Hand-checked and upstream-checked tests for the G17a predictions wrappers.

Covers alphagenome, grelu and enformer-pytorch helpers. Every wrapped function is run
against real upstream code on local inputs; no network access is needed.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.predictions import alphagenome_functions as ag
from dnakit.predictions import enformer_pytorch_functions as ef
from dnakit.predictions import grelu_functions as gr

GTF_ROWS = (
    (
        "gene",
        1001,
        2000,
        "+",
        {"gene_id": "ENSG1.1", "gene_type": "protein_coding", "gene_name": "GENEA"},
    ),
    (
        "transcript",
        1001,
        1800,
        "+",
        {
            "gene_id": "ENSG1.1",
            "transcript_id": "ENST1.1",
            "gene_type": "protein_coding",
            "gene_name": "GENEA",
            "transcript_type": "protein_coding",
            "transcript_support_level": "1",
            "tag": "basic MANE_Select",
        },
    ),
    (
        "transcript",
        1101,
        2000,
        "+",
        {
            "gene_id": "ENSG1.1",
            "transcript_id": "ENST2.1",
            "gene_type": "protein_coding",
            "gene_name": "GENEA",
            "transcript_type": "protein_coding",
            "transcript_support_level": "3",
            "tag": "basic",
        },
    ),
    ("gene", 5001, 5500, "-", {"gene_id": "ENSG2.1", "gene_type": "lncRNA", "gene_name": "GENEB"}),
    (
        "transcript",
        5001,
        5500,
        "-",
        {
            "gene_id": "ENSG2.1",
            "transcript_id": "ENST3.1",
            "gene_type": "lncRNA",
            "gene_name": "GENEB",
            "transcript_type": "lncRNA",
            "transcript_support_level": "NA",
            "tag": "basic",
        },
    ),
)


def _gtf_text() -> str:
    lines = ["##description: small annotation for G17a tests"]
    for feature, start, end, strand, attributes in GTF_ROWS:
        body = " ".join(f'{key} "{value}";' for key, value in attributes.items())
        lines.append(
            "\t".join(["chr1", "HAVANA", feature, str(start), str(end), ".", strand, ".", body])
        )
    return "\n".join(lines) + "\n"


ALL_MODULES = (ag, gr, ef)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_every_public_function_takes_keyword_only_verbose() -> None:
    for module in ALL_MODULES:
        assert module.__all__ == sorted(module.__all__)
        for name in module.__all__:
            function = getattr(module, name)
            parameter = inspect.signature(function).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is True


# ---------------------------------------------------------------- alphagenome


def test_normalize_variant_left_aligns_deletion_like_upstream(tmp_path: Path) -> None:
    genome = pytest.importorskip("alphagenome.data.genome")

    class Window:
        def extract(self, interval):  # reference "ACAAAG" starting at 0-based 0
            return "ACAAAG"[interval.start : interval.end]

    # chr1:4 AA>A sits in a run of three A's (positions 3-5); the leftmost form is 2 CA>C.
    upstream = genome.normalize_variant(genome.Variant.from_str("chr1:4:AA>A"), Window())
    result = ag.normalize_alphagenome_variant(
        "chr1:4:AA>A", reference_sequence="ACAAAG", reference_start=0
    )
    assert result["normalized"] == {
        "chromosome": "chr1",
        "position": 2,
        "reference": "CA",
        "alternate": "C",
    }
    assert result["normalized"]["position"] == upstream.position
    assert result["changed"] is True
    snv = ag.normalize_alphagenome_variant(
        "chr1:4:A>G", reference_sequence="ACAAAG", reference_start=0
    )
    assert snv["changed"] is False
    with pytest.raises(ConfigurationError):
        ag.normalize_alphagenome_variant("chr1:4:AA>A", reference_sequence="CA", reference_start=0)


def test_combine_intervals_matches_hand_calculation_and_upstream() -> None:
    genome = pytest.importorskip("alphagenome.data.genome")
    left = ["chr1:0-10", "chr1:20-30"]
    right = ["chr1:5-25"]
    inter = ag.combine_alphagenome_intervals(left, right, operation="intersect")
    assert [(item["start"], item["end"]) for item in inter["intervals"]] == [(5, 10), (20, 25)]
    union = ag.combine_alphagenome_intervals(left, right, operation="union")
    assert [(item["start"], item["end"]) for item in union["intervals"]] == [(0, 30)]
    upstream = list(
        genome.intersect_intervals(
            [genome.Interval.from_str(text) for text in left],
            [genome.Interval.from_str(text) for text in right],
        )
    )
    assert [(item.start, item.end) for item in upstream] == [(5, 10), (20, 25)]
    with pytest.raises(ConfigurationError):
        ag.combine_alphagenome_intervals(["chr1:0-10", "chr1:5-15"], right)


def test_extract_tss_points_follow_strand(tmp_path: Path) -> None:
    annotation = pytest.importorskip("alphagenome.data.gene_annotation")
    pandas = pytest.importorskip("pandas")
    gtf = _write(tmp_path / "small.gtf", _gtf_text())
    result = ag.extract_alphagenome_tss(gtf)
    assert [item["position"] for item in result["tss"]] == [1000, 1100, 5500]
    assert [item["strand"] for item in result["tss"]] == ["+", "+", "-"]
    frame = pandas.DataFrame(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [1000, 1100, 5000],
            "End": [1800, 2000, 5500],
            "Strand": ["+", "+", "-"],
            "Feature": ["transcript"] * 3,
            "gene_id": ["ENSG1.1", "ENSG1.1", "ENSG2.1"],
            "transcript_id": ["ENST1.1", "ENST2.1", "ENST3.1"],
            "gene_name": ["GENEA", "GENEA", "GENEB"],
        }
    )
    upstream = annotation.extract_tss(frame)
    assert upstream.Start.tolist() == [1000, 1100, 5500]


def test_transcript_filters_match_hand_calculation(tmp_path: Path) -> None:
    annotation = pytest.importorskip("alphagenome.data.gene_annotation")
    pandas = pytest.importorskip("pandas")
    gtf = _write(tmp_path / "small.gtf", _gtf_text())

    def ids(**options: object) -> list[str]:
        return ag.filter_alphagenome_transcripts(gtf, **options)["transcript_ids"]

    assert ids() == ["ENST1.1", "ENST2.1"]
    assert ids(method="longest") == ["ENST2.1", "ENST3.1"]
    assert ids(method="mane_select") == ["ENST1.1"]
    assert ids(method="support_level", support_levels=["1"]) == ["ENST1.1"]
    assert ids(method="transcript_type", transcript_types=["lncRNA"]) == ["ENST3.1"]
    frame = pandas.DataFrame(
        {
            "Feature": ["transcript"] * 3,
            "transcript_id": ["ENST1.1", "ENST2.1", "ENST3.1"],
            "transcript_type": ["protein_coding", "protein_coding", "lncRNA"],
        }
    )
    upstream = annotation.filter_protein_coding(frame)
    assert upstream["transcript_id"].tolist() == ["ENST1.1", "ENST2.1"]
    with pytest.raises(ConfigurationError):
        ag.filter_alphagenome_transcripts(gtf, method="transcript_type", transcript_types=["bogus"])


def test_upgrade_annotation_ids_maps_versions_and_marks_gaps() -> None:
    annotation = pytest.importorskip("alphagenome.data.gene_annotation")
    pandas = pytest.importorskip("pandas")
    result = ag.upgrade_alphagenome_annotation_ids(
        ["ENST00010.1", "ENST00011.1"], ["ENST00010.3", "ENST00012.1"]
    )
    assert result["mappings"] == [
        {"old_id": "ENST00010.1", "new_id": "ENST00010.3"},
        {"old_id": "ENST00011.1", "new_id": None},
    ]
    assert result["matched_count"] == 1
    upstream = annotation.upgrade_annotation_ids(
        pandas.Series(["ENST00010.1", "ENST00011.1"]),
        pandas.Series(["ENST00010.3", "ENST00012.1"]),
    )
    assert upstream.tolist()[0] == "ENST00010.3"
    patchless = ag.upgrade_alphagenome_annotation_ids(
        ["ENST00010", "ENST00011"], ["ENST00010.3", "ENST00011.2"], patchless=True
    )
    assert [item["new_id"] for item in patchless["mappings"]] == ["ENST00010.3", "ENST00011.2"]


def test_ism_variants_and_matrix_match_hand_calculation() -> None:
    ism = pytest.importorskip("alphagenome.interpretation.ism")
    genome = pytest.importorskip("alphagenome.data.genome")
    listed = ag.enumerate_alphagenome_ism_variants("chr1:100-103", "ACG")
    assert listed["count"] == 9
    assert listed["variants"][0] == {
        "chromosome": "chr1",
        "position": 101,
        "reference": "A",
        "alternate": "C",
    }
    upstream = ism.ism_variants(genome.Interval.from_str("chr1:100-103"), "ACG")
    assert [item.position for item in upstream] == [item["position"] for item in listed["variants"]]
    # Scores C=0.5, G=0.25, T=0.75 at one reference-A position: mean is 0.5, so A = -0.5.
    matrix = ag.build_alphagenome_ism_matrix("chr1:100-101", "A", [0.5, 0.25, 0.75])
    assert matrix["matrix"] == [[-0.5, 0.0, 0.0, 0.0]]
    full = ag.build_alphagenome_ism_matrix(
        "chr1:100-101", "A", [0.5, 0.25, 0.75], multiply_by_sequence=False
    )
    assert full["matrix"] == [[-0.5, 0.0, -0.25, 0.25]]
    with pytest.raises(ConfigurationError):
        ag.build_alphagenome_ism_matrix("chr1:100-101", "A", [0.5, 0.25])


def test_validate_sequence_length_accepts_supported_lengths_only() -> None:
    dna_client = pytest.importorskip("alphagenome.models.dna_client")
    result = ag.validate_alphagenome_sequence_length(16384)
    assert result["supported_lengths"] == [16384, 131072, 524288, 1048576]
    with pytest.raises(ValueError):
        dna_client.validate_sequence_length(1000)
    with pytest.raises(ConfigurationError):
        ag.validate_alphagenome_sequence_length(1000)


def test_fold_intervals_from_local_regions_file(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome.data.fold_intervals")
    regions = _write(
        tmp_path / "regions.bed",
        "chr1\t0\t100\tfold0\nchr2\t200\t300\tfold1\nchr3\t400\t500\tfold2\nchr1\t600\t700\tfold0\n",
    )
    valid = ag.get_alphagenome_fold_intervals(regions, model_version="fold_0", subset="valid")
    assert [(item["chromosome"], item["start"]) for item in valid["intervals"]] == [
        ("chr1", 0),
        ("chr1", 600),
    ]
    test_set = ag.get_alphagenome_fold_intervals(regions, model_version="fold_0", subset="test")
    assert [item["fold"] for item in test_set["intervals"]] == ["fold1"]
    train = ag.get_alphagenome_fold_intervals(regions, model_version="fold_0", subset="train")
    assert [item["fold"] for item in train["intervals"]] == ["fold2"]
    assert train["fold_names"] == ["fold2", "fold3", "fold4", "fold5", "fold6", "fold7"]


# ---------------------------------------------------------------- grelu


def test_convert_dna_encoding_round_trips_one_hot_and_indices() -> None:
    pytest.importorskip("grelu.sequence.format")
    torch = pytest.importorskip("torch")
    fmt = pytest.importorskip("grelu.sequence.format")
    one_hot = [
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]],
    ]
    strings = gr.convert_dna_encoding(one_hot, source="one_hot", target="strings")
    assert strings["values"] == ["ACGTN"]
    upstream = fmt.one_hot_to_strings(torch.tensor(one_hot, dtype=torch.float32))
    assert upstream == ["ACGTN"]
    expanded = gr.convert_dna_encoding([[0, 3]], source="indices", target="one_hot")
    assert expanded["values"] == [[[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]]]
    indices = gr.convert_dna_encoding(one_hot, source="one_hot", target="indices")
    assert indices["values"] == [[0, 1, 2, 3, 4]]
    with pytest.raises(ConfigurationError):
        gr.convert_dna_encoding(
            [[[1, 1], [1, 0], [0, 0], [0, 0]]], source="one_hot", target="strings"
        )
    with pytest.raises(ConfigurationError):
        gr.convert_dna_encoding([["A"]], source="strings", target="indices")


def test_check_dna_encoding_uses_grelu_checks() -> None:
    fmt = pytest.importorskip("grelu.sequence.format")
    assert fmt.check_string_dna("ACGTN") is True
    assert gr.check_dna_encoding(["ACGTN", "acgt"], encoding="strings")["valid"] is True
    assert gr.check_dna_encoding(["ACGX"], encoding="strings")["valid"] is False
    assert gr.check_dna_encoding([[0, 1, 4]], encoding="indices")["valid"] is True
    assert gr.check_dna_encoding([[0, 5]], encoding="indices")["valid"] is False
    good = [[[1, 0], [0, 1], [0, 0], [0, 0]]]
    double = [[[1, 1], [1, 0], [0, 0], [0, 0]]]
    assert gr.check_dna_encoding(good, encoding="one_hot")["valid"] is True
    assert gr.check_dna_encoding(double, encoding="one_hot")["valid"] is False
    with pytest.raises(ConfigurationError):
        gr.check_dna_encoding([1, 2], encoding="strings")


def test_gc_distribution_and_fraction_match_hand_counts() -> None:
    metrics = pytest.importorskip("grelu.sequence.metrics")
    sequences = ["ACGT", "GGGG", "AAAA"]
    hist = gr.gc_content_distribution(sequences, mode="distribution", binwidth=0.25)
    assert hist["bin_edges"] == [0.0, 0.25, 0.5, 0.75]
    assert hist["counts"] == [1.0, 0.0, 1.0, 1.0]
    assert metrics.gc_distribution(sequences, binwidth=0.25).tolist() == hist["counts"]
    fractions = gr.gc_content_distribution(sequences, mode="fraction")
    assert fractions["gc"] == [0.5, 1.0, 0.0]
    assert metrics.gc(sequences) == [0.5, 1.0, 0.0]
    with pytest.raises(ConfigurationError):
        gr.gc_content_distribution(sequences, binwidth=0.0)


def test_summarize_sequence_lengths_reports_equal_lengths() -> None:
    utils = pytest.importorskip("grelu.sequence.utils")
    equal = gr.summarize_sequence_lengths(["ACG", "TTT"])
    assert equal["lengths"] == [3, 3]
    assert equal["all_equal"] is True
    assert equal["unique_length"] == 3
    assert utils.get_lengths(["ACG", "TT"]) == [3, 2]
    ragged = gr.summarize_sequence_lengths(["ACG", "TT"])
    assert ragged["all_equal"] is False
    assert ragged["unique_length"] is None


def test_dna_sequence_differences_lists_positions_and_bases() -> None:
    mutate = pytest.importorskip("grelu.sequence.mutate")
    result = gr.dna_sequence_differences("ACGT", "ACCT")
    assert result["differences"] == [{"position": 2, "reference": "G", "alternate": "C"}]
    assert [int(item) for item in mutate.seq_differences("ACGT", "ACCT", verbose=False)] == [2]
    with pytest.raises(ConfigurationError):
        gr.dna_sequence_differences("ACG", "ACGT")


def test_extend_intervals_from_centers_uses_offset_summit() -> None:
    preprocess = pytest.importorskip("grelu.data.preprocess")
    pandas = pytest.importorskip("pandas")
    result = gr.extend_intervals_from_centers(
        [{"chrom": "chr1", "start": 100, "end": 200, "summit": 50}], seq_len=10
    )
    # start + summit - seq_len // 2 = 100 + 50 - 5
    assert result["intervals"] == [{"chrom": "chr1", "start": 145, "end": 155}]
    frame = pandas.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200], "summit": [50]})
    upstream = preprocess.extend_from_coord(frame, seq_len=10, center_col="summit")
    assert upstream["start"].tolist() == [145]
    with pytest.raises(ConfigurationError):
        gr.extend_intervals_from_centers(
            [{"chrom": "chr1", "start": 0, "end": 20, "summit": 2}], seq_len=10
        )


def test_clip_intervals_reports_rows_left_empty() -> None:
    preprocess = pytest.importorskip("grelu.data.preprocess")
    pandas = pytest.importorskip("pandas")
    rows = [
        {"chrom": "chr1", "start": 10, "end": 50},
        {"chrom": "chr2", "start": 200, "end": 300},
    ]
    result = gr.clip_intervals_to_bounds(rows, start=0, end=150)
    assert [(item["start"], item["end"]) for item in result["intervals"]] == [(10, 50), (200, 150)]
    assert result["invalid_count"] == 1
    frame = pandas.DataFrame(rows)
    clipped = preprocess.clip_intervals(frame, start=0, end=150)
    assert clipped["end"].tolist() == [50, 150]


def test_merge_intervals_by_group_spans_each_group() -> None:
    pytest.importorskip("grelu.data.preprocess")
    rows = [
        {"chrom": "chr1", "start": 0, "end": 10, "grp": "a"},
        {"chrom": "chr1", "start": 5, "end": 20, "grp": "a"},
        {"chrom": "chr2", "start": 0, "end": 5, "grp": "b"},
    ]
    result = gr.merge_intervals_by_group(rows, group_col="grp")
    assert result["intervals"] == [
        {"chrom": "chr1", "start": 0, "end": 20, "grp": "a"},
        {"chrom": "chr2", "start": 0, "end": 5, "grp": "b"},
    ]
    with pytest.raises(ConfigurationError):
        gr.merge_intervals_by_group(
            [
                {"chrom": "chr1", "start": 0, "end": 5, "grp": "x"},
                {"chrom": "chr2", "start": 0, "end": 5, "grp": "x"},
            ],
            group_col="grp",
        )


def test_overlap_filter_keeps_or_drops_overlapping_intervals() -> None:
    pytest.importorskip("grelu.data.preprocess")
    intervals = [
        {"chrom": "chr1", "start": 0, "end": 50},
        {"chrom": "chr1", "start": 100, "end": 150},
    ]
    reference = [{"chrom": "chr1", "start": 40, "end": 60}]
    kept = gr.filter_intervals_overlapping_reference(intervals, reference)
    assert kept["intervals"] == [{"chrom": "chr1", "start": 0, "end": 50}]
    dropped = gr.filter_intervals_overlapping_reference(intervals, reference, invert=True)
    assert dropped["intervals"] == [{"chrom": "chr1", "start": 100, "end": 150}]
    windowed = gr.filter_intervals_overlapping_reference(intervals, reference, window=50)
    assert windowed["count"] == 2


def test_split_by_chromosome_uses_named_lists() -> None:
    preprocess = pytest.importorskip("grelu.data.preprocess")
    pandas = pytest.importorskip("pandas")
    rows = [
        {"chrom": "chr1", "start": 0, "end": 50},
        {"chrom": "chr10", "start": 0, "end": 50},
        {"chrom": "chr11", "start": 0, "end": 50},
    ]
    result = gr.split_intervals_by_chromosome(rows)
    assert [item["chrom"] for item in result["train"]] == ["chr1"]
    assert [item["chrom"] for item in result["validation"]] == ["chr10"]
    assert [item["chrom"] for item in result["test"]] == ["chr11"]
    frame = pandas.DataFrame(rows)
    train, val, test = preprocess.split(
        frame, train_chroms=["chr1"], val_chroms=["chr10"], test_chroms=["chr11"]
    )
    assert (len(train), len(val), len(test)) == (1, 1, 1)
    with pytest.raises(ConfigurationError):
        gr.split_intervals_by_chromosome(rows, train_chroms=["chr10"], val_chroms=["chr10"])


def test_blacklist_filter_uses_bundled_hg38_regions() -> None:
    resources = pytest.importorskip("grelu.resources")
    intervals = [
        {"chrom": "chr1", "start": 0, "end": 500},
        {"chrom": "chr1", "start": 1000, "end": 1500},
        {"chrom": "chr2", "start": 5000, "end": 5500},
    ]
    result = gr.filter_intervals_blacklist(intervals, genome="hg38")
    regions = []
    for line in Path(resources.get_blacklist_file("hg38")).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and not line.startswith(("#", "track")):
            regions.append((fields[0], int(fields[1]), int(fields[2])))

    def overlaps(item: dict) -> bool:
        return any(
            chrom == item["chrom"] and start < item["end"] and item["start"] < end
            for chrom, start, end in regions
        )

    expected = [item for item in intervals if not overlaps(item)]
    assert result["intervals"] == expected
    assert result["removed_count"] == len(intervals) - len(expected)
    with pytest.raises(ConfigurationError):
        gr.filter_intervals_blacklist(intervals, genome="dm6")


def test_chromosome_filter_include_and_exclude() -> None:
    pytest.importorskip("grelu.data.preprocess")
    rows = [
        {"chrom": "chr1", "start": 0, "end": 10},
        {"chrom": "chr2", "start": 0, "end": 10},
    ]
    included = gr.filter_intervals_by_chromosome(rows, include=["chr2"])
    assert included["intervals"] == [{"chrom": "chr2", "start": 0, "end": 10}]
    excluded = gr.filter_intervals_by_chromosome(rows, exclude=["chr1"])
    assert excluded["intervals"] == [{"chrom": "chr2", "start": 0, "end": 10}]
    with pytest.raises(ConfigurationError):
        gr.filter_intervals_by_chromosome(rows)


def test_pwm_consensus_and_trim_match_hand_calculation() -> None:
    motifs = pytest.importorskip("grelu.interpret.motifs")
    np = pytest.importorskip("numpy")
    pwm = [
        [0.9, 0.05, 0.02],  # A row
        [0.05, 0.9, 0.03],  # C row
        [0.03, 0.02, 0.9],  # G row
        [0.02, 0.03, 0.05],  # T row
    ]
    assert gr.pwm_to_consensus(pwm)["consensus"] == "ACG"
    assert motifs.motifs_to_strings(np.array(pwm)) == "ACG"
    centred = np.array(pwm) - 0.25
    flanked = np.concatenate([np.zeros((4, 1)), centred, np.zeros((4, 1))], axis=1)
    trimmed = gr.trim_motif_pwm(flanked.tolist())
    assert (trimmed["start"], trimmed["end"]) == (1, 4)
    assert len(trimmed["pwm"][0]) == 3
    start, end = motifs.trim_pwm(flanked, trim_threshold=0.3, return_indices=True)
    assert (int(start), int(end)) == (1, 4)


def _motif_rows() -> list[list[float]]:
    # Strong ACG motif: column 0 favours A, column 1 favours C, column 2 favours G.
    base = [
        [0.9, 0.05, 0.03],
        [0.05, 0.9, 0.03],
        [0.03, 0.03, 0.9],
        [0.02, 0.02, 0.04],
    ]
    return base


def test_scan_motifs_reports_forward_and_reverse_hits() -> None:
    motifs = pytest.importorskip("grelu.interpret.motifs")
    np = pytest.importorskip("numpy")
    motif = np.array(_motif_rows())
    result = gr.scan_motifs_with_pwms({"s1": "TTTTACGTTT"}, {"M1": _motif_rows()}, pthresh=0.01)
    hits = [
        (item["start"], item["end"], item["strand"], item["matched_sequence"])
        for item in result["hits"]
    ]
    assert hits == [(4, 7, "+", "ACG"), (5, 8, "-", "CGT")]
    upstream = motifs.scan_sequences(["TTTTACGTTT"], {"M1": motif}, names=["M1"], pthresh=0.01)
    assert sorted(upstream["start"].tolist()) == [4, 5]
    listed = gr.scan_motifs_with_pwms(["TTTTACGTTT"], {"M1": _motif_rows()}, pthresh=0.01)
    assert listed["hits"][0]["sequence"] == "seq0"


def test_compare_motifs_scores_reference_above_mutated_allele() -> None:
    motifs = pytest.importorskip("grelu.interpret.motifs")
    np = pytest.importorskip("numpy")
    result = gr.compare_motif_pwm_alleles(
        "TTTTACGTTT", "TTTTATGTTT", {"M1": _motif_rows()}, pthresh=0.01
    )
    assert [(item["start"], item["strand"]) for item in result["sites"]] == [(4, "+"), (5, "-")]
    for site in result["sites"]:
        assert site["score_difference"] == pytest.approx(
            site["score_alternate"] - site["score_reference"]
        )
        assert site["score_reference"] > site["score_alternate"]
    frame = motifs.compare_motifs(
        "TTTTACGTTT", "TTTTATGTTT", {"M1": np.array(_motif_rows())}, names=["M1"], pthresh=0.01
    )
    upstream_diff = {
        (int(row.start), row.strand): round(float(row.fimo_score_diff), 10)
        for row in frame.itertuples(index=False)
    }
    ours = {
        (site["start"], site["strand"]): round(site["score_difference"], 10)
        for site in result["sites"]
    }
    assert ours == upstream_diff


def test_variant_filters_keep_snvs_and_drop_long_indels() -> None:
    filtering = pytest.importorskip("grelu.variant")
    pandas = pytest.importorskip("pandas")
    variants = [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"},
        {"chrom": "chr1", "pos": 20, "ref": "AC", "alt": "A"},
        {"chrom": "chr2", "pos": 30, "ref": "G", "alt": "GTT"},
    ]
    result = gr.filter_variants_by_length(variants)
    assert result["variants"] == [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"}]
    assert result["removed_count"] == 2
    relaxed = gr.filter_variants_by_length(variants, max_del_len=1, max_insert_len=2)
    assert relaxed["count"] == 3
    frame = pandas.DataFrame(variants)
    upstream = filtering.filter_variants(frame, max_insert_len=0, max_del_len=0, inplace=False)
    assert upstream["pos"].tolist() == [10]


def test_variant_intervals_follow_centering_rule() -> None:
    filtering = pytest.importorskip("grelu.variant")
    pandas = pytest.importorskip("pandas")
    variants = [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"},
        {"chrom": "chr2", "pos": 30, "ref": "G", "alt": "C"},
    ]
    result = gr.variants_to_sequence_intervals(variants, seq_len=10)
    assert [(item["start"], item["end"]) for item in result["intervals"]] == [(5, 15), (25, 35)]
    single = gr.variants_to_sequence_intervals(variants[:1], seq_len=1)
    assert (single["intervals"][0]["start"], single["intervals"][0]["end"]) == (9, 10)
    frame = pandas.DataFrame(variants)
    upstream = filtering.variants_to_intervals(frame, seq_len=10, inplace=False)
    assert upstream["start"].tolist() == [5, 25]


def test_read_fasta_returns_sequences_in_file_order(tmp_path: Path) -> None:
    pytest.importorskip("grelu.io.fasta")
    fasta = _write(tmp_path / "toy.fa", ">a\nACGTTA\n>b\nGGCC\n")
    result = gr.read_fasta_sequences(fasta)
    assert result["sequences"] == ["ACGTTA", "GGCC"]
    assert result["lengths"] == [6, 4]
    with pytest.raises(ConfigurationError):
        gr.read_fasta_sequences(fasta, max_sequences=1)


def test_sample_intervals_matches_seeded_numpy_choice() -> None:
    np = pytest.importorskip("numpy")
    pytest.importorskip("grelu.data.preprocess")
    rows = [{"chrom": "chr1", "start": i * 10, "end": i * 10 + 5} for i in range(6)]
    first = gr.sample_intervals(rows, n=3, seed=1)
    second = gr.sample_intervals(rows, n=3, seed=1)
    assert first["intervals"] == second["intervals"]
    chosen = sorted(np.random.RandomState(1).choice(6, 3, replace=False).tolist())
    assert [item["start"] for item in first["intervals"]] == [rows[i]["start"] for i in chosen]
    with pytest.raises(ConfigurationError):
        gr.sample_intervals(rows, n=7, seed=1)


# ---------------------------------------------------------------- enformer-pytorch


def test_encode_enformer_indices_matches_upstream() -> None:
    data = pytest.importorskip("enformer_pytorch.data")
    result = ef.encode_enformer_indices("ACGTN")
    assert result["indices"] == [0, 1, 2, 3, 4]
    assert data.str_to_seq_indices("ACGTN").tolist() == result["indices"]
    assert json.loads(json.dumps(result))["length"] == 5


def test_indices_to_one_hot_handles_n_and_padding() -> None:
    data = pytest.importorskip("enformer_pytorch.data")
    torch = pytest.importorskip("torch")
    result = ef.indices_to_one_hot_enformer([0, 1, 4, -1])
    assert result["one_hot"] == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25, 0.25],
    ]
    upstream = data.seq_indices_to_one_hot(torch.tensor([0, 1, 4, -1]))
    assert upstream.tolist() == result["one_hot"]
    with pytest.raises(ConfigurationError):
        ef.indices_to_one_hot_enformer([0, 7])


def test_reverse_complement_enformer_in_both_representations() -> None:
    data = pytest.importorskip("enformer_pytorch.data")
    indices = ef.reverse_complement_enformer("ACGTN")
    # reverse complement of ACGTN is NACGT
    assert indices["reverse_complement"] == [4, 0, 1, 2, 3]
    one_hot = ef.reverse_complement_enformer("ACGTN", representation="one_hot")
    expected = data.one_hot_reverse_complement(data.str_to_one_hot("ACGTN")).tolist()
    assert one_hot["reverse_complement"] == expected
    with pytest.raises(ConfigurationError):
        ef.reverse_complement_enformer("ACGT", representation="text")
