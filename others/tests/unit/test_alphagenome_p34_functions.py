"""Tests for the AlphaGenome functions of pack P34 (dnakit.predictions.alphagenome_p34_functions)."""

from __future__ import annotations

import inspect

import pytest

pytest.importorskip("alphagenome")
np = pytest.importorskip("numpy")
pytest.importorskip("pandas")

from alphagenome.data import genome

from dnakit.exceptions import ConfigurationError
from dnakit.predictions import alphagenome_p34_functions as a

Q = {"verbose": False}
TRACK = {
    "values": [[1, 2], [3, 4], [5, 6], [7, 8]],
    "names": ["a", "b"],
    "strands": ["+", "-"],
    "resolution": 2,
    "interval": "chr1:100-108:+",
}


def _span(result):
    return (result["start"], result["end"], result["strand"])


def test_adjust_interval_pad_and_boundary_shift_follow_strand():
    item = genome.Interval.from_str("chr1:100-200:-")
    assert (
        _span(a.alphagenome_adjust_interval("chr1:100-200:-", start=10, **Q)["result"])
        == _span(a._interval_dict(item.pad(10, 0)))
        == (100, 210, "-")
    )
    shifted = a.alphagenome_adjust_interval(
        {"chromosome": "chr1", "start": 100, "end": 200, "strand": "-"},
        operation="boundary_shift",
        start=10,
        end=5,
        **Q,
    )["result"]
    assert _span(shifted) == (95, 190, "-")
    plain = a.alphagenome_adjust_interval("chr1:100-200:-", start=10, use_strand=False, **Q)[
        "result"
    ]
    assert _span(plain) == (90, 200, "-")
    with pytest.raises(ConfigurationError):
        a.alphagenome_adjust_interval("chr1:100-200", start=-200, **Q)
    with pytest.raises(ConfigurationError):
        a.alphagenome_adjust_interval("chr1:100-200", operation="grow", **Q)


def test_interval_center_reference_strand_and_relations():
    assert a.alphagenome_interval_center("chr1:100-201:+", **Q)["center"] == 151
    even = "chr1:100-200:-"
    assert (
        a.alphagenome_interval_center(even, **Q)["center"]
        == genome.Interval.from_str(even).center()
    )
    assert (
        a.alphagenome_interval_within_reference("chr1:-5-100", reference_length=1000, **Q)[
            "within_reference"
        ]
        is False
    )
    assert (
        a.alphagenome_interval_within_reference("chr1:0-1000", reference_length=1000, **Q)[
            "within_reference"
        ]
        is True
    )
    assert a.alphagenome_interval_strand("chr1:1-10:+", **Q)["result"]["strand"] == "-"
    assert a.alphagenome_interval_strand("chr1:1-10:-", mode="drop", **Q)["result"]["strand"] == "."
    assert a.alphagenome_interval_relation("chr1:1-10", "chr1:9-20", **Q)["result"] is True
    assert a.alphagenome_interval_relation("chr1:1-10", "chr1:10-20", **Q)["result"] is False
    assert a.alphagenome_interval_relation("chr2:1-10", "chr1:2-5", **Q)["result"] is False
    assert (
        a.alphagenome_interval_relation("chr1:1-30", "chr1:5-20", relation="contains", **Q)[
            "result"
        ]
        is True
    )
    assert (
        a.alphagenome_interval_relation("chr1:6-30", "chr1:5-20", relation="contains", **Q)[
            "result"
        ]
        is False
    )


def test_interval_coverage_variants_match_alphagenome():
    window = genome.Interval.from_str("chr1:0-10:-")
    others = [genome.Interval.from_str(x) for x in ("chr1:2-5:+", "chr1:4-8:-", "chr2:0-10:+")]
    texts = ["chr1:2-5:+", "chr1:4-8:-", "chr2:0-10:+"]
    counts = a.alphagenome_interval_coverage("chr1:0-10:-", texts, **Q)
    assert counts["values"] == window.coverage(others).tolist() == [0, 0, 1, 1, 2, 1, 1, 1, 0, 0]
    binned = a.alphagenome_interval_coverage("chr1:0-10:-", texts, bin_size=5, **Q)
    assert binned["values"] == [4, 3] and binned["shape"] == [2]
    mask = a.alphagenome_interval_coverage("chr1:0-10:-", texts, mask=True, **Q)
    assert mask["values"] == window.binary_mask(others).tolist()
    stranded = a.alphagenome_interval_coverage("chr1:0-10:-", texts, stranded=True, **Q)
    assert stranded["values"] == window.coverage_stranded(others).tolist() and stranded[
        "shape"
    ] == [10, 2]
    both = a.alphagenome_interval_coverage("chr1:0-10:-", texts, stranded=True, mask=True, **Q)
    assert both["values"] == window.binary_mask_stranded(others).tolist()
    with pytest.raises(ConfigurationError):
        a.alphagenome_interval_coverage("chr1:0-10", texts, bin_size=3, **Q)
    empty = a.alphagenome_interval_coverage("chr1:0-4", [], **Q)
    assert empty["values"] == [0, 0, 0, 0]


def test_overlap_ranges_relative_and_clipped():
    out = a.alphagenome_overlap_ranges("chr1:3-10", ["chr1:2-5:+", "chr1:4-8:-", "chr1:20-30"], **Q)
    assert out["ranges"] == [[0, 2], [1, 5]]
    assert a.alphagenome_overlap_ranges("chr1:3-10", [], **Q)["ranges"] == []


def test_variant_overlap_and_split():
    assert a.alphagenome_variant_overlap("chr1:3:AC>T", "chr1:2-4", **Q)["overlaps"] is True
    assert a.alphagenome_variant_overlap("chr1:3:A>T", "chr1:3-6", **Q)["overlaps"] is False
    ins = "chr1:3:A>TTTT"
    assert a.alphagenome_variant_overlap(ins, "chr1:5-9", **Q)["overlaps"] is False
    assert (
        a.alphagenome_variant_overlap(ins, "chr1:5-9", allele="alternate", **Q)["overlaps"] is True
    )
    split = a.alphagenome_split_variant("chr1:3:AC>TG", anchor=3, **Q)
    assert split["upstream"]["text"] == "chr1:3:A>T" and split["downstream"]["text"] == "chr1:4:C>G"
    before = a.alphagenome_split_variant("chr1:3:AC>TG", anchor=2, **Q)
    assert before["upstream"] is None and before["downstream"]["text"] == "chr1:3:AC>TG"
    after = a.alphagenome_split_variant(
        {"chromosome": "chr1", "position": 3, "reference_bases": "AC", "alternate_bases": "TG"},
        anchor=9,
        **Q,
    )
    assert after["downstream"] is None
    with pytest.raises(ConfigurationError):
        a.alphagenome_split_variant("chr1:3:AX>T", anchor=3, **Q)


def test_junction_regions_acceptor_donor_and_dinucleotides():
    di = a.alphagenome_junction_regions("chr1:100-200:+", **Q)
    assert [(r["start"], r["end"]) for r in di["regions"]] == [(100, 102), (198, 200)]
    assert (di["donor"], di["acceptor"]) == (100, 200)
    minus = a.alphagenome_junction_regions(
        "chr1:100-200:-", region="acceptor", overhang=[5, 5], **Q
    )
    junction = genome.Junction("chr1", 100, 200, "-")
    expected = junction.acceptor_region((5, 5))
    assert (
        (minus["regions"][0]["start"], minus["regions"][0]["end"])
        == (expected.start, expected.end)
        == (95, 105)
    )
    donor = a.alphagenome_junction_regions("chr1:100-200:+", region="donor", overhang=[2, 3], **Q)
    assert (donor["regions"][0]["start"], donor["regions"][0]["end"]) == (98, 103)
    with pytest.raises(ConfigurationError):
        a.alphagenome_junction_regions("chr1:100-200:.", **Q)


TRANSCRIPT = {
    "chromosome": "chr1",
    "strand": "+",
    "exons": [[100, 200], [300, 400]],
    "cds": [[150, 200, 0], [300, 350, 2]],
    "stop_codon": [[350, 353, 0]],
    "transcript_id": "T1",
}


def test_cds_offset_both_strands_and_outside():
    assert a.alphagenome_transcript_cds_offset(TRANSCRIPT, position=150, **Q)["offset"] == 0
    assert a.alphagenome_transcript_cds_offset(TRANSCRIPT, position=310, **Q)["offset"] == 60
    assert a.alphagenome_transcript_cds_offset(TRANSCRIPT, position=250, **Q)["offset"] is None
    minus = {**TRANSCRIPT, "strand": "-"}
    assert a.alphagenome_transcript_cds_offset(minus, position=352, **Q)["offset"] == 0
    assert a.alphagenome_transcript_cds_offset(minus, position=310, **Q)["offset"] == 42
    with pytest.raises(ConfigurationError):
        a.alphagenome_transcript_cds_offset({"strand": "+"}, position=1, **Q)


def test_fix_truncation_matches_alphagenome_and_is_noop_for_tiny_cds():
    from alphagenome.data import transcript as tx

    truncated = {
        "chromosome": "chr1",
        "strand": "+",
        "exons": [[100, 200], [300, 400]],
        "cds": [[151, 200, 1], [300, 350, 0]],
        "transcript_id": "T1",
    }
    out = a.alphagenome_fix_transcript_truncation(truncated, **Q)["result"]
    reference = tx.Transcript.fix_truncation(a._transcript(truncated, "t"))
    assert out["cds"] == [[x.start, x.end, x.info["frame"]] for x in reference.cds]
    assert out["start_codon"] == [[x.start, x.end] for x in reference.start_codon]
    assert out["truncation_fixed"] is True
    tiny = {"chromosome": "chr1", "strand": "+", "exons": [[100, 200]], "cds": [[100, 104, 0]]}
    assert a.alphagenome_fix_transcript_truncation(tiny, **Q)["result"]["cds"] == [[100, 104, 0]]
    with pytest.raises(ConfigurationError):
        a.alphagenome_fix_transcript_truncation({**truncated, "cds": [[151, 200], [300, 350]]}, **Q)


GTF = (
    'chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding";\n'
    'chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding";\n'
    'chr1\tt\tCDS\t151\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\tCDS\t301\t350\t.\t+\t2\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\tstart_codon\t151\t153\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\tstop_codon\t351\t353\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t501\t600\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
)


def test_transcripts_from_gtf_reads_selected_transcripts(tmp_path):
    path = tmp_path / "t.gtf"
    path.write_text(GTF)
    out = a.alphagenome_transcripts_from_gtf(str(path), transcript_ids=["T1", "T2", "T9"], **Q)
    assert out["missing"] == ["T9"] and out["transcript_count"] == 2
    first = out["transcripts"][0]
    assert first["exons"] == [[100, 200], [300, 400]] and first["cds"] == [
        [150, 200, 0],
        [300, 350, 2],
    ]
    assert first["start_codon"] == [[150, 153]] and first["stop_codon"] == [[350, 353, 0]]
    assert out["transcripts"][1]["strand"] == "-" and out["transcripts"][1]["is_coding"] is False
    fixed = a.alphagenome_transcripts_from_gtf(
        str(path), transcript_ids=["T1"], fix_truncation=True, **Q
    )
    assert fixed["transcripts"][0]["truncation_fixed"] is True
    import gzip

    zipped = tmp_path / "t.gtf.gz"
    zipped.write_bytes(gzip.compress(GTF.encode()))
    assert (
        a.alphagenome_transcripts_from_gtf(str(zipped), transcript_ids=["T1"], **Q)[
            "transcript_count"
        ]
        == 1
    )
    with pytest.raises(ConfigurationError):
        a.alphagenome_transcripts_from_gtf(str(path), transcript_ids=[], **Q)


def test_slice_bin_index_resize_resolution_and_reverse_complement():
    by_interval = a.alphagenome_slice_track_data(TRACK, interval="chr1:102-106", **Q)["result"]
    by_position = a.alphagenome_slice_track_data(TRACK, start=2, end=6, **Q)["result"]
    assert by_interval["values"] == by_position["values"] == [[3.0, 4.0], [5.0, 6.0]]
    assert by_interval["interval"]["start"] == 102
    widened = a.alphagenome_slice_track_data(
        TRACK, interval="chr1:103-105", match_resolution=True, **Q
    )["result"]
    assert widened["shape"] == [1, 2] or widened["shape"] == [2, 2]
    with pytest.raises(ConfigurationError):
        a.alphagenome_slice_track_data(TRACK, interval="chr1:90-106", **Q)
    with pytest.raises(ConfigurationError):
        a.alphagenome_slice_track_data(TRACK, interval="chr1:102-106", start=1, end=2, **Q)
    assert a.alphagenome_track_bin_index(TRACK, relative_position=5, **Q)["bin_index"] == 2

    cropped = a.alphagenome_resize_track_data(TRACK, width=4, **Q)["result"]
    assert cropped["values"] == [[3.0, 4.0], [5.0, 6.0]]
    padded = a.alphagenome_resize_track_data(TRACK, width=12, **Q)["result"]
    assert padded["shape"] == [6, 2] and padded["values"][0] == [0.0, 0.0]
    with pytest.raises(ConfigurationError):
        a.alphagenome_resize_track_data(TRACK, width=13, **Q)

    down = a.alphagenome_change_track_resolution(TRACK, resolution=4, **Q)
    assert down["direction"] == "downsample" and down["result"]["values"] == [
        [4.0, 6.0],
        [12.0, 14.0],
    ]
    top = a.alphagenome_change_track_resolution(TRACK, resolution=4, aggregation="max", **Q)
    assert top["result"]["values"] == [[3.0, 4.0], [7.0, 8.0]]
    up = a.alphagenome_change_track_resolution(TRACK, resolution=1, **Q)
    assert up["direction"] == "upsample" and up["result"]["shape"] == [8, 2]
    assert up["result"]["values"][0] == [0.5, 1.0]
    with pytest.raises(ConfigurationError):
        a.alphagenome_change_track_resolution(TRACK, resolution=3, **Q)

    paired = {**TRACK, "names": ["a", "a"]}
    flipped = a.alphagenome_reverse_complement_track_data(paired, **Q)["result"]
    assert flipped["values"] == [[8.0, 7.0], [6.0, 5.0], [4.0, 3.0], [2.0, 1.0]]
    assert flipped["strands"] == ["-", "+"] and flipped["interval"]["strand"] == "-"
    with pytest.raises(ConfigurationError):
        a.alphagenome_reverse_complement_track_data(TRACK, **Q)


def test_junction_normalisation_and_interval_filter():
    data = {
        "junctions": ["chr1:100-200:+", "chr1:500-600:-"],
        "values": [[1, 3], [2, 6]],
        "names": ["x", "y"],
    }
    out = a.alphagenome_normalize_junction_counts(data, total_k=24.0, **Q)["result"]
    assert out["values"] == [[2.0, 6.0], [4.0, 12.0]]
    kept = a.alphagenome_junctions_in_interval(data, interval="chr1:150-160", **Q)
    assert kept["result_count"] == 1 and kept["result"]["junctions"][0]["start"] == 100
    assert kept["result"]["values"] == [[1.0, 3.0]]
    assert a.alphagenome_junctions_in_interval(data, interval="chr3:1-9", **Q)["result_count"] == 0
    with pytest.raises(ConfigurationError):
        a.alphagenome_junctions_in_interval({**data, "values": [[1, 2]]}, interval="chr1:1-9", **Q)
    with pytest.raises(ConfigurationError):
        a.alphagenome_junctions_in_interval(
            {**data, "junctions": ["chr1:1-9:.", "chr1:5-9:+"]}, interval="chr1:1-9", **Q
        )


def test_extract_fasta_interval_strand_padding_and_errors(tmp_path):
    pytest.importorskip("pyfaidx")
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGTACGTACGTAAAC\n>chr2\nGGGGCCCC\n")
    import pyfaidx

    pyfaidx.Faidx(str(fasta))
    assert (
        a.alphagenome_extract_fasta_interval(str(fasta), interval="chr1:2-8:+", **Q)["sequence"]
        == "GTACGT"
    )
    assert (
        a.alphagenome_extract_fasta_interval(str(fasta), interval="chr1:2-8:-", **Q)["sequence"]
        == "ACGTAC"
    )
    assert (
        a.alphagenome_extract_fasta_interval(str(fasta), interval="chr2:4-12", **Q)["sequence"]
        == "CCCCNNNN"
    )
    assert (
        a.alphagenome_extract_fasta_interval(str(fasta), interval="chr2:-2-2", **Q)["sequence"]
        == "NNGG"
    )
    with pytest.raises(ConfigurationError):
        a.alphagenome_extract_fasta_interval(str(fasta), interval="chr9:0-5", **Q)


def test_output_and_scorer_metadata_use_the_client_without_network(monkeypatch):
    import pandas
    from alphagenome.atlas import atlas
    from alphagenome.models import dna_client
    from alphagenome.models.v1 import dna_output

    metadata = dna_output.OutputMetadata(
        atac=pandas.DataFrame({"name": ["t1", "t2", "t3"], "strand": ["+", "-", "."]}),
        cage=pandas.DataFrame({"name": ["c1"], "strand": ["+"]}),
    )
    seen = {}

    class FakeDna:
        def output_metadata(self, organism):
            seen["organism"] = organism
            return metadata

    class FakeAtlas:
        def scorer_metadata(self):
            return {
                "S1": atlas.ScorerMetadata(
                    name="S1", is_signed=True, track_metadata=pandas.DataFrame({"name": ["z"]})
                )
            }

    monkeypatch.setattr(dna_client, "create", lambda key, **kw: seen.update(key=key) or FakeDna())
    monkeypatch.setattr(atlas, "create", lambda key, **kw: FakeAtlas())
    out = a.alphagenome_output_metadata("secret", organism="mus_musculus", max_rows=2, **Q)
    assert seen["organism"] is dna_client.Organism.MUS_MUSCULUS and seen["key"] == "secret"
    assert out["outputs"]["atac"]["row_count"] == 3 and len(out["outputs"]["atac"]["rows"]) == 2
    assert out["outputs"]["cage"]["columns"] == ["name", "strand"] and "dnase" not in out["outputs"]
    assert "secret" not in repr(out)
    scorers = a.alphagenome_scorer_metadata("secret", **Q)
    assert scorers["scorer_count"] == 1 and scorers["scorers"]["S1"]["is_signed"] is True
    with pytest.raises(ConfigurationError):
        a.alphagenome_output_metadata("secret", organism="zebrafish", **Q)
    with pytest.raises(ConfigurationError):
        a.alphagenome_output_metadata("  ", **Q)


def _is_png(path):
    with open(path, "rb") as handle:
        return handle.read(8) == b"\x89PNG\r\n\x1a\n"


def test_plots_write_png_files(tmp_path):
    pytest.importorskip("matplotlib")
    rng = np.random.RandomState(0)
    contact = {
        "values": rng.randn(8, 8, 1).tolist(),
        "names": ["c"],
        "resolution": 4,
        "interval": "chr1:0-32:.",
    }
    out = a.alphagenome_plot_contact_map_difference(contact, str(tmp_path / "cm.png"), **Q)
    assert _is_png(out["image_path"]) and out["vmin"] == -1.0
    with pytest.raises(ConfigurationError):
        a.alphagenome_plot_contact_map_difference(TRACK, str(tmp_path / "bad.png"), **Q)
    with pytest.raises(ConfigurationError):
        a.alphagenome_plot_contact_map_difference(contact, str(tmp_path / "bad.txt"), **Q)

    ref = {
        "values": [[1, 2], [3, 4], [5, 6], [7, 8]],
        "names": ["a", "b"],
        "resolution": 2,
        "interval": "chr1:100-108:.",
    }
    alt = {**ref, "values": [[2, 2], [3, 5], [5, 6], [1, 8]]}
    overlaid = a.alphagenome_plot_overlaid_tracks(
        {"REF": ref, "ALT": alt}, str(tmp_path / "ov.png"), **Q
    )
    assert _is_png(overlaid["image_path"]) and overlaid["labels"] == ["REF", "ALT"]
    with pytest.raises(ConfigurationError):
        a.alphagenome_plot_overlaid_tracks({}, str(tmp_path / "x.png"), **Q)

    both = a.alphagenome_plot_annotations(
        "chr1:100-108",
        str(tmp_path / "an.png"),
        intervals=["chr1:101-104"],
        variants=["chr1:105:A>T"],
        track=ref,
        **Q,
    )
    assert (
        _is_png(both["image_path"])
        and len(both["intervals"]) == 1
        and both["variants"][0]["position"] == 105
    )
    only = a.alphagenome_plot_annotations(
        "chr1:100-108", str(tmp_path / "an2.png"), variants=["chr1:105:A>T"], **Q
    )
    assert _is_png(only["image_path"])
    with pytest.raises(ConfigurationError):
        a.alphagenome_plot_annotations("chr1:100-108", str(tmp_path / "none.png"), **Q)


def test_all_public_functions_have_keyword_only_verbose_and_docstrings():
    for name in a.__all__:
        function = getattr(a, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
