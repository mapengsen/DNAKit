"""Tests for dnakit.sequence.biopython_p18_functions (sequence formats and Bio.Seq helpers).

The module is imported directly. Inputs are generated in tmp_path or taken from the small
Biopython test files in _biopython_p18_data.py. Nothing uses the network.
"""

from __future__ import annotations

import pytest

pytest.importorskip("Bio")

import _biopython_p18_data as data
from Bio import SeqIO
from Bio.Seq import Seq

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import biopython_p18_functions as m

QUIET = {"verbose": False}


def _records() -> list[dict[str, object]]:
    return [
        {"id": "a", "sequence": "ACGTACGTAC", "description": "first", "qualities": [30] * 10},
        {"id": "b", "sequence": "GGTTAACCGG", "qualities": [20] * 10},
    ]


@pytest.mark.parametrize("fmt", ["imgt", "pir", "seqxml", "tab"])
def test_write_then_read_multi_record_formats(tmp_path, fmt):
    out = tmp_path / f"x.{fmt}"
    written = m.write_sequence_file(_records(), out, file_format=fmt, **QUIET)
    assert written["record_count"] == 2
    back = m.read_sequence_file(out, file_format=fmt, **QUIET)
    assert [(r["id"], r["sequence"]) for r in back["records"]] == [
        ("a", "ACGTACGTAC"),
        ("b", "GGTTAACCGG"),
    ]


@pytest.mark.parametrize("fmt", ["nib", "xdna"])
def test_write_then_read_single_record_formats(tmp_path, fmt):
    out = tmp_path / f"x.{fmt}"
    m.write_sequence_file(_records()[:1], out, file_format=fmt, **QUIET)
    back = m.read_sequence_file(out, file_format=fmt, **QUIET)
    assert back["records"][0]["sequence"] == "ACGTACGTAC"
    with pytest.raises(ConfigurationError, match="exactly one"):
        m.write_sequence_file(_records(), tmp_path / f"y.{fmt}", file_format=fmt, **QUIET)


def test_write_phd_matches_biopython(tmp_path):
    out = tmp_path / "x.phd"
    m.write_sequence_file(_records(), out, file_format="phd", **QUIET)
    records = list(SeqIO.parse(str(out), "phd"))
    assert [str(r.seq) for r in records] == ["ACGTACGTAC", "GGTTAACCGG"]
    assert records[0].letter_annotations["phred_quality"] == [30] * 10
    with pytest.raises(ConfigurationError, match="qualities"):
        m.write_sequence_file(
            [{"id": "a", "sequence": "ACGT"}], tmp_path / "z.phd", file_format="phd", **QUIET
        )


def test_write_refuses_existing_file_and_bad_input(tmp_path):
    out = tmp_path / "x.tab"
    out.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        m.write_sequence_file(_records(), out, file_format="tab", **QUIET)
    assert out.read_text(encoding="utf-8") == "keep"
    with pytest.raises(ConfigurationError):
        m.write_sequence_file(_records(), tmp_path / "n.tab", file_format="fasta", **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_sequence_file(
            [{"id": "a b", "sequence": "ACGT"}], tmp_path / "n.tab", file_format="tab", **QUIET
        )


def test_read_qual_has_undefined_sequence_and_qualities(tmp_path):
    path = tmp_path / "x.qual"
    path.write_text(data.QUAL, encoding="utf-8")
    result = m.read_sequence_file(path, file_format="qual", **QUIET)
    first = result["records"][0]
    assert first["sequence"] is None
    assert len(first["letter_annotations"]["phred_quality"]) == first["length"] == 25


def test_read_gfa_gck_and_ig(tmp_path):
    gfa1 = tmp_path / "x.gfa"
    gfa1.write_text(
        "H\tVN:Z:1.0\nS\tseg1\tACGTACGT\nS\tseg2\tGGCC\nL\tseg1\t+\tseg2\t+\t0M\n", encoding="utf-8"
    )
    result = m.read_sequence_file(gfa1, file_format="gfa1", **QUIET)
    assert [(r["id"], r["sequence"]) for r in result["records"]] == [
        ("seg1", "ACGTACGT"),
        ("seg2", "GGCC"),
    ]
    gfa2 = tmp_path / "y.gfa"
    gfa2.write_text(data.GFA2, encoding="utf-8")
    assert (
        m.read_sequence_file(gfa2, file_format="gfa2", **QUIET)["records"][0]["sequence"] == "AAA"
    )
    gck = tmp_path / "x.gck"
    gck.write_bytes(data.GCK)
    sample = m.read_sequence_file(gck, file_format="gck", **QUIET)["records"][0]
    assert sample["sequence"] == "ACGTACGTACGT"
    assert sample["feature_count"] == 2
    ig = tmp_path / "x.ig"
    ig.write_text(";\n;\nseqA\nACGTACGTACGT1\n;\nseqB\nGGGTTTAAACCC1\n", encoding="utf-8")
    result = m.read_sequence_file(ig, file_format="ig", **QUIET)
    assert [r["sequence"] for r in result["records"]] == ["ACGTACGTACGT", "GGGTTTAAACCC"]


def test_read_nib_xdna_and_seqxml_samples(tmp_path):
    nib = tmp_path / "x.nib"
    nib.write_bytes(data.NIB)
    assert m.read_sequence_file(nib, file_format="nib", **QUIET)["records"][0]["length"] > 0
    xdna = tmp_path / "x.xdna"
    xdna.write_bytes(data.XDNA)
    record = m.read_sequence_file(xdna, file_format="xdna", **QUIET)["records"][0]
    assert record["length"] == 1000 and record["feature_count"] == 2
    seqxml = tmp_path / "x.xml"
    seqxml.write_text(data.SEQXML, encoding="utf-8")
    assert m.read_sequence_file(seqxml, file_format="seqxml", **QUIET)["record_count"] == 4


def test_read_cds_feature_records(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq("ATGGCCATTGTAATGGGCCGCTGA"), id="TEST1", name="TEST1", description="demo"
    )
    record.annotations["molecule_type"] = "DNA"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 24, strand=1),
            type="CDS",
            qualifiers={"gene": ["g1"], "translation": ["MAIVMGR"], "protein_id": ["P1.1"]},
        )
    )
    for fmt, cds in (("genbank", "genbank-cds"), ("embl", "embl-cds")):
        path = tmp_path / f"x.{fmt}"
        SeqIO.write(record, str(path), fmt)
        result = m.read_sequence_file(path, file_format=cds, **QUIET)
        assert result["record_count"] == 1
        assert result["records"][0]["sequence"] == "MAIVMGR"


def test_read_sequence_file_limits_and_errors(tmp_path):
    path = tmp_path / "x.tab"
    path.write_text("a\tACGT\nb\tTTTT\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="max_records"):
        m.read_sequence_file(path, file_format="tab", max_records=1, **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_sequence_file(path, file_format="fasta", **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_sequence_file(tmp_path / "missing.tab", file_format="tab", **QUIET)
    bad = tmp_path / "bad.abi"
    bad.write_bytes(b"not an abi file")
    with pytest.raises(ConfigurationError):
        m.read_sequence_file(bad, file_format="abi-trim", **QUIET)


def test_convert_sequence_format_between_text_formats(tmp_path):
    source = tmp_path / "x.tab"
    source.write_text("a\tACGTACGT\nb\tGGTTAACC\n", encoding="utf-8")
    fasta = tmp_path / "x.fasta"
    result = m.convert_sequence_format(
        source, fasta, source_format="tab", target_format="fasta", **QUIET
    )
    assert result["record_count"] == 2
    assert fasta.read_text(encoding="utf-8").splitlines()[:2] == [">a", "ACGTACGT"]
    pir = tmp_path / "x.pir"
    m.convert_sequence_format(fasta, pir, source_format="fasta", target_format="pir", **QUIET)
    assert [r.id for r in SeqIO.parse(str(pir), "pir")] == ["a", "b"]
    with pytest.raises(ConfigurationError, match="More than one sequence"):
        m.convert_sequence_format(
            fasta, tmp_path / "x.nib", source_format="fasta", target_format="nib", **QUIET
        )
    assert not (tmp_path / "x.nib").exists()


def test_sff_manifest_and_roundtrip(tmp_path):
    sff = tmp_path / "x.sff"
    sff.write_bytes(data.SFF)
    manifest = m.read_sff_manifest(sff, **QUIET)
    assert manifest["has_manifest"] and manifest["xml"].lstrip().startswith("<manifest>")
    copy = tmp_path / "copy.sff"
    result = m.convert_sequence_format(sff, copy, source_format="sff", target_format="sff", **QUIET)
    assert result["record_count"] == 10 and result["roche_manifest_copied"]
    assert m.read_sff_manifest(copy, **QUIET)["has_manifest"]
    original = [str(r.seq) for r in SeqIO.parse(str(sff), "sff")]
    assert [str(r.seq) for r in SeqIO.parse(str(copy), "sff")] == original
    fastq = tmp_path / "x.fastq"
    m.convert_sequence_format(sff, fastq, source_format="sff-trim", target_format="fastq", **QUIET)
    assert len(list(SeqIO.parse(str(fastq), "fastq"))) == 10


def test_sff_manifest_missing_and_invalid(tmp_path):
    bad = tmp_path / "bad.sff"
    bad.write_bytes(b"nothing here at all")
    with pytest.raises(ConfigurationError):
        m.read_sff_manifest(bad, **QUIET)


def test_count_subsequence_matches_biopython():
    sequence = "AAAAGAAAAGAAAA"
    assert m.count_subsequence(sequence, "AA", **QUIET)["count"] == Seq(sequence).count("AA") == 6
    overlapping = m.count_subsequence(sequence, "AA", overlap=True, **QUIET)
    assert overlapping["count"] == Seq(sequence).count_overlap("AA") == 9
    window = m.count_subsequence(sequence, "AA", overlap=True, start=2, end=8, **QUIET)
    assert window["count"] == Seq(sequence).count_overlap("AA", 2, 8) == 3
    assert m.count_subsequence(sequence, "aa", **QUIET)["count"] == 0
    with pytest.raises(ConfigurationError):
        m.count_subsequence("", "A", **QUIET)
    with pytest.raises(ConfigurationError):
        m.count_subsequence(sequence, "A", overlap="yes", **QUIET)  # type: ignore[arg-type]


def test_search_subsequences_positions():
    result = m.search_subsequences("AAAAGAAAAGAAAA", ["AAG", "GA"], **QUIET)
    assert [(h["position"], h["subsequence"]) for h in result["hits"]] == [
        (2, "AAG"),
        (4, "GA"),
        (7, "AAG"),
        (9, "GA"),
    ]
    assert m.search_subsequences("ACGT", ["TT"], **QUIET)["hit_count"] == 0
    with pytest.raises(ConfigurationError):
        m.search_subsequences("ACGT", [], **QUIET)


def test_split_sequence_left_right_and_limit():
    sequence = "AAAAGAAAAGAAAA"
    assert m.split_sequence(sequence, "G", **QUIET)["pieces"] == ["AAAA", "AAAA", "AAAA"]
    assert m.split_sequence(sequence, "G", maxsplit=1, **QUIET)["pieces"] == ["AAAA", "AAAAGAAAA"]
    right = m.split_sequence(sequence, "G", maxsplit=1, from_right=True, **QUIET)
    assert right["pieces"] == ["AAAAGAAAA", "AAAA"]
    assert m.split_sequence("AC GT", **QUIET)["pieces"] == ["AC", "GT"]


@pytest.mark.parametrize(
    ("mode", "pattern", "expected"),
    [
        ("strip", "A", "GAAAAG"),
        ("lstrip", "A", "GAAAAGAAAA"),
        ("rstrip", "A", "AAAAGAAAAG"),
        ("removeprefix", "AAAA", "GAAAAGAAAA"),
        ("removesuffix", "AAAA", "AAAAGAAAAG"),
        ("removesuffix", "CC", "AAAAGAAAAGAAAA"),
    ],
)
def test_trim_sequence_modes(mode, pattern, expected):
    assert (
        m.trim_sequence("AAAAGAAAAGAAAA", mode=mode, pattern=pattern, **QUIET)["trimmed"]
        == expected
    )


def test_trim_sequence_requires_pattern_for_affix():
    with pytest.raises(ConfigurationError):
        m.trim_sequence("ACGT", mode="removeprefix", **QUIET)
    assert m.trim_sequence("  ACGT ", **QUIET)["trimmed"] == "ACGT"


def test_parse_feature_location_compound_fuzzy_and_reference():
    joined = m.parse_feature_location("complement(join(1..10,20..30))", **QUIET)
    assert (joined["start"]["value"], joined["end"]["value"], joined["strand"]) == (0, 30, -1)
    assert joined["operator"] == "join" and joined["part_count"] == 2
    assert [p["start"]["value"] for p in joined["parts"]] == [19, 0]
    fuzzy = m.parse_feature_location("<1..>200", **QUIET)
    assert fuzzy["start"]["type"] == "BeforePosition" and fuzzy["end"]["type"] == "AfterPosition"
    assert m.parse_feature_location("AB000123.1:5..10", **QUIET)["ref"] == "AB000123.1"
    assert m.parse_feature_location("order(1..3,6..9)", **QUIET)["operator"] == "order"
    with pytest.raises(ConfigurationError):
        m.parse_feature_location("not a location", **QUIET)


def test_format_primersearch_input_text_and_file(tmp_path):
    out = tmp_path / "primers.txt"
    result = m.format_primersearch_input(
        {"p1": ["ACGTAC", "GGTACC"], "p2": ["AAAA", "TTTT"]}, output_path=out, **QUIET
    )
    assert result["text"].splitlines()[:2] == ["p1 ACGTAC GGTACC", "p2 AAAA TTTT"]
    assert out.read_text(encoding="utf-8") == result["text"]
    with pytest.raises(ConfigurationError):
        m.format_primersearch_input({"p1": ["ACGU", "GGTACC"]}, **QUIET)
    with pytest.raises(ConfigurationError):
        m.format_primersearch_input({"p1": ["ACGT"]}, **QUIET)


GAF = (
    "!gaf-version: 2.1\n"
    + "\t".join(
        [
            "UniProtKB",
            "P12345",
            "GENE1",
            "",
            "GO:0008150",
            "PMID:1",
            "IDA",
            "",
            "P",
            "gene one",
            "G1|G1b",
            "protein",
            "taxon:9606",
            "20200101",
            "UniProt",
            "",
            "",
        ]
    )
    + "\n"
    + "\t".join(
        [
            "UniProtKB",
            "Q99999",
            "GENE2",
            "NOT",
            "GO:0003674",
            "PMID:2|PMID:3",
            "IEA",
            "",
            "F",
            "gene two",
            "",
            "protein",
            "taxon:10090",
            "20200102",
            "UniProt",
            "",
            "",
        ]
    )
    + "\n"
)


def test_filter_go_annotations(tmp_path):
    path = tmp_path / "x.gaf"
    path.write_text(GAF, encoding="utf-8")
    by_symbol = m.filter_go_annotations(path, {"DB_Object_Symbol": ["GENE2"]}, **QUIET)
    assert by_symbol["records_read"] == 2 and by_symbol["kept_count"] == 1
    assert by_symbol["records"][0]["DB_Object_ID"] == "Q99999"
    by_reference = m.filter_go_annotations(path, {"DB:Reference": ["PMID:3"]}, **QUIET)
    assert by_reference["kept_count"] == 1
    either = m.filter_go_annotations(path, {"Aspect": ["P"], "Evidence": ["IEA"]}, **QUIET)
    assert either["kept_count"] == 2
    assert m.filter_go_annotations(path, {"Evidence": ["ISS"]}, **QUIET)["kept_count"] == 0
    with pytest.raises(ConfigurationError, match="not a field"):
        m.filter_go_annotations(path, {"Nope": ["x"]}, **QUIET)


def test_format_substitution_matrix_matches_biopython():
    from Bio.Align import substitution_matrices

    result = m.format_substitution_matrix("NUC.4.4", **QUIET)
    assert result["text"] == substitution_matrices.load("NUC.4.4").format("%i")
    assert result["alphabet"] == "ATGCSWRYKMBVHDN"
    assert "." in m.format_substitution_matrix("NUC.4.4", fmt="%.1f", **QUIET)["text"]
    with pytest.raises(ConfigurationError):
        m.format_substitution_matrix("NOT.A.MATRIX", **QUIET)
    with pytest.raises(ConfigurationError):
        m.format_substitution_matrix("NUC.4.4", fmt="%s%s", **QUIET)


def test_parse_entrez_xml_read_and_parse(tmp_path):
    info = tmp_path / "einfo.xml"
    info.write_text(data.ENTREZ_EINFO, encoding="utf-8")
    result = m.parse_entrez_xml(info, **QUIET)
    assert result["dtd"] == "einfo.dtd"
    assert "pubmed" in result["data"]["DbList"]
    search = tmp_path / "esearch.xml"
    search.write_text(data.ENTREZ_ESEARCH, encoding="utf-8")
    parsed = m.parse_entrez_xml(search, **QUIET)
    assert int(parsed["data"]["Count"]) == len(parsed["data"]["IdList"]) or parsed["data"]["IdList"]
    with pytest.raises(ConfigurationError, match="not valid Entrez XML"):
        m.parse_entrez_xml(search, multiple=True, **QUIET)
    summary = tmp_path / "esummary.xml"
    summary.write_text(data.ENTREZ_ESUMMARY, encoding="utf-8")
    many = m.parse_entrez_xml(summary, multiple=True, **QUIET)
    assert many["multiple"] is True and many["record_count"] >= 1


def test_parse_entrez_xml_does_not_download_missing_dtd(tmp_path):
    path = tmp_path / "unknown.xml"
    path.write_text(
        '<?xml version="1.0"?>\n<!DOCTYPE Nothing PUBLIC "-//X//DTD x//EN" '
        '"https://example.invalid/dtd/definitely_not_bundled.dtd">\n<Nothing/>\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="not available locally"):
        m.parse_entrez_xml(path, **QUIET)


def test_bgzf_virtual_offset_round_trip():
    made = m.bgzf_virtual_offset(100, 50, **QUIET)
    assert made["virtual_offset"] == (100 << 16) | 50
    split = m.bgzf_virtual_offset(virtual_offset=made["virtual_offset"], **QUIET)
    assert (split["block_start"], split["within_block"]) == (100, 50)
    with pytest.raises(ConfigurationError):
        m.bgzf_virtual_offset(**QUIET)
    with pytest.raises(ConfigurationError):
        m.bgzf_virtual_offset(1, 2, virtual_offset=3, **QUIET)
    with pytest.raises(ConfigurationError):
        m.bgzf_virtual_offset(2**48, 0, **QUIET)
    with pytest.raises(ConfigurationError):
        m.bgzf_virtual_offset(0, 2**16, **QUIET)


def test_public_functions_have_verbose_keyword_and_docstrings():
    import inspect

    for name in m.__all__:
        function = getattr(m, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
