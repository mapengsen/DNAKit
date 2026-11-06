from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.variants import pysam_p36_variant_functions as m

SPEC = dict(
    contigs=[{"id": "chr1", "length": 1000}],
    info_fields=[{"id": "DP", "number": 1, "type": "Integer", "description": "Depth"}],
    format_fields=[{"id": "GT", "number": 1, "type": "String", "description": "Genotype"}],
    filters=[{"id": "q10", "description": "Low quality"}],
    meta=[{"key": "fileDate", "value": "20260101"}],
    samples=["S1", "S2"],
)
RECORDS = [
    {
        "contig": "chr1",
        "pos": 100,
        "ref": "A",
        "alts": ["T"],
        "id": "r1",
        "qual": 30.5,
        "filter": ["PASS"],
        "info": {"DP": 5},
        "samples": {"S1": {"GT": [0, 1], "phased": True}, "S2": {"GT": [1, 1]}},
    },
    {
        "contig": "chr1",
        "pos": 200,
        "ref": "G",
        "alts": ["C", "T"],
        "qual": 1.0,
        "filter": ["q10"],
        "samples": {"S1": {"GT": [0, 2]}},
    },
]


def test_write_vcf_records_reads_back_with_pysam(tmp_path):
    out = tmp_path / "w.vcf"
    result = m.write_vcf_records(out, RECORDS, **SPEC, verbose=False)
    assert result["record_count"] == 2 and result["format"] == "VCF" and result["sample_count"] == 2
    with pysam.VariantFile(str(out)) as vcf:
        assert list(vcf.header.samples) == ["S1", "S2"]
        assert vcf.header.contigs["chr1"].length == 1000
        assert "q10" in vcf.header.filters and "DP" in vcf.header.info
        first, second = list(vcf)
    assert (first.pos, first.ref, first.alts, first.id, first.qual) == (
        100,
        "A",
        ("T",),
        "r1",
        30.5,
    )
    assert first.info["DP"] == 5 and list(first.filter) == ["PASS"]
    assert first.samples["S1"]["GT"] == (0, 1) and first.samples["S1"].phased
    assert first.samples["S2"]["GT"] == (1, 1) and not first.samples["S2"].phased
    assert (
        second.alts == ("C", "T")
        and list(second.filter) == ["q10"]
        and second.samples["S1"]["GT"] == (0, 2)
    )
    assert "##fileDate=20260101" in out.read_text()


def test_write_vcf_formats_by_suffix(tmp_path):
    gz = m.write_vcf_records(tmp_path / "w.vcf.gz", RECORDS, **SPEC, verbose=False)
    bcf = m.write_vcf_records(tmp_path / "w.bcf", RECORDS, **SPEC, verbose=False)
    assert gz["format"] == "VCF (bgzip)" and bcf["format"] == "BCF"
    for name in ("w.vcf.gz", "w.bcf"):
        with pysam.VariantFile(str(tmp_path / name)) as vcf:
            assert len(list(vcf)) == 2
    assert (tmp_path / "w.vcf.gz").read_bytes()[:2] == b"\x1f\x8b"


def test_write_vcf_header_lines_and_validation(tmp_path):
    out = tmp_path / "h.vcf"
    m.write_vcf_records(
        out, [], contigs=[{"id": "c1"}], header_lines=['##source="unit test"'], verbose=False
    )
    text = out.read_text()
    assert '##source="unit test"' in text and "##contig=<ID=c1>" in text
    with pytest.raises(ConfigurationError):
        m.write_vcf_records(out, [], verbose=False)  # exists
    with pytest.raises(ConfigurationError):
        m.write_vcf_records(tmp_path / "x.txt", [], verbose=False)
    with pytest.raises(BackendExecutionError):
        m.write_vcf_records(
            tmp_path / "bad.vcf",
            [{"contig": "nope", "pos": 1, "ref": "A", "alts": ["T"]}],
            **SPEC,
            verbose=False,
        )
    assert not (tmp_path / "bad.vcf").exists()


def test_merge_headers_adds_missing_lines_only(tmp_path):
    first, second = tmp_path / "a.vcf", tmp_path / "b.vcf"
    m.write_vcf_records(first, RECORDS, **SPEC, verbose=False)
    m.write_vcf_records(
        second,
        [],
        contigs=[{"id": "chr1", "length": 1000}],
        info_fields=[
            {"id": "AF", "number": "A", "type": "Float", "description": "Allele frequency"}
        ],
        samples=["S9"],
        verbose=False,
    )
    result = m.merge_vcf_headers(first, second, tmp_path / "merged.vcf", verbose=False)
    assert result["added_header_lines"] == [
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    ]
    assert result["samples"] == ["S1", "S2"]  # pysam does not merge sample columns
    with pysam.VariantFile(str(tmp_path / "merged.vcf")) as vcf:
        assert {"DP", "AF"} <= set(vcf.header.info)
        assert list(vcf) == []


def test_translate_records_to_merged_header(tmp_path):
    first, second = tmp_path / "a.vcf", tmp_path / "b.vcf"
    m.write_vcf_records(first, RECORDS, **SPEC, verbose=False)
    m.write_vcf_records(
        second,
        [],
        **{
            **SPEC,
            "info_fields": SPEC["info_fields"]
            + [{"id": "AF", "number": "A", "type": "Float", "description": "AF"}],
        },
        verbose=False,
    )
    result = m.translate_vcf_records(first, second, tmp_path / "t.vcf", verbose=False)
    assert result["record_count"] == 2
    with pysam.VariantFile(str(tmp_path / "t.vcf")) as vcf, pysam.VariantFile(str(first)) as ref:
        assert "AF" in vcf.header.info
        assert [(r.pos, r.alleles, r.qual) for r in vcf] == [
            (r.pos, r.alleles, r.qual) for r in ref
        ]
    # a header that lacks the contig of the source cannot take the records
    third = tmp_path / "c.vcf"
    m.write_vcf_records(
        third, [], contigs=[{"id": "chr9", "length": 5}], samples=["S1", "S2"], verbose=False
    )
    with pytest.raises(BackendExecutionError):
        m.translate_vcf_records(first, third, tmp_path / "bad.vcf", verbose=False)
    assert not (tmp_path / "bad.vcf").exists()
    with pytest.raises(ConfigurationError):
        m.translate_vcf_records(first, second, tmp_path / "t2.vcf", max_records=0, verbose=False)


@pytest.mark.parametrize(
    ("call_a", "call_b", "expected"),
    [
        ((10, "ATT", "AT"), (10, "AT", "A"), True),
        ((10, "A", "C"), (10, "A", "G"), False),
        ((10, "A", "C"), (11, "A", "C"), False),
        ((5, "GA", "GAA"), (5, "G", "GA"), True),
    ],
)
def test_compare_variant_calls_equals_pysam(call_a, call_b, expected):
    result = m.compare_variant_calls(*call_a, *call_b, verbose=False)
    assert result["equal"] is expected
    assert result["equal"] == pysam.VCF().compare_calls(*call_a, *call_b)


def test_compare_variant_calls_validation():
    with pytest.raises(ConfigurationError):
        m.compare_variant_calls(1, "", "A", 1, "A", "A", verbose=False)


def test_iupac_genotype_codes_both_ways():
    encoded = m.convert_iupac_genotype_codes(["AG", "GA", "CT", "GG", "ZZ"], verbose=False)
    assert encoded["results"] == {"AG": "r", "GA": "R", "CT": "y", "GG": "G", "ZZ": None}
    assert encoded["errors"] == ["ZZ"]
    decoded = m.convert_iupac_genotype_codes(
        ["r", "R", "m", "G"], direction="decode", verbose=False
    )
    assert decoded["results"] == {"r": "AG", "R": "AG", "m": "AC", "G": "GG"}
    with pytest.raises(ConfigurationError):
        m.convert_iupac_genotype_codes([], verbose=False)
    with pytest.raises(ConfigurationError):
        m.convert_iupac_genotype_codes(["AG"], direction="sideways", verbose=False)


def test_gtf_and_gff3_attribute_conversion():
    gtf = m.convert_gff_attributes(
        attribute_string='gene_id "A"; gene_name "B"; z "1,2";', verbose=False
    )
    assert gtf["attributes"] == {"gene_id": "A", "gene_name": "B", "z": "1,2"}
    back = m.convert_gff_attributes(attributes={"gene_id": "X", "n": 5, "s": "a b"}, verbose=False)
    assert back["attribute_string"] == 'gene_id "X"; n 5; s "a b";'
    gff = m.convert_gff_attributes(attribute_string="ID=a;Name=b%2Cc", flavor="gff3", verbose=False)
    assert gff["attributes"] == {"ID": "a", "Name": "b%2Cc"}
    text = m.convert_gff_attributes(
        attributes={"ID": "a", "Name": "b"}, flavor="gff3", verbose=False
    )
    assert text["attribute_string"] == "ID=a;Name=b"
    with pytest.raises(ConfigurationError):
        m.convert_gff_attributes(verbose=False)
    with pytest.raises(ConfigurationError):
        m.convert_gff_attributes(attribute_string="a=1", attributes={"a": 1}, verbose=False)
    with pytest.raises(ConfigurationError):
        m.convert_gff_attributes(attribute_string="a=1", flavor="bed", verbose=False)
    with pytest.raises(ConfigurationError):
        m.convert_gff_attributes(attributes={}, verbose=False)
