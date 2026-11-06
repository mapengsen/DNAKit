"""Tests for dnakit.variants.hgvs_p35_functions (parsing, validation, mapping, projection, data queries).

Offline: the mapping tests run hgvs' own mappers on a synthetic local data snapshot (two toy transcripts on a
toy chromosome named like GRCh38 chromosome 1, so the AssemblyMapper accepts it). The UTA connection is only
checked for its guard because the postgres driver and the database are not available here.
"""

from __future__ import annotations

import copy
import inspect
import random
import socket

import pytest

pytest.importorskip("hgvs")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import hgvs_p35_functions as h

CHROM = "NC_000001.11"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class _NoConnectSocket(socket.socket):
    def connect(self, *args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    connect_ex = connect  # type: ignore[assignment]


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", _NoConnectSocket)


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT", "TGCA"))


def build_data() -> dict:
    rng = random.Random(5)
    stops = {"TAA", "TAG", "TGA", "ATG"}
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in stops]
    genome = [rng.choice("ACGT") for _ in range(1500)]
    plus_exons = [(100, 200), (300, 400), (500, 560)]
    cds = "ATG" + "".join(rng.choice(codons) for _ in range(58)) + "TAA"
    utr5 = "".join(rng.choice("ACGT") for _ in range(20))
    utr3 = "".join(rng.choice("ACGT") for _ in range(60))
    tx = utr5 + cds + utr3
    position = 0
    for start, end in plus_exons:
        genome[start:end] = list(tx[position : position + end - start])
        position += end - start
    minus_tx = (
        utr5[:10] + "ATG" + "".join(rng.choice(codons) for _ in range(40)) + "TGA" + utr3[:20]
    )
    first = 80
    e1 = (1000 - first, 1000)
    e2 = (e1[0] - 30 - (len(minus_tx) - first), e1[0] - 30)
    genome[e1[0] : e1[1]] = list(revcomp(minus_tx[:first]))
    genome[e2[0] : e2[1]] = list(revcomp(minus_tx[first:]))
    # transcript version 2: the same transcript with ten more bases at the 5' end (exon 1 starts at 90)
    genome[90:100] = list("ACGTACGTAC")
    tx2 = "ACGTACGTAC" + tx
    genome_text = "".join(genome)

    def exons(spans: list[tuple[int, int]]) -> list[dict]:
        rows, used = [], 0
        for index, (start, end) in enumerate(spans):
            rows.append(
                {
                    "ord": index,
                    "tx_start_i": used,
                    "tx_end_i": used + end - start,
                    "alt_start_i": start,
                    "alt_end_i": end,
                    "cigar": f"{end - start}=",
                }
            )
            used += end - start
        return rows

    return {
        "sequences": {
            CHROM: genome_text,
            "NM_999999.1": tx,
            "NM_999999.2": tx2,
            "NM_999998.1": minus_tx,
        },
        "transcripts": {
            "NM_999999.1": {
                "hgnc": "TOYA",
                "cds_start_i": 20,
                "cds_end_i": 200,
                "pro_ac": "NP_999999.1",
            },
            "NM_999999.2": {
                "hgnc": "TOYA",
                "cds_start_i": 30,
                "cds_end_i": 210,
                "pro_ac": "NP_999999.2",
            },
            "NM_999998.1": {"hgnc": "TOYB", "cds_start_i": 10, "cds_end_i": 136},
        },
        "alignments": [
            {
                "tx_ac": "NM_999999.1",
                "alt_ac": CHROM,
                "alt_aln_method": "splign",
                "alt_strand": 1,
                "exons": exons(plus_exons),
            },
            {
                "tx_ac": "NM_999999.2",
                "alt_ac": CHROM,
                "alt_aln_method": "splign",
                "alt_strand": 1,
                "exons": exons([(90, 200), (300, 400), (500, 560)]),
            },
            {
                "tx_ac": "NM_999998.1",
                "alt_ac": CHROM,
                "alt_aln_method": "splign",
                "alt_strand": -1,
                "exons": [
                    {
                        "ord": 0,
                        "tx_start_i": 0,
                        "tx_end_i": first,
                        "alt_start_i": e1[0],
                        "alt_end_i": e1[1],
                        "cigar": f"{first}=",
                    },
                    {
                        "ord": 1,
                        "tx_start_i": first,
                        "tx_end_i": len(minus_tx),
                        "alt_start_i": e2[0],
                        "alt_end_i": e2[1],
                        "cigar": f"{len(minus_tx) - first}=",
                    },
                ],
            },
        ],
        "genes": {"TOYA": {"maploc": "1p1", "descr": "toy gene A", "summary": "", "aliases": "A1"}},
    }


DATA = build_data()
GENOME = DATA["sequences"][CHROM]
KW = {"data_source": "local", "local_data": DATA, "verbose": False}


def snv(position: int) -> str:
    ref = GENOME[position - 1]
    return f"{CHROM}:g.{position}{ref}>{COMPLEMENT[ref]}"


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(h.__all__) == 14
    for name in h.__all__:
        parameter = inspect.signature(getattr(h, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


# Parsing.


def test_parse_posedit_and_position() -> None:
    coding = h.parse_hgvs_posedit("22+1A>T", verbose=False)
    assert coding["posedit"] == "22+1A>T" and coding["position"] == {"start": "22+1", "end": "22+1"}
    assert coding["edit"] == {"type": "sub", "ref": "A", "alt": "T"}
    protein = h.parse_hgvs_posedit("Ala22Trp", variant_type="p", verbose=False)
    assert protein["posedit"] == "Ala22Trp" and protein["edit"]["type"] == "sub"
    genomic = h.parse_hgvs_posedit("12_14del", variant_type="g", verbose=False)
    assert genomic["edit"]["type"] == "del" and genomic["position"]["end"] == "14"
    assert coding["uncertain"] is False
    with pytest.raises(ConfigurationError, match="posedit"):
        h.parse_hgvs_posedit("not a posedit", verbose=False)
    with pytest.raises(ConfigurationError, match="variant_type"):
        h.parse_hgvs_posedit("1A>T", variant_type="x", verbose=False)

    bare = h.parse_hgvs_position("26+2_57-3", variant_type="c", verbose=False)
    # Coding coordinates are measured from CDS_START, as the upstream parser reports.
    from hgvs.parser import Parser

    direct = Parser().parse_c_interval("26+2_57-3")
    assert bare["start"] == {
        "text": "26+2",
        "base": 26,
        "offset": 2,
        "datum": direct.start.datum.name,
    }
    assert bare["end"]["offset"] == -3
    prefixed = h.parse_hgvs_position("NM_000001.1:c.26+2_57-3", verbose=False)
    assert prefixed["accession"] == "NM_000001.1" and prefixed["type"] == "c"
    typed = h.parse_hgvs_position("NC_000001.11:g.10_20", variant_type="g", verbose=False)
    assert typed["start"]["base"] == 10 and typed["end"]["base"] == 20
    with pytest.raises(ConfigurationError, match="valid HGVS position"):
        h.parse_hgvs_position("26+2_57-3", verbose=False)


# Validation and small variant operations.


def test_validate_position_and_variant() -> None:
    assert h.validate_hgvs_position("26_57", variant_type="c", verbose=False)["valid"] is True
    reversed_interval = h.validate_hgvs_position("57_26", variant_type="c", verbose=False)
    assert reversed_interval["valid"] is False and reversed_interval["level"] == "ERROR"
    assert h.validate_hgvs_position("0_5", variant_type="g", verbose=False)["valid"] is False

    assert h.validate_hgvs_variant("NM_999999.1:c.5A>T", verbose=False)["valid"] is True
    bad_type = h.validate_hgvs_variant("NM_999999.1:g.5A>T", verbose=False)
    assert bad_type["valid"] is False and "NM_" in bad_type["message"]
    bad_length = h.validate_hgvs_variant("NM_999999.1:c.5_7del2", verbose=False)
    assert bad_length["valid"] is False

    start = DATA["sequences"]["NM_999999.1"][20:23]
    assert start == "ATG"
    good = h.validate_hgvs_variant("NM_999999.1:c.1A>T", level="full", **KW)
    assert good["valid"] is True and good["level"] == "full"
    wrong_ref = h.validate_hgvs_variant("NM_999999.1:c.1C>T", level="extrinsic", **KW)
    assert wrong_ref["valid"] is False and "reference" in wrong_ref["message"]
    outside = h.validate_hgvs_variant("NM_999999.1:c.500A>T", level="extrinsic", **KW)
    assert outside["valid"] is False
    with pytest.raises(ConfigurationError, match="local_data"):
        h.validate_hgvs_variant("NM_999999.1:c.1A>T", level="extrinsic", verbose=False)


def test_fill_reference_and_length_change() -> None:
    tx = DATA["sequences"]["NM_999999.1"]
    filled = h.fill_hgvs_reference("NM_999999.1:c.5_7del", **KW)
    assert (
        filled["output"]["hgvs"] == f"NM_999999.1:c.5_7del{tx[24:27]}"
        and filled["reference"] == tx[24:27]
    )
    unchanged = h.fill_hgvs_reference("NM_999999.1:c.5A>T", **KW)
    assert unchanged["output"]["hgvs"] == "NM_999999.1:c.5A>T"
    changes = {
        "NM_1.1:c.5_7delACG": -3,
        "NM_1.1:c.5A>T": 0,
        "NM_1.1:c.5_6insTT": 2,
        "NM_1.1:c.5_7dup": 3,
        "NM_1.1:c.5_7delinsAA": -1,
    }
    for text, expected in changes.items():
        assert h.hgvs_length_change(text, verbose=False)["length_change"] == expected, text
    # Explicit deleted bases are optional: HGVS can infer the length from the interval.
    from hgvs.parser import Parser

    direct = Parser().parse_hgvs_variant("NM_1.1:c.5_7del")
    assert h.hgvs_length_change("NM_1.1:c.5_7del", verbose=False)[
        "length_change"
    ] == direct.posedit.length_change(on_error_raise=True)
    with pytest.raises(ConfigurationError, match="not a valid HGVS variant"):
        h.hgvs_length_change("garbage", verbose=False)


# Mapping.


def test_map_variant_plus_strand_round_trip() -> None:
    genomic = snv(120)  # first base before the start codon
    to_c = h.map_hgvs_variant(genomic, to_type="c", transcript_accession="NM_999999.1", **KW)
    assert to_c["output"]["hgvs"] == "NM_999999.1:c.-1C>G" or to_c["output"]["type"] == "c"
    coding = to_c["output"]["hgvs"]
    assert to_c["mapper_method"] == "g_to_c" and to_c["mapper"] == "VariantMapper"
    assert (
        h.map_hgvs_variant(coding, to_type="g", genomic_accession=CHROM, **KW)["output"]["hgvs"]
        == genomic
    )
    non_coding = h.map_hgvs_variant(genomic, to_type="n", transcript_accession="NM_999999.1", **KW)[
        "output"
    ]["hgvs"]
    assert non_coding == "NM_999999.1:n.20C>G"
    assert (
        h.map_hgvs_variant(non_coding, to_type="g", genomic_accession=CHROM, **KW)["output"]["hgvs"]
        == genomic
    )
    assert h.map_hgvs_variant(non_coding, to_type="c", **KW)["output"]["hgvs"] == coding
    assert h.map_hgvs_variant(coding, to_type="n", **KW)["output"]["hgvs"] == non_coding
    automatic = h.map_hgvs_variant(genomic, to_type="t", transcript_accession="NM_999999.1", **KW)
    assert automatic["output"]["hgvs"] == coding and automatic["mapper_method"] == "g_to_t"
    back = h.map_hgvs_variant(coding, to_type="t", genomic_accession=CHROM, **KW)
    assert back["output"]["hgvs"] == genomic and back["mapper_method"] == "t_to_g"


def test_map_variant_with_assembly_mapper_and_minus_strand() -> None:
    genomic = snv(950)
    via_assembly = h.map_hgvs_variant(
        genomic, to_type="c", transcript_accession="NM_999998.1", assembly_name="GRCh38", **KW
    )
    assert via_assembly["mapper"] == "AssemblyMapper"
    explicit = h.map_hgvs_variant(genomic, to_type="c", transcript_accession="NM_999998.1", **KW)
    assert via_assembly["output"]["hgvs"] == explicit["output"]["hgvs"]
    coding = explicit["output"]["hgvs"]
    ref = GENOME[949]
    assert coding.endswith(
        f"{COMPLEMENT[ref]}>{COMPLEMENT[COMPLEMENT[ref]]}"
    )  # the minus strand flips the bases
    assert (
        h.map_hgvs_variant(coding, to_type="g", assembly_name="GRCh38", **KW)["output"]["hgvs"]
        == genomic
    )


def test_map_variant_argument_errors() -> None:
    with pytest.raises(ConfigurationError, match="transcript_accession"):
        h.map_hgvs_variant(snv(120), to_type="c", **KW)
    with pytest.raises(ConfigurationError, match="genomic_accession or assembly_name"):
        h.map_hgvs_variant("NM_999999.1:c.5A>T", to_type="g", **KW)
    with pytest.raises(ConfigurationError, match="cannot map"):
        h.map_hgvs_variant(snv(120), to_type="g", transcript_accession="NM_999999.1", **KW)
    with pytest.raises(ConfigurationError, match="reference"):
        h.map_hgvs_variant(
            f"{CHROM}:g.120A>G", to_type="c", transcript_accession="NM_999999.1", **KW
        )
    unchecked = h.map_hgvs_variant(
        f"{CHROM}:g.120A>G",
        to_type="c",
        transcript_accession="NM_999999.1",
        prevalidation_level="NONE",
        replace_reference=False,
        **KW,
    )
    assert unchecked["output"]["posedit"].startswith("-1")
    with pytest.raises(ConfigurationError, match=r"No transcript definition|NM_000000"):
        h.map_hgvs_variant(snv(120), to_type="c", transcript_accession="NM_000000.1", **KW)


def test_map_interval() -> None:
    to_genome = h.map_hgvs_interval(
        "1_10",
        from_type="c",
        to_type="g",
        transcript_accession="NM_999999.1",
        genomic_accession=CHROM,
        **KW,
    )
    assert to_genome["output"] == "121_130"
    back = h.map_hgvs_interval(
        "121_130",
        from_type="g",
        to_type="c",
        transcript_accession="NM_999999.1",
        genomic_accession=CHROM,
        **KW,
    )
    assert back["output"] == "1_10"
    non_coding = h.map_hgvs_interval(
        "1_10", from_type="c", to_type="n", transcript_accession="NM_999999.1", **KW
    )
    assert non_coding["output"] == "21_30" and non_coding["start"] == "21"
    across_exons = h.map_hgvs_interval(
        "70_130",
        from_type="c",
        to_type="g",
        transcript_accession="NM_999999.1",
        genomic_accession=CHROM,
        **KW,
    )
    assert (
        across_exons["output"] == "190_350"
    )  # c.70 -> transcript index 89 in exon 1, c.130 -> index 149 in exon 2
    with pytest.raises(ConfigurationError, match="differ"):
        h.map_hgvs_interval(
            "1_10", from_type="c", to_type="c", transcript_accession="NM_999999.1", **KW
        )
    with pytest.raises(ConfigurationError, match="genomic_accession"):
        h.map_hgvs_interval(
            "1_10", from_type="c", to_type="g", transcript_accession="NM_999999.1", **KW
        )
    with pytest.raises(ConfigurationError, match=r"out of bounds|bounds"):
        h.map_hgvs_interval(
            "1_1000", from_type="c", to_type="n", transcript_accession="NM_999999.1", **KW
        )


def test_predict_protein_consequence() -> None:
    tx = DATA["sequences"]["NM_999999.1"]
    codon = tx[23:26]
    alt = COMPLEMENT[codon[0]]
    result = h.predict_hgvs_protein_consequence(f"NM_999999.1:c.4{codon[0]}>{alt}", **KW)
    assert result["output"]["accession"] == "NP_999999.1" and result["output"]["type"] == "p"
    assert (
        result["output"]["posedit"].startswith("(") and result["output"]["posedit"].count("2") == 1
    )
    assembly = h.predict_hgvs_protein_consequence(
        f"NM_999999.1:c.4{codon[0]}>{alt}", assembly_name="GRCh38", **KW
    )
    assert assembly["output"]["hgvs"] == result["output"]["hgvs"]
    start_loss = h.predict_hgvs_protein_consequence("NM_999999.1:c.1A>G", **KW)
    assert start_loss["output"]["posedit"] == "Met1?"
    assert h.predict_hgvs_protein_consequence("NM_999999.1:n.25A>G", **KW)["output"] == "non-coding"
    with pytest.raises(ConfigurationError, match="coding"):
        h.predict_hgvs_protein_consequence(snv(120), **KW)


def test_project_variant_between_transcripts() -> None:
    forward = h.project_hgvs_variant(
        "NM_999999.1:c.5A>T".replace(
            "A>T",
            f"{DATA['sequences']['NM_999999.1'][24]}>{COMPLEMENT[DATA['sequences']['NM_999999.1'][24]]}",
        ),
        source_accession="NM_999999.1",
        destination_accession="NM_999999.2",
        genomic_accession=CHROM,
        **KW,
    )
    assert forward["kind"] == "variant" and forward["output"]["accession"] == "NM_999999.2"
    assert forward["output"]["posedit"].startswith(
        "5"
    )  # same CDS-relative position on the longer transcript
    interval = h.project_hgvs_variant(
        "-5_10",
        source_accession="NM_999999.1",
        destination_accession="NM_999999.2",
        genomic_accession=CHROM,
        **KW,
    )
    assert interval["kind"] == "interval" and interval["output"] == "-5_10"
    backward = h.project_hgvs_variant(
        forward["output"]["hgvs"],
        source_accession="NM_999999.1",
        destination_accession="NM_999999.2",
        genomic_accession=CHROM,
        direction="backward",
        **KW,
    )
    assert backward["output"]["accession"] == "NM_999999.1"
    with pytest.raises(ConfigurationError, match="direction"):
        h.project_hgvs_variant(
            "-5_10",
            source_accession="NM_999999.1",
            destination_accession="NM_999999.2",
            genomic_accession=CHROM,
            direction="up",
            **KW,
        )


def test_relevant_transcripts_and_vcf_allele() -> None:
    relevant = h.list_hgvs_relevant_transcripts(snv(120), **KW)
    assert relevant["transcripts"] == ["NM_999999.1", "NM_999999.2"] or sorted(
        relevant["transcripts"]
    ) == [
        "NM_999999.1",
        "NM_999999.2",
    ]
    assert h.list_hgvs_relevant_transcripts(snv(1200), **KW)["transcript_count"] == 0
    with pytest.raises(ConfigurationError, match="genomic"):
        h.list_hgvs_relevant_transcripts("NM_999999.1:c.5A>T", **KW)

    allele = h.hgvs_to_vcf_allele(snv(120), **KW)
    assert (allele["chrom"], allele["pos"], allele["ref"], allele["edit_type"]) == (
        "1",
        120,
        GENOME[119],
        "sub",
    )
    deletion = h.hgvs_to_vcf_allele(f"{CHROM}:g.120_122del", **KW)
    assert (deletion["edit_type"] == "del" and deletion["ref"] == GENOME[117:122]) or deletion[
        "ref"
    ].startswith(GENOME[117])
    assert deletion["alt"] == deletion["ref"][0] and len(deletion["ref"]) == 4
    duplication = h.hgvs_to_vcf_allele(f"{CHROM}:g.120_122dup", **KW)
    assert (
        duplication["edit_type"] == "dup" and len(duplication["alt"]) == len(duplication["ref"]) + 3
    )
    with pytest.raises(ConfigurationError, match="g\\. variant"):
        h.hgvs_to_vcf_allele("NM_999999.1:c.5A>T", **KW)


# Data-provider queries.


def test_fetch_reference_sequence() -> None:
    result = h.fetch_hgvs_reference_sequence(CHROM, start_i=100, end_i=110, **KW)
    assert result["sequence"] == GENOME[100:110] and result["length"] == 10
    assert h.fetch_hgvs_reference_sequence("NM_999999.1", **KW)["length"] == 260
    with pytest.raises(ConfigurationError, match="No sequence"):
        h.fetch_hgvs_reference_sequence("NC_404", **KW)


def test_transcript_database_queries() -> None:
    query = h.query_hgvs_transcript_database
    gene = query("gene_info", gene="TOYA", **KW)["rows"]
    assert gene[0]["hgnc"] == "TOYA" and gene[0]["maploc"] == "1p1"
    assert query("gene_info", gene="NOPE", **KW)["rows"] == []
    for_gene = query("tx_for_gene", gene="TOYA", **KW)["rows"]
    assert [row["tx_ac"] for row in for_gene] == ["NM_999999.1", "NM_999999.2"]
    exons = query(
        "tx_exons",
        transcript_accession="NM_999998.1",
        genomic_accession=CHROM,
        alt_aln_method="splign",
        **KW,
    )
    assert [row["alt_start_i"] for row in exons["rows"]] == sorted(
        row["alt_start_i"] for row in exons["rows"]
    )
    assert exons["rows"][0]["alt_strand"] == -1 and exons["row_count"] == 2
    identity = query("tx_identity_info", transcript_accession="NM_999999.1", **KW)["rows"][0]
    assert identity["lengths"] == [100, 100, 60] and (
        identity["cds_start_i"],
        identity["cds_end_i"],
    ) == (20, 200)
    info = query(
        "tx_info",
        transcript_accession="NM_999999.1",
        genomic_accession=CHROM,
        alt_aln_method="splign",
        **KW,
    )["rows"][0]
    assert info["hgnc"] == "TOYA" and info["alt_ac"] == CHROM
    options = query("tx_mapping_options", transcript_accession="NM_999999.1", **KW)["rows"]
    assert options == [{"tx_ac": "NM_999999.1", "alt_ac": CHROM, "alt_aln_method": "splign"}]
    similar = query("similar_transcripts", transcript_accession="NM_999999.1", **KW)["rows"]
    assert (
        similar
        == [
            {
                "tx_ac1": "NM_999999.1",
                "tx_ac2": "NM_999999.2",
                "cds_eq": True,
                "es_fp_eq": False,
                "cds_es_fp_eq": False,
            }
        ]
        or similar[0]["cds_eq"] is True
    )
    assert similar[0]["es_fp_eq"] is False  # exon 1 starts at a different genomic position
    assert query("protein_accession", transcript_accession="NM_999999.1", **KW)["rows"] == [
        {"tx_ac": "NM_999999.1", "pro_ac": "NP_999999.1"}
    ]
    assert query("protein_accession", transcript_accession="NM_999998.1", **KW)["rows"] == []
    assembly = query("assembly_map", assembly_name="GRCh38", **KW)
    assert {"accession": CHROM, "name": "1"} in assembly["rows"]
    with pytest.raises(ConfigurationError, match="gene is required"):
        query("gene_info", **KW)
    with pytest.raises(ConfigurationError, match="query"):
        query("everything", **KW)
    with pytest.raises(ConfigurationError, match="No tx_exons"):
        query(
            "tx_exons",
            transcript_accession="NM_999998.1",
            genomic_accession="NC_X",
            alt_aln_method="splign",
            **KW,
        )


def test_local_data_is_validated() -> None:
    bad = copy.deepcopy(DATA)
    bad["alignments"][0]["alt_strand"] = 0
    with pytest.raises(ConfigurationError, match="alt_strand"):
        h.fetch_hgvs_reference_sequence(CHROM, data_source="local", local_data=bad, verbose=False)
    with pytest.raises(ConfigurationError, match="unknown keys"):
        h.fetch_hgvs_reference_sequence(
            CHROM, data_source="local", local_data={"x": 1}, verbose=False
        )
    with pytest.raises(ConfigurationError, match="needs local_data"):
        h.fetch_hgvs_reference_sequence(CHROM, verbose=False)
    with pytest.raises(ConfigurationError, match="data_source"):
        h.fetch_hgvs_reference_sequence(CHROM, data_source="nowhere", verbose=False)


def test_uta_source_needs_the_postgres_driver() -> None:
    import importlib.util

    if importlib.util.find_spec("psycopg2") is not None:
        pytest.skip("psycopg2 is installed; the live UTA path is not exercised offline")
    with pytest.raises(BackendUnavailableError, match="psycopg2"):
        h.query_hgvs_transcript_database("gene_info", gene="TP53", data_source="uta", verbose=False)
