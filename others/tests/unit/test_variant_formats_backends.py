"""Contract tests for PLINK, BGEN and varcode adapters in dnakit.variants."""

from __future__ import annotations

import math
import struct
import zlib
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants.genotype_formats import (
    calculate_bgen_allele_frequency,
    read_bgen_probabilities,
    read_bgen_variants,
    read_plink_bed,
    read_plink_genotypes,
    read_plink_metadata,
    summarize_plink_bed,
    write_plink_bed,
    write_plink_genotypes,
)
from dnakit.variants.nomenclature import (
    check_hgvs_syntax,
    describe_variant,
    format_hgvs_variant,
    load_vcf_variants,
    parse_hgvs_variant,
)

SAMPLES = ["S1", "S2", "S3"]
VARIANT_IDS = ["rs1", "rs2", "rs3"]
GENOTYPES = [[0, 1, 2], [2, None, 1], [1, 1, 0]]
PLINK_KWARGS = {
    "chromosomes": ["1", "1", "2"],
    "positions": [100, 200, 300],
    "allele_1": ["A", "C", "G"],
    "allele_2": ["G", "T", "A"],
}


def _write_plink(tmp_path: Path, genotypes=GENOTYPES, name: str = "panel.bed") -> Path:
    path = tmp_path / name
    write_plink_bed(
        path,
        genotypes,
        sample_ids=SAMPLES,
        variant_ids=VARIANT_IDS,
        verbose=False,
        **PLINK_KWARGS,
    )
    return path


def _bgen_bytes(samples: list[str], variants: list[tuple], *, phased: bool = False) -> bytes:
    """BGEN 1.2 layout 2 with zlib compression, one tuple per variant.

    A variant tuple is (id, rsid, chromosome, position, alleles, probabilities, missing_samples).
    Unphased probabilities are integers in [0, 255] for the first two genotype combinations; the
    last combination is implied. Phased blocks store two values per sample for diploid biallelic
    data. For compressed blocks C counts D plus the zlib bytes.
    """

    n = len(samples)
    flags = 0x80000000 | (2 << 2) | 1  # sample identifiers, layout 2, zlib
    header = struct.pack("<III", 20, len(variants), n) + b"bgen" + struct.pack("<I", flags)
    names = b"".join(struct.pack("<H", len(s)) + s.encode() for s in samples)
    sample_block = struct.pack("<II", 8 + len(names), n) + names
    body = b""
    for ident, rsid, chrom, pos, alleles, probs, missing in variants:
        body += struct.pack("<H", len(ident)) + ident.encode()
        body += struct.pack("<H", len(rsid)) + rsid.encode()
        body += struct.pack("<H", len(chrom)) + chrom.encode()
        body += struct.pack("<IH", pos, len(alleles))
        for allele in alleles:
            body += struct.pack("<I", len(allele)) + allele.encode()
        ploidy = [2 | (0x80 if index in missing else 0) for index in range(n)]
        raw = struct.pack("<IHBB", n, len(alleles), 2, 2) + bytes(ploidy)
        raw += bytes([1 if phased else 0, 8])  # phasing flag, 8 bits per probability
        for row in probs:
            raw += bytes(row)
        packed = zlib.compress(raw)
        body += struct.pack("<II", 4 + len(packed), len(raw)) + packed
    offset = len(header) + len(sample_block)
    return struct.pack("<I", offset) + header + sample_block + body


def _write_bgen(
    tmp_path: Path,
    name: str,
    samples: list[str],
    variants: list[tuple],
    *,
    phased: bool = False,
) -> Path:
    path = tmp_path / name
    path.write_bytes(_bgen_bytes(samples, variants, phased=phased))
    return path


def _two_variant_bgen(tmp_path: Path) -> Path:
    # S2 is missing at rs1 only; probabilities are AA, AB (BB implied) in 1/255 steps.
    return _write_bgen(
        tmp_path,
        "two.bgen",
        ["S1", "S2", "S3"],
        [
            ("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 255], [128, 127]], ()),
            ("v2", "rs2", "1", 200, ["C", "T"], [[0, 0], [255, 0], [0, 0]], ()),
        ],
    )


def _write_vcf(tmp_path: Path) -> Path:
    text = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t101\trs1\tA\tG\t50\tPASS\tDP=5\n"
        "1\t200\t.\tACG\tA\t40\tPASS\tDP=6\n"
        "1\t300\t.\tA\tAGT\t.\tq10\tDP=3\n"
        "1\t400\t.\tAC\tTG\t.\t.\tDP=3\n"
        "1\t500\t.\tA\tG\t.\tPASS\tDP=3\n"
    )
    path = tmp_path / "records.vcf"
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------- bed-reader and PLINK


def test_plink_round_trip_keeps_dosages_missing_calls_and_ids(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    path = _write_plink(tmp_path)

    result = read_plink_bed(path, verbose=False)

    assert result["samples"] == SAMPLES
    assert [variant["id"] for variant in result["variants"]] == VARIANT_IDS
    assert result["variants"][1] == {
        "id": "rs2",
        "chromosome": "1",
        "position": 200,
        "allele_1": "C",
        "allele_2": "T",
    }
    assert result["genotypes"] == GENOTYPES


def test_read_plink_bed_matches_direct_bed_reader_call(tmp_path: Path) -> None:
    bed_reader = pytest.importorskip("bed_reader")
    path = _write_plink(tmp_path)

    direct = bed_reader.open_bed(str(path), count_A1=True).read()
    result = read_plink_bed(path, verbose=False)

    for row, values in zip(result["genotypes"], direct.tolist(), strict=True):
        expected = [None if math.isnan(value) else int(value) for value in values]
        assert row == expected


def test_read_plink_bed_selects_samples_variants_and_counted_allele(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    path = _write_plink(tmp_path)

    selected = read_plink_bed(
        path, samples=["S3", "S1"], variants=["rs3"], count_allele="allele_2", verbose=False
    )
    everything = read_plink_bed(path, count_allele="allele_2", verbose=False)

    # Counting allele_2 is the complement of allele_1 where both alleles are called.
    assert selected["samples"] == ["S3", "S1"]
    assert selected["genotypes"] == [[2], [0]]
    assert selected["variants"][0]["id"] == "rs3"
    assert everything["genotypes"] == [[2, 1, 0], [0, None, 1], [1, 1, 2]]
    assert selected["count_allele"] == "allele_2"


def test_summarize_plink_bed_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    path = _write_plink(
        tmp_path,
        genotypes=[[0, 2, 1], [0, 2, None], [0, None, 2]],
        name="summary.bed",
    )

    result = summarize_plink_bed(path, verbose=False)

    # rs1: all zero -> frequency 0. rs2: sum 4 over 2 calls -> 1.0. rs3: sum 3 over 2 calls -> 0.75.
    assert result["samples"] == 3
    assert [item["allele_frequency"] for item in result["variants"]] == [0.0, 1.0, 0.75]
    assert [item["missing"] for item in result["variants"]] == [0, 1, 1]
    assert result["missing_genotypes"] == 2


def test_write_plink_bed_refuses_bad_values_and_overwrites(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    path = _write_plink(tmp_path)

    with pytest.raises(ConfigurationError, match="Refusing to overwrite"):
        _write_plink(tmp_path)
    with pytest.raises(ConfigurationError, match="0, 1, 2 or None"):
        write_plink_bed(
            tmp_path / "bad.bed",
            [[0, 3, 1], [0, 1, 1], [0, 1, 1]],
            sample_ids=SAMPLES,
            variant_ids=VARIANT_IDS,
            verbose=False,
            **PLINK_KWARGS,
        )
    with pytest.raises(ConfigurationError, match="must have 3 entries"):
        write_plink_bed(
            tmp_path / "short.bed",
            GENOTYPES,
            sample_ids=SAMPLES,
            variant_ids=VARIANT_IDS,
            verbose=False,
            **{**PLINK_KWARGS, "chromosomes": ["1", "1"]},
        )
    with pytest.raises(ConfigurationError, match="must be unique"):
        write_plink_bed(
            tmp_path / "dup.bed",
            GENOTYPES,
            sample_ids=["S1", "S1", "S3"],
            variant_ids=VARIANT_IDS,
            verbose=False,
            **PLINK_KWARGS,
        )
    assert path.exists()


def test_read_plink_bed_enforces_record_and_cell_limits(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    path = _write_plink(tmp_path)

    with pytest.raises(ConfigurationError, match="more than 2 records"):
        read_plink_bed(path, max_records=2, verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells=8"):
        read_plink_bed(path, max_cells=8, verbose=False)
    with pytest.raises(ConfigurationError, match="not in the file"):
        read_plink_bed(path, samples=["S9"], verbose=False)


def test_read_plink_metadata_matches_pandas_plink_tables(tmp_path: Path) -> None:
    pandas_plink = pytest.importorskip("pandas_plink")
    path = _write_plink(tmp_path)

    result = read_plink_metadata(path, verbose=False)
    bim, fam, _bed = pandas_plink.read_plink(str(path.with_suffix("")), verbose=False)

    assert [row["sample_id"] for row in result["samples"]] == fam["iid"].tolist()
    assert [row["id"] for row in result["variants"]] == bim["snp"].tolist()
    assert [row["position"] for row in result["variants"]] == bim["pos"].tolist()
    assert result["variants"][2]["allele_1"] == "G"
    assert result["variants"][2]["allele_2"] == "A"


def test_read_plink_genotypes_agrees_with_bed_reader_for_both_alleles(tmp_path: Path) -> None:
    bed_reader = pytest.importorskip("bed_reader")
    pytest.importorskip("pandas_plink")
    path = _write_plink(tmp_path)

    first = read_plink_genotypes(path, verbose=False)
    second = read_plink_genotypes(path, count_allele="allele_2", verbose=False)
    bed_first = read_plink_bed(path, verbose=False)
    bed_second = read_plink_bed(path, count_allele="allele_2", verbose=False)
    direct = bed_reader.open_bed(str(path), count_A1=False).read().tolist()

    assert first["genotypes"] == bed_first["genotypes"] == GENOTYPES
    assert second["genotypes"] == bed_second["genotypes"]
    assert second["genotypes"][0] == [2, 1, 0]
    assert direct[0] == [2.0, 1.0, 0.0]


def test_write_plink_genotypes_round_trips_through_bed_reader(tmp_path: Path) -> None:
    pytest.importorskip("pandas_plink")
    pytest.importorskip("bed_reader")

    for count_allele in ("allele_1", "allele_2"):
        path = tmp_path / f"{count_allele}.bed"
        result = write_plink_genotypes(
            path,
            GENOTYPES,
            sample_ids=SAMPLES,
            variant_ids=VARIANT_IDS,
            count_allele=count_allele,
            verbose=False,
            **PLINK_KWARGS,
        )
        resolved = path.resolve()
        assert result["files"] == [
            str(resolved.with_suffix(suffix)) for suffix in (".bed", ".bim", ".fam")
        ]
        assert read_plink_bed(path, count_allele=count_allele, verbose=False)["genotypes"] == (
            GENOTYPES
        )


def test_plink_adapters_report_missing_extra(tmp_path: Path, monkeypatch) -> None:
    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"missing{suffix}").write_bytes(b"")
    monkeypatch.setitem(__import__("sys").modules, "bed_reader", None)
    with pytest.raises(BackendUnavailableError, match="bio-bed-reader"):
        read_plink_bed(tmp_path / "missing.bed", verbose=False)


# ---------------------------------------------------------------- bgen-reader


def test_read_bgen_variants_reports_hand_written_metadata(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _two_variant_bgen(tmp_path)

    result = read_bgen_variants(path, verbose=False)

    assert result["sample_count"] == 3
    assert result["variants"][0] == {
        "id": "v1",
        "rsid": "rs1",
        "chromosome": "1",
        "position": 100,
        "alleles": ["A", "G"],
        "phased": False,
        "genotype_combinations": 3,
    }
    assert result["variants"][1]["alleles"] == ["C", "T"]


def test_read_bgen_probabilities_follows_requested_variant_order(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _two_variant_bgen(tmp_path)

    result = read_bgen_probabilities(path, variant_indices=[1, 0], verbose=False)

    assert [variant["id"] for variant in result["variants"]] == ["v2", "v1"]
    assert result["probabilities"][0] == [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    # Third sample at v1: AA = 128/255, AB = 127/255, BB = 0.
    assert result["probabilities"][1][2] == pytest.approx([128 / 255, 127 / 255, 0.0])


def test_read_bgen_probabilities_returns_none_for_missing_samples(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _write_bgen(
        tmp_path,
        "missing.bgen",
        ["S1", "S2", "S3"],
        [("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 0], [128, 127]], (1,))],
    )

    result = read_bgen_probabilities(path, verbose=False)

    assert result["probabilities"][0][1] is None
    assert result["probabilities"][0][0] == [1.0, 0.0, 0.0]


def test_calculate_bgen_allele_frequency_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _two_variant_bgen(tmp_path)

    result = calculate_bgen_allele_frequency(path, verbose=False)
    first, second = result["variants"]

    # rs1 expected A copies: 2 + 1 + 2 * 128/255 + 127/255 = 1 + 1 + 383/255 + 1.
    a_copies = 2 + 1 + 383 / 255
    g_copies = 1 + 127 / 255
    assert first["allele_copies"] == pytest.approx([a_copies, g_copies])
    assert first["allele_frequencies"] == pytest.approx([a_copies / 6, g_copies / 6])
    # rs2: C in two homozygous AA calls, T in four copies across three samples.
    assert second["allele_frequencies"] == pytest.approx([1 / 3, 2 / 3])
    assert first["called_samples"] == 3


def test_calculate_bgen_allele_frequency_excludes_missing_samples(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _write_bgen(
        tmp_path,
        "missing.bgen",
        ["S1", "S2", "S3"],
        [("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 0], [128, 127]], (1,))],
    )

    variant = calculate_bgen_allele_frequency(path, verbose=False)["variants"][0]

    assert variant["called_samples"] == 2
    assert variant["allele_frequencies"] == pytest.approx([(2 + 383 / 255) / 4, (127 / 255) / 4])


def test_calculate_bgen_allele_frequency_matches_upstream_expectation(tmp_path: Path) -> None:
    bgen_reader = pytest.importorskip("bgen_reader")
    path = _two_variant_bgen(tmp_path)

    with bgen_reader.open_bgen(str(path), verbose=False) as handle:
        expectation = handle.allele_expectation(None)
    copies = expectation.sum(axis=0)
    result = calculate_bgen_allele_frequency(path, verbose=False)

    for index, variant in enumerate(result["variants"]):
        assert variant["allele_copies"] == pytest.approx(copies[index].tolist())


def test_calculate_bgen_allele_frequency_rejects_phased_variants(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    phased = _write_bgen(
        tmp_path,
        "phased.bgen",
        ["S1", "S2"],
        [("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 255]], ())],
        phased=True,
    )

    with pytest.raises(ConfigurationError, match="unphased variants only"):
        calculate_bgen_allele_frequency(phased, verbose=False)


def test_bgen_adapters_validate_indices_and_limits(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = _two_variant_bgen(tmp_path)

    with pytest.raises(ConfigurationError, match="integers in"):
        read_bgen_probabilities(path, variant_indices=[2], verbose=False)
    with pytest.raises(ConfigurationError, match="integers in"):
        read_bgen_probabilities(path, variant_indices=[True], verbose=False)
    with pytest.raises(ConfigurationError, match="must not be empty"):
        read_bgen_probabilities(path, variant_indices=[], verbose=False)
    with pytest.raises(ConfigurationError, match="max_variants=1"):
        read_bgen_variants(path, max_variants=1, verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells=4"):
        read_bgen_probabilities(path, max_cells=4, verbose=False)


def test_bgen_adapters_report_missing_extra(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(__import__("sys").modules, "bgen_reader", None)
    path = tmp_path / "missing.bgen"
    path.write_bytes(b"\x00" * 16)
    with pytest.raises(BackendUnavailableError, match="bio-bgen-reader"):
        read_bgen_variants(path, verbose=False)


# ---------------------------------------------------------------- varcode


def test_describe_variant_reports_trimmed_alleles_and_kind() -> None:
    pytest.importorskip("varcode")

    deletion = describe_variant("1", 200, "ACG", "A", verbose=False)
    insertion = describe_variant("1", 300, "A", "AGT", verbose=False)
    snv = describe_variant("1", 101, "a", "g", verbose=False)
    substitution = describe_variant("1", 400, "AC", "TG", verbose=False)
    complex_change = describe_variant("1", 500, "AT", "G", verbose=False)

    assert (deletion["kind"], deletion["trimmed_start"], deletion["trimmed_end"]) == (
        "deletion",
        201,
        202,
    )
    assert (deletion["trimmed_ref"], deletion["trimmed_alt"]) == ("CG", "")
    assert (insertion["kind"], insertion["trimmed_start"], insertion["trimmed_end"]) == (
        "insertion",
        300,
        300,
    )
    assert insertion["trimmed_alt"] == "GT"
    assert (snv["kind"], snv["ref"], snv["alt"], snv["is_transition"]) == ("snv", "A", "G", True)
    assert substitution["kind"] == "substitution"
    assert complex_change["kind"] == "complex"
    assert complex_change["is_transition"] is None


def test_describe_variant_agrees_with_varcode_objects() -> None:
    varcode = pytest.importorskip("varcode")

    variant = varcode.Variant("1", 200, "ACG", "A")
    result = describe_variant("1", 200, "ACG", "A", verbose=False)

    assert result["trimmed_start"] == variant.trimmed_base1_start
    assert result["trimmed_end"] == variant.trimmed_base1_end
    assert result["trimmed_ref"] == variant.trimmed_ref
    assert result["trimmed_alt"] == variant.trimmed_alt


def test_describe_variant_rejects_invalid_alleles_and_positions() -> None:
    pytest.importorskip("varcode")

    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        describe_variant("1", 10, "AX", "G", verbose=False)
    with pytest.raises(ConfigurationError, match="cannot both be empty"):
        describe_variant("1", 10, "", "", verbose=False)
    with pytest.raises(ConfigurationError, match="must differ"):
        describe_variant("1", 10, "A", "a", verbose=False)
    with pytest.raises(ConfigurationError, match="1-based integer"):
        describe_variant("1", 0, "A", "G", verbose=False)
    with pytest.raises(ConfigurationError, match="1-based integer"):
        describe_variant("1", True, "A", "G", verbose=False)
    with pytest.raises(ConfigurationError, match="contig must"):
        describe_variant("chr 1", 10, "A", "G", verbose=False)


def test_load_vcf_variants_filters_and_describes_records(tmp_path: Path) -> None:
    varcode = pytest.importorskip("varcode")
    path = _write_vcf(tmp_path)

    passing = load_vcf_variants(path, verbose=False)
    everything = load_vcf_variants(path, only_passing=False, verbose=False)
    direct = varcode.load_vcf(str(path), genome="GRCh38", only_passing=False, distinct=False)

    assert passing["count"] == 4
    assert everything["count"] == 5
    assert [(v.trimmed_base1_start, v.trimmed_ref, v.trimmed_alt) for v in direct] == [
        (item["trimmed_start"], item["trimmed_ref"], item["trimmed_alt"])
        for item in everything["variants"]
    ]
    kinds = [item["kind"] for item in everything["variants"]]
    assert kinds == ["snv", "deletion", "insertion", "substitution", "snv"]


def test_load_vcf_variants_enforces_max_variants(tmp_path: Path) -> None:
    pytest.importorskip("varcode")
    path = _write_vcf(tmp_path)

    with pytest.raises(ConfigurationError, match="max_variants=2"):
        load_vcf_variants(path, max_variants=2, only_passing=False, verbose=False)
    with pytest.raises(ConfigurationError, match="max_variants must be"):
        load_vcf_variants(path, max_variants=0, verbose=False)


def test_varcode_adapters_report_missing_extra(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(__import__("sys").modules, "varcode", None)
    with pytest.raises(BackendUnavailableError, match="bio-varcode"):
        describe_variant("1", 10, "A", "G", verbose=False)


# ---------------------------------------------------------------- hgvs (parse and format only)


def test_parse_hgvs_variant_reports_coordinates_and_edit() -> None:
    pytest.importorskip("hgvs.parser")

    substitution = parse_hgvs_variant("NM_004006.2:c.4375C>T", verbose=False)
    duplication = parse_hgvs_variant("NM_004006.2:c.-5_10+3dupA", verbose=False)
    protein = parse_hgvs_variant("NP_003997.1:p.Val600Glu", verbose=False)
    genomic = parse_hgvs_variant("NC_000023.10:g.33038255C>A", verbose=False)

    assert substitution["reference"] == "NM_004006.2"
    assert substitution["coordinate_type"] == "c"
    assert substitution["start"] == {"base": 4375, "offset": 0}
    assert substitution["end"] == {"base": 4375, "offset": 0}
    assert (substitution["edit_type"], substitution["reference_bases"]) == ("sub", "C")
    assert substitution["alternate_bases"] == "T"
    assert substitution["formatted"] == "NM_004006.2:c.4375C>T"

    assert duplication["start"] == {"base": -5, "offset": 0}
    assert duplication["end"] == {"base": 10, "offset": 3}
    assert (duplication["edit_type"], duplication["reference_bases"]) == ("dup", "A")
    assert duplication["formatted"] == "NM_004006.2:c.-5_10+3dup"

    assert protein["coordinate_type"] == "p"
    assert protein["start"] == {"base": 600, "offset": None}
    assert protein["alternate_bases"] == "E"
    assert genomic["coordinate_type"] == "g"
    assert genomic["start"]["base"] == 33038255


def test_parse_hgvs_variant_matches_direct_hgvs_parser() -> None:
    parser_module = pytest.importorskip("hgvs.parser")
    text = "NM_004006.2:c.4375_4376delinsTT"

    direct = parser_module.Parser().parse_hgvs_variant(text)
    result = parse_hgvs_variant(text, verbose=False)

    assert result["formatted"] == direct.format()
    assert result["edit_type"] == direct.posedit.edit.type == "delins"
    assert result["alternate_bases"] == "TT"
    assert result["reference_bases"] is None


def test_format_hgvs_variant_reports_canonical_text() -> None:
    pytest.importorskip("hgvs.parser")

    deletion = format_hgvs_variant("NM_004006.2:c.123delA", verbose=False)
    unchanged = format_hgvs_variant("NM_004006.2:c.4375C>T", verbose=False)

    assert deletion["formatted"] == "NM_004006.2:c.123del"
    assert deletion["changed"] is True
    assert unchanged["formatted"] == "NM_004006.2:c.4375C>T"
    assert unchanged["changed"] is False


def test_parse_hgvs_variant_rejects_invalid_strings() -> None:
    pytest.importorskip("hgvs.parser")

    with pytest.raises(ConfigurationError, match="Invalid HGVS string"):
        parse_hgvs_variant("NM_000088.3:c.[589G>T;590C>A]", verbose=False)
    with pytest.raises(ConfigurationError, match="start is after its end"):
        parse_hgvs_variant("NM_004006.2:c.123_122del", verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty HGVS string"):
        parse_hgvs_variant("NM_004006.2:c.4375 C>T", verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty HGVS string"):
        parse_hgvs_variant("", verbose=False)


def test_check_hgvs_syntax_reports_each_entry() -> None:
    pytest.importorskip("hgvs.parser")
    texts = [
        "NM_004006.2:c.4375C>T",
        "NM_004006.2:c.123_122del",
        "not an hgvs string",
        "NP_003997.1:p.Val600Glu",
    ]

    result = check_hgvs_syntax(texts, verbose=False)

    assert (result["count"], result["valid"], result["invalid"]) == (4, 2, 2)
    assert [item["valid"] for item in result["items"]] == [True, False, False, True]
    assert result["items"][0]["formatted"] == "NM_004006.2:c.4375C>T"
    assert "start is after its end" in result["items"][1]["error"]
    assert result["items"][2]["formatted"] is None


def test_check_hgvs_syntax_enforces_batch_limits_and_types() -> None:
    pytest.importorskip("hgvs.parser")

    with pytest.raises(ConfigurationError, match="max_variants=1"):
        check_hgvs_syntax(
            ["NM_004006.2:c.1A>G", "NM_004006.2:c.2A>G"], max_variants=1, verbose=False
        )
    with pytest.raises(ConfigurationError, match="every entry"):
        check_hgvs_syntax(["NM_004006.2:c.1A>G", 7], verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty list"):
        check_hgvs_syntax([], verbose=False)


def test_hgvs_adapters_report_missing_extra(monkeypatch) -> None:
    monkeypatch.setitem(__import__("sys").modules, "hgvs", None)
    monkeypatch.setitem(__import__("sys").modules, "hgvs.parser", None)
    with pytest.raises(BackendUnavailableError, match="bio-hgvs"):
        parse_hgvs_variant("NM_004006.2:c.4375C>T", verbose=False)
