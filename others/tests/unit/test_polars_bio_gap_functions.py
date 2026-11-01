"""Tests for dnakit.annotation.polars_bio_gap_functions (BGEN, PGEN, VCF Zarr, CRAM, interval plot).

Offline: every input is a small file written into tmp_path. The BGEN fixture is a hand-written BGEN
v1.2 file (uncompressed, layout 2) read by polars-bio. The VCF Zarr fixture is a zarr v2 store
written with zarr. The CRAM fixture is written with pysam from unmapped reads. No PGEN writer is
installed, so the PGEN happy-path test is skipped; its guard paths still run.
"""

from __future__ import annotations

import inspect
import socket
import struct
import sys
from pathlib import Path

import pytest

pytest.importorskip("polars_bio")

from dnakit.annotation import polars_bio_gap_functions as gap
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
PGEN_SKIP = (
    "no PGEN writer in this environment: pgenlib and plink2 are not installed, "
    "so no PGEN fileset is written"
)


class _NoConnectSocket(socket.socket):
    """Refuse outgoing connections; socketpair, which asyncio uses inside zarr, still works."""

    def connect(self, *args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    def connect_ex(self, *args: object, **kwargs: object) -> int:
        raise AssertionError("network access is not allowed in this test")


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", _NoConnectSocket)


def _short(text: str) -> bytes:
    raw = text.encode("ascii")
    return struct.pack("<H", len(raw)) + raw


def write_bgen(path: Path, samples: list[str], variants: list[tuple]) -> Path:
    """Write a BGEN v1.2 file: uncompressed, layout 2, 8-bit probabilities, with sample identifiers.

    Each variant is (rsid, chrom, pos, allele_a, allele_b, stored) where stored gives, per sample,
    the integers (p_AA, p_AB) out of 255; p_BB is 255 minus their sum.
    """

    n_samples = len(samples)
    header_length = 20
    flags = 0x80000000 | (2 << 2)
    header = struct.pack("<I", header_length) + struct.pack("<II", len(variants), n_samples)
    header += b"bgen" + struct.pack("<I", flags)
    ids = struct.pack("<I", n_samples) + b"".join(_short(sample) for sample in samples)
    sample_block = struct.pack("<I", 4 + len(ids)) + ids
    body = b""
    for index, (rsid, chrom, pos, allele_a, allele_b, stored) in enumerate(variants):
        block = _short(f"var{index}") + _short(rsid) + _short(chrom) + struct.pack("<I", pos)
        block += struct.pack("<H", 2)
        block += struct.pack("<I", len(allele_a)) + allele_a.encode("ascii")
        block += struct.pack("<I", len(allele_b)) + allele_b.encode("ascii")
        geno = struct.pack("<IHBB", n_samples, 2, 2, 2) + bytes([2] * n_samples) + bytes([0, 8])
        for p_aa, p_ab in stored:
            geno += bytes([p_aa, p_ab])
        block += struct.pack("<I", len(geno)) + geno
        body += block
    data = struct.pack("<I", header_length + len(sample_block)) + header + sample_block + body
    path.write_bytes(data)
    return path


@pytest.fixture()
def bgen(tmp_path: Path) -> Path:
    # var0: S1 AA, S2 AB, S3 BB -> dosages 0, 1, 2. var1: S1 BB, S2 AA, S3 AB -> dosages 2, 0, 1.
    return write_bgen(
        tmp_path / "cohort.bgen",
        ["S1", "S2", "S3"],
        [
            ("rs1", "1", 100, "A", "G", [(255, 0), (0, 255), (0, 0)]),
            ("rs2", "1", 250, "C", "T", [(0, 0), (255, 0), (0, 255)]),
        ],
    )


def write_vcf_zarr(path: Path) -> Path:
    """Write a two-variant, two-sample VCF Zarr store (zarr format 2 arrays)."""

    zarr = pytest.importorskip("zarr")
    import numpy as np
    from zarr.core.dtype import VariableLengthUTF8

    root = zarr.open_group(str(path), mode="w", zarr_format=2)
    root.attrs["vcf_zarr_version"] = "0.4"

    def array(name: str, data: object, dtype: object = None) -> None:
        values = np.asarray(data, dtype=dtype)
        zdtype = VariableLengthUTF8() if values.dtype == object else values.dtype
        stored = root.create_array(name, shape=values.shape, dtype=zdtype, chunks=values.shape)
        stored[...] = values

    array("contig_id", np.array(["chr1"], dtype=object), object)
    array("sample_id", np.array(["S1", "S2"], dtype=object), object)
    array("variant_contig", np.array([0, 0], dtype=np.int32))
    array("variant_position", np.array([100, 250], dtype=np.int32))
    array("variant_id", np.array(["rs1", "rs2"], dtype=object), object)
    array("variant_allele", np.array([["A", "G"], ["C", "T"]], dtype=object), object)
    array("variant_quality", np.array([50.0, 60.0], dtype=np.float32))
    array("variant_filter", np.zeros((2, 1), dtype=bool))
    array("filter_id", np.array(["PASS"], dtype=object), object)
    genotypes = np.array([[[0, 0], [0, 1]], [[1, 1], [0, 0]]], dtype=np.int8)
    array("call_genotype", genotypes)
    array("call_genotype_mask", genotypes < 0)
    array("call_genotype_phased", np.zeros((2, 2), dtype=bool))
    return path


@pytest.fixture()
def vcf_zarr(tmp_path: Path) -> Path:
    return write_vcf_zarr(tmp_path / "cohort.vcz")


def write_unmapped_cram(path: Path) -> Path:
    """Write two unmapped reads as CRAM; htslib embeds the sequences, so no reference is needed."""

    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wc", header=header) as out:
        for index, sequence in enumerate(["ACGTACGTAC", "GGGGCCCCAA"]):
            read = pysam.AlignedSegment()
            read.query_name = f"r{index}"
            read.query_sequence = sequence
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            read.flag = 4
            read.reference_id = -1
            read.reference_start = -1
            read.mapping_quality = 0
            out.write(read)
    return path


def test_public_functions_take_keyword_only_verbose() -> None:
    assert gap.__all__
    for name in gap.__all__:
        function = getattr(gap, name)
        assert callable(function), name
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name
        assert list(inspect.signature(function).parameters)[-1] == "verbose", name


# BGEN.


def test_read_bgen_returns_probabilities_per_variant(bgen: Path) -> None:
    result = gap.polars_bio_read_bgen(str(bgen), verbose=False)
    assert result["variant_count"] == 2
    first, second = result["variants"]
    assert first["id"] == "var0" and first["rsid"] == "rs1"
    assert first["start"] == 100 and first["alleles"] == ["A", "G"]
    assert first["genotypes"]["GP"] == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    assert second["genotypes"]["GP"][0] == [0.0, 0.0, 1.0]


def test_read_bgen_dosage_output(bgen: Path) -> None:
    result = gap.polars_bio_read_bgen(str(bgen), genotype_output="dosage", verbose=False)
    assert result["genotype_output"] == "dosage"
    assert result["variants"][0]["genotypes"]["DS"] == [0.0, 1.0, 2.0]
    assert result["variants"][1]["genotypes"]["DS"] == [2.0, 0.0, 1.0]


def test_read_bgen_matrix_gives_alt_dosages(bgen: Path) -> None:
    result = gap.polars_bio_read_bgen_matrix(str(bgen), verbose=False)
    assert result["positions"] == [100, 250]
    assert result["sample_names"] == ["S1", "S2", "S3"]
    assert result["values"] == [[0.0, 1.0, 2.0], [2.0, 0.0, 1.0]]
    assert (result["variant_count"], result["sample_count"]) == (2, 3)


def test_describe_bgen_lists_the_output_columns(bgen: Path) -> None:
    result = gap.polars_bio_describe_bgen(str(bgen), verbose=False)
    names = [entry["name"] for entry in result["fields"]]
    assert names == [
        "chrom",
        "start",
        "end",
        "id",
        "rsid",
        "alleles",
        "phased",
        "bits",
        "genotypes",
    ]
    assert result["field_count"] == 9


def test_scan_bgen_counts_all_variants_and_materializes_the_head(bgen: Path) -> None:
    result = gap.polars_bio_scan_bgen(str(bgen), max_rows=1, verbose=False)
    assert result["variant_count"] == 2
    assert result["returned_count"] == 1
    assert result["truncated"] is True
    assert result["variants"][0]["id"] == "var0"
    assert [column["name"] for column in result["columns"]][:2] == ["chrom", "start"]


def test_register_bgen_queries_the_table_and_removes_it(bgen: Path) -> None:
    pb = pytest.importorskip("polars_bio")
    result = gap.polars_bio_register_bgen(str(bgen), "gap_bgen_table", verbose=False)
    assert result["table_name"] == "gap_bgen_table"
    assert result["row_count"] == 2
    assert "chrom" in result["columns"]
    assert result["rows"][0]["id"] == "var0"
    with pytest.raises(KeyError):
        pb.ctx.table("gap_bgen_table")


def test_register_bgen_refuses_a_table_name_already_in_the_session(bgen: Path) -> None:
    pb = pytest.importorskip("polars_bio")
    pb.register_bgen(str(bgen), "gap_clash_table")
    try:
        with pytest.raises(ConfigurationError, match=r"already registered"):
            gap.polars_bio_register_bgen(str(bgen), "gap_clash_table", verbose=False)
        assert pb.ctx.table("gap_clash_table") is not None
    finally:
        pb.ctx.deregister_table("gap_clash_table")


def test_bgen_guards(tmp_path: Path, bgen: Path) -> None:
    with pytest.raises(ConfigurationError):
        gap.polars_bio_read_bgen(str(tmp_path / "missing.bgen"), verbose=False)
    renamed = tmp_path / "cohort.bgn"
    renamed.write_bytes(bgen.read_bytes())
    with pytest.raises(ConfigurationError, match=r"ending in \.bgen"):
        gap.polars_bio_read_bgen(str(renamed), verbose=False)
    with pytest.raises(ConfigurationError, match=r"genotype_output"):
        gap.polars_bio_read_bgen(str(bgen), genotype_output="counts", verbose=False)
    with pytest.raises(ConfigurationError, match=r"max_rows"):
        gap.polars_bio_read_bgen(str(bgen), max_rows=1, verbose=False)
    with pytest.raises(ConfigurationError, match=r"table_name"):
        gap.polars_bio_register_bgen(str(bgen), "Bad-Name", verbose=False)
    with pytest.raises(ConfigurationError, match=r"remote URLs"):
        gap.polars_bio_read_bgen("https://example.org/cohort.bgen", verbose=False)


def test_missing_polars_bio_is_a_backend_error(monkeypatch: pytest.MonkeyPatch, bgen: Path) -> None:
    monkeypatch.setitem(sys.modules, "polars_bio", None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        gap.polars_bio_read_bgen(str(bgen), verbose=False)
    assert "bio-polarsbio" in excinfo.value.hint


# PGEN: guard paths always run; the happy path needs a PGEN writer, which is not installed.


def test_pgen_guards(tmp_path: Path, bgen: Path) -> None:
    placeholder = tmp_path / "cohort.pgen"
    placeholder.write_bytes(b"")
    with pytest.raises(ConfigurationError, match=r"ending in \.pgen"):
        gap.polars_bio_read_pgen(str(bgen), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.polars_bio_describe_pgen(str(tmp_path / "missing.pgen"), verbose=False)
    with pytest.raises(ConfigurationError, match=r"psam_id_mode"):
        gap.polars_bio_read_pgen(str(placeholder), psam_id_mode="sample", verbose=False)
    with pytest.raises(ConfigurationError, match=r"field"):
        gap.polars_bio_read_pgen_matrix(str(placeholder), field="GT", verbose=False)


@pytest.mark.skip(reason=PGEN_SKIP)
def test_pgen_happy_path_reads_a_pgen_fileset() -> None:  # pragma: no cover - needs a PGEN fileset
    """Would check the five PGEN readers on a real .pgen/.pvar/.psam fileset (two variants)."""
    pytest.fail("a PGEN fileset is required for this check")


# VCF Zarr.


def test_read_vcf_zarr_returns_one_record_per_variant(vcf_zarr: Path) -> None:
    result = gap.polars_bio_read_vcf_zarr(str(vcf_zarr), verbose=False)
    assert result["variant_count"] == 2
    assert [row["id"] for row in result["variants"]] == ["rs1", "rs2"]
    assert [row["alt"] for row in result["variants"]] == ["G", "T"]
    assert [row["start"] for row in result["variants"]] == [100, 250]


def test_describe_vcf_zarr_lists_the_schema(vcf_zarr: Path) -> None:
    result = gap.polars_bio_describe_vcf_zarr(str(vcf_zarr), verbose=False)
    assert result["field_count"] > 0
    assert result["fields"]


def test_scan_vcf_zarr_counts_and_truncates(vcf_zarr: Path) -> None:
    result = gap.polars_bio_scan_vcf_zarr(str(vcf_zarr), max_rows=1, verbose=False)
    assert result["variant_count"] == 2
    assert result["returned_count"] == 1
    assert result["truncated"] is True


def test_register_vcf_zarr_queries_the_table(vcf_zarr: Path) -> None:
    result = gap.polars_bio_register_vcf_zarr(str(vcf_zarr), "gap_vcz_table", verbose=False)
    assert result["row_count"] == 2
    assert result["rows"][1]["id"] == "rs2"


def test_vcf_zarr_guards(tmp_path: Path, bgen: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"directory"):
        gap.polars_bio_read_vcf_zarr(str(bgen), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.polars_bio_read_vcf_zarr(str(tmp_path / "absent.vcz"), verbose=False)
    with pytest.raises(ConfigurationError, match=r"max_rows"):
        gap.polars_bio_read_vcf_zarr(str(tmp_path), max_rows=0, verbose=False)


# CRAM.


def test_register_cram_reads_embedded_sequences(tmp_path: Path) -> None:
    cram = write_unmapped_cram(tmp_path / "reads.cram")
    result = gap.polars_bio_register_cram(str(cram), "gap_cram_table", verbose=False)
    assert result["row_count"] == 2
    assert [row["name"] for row in result["rows"]] == ["r0", "r1"]
    assert result["rows"][1]["sequence"] == "GGGGCCCCAA"


def test_register_cram_guards(tmp_path: Path) -> None:
    cram = write_unmapped_cram(tmp_path / "reads.cram")
    with pytest.raises(ConfigurationError):
        gap.polars_bio_register_cram(
            str(tmp_path / "absent.cram"), "gap_cram_missing", verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"table_name"):
        gap.polars_bio_register_cram(str(cram), "1bad", verbose=False)


# Interval plot.


def _write_bed(path: Path, rows: list[str]) -> Path:
    path.write_text("".join(f"{row}\n" for row in rows), encoding="utf-8")
    return path


def _colour_present(image: object, rgb: tuple[int, int, int], tolerance: float = 0.06) -> bool:
    """True when a pixel of the RGB image is within ``tolerance`` of the colour ``rgb`` (0-255)."""

    import numpy as np

    target = np.array(rgb, dtype=float) / 255.0
    pixels = np.asarray(image)[..., :3]
    return bool((np.abs(pixels - target).max(axis=-1) < tolerance).any())


def test_visualize_intervals_draws_each_overlapping_pair(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    left = _write_bed(tmp_path / "left.bed", ["chr1\t10\t20\ta", "chr1\t40\t60\tb"])
    right = _write_bed(tmp_path / "right.bed", ["chr1\t15\t30\tc", "chr1\t50\t55\td"])
    target = tmp_path / "plots" / "overlaps.png"
    result = gap.polars_bio_visualize_intervals(str(left), str(right), str(target), verbose=False)
    assert result["path"] == str(target.resolve())
    assert target.read_bytes()[:8] == PNG_SIGNATURE
    assert plt.get_fignums() == []
    # polars-bio colours the left interval of each pair skyblue and the right one lightpink.
    image = imread(str(target))
    assert _colour_present(image, (135, 206, 235))
    assert _colour_present(image, (255, 182, 193))
    data = result["data"]
    assert data["pair_count"] == 2 and data["plotted_count"] == 2
    assert [
        (pair["start_1"], pair["end_1"], pair["start_2"], pair["end_2"]) for pair in data["pairs"]
    ] == [
        (10, 20, 15, 30),
        (40, 60, 50, 55),
    ]


def test_visualize_intervals_respects_max_pairs(tmp_path: Path) -> None:
    left = _write_bed(tmp_path / "left.bed", ["chr1\t10\t20\ta", "chr1\t40\t60\tb"])
    right = _write_bed(tmp_path / "right.bed", ["chr1\t15\t30\tc", "chr1\t50\t55\td"])
    result = gap.polars_bio_visualize_intervals(
        str(left), str(right), str(tmp_path / "one.png"), label="pair", max_pairs=1, verbose=False
    )
    assert result["data"]["pair_count"] == 2
    assert result["data"]["plotted_count"] == 1
    assert result["data"]["label"] == "pair"


def test_visualize_intervals_guards(tmp_path: Path) -> None:
    left = _write_bed(tmp_path / "left.bed", ["chr1\t10\t20\ta"])
    far = _write_bed(tmp_path / "far.bed", ["chr1\t100\t120\tz"])
    with pytest.raises(ConfigurationError, match=r"no overlapping intervals"):
        gap.polars_bio_visualize_intervals(
            str(left), str(far), str(tmp_path / "none.png"), verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"\.png"):
        gap.polars_bio_visualize_intervals(
            str(left), str(left), str(tmp_path / "plot.pdf"), verbose=False
        )
    existing = tmp_path / "exists.png"
    existing.write_bytes(PNG_SIGNATURE)
    with pytest.raises(ConfigurationError, match=r"already exists"):
        gap.polars_bio_visualize_intervals(str(left), str(left), str(existing), verbose=False)
