"""Contract tests for G6b wrappers: pandas-plink relationship matrices, bed-reader genotype
matrices and bgen-reader allele expectations, dosages and metafiles.

Every BGEN and PLINK fixture is written inside tmp_path. Each wrapper is checked against hand
calculations or a direct call to the upstream package. No test touches the network.
"""

from __future__ import annotations

import math
import struct
import sys
import zlib
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants.bed_reader_functions import (
    calculate_genotype_relationship_matrix,
    standardize_genotypes,
    subset_genotype_matrix,
)
from dnakit.variants.bgen_reader_functions import (
    create_bgen_metafile,
    read_bgen_allele_expectation,
    read_bgen_dosage,
)
from dnakit.variants.pandas_plink_functions import (
    read_gcta_relationship_matrix,
    read_plink_relationship_matrix,
)

REL_MATRIX = [
    [1.0, 0.5, 0.25],
    [0.5, 2.0, -0.5],
    [0.25, -0.5, 1.5],
]
SAMPLES = [("fam1", "ind1"), ("fam1", "ind2"), ("fam2", "ind3")]


def _write_id_file(path: Path) -> None:
    path.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in SAMPLES), encoding="utf-8")


def _write_bgen(path: Path, samples: list[str], variants: list[dict], bits: int = 8) -> None:
    """Write a zlib-compressed BGEN 1.2 file (layout 2, diploid, biallelic) for tests.

    Each variant gives one entry per sample: None for a missing sample, otherwise the unphased
    probabilities of genotypes AA, AB, BB, or the two phased haplotype probabilities when the
    variant is marked phased.
    """

    scale = (1 << bits) - 1
    sample_block = struct.pack("<I", len(samples)) + b"".join(
        struct.pack("<H", len(sample)) + sample.encode() for sample in samples
    )
    sample_block = struct.pack("<I", len(sample_block) + 4) + sample_block
    header = struct.pack("<I", 20) + struct.pack("<II", len(variants), len(samples)) + b"bgen"
    header += struct.pack("<I", 1 | (2 << 2) | (1 << 31))  # zlib, layout 2, sample IDs present
    out = bytearray(struct.pack("<I", len(header) + len(sample_block)))
    out += header + sample_block
    for variant in variants:
        block = bytearray()
        for text in (variant["id"], variant["rsid"], variant["chrom"]):
            block += struct.pack("<H", len(text)) + text.encode()
        block += struct.pack("<IH", variant["pos"], 2)
        for allele in variant["alleles"]:
            block += struct.pack("<I", len(allele)) + allele.encode()
        data = bytearray(struct.pack("<IHBB", len(samples), 2, 2, 2))
        for probs in variant["probs"]:
            data.append(0x80 | 2 if probs is None else 2)
        data += bytes([1 if variant.get("phased") else 0, bits])
        for probs in variant["probs"]:
            values = [0.0, 0.0] if probs is None else list(probs[:2])
            data += bytes(round(value * scale) for value in values)
        packed = zlib.compress(bytes(data))
        block += struct.pack("<II", len(packed) + 4, len(data)) + packed
        out += block
    path.write_bytes(bytes(out))


# Probabilities are k/255 so that 8-bit storage is exact and hand calculations are exact.
BGEN_VARIANTS = [
    {
        "id": "v1",
        "rsid": "rs1",
        "chrom": "1",
        "pos": 100,
        "alleles": ["A", "G"],
        "probs": [[51 / 255, 76 / 255, 128 / 255], [1.0, 0.0, 0.0], None],
    },
    {
        "id": "v2",
        "rsid": "rs2",
        "chrom": "1",
        "pos": 200,
        "alleles": ["C", "T"],
        "probs": [[0.0, 0.0, 1.0], [64 / 255, 128 / 255, 63 / 255], None],
    },
]
BGEN_SAMPLES = ["s1", "s2", "s3"]


@pytest.fixture()
def bgen_file(tmp_path: Path) -> Path:
    path = tmp_path / "small.bgen"
    _write_bgen(path, BGEN_SAMPLES, BGEN_VARIANTS)
    return path


def _write_bed(path: Path, values: list[list[float]], names: list[str]) -> None:
    import numpy as np
    from bed_reader import to_bed

    count = len(names)
    to_bed(
        str(path),
        np.array(values, dtype=np.float64),
        properties={
            "iid": [f"i{index}" for index in range(len(values))],
            "sid": [f"s{index}" for index in range(count)],
            "chromosome": ["1"] * count,
            "bp_position": list(range(1, count + 1)),
            "allele_1": ["A"] * count,
            "allele_2": ["G"] * count,
        },
        count_A1=True,
    )


# ---------------------------------------------------------------------------- pandas-plink


def _gcta_binary(tmp_path: Path) -> Path:
    import numpy as np

    lower = [REL_MATRIX[row][column] for row in range(3) for column in range(row + 1)]
    prefix = tmp_path / "study.grm"
    np.array(lower, dtype=np.float32).tofile(prefix.with_suffix(".grm.bin"))
    np.array([10] * len(lower), dtype=np.float32).tofile(tmp_path / "study.grm.N.bin")
    _write_id_file(tmp_path / "study.grm.id")
    return prefix.with_suffix(".grm.bin")


def test_gcta_binary_grm_matches_hand_values_and_upstream(tmp_path: Path) -> None:
    pandas_plink = pytest.importorskip("pandas_plink")
    path = _gcta_binary(tmp_path)

    result = read_gcta_relationship_matrix(path, verbose=False)

    assert result["matrix"] == REL_MATRIX
    assert result["samples"] == [{"family_id": fid, "sample_id": iid} for fid, iid in SAMPLES]
    assert result["n_snps"] == [10] * 6
    kinship, n_snps = pandas_plink.read_grm(str(path))
    assert kinship.values.tolist() == result["matrix"]
    assert list(n_snps) == result["n_snps"]


def test_gcta_text_grm_matches_binary_and_infers_id_file(tmp_path: Path) -> None:
    pytest.importorskip("pandas_plink")
    binary = _gcta_binary(tmp_path)
    text = tmp_path / "text.grm"
    lines = [
        "1\t1\t10\t1.0",
        "2\t1\t10\t0.5",
        "2\t2\t10\t2.0",
        "3\t1\t10\t0.25",
        "3\t2\t10\t-0.5",
        "3\t3\t10\t1.5",
    ]
    text.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _write_id_file(tmp_path / "text.grm.id")

    result = read_gcta_relationship_matrix(text, verbose=False)
    reference = read_gcta_relationship_matrix(binary, verbose=False)

    assert result["matrix"] == REL_MATRIX
    assert result["matrix"] == reference["matrix"]
    assert result["n_snps"] == [10] * 6


def test_plink_rel_binary_and_text_match_upstream(tmp_path: Path) -> None:
    pandas_plink = pytest.importorskip("pandas_plink")
    import numpy as np

    binary = tmp_path / "study.rel.bin"
    np.array(REL_MATRIX, dtype=np.float64).tofile(binary)
    _write_id_file(tmp_path / "study.rel.id")
    text = tmp_path / "study_text.rel"
    text.write_text(
        "1.0\n0.5\t2.0\n0.25\t-0.5\t1.5\n",
        encoding="utf-8",
    )
    _write_id_file(tmp_path / "study_text.rel.id")

    from_binary = read_plink_relationship_matrix(binary, verbose=False)
    from_text = read_plink_relationship_matrix(text, verbose=False)

    assert from_binary["matrix"] == REL_MATRIX
    assert from_text["matrix"] == REL_MATRIX
    assert from_binary["samples"] == from_text["samples"]
    upstream = pandas_plink.read_rel(str(binary))
    assert upstream.values.tolist() == from_binary["matrix"]


def test_relationship_readers_bound_samples_and_check_inputs(tmp_path: Path) -> None:
    pytest.importorskip("pandas_plink")
    path = _gcta_binary(tmp_path)

    # max_cells=4 allows at most a 2 x 2 matrix; the file holds 3 samples.
    with pytest.raises(ConfigurationError, match="max_cells"):
        read_gcta_relationship_matrix(path, max_cells=4, verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells"):
        read_gcta_relationship_matrix(path, max_cells=0, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        read_gcta_relationship_matrix(tmp_path / "missing.grm.bin", verbose=False)
    (tmp_path / "study.grm.id").unlink()
    with pytest.raises(ConfigurationError, match="sample file does not exist"):
        read_gcta_relationship_matrix(path, verbose=False)


def test_relationship_reader_reports_missing_pandas_plink(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import numpy as np

    binary = tmp_path / "study.rel.bin"
    np.array(REL_MATRIX, dtype=np.float64).tofile(binary)
    _write_id_file(tmp_path / "study.rel.id")
    monkeypatch.setitem(sys.modules, "pandas_plink", None)
    with pytest.raises(BackendUnavailableError, match="bio-pandas-plink"):
        read_plink_relationship_matrix(binary, verbose=False)


# ---------------------------------------------------------------------------- bed-reader


def test_standardize_matches_hand_calculation_and_upstream(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    import numpy as np
    from bed_reader import standardize_f64

    rows = [[0, 1], [1, 1], [2, 1], [None, 2]]
    # Column 0 observed values 0, 1, 2: mean 1, population SD sqrt(2/3).
    # Column 1 observed values 1, 1, 1, 2: mean 1.25, population SD sqrt(3/16).
    sd0 = math.sqrt(2 / 3)
    sd1 = math.sqrt(3 / 16)
    expected = [
        [-1 / sd0, -0.25 / sd1],
        [0.0, -0.25 / sd1],
        [1 / sd0, -0.25 / sd1],
        [0.0, 0.75 / sd1],
    ]

    result = standardize_genotypes(rows, verbose=False)

    assert result["dtype"] == "float64"
    assert result["sample_count"] == 4 and result["variant_count"] == 2
    for got_row, want_row in zip(result["standardized"], expected, strict=True):
        assert got_row == pytest.approx(want_row, abs=1e-12)
    assert result["variant_statistics"][0]["mean"] == pytest.approx(1.0)
    assert result["variant_statistics"][0]["std"] == pytest.approx(sd0, abs=1e-12)
    assert result["variant_statistics"][1]["mean"] == pytest.approx(1.25)
    assert result["variant_statistics"][1]["std"] == pytest.approx(sd1, abs=1e-12)

    # Direct upstream call on the same matrix: in-place kernel with stats output.
    values = np.asfortranarray(np.array([[0, 1], [1, 1], [2, 1], [np.nan, 2]], dtype=np.float64))
    stats = np.zeros((2, 2), order="F")
    standardize_f64(values, False, 1.0, 25.0, True, False, stats, 1)
    assert values.tolist() == result["standardized"]


def test_standardize_float32_variant_matches_float64(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    rows = [[0, 1], [1, 1], [2, 1], [None, 2]]
    single = standardize_genotypes(rows, dtype="float32", verbose=False)
    double = standardize_genotypes(rows, dtype="float64", verbose=False)

    assert single["dtype"] == "float32"
    for got_row, want_row in zip(single["standardized"], double["standardized"], strict=True):
        assert got_row == pytest.approx(want_row, abs=1e-6)


def test_standardize_rejects_invalid_genotypes(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("bed_reader")
    with pytest.raises(ConfigurationError, match="no observed genotypes"):
        standardize_genotypes([[1, None], [2, None]], verbose=False)
    with pytest.raises(ConfigurationError, match="dtype"):
        standardize_genotypes([[1, 2]], dtype="int8", verbose=False)
    with pytest.raises(ConfigurationError, match="finite numbers"):
        standardize_genotypes([[True, 2]], verbose=False)
    with pytest.raises(ConfigurationError, match="variant columns"):
        standardize_genotypes([[1, 2], [1]], verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells"):
        standardize_genotypes([[1, 2], [2, 1]], max_cells=3, verbose=False)
    with pytest.raises(ConfigurationError, match="finite numbers"):
        standardize_genotypes([[float("inf"), 1.0]], verbose=False)


def test_subset_matches_hand_selection_and_upstream() -> None:
    pytest.importorskip("bed_reader")
    import numpy as np
    from bed_reader import subset_f64_f64

    matrix = [[0, 1, 2], [None, 1, 0], [2, 2, 2]]

    result = subset_genotype_matrix(
        matrix, sample_indices=[0, 2], variant_indices=[1, 2], verbose=False
    )
    assert result["genotypes"] == [[1.0, 2.0], [2.0, 2.0]]
    assert result["sample_indices"] == [0, 2] and result["variant_indices"] == [1, 2]

    # Rows and columns follow the requested order and indices may repeat.
    repeated = subset_genotype_matrix(
        matrix, sample_indices=[2, 2], variant_indices=[0, 2], verbose=False
    )
    assert repeated["genotypes"] == [[2.0, 2.0], [2.0, 2.0]]
    missing = subset_genotype_matrix(
        matrix, sample_indices=[1], variant_indices=[0, 1], verbose=False
    )
    assert missing["genotypes"] == [[None, 1.0]]

    source = np.asfortranarray(np.array(matrix, dtype=np.float64)[:, :, None])
    target = np.zeros((2, 2, 1), order="F")
    subset_f64_f64(
        source,
        np.array([0, 2], dtype=np.uintp),
        np.array([1, 2], dtype=np.uintp),
        target,
        1,
    )
    assert target[:, :, 0].tolist() == result["genotypes"]

    single = subset_genotype_matrix(
        matrix, sample_indices=[0, 2], variant_indices=[1, 2], dtype="float32", verbose=False
    )
    assert single["genotypes"] == result["genotypes"]


def test_subset_rejects_bad_indices() -> None:
    pytest.importorskip("bed_reader")
    matrix = [[0, 1, 2], [1, 1, 0], [2, 2, 2]]
    with pytest.raises(ConfigurationError, match="sample_indices"):
        subset_genotype_matrix(matrix, sample_indices=[3], variant_indices=[0], verbose=False)
    with pytest.raises(ConfigurationError, match="variant_indices"):
        subset_genotype_matrix(matrix, sample_indices=[0], variant_indices=[True], verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty"):
        subset_genotype_matrix(matrix, sample_indices=[], variant_indices=[0], verbose=False)


def test_relationship_matrix_from_bed_matches_manual_formula(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    bed = tmp_path / "panel.bed"
    columns = [
        [0, 1, 2, None],  # polymorphic with one missing call
        [1, 1, 0, 2],
        [2, None, 1, 0],
        [0, 0, 2, 1],
        [1, 1, 1, 1],  # monomorphic: observed, but contributes nothing
        [None, None, None, None],  # unobserved: skipped
    ]
    values = [[columns[c][r] for c in range(len(columns))] for r in range(4)]
    values = [[math.nan if v is None else float(v) for v in row] for row in values]
    _write_bed(bed, values, names=["a", "b", "c", "d", "e", "f"])

    result = calculate_genotype_relationship_matrix(bed, verbose=False)

    # Manual standardization: population SD, missing and monomorphic columns become zero.
    z = [[0.0] * 5 for _ in range(4)]
    polymorphic = 0
    for c in range(5):
        observed = [values[r][c] for r in range(4) if not math.isnan(values[r][c])]
        mean = sum(observed) / len(observed)
        sd = math.sqrt(sum((x - mean) ** 2 for x in observed) / len(observed))
        if sd > 0:
            polymorphic += 1
            for r in range(4):
                if not math.isnan(values[r][c]):
                    z[r][c] = (values[r][c] - mean) / sd
    expected = [
        [sum(z[i][k] * z[j][k] for k in range(5)) / polymorphic for j in range(4)] for i in range(4)
    ]

    assert result["variant_count"] == 6
    assert result["observed_variants"] == 5
    assert result["polymorphic_variants"] == polymorphic == 4
    assert len(result["samples"]) == 4
    for got_row, want_row in zip(result["matrix"], expected, strict=True):
        assert got_row == pytest.approx(want_row, abs=1e-12)
    matrix = result["matrix"]
    assert all(matrix[i][j] == pytest.approx(matrix[j][i]) for i in range(4) for j in range(4))


def test_relationship_matrix_bounds_and_missing_files(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    bed = tmp_path / "small.bed"
    _write_bed(bed, [[0.0, 1.0], [1.0, 2.0], [2.0, 0.0]], names=["a", "b"])
    with pytest.raises(ConfigurationError, match="max_cells"):
        calculate_genotype_relationship_matrix(bed, max_cells=3, verbose=False)
    (tmp_path / "small.fam").unlink()
    with pytest.raises(ConfigurationError, match="does not exist"):
        calculate_genotype_relationship_matrix(bed, verbose=False)


def test_bed_wrappers_report_missing_bed_reader(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "bed_reader", None)
    with pytest.raises(BackendUnavailableError, match="bio-bed-reader"):
        standardize_genotypes([[0, 1], [2, 1]], verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-bed-reader"):
        subset_genotype_matrix([[0, 1]], sample_indices=[0], variant_indices=[0], verbose=False)


# ---------------------------------------------------------------------------- bgen-reader


def test_allele_expectation_matches_hand_values(bgen_file: Path) -> None:
    pytest.importorskip("bgen_reader")
    result = read_bgen_allele_expectation(bgen_file, verbose=False)

    assert result["samples"] == BGEN_SAMPLES
    assert [variant["id"] for variant in result["variants"]] == ["v1", "v2"]
    v1, v2 = result["expectations"]
    # v1, sample s1: AA=51/255, AG=76/255, GG=128/255 -> A copies = 2*51/255 + 76/255.
    assert v1[0] == pytest.approx([178 / 255, 332 / 255], abs=1e-12)
    assert v1[1] == pytest.approx([2.0, 0.0], abs=1e-12)
    assert v1[2] is None  # missing sample
    assert v2[0] == pytest.approx([0.0, 2.0], abs=1e-12)
    assert v2[1] == pytest.approx([256 / 255, 254 / 255], abs=1e-12)


def test_allele_expectation_and_dosage_match_upstream_functions(bgen_file: Path) -> None:
    bgen_reader = pytest.importorskip("bgen_reader")
    import numpy as np

    expectation = read_bgen_allele_expectation(bgen_file, verbose=False)
    dosage = read_bgen_dosage(bgen_file, verbose=False)

    upstream = bgen_reader.read_bgen(str(bgen_file), verbose=False)
    for column in range(2):
        direct = np.asarray(bgen_reader.allele_expectation(upstream, column), dtype=float)
        ours = np.array(
            [[np.nan] * 2 if row is None else row for row in expectation["expectations"][column]],
            dtype=float,
        )
        assert np.allclose(ours, direct, equal_nan=True)
        direct_dosage = np.asarray(bgen_reader.compute_dosage(direct, alt=1), dtype=float)
        ours_dosage = np.array(
            [np.nan if value is None else value for value in dosage["dosages"][column]],
            dtype=float,
        )
        assert np.allclose(ours_dosage, direct_dosage, equal_nan=True)


def test_dosage_default_and_explicit_alternative_allele(bgen_file: Path) -> None:
    pytest.importorskip("bgen_reader")
    default = read_bgen_dosage(bgen_file, verbose=False)
    first = read_bgen_dosage(bgen_file, alt_allele_index=0, verbose=False)

    # Default is the last allele column (G for v1, T for v2).
    assert [variant["alt_allele"] for variant in default["variants"]] == ["G", "T"]
    assert default["dosages"][0][0] == pytest.approx(332 / 255, abs=1e-12)
    assert default["dosages"][0][2] is None
    assert first["variants"][0]["alt_allele_index"] == 0
    assert first["dosages"][0][0] == pytest.approx(178 / 255, abs=1e-12)
    assert first["dosages"][0][1] == pytest.approx(2.0, abs=1e-12)

    selected = read_bgen_dosage(bgen_file, variant_indices=[1], verbose=False)
    assert [variant["id"] for variant in selected["variants"]] == ["v2"]
    assert selected["dosages"][0][0] == pytest.approx(2.0, abs=1e-12)


def test_bgen_rejects_bad_alleles_bounds_and_suffixes(bgen_file: Path, tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    with pytest.raises(ConfigurationError, match="out of range"):
        read_bgen_dosage(bgen_file, alt_allele_index=2, verbose=False)
    with pytest.raises(ConfigurationError, match="alt_allele_index"):
        read_bgen_dosage(bgen_file, alt_allele_index=True, verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells"):
        read_bgen_allele_expectation(bgen_file, max_cells=3, verbose=False)
    with pytest.raises(ConfigurationError, match="above max_variants"):
        read_bgen_allele_expectation(bgen_file, max_variants=1, verbose=False)
    with pytest.raises(ConfigurationError, match="variant_indices"):
        read_bgen_allele_expectation(bgen_file, variant_indices=[5], verbose=False)
    renamed = tmp_path / "small.bgn"
    renamed.write_bytes(bgen_file.read_bytes())
    with pytest.raises(ConfigurationError, match=r"\.bgen"):
        read_bgen_allele_expectation(renamed, verbose=False)


def test_phased_variants_are_rejected(tmp_path: Path) -> None:
    pytest.importorskip("bgen_reader")
    path = tmp_path / "phased.bgen"
    _write_bgen(
        path,
        BGEN_SAMPLES,
        [
            {
                "id": "p1",
                "rsid": "rsp",
                "chrom": "1",
                "pos": 300,
                "alleles": ["A", "G"],
                "phased": True,
                "probs": [[0.5, 0.0], [1.0, 1.0], None],
            }
        ],
    )
    with pytest.raises(ConfigurationError, match="unphased variants only"):
        read_bgen_allele_expectation(path, verbose=False)
    with pytest.raises(ConfigurationError, match="unphased variants only"):
        read_bgen_dosage(path, verbose=False)


def test_create_bgen_metafile_round_trip(bgen_file: Path, tmp_path: Path) -> None:
    bgen_reader = pytest.importorskip("bgen_reader")
    metafile = tmp_path / "small.bgen.metafile"

    created = create_bgen_metafile(bgen_file, metafile, verbose=False)

    assert created["variants"] == 2
    assert created["metafile"] == str(metafile.resolve())
    assert created["bytes"] == metafile.stat().st_size > 0
    reread = bgen_reader.read_bgen(str(bgen_file), metafile_filepath=str(metafile), verbose=False)
    assert len(reread["variants"]) == 2

    with pytest.raises(ConfigurationError, match="Refusing to overwrite"):
        create_bgen_metafile(bgen_file, metafile, verbose=False)
    replaced = create_bgen_metafile(bgen_file, metafile, overwrite=True, verbose=False)
    assert replaced["variants"] == 2
    with pytest.raises(ConfigurationError, match="differ from the BGEN file"):
        create_bgen_metafile(bgen_file, bgen_file, overwrite=True, verbose=False)


def test_bgen_wrappers_report_missing_bgen_reader(
    bgen_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "bgen_reader", None)
    with pytest.raises(BackendUnavailableError, match="bio-bgen-reader"):
        read_bgen_dosage(bgen_file, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-bgen-reader"):
        read_bgen_allele_expectation(bgen_file, verbose=False)
