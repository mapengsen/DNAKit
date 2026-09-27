"""Check bedMethyl units, denominators, grouping, and bounded input failures."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.genomics.methylation import read_bedmethyl, summarize_methylation

# Same mixed whitespace and motif-name format as the modkit documentation:
# https://nanoporetech.github.io/modkit/intro_pileup.html
_DOCUMENTED_ROW = "oligo_741_adapters 39 40 m,CG,0 4\t-\t39\t40\t255,0,0\t4 100.00 4 0 0 0 0 0 0\n"


def _row(
    *,
    chrom: str = "chr1",
    start: int = 10,
    code: str = "m",
    strand: str = "+",
    modified: int = 1,
    canonical: int = 3,
    other: int = 0,
) -> str:
    coverage = modified + canonical + other
    percent = 100 * modified / coverage if coverage else 0
    return (
        "\t".join(
            map(
                str,
                (
                    chrom,
                    start,
                    start + 1,
                    code,
                    coverage,
                    strand,
                    start,
                    start + 1,
                    "255,0,0",
                    coverage,
                    f"{percent:.2f}",
                    modified,
                    canonical,
                    other,
                    2,
                    3,
                    4,
                    5,
                ),
            )
        )
        + "\n"
    )


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "input.bed"
    path.write_text(text)
    return path


def test_documented_format_and_gzip_magic_detection(tmp_path: Path) -> None:
    text = "# example\ntrack name=methylation\nbrowser position chr1:1-2\n\n" + _DOCUMENTED_ROW
    path = _write(tmp_path, text)
    plain = read_bedmethyl(path, verbose=False)
    compressed = tmp_path / "compressed_without_suffix"
    compressed.write_bytes(gzip.compress(text.encode()))
    assert read_bedmethyl(compressed, verbose=False) == plain
    row = plain[0]
    assert (row["start"], row["end"], row["modification_code"]) == (39, 40, "m")
    assert (row["motif"], row["motif_offset"]) == ("CG", 0)
    assert (row["percent_modified"], row["fraction_modified"]) == (100, 1)
    json.dumps(plain, allow_nan=False)


def test_weighted_denominator_includes_other_modifications(tmp_path: Path) -> None:
    # 1/4 and 9/12 => (1+9)/(4+12)=0.625; unweighted mean is 0.5.
    rows = read_bedmethyl(
        _write(tmp_path, _row() + _row(start=20, modified=9, canonical=1, other=2)), verbose=False
    )
    result = summarize_methylation(rows, verbose=False)
    group = result["groups"][0]
    assert group["fraction_modified"] == 0.625
    assert group["valid_coverage"] == 16
    assert group["other_modified_count"] == 2
    assert group["deletion_count"] == 4  # excluded from valid coverage
    assert group["record_count"] == 2
    assert result["coordinate_system"] == "0-based-half-open"
    json.dumps(result, allow_nan=False)


def test_modification_strand_chromosome_and_motif_are_not_mixed(tmp_path: Path) -> None:
    text = (
        _row()
        + _row(code="h")
        + _row(strand="-")
        + _row(strand=".")
        + _row(chrom="chr2")
        + _row(code="m,CG,0")
        + _row(code="m,CGCG,2")
        + _row(code="21839")
    )
    rows = read_bedmethyl(_write(tmp_path, text), verbose=False)
    assert len(summarize_methylation(rows, verbose=False)["groups"]) == 8
    selected = summarize_methylation(
        rows,
        modification_codes=["m"],
        chromosomes=["chr1"],
        strands=["+"],
        min_coverage=4,
        verbose=False,
    )
    assert selected["records_selected"] == 3
    assert selected["records_filtered"] == 5
    assert summarize_methylation(rows, modification_codes=[], verbose=False)["groups"] == []


def test_zero_coverage_is_undefined_and_empty_input_is_valid(tmp_path: Path) -> None:
    rows = read_bedmethyl(_write(tmp_path, _row(modified=0, canonical=0)), verbose=False)
    assert rows[0]["fraction_modified"] is None
    group = summarize_methylation(rows, verbose=False)["groups"][0]
    assert group["fraction_modified"] is None
    assert group["percent_modified"] is None
    assert group["zero_coverage_records"] == 1
    assert summarize_methylation(rows, min_coverage=1, verbose=False)["groups"] == []
    empty = summarize_methylation([], verbose=False)
    assert empty["records_total"] == 0 and empty["groups"] == []
    json.dumps(group, allow_nan=False)


@pytest.mark.parametrize(
    "index,value",
    [
        (1, "-1"),
        (2, "10"),
        (4, "5"),
        (5, "?"),
        (7, "12"),
        (8, "256,0,0"),
        (9, "5"),
        (10, "NaN"),
        (10, "0.25"),
        (10, "101"),
        (11, "1.0"),
        (12, "-1"),
        (3, "m,CG,2"),
        (3, "5mC"),
        (3, "m,,0"),
        (3, "m,CG"),
        (3, "m,UU,0"),
    ],
)
def test_bad_fields_are_rejected_with_line_number(tmp_path: Path, index: int, value: str) -> None:
    fields = _row().split()
    fields[index] = value
    with pytest.raises(InputFormatError) as error:
        read_bedmethyl(_write(tmp_path, "# header\n" + "\t".join(fields)), verbose=False)
    assert error.value.context["line"] == 2


def test_invalid_rows_are_not_hidden_by_filters(tmp_path: Path) -> None:
    rows = read_bedmethyl(_write(tmp_path, _row()), verbose=False)
    rows[0]["modified_count"] = True
    with pytest.raises(InputFormatError):
        summarize_methylation(rows, modification_codes=[], verbose=False)


def test_record_line_and_decompressed_byte_limits_reject_partial_results(tmp_path: Path) -> None:
    text = _row() + _row(start=20)
    path = _write(tmp_path, text)
    with pytest.raises(InputFormatError, match="max_records"):
        read_bedmethyl(path, max_records=1, verbose=False)
    with pytest.raises(InputFormatError, match="limit"):
        read_bedmethyl(path, max_line_bytes=10, verbose=False)
    compressed = tmp_path / "input.gz"
    compressed.write_bytes(gzip.compress(text.encode()))
    with pytest.raises(InputFormatError, match="limit"):
        read_bedmethyl(compressed, max_bytes=len(text.encode()) - 1, verbose=False)
    rows = read_bedmethyl(path, max_bytes=len(text.encode()), verbose=False)
    with pytest.raises(InputFormatError, match="max_records"):
        summarize_methylation(iter(rows), max_records=1, verbose=False)


def test_truncated_gzip_columns_and_bedrmod_rejected(tmp_path: Path) -> None:
    path = tmp_path / "truncated.gz"
    path.write_bytes(gzip.compress(_row().encode())[:-6])
    with pytest.raises(InputFormatError):
        read_bedmethyl(path, verbose=False)
    for text in ("#fileformat=bedRMod\n" + _row(), "chr1 1 2", _row() + "extra"):
        with pytest.raises(InputFormatError):
            read_bedmethyl(_write(tmp_path, text), verbose=False)


def test_invalid_configuration_and_quiet_progress(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    path = _write(tmp_path, _row())
    for value in (0, -1, True, 1.5):
        with pytest.raises(ConfigurationError):
            read_bedmethyl(path, max_records=value, verbose=False)
    rows = read_bedmethyl(path, verbose=False)
    for kwargs in ({"strands": ["?"]}, {"modification_codes": "m"}, {"min_coverage": True}):
        with pytest.raises(ConfigurationError):
            summarize_methylation(rows, **kwargs, verbose=False)
    summarize_methylation(rows, verbose=False)
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("records", [None, 42, "", "not records", b"", {}])
def test_invalid_json_record_container_cannot_return_empty_success(records: object) -> None:
    with pytest.raises(ConfigurationError, match="iterable"):
        summarize_methylation(records, verbose=False)  # type: ignore[arg-type]


def test_json_round_trip_recomputes_fraction_from_counts(tmp_path: Path) -> None:
    rows = json.loads(json.dumps(read_bedmethyl(_write(tmp_path, _row()), verbose=False)))
    rows[0]["fraction_modified"] = None
    # The counts, not cached fraction metadata, are authoritative.
    result = summarize_methylation(rows, verbose=False)
    assert result["groups"][0]["fraction_modified"] == 0.25
    json.dumps(result, allow_nan=False)


def test_oversized_json_percent_raises_input_error_not_overflow(tmp_path: Path) -> None:
    rows = read_bedmethyl(_write(tmp_path, _row()), verbose=False)
    rows[0]["percent_modified"] = 10**400
    with pytest.raises(InputFormatError, match="percent_modified"):
        summarize_methylation(rows, verbose=False)


@pytest.mark.parametrize("path", [None, True, 42, ["input.bed"]])
def test_invalid_json_path_has_configuration_error(path: object) -> None:
    with pytest.raises(ConfigurationError, match="path"):
        read_bedmethyl(path, verbose=False)  # type: ignore[arg-type]
