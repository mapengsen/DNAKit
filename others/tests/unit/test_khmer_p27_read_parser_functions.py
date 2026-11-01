"""Compare the public read/pair output with the actual khmer parser."""

import json

import pytest

khmer = pytest.importorskip("khmer")

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import khmer_p27_read_parser_functions as m


@pytest.mark.parametrize("text", [">one\nACGT\n", ">one\nacgtn\n>two\nN\n", "@one\nACGT\n+\nIIII\n"])
def test_read_records_match_direct_parser(tmp_path, text):
    path = tmp_path / "reads.fa"
    path.write_text(text)
    expected = list(khmer.ReadParser(str(path)))
    actual = m.khmer_parse_reads(path, verbose=False)
    assert [r["sequence"] for r in actual["records"]] == [r.sequence for r in expected]
    assert [r["name"] for r in actual["records"]] == [r.name for r in expected]
    assert actual["truncated"] is False
    json.dumps(actual, allow_nan=False)


@pytest.mark.parametrize("length", [1, 4, 100])
def test_pairs_match_direct_parser(tmp_path, length):
    path = tmp_path / "pairs.fa"
    path.write_text(f">r/1\n{'A' * length}\n>r/2\n{'T' * length}\n")
    expected = list(khmer.ReadParser(str(path)).iter_read_pairs())
    actual = m.khmer_parse_read_pairs(path, verbose=False)
    assert [[r["name"] for r in pair] for pair in actual["records"]] == [[r.name for r in pair] for pair in expected]
    assert actual["base_count"] == 2 * length


def test_limits_report_truncation_and_invalid_pairs(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text(">a\nACGT\n>b\nACGT\n")
    assert m.khmer_parse_reads(path, max_records=1, verbose=False)["truncated"]
    actual = m.khmer_parse_reads(path, max_bases=3, verbose=False)
    assert actual["records"] == [] and actual["truncated"]
    with pytest.raises(ConfigurationError):
        m.khmer_parse_read_pairs(path, verbose=False)
    with pytest.raises(ConfigurationError):
        m.khmer_parse_reads(path, max_records=True, verbose=False)
