from dnakit.assembly import (
    calculate_l50,
    calculate_n50,
    calculate_ng50,
    calculate_nx,
    calculate_stats,
)


def test_n50_and_l50_use_descending_contigs():
    lengths = [100, 60, 40]
    assert calculate_n50(lengths)["result"]["n50"] == 100
    assert calculate_l50(lengths)["result"]["l50"] == 1


def test_nx_and_ng50_use_explicit_denominators():
    lengths = {"lengths": [100, 60, 40]}
    assert calculate_nx(lengths, options={"x": 90})["result"]["nx"] == 40
    assert calculate_ng50(lengths, options={"genome_size": 250})["result"]["ng50"] == 60


def test_stats_contains_common_assembly_metrics():
    result = calculate_stats([100, 60, 40], options={"genome_size": 250})["result"]
    assert result["contig_count"] == 3
    assert result["total_length"] == 200
    assert result["n50"] == 100
    assert result["ng50"] == 60
