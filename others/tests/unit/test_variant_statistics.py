"""Contract tests for dependency-free VCF-like variant statistics."""

from dnakit.variants import (
    calculate_allele_frequency,
    calculate_heterozygosity,
    calculate_missingness,
    calculate_transition_transversion,
    calculate_variant_density,
)


def _variants() -> dict:
    return {
        "variants": [
            {
                "chrom": "chr1",
                "pos": 10,
                "ref": "A",
                "alt": "G",
                "genotypes": ["0/1", "1/1", "0/0"],
            }
        ]
    }


def test_variant_statistics_are_reproducible() -> None:
    assert calculate_allele_frequency(_variants(), verbose=False)["result"]["aggregate_frequency"] == 0.5
    assert calculate_heterozygosity(_variants(), verbose=False)["result"]["observed_heterozygosity"] == 1 / 3
    assert calculate_missingness(_variants(), verbose=False)["result"]["missingness"] == 0.0
    transition = calculate_transition_transversion(_variants(), verbose=False)["result"]
    assert transition["transitions"] == 1
    assert transition["transversions"] == 0
    density = calculate_variant_density(_variants(), options={"region_length": 100}, verbose=False)
    assert density["result"]["density_per_base"] == 0.01
