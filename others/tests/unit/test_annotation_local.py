from __future__ import annotations

import pytest

from dnakit.annotation import (
    convert_annotation_coordinates,
    make_intergenic_regions,
    make_promoter_regions,
    merge_gene_models,
)
from dnakit.exceptions import ConfigurationError


def test_convert_annotation_coordinates() -> None:
    result = convert_annotation_coordinates(
        {"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
        options={"coordinate_system": "1-based-closed", "target_coordinate_system": "0-based-half-open"},
        verbose=False,
    )
    assert result["result"]["records"][0]["start"] == 1
    assert result["result"]["records"][0]["end"] == 5


def test_make_intergenic_regions_respects_one_based_gaps() -> None:
    result = make_intergenic_regions(
        {
            "genes": [
                {"chrom": "chr1", "start": 2, "end": 4},
                {"chrom": "chr1", "start": 7, "end": 8},
            ],
            "contig_lengths": {"chr1": 10},
        },
        options={"coordinate_system": "1-based-closed"},
        verbose=False,
    )
    assert [(row["start"], row["end"]) for row in result["result"]["regions"]] == [(1, 1), (5, 6), (9, 10)]


def test_make_promoter_regions_is_strand_aware() -> None:
    result = make_promoter_regions(
        {
            "genes": [
                {"id": "plus", "chrom": "chr1", "start": 4, "end": 7, "strand": "+"},
                {"id": "minus", "chrom": "chr1", "start": 4, "end": 7, "strand": "-"},
            ],
            "contig_lengths": {"chr1": 12},
        },
        options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1},
        verbose=False,
    )
    promoters = result["result"]["promoters"]
    assert (promoters[0]["start"], promoters[0]["end"]) == (2, 4)
    assert (promoters[1]["start"], promoters[1]["end"]) == (7, 9)


def test_merge_gene_models_and_validation() -> None:
    result = merge_gene_models(
        {
            "segments": [
                {"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3},
                {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5},
                {"gene_id": "g2", "chrom": "chr1", "start": 10, "end": 11},
            ]
        },
        options={"coordinate_system": "1-based-closed", "merge_adjacent": True},
        verbose=False,
    )
    assert result["result"]["models"][0]["segments"] == [{"start": 1, "end": 5}]
    with pytest.raises(ConfigurationError):
        convert_annotation_coordinates({"records": [{"chrom": "chr1", "start": 0, "end": 2}]})
