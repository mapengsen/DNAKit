"""Upstream-checked AlphaGenome variant-score wrappers; selene and enformer review (G18q)."""

from __future__ import annotations

import inspect
import json
import sys
from typing import Any

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import alphagenome_review_functions as ag
from dnakit.predictions import enformer_pytorch_review_functions as enformer_api
from dnakit.predictions import reverse_complement_one_hot
from dnakit.predictions import selene_sdk_review_functions as selene_api

ALPHAGENOME_WRAPPERS = [
    ag.convert_alphagenome_variant_scores,
    ag.list_alphagenome_recommended_scorers,
    ag.merge_alphagenome_stranded_gene_tracks,
    ag.merge_alphagenome_track_metadata,
    ag.tidy_alphagenome_anndata,
    ag.tidy_alphagenome_scores,
    ag.unmerge_alphagenome_stranded_gene_tracks,
]
TRACKS = [{"name": "T1", "strand": "."}, {"name": "T2", "strand": "."}]


def _hide_package(monkeypatch: pytest.MonkeyPatch, package: str) -> None:
    """Make ``package`` and every loaded submodule unimportable for the duration of a test."""

    for name in [key for key in sys.modules if key == package or key.startswith(package + ".")]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, package, None)


def _score_bytes(scorer: str, values: list[float], shape: list[int]) -> bytes:
    """A serialized DenseVariantScores message for one variant and one scorer."""

    import numpy as np
    from alphagenome.protos import atlas_service_pb2 as atlas
    from alphagenome.protos import dna_model_pb2 as model

    message = atlas.DenseVariantScores()
    message.variant.chromosome = "chr22"
    message.variant.position = 36201698
    message.variant.reference_bases = "A"
    message.variant.alternate_bases = "C"
    message.interval.chromosome = "chr22"
    message.interval.start = 36100000
    message.interval.end = 36300000
    message.interval.strand = model.STRAND_POSITIVE
    score = message.scores.add()
    score.variant_scorer.name = scorer
    score.shape.extend(shape)
    score.scores = np.asarray(values, dtype=np.float32).tobytes()
    return message.SerializeToString()


def _variant_anndata(scorer: str, values: list[float]) -> Any:
    """One-variant score AnnData with empty obs, as produced by the variant scoring flow."""

    import anndata
    import numpy as np
    import pandas as pd
    from alphagenome.data import genome
    from alphagenome.models.v1 import variant_scorers

    var = pd.DataFrame(
        {"name": ["T1", "T2"], "strand": [".", "."]},
        index=["0", "1"],
    )
    adata = anndata.AnnData(X=np.array([values], dtype=np.float32), var=var)
    adata.uns["interval"] = genome.Interval("chr22", 36100000, 36300000, "+")
    adata.uns["variant"] = genome.Variant("chr22", 36201698, "A", "C")
    adata.uns["variant_scorer"] = variant_scorers.RECOMMENDED_VARIANT_SCORERS[scorer]
    return adata


def _gene_anndata(values: list[list[float]], var_names: list[str], strands: list[str]) -> Any:
    import anndata
    import numpy as np
    import pandas as pd

    obs = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "strand": ["+", "-"]}, index=["0", "1"])
    var = pd.DataFrame(
        {"name": var_names, "strand": strands},
        index=[str(index) for index in range(len(var_names))],
    )
    return anndata.AnnData(X=np.array(values, dtype=float), obs=obs, var=var)


def test_every_alphagenome_wrapper_has_native_keyword_only_verbose() -> None:
    for function in ALPHAGENOME_WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_selene_and_enformer_modules_export_no_functions() -> None:
    assert list(selene_api.__all__) == []
    assert list(enformer_api.__all__) == []


def test_selene_get_reverse_encoding_raises_on_one_hot_and_rc_wrapper_matches_hand_values() -> None:
    pytest.importorskip("selene_sdk")
    import numpy as np
    from selene_sdk.sequences import sequence as upstream

    encoding = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0.25, 0.25, 0.25, 0.25]],
        dtype=np.float32,
    )
    bases = ["A", "C", "G", "T"]
    base_to_index = {base: index for index, base in enumerate(bases)}
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    with pytest.raises(IndexError):
        upstream.get_reverse_encoding(encoding, bases, base_to_index, complement)

    # ACGN reverse-complements to NCGT; the N row stays uniform.
    rc = reverse_complement_one_hot(encoding.tolist(), verbose=False)["encoding"]
    assert rc == [[0.25] * 4, [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


def test_convert_variant_scores_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("alphagenome")
    import pandas as pd
    from alphagenome.atlas import atlas as upstream
    from alphagenome.protos import atlas_service_pb2

    payload = _score_bytes("ATAC", [0.5, -1.0], [1, 2])
    result = ag.convert_alphagenome_variant_scores(
        [payload], track_metadata={"ATAC": TRACKS}, verbose=False
    )
    assert result["operation"] == "convert_alphagenome_variant_scores"
    assert result["variant_count"] == 1
    assert [item["scorer"] for item in result["scorers"]] == ["ATAC"]
    record = result["scorers"][0]["anndata"]
    assert record["shape"] == [1, 2]
    assert record["scores"] == [[0.5, -1.0]]
    assert record["var"] == TRACKS
    assert record["obs"] == [
        {"variant": "chr22:36201698:A>C", "interval": "chr22:36100000-36300000:+"}
    ]

    direct = upstream.convert_variant_scores_to_anndata(
        [atlas_service_pb2.DenseVariantScores.FromString(payload)],
        {"ATAC": pd.DataFrame(TRACKS)},
    )["ATAC"]
    assert direct.X.tolist() == record["scores"]
    json.dumps(result)


def test_convert_variant_scores_input_validation() -> None:
    pytest.importorskip("alphagenome")
    payload = _score_bytes("ATAC", [0.5, -1.0], [1, 2])
    with pytest.raises(ConfigurationError):
        ag.convert_alphagenome_variant_scores([], track_metadata={"ATAC": TRACKS}, verbose=False)
    with pytest.raises(ConfigurationError):
        ag.convert_alphagenome_variant_scores(["not bytes"], track_metadata={"ATAC": TRACKS})
    with pytest.raises(ConfigurationError):
        ag.convert_alphagenome_variant_scores([b"\x0a\x03abc"], track_metadata={"ATAC": TRACKS})
    with pytest.raises(ConfigurationError, match="no tracks for scorers"):
        ag.convert_alphagenome_variant_scores([payload], track_metadata={"CAGE": TRACKS})
    with pytest.raises(ConfigurationError):
        ag.convert_alphagenome_variant_scores([payload], track_metadata={})


def test_tidy_anndata_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.models.v1 import variant_scorers

    adata = _variant_anndata("ATAC", [0.5, -1.0])
    result = ag.tidy_alphagenome_anndata(adata, verbose=False)
    records = result["records"]
    assert result["operation"] == "tidy_alphagenome_anndata"
    assert result["row_count"] == 2
    assert [row["track_name"] for row in records] == ["T1", "T2"]
    assert [row["raw_score"] for row in records] == [0.5, -1.0]
    assert {row["output_type"] for row in records} == {"ATAC"}
    assert records[0]["variant_id"] == "chr22:36201698:A>C"
    assert records[0]["scored_interval"] == "chr22:36100000-36300000:+"

    direct = variant_scorers.tidy_anndata(adata)
    assert direct["raw_score"].tolist() == [row["raw_score"] for row in records]
    json.dumps(result)

    flags = ag.tidy_alphagenome_anndata(
        adata, match_gene_strand=False, include_extended_metadata=False, verbose=False
    )
    assert flags["row_count"] == 2


def test_tidy_scores_treats_flat_and_nested_inputs_as_the_same_groups() -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.models.v1 import variant_scorers

    first = _variant_anndata("ATAC", [0.5, -1.0])
    second = _variant_anndata("DNASE", [2.0, 3.0])
    nested = ag.tidy_alphagenome_scores([[first], [second]], verbose=False)
    flat = ag.tidy_alphagenome_scores([first, second], verbose=False)
    assert nested["records"] == flat["records"]
    assert flat["group_count"] == 2
    assert [
        (row["output_type"], row["track_name"], row["raw_score"]) for row in flat["records"]
    ] == [
        ("ATAC", "T1", 0.5),
        ("ATAC", "T2", -1.0),
        ("DNASE", "T1", 2.0),
        ("DNASE", "T2", 3.0),
    ]
    direct = variant_scorers.tidy_scores([[first], [second]])
    assert direct["raw_score"].tolist() == [0.5, -1.0, 2.0, 3.0]
    json.dumps(flat)


def test_tidy_inputs_are_validated() -> None:
    pytest.importorskip("alphagenome")
    adata = _variant_anndata("ATAC", [0.5, -1.0])
    with pytest.raises(ConfigurationError):
        ag.tidy_alphagenome_anndata([adata], verbose=False)
    with pytest.raises(ConfigurationError):
        ag.tidy_alphagenome_anndata(adata, match_gene_strand="yes", verbose=False)
    with pytest.raises(ConfigurationError):
        ag.tidy_alphagenome_scores([], verbose=False)
    with pytest.raises(ConfigurationError):
        ag.tidy_alphagenome_scores([adata, [1]], verbose=False)


def test_merge_track_metadata_hand_values_and_duplicate_removal() -> None:
    pytest.importorskip("alphagenome")
    rows = [
        {"name": "A", "strand": "+"},
        {"name": "A", "strand": "-"},
        {"name": "B", "strand": "+"},
        {"name": "B", "strand": "-"},
        {"name": "C", "strand": "."},
    ]
    result = ag.merge_alphagenome_track_metadata(rows, verbose=False)
    assert result["operation"] == "merge_alphagenome_track_metadata"
    assert result["track_count"] == 3
    assert [(track["name"], track["strand"]) for track in result["tracks"]] == [
        ("A", "."),
        ("B", "."),
        ("C", "."),
    ]

    # An unstranded track with the same name as a stranded pair is dropped.
    duplicate = [
        {"name": "A", "strand": "+"},
        {"name": "A", "strand": "-"},
        {"name": "A", "strand": "."},
    ]
    merged = ag.merge_alphagenome_track_metadata(duplicate, verbose=False)
    assert merged["tracks"] == [{"name": "A", "strand": "."}]
    json.dumps(result)


def test_merge_track_metadata_validation() -> None:
    pytest.importorskip("alphagenome")
    with pytest.raises(ConfigurationError):
        ag.merge_alphagenome_track_metadata([{"name": "A", "strand": "x"}], verbose=False)
    with pytest.raises(ConfigurationError):
        ag.merge_alphagenome_track_metadata([{"name": "A"}], verbose=False)
    with pytest.raises(ConfigurationError):
        # Positive track A has no negative partner, so the strand lists differ.
        ag.merge_alphagenome_track_metadata(
            [{"name": "A", "strand": "+"}, {"name": "B", "strand": "-"}], verbose=False
        )


def test_merge_stranded_gene_tracks_matches_upstream_and_hand_values() -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.models.v1 import variant_scoring_utils

    scores = _gene_anndata(
        [[1, 2, 3, 4], [5, 6, 7, 8]],
        var_names=["A", "A", "B", "B"],
        strands=["+", "-", "+", "-"],
    )
    result = ag.merge_alphagenome_stranded_gene_tracks(scores, verbose=False)
    record = result["anndata"]
    # Gene 1 is on +, so it takes the positive columns (A+, B+); gene 2 is on -, so the
    # negative columns (A-, B-).
    assert record["shape"] == [2, 2]
    assert record["scores"] == [[1.0, 3.0], [6.0, 8.0]]
    assert record["var"] == [{"name": "A", "strand": "."}, {"name": "B", "strand": "."}]
    assert record["obs"] == [
        {"gene_id": "ENSG1", "strand": "+"},
        {"gene_id": "ENSG2", "strand": "-"},
    ]
    direct = variant_scoring_utils.merge_stranded_gene_tracks(scores)
    assert direct.X.tolist() == record["scores"]
    json.dumps(result)


def test_unmerge_stranded_gene_tracks_matches_hand_values_with_nulls() -> None:
    pytest.importorskip("alphagenome")
    merged = _gene_anndata([[1.0, 3.0], [6.0, 8.0]], var_names=["A", "B"], strands=[".", "."])
    tracks = [
        {"name": "A", "strand": "+"},
        {"name": "A", "strand": "-"},
        {"name": "B", "strand": "+"},
        {"name": "B", "strand": "-"},
    ]
    result = ag.unmerge_alphagenome_stranded_gene_tracks(
        merged, track_metadata=tracks, verbose=False
    )
    record = result["anndata"]
    # Gene 1 (+) keeps A+ and B+; gene 2 (-) keeps A- and B-; the other columns are null.
    assert record["scores"] == [[1.0, None, 3.0, None], [None, 6.0, None, 8.0]]
    assert record["var"] == tracks
    assert record["var_names"] == ["0", "1", "2", "3"]
    json.dumps(result)

    with pytest.raises(ConfigurationError):
        ag.unmerge_alphagenome_stranded_gene_tracks(
            merged, track_metadata=[{"name": "A"}], verbose=False
        )


def test_list_recommended_scorers_matches_upstream_registry() -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.models.v1 import dna_model, variant_scorers

    for organism, species in (
        ("homo_sapiens", dna_model.Organism.HOMO_SAPIENS),
        ("mus_musculus", dna_model.Organism.MUS_MUSCULUS),
    ):
        result = ag.list_alphagenome_recommended_scorers(organism, verbose=False)
        recommended = variant_scorers.get_recommended_scorers(species.value)
        expected = sorted(
            key
            for key, value in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()
            if any(value is item for item in recommended)
        )
        assert [item["name"] for item in result["scorers"]] == expected
        assert result["scorer_count"] == len(expected)
        assert result["organism"] == organism
    assert ag.list_alphagenome_recommended_scorers(verbose=False)["scorer_count"] == 19
    assert (
        ag.list_alphagenome_recommended_scorers("mus_musculus", verbose=False)["scorer_count"] == 18
    )
    with pytest.raises(ConfigurationError):
        ag.list_alphagenome_recommended_scorers("danio_rerio", verbose=False)


def test_alphagenome_wrappers_raise_backend_unavailable_without_alphagenome(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _hide_package(monkeypatch, "alphagenome")
    with pytest.raises(BackendUnavailableError):
        ag.list_alphagenome_recommended_scorers(verbose=False)
    with pytest.raises(BackendUnavailableError):
        ag.merge_alphagenome_track_metadata(TRACKS, verbose=False)
