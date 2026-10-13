"""Family-variant verification for group H5: new wrappers and the checks behind the ledger rows."""

from __future__ import annotations

import inspect
import random
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError


# ---------------------------------------------------------------- sourmash
def _sm_seq(rng: random.Random, size: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(size))


def _sm_minhash(sequence: str, *, track: bool = False, **kwargs: Any) -> Any:
    import sourmash

    minhash = sourmash.MinHash(track_abundance=track, **kwargs)
    minhash.add_sequence(sequence, True)
    minhash.add_sequence(sequence[:700], True)
    return minhash


def _sm_sketch(sequence: str, name: str, *, track: bool = False) -> dict[str, Any]:
    minhash = _sm_minhash(sequence, track=track, n=0, ksize=21, scaled=5)
    sketch: dict[str, Any] = {
        "id": name,
        "ksize": 21,
        "scaled": 5,
        "hashes": sorted(minhash.hashes),
    }
    if track:
        sketch["abundances"] = [minhash.hashes[key] for key in sketch["hashes"]]
    return sketch


@pytest.fixture(scope="module")
def sm_data() -> dict[str, Any]:
    pytest.importorskip("sourmash")
    rng = random.Random(1)
    first = _sm_seq(rng, 3000)
    second = first[:1800] + _sm_seq(rng, 1200)
    third = _sm_seq(rng, 3000)
    named = [(first, "a"), (second, "b"), (third, "c")]
    return {
        "sequences": named,
        "plain": [_sm_sketch(seq, name) for seq, name in named],
        "tracked": [_sm_sketch(seq, name, track=True) for seq, name in named],
    }


def _sm_signatures(sketches: list[dict[str, Any]]) -> list[Any]:
    import sourmash

    out = []
    for sketch in sketches:
        minhash = sourmash.MinHash(n=0, ksize=sketch["ksize"], scaled=sketch["scaled"])
        minhash.add_many(sketch["hashes"])
        out.append(sourmash.SourmashSignature(minhash, name=sketch["id"]))
    return out


def test_sm_archive_round_trip(tmp_path: Path, sm_data: dict[str, Any]) -> None:
    from dnakit.microbiome import read_sourmash_signatures
    from dnakit.microbiome.sourmash_family_functions import write_sourmash_signature_archive

    for name, expected in (("x.sig.zip", "zip"), ("y.zip", "zip"), ("z.sqldb", "sqlite")):
        result = write_sourmash_signature_archive(sm_data["plain"], tmp_path / name, verbose=False)
        assert result["format"] == expected and result["signatures"] == 3
        back = read_sourmash_signatures(tmp_path / name, verbose=False)
        assert [item["id"] for item in back["sketches"]] == ["a", "b", "c"]
        assert [item["hashes"] for item in back["sketches"]] == [
            item["hashes"] for item in sm_data["plain"]
        ]
    with pytest.raises(ConfigurationError):
        write_sourmash_signature_archive(sm_data["plain"], tmp_path / "x.sig.zip", verbose=False)
    with pytest.raises(ConfigurationError):
        write_sourmash_signature_archive(sm_data["plain"], tmp_path / "x.sig", verbose=False)


def test_sm_archive_keeps_abundance(tmp_path: Path, sm_data: dict[str, Any]) -> None:
    import sourmash

    from dnakit.microbiome.sourmash_family_functions import write_sourmash_signature_archive

    path = tmp_path / "t.zip"
    write_sourmash_signature_archive(sm_data["tracked"], path, verbose=False)
    loaded = list(sourmash.load_file_as_signatures(str(path)))
    expected = sm_data["tracked"][0]
    assert loaded[0].minhash.track_abundance
    assert [loaded[0].minhash.hashes[key] for key in expected["hashes"]] == expected["abundances"]


def _sm_json(sm_data: dict[str, Any], tmp_path: Path) -> str:
    import sourmash

    first, second, _third = (seq for seq, _name in sm_data["sequences"])
    signatures = [
        sourmash.SourmashSignature(_sm_minhash(first, n=0, ksize=21, scaled=5), name="a"),
        sourmash.SourmashSignature(_sm_minhash(second, n=0, ksize=21, scaled=5), name="b"),
        sourmash.SourmashSignature(_sm_minhash(first, n=0, ksize=31, scaled=5), name="k31"),
    ]
    path = tmp_path / "json.sig"
    with sourmash.save_load.SaveSignaturesToLocation(str(path)) as saver:
        for signature in signatures:
            saver.add(signature)
    return path.read_text()


def test_sm_parse_signature_json(tmp_path: Path, sm_data: dict[str, Any]) -> None:
    from sourmash.signature import load_signatures_from_json

    from dnakit.microbiome.sourmash_family_functions import parse_sourmash_signature_json

    text = _sm_json(sm_data, tmp_path)
    result = parse_sourmash_signature_json(text, verbose=False)
    direct = sorted(load_signatures_from_json(text), key=lambda item: item.name)
    assert [item["id"] for item in result["sketches"]] == [item.name for item in direct]
    assert result["sketches"][0]["hashes"] == sorted(direct[0].minhash.hashes)
    assert result["sketches"][0]["hashes"] == sm_data["plain"][0]["hashes"]
    only31 = parse_sourmash_signature_json(text, ksize=31, require_single=True, verbose=False)
    assert [item["id"] for item in only31["sketches"]] == ["k31"]
    with pytest.raises(ConfigurationError):
        parse_sourmash_signature_json(text, require_single=True, verbose=False)
    with pytest.raises(ConfigurationError):
        parse_sourmash_signature_json("not json", verbose=False)
    with pytest.raises(ConfigurationError):
        parse_sourmash_signature_json(text, max_signatures=1, verbose=False)


def test_sm_select_signatures(tmp_path: Path, sm_data: dict[str, Any]) -> None:
    from dnakit.microbiome.sourmash_family_functions import select_sourmash_signatures

    (tmp_path / "s.sig").write_text(_sm_json(sm_data, tmp_path))
    path = tmp_path / "s.sig"

    def ids(**kwargs: Any) -> list[str]:
        result = select_sourmash_signatures([path], verbose=False, **kwargs)
        return [item["id"] for item in result["sketches"]]

    assert ids() == ["a", "b", "k31"]
    assert ids(ksize=21) == ["a", "b"]
    assert ids(name_pattern="^B") == ["b"]
    assert ids(ksize=99) == []
    assert ids(expect="one", ksize=31) == ["k31"]
    assert ids(expect="query", ksize=31) == ["k31"]
    first = select_sourmash_signatures([path], ksize=21, verbose=False)["sketches"][0]
    assert ids(expect="query", ksize=21, md5_prefix=first["md5sum"][:8]) == ["a"]
    for bad in (
        {"expect": "one"},
        {"expect": "query", "name_pattern": "a"},
        {"ksize": 0},
        {"expect": "other"},
    ):
        with pytest.raises(ConfigurationError):
            select_sourmash_signatures([path], verbose=False, **bad)


def test_sm_abundance_matrix_matches_upstream(sm_data: dict[str, Any]) -> None:
    import numpy as np
    import sourmash
    from sourmash.compare import compare_serial

    from dnakit.microbiome.sourmash_family_functions import compare_minhash_abundance_matrix

    signatures = []
    for sketch in sm_data["tracked"]:
        minhash = sourmash.MinHash(n=0, ksize=21, scaled=5, track_abundance=True)
        minhash.set_abundances(dict(zip(sketch["hashes"], sketch["abundances"], strict=True)))
        signatures.append(sourmash.SourmashSignature(minhash, name=sketch["id"]))
    direct = compare_serial(signatures, False, downsample=False)
    result = compare_minhash_abundance_matrix(sm_data["tracked"], verbose=False)
    assert result["measure"] == "angular_similarity" and result["abundance_weighted"]
    assert np.allclose(result["matrix"], direct)
    parallel = compare_minhash_abundance_matrix(sm_data["tracked"], n_jobs=2, verbose=False)
    assert parallel["matrix"] == result["matrix"]
    jaccard = compare_minhash_abundance_matrix(
        sm_data["tracked"], ignore_abundance=True, verbose=False
    )
    assert jaccard["measure"] == "jaccard"
    assert jaccard["matrix"][0][1] != result["matrix"][0][1]
    with pytest.raises(ConfigurationError):
        compare_minhash_abundance_matrix(
            sm_data["tracked"][:2] + sm_data["plain"][2:], verbose=False
        )


def test_sm_abundance_matrix_plain_equals_jaccard_wrapper(sm_data: dict[str, Any]) -> None:
    from dnakit.microbiome import compare_minhash_matrix
    from dnakit.microbiome.sourmash_family_functions import compare_minhash_abundance_matrix

    plain = compare_minhash_abundance_matrix(sm_data["plain"], verbose=False)
    assert plain["measure"] == "jaccard"
    assert plain["matrix"] == compare_minhash_matrix(sm_data["plain"], verbose=False)["matrix"]
    ani = compare_minhash_abundance_matrix(sm_data["plain"], return_ani=True, verbose=False)
    assert ani["measure"] == "ani" and ani["matrix"][0][0] == 1.0


def test_sm_search_ranked_by_hand(sm_data: dict[str, Any]) -> None:
    from dnakit.microbiome.sourmash_family_functions import search_minhash_database_ranked

    plain = sm_data["plain"]
    query = set(plain[0]["hashes"])
    other = set(plain[1]["hashes"])
    jaccard = len(query & other) / len(query | other)
    contain = len(query & other) / len(query)
    maximum = len(query & other) / min(len(query), len(other))
    for mode, expected in (
        ("jaccard", jaccard),
        ("containment", contain),
        ("max_containment", maximum),
    ):
        result = search_minhash_database_ranked(plain[0], plain, mode=mode, verbose=False)
        assert [m["id"] for m in result["matches"]] == ["a", "b"]
        assert result["matches"][0]["score"] == pytest.approx(1.0)
        assert result["matches"][1]["score"] == pytest.approx(expected)


def test_sm_search_ranked_abundance_matches_upstream(sm_data: dict[str, Any]) -> None:
    import numpy as np
    import sourmash

    from dnakit.microbiome.sourmash_family_functions import (
        compare_minhash_abundance_matrix,
        search_minhash_database_ranked,
    )

    tracked = sm_data["tracked"]
    result = search_minhash_database_ranked(tracked[0], tracked, mode="abundance", verbose=False)
    matrix = compare_minhash_abundance_matrix(tracked, verbose=False)["matrix"]
    assert [m["id"] for m in result["matches"]] == ["a", "b"]
    assert np.isclose(result["matches"][1]["score"], matrix[0][1], atol=1e-3)
    del sourmash
    with pytest.raises(ConfigurationError):
        search_minhash_database_ranked(
            tracked[0], sm_data["plain"], mode="abundance", verbose=False
        )
    with pytest.raises(ConfigurationError):
        search_minhash_database_ranked(tracked[0], tracked, mode="jaccard", verbose=False)


def test_sm_sbt_search_equals_linear(tmp_path: Path, sm_data: dict[str, Any]) -> None:
    from sourmash.sbtmh import create_sbt_index

    from dnakit.microbiome.sourmash_family_functions import (
        search_minhash_database_ranked,
        search_sourmash_sbt_index,
    )

    tree = create_sbt_index()
    for signature in _sm_signatures(sm_data["plain"]):
        tree.insert(signature)
    path = tmp_path / "x.sbt.zip"
    tree.save(str(path))
    result = search_sourmash_sbt_index(path, sm_data["plain"][0], threshold=0.01, verbose=False)
    linear = search_minhash_database_ranked(
        sm_data["plain"][0], sm_data["plain"], threshold=0.01, verbose=False
    )
    assert [(m["id"], m["jaccard"]) for m in result["matches"]] == [
        (m["id"], m["score"]) for m in linear["matches"]
    ]
    with pytest.raises(ConfigurationError):
        search_sourmash_sbt_index(tmp_path / "missing.sbt.zip", sm_data["plain"][0], verbose=False)


def test_sm_downsample_num_sketch() -> None:
    pytest.importorskip("sourmash")
    from dnakit.microbiome.sourmash_family_functions import downsample_minhash_num_sketch

    rng = random.Random(5)
    minhash = _sm_minhash(_sm_seq(rng, 2000), n=50, ksize=21)
    sketch = {"id": "n", "ksize": 21, "num": 50, "hashes": sorted(minhash.hashes)}
    result = downsample_minhash_num_sketch(sketch, num=20, verbose=False)
    assert result["num"] == 20 and result["hashes"] == sketch["hashes"][:20]
    with pytest.raises(ConfigurationError):
        downsample_minhash_num_sketch(sketch, num=60, verbose=False)


def test_sm_family_compare_metric_reaches_upstream(sm_data: dict[str, Any]) -> None:
    import numpy as np
    from sourmash import compare

    from dnakit.microbiome import compare_minhash_matrix

    assert "metric" in inspect.signature(compare_minhash_matrix).parameters
    signatures = _sm_signatures(sm_data["plain"])
    for metric, function in (
        ("containment", compare.compare_serial_containment),
        ("max_containment", compare.compare_serial_max_containment),
        ("avg_containment", compare.compare_serial_avg_containment),
    ):
        wrapped = compare_minhash_matrix(sm_data["plain"], metric=metric, verbose=False)["matrix"]
        assert np.allclose(wrapped, function(signatures, downsample=False).T)


def test_sm_family_mutation_rate_containment_reaches_upstream() -> None:
    pytest.importorskip("sourmash")
    from sourmash.minhash import containment_to_distance

    from dnakit.microbiome import estimate_mutation_rate

    assert "measure" in inspect.signature(estimate_mutation_rate).parameters
    direct = containment_to_distance(
        0.6, 21, 1000, sequence_len_bp=5_000_000, confidence=0.95, estimate_ci=True
    )
    wrapped = estimate_mutation_rate(
        0.6, measure="containment", ksize=21, scaled=1000, sequence_len_bp=5_000_000, verbose=False
    )
    assert wrapped["distance"] == pytest.approx(direct.dist)
    assert wrapped["distance_low"] == pytest.approx(direct.dist_low)


def test_sm_family_read_reaches_zip_sqlite_and_nodegraph(
    tmp_path: Path, sm_data: dict[str, Any]
) -> None:
    import sourmash
    from sourmash.nodegraph import Nodegraph, extract_nodegraph_info

    from dnakit.microbiome import inspect_khmer_graph, read_sourmash_signatures

    assert "source" in inspect.signature(read_sourmash_signatures).parameters
    for name in ("p.sig.zip", "p.sqldb"):
        with sourmash.save_load.SaveSignaturesToLocation(str(tmp_path / name)) as saver:
            for signature in _sm_signatures(sm_data["plain"]):
                saver.add(signature)
        loaded = read_sourmash_signatures(tmp_path / name, verbose=False)
        assert [item["id"] for item in loaded["sketches"]] == ["a", "b", "c"]
    graph = Nodegraph(21, 1000, 4)
    for sequence, _name in sm_data["sequences"][:1]:
        for start in range(len(sequence) - 20):
            graph.count(sequence[start : start + 21])
    graph.save(str(tmp_path / "g.ng"))
    pytest.importorskip("khmer")
    info = inspect_khmer_graph(tmp_path / "g.ng", kind="nodegraph", verbose=False)
    upstream = extract_nodegraph_info(str(tmp_path / "g.ng"))
    assert (info["ksize"], info["table_size"], info["n_tables"], info["n_occupied"]) == (
        upstream[0],
        upstream[1],
        upstream[2],
        upstream[5],
    )


def test_sm_family_nodegraph_collisions_equal_khmer(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    import khmer
    from sourmash.nodegraph import Nodegraph, calc_expected_collisions

    from dnakit.microbiome import expected_kmer_table_collisions

    rng = random.Random(3)
    sequence = _sm_seq(rng, 3000)
    fasta = tmp_path / "a.fa"
    fasta.write_text(f">r1\n{sequence}\n")
    sourmash_graph = Nodegraph(21, 5000, 4)
    for start in range(len(sequence) - 20):
        sourmash_graph.count(sequence[start : start + 21])
    rate = calc_expected_collisions(sourmash_graph, force=True, max_false_pos=1.0)
    wrapped = expected_kmer_table_collisions(
        fasta, ksize=21, table_size=5000, n_tables=4, verbose=False
    )
    assert wrapped["expected_false_positive_rate"] == pytest.approx(rate)
    assert wrapped["occupied_bins"] == sourmash_graph.n_occupied()
    del khmer


def test_sm_aliases_are_same_objects() -> None:
    pytest.importorskip("sourmash")
    import sourmash
    from sourmash import (
        commands,
        compare,
        distance_utils,
        index,
        minhash,
        nodegraph,
        save_load,
        sbt,
        sbtmh,
        signature,
        sourmash_args,
    )
    from sourmash.tax import tax_utils

    assert commands.SaveSignaturesToLocation is save_load.SaveSignaturesToLocation
    assert sourmash_args.SaveSignaturesToLocation is save_load.SaveSignaturesToLocation
    assert index.load_signatures_from_json is signature.load_signatures_from_json
    assert sourmash.load_signatures_from_json is signature.load_signatures_from_json
    assert sourmash.load_one_signature_from_json is signature.load_one_signature_from_json
    assert sourmash.load_file_as_index is save_load.load_file_as_index
    assert sourmash_args.load_file_as_index is save_load.load_file_as_index
    for kind in ("containment", "max_containment", "avg_containment"):
        assert getattr(commands, f"compare_serial_{kind}") is getattr(
            compare, f"compare_serial_{kind}"
        )
    assert commands.compare_all_pairs is compare.compare_all_pairs
    assert distance_utils.containment_to_distance is minhash.containment_to_distance
    assert tax_utils.containment_to_distance is minhash.containment_to_distance
    assert index.flatten_and_downsample_num is minhash.flatten_and_downsample_num
    assert sbt.calc_expected_collisions is nodegraph.calc_expected_collisions
    assert sbt.extract_nodegraph_info is nodegraph.extract_nodegraph_info
    assert sourmash.search_sbt_index_private is sbtmh.search_sbt_index


# ------------------------------------------------------------- gimmemotifs
def test_gm_normalize_peak_widths_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from dnakit.regulation.gimmemotifs_family_functions import normalize_peak_widths

    path = tmp_path / "peaks.bed"
    path.write_text(
        "chr1\t100\t300\tp1\t5\t+\nchr1\t1000\t1100\tp2\t7\t-\nchr2\t10\t50\n\n# note\n"
        "chr3\t500\t700\tp4\n"
    )
    result = normalize_peak_widths(path, window=100, verbose=False)
    rows = [
        (i["chrom"], i["start"], i["end"], i["name"], i["score"], i["strand"])
        for i in result["intervals"]
    ]
    # midpoint 200 -> 150..250; minus strand midpoint 1050 -> start 1100, end 1000 (swapped);
    # chr2 midpoint 30 gives start -20 and is dropped; chr3 midpoint 600 -> 550..650 with defaults.
    assert rows == [
        ("chr1", 150, 250, "p1", "5", "+"),
        ("chr1", 1100, 1000, "p2", "7", "-"),
        ("chr3", 550, 650, "p4", "0", "+"),
    ]
    assert [i["start_after_end"] for i in result["intervals"]] == [False, True, False]
    assert result["dropped_rows"] == 1 and result["input_rows"] == 4
    with pytest.raises(ConfigurationError):
        normalize_peak_widths(path, window=0, verbose=False)
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\tx\t5\n")
    with pytest.raises(ConfigurationError):
        normalize_peak_widths(bad, verbose=False)


def test_gm_normalize_differs_from_equal_size(tmp_path: Path) -> None:
    pytest.importorskip("gimmemotifs")
    from dnakit.regulation import equal_size_intervals
    from dnakit.regulation.gimmemotifs_family_functions import normalize_peak_widths

    path = tmp_path / "p.bed"
    path.write_text("chr2\t10\t50\tx\t1\t+\nchr1\t1000\t1100\ty\t2\t-\n")
    equal = equal_size_intervals(path, size=100, verbose=False)["intervals"]
    strand = normalize_peak_widths(path, window=100, verbose=False)["intervals"]
    assert [i["start"] for i in equal] == [0, 1000]  # clamped, strand ignored
    assert [i["start"] for i in strand] == [1100]  # row with negative start dropped


def test_gm_family_metric_values_reach_rocmetrics() -> None:
    pytest.importorskip("gimmemotifs")
    import numpy as np
    from gimmemotifs import rocmetrics

    from dnakit.regulation import motif_score_metric

    assert {"metric", "fpr", "fdr", "minbg"} <= set(
        inspect.signature(motif_score_metric).parameters
    )
    rng = np.random.default_rng(7)
    foreground = rng.normal(1.0, 1.0, 80)
    background = rng.normal(0.0, 1.0, 300)
    plain = {
        "pr_auc": rocmetrics.pr_auc,
        "mncp": rocmetrics.mncp,
        "max_fmeasure": rocmetrics.max_fmeasure,
        "max_enrichment": rocmetrics.max_enrichment,
    }
    at_fpr = {
        "phyper_at_fpr": rocmetrics.phyper_at_fpr,
        "fraction_fpr": rocmetrics.fraction_fpr,
        "score_at_fpr": rocmetrics.score_at_fpr,
        "enr_at_fpr": rocmetrics.enr_at_fpr,
    }
    for name, function in plain.items():
        value = motif_score_metric(list(foreground), list(background), metric=name, verbose=False)[
            "value"
        ]
        assert value == pytest.approx(float(function(foreground, background)))
    for name, function in at_fpr.items():
        value = motif_score_metric(
            list(foreground), list(background), metric=name, fpr=0.05, verbose=False
        )["value"]
        assert value == pytest.approx(float(function(foreground, background, fpr=0.05)))
    recall = motif_score_metric(
        list(foreground), list(background), metric="recall_at_fdr", fdr=0.2, verbose=False
    )["value"]
    assert recall == pytest.approx(
        float(rocmetrics.recall_at_fdr(foreground, background, fdr_cutoff=0.2))
    )
    xs, ys = rocmetrics.roc_values(foreground, background)
    curve = motif_score_metric(
        list(foreground), list(background), metric="roc_values", verbose=False
    )["curve"]
    assert curve["fpr"] == [float(v) for v in xs] and curve["tpr"] == [float(v) for v in ys]


def test_gm_family_ks_significance_and_position_metrics() -> None:
    pytest.importorskip("gimmemotifs")
    import numpy as np
    from gimmemotifs import comparison, rocmetrics

    from dnakit.regulation import compare_motif_positions, motif_position_ks

    assert "output" in inspect.signature(motif_position_ks).parameters
    positions = [10.0, 12.0, 15.0, 40.0, 41.0, 90.0]
    background = [1.0, 30.0, 55.0, 60.0, 80.0, 99.0, 70.0]
    value = motif_position_ks(positions, background, output="significance", verbose=False)["value"]
    assert value == pytest.approx(
        float(rocmetrics.ks_significance(np.array(positions), np.array(background)))
    )
    pfm_a = [[8, 1, 1, 0], [0, 9, 1, 0], [1, 1, 7, 1]]
    pfm_b = [[7, 1, 1, 1], [1, 7, 1, 1], [2, 2, 5, 1]]
    for metric in ("chisq", "ssd"):
        result = compare_motif_positions(pfm_a, pfm_b, metric=metric, verbose=False)
        left = np.asarray(pfm_a, dtype=float)
        right = np.asarray(pfm_b, dtype=float)
        from gimmemotifs.utils import ppm_pseudocount

        left = ppm_pseudocount(left / left.sum(axis=1, keepdims=True))
        right = ppm_pseudocount(right / right.sum(axis=1, keepdims=True))
        direct = [
            float(getattr(comparison, metric)(a, b)) for a, b in zip(left, right, strict=True)
        ]
        assert result["position_scores"] == pytest.approx(direct)


def test_gm_count_motif_matches_at_fpr_by_hand() -> None:
    pytest.importorskip("gimmemotifs")
    from dnakit.regulation import motif_score_metric
    from dnakit.regulation.gimmemotifs_family_functions import count_motif_matches_at_fpr

    background = [float(value) for value in range(100)]  # 0..99
    foreground = [95.0, 99.0, 50.0, 80.0, 98.5]
    result = count_motif_matches_at_fpr(foreground, background, fpr=0.1, verbose=False)
    # 90th percentile of 0..99 is 89.1; background values 90..99 pass, foreground 95, 99, 98.5 pass.
    assert result["background_matches"] == 10 and result["foreground_matches"] == 3
    assert result["foreground_fraction"] == pytest.approx(0.6)
    # the generic metric wrapper cannot return the two counts and fails on this upstream function
    with pytest.raises(TypeError):
        motif_score_metric(foreground, background, metric="matches_at_fpr", fpr=0.1, verbose=False)
    for bad in (0.0, 1.5, True):
        with pytest.raises(ConfigurationError):
            count_motif_matches_at_fpr(foreground, background, fpr=bad, verbose=False)
    with pytest.raises(ConfigurationError):
        count_motif_matches_at_fpr([], background, verbose=False)


def test_gm_partial_roc_auc_by_hand_and_broken_generic_path() -> None:
    pytest.importorskip("gimmemotifs")
    from dnakit.exceptions import BackendExecutionError
    from dnakit.regulation import motif_score_metric
    from dnakit.regulation.gimmemotifs_family_functions import calculate_partial_motif_roc_auc

    foreground = [3.0, 4.0, 5.0]
    background = [1.0, 2.0, 3.5, 4.5]
    # ROC points (fpr, tpr): (0, 0) (0, 1/3) (1/4, 1/3) (1/4, 2/3) (1/2, 2/3) -> area to 0.5:
    # 0.25 * 1/3 + 0.25 * 2/3 = 0.25
    result = calculate_partial_motif_roc_auc(foreground, background, fpr_limit=0.5, verbose=False)
    assert result["area"] == pytest.approx(0.25)
    assert result["area_fraction_of_limit"] == pytest.approx(0.5)
    with pytest.raises(BackendExecutionError):
        motif_score_metric(foreground, background, metric="roc_auc_xlim", fpr=0.5, verbose=False)
    for bad in (0.0, 2.0):
        with pytest.raises(ConfigurationError):
            calculate_partial_motif_roc_auc(foreground, background, fpr_limit=bad, verbose=False)
    with pytest.raises(ConfigurationError):
        calculate_partial_motif_roc_auc([float("nan")], background, verbose=False)


# ------------------------------------------------------------------- MACS3
_MA_LEFT = [10, 20, 40, 100]
_MA_RIGHT = [50, 60, 90, 150]
_MA_COUNTS = [2, 1, 3, 1]


def test_ma_counted_fragment_pileup_by_hand() -> None:
    pytest.importorskip("MACS3")
    from dnakit.epigenome.MACS3_family_functions import build_counted_fragment_pileup

    result = build_counted_fragment_pileup(_MA_LEFT, _MA_RIGHT, _MA_COUNTS, verbose=False)
    # segments (10,20] 2, (20,40] 2+1, (40,50] 2+1+3, (50,60] 1+3, (60,90] 3, (90,100] 0,
    # (100,150] 1
    assert result["positions"] == [10, 20, 40, 50, 60, 90, 100, 150]
    assert result["values"] == [0.0, 2.0, 3.0, 6.0, 4.0, 3.0, 0.0, 1.0]
    assert result["total_count"] == 7
    scaled = build_counted_fragment_pileup(
        _MA_LEFT, _MA_RIGHT, _MA_COUNTS, scale_factor=0.5, left_sorted=True, verbose=False
    )
    assert scaled["values"] == [0.0, 1.0, 1.5, 3.0, 2.0, 1.5, 0.0, 0.5]
    with pytest.raises(ConfigurationError):
        build_counted_fragment_pileup([20, 10], [30, 40], [1, 1], left_sorted=True, verbose=False)
    with pytest.raises(ConfigurationError):
        build_counted_fragment_pileup(_MA_LEFT, _MA_RIGHT, [1, 2, 3], verbose=False)
    with pytest.raises(ConfigurationError):
        build_counted_fragment_pileup(_MA_LEFT, _MA_RIGHT, [70000, 1, 1, 1], verbose=False)


def test_ma_fragment_end_pileup_by_hand() -> None:
    pytest.importorskip("MACS3")
    from dnakit.epigenome.MACS3_family_functions import build_fragment_end_pileup

    result = build_fragment_end_pileup([100], [200], [2], window=20, verbose=False)
    # windows [90, 110) and [190, 210), each of weight 2
    assert result["positions"] == [90, 110, 190, 210]
    assert result["values"] == [0.0, 2.0, 0.0, 2.0]
    with pytest.raises(ConfigurationError):
        build_fragment_end_pileup([100], [200], [2], window=0, verbose=False)


def test_ma_hmmratac_fragment_pileup_by_hand() -> None:
    pytest.importorskip("MACS3")
    from dnakit.epigenome.MACS3_family_functions import build_hmmratac_fragment_pileup

    result = build_hmmratac_fragment_pileup(
        _MA_LEFT, _MA_RIGHT, {"40": 0.5, "50": 2.0}, verbose=False
    )
    # fragments of length 40 weigh 0.5 and of length 50 weigh 2.0
    assert result["positions"] == [10, 20, 40, 50, 60, 90, 100, 150]
    assert result["values"] == [0.0, 0.5, 1.0, 3.0, 2.5, 2.0, 0.0, 2.0]
    with pytest.raises(ConfigurationError):
        build_hmmratac_fragment_pileup(_MA_LEFT, _MA_RIGHT, {"40": 0.5}, verbose=False)
    with pytest.raises(ConfigurationError):
        build_hmmratac_fragment_pileup(_MA_LEFT, _MA_RIGHT, {"x": 0.5}, verbose=False)


def test_ma_allele_specific_model_matches_upstream_and_hand_gq() -> None:
    pytest.importorskip("MACS3")
    import math

    import numpy as np
    from MACS3.Signal import VariantStat

    from dnakit.epigenome.MACS3_family_functions import (
        calculate_allele_specific_genotype_quality,
        calculate_allele_specific_model_likelihood,
    )

    reads = ([30, 32, 35, 28, 30, 33], [2, 3, 1, 2, 3, 2], [0] * 6, [0] * 6)
    arrays = [np.array(item, dtype=np.int32) for item in reads]
    direct = VariantStat.CalModel_Heter_AS(*arrays, 0.9)
    result = calculate_allele_specific_model_likelihood(
        *reads, max_allowed_ratio=0.9, verbose=False
    )
    assert (result["log_likelihood"], result["bic"]) == pytest.approx(direct)
    free = calculate_allele_specific_model_likelihood(*reads, verbose=False)
    assert free["log_likelihood"] != result["log_likelihood"]
    for first, second in ((-10.0, -15.0), (-15.0, -10.0), (-10.0, -10.0), (-5.0, -400.0)):
        ratio = max(1e-110, min(1.0, math.exp(second - first)))
        share = ratio / (1 + ratio)
        expected = int(-4.34294 * math.log(share)) if share > 1e-110 else 255
        score = calculate_allele_specific_genotype_quality(first, second, verbose=False)["score"]
        assert score == expected
    assert calculate_allele_specific_genotype_quality(-10, -300, verbose=False)["score"] == 255
    with pytest.raises(ConfigurationError):
        calculate_allele_specific_model_likelihood(*reads, max_allowed_ratio=2.0, verbose=False)


def test_ma_family_heterozygous_model_reaches_noas() -> None:
    pytest.importorskip("MACS3")
    import numpy as np
    from MACS3.Signal import VariantStat

    from dnakit.epigenome import calculate_allele_model_likelihood

    assert "model" in inspect.signature(calculate_allele_model_likelihood).parameters
    reads = (
        [20, 22, 25, 18, 20, 23],
        [15, 13, 17, 14, 12, 16],
        [2, 1, 2, 1, 1, 2],
        [1, 1, 1, 0, 1, 1],
    )
    direct = VariantStat.CalModel_Heter_noAS(*[np.array(item, dtype=np.int32) for item in reads])
    wrapped = calculate_allele_model_likelihood(*reads, model="heterozygous", verbose=False)
    assert (wrapped["log_likelihood"], wrapped["bic"]) == pytest.approx(direct)


def test_ma_family_pileup_variants_equal_existing_wrappers() -> None:
    pytest.importorskip("MACS3")
    import numpy as np
    from MACS3.Signal import PileupV2
    from MACS3.Signal.Pileup import se_all_in_one_pileup

    from dnakit.epigenome import build_paired_fragment_pileup, build_shifted_read_pileup
    from dnakit.epigenome.MACS3_family_functions import build_counted_fragment_pileup

    plus, minus = [10, 50, 100], [30, 80, 300]
    pa, na = np.array(plus, dtype=np.int32), np.array(minus, dtype=np.int32)
    # se_all_in_one_pileup gives the same pileup as the shifted wrapper
    shifted = build_shifted_read_pileup(
        plus, minus, five_shift=-25, three_shift=75, chrom_length=1000, verbose=False
    )
    positions, values = se_all_in_one_pileup(pa, na, -25, 75, 1000, 1.0, 0.0)
    assert shifted["positions"] == positions.tolist() and shifted["values"] == values.tolist()
    # pileup_from_PN(extsize) equals five_shift=0, three_shift=extsize apart from an unclipped
    # leading zero point that MACS3 adds when a minus read reaches below 0
    ext = PileupV2.pileup_from_PN(pa, na, 60)
    fixed = build_shifted_read_pileup(
        plus, minus, five_shift=0, three_shift=60, chrom_length=1_000_000, verbose=False
    )
    upstream = [(int(p), float(v)) for p, v in ext]
    assert upstream[0] == (-30, 0.0)
    assert upstream[1:] == list(zip(fixed["positions"], fixed["values"], strict=True))
    # pileup_from_LR returns the same pileup as an array of (p, v) records
    array = np.array(
        list(zip(_MA_LEFT, _MA_RIGHT, strict=True)), dtype=[("l", "<i4"), ("r", "<i4")]
    )
    track = PileupV2.pileup_from_LR(array)
    paired = build_paired_fragment_pileup(_MA_LEFT, _MA_RIGHT, verbose=False)
    assert [int(p) for p, _v in track] == paired["positions"]
    assert [float(v) for _p, v in track] == paired["values"]
    # pileup_from_LRC likewise equals the counted wrapper
    counted_array = np.array(
        list(zip(_MA_LEFT, _MA_RIGHT, _MA_COUNTS, strict=True)),
        dtype=[("l", "<i4"), ("r", "<i4"), ("c", "<u2")],
    )
    counted = build_counted_fragment_pileup(_MA_LEFT, _MA_RIGHT, _MA_COUNTS, verbose=False)
    track = PileupV2.pileup_from_LRC(counted_array)
    assert [int(p) for p, _v in track] == counted["positions"]
    assert [float(v) for _p, v in track] == counted["values"]


def test_ma_aliases_are_same_objects() -> None:
    pytest.importorskip("MACS3")
    from MACS3.Signal import PairedEndTrack, PileupV2

    for name in (
        "pileup_from_LR_hmmratac",
        "pileup_from_LRC",
        "pileup_from_LRC_as_list",
        "pileup_from_LRC_centers_as_list",
    ):
        assert getattr(PairedEndTrack, name) is getattr(PileupV2, name)


# --------------------------------------------------------------- cooltools
def _floats(values: list[Any]) -> Any:
    import numpy as np

    return np.array([np.nan if value is None else value for value in values], dtype=float)


def _cooler_file(tmp_path: Path, name: str = "synthetic.cool") -> Path:
    """Two chromosomes (1500 and 2000 bp) at 100 bp, balanced."""

    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(0)
    chromsizes = pd.Series({"chr1": 1500, "chr2": 2000}, name="length")
    bins = cooler.util.binnify(chromsizes, 100)
    rows = []
    for i in range(len(bins)):
        for j in range(i, len(bins)):
            distance = abs(i - j)
            same = bins.chrom[i] == bins.chrom[j]
            count = rng.poisson(6.0 / (1 + distance) if same else 0.5)
            if distance == 0:
                count += 3
            if count > 0:
                rows.append((i, j, int(count)))
    pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
    path = tmp_path / name
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
    cooler.balance_cooler(cooler.Cooler(str(path)), store=True)
    return path


def test_ct_pairwise_region_sums_has_diagonals_and_block_rows(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    from dnakit.epigenome import pairwise_diagonal_sums
    from dnakit.epigenome.cooltools_family_functions import pairwise_region_sums

    path = _cooler_file(tmp_path)
    result = pairwise_region_sums(path, verbose=False)
    diagonal = pairwise_diagonal_sums(path, verbose=False)
    intra = [row for row in result["rows"] if row["dist"] >= 0]
    assert intra == diagonal["rows"]
    blocks = [row for row in result["rows"] if row["dist"] == -1]
    assert result["block_row_count"] == len(blocks) == 1
    # the block row is the sum of all chr1-chr2 pixel counts and the number of valid bin pairs
    clr = cooler.Cooler(str(path))
    pixels = clr.pixels()[:]
    chrom = clr.bins()[:]["chrom"].to_numpy()
    cross = pixels[
        (chrom[pixels["bin1_id"].to_numpy()] == "chr1")
        & (chrom[pixels["bin2_id"].to_numpy()] == "chr2")
    ]
    assert blocks[0]["region1"] == "chr1" and blocks[0]["region2"] == "chr2"
    assert blocks[0]["count.sum"] == int(cross["count"].sum())
    assert blocks[0]["n_valid"] == 15 * 20
    assert blocks[0]["n_total"] is None
    with pytest.raises(ConfigurationError):
        pairwise_region_sums(path, weight="nope", verbose=False)


def test_ct_region_contact_scaling_by_hand() -> None:
    pytest.importorskip("cooltools")
    pytest.importorskip("pandas")
    from dnakit.epigenome.cooltools_family_functions import compute_region_contact_scaling

    result = compute_region_contact_scaling(
        [100, 150, 200, 1000, 5000, 200],
        [130, 400, 1000, 1500, 9000, 260],
        region1=[0, 10000],
        min_distance=10,
        max_distance=10000,
        n_bins=5,
        verbose=False,
    )
    # distances 30, 250, 800, 500, 4000 and 60 fall in bins [10,56) [56,316) [316,1778) [1778,10000)
    edges = [10, 56, 316, 1778, 10000]
    assert [item["min_dist"] for item in result["bins"]] == edges[:-1]
    assert [item["observed"] for item in result["bins"]] == [1, 2, 2, 1]
    assert result["contacts_counted"] == 6
    length = 10000
    for item in result["bins"]:
        low, high = item["min_dist"], item["max_dist"]
        area = (high - low) * length - (high**2 - low**2) / 2
        assert item["area"] == pytest.approx(area)
        assert item["frequency"] == pytest.approx(item["observed"] / area)
    restricted = compute_region_contact_scaling(
        [100, 150, 200, 1000, 5000, 200],
        [130, 400, 1000, 1500, 9000, 260],
        region1=[0, 2000],
        region2=[0, 10000],
        min_distance=10,
        max_distance=10000,
        n_bins=5,
        verbose=False,
    )
    assert [item["observed"] for item in restricted["bins"]] == [1, 2, 2, 0]
    for bad in ({"region1": [5, 5]}, {"region1": [0]}, {"n_bins": 1}, {"max_distance": 10}):
        arguments = {"region1": [0, 100], "min_distance": 10, "max_distance": 1000, **bad}
        with pytest.raises(ConfigurationError):
            compute_region_contact_scaling([1, 2], [3, 4], verbose=False, **arguments)


def test_ct_family_smoothing_sampling_balancing_and_expected_reach_upstream(
    tmp_path: Path,
) -> None:
    pytest.importorskip("cooltools")
    import numpy as np
    from cooltools.api import expected, sample
    from cooltools.lib import numutils

    from dnakit.epigenome import (
        expected_full_table,
        iterative_correction_dense,
        sample_contacts,
        smooth_contact_decay,
    )

    path = _cooler_file(tmp_path)
    assert "mode" in inspect.signature(smooth_contact_decay).parameters
    assert "fast" in inspect.signature(expected_full_table).parameters
    assert {"count", "fraction"} <= set(inspect.signature(sample_contacts).parameters)
    assert "symmetric" in inspect.signature(iterative_correction_dense).parameters

    cooler = pytest.importorskip("cooler")
    from cooltools.lib import common

    clr = cooler.Cooler(str(path))
    view = common.make_cooler_view(clr)
    cvd = expected.expected_cis(clr, view_df=view, intra_only=True, smooth=False)
    direct = expected.per_region_smooth_cvd(
        cvd, cols=expected.DEFAULT_BALANCED_CVD_COLS, suffix=".smoothed"
    )
    wrapped = smooth_contact_decay(path, mode="per_region", verbose=False)
    column = "balanced.avg.smoothed"
    direct = direct.sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    assert np.allclose(
        _floats([row[column] for row in wrapped["rows"]]),
        direct[column].to_numpy(),
        equal_nan=True,
    )
    assert np.isfinite(direct[column].to_numpy()).any()

    # balancing without symmetry
    matrix = np.array([[2.0, 1.0, 0.5], [1.0, 3.0, 1.0], [0.5, 2.0, 4.0]])
    direct_corrected, row_bias, _col_bias, _ = numutils.iterative_correction_asymmetric(
        matrix, max_iter=100, tol=1e-6
    )
    result = iterative_correction_dense(
        matrix.tolist(), symmetric=False, max_iter=100, tol=1e-6, verbose=False
    )
    assert np.allclose(np.array(result["corrected"]), direct_corrected, equal_nan=True)
    assert np.allclose(_floats(result["bias"]["row"]), row_bias, equal_nan=True)

    # sampling: seeded wrapper output equals the pixel-level function with the same seed
    pixels = clr.pixels()[:]
    np.random.seed(5)
    exact = sample.sample_pixels_exact(pixels.copy(), 100)
    out = tmp_path / "exact.cool"
    sample_contacts(path, out, count=100, seed=5, verbose=False)
    sampled = cooler.Cooler(str(out)).pixels()[:]
    assert sampled["count"].tolist() == exact["count"].tolist()
    np.random.seed(5)
    approx = sample.sample_pixels_approx(pixels.copy(), 0.5)
    out2 = tmp_path / "approx.cool"
    sample_contacts(path, out2, fraction=0.5, seed=5, verbose=False)
    sampled2 = cooler.Cooler(str(out2)).pixels()[:]
    assert sampled2["count"].tolist() == approx["count"].tolist()


def test_ct_family_fast_expected_and_oe_chunks_equal_wrappers(tmp_path: Path) -> None:
    cooler = pytest.importorskip("cooler")
    pytest.importorskip("cooltools")
    import numpy as np
    from cooltools.lib import common
    from cooltools.sandbox import obs_over_exp_cooler as sandbox

    from dnakit.epigenome import expected_full_table, observed_over_expected_pixels

    path = _cooler_file(tmp_path)
    clr = cooler.Cooler(str(path))
    view = common.make_cooler_view(clr)
    fast = expected_full_table(path, fast=True, verbose=False)
    direct = sandbox.expected_full_fast(clr, view_df=view, clr_weight_name="weight", nproc=1)
    direct = direct.sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    assert np.allclose(
        _floats([row["expected"] for row in fast["rows"]]),
        direct["expected"].to_numpy(),
        equal_nan=True,
    )
    expected = sandbox.expected_full(clr, view_df=view, clr_weight_name="weight", nproc=1)
    chunks = list(
        sandbox.obs_over_exp_generator(
            clr, expected, view_df=view, chunksize=50, oe_column_name="oe"
        )
    )
    assert len(chunks) > 1
    joined = [value for chunk in chunks for value in chunk["oe"].tolist()]
    wrapped = [row["oe"] for row in observed_over_expected_pixels(path, verbose=False)["pixels"]]
    assert np.allclose(np.array(joined, dtype=float), _floats(wrapped), equal_nan=True)


# ---------------------------------------------------------------- pairtools
_PT_HEADER = (
    "## pairs format v1.0.0\n#shape: upper triangle\n#chromsize: chr1 10000\n"
    "#chromsize: chr2 10000\n"
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
)


def _pairs_file(tmp_path: Path) -> Path:
    rows = [
        ("r1", "chr1", 100, "chr2", 500, "+", "-"),
        ("r2", "chr1", 102, "chr2", 500, "+", "-"),
        ("r3", "chr1", 104, "chr2", 500, "+", "-"),
        ("r4", "chr1", 106, "chr2", 500, "+", "-"),
        ("r5", "chr1", 100, "chr2", 502, "+", "-"),
        ("r6", "chr1", 3000, "chr1", 4000, "+", "-"),
        ("r7", "chr1", 3001, "chr1", 4002, "+", "-"),
        ("r8", "chr1", 3000, "chr1", 4000, "-", "+"),
    ]
    rows.sort(key=lambda row: (row[1], row[2]))
    path = tmp_path / "a.pairs"
    path.write_text(_PT_HEADER + "".join("\t".join(map(str, row)) + "\tUU\n" for row in rows))
    return path


def test_pt_deduplicate_online_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    from dnakit.epigenome import deduplicate_pairs
    from dnakit.epigenome.pairtools_family_functions import deduplicate_pairs_online

    path = _pairs_file(tmp_path)
    result = deduplicate_pairs_online(
        path,
        max_mismatch=2,
        output_path=tmp_path / "unique.pairs",
        duplicates_path=tmp_path / "dups.pairs",
        mark_duplicates=True,
        verbose=False,
    )
    # kept: r1 (100), r3 (104; 4 bp from r1), r6, and r8 (other strands); r2 is 2 bp from r1,
    # r4 is 2 bp from r3, r5 differs by 2 bp on the second side, r7 by 1 and 2 bp from r6
    assert result["unique_read_ids"] == ["r1", "r3", "r6", "r8"]
    assert sorted(result["duplicate_read_ids"]) == ["r2", "r4", "r5", "r7"]
    assert (result["total_nodups"], result["total_dups"]) == (4, 4)
    dup_lines = (tmp_path / "dups.pairs").read_text().splitlines()
    assert all(line.endswith("\tDD") for line in dup_lines if not line.startswith("#"))
    assert (tmp_path / "unique.pairs").read_text().startswith("## pairs format")
    # the default clustering joins chains of close pairs, so it keeps fewer pairs
    chained = deduplicate_pairs(path, max_mismatch=2, verbose=False)
    assert (chained["total_nodups"], chained["total_dups"]) == (3, 5)
    with pytest.raises(ConfigurationError):
        deduplicate_pairs_online(path, output_path=tmp_path / "unique.pairs", verbose=False)
    with pytest.raises(ConfigurationError):
        deduplicate_pairs_online(path, method="min", verbose=False)


def test_pt_deduplicate_online_equals_direct_call(tmp_path: Path) -> None:
    pytest.importorskip("pairtools")
    import io

    from pairtools.lib import dedup, pairsam_format
    from pairtools.lib.stats import PairCounter

    from dnakit.epigenome.pairtools_family_functions import deduplicate_pairs_online

    path = _pairs_file(tmp_path)
    lines = [line for line in path.read_text().splitlines(keepends=True) if line[0] != "#"]
    for method, mismatch in (("max", 3), ("sum", 2)):
        counter = PairCounter(
            bytile_dups=False, filters=None, startup_code="", type_cast=(), engine="pandas"
        )
        dedup.streaming_dedup_cython(
            method, mismatch, "\t", 1, 3, 2, 4, 5, 6, [], [], pairsam_format.UNMAPPED_CHROM,
            iter(lines), io.StringIO(), io.StringIO(), io.StringIO(), counter, False,
        )  # fmt: skip
        result = deduplicate_pairs_online(path, max_mismatch=mismatch, method=method, verbose=False)
        assert (result["total_nodups"], result["total_dups"]) == (
            counter["total_nodups"],
            counter["total_dups"],
        )


def test_pt_unity_normalisation_by_hand() -> None:
    pytest.importorskip("pairtools")
    from dnakit.epigenome import normalize_scaling_curve
    from dnakit.epigenome.pairtools_family_functions import normalize_scaling_curve_unity

    edges = [10, 100, 1000, 10000]
    frequencies = [0.5, 0.1, 0.01]
    result = normalize_scaling_curve_unity(edges, frequencies, norm_range=None, verbose=False)
    area = 90 * 0.5 + 900 * 0.1 + 9000 * 0.01  # 225
    assert [item["frequency"] for item in result["bins"]] == pytest.approx(
        [value / area for value in frequencies]
    )
    # only bins whose midpoint (31.6, 316, 3162) is inside (100, 5000) count: the last two
    windowed = normalize_scaling_curve_unity(
        edges, frequencies, norm_range=[100, 5000], verbose=False
    )
    area2 = 900 * 0.1 + 9000 * 0.01
    assert windowed["bins"][0]["frequency"] == pytest.approx(0.5 / area2)
    mean_scaled = normalize_scaling_curve(
        edges[1:], frequencies, norm_window=[100, 10000], verbose=False
    )
    assert mean_scaled["points"][0]["frequency"] != pytest.approx(0.5 / area)
    for bad in (
        {"bin_edges": [10, 100], "frequencies": [1.0, 2.0]},
        {"bin_edges": [100, 10, 5], "frequencies": [1.0, 2.0]},
    ):
        with pytest.raises(ConfigurationError):
            normalize_scaling_curve_unity(verbose=False, **bad)


def test_pt_family_xa_phase_reaches_upstream() -> None:
    pytest.importorskip("pairtools")
    from pairtools.lib import phase

    from dnakit.epigenome import assign_haplotype_phase

    assert {"variant", "mismatches"} <= set(inspect.signature(assign_haplotype_phase).parameters)
    tag = "chr1_2,+100,50M,2;chr1_1,+100,50M,3;"
    direct = phase.phase_side_XA("chr1_1", tag, 40, 41, 1, ("_1", "_2"))
    result = assign_haplotype_phase(
        "chr1_1",
        tag_value=tag,
        primary_score=40,
        secondary_score=41,
        variant="XA",
        mismatches=1,
        verbose=False,
    )
    assert (result["phase"], result["chromosome_base"], *result["scores"]) == tuple(direct)


def test_pt_aliases_are_same_objects() -> None:
    pytest.importorskip("pairtools")
    from pairtools.cli import dedup as cli_dedup
    from pairtools.cli import phase as cli_phase
    from pairtools.lib import dedup, phase

    assert cli_dedup.streaming_dedup_cython is dedup.streaming_dedup_cython
    assert cli_phase.phase_side_XA is phase.phase_side_XA


# -------------------------------------------------------------- tangermeme
def _tm_linear_model(length: int = 20) -> Any:
    torch = pytest.importorskip("torch")
    weights = torch.tensor(
        [[((base * 7 + pos * 3) % 11) / 4.0 - 1.0 for pos in range(length)] for base in range(4)],
        dtype=torch.float32,
    )

    class Linear(torch.nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.register_buffer("weights", weights)

        def forward(self, x: Any, *args: Any) -> Any:
            return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

    return Linear().eval()


_TM_SEQUENCE = "ACGTTGATAAGCCGTAGGAT"


def test_tm_insert_motif_by_hand() -> None:
    pytest.importorskip("tangermeme")
    from dnakit.regulation import place_motif
    from dnakit.regulation.tangermeme_family_functions import insert_motif

    result = insert_motif({"a": _TM_SEQUENCE, "b": "T" * 20}, "GG", start=3, verbose=False)
    by_name = {item["name"]: item["sequence"] for item in result["sequences"]}
    assert by_name["a"] == "ACG" + "GG" + _TM_SEQUENCE[3:]
    assert by_name["b"] == "TTT" + "GG" + "T" * 17
    assert (result["input_length"], result["output_length"]) == (20, 22)
    trimmed = place_motif({"a": _TM_SEQUENCE}, "GG", start=3, mode="insert", verbose=False)
    assert trimmed["sequences"][0]["sequence"] == by_name["a"][:20]
    centred = insert_motif({"a": _TM_SEQUENCE}, "GG", verbose=False)
    assert centred["sequences"][0]["sequence"] == _TM_SEQUENCE[:10] + "GG" + _TM_SEQUENCE[10:]
    with pytest.raises(ConfigurationError):
        insert_motif({"a": _TM_SEQUENCE}, "GG", start=21, verbose=False)
    with pytest.raises(ConfigurationError):
        insert_motif({"a": _TM_SEQUENCE}, "GX", verbose=False)


def test_tm_design_beam_substitution_matches_upstream_and_greedy() -> None:
    pytest.importorskip("tangermeme")
    import torch
    from tangermeme.design import beam_substitution
    from tangermeme.utils import characters, one_hot_encode

    from dnakit.regulation import design_greedy_substitution
    from dnakit.regulation.tangermeme_family_functions import design_beam_substitution

    model = _tm_linear_model()
    result = design_beam_substitution(
        model,
        _TM_SEQUENCE,
        [5.0],
        motifs=["GGG", "AC"],
        beam_size=4,
        n_best=3,
        max_iter=3,
        verbose=False,
    )
    direct = beam_substitution(
        model,
        one_hot_encode(_TM_SEQUENCE)[None].float(),
        y=torch.tensor([[5.0]]),
        motifs=["GGG", "AC"],
        beam_size=4,
        n_best=3,
        max_iter=3,
        tol=0.001,
    )
    assert result["designed_sequences"] == [characters(item) for item in direct]
    assert result["designed_sequence"] == result["designed_sequences"][0]
    # a beam of one is the greedy search
    single = design_beam_substitution(
        model, _TM_SEQUENCE, [5.0], motifs=["GGG"], beam_size=1, max_iter=2, verbose=False
    )
    greedy = design_greedy_substitution(
        model, _TM_SEQUENCE, [5.0], motifs=["GGG"], max_iter=2, verbose=False
    )
    assert single["designed_sequence"] == greedy["designed_sequence"]
    with pytest.raises(ConfigurationError):
        design_beam_substitution(model, _TM_SEQUENCE, [5.0], beam_size=2, n_best=3, verbose=False)
    with pytest.raises(ConfigurationError):
        design_beam_substitution(model, _TM_SEQUENCE, [], verbose=False)


def test_tm_family_spacing_reaches_upstream_and_aliases() -> None:
    pytest.importorskip("tangermeme")
    import tangermeme.design as design
    import torch
    from tangermeme import annotate, ersatz, variant_effect

    from dnakit.regulation import pairwise_annotation_counts

    assert "max_distance" in inspect.signature(pairwise_annotation_counts).parameters
    rows = [
        {"example_idx": 0, "annotation": 0, "start": 0, "end": 5},
        {"example_idx": 0, "annotation": 1, "start": 8, "end": 12},
        {"example_idx": 0, "annotation": 1, "start": 20, "end": 25},
        {"example_idx": 1, "annotation": 0, "start": 3, "end": 6},
        {"example_idx": 1, "annotation": 1, "start": 6, "end": 9},
    ]
    result = pairwise_annotation_counts(rows, max_distance=12, verbose=False)
    table = torch.tensor(
        [[row["example_idx"], row["annotation"], row["start"], row["end"]] for row in rows]
    )
    direct = annotate.pairwise_annotations_spacing(
        table, max_distance=12, dtype=torch.int64, symmetric=True
    )
    assert result["mode"] == "spacing" and result["counts"] == direct.tolist()
    assert variant_effect.insert is ersatz.insert
    import importlib

    module = importlib.import_module("tangermeme.design.beam_substitution")
    assert module.beam_substitution is design.beam_substitution


# ---------------------------------------------------------------- deeptools
def test_dt_chop_regions_from_middle_by_hand() -> None:
    pytest.importorskip("deeptools")
    from deeptools.heatmapper import chopRegionsFromMiddle

    from dnakit.epigenome.deeptools_family_functions import chop_regions_from_middle

    exons = [[0, 200], [300, 400], [800, 900]]
    # 400 exonic bases: the centre falls after 200 bases, at the end of the first exon
    result = chop_regions_from_middle(exons, left=150, right=80, verbose=False)
    assert result["left_bins"] == [[50, 200]]
    assert result["right_bins"] == [[300, 380]]
    assert (result["pad_left"], result["pad_right"]) == (0, 0)
    short = chop_regions_from_middle([[0, 50]], left=100, right=10, verbose=False)
    # centre after 25 bases: 25 bases on the left (75 missing), right side keeps 10
    assert short["left_bins"] == [[0, 25]] and short["pad_left"] == 75
    assert short["right_bins"] == [[25, 35]]
    direct = chopRegionsFromMiddle([(0, 200), (300, 400), (800, 900)], 150, 80)
    assert [[list(item) for item in direct[0]], [list(item) for item in direct[1]]] == [
        result["left_bins"],
        result["right_bins"],
    ]
    with pytest.raises(ConfigurationError):
        chop_regions_from_middle([[10, 5]], verbose=False)


def test_dt_aliases_are_same_objects() -> None:
    pytest.importorskip("deeptools")
    from deeptools import computeGCBias, correctGCBias, utilities

    from dnakit.epigenome import gc_fraction_from_twobit, map_twobit_to_bam_chromosomes

    assert "fraction" in inspect.signature(gc_fraction_from_twobit).parameters
    assert map_twobit_to_bam_chromosomes is not None
    for module in (computeGCBias, correctGCBias):
        assert module.tbitToBamChrName is utilities.tbitToBamChrName
        assert module.getGC_content is utilities.getGC_content


# ------------------------------------------------------------------ varcode
_VC_MAF_COLUMNS = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2",
]


def test_vc_maf_and_vcf_tables(tmp_path: Path) -> None:
    pytest.importorskip("varcode")
    from dnakit.variants.varcode_family_functions import load_maf_table, read_vcf_table

    row = dict.fromkeys(_VC_MAF_COLUMNS, "")
    row.update(
        Hugo_Symbol="TP53",
        NCBI_Build="GRCh38",
        Chromosome="17",
        Start_Position="7676154",
        End_Position="7676154",
        Strand="+",
        Variant_Classification="Missense_Mutation",
        Variant_Type="SNP",
        Reference_Allele="C",
        Tumor_Seq_Allele1="C",
        Tumor_Seq_Allele2="T",
        Tumor_Sample_Barcode="TCGA-AA-0001-01A",
    )
    maf = tmp_path / "t.maf"
    maf.write_text(
        "\t".join(_VC_MAF_COLUMNS) + "\n" + "\t".join(row[c] for c in _VC_MAF_COLUMNS) + "\n"
    )
    table = load_maf_table(maf, verbose=False)
    assert table["columns"] == _VC_MAF_COLUMNS and table["row_count"] == 1
    first = table["rows"][0]
    assert first["Hugo_Symbol"] == "TP53" and first["Variant_Classification"] == "Missense_Mutation"
    assert first["Start_Position"] == 7676154 and first["Entrez_Gene_Id"] is None
    short = tmp_path / "short.maf"
    short.write_text("Hugo_Symbol\tChromosome\nTP53\t17\n")
    with pytest.raises(ConfigurationError):
        load_maf_table(short, verbose=False)

    vcf = tmp_path / "t.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "1\t100\t.\tA\tG,T\t50\tPASS\tDP=10;AF=0.1,0.2\tGT:GQ\t0/1:35\t1/2:40\n"
        "1\t200\trs1\tCT\tC\t60\tPASS\tDP=20\tGT:GQ\t0/0:99\t0/1:12\n"
    )
    plain = read_vcf_table(vcf, verbose=False)
    assert plain["columns"] == ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    assert plain["rows"][0]["ALT"] == "G,T" and plain["rows"][0]["POS"] == 100
    full = read_vcf_table(vcf, include_info=True, sample_names=["S1", "S2"], verbose=False)
    assert full["rows"][0]["S2"] == "1/2:40" and full["rows"][1]["S1"] == "0/0:99"
    assert full["rows"][0]["INFO"] == "DP=10;AF=0.1,0.2"
    one = read_vcf_table(vcf, include_info=True, sample_names=["S1"], verbose=False)
    assert one["columns"][-1] == "S1" and one["rows"][1]["S1"] == "0/0:99"
    for bad in ({"sample_names": ["S2"], "include_info": True}, {"sample_names": ["S1"]}):
        with pytest.raises(ConfigurationError):
            read_vcf_table(vcf, verbose=False, **bad)


def test_vc_family_acceptor_and_prefix_trim_reach_upstream() -> None:
    pytest.importorskip("varcode")
    from varcode import cryptic_exons
    from varcode.string_helpers import trim_shared_flanking_strings, trim_shared_prefix

    from dnakit.variants import trim_shared_allele_flanks
    from dnakit.variants.varcode_functions import score_splice_site

    assert "site" in inspect.signature(score_splice_site).parameters
    result = score_splice_site("CAGG", site="acceptor", verbose=False)
    assert result["score"] == pytest.approx(cryptic_exons.score_acceptor("CAGG"), abs=1e-6)
    # prefix-only trimming equals the flank-trimmed alleles with the shared suffix put back
    trimmed = trim_shared_allele_flanks("ACGT", "ACTT", verbose=False)
    ref_rest, alt_rest, prefix = trim_shared_prefix("ACGT", "ACTT")
    assert (trimmed["trimmed_ref"], trimmed["trimmed_alt"]) == ("G", "T")
    assert (ref_rest, alt_rest, prefix) == (
        trimmed["trimmed_ref"] + trimmed["suffix"],
        trimmed["trimmed_alt"] + trimmed["suffix"],
        trimmed["prefix"],
    )
    assert trim_shared_flanking_strings("ACGT", "ACTT")[:2] == ("G", "T")


# --------------------------------------------------------------- bed_reader
def test_br_family_float32_kernels_reach_upstream() -> None:
    pytest.importorskip("bed_reader")
    import bed_reader
    import numpy as np

    from dnakit.variants import standardize_genotypes, subset_genotype_matrix

    assert "dtype" in inspect.signature(standardize_genotypes).parameters
    assert "dtype" in inspect.signature(subset_genotype_matrix).parameters
    genotypes = [[0.0, 1.0, 2.0], [1.0, 1.0, 0.0], [2.0, 0.0, 1.0], [0.0, 2.0, 2.0]]
    values = np.asfortranarray(np.array(genotypes, dtype=np.float32))
    stats = np.zeros((3, 2), dtype=np.float32, order="F")
    bed_reader.standardize_f32(values, False, 1.0, 25.0, True, False, stats, 1)
    wrapped = standardize_genotypes(genotypes, dtype="float32", verbose=False)
    assert np.allclose(np.array(wrapped["standardized"]), values, atol=1e-6)
    assert wrapped["dtype"] == "float32"
    # column 0 is [0, 1, 2, 0]: mean 0.75, population SD sqrt(0.6875)
    assert wrapped["variant_statistics"][0]["mean"] == pytest.approx(0.75)
    assert wrapped["variant_statistics"][0]["std"] == pytest.approx(0.6875**0.5, rel=1e-6)

    source = np.asfortranarray(np.array(genotypes, dtype=np.float32)[:, :, None])
    rows, columns = [2, 0], [1, 1, 2]
    target32 = np.zeros((2, 3, 1), dtype=np.float32, order="F")
    bed_reader.subset_f32_f32(
        source, np.asarray(rows, dtype=np.uintp), np.asarray(columns, dtype=np.uintp), target32, 1
    )
    target64 = np.zeros((2, 3, 1), dtype=np.float64, order="F")
    bed_reader.subset_f32_f64(
        source, np.asarray(rows, dtype=np.uintp), np.asarray(columns, dtype=np.uintp), target64, 1
    )
    wrapped_subset = subset_genotype_matrix(
        genotypes, sample_indices=rows, variant_indices=columns, dtype="float32", verbose=False
    )
    assert (
        wrapped_subset["genotypes"]
        == target32[:, :, 0].tolist()
        == [[0.0, 0.0, 1.0], [1.0, 1.0, 2.0]]
    )
    # the mixed-precision kernel gives exactly the same numbers as plain Python floats
    assert target64[:, :, 0].tolist() == wrapped_subset["genotypes"]


# ----------------------------------------------------------------- methylpy
def test_mp_paired_fastq_conversion_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    from dnakit.epigenome import convert_fastq_to_bisulfite
    from dnakit.epigenome.methylpy_family_functions import convert_paired_fastq_to_bisulfite

    fastq = tmp_path / "r.fq"
    fastq.write_text("@r1 x\nACCGTTCA\n+\nIIIIIIII\n@r2\nGGATCC\n+\nIIIIII\n")
    one = convert_paired_fastq_to_bisulfite(fastq, tmp_path / "one.fq", verbose=False)
    two = convert_paired_fastq_to_bisulfite(
        fastq, tmp_path / "two.fq", is_read2=True, verbose=False
    )
    lines_one = (tmp_path / "one.fq").read_text().splitlines()
    lines_two = (tmp_path / "two.fq").read_text().splitlines()
    assert (one["reads"], two["reads"]) == (2, 2)
    assert (one["converted_base"], two["converted_base"]) == ("C to T", "G to A")
    # read 1: C -> T; read 2: G -> A. Header spaces become "!" and an encoding is appended.
    assert (lines_one[1], lines_one[5]) == ("ATTGTTTA", "GGATTT")
    assert (lines_two[1], lines_two[5]) == ("ACCATTCA", "AAATCC")
    assert lines_one[0].startswith("@r1!x!") and lines_two[0].startswith("@r1!x!")
    # read 1 conversion equals the single-end converter
    convert_fastq_to_bisulfite(fastq, tmp_path / "se.fq", verbose=False)
    assert (tmp_path / "se.fq").read_text() == (tmp_path / "one.fq").read_text()
    with pytest.raises(ConfigurationError):
        convert_paired_fastq_to_bisulfite(fastq, tmp_path / "one.fq", verbose=False)
    with pytest.raises(ConfigurationError):
        convert_paired_fastq_to_bisulfite(fastq, tmp_path / "x.fq", is_read2=1, verbose=False)


def test_mp_family_mini_batch_merge_reaches_upstream(tmp_path: Path) -> None:
    pytest.importorskip("methylpy")
    from methylpy import utilities

    from dnakit.epigenome import merge_allc_tables

    assert "mini_batch" in inspect.signature(merge_allc_tables).parameters
    tables = {
        "a.tsv": "chr1\t10\t+\tCGA\t3\t5\t1\nchr1\t20\t+\tCGT\t1\t4\t0\n",
        "b.tsv": "chr1\t10\t+\tCGA\t2\t6\t1\n",
        "c.tsv": "chr1\t10\t+\tCGA\t1\t1\t1\nchr1\t30\t+\tCHH\t0\t3\t0\n",
    }
    paths = []
    for name, text in tables.items():
        (tmp_path / name).write_text(text)
        paths.append(tmp_path / name)
    result = merge_allc_tables(
        paths, tmp_path / "merged.tsv", mini_batch=2, compress_output=False, verbose=False
    )
    direct = tmp_path / "direct.tsv"
    utilities.merge_allc_files_minibatch(
        [str(path) for path in paths], str(direct), mini_batch=2, compress_output=False
    )
    merged = Path(result["output"]).read_text()
    # methylated 3 + 2 + 1 = 6 of 5 + 6 + 1 = 12 at position 10
    assert merged.splitlines()[0] == "chr1\t10\t+\tCGA\t6\t12\t1"
    assert merged == direct.read_text()


# ------------------------------------------------- pysam, hgvs, vcfpy, dendropy
def test_ps_export_vcf_to_phasing_formats_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    from dnakit.variants.pysam_family_functions import export_vcf_to_phasing_format

    vcf = tmp_path / "ph.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|1\n"
        "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0|0\t1|0\n"
    )
    gen = export_vcf_to_phasing_format(vcf, tmp_path / "g", fmt="gensample", verbose=False)
    # per sample the three columns are P(hom ref), P(het), P(hom alt)
    assert gen["rows"] == [
        "1:100_A_G 1:100_A_G 100 A G 0 1 0 0 0 1",
        "1:200_C_T 1:200_C_T 200 C T 1 0 0 0 1 0",
    ]
    assert gen["sample_ids"] == ["S1", "S2"] and gen["site_count"] == 2
    assert (tmp_path / "g.samples").is_file()
    hap = export_vcf_to_phasing_format(vcf, tmp_path / "h", fmt="hapsample", verbose=False)
    # two haplotype alleles per sample: S1 0|1 and S2 1|1 at the first site
    assert hap["rows"][0] == "1 1:100_A_G 100 A G 0 1 1 1"
    legend = export_vcf_to_phasing_format(vcf, tmp_path / "l", fmt="haplegendsample", verbose=False)
    assert legend["rows"] == ["0 1 1 1", "0 0 1 0"] and len(legend["files"]) == 3
    with pytest.raises(ConfigurationError):
        export_vcf_to_phasing_format(vcf, tmp_path / "g", fmt="gensample", verbose=False)
    with pytest.raises(ConfigurationError):
        export_vcf_to_phasing_format(vcf, tmp_path / "z", fmt="other", verbose=False)
    with pytest.raises(ConfigurationError):
        export_vcf_to_phasing_format(vcf, tmp_path / "missing_dir" / "z", verbose=False)


def test_hg_family_left_and_right_reach_upstream() -> None:
    pytest.importorskip("hgvs")
    from hgvs.utils import norm

    from dnakit.variants import normalize_indel_alleles

    assert "left" in inspect.signature(normalize_indel_alleles).parameters
    reference = "GGCACACACATT"
    # the CACACACA repeat covers bases 2 to 10; deleting one CA (bases 2 to 4) slides to 8 to 10
    for left in (True, False):
        wrapped = normalize_indel_alleles(reference, 2, 4, ["CA", ""], left=left, verbose=False)
        function = norm.normalize_alleles_left if left else norm.normalize_alleles_right
        bound = 0 if left else len(reference)
        direct = function(reference, 2, 4, ("CA", ""), bound, 1, True)
        assert (wrapped["start"], wrapped["stop"], wrapped["alleles"]) == (
            direct.start,
            direct.stop,
            list(direct.alleles),
        )
    left_result = normalize_indel_alleles(reference, 2, 4, ["CA", ""], left=True, verbose=False)
    right_result = normalize_indel_alleles(reference, 2, 4, ["CA", ""], left=False, verbose=False)
    assert (left_result["start"], right_result["start"]) == (2, 8)


def test_vf_family_substitution_branches_reach_upstream() -> None:
    pytest.importorskip("vcfpy")
    from vcfpy import parser

    from dnakit.variants import classify_vcf_substitution

    # process_sub sends a longer REF to process_sub_grow and a longer ALT to process_sub_shrink
    cases = [("CTT", "C", parser.process_sub_grow), ("C", "CTT", parser.process_sub_shrink)]
    for ref, alt, function in cases:
        direct = function(ref, alt)
        wrapped = classify_vcf_substitution(ref, alt, verbose=False)
        assert (wrapped["subtype"], wrapped["value"]) == (direct.type, direct.value)
    assert classify_vcf_substitution("CTT", "C", verbose=False)["subtype"] == "DEL"
    assert classify_vcf_substitution("C", "CTT", verbose=False)["subtype"] == "INS"


def test_dp_parsimony_alias_is_same_object() -> None:
    pytest.importorskip("dendropy")
    from dendropy.calculate import treescore
    from dendropy.model import parsimony

    assert parsimony.parsimony_score is treescore.parsimony_score


# -------------------------------------------------------------- predictions
def test_gr_family_encoding_rows_reach_upstream() -> None:
    pytest.importorskip("grelu")
    import numpy as np
    import torch
    from grelu.data import dataset
    from grelu.sequence import format as fmt
    from grelu.sequence import metrics, utils

    from dnakit.predictions import (
        check_dna_encoding,
        convert_dna_encoding,
        gc_content_distribution,
        summarize_sequence_lengths,
    )

    indices = np.array([[0, 1, 2, 3, 4], [3, 3, 0, 0, 1]], dtype=np.int8)
    one_hot = fmt.indices_to_one_hot(indices)
    wrapped = convert_dna_encoding(
        indices.tolist(), source="indices", target="one_hot", verbose=False
    )
    assert np.allclose(np.array(wrapped["values"]), one_hot.numpy())
    assert wrapped["values"][0][0] == [
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]  # A at position 0, N column is zeros
    back = convert_dna_encoding(one_hot.tolist(), source="one_hot", target="indices", verbose=False)
    assert back["values"] == fmt.one_hot_to_indices(one_hot).tolist() == indices.tolist()
    strings = convert_dna_encoding(
        indices.tolist(), source="indices", target="strings", verbose=False
    )
    assert strings["values"] == fmt.indices_to_strings(indices) == ["ACGTN", "TTAAC"]
    assert check_dna_encoding(indices.tolist(), encoding="indices", verbose=False)["valid"] is bool(
        fmt.check_indices(indices)
    )
    assert check_dna_encoding(one_hot.tolist(), encoding="one_hot", verbose=False)["valid"] is bool(
        fmt.check_one_hot(torch.tensor(one_hot.tolist()))
    )
    assert check_dna_encoding([[0, 9]], encoding="indices", verbose=False)["valid"] is False
    text = ["ACGT", "GGCC", "AATT"]
    gc = gc_content_distribution(text, mode="fraction", verbose=False)["gc"]
    assert gc == pytest.approx([0.5, 1.0, 0.0])
    assert gc == pytest.approx(list(metrics.gc(text)))
    lengths = summarize_sequence_lengths(text, verbose=False)
    assert lengths["all_equal"] is bool(utils.check_equal_lengths(text)) is True
    assert lengths["unique_length"] == utils.get_unique_length(text) == 4
    assert dataset.get_unique_length is utils.get_unique_length
    assert dataset.indices_to_one_hot is fmt.indices_to_one_hot
    assert summarize_sequence_lengths(["AC", "ACG"], verbose=False)["unique_length"] is None
    parameters = inspect.signature(convert_dna_encoding).parameters
    assert {"source", "target"} <= set(parameters)
    assert "encoding" in inspect.signature(check_dna_encoding).parameters
    assert "mode" in inspect.signature(gc_content_distribution).parameters


_AG_GTF = (
    "##small annotation\n"
    'chr1\tHAVANA\tgene\t1001\t2000\t.\t+\t.\tgene_id "ENSG1.1"; gene_type "protein_coding"; '
    'gene_name "GENEA";\n'
    'chr1\tHAVANA\ttranscript\t1001\t1800\t.\t+\t.\tgene_id "ENSG1.1"; transcript_id "ENST1.1"; '
    'gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; '
    'transcript_support_level "1"; tag "basic MANE_Select";\n'
    'chr1\tHAVANA\ttranscript\t1101\t2000\t.\t+\t.\tgene_id "ENSG1.1"; transcript_id "ENST2.1"; '
    'gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; '
    'transcript_support_level "3"; tag "basic";\n'
    'chr1\tHAVANA\tgene\t5001\t5500\t.\t-\t.\tgene_id "ENSG2.1"; gene_type "lncRNA"; '
    'gene_name "GENEB";\n'
    'chr1\tHAVANA\ttranscript\t5001\t5500\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST3.1"; '
    'gene_type "lncRNA"; gene_name "GENEB"; transcript_type "lncRNA"; '
    'transcript_support_level "NA"; tag "basic";\n'
)


def test_ag_family_filters_gene_interval_and_union_reach_upstream(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.data import gene_annotation as annotation
    from alphagenome.data import genome

    from dnakit.predictions import (
        combine_alphagenome_intervals,
        filter_alphagenome_transcripts,
        read_gencode_gene_intervals,
    )
    from dnakit.predictions.alphagenome_functions import _annotation_frame

    gtf = tmp_path / "small.gtf"
    gtf.write_text(_AG_GTF)
    assert {"method", "transcript_types", "support_levels"} <= set(
        inspect.signature(filter_alphagenome_transcripts).parameters
    )
    frame = _annotation_frame(gtf, 1000)

    def direct_ids(result: Any) -> list[str]:
        return sorted(str(item) for item in result["transcript_id"].dropna())

    longest = annotation.filter_to_longest_transcript(frame)
    assert filter_alphagenome_transcripts(gtf, method="longest", verbose=False)[
        "transcript_ids"
    ] == direct_ids(longest)
    mane = annotation.filter_to_mane_select_transcript(frame)
    assert (
        filter_alphagenome_transcripts(gtf, method="mane_select", verbose=False)["transcript_ids"]
        == direct_ids(mane)
        == ["ENST1.1"]
    )
    support = annotation.filter_transcript_support_level(frame, ["1"])
    assert filter_alphagenome_transcripts(
        gtf, method="support_level", support_levels=["1"], verbose=False
    )["transcript_ids"] == direct_ids(support)
    typed = annotation.filter_transcript_type(frame, (annotation.TranscriptType.LNCRNA,))
    assert (
        filter_alphagenome_transcripts(
            gtf, method="transcript_type", transcript_types=["lncRNA"], verbose=False
        )["transcript_ids"]
        == direct_ids(typed)
        == ["ENST3.1"]
    )

    single = annotation.get_gene_interval(frame, gene_symbol="GENEB")
    wrapped = read_gencode_gene_intervals(gtf, gene_symbols=["GENEB"], verbose=False)["genes"][0]
    assert (wrapped["chromosome"], wrapped["start"], wrapped["end"]) == (
        single.chromosome,
        single.start,
        single.end,
    )
    by_id = annotation.get_gene_interval(frame, gene_id="ENSG1")
    wrapped_id = read_gencode_gene_intervals(gtf, gene_ids=["ENSG1"], verbose=False)["genes"][0]
    assert (wrapped_id["start"], wrapped_id["end"]) == (by_id.start, by_id.end)

    left = [genome.Interval.from_str(text) for text in ("chr1:0-10", "chr1:20-30")]
    right = [genome.Interval.from_str("chr1:5-25")]
    upstream = [(item.start, item.end) for item in genome.union_intervals(left, right)]
    union = combine_alphagenome_intervals(
        ["chr1:0-10", "chr1:20-30"], ["chr1:5-25"], operation="union", verbose=False
    )
    assert [(item["start"], item["end"]) for item in union["intervals"]] == upstream == [(0, 30)]


def test_ag_family_fold_names_reach_upstream(tmp_path: Path) -> None:
    pytest.importorskip("alphagenome")
    from alphagenome.data import fold_intervals
    from alphagenome.models.v1 import dna_model

    from dnakit.predictions import get_alphagenome_fold_intervals

    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t0\t1000\ttest\nchr1\t1000\t2000\tfold0\n")
    assert {"model_version", "subset"} <= set(
        inspect.signature(get_alphagenome_fold_intervals).parameters
    )
    version = dna_model.ModelVersion.FOLD_0
    for subset in ("train", "valid", "test"):
        result = get_alphagenome_fold_intervals(
            regions, model_version="fold_0", subset=subset, verbose=False
        )
        direct = fold_intervals.get_fold_names(version, fold_intervals.Subset[subset.upper()])
        assert result["fold_names"] == sorted(direct)


def test_ef_family_one_hot_reverse_complement_reaches_upstream() -> None:
    pytest.importorskip("enformer_pytorch")
    from enformer_pytorch import data

    from dnakit.predictions import reverse_complement_enformer

    assert "representation" in inspect.signature(reverse_complement_enformer).parameters
    result = reverse_complement_enformer("ACGTN", representation="one_hot", verbose=False)
    direct = data.one_hot_reverse_complement(data.str_to_one_hot("ACGTN"))
    assert result["reverse_complement"] == direct.tolist()
    # reverse complement of ACGTN is NACGT: N is all zeros first, then A, C, G, T
    assert result["reverse_complement"][0] == [0.0, 0.0, 0.0, 0.0]
    assert result["reverse_complement"][1] == [1.0, 0.0, 0.0, 0.0]
