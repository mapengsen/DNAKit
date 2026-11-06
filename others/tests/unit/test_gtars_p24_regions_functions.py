"""Tests for dnakit.annotation.gtars_p24_regions_functions (region sets, annotation, tokenizers,
LOLA databases)."""

from __future__ import annotations

import inspect

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("gtars")
models = pytest.importorskip("gtars.models")
pytest.importorskip("gtars.tokenizers")
pytest.importorskip("gtars.lola")

from dnakit.annotation import gtars_p24_regions_functions as m

A_ROWS = [("chr1", 100, 200), ("chr1", 150, 300), ("chr1", 500, 600), ("chr2", 10, 50)]
B_ROWS = [("chr1", 180, 520), ("chr2", 40, 90), ("chr3", 5, 15)]
SIZES = {"chr1": 1000, "chr2": 100}


def _bed(path, rows):
    path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows), encoding="utf-8")
    return path


@pytest.fixture()
def beds(tmp_path):
    return {
        "a": _bed(tmp_path / "a.bed", A_ROWS),
        "b": _bed(tmp_path / "b.bed", B_ROWS),
        "dir": tmp_path,
    }


def _tuples(result):
    return [(r["chrom"], r["start"], r["end"]) for r in result["regions"]]


def test_public_functions_are_keyword_only_verbose():
    for name in m.__all__:
        parameters = inspect.signature(getattr(m, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_summary_matches_hand_computation(beds):
    result = m.region_set_summary(beds["a"], include_widths=True, verbose=False)
    assert result["widths"] == [100, 150, 100, 40]
    assert result["total_nucleotides"] == 390 and result["region_count"] == 4
    assert result["max_end_per_chromosome"] == {"chr1": 600, "chr2": 50}
    assert result["mean_region_width"] == pytest.approx(97.5)
    assert result["chromosomes"][0]["number_of_regions"] == 3


def test_neighbor_distances_and_bins(beds):
    gaps = m.region_neighbor_distances(beds["a"], verbose=False)["distances"]
    assert gaps == [200]  # gtars merges the overlapping first two regions before measuring
    nearest = m.region_neighbor_distances(beds["a"], kind="nearest", verbose=False)["distances"]
    assert nearest[0] == 0 and nearest[2] == 200
    bins = m.region_distribution_bins(beds["a"], n_bins=2, chrom_sizes=SIZES, verbose=False)
    assert sum(row["count"] for row in bins["bins"]) == 4
    assert {row["chrom"] for row in bins["bins"]} == {"chr1", "chr2"}
    with pytest.raises(ConfigurationError):
        m.region_neighbor_distances(beds["a"], kind="other", verbose=False)


def test_promoters_and_output_file(beds):
    out = beds["dir"] / "prom.bed"
    result = m.region_promoters(beds["a"], 50, 10, output=out, verbose=False)
    assert _tuples(result)[0] == ("chr1", 50, 110)
    assert out.read_text(encoding="utf-8").splitlines()[0].startswith("chr1\t50\t110")
    with pytest.raises(ConfigurationError):
        m.region_promoters(beds["a"], 50, 10, output=out, verbose=False)


def test_pairwise_operations_overlap_and_subset(beds):
    all_vs_all = m.region_pairwise_operation(
        beds["a"], beds["b"], operation="intersect_all", verbose=False
    )
    assert _tuples(all_vs_all) == [("chr1", 180, 300), ("chr1", 500, 520), ("chr2", 40, 50)]
    joined = m.region_pairwise_operation(beds["a"], beds["b"], operation="concat", verbose=False)
    assert joined["region_count"] == 7
    paired = m.region_pairwise_operation(
        beds["a"], beds["b"], operation="pintersect", verbose=False
    )
    assert paired["region_count"] == 3
    kept = m.region_subset_by_overlaps(beds["a"], beds["b"], verbose=False)
    assert _tuples(kept) == [
        ("chr1", 100, 200),
        ("chr1", 150, 300),
        ("chr1", 500, 600),
        ("chr2", 10, 50),
    ]
    found = m.region_overlap_indices(beds["a"], beds["b"], verbose=False)
    assert found["overlaps"] == [[0], [0], [0], [1]]
    assert found["overlapping_region_count"] == 4
    value = m.region_overlap_coefficient(beds["a"], beds["b"], verbose=False)["overlap_coefficient"]
    assert 0.0 <= value <= 1.0
    with pytest.raises(ConfigurationError):
        m.region_pairwise_operation(beds["a"], beds["b"], operation="union", verbose=False)


def test_gaps_and_disjoin(beds):
    gaps = m.region_gaps(beds["a"], SIZES, verbose=False)
    assert ("chr1", 0, 100) in _tuples(gaps) and ("chr1", 600, 1000) in _tuples(gaps)
    assert ("chr2", 0, 10) in _tuples(gaps) and ("chr2", 50, 100) in _tuples(gaps)
    pieces = m.region_disjoin(beds["a"], verbose=False)
    assert _tuples(pieces)[:3] == [("chr1", 100, 150), ("chr1", 150, 200), ("chr1", 200, 300)]
    with pytest.raises(ConfigurationError):
        m.region_gaps(beds["a"], {}, verbose=False)


def _gtf(path):
    lines = []
    for i in range(3):
        s = 1000 + i * 5000
        e = s + 3000
        strand = "+" if i % 2 == 0 else "-"
        attrs = (
            f'gene_id "G{i}"; transcript_id "T{i}"; gene_biotype "protein_coding"; '
            f'gene_name "N{i}";'
        )
        lines.append(f"chr1\tdemo\tgene\t{s + 1}\t{e}\t.\t{strand}\t.\t{attrs}")
        lines.append(f"chr1\tdemo\ttranscript\t{s + 1}\t{e}\t.\t{strand}\t.\t{attrs}")
        lines.append(f"chr1\tdemo\texon\t{s + 1}\t{s + 500}\t.\t{strand}\t.\t{attrs}")
        lines.append(f"chr1\tdemo\texon\t{s + 1500}\t{e}\t.\t{strand}\t.\t{attrs}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_gene_annotation_and_tss_distances(beds):
    gtf = _gtf(beds["dir"] / "g.gtf")
    summary = m.gene_annotation_summary(gtf=gtf, verbose=False)
    assert summary["n_genes"] == 3 and summary["n_exons"] == 6 and summary["tss_count"] == 3
    assert summary["partition_names"] == ["promoterCore", "promoterProx", "exon", "intron"]
    query = _bed(beds["dir"] / "q.bed", [("chr1", 900, 1000), ("chr1", 20000, 20100)])
    unsigned = m.tss_distances(query, gtf=gtf, verbose=False)["distances"]
    signed = m.tss_distances(query, gtf=gtf, signed=True, verbose=False)["distances"]
    assert unsigned[0]["distance"] == pytest.approx(abs(signed[0]["distance"]))
    tss = _bed(beds["dir"] / "tss.bed", [("chr1", 1000, 1001)])
    from_bed = m.tss_distances(query, tss_regions=tss, verbose=False)
    assert (
        from_bed["tss_count"] == 1 and from_bed["distances"][0]["distance"] == 50
    )  # from the region midpoint (950) to the TSS (1000)
    with pytest.raises(ConfigurationError):
        m.gene_annotation_summary(verbose=False)
    with pytest.raises(ConfigurationError):
        m.tss_distances(query, verbose=False)


def test_signal_matrix_info(beds):
    table = beds["dir"] / "sm.tsv"
    table.write_text(
        "label\tc1\tc2\nchr1_100_500\t1.5\t2\nchr1_1000_2000\t0.5\t3\n", encoding="utf-8"
    )
    info = m.signal_matrix_info(table, verbose=False)
    assert info["condition_names"] == ["c1", "c2"] and info["n_regions"] == 2
    bad = beds["dir"] / "bad.tsv"
    bad.write_text("nothing\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.signal_matrix_info(bad, verbose=False)


def test_tokenizer_vocabulary_ids_and_mask(beds):
    universe = _bed(
        beds["dir"] / "u.bed", [("chr1", 0, 1000), ("chr1", 1000, 2000), ("chr2", 0, 500)]
    )
    vocab = m.region_tokenizer_vocabulary(universe, verbose=False)
    assert vocab["vocab_size"] == 10 and vocab["vocabulary"]["chr2:0-500"] == 2
    assert vocab["special_tokens"]["cls_token"] == "<cls>"
    ids = m.region_token_ids(universe, ["chr1:1000-2000", "unknown"], verbose=False)["outputs"]
    assert ids[0] == 1 and ids[1] == vocab["vocabulary"]["<unk>"]
    back = m.region_token_ids(
        universe, [0, vocab["vocabulary"]["<cls>"]], direction="to_tokens", verbose=False
    )
    assert back["outputs"] == ["chr1:0-1000", "<cls>"]
    mask = m.region_special_tokens_mask(universe, ["<cls>", "chr1:0-1000"], verbose=False)
    assert mask["mask"] == [True, False]
    config = beds["dir"] / "c.toml"
    config.write_text(f'universe = "{universe}"\ntokenizer_type = "ailist"\n', encoding="utf-8")
    from_config = m.region_tokenizer_vocabulary(
        config, source_type="config", include_vocabulary=False, verbose=False
    )
    assert from_config["vocab_size"] == 10 and "vocabulary" not in from_config
    with pytest.raises(ConfigurationError):
        m.region_token_ids(universe, [], verbose=False)
    with pytest.raises(ConfigurationError):
        m.region_tokenizer_vocabulary(universe, source_type="remote", verbose=False)


def test_region_database(beds):
    listing = m.region_database_list_sets(bed_files=[beds["a"], beds["b"]], verbose=False)
    assert [row["name"] for row in listing["sets"]] == ["a.bed", "b.bed"]
    picked = m.region_database_get_sets(
        bed_files=[beds["a"], beds["b"]], indices=[1], verbose=False
    )
    assert [(r["chrom"], r["start"], r["end"]) for r in picked["sets"][0]["regions"]] == B_ROWS
    everything = m.region_database_get_sets(bed_files=[beds["a"], beds["b"]], verbose=False)
    assert len(everything["sets"]) == 2
    with pytest.raises(ConfigurationError):
        m.region_database_get_sets(bed_files=[beds["a"]], indices=[3], verbose=False)
    with pytest.raises(ConfigurationError):
        m.region_database_list_sets(verbose=False)
