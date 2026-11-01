"""Tests for dnakit.metagenomics_verify_sourmash_functions.

Each wrapper is compared with the upstream sourmash command it replaces (sourmash 4.9.4) on the same small
synthetic genomes: the CSV files, the printed rank table, the compare-csv verdicts and the index report must
match. Tests skip when sourmash is not installed and never use the network.
"""

from __future__ import annotations

import csv
import inspect
import random
from argparse import Namespace
from pathlib import Path

import pytest

from dnakit import metagenomics_verify_sourmash_functions as mod
from dnakit.exceptions import ConfigurationError

PUBLIC = [
    "classify_lca_signatures",
    "compare_lineage_assignment_csvs",
    "lca_rank_info",
    "summarize_lca_signatures",
    "write_lca_index_report",
]

LINEAGE_HEADER = ["ident", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
LINEAGE_A = ["Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus", "Bacillus alpha", ""]
LINEAGE_B = ["Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides", "Bacteroides beta", ""]


def _random_dna(seed: int, length: int = 3000) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _write_signature(path: Path, name: str, sequence: str, scaled: int = 100) -> None:
    sourmash = pytest.importorskip("sourmash")
    from sourmash.signature import SourmashSignature, save_signatures_to_json

    minhash = sourmash.MinHash(n=0, ksize=21, scaled=scaled)
    minhash.add_sequence(sequence, False)
    path.write_bytes(save_signatures_to_json([SourmashSignature(minhash, name=name)]))


@pytest.fixture(scope="module")
def lca_fixture(tmp_path_factory):
    """Two reference genomes indexed as an LCA database, plus query signatures.

    ``genA`` is a reference with lineage LINEAGE_A, ``genB`` one with LINEAGE_B. ``qA`` is a copy of genA,
    ``qC`` an unrelated genome (no match in the database).
    """
    pytest.importorskip("sourmash.lca.command_index")
    root = tmp_path_factory.mktemp("sourmash_verify")
    seq_a = _random_dna(1)
    seq_b = _random_dna(2)
    seq_c = _random_dna(3)
    _write_signature(root / "genA.sig", "genA", seq_a)
    _write_signature(root / "genB.sig", "genB", seq_b)
    _write_signature(root / "qA.sig", "qA", seq_a)
    _write_signature(root / "qC.sig", "qC", seq_c)
    _write_signature(root / "qA_fine.sig", "qA_fine", seq_a, scaled=50)
    _write_signature(root / "qA_coarse.sig", "qA_coarse", seq_a, scaled=200)
    lineages = root / "lineages.csv"
    with lineages.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(LINEAGE_HEADER)
        writer.writerow(["genA"] + LINEAGE_A)
        writer.writerow(["genB"] + LINEAGE_B)
    db = root / "lca.db"
    import sourmash

    args = sourmash.cli.get_parser().parse_args(
        ["lca", "index", str(lineages), str(db), str(root / "genA.sig"), str(root / "genB.sig"),
         "-k", "21", "--scaled", "100", "--dna"]
    )
    sourmash.lca.command_index.index(args)
    # sourmash writes the database as <output>.lca.json
    return {
        "root": root, "db": root / "lca.db.lca.json", "qA": root / "qA.sig", "qC": root / "qC.sig",
        "qA_fine": root / "qA_fine.sig", "qA_coarse": root / "qA_coarse.sig",
    }


def _rows(path: Path) -> list[list[str]]:
    with path.open(newline="") as fp:
        return [row for row in csv.reader(fp)]


def _upstream_args(**values):
    base = {"quiet": True, "debug": False}
    base.update(values)
    return Namespace(**base)


def test_public_functions_are_exported_and_keyword_verbose():
    for name in PUBLIC:
        function = getattr(mod, name)
        assert name in mod.__all__
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_classify_matches_upstream_csv(lca_fixture, tmp_path, capsys):
    pytest.importorskip("sourmash.lca.command_classify")
    from sourmash.lca.command_classify import classify

    for majority in (False, True):
        upstream = tmp_path / f"upstream_{majority}.csv"
        classify(_upstream_args(
            db=[[str(lca_fixture["db"])]], query=[[str(lca_fixture["qA"]), str(lca_fixture["qC"])]],
            query_from_file=None, output=str(upstream), scaled=None, threshold=5, majority=majority,
        ))
        ours = tmp_path / f"ours_{majority}.csv"
        result = mod.classify_lca_signatures(
            [str(lca_fixture["qA"]), str(lca_fixture["qC"])], [str(lca_fixture["db"])], str(ours),
            threshold=5, majority=majority,
        )
        assert _rows(ours) == _rows(upstream)
        assert result["classification_count"] == 2
        statuses = {row["id"]: row["status"] for row in result["classifications"]}
        assert statuses == {"qA": "found", "qC": "nomatch"}


def test_classify_downsamples_finer_query_like_upstream(lca_fixture, tmp_path):
    pytest.importorskip("sourmash.lca.command_classify")
    from sourmash.lca.command_classify import classify

    upstream = tmp_path / "upstream_fine.csv"
    classify(_upstream_args(
        db=[[str(lca_fixture["db"])]], query=[[str(lca_fixture["qA_fine"])]],
        query_from_file=None, output=str(upstream), scaled=None, threshold=5, majority=False,
    ))
    ours = tmp_path / "ours_fine.csv"
    result = mod.classify_lca_signatures([str(lca_fixture["qA_fine"])], [str(lca_fixture["db"])], str(ours))
    assert _rows(ours) == _rows(upstream)
    assert result["classifications"][0]["status"] == "found"


def test_classify_rejects_query_coarser_than_database(lca_fixture, tmp_path):
    with pytest.raises(ConfigurationError):
        mod.classify_lca_signatures([str(lca_fixture["qA_coarse"])], [str(lca_fixture["db"])],
                                    str(tmp_path / "coarse.csv"))


def test_classify_refuses_existing_output_and_empty_inputs(lca_fixture, tmp_path):
    existing = tmp_path / "exists.csv"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        mod.classify_lca_signatures([str(lca_fixture["qA"])], [str(lca_fixture["db"])], str(existing))
    with pytest.raises(ConfigurationError):
        mod.classify_lca_signatures([], [str(lca_fixture["db"])], str(tmp_path / "a.csv"))
    with pytest.raises(ConfigurationError):
        mod.classify_lca_signatures([str(lca_fixture["qA"])], [str(lca_fixture["db"])],
                                    str(tmp_path / "b.csv"), threshold=-1)


def test_summarize_csv_matches_upstream(lca_fixture, tmp_path):
    pytest.importorskip("sourmash.lca.command_summarize")
    from sourmash.lca.command_summarize import summarize_main

    for ignore in (False, True):
        upstream = tmp_path / f"upstream_{ignore}.csv"
        summarize_main(_upstream_args(
            db=[[str(lca_fixture["db"])]], query=[[str(lca_fixture["qA"]), str(lca_fixture["qC"])]],
            query_from_file=None, output=str(upstream), scaled=None, threshold=5, ignore_abundance=ignore,
        ))
        ours = tmp_path / f"ours_{ignore}.csv"
        result = mod.summarize_lca_signatures(
            [str(lca_fixture["qA"]), str(lca_fixture["qC"])], [str(lca_fixture["db"])],
            output_csv=str(ours), threshold=5, ignore_abundance=ignore,
        )
        assert _rows(ours) == _rows(upstream)
        assert result["signature_count"] == 2
        assert result["output_csv"] == str(ours)


def test_summarize_counts_match_summarize_facade(lca_fixture):
    from dnakit.metagenomics_sourmash_facades3 import summarize_lineage_counts

    result = mod.summarize_lca_signatures([str(lca_fixture["qA"])], [str(lca_fixture["db"])])
    assert result["output_csv"] is None
    facade = summarize_lineage_counts(str(lca_fixture["qA"]), [str(lca_fixture["db"])], verbose=False)
    ours = {tuple((item["rank"], item["name"]) for item in row["lineage"]): row["count"]
            for row in result["signatures"][0]["counts"]}
    theirs = {tuple((item["rank"], item["name"]) for item in row["lineage"]): row["count"]
              for row in facade["counts"]}
    assert ours == theirs
    assert any(count > 0 for count in ours.values())


def test_rank_info_matches_upstream_printed_table(lca_fixture, capsys):
    pytest.importorskip("sourmash.lca.command_rankinfo")
    from sourmash.lca.command_rankinfo import rankinfo_main

    for min_num in (0, 3):
        capsys.readouterr()
        rankinfo_main(_upstream_args(db=[str(lca_fixture["db"])], scaled=None, minimum_num=min_num))
        upstream_lines = capsys.readouterr().out.strip().splitlines()
        result = mod.lca_rank_info([str(lca_fixture["db"])], min_num=min_num)
        ours_lines = [f"{row['rank']}: {row['count']} ({row['percent']:.1f}%)" for row in result["ranks"]]
        assert ours_lines == upstream_lines
        assert result["total"] == sum(row["count"] for row in result["ranks"]) or result["total"] == 0


def test_rank_info_rejects_negative_minimum(lca_fixture):
    with pytest.raises(ConfigurationError):
        mod.lca_rank_info([str(lca_fixture["db"])], min_num=-1)


def _custom_lineage_table(path: Path) -> None:
    rows = [
        ["ident", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"],
        ["s1", *LINEAGE_A],
        ["s2", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "", "", "", ""],
        ["s3", "Bacteria", "Firmicutes", "", "", "", "", "", ""],
        ["s5", *LINEAGE_B],
    ]
    with path.open("w", newline="") as fp:
        csv.writer(fp).writerows(rows)


def _classify_style_table(path: Path) -> None:
    rows = [
        ["ID", "status", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"],
        ["s1", "found", *LINEAGE_A],
        ["s2", "found", "Bacteria", "Firmicutes", "", "", "", "", "", ""],
        ["s3", "found", "Bacteria", "Bacteroidota", "", "", "", "", "", ""],
        ["s4", "found", "Bacteria", "Firmicutes", "Bacilli", "", "", "", "", ""],
    ]
    with path.open("w", newline="") as fp:
        csv.writer(fp).writerows(rows)


def test_compare_csv_matches_upstream_verdicts(tmp_path, capsys):
    pytest.importorskip("sourmash.lca.command_compare_csv")
    from sourmash.lca.command_compare_csv import compare_csv

    classify_csv = tmp_path / "classify.csv"
    custom_csv = tmp_path / "custom.csv"
    _classify_style_table(classify_csv)
    _custom_lineage_table(custom_csv)
    capsys.readouterr()
    compare_csv(_upstream_args(
        csv1=str(classify_csv), csv2=str(custom_csv), start_column=2, tabs=False, no_headers=False, force=False,
        quiet=False,
    ))
    upstream = {}
    for line in capsys.readouterr().out.strip().splitlines():
        identifier, verdict, _ = line.split(",", 2)
        upstream[identifier] = verdict
    result = mod.compare_lineage_assignment_csvs(str(classify_csv), str(custom_csv))
    ours = {row["identifier"]: row["verdict"] for row in result["differences"]}
    assert ours == upstream
    assert ours == {"s2": "compatible", "s3": "incompatible"}
    assert result["missing_in_custom"] == ["s4"]
    assert result["missing_in_classify"] == ["s5"]
    assert result["common_count"] == 3
    assert result["incompatible_count"] == 1
    by_rank = {row["rank"]: row["count"] for row in result["incompatible_by_rank"]}
    assert by_rank["superkingdom"] == 1 and sum(by_rank.values()) == 1


def test_compare_csv_rejects_start_column_below_two(tmp_path):
    table = tmp_path / "t.csv"
    _classify_style_table(table)
    with pytest.raises(ConfigurationError):
        mod.compare_lineage_assignment_csvs(str(table), str(table), start_column=1)


def test_write_report_matches_upstream_file(tmp_path):
    pytest.importorskip("sourmash.lca.command_index")
    from sourmash.lca.command_index import generate_report
    from sourmash.lca.lca_utils import LineagePair

    duplicates = ["dupA", "dupB"]
    no_lineage = ["idX"]
    no_signature = ["idY", "idZ"]
    unused_ids = ["idU"]
    unused_lineage_dicts = [
        [{"rank": "superkingdom", "name": "Archaea"}, {"rank": "phylum", "name": "Euryarchaeota"}],
        [],
    ]
    unused_pairs = [
        tuple(LineagePair(item["rank"], item["name"]) for item in lineage) for lineage in unused_lineage_dicts
    ]
    upstream = tmp_path / "upstream_report.txt"
    generate_report(duplicates, no_lineage, no_signature, unused_pairs, unused_ids, str(upstream))
    ours = tmp_path / "ours_report.txt"
    result = mod.write_lca_index_report(
        duplicates, no_lineage, no_signature, unused_ids, unused_lineage_dicts, str(ours)
    )
    assert ours.read_text() == upstream.read_text()
    assert result["unused_lineage_count"] == 2
    assert result["no_signature_count"] == 2


def test_write_report_refuses_existing_file(tmp_path):
    existing = tmp_path / "report.txt"
    existing.write_text("old")
    with pytest.raises(ConfigurationError):
        mod.write_lca_index_report([], [], [], [], [], str(existing))


def test_verbose_flag_does_not_change_results(lca_fixture, tmp_path):
    first = mod.lca_rank_info([str(lca_fixture["db"])], verbose=False)
    second = mod.lca_rank_info([str(lca_fixture["db"])], verbose=True)
    assert first == second
