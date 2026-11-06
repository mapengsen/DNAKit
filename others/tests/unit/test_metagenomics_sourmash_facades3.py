"""Tests for dnakit.metagenomics_sourmash_facades3 (sourmash plot data and needs-data operations).

All inputs are built locally from the shared consistency data; nothing is downloaded.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

sourmash = pytest.importorskip("sourmash")
numpy = pytest.importorskip("numpy")

from dnakit import metagenomics_sourmash_facades3 as facades

RANDOM_FASTA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data/random.fa"
)


def _sequence(length: int = 3000) -> str:
    if not RANDOM_FASTA.is_file():
        pytest.skip("shared consistency data is not available")
    lines = RANDOM_FASTA.read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))[:length]


@pytest.fixture(scope="module")
def signatures(tmp_path_factory):
    from sourmash import MinHash
    from sourmash.signature import SourmashSignature, save_signatures_to_json

    seq = _sequence()

    def make(text: str, name: str) -> SourmashSignature:
        minhash = MinHash(n=0, ksize=21, scaled=1)
        minhash.add_sequence(text, True)
        return SourmashSignature(minhash, name=name)

    root = tmp_path_factory.mktemp("sourmash_inputs")
    sig_a = make(seq[:2000], "a")
    sig_b = make(seq[1000:3000], "b")
    sig_file = root / "query.sig"
    with sig_file.open("wt", encoding="utf-8") as handle:
        save_signatures_to_json([sig_a], handle)
    return {"root": root, "a": sig_a, "b": sig_b, "query": sig_file}


@pytest.fixture(scope="module")
def lca_databases(signatures):
    from sourmash.lca.lca_db import LCA_Database
    from sourmash.tax.tax_utils import LineagePair

    root = signatures["root"]
    db = LCA_Database(ksize=21, scaled=1, moltype="DNA")
    db.insert(
        signatures["a"],
        ident="gA",
        lineage=(
            LineagePair("superkingdom", "Bacteria", None),
            LineagePair("phylum", "Firmicutes", None),
        ),
    )
    db.insert(
        signatures["b"],
        ident="gB",
        lineage=(
            LineagePair("superkingdom", "Bacteria", None),
            LineagePair("phylum", "Bacteroidota", None),
        ),
    )
    path = root / "test.lca.json.gz"
    db.save(str(path))
    return {"path": str(path), "db": db}


@pytest.fixture(scope="module")
def taxonomy(signatures):
    path = signatures["root"] / "tax.csv"
    path.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species,strain\n"
        "gA,Bacteria,Firmicutes,,,,,,\n"
        "gB,Bacteria,Bacteroidota,,,,,,\n",
        encoding="utf-8",
    )
    return str(path)


@pytest.fixture(scope="module")
def sbt_index(signatures):
    from sourmash.sbtmh import create_sbt_index

    tree = create_sbt_index()
    tree.insert(signatures["a"])
    tree.insert(signatures["b"])
    tree.save(str(signatures["root"] / "index"))
    return str(signatures["root"] / "index.sbt.zip")


@pytest.fixture(scope="module")
def sqlite_index(signatures):
    from sourmash.index.sqlite_index import SqliteIndex

    path = signatures["root"] / "index.sqlite"
    index = SqliteIndex.create(str(path))
    index.insert(signatures["a"])
    index.insert(signatures["b"])
    index.close()
    return str(path)


@pytest.fixture(scope="module")
def comparison_base(tmp_path_factory):
    root = tmp_path_factory.mktemp("compare")
    labels = ["genome_one", "genome_two", "genome_three", "genome_four"]
    rng = numpy.random.default_rng(7)
    raw = rng.uniform(0.05, 0.95, size=(4, 4))
    matrix = (raw + raw.T) / 2
    numpy.fill_diagonal(matrix, 1.0)
    base = root / "compare_out"
    with base.open("wb") as handle:
        numpy.save(handle, matrix)
    (root / "compare_out.labels.txt").write_text("\n".join(labels) + "\n", encoding="utf-8")
    return str(base)


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert "verbose" in parameters, name
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_load_lca_database_describes_the_database(lca_databases):
    result = facades.load_lca_database(lca_databases["path"], verbose=False)
    assert result["status"] == "completed"
    assert result["ksize"] == 21
    assert result["scaled"] == 1
    assert result["signature_count"] == 2
    assert result["hash_count"] > 0


def test_load_lca_databases_returns_one_row_per_file(lca_databases):
    result = facades.load_lca_databases(
        [lca_databases["path"], lca_databases["path"]], verbose=False
    )
    assert len(result["databases"]) == 2
    assert result["databases"][0]["path"] == lca_databases["path"]


def test_classify_signature_lineage_finds_the_shared_lineage(signatures, lca_databases):
    result = facades.classify_signature_lineage(
        str(signatures["query"]), [lca_databases["path"]], threshold=1, verbose=False
    )
    [row] = result["classifications"]
    assert row["status"] == "found"
    assert row["lineage"] == [
        {"rank": "superkingdom", "name": "Bacteria"},
        {"rank": "phylum", "name": "Firmicutes"},
    ]


def test_classify_signature_lineage_majority_and_threshold(signatures, lca_databases):
    majority = facades.classify_signature_lineage(
        str(signatures["query"]), [lca_databases["path"]], threshold=1, majority=True, verbose=False
    )
    assert majority["majority"] is True
    assert majority["classifications"][0]["status"] == "found"
    high = facades.classify_signature_lineage(
        str(signatures["query"]), [lca_databases["path"]], threshold=10**9, verbose=False
    )
    assert high["classifications"][0]["status"] == "nomatch"
    assert high["classifications"][0]["lineage"] == []


def test_lca_rank_counts_counts_every_lineage(lca_databases):
    result = facades.lca_rank_counts([lca_databases["path"]], verbose=False)
    assert result["counts"]
    total = sum(row["count"] for row in result["counts"])
    assert total > 0
    names = {item["name"] for row in result["counts"] for item in row["lineage"]}
    assert {"Bacteria", "Firmicutes", "Bacteroidota"} <= names


def test_summarize_lineage_counts_aggregates_ancestors(signatures, lca_databases):
    result = facades.summarize_lineage_counts(
        str(signatures["query"]), [lca_databases["path"]], threshold=1, verbose=False
    )
    by_lineage = {
        tuple(item["name"] for item in row["lineage"]): row["count"] for row in result["counts"]
    }
    assert by_lineage[("Bacteria",)] >= by_lineage[("Bacteria", "Firmicutes")]
    assert by_lineage[("Bacteria", "Firmicutes")] > 0


def test_gather_lineage_assignments_lists_lineages_per_hash(signatures, lca_databases):
    result = facades.gather_lineage_assignments(
        str(signatures["query"]), [lca_databases["path"]], verbose=False
    )
    assert result["assigned_hash_count"] == len(result["assignments"]) > 0
    first = result["assignments"][0]
    assert first["hashval"].isdigit()
    assert all(isinstance(lineage, list) for lineage in first["lineages"])


def test_gather_lineage_assignments_requires_one_signature(signatures, lca_databases, tmp_path):
    from sourmash.signature import save_signatures_to_json

    two = tmp_path / "two.sig"
    with two.open("wt", encoding="utf-8") as handle:
        save_signatures_to_json([signatures["a"], signatures["b"]], handle)
    with pytest.raises(ConfigurationError):
        facades.gather_lineage_assignments(str(two), [lca_databases["path"]], verbose=False)


def test_load_sbt_index_summary_lists_signatures(sbt_index):
    result = facades.load_sbt_index_summary(sbt_index, verbose=False)
    assert result["index_type"] == "SBT"
    assert result["signature_count"] == 2
    assert sorted(result["signature_names"]) == ["a", "b"]
    assert result["ksizes"] == [21]


def test_load_sqlite_index_summary_lists_signatures(sqlite_index):
    result = facades.load_sqlite_index_summary(sqlite_index, verbose=False)
    assert result["index_type"] == "SqliteIndex"
    assert result["signature_count"] == 2
    assert sorted(result["signature_names"]) == ["a", "b"]


def test_load_comparison_matrix_returns_labels_and_matrix(comparison_base):
    result = facades.load_comparison_matrix(comparison_base, verbose=False)
    assert result["size"] == 4
    assert len(result["matrix"]) == 4 and len(result["matrix"][0]) == 4
    assert result["labels"] == ["genome_one", "genome_two", "genome_three", "genome_four"]


def test_plot_comparison_matrix_writes_image_and_data(comparison_base, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    pytest.importorskip("scipy")
    monkeypatch.setenv("MPLBACKEND", "Agg")
    out = tmp_path / "composite.png"
    result = facades.plot_comparison_matrix(comparison_base, str(out), verbose=False)
    assert Path(result["image"]).is_file()
    assert Path(result["image"]).stat().st_size > 0
    assert sorted(result["labels"]) == ["genome_four", "genome_one", "genome_three", "genome_two"]
    assert len(result["matrix"]) == 4


def test_plot_comparison_matrix_rejects_bad_inputs(comparison_base, tmp_path):
    with pytest.raises(ConfigurationError):
        facades.plot_comparison_matrix(comparison_base, str(tmp_path / "x.gif"), verbose=False)
    with pytest.raises(ConfigurationError):
        facades.plot_comparison_matrix(
            comparison_base, str(tmp_path / "x.png"), vmin=1.0, vmax=0.0, verbose=False
        )
    existing = tmp_path / "exists.png"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        facades.plot_comparison_matrix(comparison_base, str(existing), verbose=False)


def test_load_functions_reject_missing_files(tmp_path):
    missing = str(tmp_path / "missing.lca.json.gz")
    with pytest.raises(ConfigurationError):
        facades.load_lca_database(missing, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.load_sbt_index_summary(missing, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.load_sqlite_index_summary(missing, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.load_lca_databases([], verbose=False)
    with pytest.raises(ConfigurationError):
        facades.load_lca_databases([missing], scaled=0, verbose=False)


def test_missing_sourmash_gives_install_hint(monkeypatch, lca_databases):
    import importlib

    real = importlib.import_module

    def blocked(name, package=None):
        if name.startswith("sourmash"):
            raise ImportError("blocked for test")
        return real(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", blocked)
    with pytest.raises(BackendUnavailableError) as info:
        facades.load_lca_database(lca_databases["path"], verbose=False)
    assert "bio-sourmash" in (info.value.hint or "")


def test_gather_summaries_read_gather_csv(tmp_path, taxonomy):
    gather = tmp_path / "gather.csv"
    gather.write_text(
        "query_name,query_md5,query_filename,name,f_unique_weighted,f_unique_to_query,"
        "unique_intersect_bp,remaining_bp,query_bp,ksize,scaled\n"
        "q1,abc,q1.fa,gA,0.5,0.5,1000,0,2000,21,1\n"
        "q1,abc,q1.fa,gB,0.25,0.25,500,500,2000,21,1\n",
        encoding="utf-8",
    )
    loaded = facades.load_gather_results_summary(str(gather), [taxonomy], verbose=False)
    [query] = loaded["queries"]
    assert query["query_name"] == "q1"
    assert query["match_count"] == 2
    assert "f_unique_weighted" in loaded["columns"]

    checked = facades.check_and_load_gather_csvs_summary([str(gather)], [taxonomy], verbose=False)
    assert [q["query_name"] for q in checked["queries"]] == ["q1"]
