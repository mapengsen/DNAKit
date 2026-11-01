"""Tests for dnakit.annotation.gtars_p24_refget_functions (local refget stores built in tmp_path; no
network)."""

from __future__ import annotations

import inspect
import random

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("gtars")
refget = pytest.importorskip("gtars.refget")

from dnakit.annotation import gtars_p24_refget_functions as m  # noqa: E402


def _fasta(path, seed, names=("chrA", "chrB"), length=300):
    random.seed(seed)
    text = "".join(
        f">{n}\n{''.join(random.choice('ACGT') for _ in range(length))}\n" for n in names
    )
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def env(tmp_path_factory):
    root = tmp_path_factory.mktemp("p24_refget")
    one = _fasta(root / "one.fa", 1)
    two = _fasta(root / "two.fa", 2, names=("chrA", "chrC"))
    store = root / "store"
    first = m.refget_store_add_fasta(store, [one], verbose=False)
    second = m.refget_store_add_fasta(store, [two], verbose=False)
    return {
        "root": root,
        "one": one,
        "two": two,
        "store": store,
        "d1": first["collections"][0]["digest"],
        "d2": second["collections"][0]["digest"],
    }


def test_public_functions_are_keyword_only_verbose():
    for name in m.__all__:
        parameters = inspect.signature(getattr(m, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_add_fasta_matches_gtars_digest_and_is_idempotent(env):
    expected = str(refget.digest_fasta(str(env["one"])).digest)
    assert env["d1"] == expected
    again = m.refget_store_add_fasta(env["store"], [env["one"]], verbose=False)
    assert again["collections"][0]["added"] is False
    assert again["collection_count"] == 2


def test_add_fasta_many_and_via_collection_object(env, tmp_path):
    result = m.refget_store_add_fasta(tmp_path / "s", [env["one"], env["two"]], verbose=False)
    assert [row["digest"] for row in result["collections"]] == [env["d1"], env["d2"]]
    other = m.refget_store_add_fasta(
        tmp_path / "t", [env["one"]], use_collection_object=True, verbose=False
    )
    assert other["collections"][0]["digest"] == env["d1"]


def test_add_sequence_digest_and_edge_cases(tmp_path):
    result = m.refget_store_add_sequence(tmp_path / "s", "acgtnn", name="x", verbose=False)
    assert result["sha512t24u"] == refget.sha512t24u_digest("ACGTNN")
    assert result["length"] == 6
    fetched = m.refget_sequence_retrieve(
        tmp_path / "s", digest=result["sha512t24u"], include_sequence=True, verbose=False
    )
    assert fetched["sequence"] == "ACGTNN"
    with pytest.raises(ConfigurationError):
        m.refget_store_add_sequence(tmp_path / "s", "", verbose=False)


def test_collection_levels_and_retrieval(env):
    level1 = m.refget_collection_digests(env["store"], env["d1"], verbose=False)["representation"]
    level2 = m.refget_collection_digests(env["store"], env["d1"], level=2, verbose=False)[
        "representation"
    ]
    assert level2["names"] == ["chrA", "chrB"] and level2["lengths"] == [300, 300]
    assert set(level1) >= {"names", "lengths", "sequences"}
    meta = m.refget_collection_retrieve(env["store"], digest=env["d1"], verbose=False)
    full = m.refget_collection_retrieve(
        env["store"], digest=env["d1"], include_sequences=True, verbose=False
    )
    assert meta["total_bases"] == 600 and "sequence" not in meta["sequences"][0]
    assert len(full["sequences"][0]["sequence"]) == 300
    with pytest.raises(ConfigurationError):
        m.refget_collection_retrieve(env["store"], digest="nope", verbose=False)
    with pytest.raises(ConfigurationError):
        m.refget_collection_retrieve(
            env["store"], include_sequences=True, max_bases=10, digest=env["d1"], verbose=False
        )


def test_list_compare_and_attribute_search(env):
    listing = m.refget_collections_list(env["store"], verbose=False)
    assert {row["digest"] for row in listing["collections"]} == {env["d1"], env["d2"]}
    result = m.refget_collection_compare(env["store"], env["d1"], env["d2"], verbose=False)[
        "comparison"
    ]
    assert "names" in result["attributes"]["a_and_b"] or "names" in result["attributes"]["a_only"]
    names = m.refget_collection_digests(env["store"], env["d1"], verbose=False)["representation"][
        "names"
    ]
    assert m.refget_collections_by_attribute(env["store"], "names", names, verbose=False)[
        "collections"
    ] == [env["d1"]]
    assert m.refget_collection_attribute_array(env["store"], "names", names, verbose=False)[
        "values"
    ] == ["chrA", "chrB"]
    filtered = m.refget_collections_list(
        env["store"], attribute="names", attribute_digest=names, verbose=False
    )
    assert [row["digest"] for row in filtered["collections"]] == [env["d1"]]
    with pytest.raises(ConfigurationError):
        m.refget_collections_by_attribute(env["store"], "bogus", names, verbose=False)


def test_aliases_register_list_lookup_remove(env, tmp_path):
    store = tmp_path / "alias_store"
    added = m.refget_store_add_fasta(store, [env["one"]], verbose=False)
    digest = added["collections"][0]["digest"]
    seq = m.refget_collection_digests(store, digest, level=2, verbose=False)["representation"][
        "sequences"
    ][0]
    m.refget_alias_register(store, "collection", "demo", alias="v1", digest=digest, verbose=False)
    m.refget_alias_register(store, "sequence", "demo", alias="chrA", digest=seq, verbose=False)
    table = tmp_path / "aliases.tsv"
    table.write_text(f"second\t{seq[3:]}\n", encoding="utf-8")
    assert (
        m.refget_alias_register(store, "sequence", "file_ns", alias_file=table, verbose=False)[
            "aliases_registered"
        ]
        == 1
    )
    assert m.refget_alias_list(store, "sequence", verbose=False)["namespaces"] == [
        "demo",
        "file_ns",
    ]
    assert m.refget_alias_list(store, "collection", namespace="demo", verbose=False)["aliases"] == [
        "v1"
    ]
    assert (
        m.refget_alias_list(store, "collection", namespace="zzz", verbose=False)["aliases"] is None
    )
    back = m.refget_alias_reverse_lookup(store, "sequence", seq, verbose=False)["aliases"]
    assert back == [
        {"namespace": "demo", "alias": "chrA"},
        {"namespace": "file_ns", "alias": "second"},
    ]
    by_alias = m.refget_collection_retrieve(store, namespace="demo", alias="v1", verbose=False)
    assert by_alias["digest"] == digest
    sequence = m.refget_sequence_retrieve(
        store, namespace="demo", alias="chrA", include_sequence=True, verbose=False
    )
    assert sequence["name"] == "chrA" and len(sequence["sequence"]) == 300
    assert (
        m.refget_alias_remove(store, "sequence", "demo", "chrA", verbose=False)["removed"] is True
    )
    assert (
        m.refget_alias_remove(store, "sequence", "demo", "chrA", verbose=False)["removed"] is False
    )
    with pytest.raises(ConfigurationError):
        m.refget_alias_register(store, "sequence", "demo", verbose=False)
    with pytest.raises(ConfigurationError):
        m.refget_collection_retrieve(store, namespace="demo", alias="missing", verbose=False)


def test_sequence_retrieve_by_md5_and_substrings(env):
    level2 = m.refget_collection_digests(env["store"], env["d1"], level=2, verbose=False)[
        "representation"
    ]
    seq = m.refget_sequence_retrieve(
        env["store"], digest=level2["sequences"][0], include_sequence=True, verbose=False
    )
    by_md5 = m.refget_sequence_retrieve(env["store"], digest=seq["md5"], verbose=False)
    assert by_md5["sha512t24u"] == seq["sha512t24u"]
    pieces = m.refget_sequence_substring(
        env["store"],
        [[0, 10], [295, 300], [5, 5]],
        collection=env["d1"],
        name="chrA",
        verbose=False,
    )
    assert [p["sequence"] for p in pieces["substrings"]] == [
        seq["sequence"][0:10],
        seq["sequence"][295:300],
        "",
    ]
    by_digest = m.refget_sequence_substring(
        env["store"], [[2, 4]], digest=level2["sequences"][0], verbose=False
    )
    assert by_digest["substrings"][0]["sequence"] == seq["sequence"][2:4]
    with pytest.raises(ConfigurationError):
        m.refget_sequence_substring(
            env["store"], [[4, 2]], digest=level2["sequences"][0], verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.refget_sequence_substring(
            env["store"], [[0, 2]], collection=env["d1"], name="nope", verbose=False
        )


def test_import_plan_and_remove(env, tmp_path):
    destination = tmp_path / "dest"
    imported = m.refget_collection_import(destination, env["store"], env["d1"], verbose=False)
    assert imported["imported"]["digest"] == env["d1"]
    plan = m.refget_orphan_removal_plan(destination, env["d1"], verbose=False)
    assert len(plan["orphan_sequences"]) == 2
    removed = m.refget_collection_remove(
        destination, env["d1"], remove_orphan_sequences=True, verbose=False
    )
    assert removed["removed"] is True
    assert m.refget_collection_remove(destination, env["d1"], verbose=False)["removed"] is False
    assert m.refget_collections_list(destination, verbose=False)["collections"] == []


def test_write_rgsi(env, tmp_path):
    out = tmp_path / "x.rgsi"
    result = m.refget_write_rgsi(env["one"], out, verbose=False)
    text = out.read_text(encoding="utf-8")
    assert text.startswith("##seqcol_digest=" + env["d1"])
    assert result["sequence_count"] == 2
    with pytest.raises(ConfigurationError):
        m.refget_write_rgsi(env["one"], out, verbose=False)


def test_input_checks(env, tmp_path):
    with pytest.raises(ConfigurationError):
        m.refget_collections_list(tmp_path / "missing_store", verbose=False)
    with pytest.raises(ConfigurationError):
        m.refget_store_add_fasta(tmp_path / "s", [], verbose=False)
    with pytest.raises(ConfigurationError):
        m.refget_alias_list(env["store"], "bogus", verbose=False)
    with pytest.raises(ConfigurationError):
        m.refget_alias_pull("ftp://x", tmp_path / "c", verbose=False)
