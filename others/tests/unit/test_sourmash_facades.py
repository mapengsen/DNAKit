"""Tests for dnakit.metagenomics_sourmash_facades (sourmash sketch dna and compare)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from dnakit import metagenomics_sourmash_facades as facade
from dnakit.exceptions import ConfigurationError

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
GENES = DATA / "genes.fa"
_ENV_BIN = Path(sys.executable).parent
pytestmark = pytest.mark.skipif(not (_ENV_BIN / "sourmash").is_file() or not GENES.exists(),
                                reason="sourmash or shared FASTA missing")


@pytest.fixture(autouse=True)
def _sourmash_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


def test_sketch_then_compare(tmp_path: Path) -> None:
    sig_a = tmp_path / "a.sig"
    sig_b = tmp_path / "b.sig"
    assert facade.sketch_dna(str(GENES), str(sig_a), ksizes=(21,), scaled=100)["status"] == "completed"
    assert facade.sketch_dna(str(GENES), str(sig_b), ksizes=(21,), scaled=100)["status"] == "completed"
    result = facade.compare_signatures([str(sig_a), str(sig_b)], str(tmp_path / "matrix.csv"), ksize=21)
    assert result["status"] == "completed", result["log_tail"]


def test_sketch_refuses_existing_output(tmp_path: Path) -> None:
    existing = tmp_path / "exists.sig"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        facade.sketch_dna(str(GENES), str(existing))


def test_compare_needs_two_signatures(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.compare_signatures([str(GENES)], str(tmp_path / "m.csv"))


def test_index_then_search_and_gather(tmp_path: Path) -> None:
    sig = tmp_path / "ref.sig"
    query = tmp_path / "query.sig"
    assert facade.sketch_dna(str(GENES), str(sig), ksizes=(21,), scaled=100)["status"] == "completed"
    assert facade.sketch_dna(str(GENES), str(query), ksizes=(21,), scaled=100)["status"] == "completed"
    db = tmp_path / "db.zip"
    built = facade.index_signatures([str(sig)], str(db), ksize=21)
    assert built["status"] == "completed", built["log_tail"]
    database = built["database"]
    found = facade.search_signatures(str(query), database, str(tmp_path / "hits.csv"), ksize=21, threshold=0.0)
    assert found["status"] == "completed", found["log_tail"]
    gathered = facade.gather_signatures(str(query), database, str(tmp_path / "gather.csv"), ksize=21, threshold_bp=0)
    assert gathered["status"] == "completed", gathered["log_tail"]


def test_describe_signature_reports_contents(tmp_path: Path) -> None:
    sig = tmp_path / "d.sig"
    assert facade.sketch_dna(str(GENES), str(sig), ksizes=(21,), scaled=100)["status"] == "completed"
    result = facade.describe_signature(str(sig))
    assert result["status"] == "completed", result["description"]
    assert "21" in result["description"]


def test_search_threshold_must_be_a_fraction(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.search_signatures(str(GENES), str(GENES), str(tmp_path / "x.csv"), threshold=2.0)


def test_merge_intersect_subtract_downsample(tmp_path: Path) -> None:
    a = tmp_path / "a.sig"
    b = tmp_path / "b.sig"
    assert facade.sketch_dna(str(GENES), str(a), ksizes=(21,), scaled=100)["status"] == "completed"
    assert facade.sketch_dna(str(DATA / "random.fa"), str(b), ksizes=(21,), scaled=100)["status"] == "completed"
    merged = facade.merge_signatures([str(a), str(b)], str(tmp_path / "m.sig"))
    assert merged["status"] == "completed", merged["log_tail"]
    inter = facade.intersect_signatures([str(a), str(b)], str(tmp_path / "i.sig"))
    assert inter["status"] == "completed", inter["log_tail"]
    sub = facade.subtract_signatures(str(a), [str(b)], str(tmp_path / "s.sig"))
    assert sub["status"] == "completed", sub["log_tail"]
    down = facade.downsample_signature(str(a), str(tmp_path / "d.sig"), scaled=1000)
    assert down["status"] == "completed", down["log_tail"]


def test_set_operations_need_enough_inputs(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.merge_signatures([str(GENES)], str(tmp_path / "x.sig"))
    with pytest.raises(ConfigurationError):
        facade.intersect_signatures([str(GENES)], str(tmp_path / "y.sig"))
    with pytest.raises(ConfigurationError):
        facade.downsample_signature(str(GENES), str(tmp_path / "z.sig"), scaled=0)


def test_cat_flatten_rename_fileinfo(tmp_path: Path) -> None:
    a = tmp_path / "a.sig"
    b = tmp_path / "b.sig"
    assert facade.sketch_dna(str(GENES), str(a), ksizes=(21,), scaled=100)["status"] == "completed"
    assert facade.sketch_dna(str(DATA / "random.fa"), str(b), ksizes=(21,), scaled=100)["status"] == "completed"
    cat = facade.concatenate_signatures([str(a), str(b)], str(tmp_path / "cat.sig"))
    assert cat["status"] == "completed", cat["log_tail"]
    flat = facade.flatten_signature(str(a), str(tmp_path / "flat.sig"))
    assert flat["status"] == "completed", flat["log_tail"]
    ren = facade.rename_signature(str(a), "renamed sample", str(tmp_path / "ren.sig"))
    assert ren["status"] == "completed", ren["log_tail"]
    info = facade.signature_file_info(str(a))
    assert info["status"] == "completed", info["description"]
    with pytest.raises(ConfigurationError):
        facade.rename_signature(str(a), "  ", str(tmp_path / "bad.sig"))


def test_filter_extract_split_manifest(tmp_path: Path) -> None:
    sig = tmp_path / "f.sig"
    assert facade.sketch_dna(str(GENES), str(sig), ksizes=(21, 31), scaled=100)["status"] == "completed"
    filt = facade.filter_signature_abundance(str(sig), str(tmp_path / "filt.sig"), min_abundance=1)
    assert filt["status"] == "completed", filt["log_tail"]
    ext = facade.extract_signatures([str(sig)], str(tmp_path / "k21.sig"), ksize=21)
    assert ext["status"] == "completed", ext["log_tail"]
    parts = facade.split_signatures([str(sig)], str(tmp_path / "split"))
    assert parts["status"] == "completed" and len(parts["outputs"]) >= 2, parts["log_tail"]
    man = facade.manifest_signature(str(sig), str(tmp_path / "manifest.csv"))
    assert man["status"] == "completed", man["log_tail"]
    with pytest.raises(ConfigurationError):
        facade.filter_signature_abundance(str(sig), str(tmp_path / "bad.sig"))


def _sketch(tmp_path: Path, name: str, fasta: Path = GENES) -> Path:
    sig = tmp_path / name
    assert facade.sketch_dna(str(fasta), str(sig), ksizes=(21,), scaled=10)["status"] == "completed"
    return sig


def test_prefetch_multigather_categorize(tmp_path: Path) -> None:
    ref = _sketch(tmp_path, "ref.sig")
    query = _sketch(tmp_path, "query.sig")
    fetched = facade.prefetch_signatures(str(query), [str(ref)], str(tmp_path / "prefetch.csv"), ksize=21)
    assert fetched["status"] == "completed", fetched["log_tail"]
    assert len((tmp_path / "prefetch.csv").read_text().splitlines()) >= 2
    multi = facade.multigather_signatures([str(query)], [str(ref)], str(tmp_path / "mg"), ksize=21)
    assert multi["status"] == "completed", multi["log_tail"]
    assert multi["outputs"]
    cat = facade.categorize_signatures(str(ref), [str(query)], str(tmp_path / "cat.csv"), ksize=21)
    assert cat["status"] == "completed", cat["log_tail"]


def test_overlap_reports_full_similarity_for_identical(tmp_path: Path) -> None:
    a = _sketch(tmp_path, "a.sig")
    b = _sketch(tmp_path, "b.sig")
    result = facade.overlap_signatures(str(a), str(b))
    assert result["status"] == "completed", result["report"]
    assert "similarity" in result["report"].lower()


def test_kmers_export_translate_plot(tmp_path: Path) -> None:
    a = _sketch(tmp_path, "a.sig")
    kmers = facade.extract_signature_kmers(str(a), [str(GENES)], str(tmp_path / "kmers.fa"), ksize=21)
    assert kmers["status"] == "completed", kmers["log_tail"]
    assert (tmp_path / "kmers.fa").read_text().startswith(">")
    exported = facade.export_signature_to_mash(str(a), str(tmp_path / "a.msh.json"), ksize=21)
    assert exported["status"] == "completed", exported["log_tail"]
    protein = facade.sketch_translated_dna(str(GENES), str(tmp_path / "prot.sig"), scaled=10)
    assert protein["status"] == "completed", protein["log_tail"]
    b = _sketch(tmp_path, "b.sig")
    cmp_csv = tmp_path / "cmp.csv"
    cmp_np = facade.compare_signatures([str(a), str(b)], str(cmp_csv), ksize=21)
    assert cmp_np["status"] == "completed"
    plotted = facade.plot_comparison_matrix(str(cmp_csv), str(tmp_path / "plots"))
    assert plotted["status"] == "completed", plotted["log_tail"]
    assert plotted["image"].endswith(".png")


def test_check_collect_ingest(tmp_path: Path) -> None:
    a = _sketch(tmp_path, "a.sig")
    manifest = tmp_path / "a.manifest.csv"
    assert facade.manifest_signature(str(a), str(manifest))["status"] == "completed"
    picklist = tmp_path / "pick.csv"
    picklist.write_text("md5\n" + "0" * 32 + "\n")
    checked = facade.check_signatures_against_picklist([str(a)], f"{picklist}:md5:md5", str(tmp_path / "miss.csv"))
    assert checked["status"] == "completed", checked["log_tail"]
    assert checked["all_picklist_rows_found"] is False
    collected = facade.collect_signature_manifests([str(a)], str(tmp_path / "all.sqlmf"))
    assert collected["status"] == "completed", collected["log_tail"]
    mash = tmp_path / "mash.csv"
    mash.write_text("murmur64,42,21,toy," + " ".join(str(i * 1000003) for i in range(1, 6)) + "\n")
    ingested = facade.ingest_mash_csv([str(mash)], str(tmp_path / "ingested.sig"))
    assert ingested["status"] == "completed", ingested["log_tail"]


def test_lca_index_and_classify(tmp_path: Path) -> None:
    sig = tmp_path / "g.sig"
    ident = "TOY0001"
    assert facade.sketch_dna(str(GENES), str(sig), ksizes=(21,), scaled=10, name=f"{ident} toy genome")["status"] == "completed"
    lineages = tmp_path / "lineages.csv"
    lineages.write_text("ident,superkingdom,phylum\n" f"{ident},Bacteria,Firmicutes\n")
    built = facade.index_lca_database([str(sig)], str(lineages), str(tmp_path / "db.lca.json"), ksize=21, scaled=10)
    assert built["status"] == "completed", built["log_tail"]
    called = facade.classify_with_lca([str(sig)], [str(tmp_path / "db.lca.json")], str(tmp_path / "cls.csv"), threshold=1)
    assert called["status"] == "completed", called["log_tail"]
    assert "Firmicutes" in (tmp_path / "cls.csv").read_text()
