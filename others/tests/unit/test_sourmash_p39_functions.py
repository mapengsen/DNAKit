"""Tests for dnakit.microbiome.sourmash_p39_functions (sourmash tax, index, watch and sketch operations)."""

from __future__ import annotations

import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("sourmash")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.microbiome import sourmash_p39_functions as sm  # noqa: E402

PROGRAM = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
pytestmark = pytest.mark.skipif(not Path(PROGRAM).exists(), reason="sourmash program not found")


def _cli(*args: str, cwd: Path) -> None:
    subprocess.run([PROGRAM, *args], cwd=cwd, check=True, capture_output=True)


@pytest.fixture(scope="module")
def workspace(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("sm_p39")
    rng = random.Random(5)
    genomes = {f"g{i}": "".join(rng.choice("ACGT") for _ in range(20000)) for i in range(3)}
    for name, seq in genomes.items():
        (root / f"{name}.fa").write_text(f">{name}\n{seq}\n")
        _cli("sketch", "dna", "-p", "k=31,scaled=10", f"{name}.fa", "--name", name, "-o", f"{name}.sig", cwd=root)
        _cli("sketch", "dna", "-p", "k=31,num=500", f"{name}.fa", "--name", name, "-o", f"n{name}.sig", cwd=root)
    (root / "q.fa").write_text(">q\n" + genomes["g0"] + genomes["g1"][:8000] + "\n")
    _cli("sketch", "dna", "-p", "k=31,scaled=10", "q.fa", "--name", "q", "-o", "q.sig", cwd=root)
    _cli("index", "db", "g0.sig", "g1.sig", "g2.sig", "-k", "31", cwd=root)
    _cli("gather", "q.sig", "db.sbt.zip", "-o", "gather.csv", "--threshold-bp", "100", cwd=root)
    _cli("index", "ndb", "ng0.sig", "ng1.sig", "ng2.sig", "-k", "31", cwd=root)
    (root / "tax.csv").write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species\n"
        + "".join(f"g{i},Bacteria,P{i // 2},C,O,F{i // 2},G{i // 2},S{i}\n" for i in range(3))
    )
    reads = []
    for i in range(200):
        p = (i * 97) % 19850
        reads.append(f"@r{i}\n{genomes['g1'][p:p + 150]}\n+\n{'I' * 150}\n")
    (root / "r.fq").write_text("".join(reads))
    return root


def test_tax_genome_formats(workspace: Path) -> None:
    result = sm.tax_genome_classification(
        [str(workspace / "gather.csv")], [str(workspace / "tax.csv")], "out", str(workspace / "genome_out"),
        output_formats=["csv_summary", "krona", "lineage_csv"], rank="species", verbose=False,
    )
    assert result["status"] == "completed"
    names = sorted(Path(p).name for p in result["outputs"])
    assert names == ["out.classifications.csv", "out.krona.tsv", "out.lineage.csv"]
    assert "Bacteria;P0;C;O;F0;G0;S0" in (workspace / "genome_out" / "out.classifications.csv").read_text()


def test_tax_genome_validation(workspace: Path) -> None:
    with pytest.raises(ConfigurationError):
        sm.tax_genome_classification(
            [str(workspace / "gather.csv")], [str(workspace / "tax.csv")], "o", str(workspace / "x"),
            output_formats=["krona"], verbose=False,
        )


def test_tax_annotate_prepare_grep_summarize(workspace: Path) -> None:
    annotated = sm.tax_annotate_gather_results(
        [str(workspace / "gather.csv")], [str(workspace / "tax.csv")], str(workspace / "ann"), verbose=False
    )
    assert annotated["status"] == "completed"
    assert "lineage" in Path(annotated["outputs"][0]).read_text().splitlines()[0]
    prepared = sm.tax_prepare_taxonomy_database(
        [str(workspace / "tax.csv")], str(workspace / "tax.db"), verbose=False
    )
    assert prepared["status"] == "completed"
    grep = sm.tax_grep_lineages("S1", [str(workspace / "tax.csv")], output=str(workspace / "grep.csv"), verbose=False)
    assert grep["match_count"] == 1
    assert "g1," in (workspace / "grep.csv").read_text()
    counted = sm.tax_grep_lineages("Bacteria", [str(workspace / "tax.csv")], verbose=False)
    assert counted["match_count"] == 3
    summary = sm.tax_summarize_lineages(
        [str(workspace / "tax.csv")], output_lineage_counts=str(workspace / "lin.csv"), verbose=False
    )
    assert summary["distinct_lineages"] == 3
    assert summary["rank_counts"]["species"] == 3


def test_rocksdb_index_and_existing_output(workspace: Path) -> None:
    result = sm.build_signature_index(
        [str(workspace / f"g{i}.sig") for i in range(3)], str(workspace / "rdb"), verbose=False
    )
    assert result["status"] == "completed"
    assert (workspace / "rdb").is_dir()
    with pytest.raises(ConfigurationError):
        sm.build_signature_index([str(workspace / "g0.sig")], str(workspace / "rdb"), verbose=False)


def test_watch_finds_source_genome(workspace: Path) -> None:
    result = sm.watch_reads_against_index(
        str(workspace / "ndb.sbt.zip"), str(workspace / "r.fq"), str(workspace / "w.sig"),
        ksize=31, threshold=0.001, verbose=False,
    )
    assert result["status"] == "completed"
    assert result["match"]["name"] == "g1"


def _mh(hashes: list[int], abundances: list[int] | None = None) -> dict:
    sketch: dict = {"ksize": 31, "scaled": 1, "hashes": hashes}
    if abundances is not None:
        sketch["abundances"] = abundances
    return sketch


def test_weighted_containment_and_union() -> None:
    weighted = sm.sketch_weighted_containment(_mh([1, 2, 3, 4], [1, 1, 2, 4]), _mh([2, 4, 9]), verbose=False)
    assert weighted["weighted_containment"] == pytest.approx(5 / 8)
    sizes = sm.sketch_intersection_union_size(_mh([1, 2, 3]), _mh([2, 3, 4, 5]), verbose=False)
    assert (sizes["shared_hashes"], sizes["intersection"], sizes["union"]) == (2, 2, 5)
    with pytest.raises(ConfigurationError):
        sm.sketch_weighted_containment(_mh([1, 2]), _mh([1]), verbose=False)


def test_kmer_hashes_match_api() -> None:
    import sourmash

    plain = sm.sequence_kmer_hashes("ACGTACGTAC", 5, verbose=False)
    reference = sourmash.MinHash(n=0, ksize=5, scaled=1).seq_to_hashes("ACGTACGTAC")
    assert plain["hashes"] == [int(h) for h in reference]
    paired = sm.sequence_kmer_hashes("acgtacgtac", 5, include_kmers=True, verbose=False)
    assert paired["kmers"][0][0] == "ACGTA"
    assert paired["kmer_count"] == 6


def test_size_accuracy_and_murmur() -> None:
    small = sm.sketch_size_accuracy({"ksize": 31, "scaled": 1000, "hashes": list(range(5))}, verbose=False)
    assert small["is_accurate"] is False
    large = sm.sketch_size_accuracy({"ksize": 31, "scaled": 10, "hashes": list(range(5000))}, verbose=False)
    assert large["is_accurate"] is True
    from sourmash.minhash import hash_murmur

    assert sm.murmur_hash_string("ACGTA", verbose=False)["hash"] == hash_murmur("ACGTA", 42)
    with pytest.raises(ConfigurationError):
        sm.murmur_hash_string("", verbose=False)


def test_ingest_mash_json_round_trip(workspace: Path) -> None:
    exported = subprocess.run(
        [PROGRAM, "sig", "export", "ng0.sig", "-k", "31"], cwd=workspace, check=True, capture_output=True, text=True
    ).stdout
    (workspace / "mash.json").write_text(exported)
    result = sm.ingest_mash_json_signatures([str(workspace / "mash.json")], str(workspace / "mash.sig"), verbose=False)
    assert result["status"] == "completed"
    described = subprocess.run(
        [PROGRAM, "sig", "describe", "mash.sig"], cwd=workspace, check=True, capture_output=True, text=True
    ).stdout
    assert "num=500" in described
