"""Tests for the Kraken 2 facades; program-dependent tests skip when Kraken 2 is not installed.
No database is downloaded: a toy taxonomy and two random genomes are built locally."""

from __future__ import annotations

import random
import shutil

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.microbiome import kraken2_p14_cli_functions as k2

HAS = all(shutil.which(name) for name in ("kraken2", "kraken2-build", "kraken2-inspect", "k2"))
needs_kraken = pytest.mark.skipif(not HAS, reason="Kraken 2 (with k2) is not installed")

NODES = "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tgenus\t|\n11\t|\t10\t|\tspecies\t|\n12\t|\t10\t|\tspecies\t|\n"
NAMES = "".join(f"{i}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for i, n in ((1, "root"), (2, "Bacteria"), (10, "Genusa"), (11, "Genusa alpha"), (12, "Genusa beta")))


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    folder = tmp_path_factory.mktemp("kraken")
    rng = random.Random(3)
    genomes = {11: "".join(rng.choice("ACGT") for _ in range(4000)), 12: "".join(rng.choice("ACGT") for _ in range(4000))}
    files = []
    for taxid, sequence in genomes.items():
        path = folder / f"g{taxid}.fa"
        path.write_text(f">seq{taxid}|kraken:taxid|{taxid}\n{sequence}\n")
        files.append(str(path))
    database = folder / "db"
    (database / "taxonomy").mkdir(parents=True)
    (database / "taxonomy" / "nodes.dmp").write_text(NODES)
    (database / "taxonomy" / "names.dmp").write_text(NAMES)
    reads = "".join(f">a{i}\n{genomes[11][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
    reads += "".join(f">b{i}\n{genomes[12][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
    reads += "".join(f">u{i}\n" + "".join(rng.choice("ACGT") for _ in range(150)) + "\n" for i in range(5))
    read_file = folder / "reads.fa"
    read_file.write_text(reads)
    return folder, database, files, read_file


def test_missing_program_and_validation(tmp_path):
    with pytest.raises(ConfigurationError):
        k2.kraken2_classify_reads([], str(tmp_path), str(tmp_path / "o"))
    with pytest.raises(ConfigurationError):
        k2.kraken2_classify_reads([str(tmp_path / "missing.fa")], str(tmp_path), str(tmp_path / "o"))
    with pytest.raises(ConfigurationError):
        k2.kraken2_build_special_database(str(tmp_path / "d"), "gtdb")  # GTDB needs backend='k2'
    with pytest.raises(ConfigurationError):
        k2.kraken2_download_library(str(tmp_path / "d"), "bacteria", assembly_source="refseq")  # k2-only option
    reads = tmp_path / "r.fa"
    reads.write_text(">r\nACGT\n")
    fasta = tmp_path / "x.fa"
    fasta.write_text(">s|kraken:taxid|1\nACGT\n")
    with pytest.raises(BackendUnavailableError):
        k2.kraken2_add_to_library(str(tmp_path / "d2"), [str(fasta)], kraken2_path=str(tmp_path))


@needs_kraken
def test_build_classify_inspect_merge(toy, tmp_path):
    folder, database, files, read_file = toy
    added = k2.kraken2_add_to_library(str(database), files, no_masking=True, verbose=False)
    assert added["status"] == "completed" and len(added["runs"]) == 2
    built = k2.kraken2_build_database(str(database), kmer_length=31, minimizer_length=15, minimizer_spaces=3, threads=2, verbose=False)
    assert built["status"] == "completed" and len(built["outputs"]) == 3
    with pytest.raises(ConfigurationError):
        k2.kraken2_build_database(str(database), verbose=False)  # already built
    first = k2.kraken2_classify_reads([str(read_file)], str(database), str(tmp_path / "c1.txt"), report=str(tmp_path / "c1.rep"), verbose=False)
    assert first["status"] == "completed"
    assert (first["sequences_processed"], first["sequences_classified"], first["sequences_unclassified"]) == (45, 40, 5)
    by_taxid = {row["taxid"]: int(row["taxon_reads"]) for row in first["report_rows"]}
    assert by_taxid["11"] == 20 and by_taxid["12"] == 20
    quick = k2.kraken2_classify_reads([str(read_file)], str(database), str(tmp_path / "c2.txt"), report=str(tmp_path / "c2.rep"), quick=True, report_minimizer_data=True, verbose=False)
    assert quick["report_rows"][1]["distinct_minimizers"]
    multi = k2.kraken2_classify_reads([str(read_file)], [str(database), str(database)], str(tmp_path / "c3.txt"), merge_policy="resolve-lca", verbose=False)
    assert multi["status"] == "completed" and multi["command"][0] == "k2" and multi["sequences_classified"] == 40
    inspected = k2.kraken2_inspect_database(str(database), str(tmp_path / "ins.txt"), verbose=False)
    assert any("k = 31" in line for line in inspected["header"])
    assert {row["taxid"] for row in inspected["rows"]} >= {"1", "2", "10", "11", "12"}
    merged = k2.kraken2_merge_classifications([str(tmp_path / "c1.txt"), str(tmp_path / "c2.txt")], [str(database), str(database)], str(tmp_path / "m.txt"), report=str(tmp_path / "m.rep"), verbose=False)
    assert merged["status"] == "completed"
    assert len((tmp_path / "m.txt").read_text().splitlines()) == 45


@needs_kraken
@pytest.mark.skipif(not __import__("os").environ.get("DNAKIT_LIVE_TESTS"), reason="needs the network (DNAKIT_LIVE_TESTS=1)")
def test_live_download_univec_core(tmp_path):
    result = k2.kraken2_download_library(str(tmp_path / "db"), "UniVec_Core", use_ftp=True, no_masking=True, verbose=False)
    assert result["status"] == "completed"
