"""Tests for the VSEARCH command-line facades (pack P09)."""

from __future__ import annotations

import random
import shutil
import struct
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import vsearch_p09_cli_functions as vs

needs_vsearch = pytest.mark.skipif(
    shutil.which("vsearch") is None, reason="vsearch is not installed"
)

_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(text: str) -> str:
    return text.translate(_COMP)[::-1]


def _random(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


@pytest.fixture()
def toy(tmp_path: Path) -> dict[str, Path]:
    rng = random.Random(1)
    a, b = _random(rng, 300), _random(rng, 300)
    variant = a[:290] + ("T" if a[290] != "T" else "A")
    fasta = tmp_path / "toy.fa"
    fasta.write_text(
        f">a;size=50\n{a}\n>b;size=30\n{b}\n>a2;size=5\n{variant}\n>chim;size=3\n{a[:150] + b[150:]}\n>dup1;size=2\n{a}\n"
    )
    site = tmp_path / "site.fa"
    site.write_text(f">withsite\n{a[:50]}GAATTC{b[:60]}\n>nosite\n{b[:40]}\n")
    taxonomy = tmp_path / "tax.fa"
    taxonomy.write_text(f">t1;tax=d:Bac,p:Fir,c:Cl\n{a}\n>t2;tax=d:Bac,p:Pro,c:Gam\n{b}\n")
    forward, reverse, orphan = [], [], []
    for index in range(5):
        fragment = _random(rng, 120)
        quality = "I" * 80
        forward.append(f"@p{index}/1\n{fragment[:80]}\n+\n{quality}\n")
        reverse.append(f"@p{index}/2\n{_revcomp(fragment[-80:])}\n+\n{quality}\n")
    orphan.append(f"@orph/1\n{_random(rng, 80)}\n+\n{'I' * 80}\n")
    paths = {
        "fasta": fasta,
        "site": site,
        "taxonomy": taxonomy,
        "forward": tmp_path / "r1.fq",
        "reverse": tmp_path / "r2.fq",
        "forward_orphan": tmp_path / "r1_orphan.fq",
    }
    paths["forward"].write_text("".join(forward))
    paths["reverse"].write_text("".join(reverse))
    paths["forward_orphan"].write_text("".join(forward + orphan))
    reads = tmp_path / "reads.fq"
    lines = []
    for index in range(40):
        length = 60 + index
        qualities = "".join(chr(33 + (38 if position < 40 else 12)) for position in range(length))
        lines.append(f"@read{index}\n{_random(rng, length)}\n+\n{qualities}\n")
    reads.write_text("".join(lines))
    paths["reads"] = reads
    return paths


def _sff(path: Path) -> Path:
    flows, key = "TACG" * 3, "TCAG"

    def pad(count: int) -> int:
        return (-count) % 8

    header_length = 31 + len(flows) + len(key)
    header_length += pad(header_length)
    header = struct.pack(
        ">IIQIIHHHB", 0x2E736666, 1, 0, 0, 2, header_length, len(key), len(flows), 1
    )
    header += flows.encode() + key.encode()
    header += b"\0" * (header_length - len(header))
    blocks = [header]
    for name, bases, qualities in (
        ("r1", "ACGTACGTAC", [30] * 10),
        ("r2", "GGGTTTAAAC", [20] * 5 + [35] * 5),
    ):
        read_header_length = 16 + len(name)
        read_header_length += pad(read_header_length)
        read_header = struct.pack(
            ">HHIHHHH", read_header_length, len(name), len(bases), 1, len(bases), 0, 0
        )
        read_header += name.encode()
        read_header += b"\0" * (read_header_length - len(read_header))
        data = (
            struct.pack(">12H", *([100] * 12))
            + bytes([1] * len(bases))
            + bases.encode()
            + bytes(qualities)
        )
        data += b"\0" * pad(len(data))
        blocks += [read_header, data]
    path.write_bytes(b"".join(blocks))
    return path


# ---------------------------------------------------------------------------
# Tests that need no VSEARCH program
# ---------------------------------------------------------------------------


def test_missing_program_raises_before_creating_output(
    tmp_path: Path, toy: dict[str, Path]
) -> None:
    target = tmp_path / "out"
    with pytest.raises(BackendUnavailableError):
        vs.vsearch_sort_by_size(
            toy["fasta"], target, vsearch_path=str(tmp_path / "no-such-vsearch"), verbose=False
        )
    assert not target.exists()


def test_validation_errors(tmp_path: Path, toy: dict[str, Path]) -> None:
    fake = tmp_path / "fakevsearch"
    fake.write_text("#!/bin/sh\nexit 0\n")
    fake.chmod(0o755)
    with pytest.raises(ConfigurationError):
        vs.vsearch_detect_chimeras(
            toy["fasta"], tmp_path / "o1", method="uchime_ref", vsearch_path=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vs.vsearch_detect_chimeras(
            toy["fasta"], tmp_path / "o2", method="uchime_denovo", chimeras_parts=4, vsearch_path=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vs.vsearch_cut_restriction_sites(
            toy["fasta"], tmp_path / "o3", cut_pattern="GAATTC", vsearch_path=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vs.vsearch_extract_records(toy["fasta"], tmp_path / "o4", vsearch_path=str(fake), verbose=False)
    with pytest.raises(ConfigurationError):
        vs.vsearch_fastq_eestats2(
            toy["reads"], tmp_path / "o5", length_cutoffs="50;100", vsearch_path=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vs.vsearch_dereplicate_fulllength(
            toy["fasta"], tmp_path / "o6", method="bogus", vsearch_path=str(fake), verbose=False
        )
    with pytest.raises(ConfigurationError):
        vs.vsearch_sort_by_size(tmp_path / "missing.fa", tmp_path / "o7", vsearch_path=str(fake), verbose=False)


def test_public_names_have_verbose() -> None:
    import inspect

    for name in vs.__all__:
        signature = inspect.signature(getattr(vs, name))
        assert "verbose" in signature.parameters, name
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY


# ---------------------------------------------------------------------------
# Tests that run VSEARCH
# ---------------------------------------------------------------------------


@needs_vsearch
def test_allpairs_and_search(tmp_path: Path, toy: dict[str, Path]) -> None:
    pairs = vs.vsearch_allpairs_global(toy["fasta"], tmp_path / "pairs", verbose=False)
    assert pairs["status"] == "completed" and pairs["record_count"] == 10
    assert {"query", "target", "identity"} <= set(pairs["records"][0])
    found = vs.vsearch_search_global(
        toy["fasta"],
        toy["fasta"],
        tmp_path / "search",
        min_identity=0.9,
        abundance_in=True,
        write_alignments=True,
        write_uc=True,
        write_otu_table=True,
        write_matched_queries=True,
        verbose=False,
    )
    assert found["record_count"] == 5
    assert all(Path(path).is_file() for path in found["outputs"].values())
    top = {row["query"]: row["target"] for row in found["records"]}
    assert top["a;size=50"] == "a;size=50"


@needs_vsearch
def test_chimera_methods(tmp_path: Path, toy: dict[str, Path]) -> None:
    for method in ("uchime_denovo", "uchime2_denovo", "chimeras_denovo"):
        result = vs.vsearch_detect_chimeras(
            toy["fasta"], tmp_path / method, method=method, abundance_in=True, verbose=False
        )
        assert result["method"] == method
        assert result["chimera_count"] + result["nonchimera_count"] == 5
    de_novo = vs.vsearch_detect_chimeras(
        toy["fasta"], tmp_path / "check", method="uchime_denovo", abundance_in=True, verbose=False
    )
    flagged = [row["query"] for row in de_novo["records"] if row["chimeric"] == "Y"]
    assert flagged == ["chim;size=3"]
    reference = vs.vsearch_detect_chimeras(
        toy["fasta"],
        tmp_path / "ref",
        method="uchime_ref",
        database=toy["fasta"],
        exclude_same_label=True,
        verbose=False,
    )
    assert reference["chimera_count"] == 1


@needs_vsearch
def test_unoise_and_dereplication(tmp_path: Path, toy: dict[str, Path]) -> None:
    unoise = vs.vsearch_cluster_unoise(
        toy["fasta"],
        tmp_path / "unoise",
        min_abundance=2,
        abundance_in=True,
        abundance_out=True,
        verbose=False,
    )
    assert unoise["centroid_count"] == 3
    full = vs.vsearch_dereplicate_fulllength(
        toy["fasta"], tmp_path / "full", abundance_in=True, abundance_out=True, verbose=False
    )
    assert full["unique_count"] == 4
    assert ">a;size=52" in Path(full["outputs"]["unique"]).read_text()
    small = vs.vsearch_dereplicate_fulllength(
        toy["fasta"], tmp_path / "small", method="smallmem", verbose=False
    )
    assert small["unique_count"] == 4 and "clusters" not in small["outputs"]
    fastq = vs.vsearch_dereplicate_fulllength(
        toy["reads"], tmp_path / "uniques", method="fastx_uniques", fastq_output=True, verbose=False
    )
    assert fastq["unique_count"] == 40
    assert (
        vs.vsearch_dereplicate_by_id(toy["fasta"], tmp_path / "byid", verbose=False)["unique_count"]
        == 5
    )
    assert (
        vs.vsearch_dereplicate_prefix(toy["fasta"], tmp_path / "prefix", verbose=False)[
            "unique_count"
        ]
        == 4
    )
    assert (
        vs.vsearch_rereplicate(toy["fasta"], tmp_path / "rerep", verbose=False)["record_count"]
        == 90
    )
    ranked = vs.vsearch_sort_by_size(toy["fasta"], tmp_path / "sorted", min_size=3, verbose=False)
    assert ranked["record_count"] == 4


@needs_vsearch
def test_fastq_statistics_and_filtering(tmp_path: Path, toy: dict[str, Path]) -> None:
    ee = vs.vsearch_fastq_eestats(toy["reads"], tmp_path / "ee", verbose=False)
    assert ee["position_count"] == 99 and ee["records"][0]["recs"] == 40
    ee2 = vs.vsearch_fastq_eestats2(
        toy["reads"],
        tmp_path / "ee2",
        ee_cutoffs=[0.5, 1],
        length_cutoffs="50,100,25",
        verbose=False,
    )
    assert ee2["reads"] == 40 and [row["length"] for row in ee2["records"]] == [50, 75]
    stats = vs.vsearch_fastq_stats(toy["reads"], tmp_path / "stats", verbose=False)
    assert stats["read_count"] == 40 and stats["read_length_distribution"][0]["length"] == 99
    assert sum(row["bases"] for row in stats["quality_distribution"]) == sum(
        60 + i for i in range(40)
    )
    kept = vs.vsearch_filter_reads(
        toy["reads"], tmp_path / "filter", min_length=80, write_discarded=True, verbose=False
    )
    assert kept["passed_count"] == 20 and kept["discarded_count"] == 20
    truncated = vs.vsearch_filter_reads(
        toy["reads"], tmp_path / "trunc", truncate_quality=20, verbose=False
    )
    lengths = {
        len(line) for line in Path(truncated["outputs"]["filtered"]).read_text().splitlines()[1::4]
    }
    assert lengths == {40}
    by_abundance = vs.vsearch_filter_reads(
        toy["fasta"], tmp_path / "abund", min_size=5, verbose=False
    )
    assert by_abundance["passed_count"] == 3


@needs_vsearch
def test_pair_operations(tmp_path: Path, toy: dict[str, Path]) -> None:
    merged = vs.vsearch_merge_pairs(
        toy["forward"], toy["reverse"], tmp_path / "merge", write_error_table=True, verbose=False
    )
    assert merged["merged_count"] == 5 and merged["not_merged_count"] == 0
    first = Path(merged["outputs"]["merged"]).read_text().splitlines()[1]
    assert len(first) == 120
    joined = vs.vsearch_join_pairs(
        toy["forward"], toy["reverse"], tmp_path / "join", write_fasta=True, verbose=False
    )
    sequence = Path(joined["outputs"]["joined"]).read_text().splitlines()[1]
    assert len(sequence) == 80 + 8 + 80 and "NNNNNNNN" in sequence
    synced = vs.vsearch_sync_pairs(
        toy["forward_orphan"], toy["reverse"], tmp_path / "sync", verbose=False
    )
    assert (synced["synced_pairs"], synced["forward_orphans"], synced["reverse_orphans"]) == (
        5,
        1,
        0,
    )


@needs_vsearch
def test_retrieval_orientation_and_editing(tmp_path: Path, toy: dict[str, Path]) -> None:
    one = vs.vsearch_extract_records(
        toy["fasta"], tmp_path / "get", label="a;size=50", write_unmatched=True, verbose=False
    )
    assert one["extracted_count"] == 1
    assert _count_headers(one["outputs"]["unmatched"]) == 4
    sub = vs.vsearch_extract_subsequence(
        toy["fasta"], tmp_path / "sub", label="a;size=50", start=1, end=20, verbose=False
    )
    assert len(Path(sub["outputs"]["subsequence"]).read_text().splitlines()[1]) == 20
    reverse_input = tmp_path / "reverse.fa"
    original = Path(toy["fasta"]).read_text().splitlines()[1]
    reverse_input.write_text(f">flip\n{_revcomp(original)}\n>same\n{original}\n")
    oriented = vs.vsearch_orient_sequences(
        reverse_input, toy["fasta"], tmp_path / "orient", verbose=False
    )
    assert [row["strand"] for row in oriented["records"]] == ["-", "+"]
    assert _fasta(oriented["outputs"]["oriented"]) == [original, original]
    cut = vs.vsearch_cut_restriction_sites(
        toy["site"], tmp_path / "cut", cut_pattern="G^AATT_C", verbose=False
    )
    assert cut["forward_fragment_count"] == 2
    assert _count_headers(cut["outputs"]["discarded_forward"]) == 1
    scrambled = vs.vsearch_scramble_sequences(
        toy["fasta"], tmp_path / "scramble", random_seed=3, verbose=False
    )
    new = _fasta(scrambled["outputs"]["scrambled"])[0]
    assert sorted(new) == sorted(original) and new != original
    converted = vs.vsearch_convert_sff(_sff(tmp_path / "toy.sff"), tmp_path / "sff", verbose=False)
    lines = Path(converted["outputs"]["converted"]).read_text().splitlines()
    assert lines[:4] == ["@r1", "ACGTACGTAC", "+", "??????????"]


@needs_vsearch
def test_sintax_and_udb(tmp_path: Path, toy: dict[str, Path]) -> None:
    result = vs.vsearch_sintax(
        toy["fasta"], toy["taxonomy"], tmp_path / "sintax", cutoff=0.5, random_seed=1, verbose=False
    )
    by_query = {row["query"]: row["taxonomy"] for row in result["records"]}
    assert by_query["a;size=50"] == "d:Bac,p:Fir,c:Cl"
    assert by_query["b;size=30"] == "d:Bac,p:Pro,c:Gam"
    built = vs.vsearch_make_udb(toy["fasta"], tmp_path / "udb", word_length=6, verbose=False)
    database = built["outputs"]["database"]
    exported = vs.vsearch_udb_to_fasta(database, tmp_path / "export", verbose=False)
    assert exported["record_count"] == 5
    info = vs.vsearch_udb_report(database, tmp_path / "info", verbose=False)
    assert info["fields"]["Seqs"] == 5 and info["fields"]["Word width"] == 6
    stats = vs.vsearch_udb_report(database, tmp_path / "stats", report="stats", verbose=False)
    assert "Word width" in stats["report_text"] and "DB size" in stats["report_text"]
    hits = vs.vsearch_search_global(
        toy["fasta"], database, tmp_path / "udbsearch", min_identity=0.9, verbose=False
    )
    assert hits["record_count"] == 5


@needs_vsearch
def test_vsearch_failure_raises(tmp_path: Path, toy: dict[str, Path]) -> None:
    with pytest.raises(BackendExecutionError):
        vs.vsearch_fastq_eestats(toy["fasta"], tmp_path / "bad", verbose=False)


def _count_headers(path: str) -> int:
    return sum(1 for line in Path(path).read_text().splitlines() if line.startswith(">"))


def _fasta(path: str) -> list[str]:
    """Sequences of a FASTA file with the 80-column wrapping removed."""
    sequences: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            sequences.append("")
        elif sequences:
            sequences[-1] += line
    return sequences
