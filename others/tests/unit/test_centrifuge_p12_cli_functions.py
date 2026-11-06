"""Tests of the Centrifuge facades: validation, argument vectors (stand-in programs) and real Centrifuge runs."""

from __future__ import annotations

import json
import random
import shutil
import stat
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.microbiome import centrifuge_p12_cli_functions as cf

has_build = (
    shutil.which("centrifuge-build") is not None and shutil.which("centrifuge-inspect") is not None
)
has_perl_scripts = (
    shutil.which("perl") is not None and shutil.which("centrifuge-RemoveN.pl") is not None
)
has_compress = (
    shutil.which("centrifuge-compress.pl") is not None
    and shutil.which("jellyfish") is not None
    and shutil.which("nucmer") is not None
)
needs_build = pytest.mark.skipif(not has_build, reason="Centrifuge is not installed")
needs_scripts = pytest.mark.skipif(
    not has_perl_scripts, reason="Centrifuge Perl helpers are not installed"
)
needs_compress = pytest.mark.skipif(
    not has_compress, reason="Centrifuge, Jellyfish or MUMmer is missing"
)


def _genome(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _wrap(sequence: str, width: int = 70) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


@pytest.fixture()
def reference(tmp_path: Path) -> dict[str, Path]:
    sequences = {"seqA": _genome(1, 2500), "seqB": _genome(2, 1800), "seqC": _genome(3, 2200)}
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name} test\n{_wrap(seq)}\n" for name, seq in sequences.items()))
    conversion = tmp_path / "conv.tsv"
    conversion.write_text("seqA\t3\nseqB\t3\nseqC\t4\n")
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n3\t|\t2\t|\tspecies\t|\n4\t|\t2\t|\tspecies\t|\n"
    )
    names = tmp_path / "names.dmp"
    names.write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "3\t|\tSpecies A\t|\t\t|\tscientific name\t|\n4\t|\tSpecies B\t|\t\t|\tscientific name\t|\n"
    )
    return {"fa": fasta, "conv": conversion, "nodes": nodes, "names": names, "dir": tmp_path}


def _stand_in(folder: Path, name: str, stdout: str = "", create: tuple[str, ...] = ()) -> Path:
    """A fake Centrifuge program: records its arguments, prints ``stdout`` and creates the files in ``create``."""
    record = folder / "record.jsonl"
    script = folder / "bin" / name
    script.parent.mkdir(exist_ok=True)
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, pathlib, sys\n"
        f"open({str(record)!r}, 'a').write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:] == ['--version']:\n"
        "    print('centrifuge-class version 9.9.9'); sys.exit(0)\n"
        f"for path in {list(create)!r}:\n"
        "    target = pathlib.Path(path); target.parent.mkdir(parents=True, exist_ok=True); target.write_text('x')\n"
        f"sys.stdout.write({stdout!r})\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _calls(folder: Path) -> list[list[str]]:
    return [json.loads(line) for line in (folder / "record.jsonl").read_text().splitlines()]


def test_public_names_have_verbose() -> None:
    import inspect

    assert len(cf.__all__) == 9
    for name in cf.__all__:
        assert inspect.signature(getattr(cf, name)).parameters["verbose"].default is True


def test_missing_programs_raise(reference: dict[str, Path]) -> None:
    with pytest.raises(ConfigurationError):
        cf.centrifuge_remove_n(
            str(reference["fa"]),
            str(reference["dir"] / "o.fa"),
            centrifuge_dir="/nonexistent",
            verbose=False,
        )
    empty = reference["dir"] / "emptybin"
    empty.mkdir()
    with pytest.raises(BackendUnavailableError):
        cf.centrifuge_remove_n(
            str(reference["fa"]),
            str(reference["dir"] / "o.fa"),
            centrifuge_dir=str(empty),
            verbose=False,
        )


def test_build_validates_arguments(reference: dict[str, Path]) -> None:
    common = (str(reference["conv"]), str(reference["nodes"]))
    with pytest.raises(ConfigurationError):
        cf.centrifuge_build([], str(reference["dir"] / "idx"), *common, verbose=False)
    with pytest.raises(ConfigurationError):
        cf.centrifuge_build(
            [str(reference["fa"])],
            str(reference["dir"] / "idx"),
            *common,
            kmer_count=65,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        cf.centrifuge_build(
            [str(reference["fa"])], str(reference["dir"] / "idx"), *common, threads=0, verbose=False
        )
    with pytest.raises(ConfigurationError, match="folder"):
        cf.centrifuge_build(
            [str(reference["fa"])], str(reference["dir"] / "no" / "idx"), *common, verbose=False
        )
    for suffix in cf._INDEX_SUFFIXES:
        (reference["dir"] / f"done{suffix}").write_text("x")
    with pytest.raises(ConfigurationError, match="already exists"):
        cf.centrifuge_build(
            [str(reference["fa"])], str(reference["dir"] / "done"), *common, verbose=False
        )


def test_build_argument_vector_and_kmer_parsing(reference: dict[str, Path]) -> None:
    prefix = reference["dir"] / "idx"
    outputs = tuple(str(prefix) + suffix for suffix in cf._INDEX_SUFFIXES)
    stdout = "Number of distinct 1-mers is 4\nNumber of distinct 2-mers is 16\n"
    script = _stand_in(reference["dir"], "centrifuge-build", stdout, outputs)
    result = cf.centrifuge_build(
        [str(reference["fa"])],
        str(prefix),
        str(reference["conv"]),
        str(reference["nodes"]),
        name_table=str(reference["names"]),
        threads=2,
        offrate=4,
        seed=7,
        kmer_count=2,
        centrifuge_dir=str(script.parent),
        verbose=False,
    )
    call = _calls(reference["dir"])[0]
    assert call[:2] == ["-p", "2"]
    assert call[-2:] == [str(reference["fa"]), str(prefix)]
    assert call[call.index("--kmer-count") : call.index("--kmer-count") + 2] == [
        "--kmer-count",
        "2",
    ]
    assert ["--name-table", str(reference["names"])] == call[
        call.index("--name-table") : call.index("--name-table") + 2
    ]
    assert result["distinct_kmers"] == {"1": 4, "2": 16}
    assert result["version"] == "9.9.9" and len(result["index_files"]) == 4
    json.dumps(result)


def test_download_validates_and_builds_arguments(reference: dict[str, Path]) -> None:
    out = reference["dir"] / "dl"
    with pytest.raises(ConfigurationError):
        cf.centrifuge_download(str(out), "ncbi", verbose=False)
    with pytest.raises(ConfigurationError, match="domains"):
        cf.centrifuge_download(
            str(reference["dir"] / "dl1"), "refseq", domains=["mars"], verbose=False
        )
    with pytest.raises(ConfigurationError, match="genome options"):
        cf.centrifuge_download(
            str(reference["dir"] / "dl2"), "contaminants", include_rna=True, verbose=False
        )
    with pytest.raises(ConfigurationError):
        cf.centrifuge_download(
            str(reference["dir"] / "dl3"), "refseq", assembly_level="Draft", verbose=False
        )
    with pytest.raises(ConfigurationError):
        cf.centrifuge_download(
            str(reference["dir"] / "dl4"), "refseq", refseq_category="bad;rm", verbose=False
        )
    script = _stand_in(
        reference["dir"],
        "centrifuge-download",
        stdout="seq1\t562\nseq2\t562\n",
        create=(str(out / "bacteria" / "g.fna"),),
    )
    result = cf.centrifuge_download(
        str(out),
        "refseq",
        domains=["bacteria", "viral"],
        assembly_level="Chromosome",
        refseq_category="reference genome",
        taxids=[562, 9606],
        include_rna=True,
        mask_low_complexity=True,
        taxid_in_headers=True,
        threads=4,
        download_program="curl",
        centrifuge_dir=str(script.parent),
        verbose=False,
    )
    call = _calls(reference["dir"])[0]
    assert call == [
        "-o", str(out), "-P", "4", "-d", "bacteria,viral", "-a", "Chromosome", "-c", "reference genome",
        "-t", "562,9606", "-r", "-m", "-l", "-g", "curl", "refseq",
    ]  # fmt: skip
    assert result["mapping_rows"] == 2 and result["files"] == ["bacteria/g.fna"]
    assert result["mapping_file"] == str(out / "seqid_to_taxid.map")
    # an existing output directory is refused
    with pytest.raises(ConfigurationError, match="must not exist"):
        cf.centrifuge_download(
            str(out), "taxonomy", centrifuge_dir=str(script.parent), verbose=False
        )


def test_failing_program_raises_with_log(reference: dict[str, Path]) -> None:
    script = reference["dir"] / "centrifuge-RemoveN.pl"
    script.write_text(
        f"#!{sys.executable}\nimport sys\nsys.stderr.write('bad input')\nsys.exit(2)\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    with pytest.raises(BackendExecutionError) as caught:
        cf.centrifuge_remove_n(
            str(reference["fa"]),
            str(reference["dir"] / "o.fa"),
            centrifuge_dir=str(reference["dir"]),
            verbose=False,
        )
    assert "bad input" in caught.value.context["log_tail"]


@needs_build
def test_build_inspect_and_kreport_roundtrip(reference: dict[str, Path]) -> None:
    prefix = reference["dir"] / "idx"
    built = cf.centrifuge_build(
        [str(reference["fa"])],
        str(prefix),
        str(reference["conv"]),
        str(reference["nodes"]),
        name_table=str(reference["names"]),
        kmer_count=8,
        verbose=False,
    )
    assert built["status"] == "completed" and built["version"].startswith("1.")
    assert set(built["distinct_kmers"]) == {str(k) for k in range(1, 9)}
    assert built["distinct_kmers"]["1"] == 4 and built["distinct_kmers"]["2"] == 16
    fasta = cf.centrifuge_inspect(str(prefix), str(reference["dir"] / "i.fa"), verbose=False)
    assert (fasta["record_count"], fasta["base_count"]) == (3, 2500 + 1800 + 2200)
    names = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "n.txt"), output="names", verbose=False
    )
    assert names["records"] == ["seqA test", "seqB test", "seqC test"]  # the header text is kept
    summary = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "s.txt"), output="summary", verbose=False
    )
    assert summary["records"][1] == {"index": 2, "name": "seqB test", "length": 1800}
    assert summary["parameters"]["Flags"] == 1
    conv = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "c.txt"), output="conversion_table", verbose=False
    )
    assert conv["records"] == [
        {"sequence": "seqA", "taxid": 3},
        {"sequence": "seqB", "taxid": 3},
        {"sequence": "seqC", "taxid": 4},
    ]
    sizes = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "z.txt"), output="size_table", verbose=False
    )
    assert sizes["records"] == [
        {"taxid": 3, "total_length": 4300},
        {"taxid": 4, "total_length": 2200},
    ]
    tree = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "t.txt"), output="taxonomy_tree", verbose=False
    )
    assert tree["record_count"] == 4 and tree["records"][2]["rank"] == "species"
    table = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "m.txt"), output="name_table", verbose=False
    )
    assert table["records"][2] == {"taxid": 3, "name": "Species A"}
    estimate = cf.centrifuge_inspect(
        str(prefix), str(reference["dir"] / "e.txt"), output="estimate_n_kmers", verbose=False
    )
    assert 0.8 * 6500 < estimate["estimated_kmers"] < 1.2 * 6500  # about 6.5 kb of random sequence
    counts = reference["dir"] / "counts.tsv"
    counts.write_text("3\t6\n4\t2\n")
    report = cf.centrifuge_kreport(
        str(prefix),
        [str(counts)],
        str(reference["dir"] / "report.txt"),
        is_count_table=True,
        verbose=False,
    )
    by_taxid = {row["taxid"]: row for row in report["records"]}
    assert by_taxid[3]["clade_reads"] == 6 and by_taxid[4]["percent"] == 25.0
    assert (
        by_taxid[2]["name"] == "Bacteria"
        and by_taxid[2]["clade_reads"] == 8
        and by_taxid[3]["depth"] == 2
    )
    with pytest.raises(ConfigurationError, match=r"complete Centrifuge index|No complete"):
        cf.centrifuge_inspect(
            str(reference["dir"] / "nope"), str(reference["dir"] / "x.txt"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        cf.centrifuge_inspect(
            str(prefix), str(reference["dir"] / "y.txt"), output="bogus", verbose=False
        )


@needs_scripts
def test_remove_n_and_remove_empty_sequences(tmp_path: Path) -> None:
    source = tmp_path / "n.fa"
    source.write_text(">r1 d\nACGTNNacgnnTTAA\nNNNN\n>r2\nNNNN\n>r3\nGGCCnnTT\n")
    cleaned = cf.centrifuge_remove_n(str(source), str(tmp_path / "noN.fa"), verbose=False)
    assert (tmp_path / "noN.fa").read_text() == ">r1 d\nACGTacgTTAA\n>r2\n>r3\nGGCCTT\n"
    assert (cleaned["record_count"], cleaned["base_count"]) == (3, 17)
    wide = tmp_path / "wide.fa"
    wide.write_text(">w\n" + "A" * 200 + "N" * 5 + "\n")
    cf.centrifuge_remove_n(str(wide), str(tmp_path / "wide80.fa"), verbose=False)
    assert [len(line) for line in (tmp_path / "wide80.fa").read_text().splitlines()[1:]] == [
        80,
        80,
        40,
    ]
    result = cf.centrifuge_remove_empty_sequences(
        str(tmp_path / "noN.fa"), str(tmp_path / "nonempty.fa"), verbose=False
    )
    assert (tmp_path / "nonempty.fa").read_text() == ">r1 d\nACGTacgTTAA\n>r3\nGGCCTT\n"
    assert result["record_count"] == 2 and result["version"]
    every = tmp_path / "every.fa"
    every.write_text(">e1\n>e2\nACGT\n>e3\n>e4\n>e5\nTT\n")
    cf.centrifuge_remove_empty_sequences(
        str(every), str(tmp_path / "every_clean.fa"), verbose=False
    )
    assert (tmp_path / "every_clean.fa").read_text() == ">e2\nACGT\n>e5\nTT\n"


@needs_scripts
def test_sort_nt_orders_records_by_taxid(tmp_path: Path) -> None:
    nt = tmp_path / "nt.fa"
    nt.write_text(
        ">AC1.1 first\nACGT\nACGT\n>AC2.1 second\nGGGG\n>AC3.1 third\nTTTT\n>AC4.1 unmapped\nCCCC\n"
    )
    mapping = tmp_path / "a2t.txt"
    mapping.write_text(
        "accession\taccession.version\ttaxid\tgi\nAC1\tAC1.1\t30\t0\nAC2\tAC2.1\t10\t0\nAC3\tAC3.1\t30\t0\n"
    )
    result = cf.centrifuge_sort_nt(
        str(nt),
        [str(mapping)],
        str(tmp_path / "sorted.fa"),
        taxid_map_output=str(tmp_path / "map.txt"),
        verbose=False,
    )
    assert (
        tmp_path / "sorted.fa"
    ).read_text() == ">AC2.1 second\nGGGG\n>AC1.1 first\nACGT\nACGT\n>AC3.1 third\nTTTT\n"
    assert result["record_count"] == 3
    rows = sorted(line.split("\t") for line in (tmp_path / "map.txt").read_text().splitlines())
    assert rows == [["AC1.1", "30"], ["AC2.1", "10"], ["AC3.1", "30"]]


@needs_compress
def test_compress_genomes_and_shared_sequence(tmp_path: Path) -> None:
    rng = random.Random(5)
    base = _genome(11, 6000)
    variant = list(base)
    for position in rng.sample(range(6000), 10):
        variant[position] = rng.choice("ACGT")
    other = _genome(12, 6000)
    download = tmp_path / "dl"
    download.mkdir()
    for name, sequence in (("gA1", base), ("gA2", "".join(variant)), ("gB1", other)):
        (download / f"{name}.fna").write_text(f">{name}_chr1 test\n{_wrap(sequence)}\n")
    taxonomy = tmp_path / "tax"
    taxonomy.mkdir()
    (taxonomy / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tgenus\t|\n3\t|\t2\t|\tspecies\t|\n4\t|\t2\t|\tspecies\t|\n"
    )
    (taxonomy / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tGenusX\t|\t\t|\tscientific name\t|\n"
        "3\t|\tSpecies A\t|\t\t|\tscientific name\t|\n4\t|\tSpecies B\t|\t\t|\tscientific name\t|\n"
    )
    mapping = tmp_path / "map.txt"
    mapping.write_text("gA1_chr1\t3\ngA2_chr1\t3\ngB1_chr1\t4\n")
    result = cf.centrifuge_compress_genomes(
        str(download),
        str(taxonomy),
        str(tmp_path / "comp"),
        mapping_file=str(mapping),
        no_dustmasker=True,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert (result["sequence_count"], result["base_count"]) == (
        2,
        12000,
    )  # the two species A genomes are merged
    sizes = dict(line.split() for line in Path(result["size_table"]).read_text().splitlines())
    assert sizes == {"3": "12000", "4": "6000"}
    shared = cf.centrifuge_shared_sequence(
        [str(download / "gA1.fna"), str(download / "gA2.fna"), str(download / "gB1.fna")],
        str(tmp_path / "shared"),
        verbose=False,
    )
    assert shared["status"] == "completed" and shared["files"]
    assert all(item["records"] >= 1 for item in shared["files"])
