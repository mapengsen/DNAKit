"""Tests for the BLAST+ facades; program-dependent tests skip when BLAST+ is not installed."""

from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from dnakit.alignment import blast_p14_cli_functions as blast
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

HAS_BLAST = shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None
needs_blast = pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ is not installed")


def _genes(path, count=6, length=900, seed=1):
    rng = random.Random(seed)
    records = ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(count)]
    path.write_text("".join(f">cds{i}\n{seq}\n" for i, seq in enumerate(records)))
    return records


@pytest.fixture()
def database(tmp_path):
    fasta = tmp_path / "genes.fa"
    records = _genes(fasta)
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(tmp_path / "genes"),
            "-parse_seqids",
        ],
        check=True,
        capture_output=True,
    )
    query = tmp_path / "q.fa"
    query.write_text(f">q0\n{records[2][100:300]}\n")
    return tmp_path, records, query


def test_missing_program_reports_install_hint(tmp_path):
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGTACGTACGTACGTACGT\n")
    (tmp_path / "db.nin").write_text("x")
    with pytest.raises(BackendUnavailableError):
        blast.blastn_search(
            str(query),
            str(tmp_path / "o.tsv"),
            database=str(tmp_path / "db"),
            blast_path=str(tmp_path),
        )


def test_argument_validation(tmp_path):
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    with pytest.raises(ConfigurationError):
        blast.blastn_search(
            str(query), str(tmp_path / "o"), task="blastn"
        )  # neither database nor subject
    with pytest.raises(ConfigurationError):
        blast.blastn_search(str(query), str(tmp_path / "o"), subject=str(query), task="nope")
    with pytest.raises(ConfigurationError):
        blast.blastn_search(str(query), str(tmp_path / "o"), subject=str(query), outfmt="6; rm")
    with pytest.raises(ConfigurationError):
        blast.blastdbcmd_fetch_sequences(
            str(tmp_path / "missing"), str(tmp_path / "o"), entries=["a"]
        )


@needs_blast
@pytest.mark.parametrize("task", ["blastn", "megablast", "dc-megablast"])
def test_blastn_tasks_find_the_source_gene(database, task):
    folder, _, query = database
    result = blast.blastn_search(
        str(query),
        str(folder / f"{task}.tsv"),
        database=str(folder / "genes"),
        task=task,
        outfmt=6,
        verbose=False,
    )
    assert result["status"] == "completed" and result["version"]
    assert result["rows"][0]["sseqid"] == "cds2"
    assert float(result["rows"][0]["pident"]) == 100.0


@needs_blast
def test_blastn_short_and_subject_mode(database):
    folder, records, _ = database
    short = folder / "primer.fa"
    short.write_text(f">p\n{records[1][50:72]}\n")
    result = blast.blastn_search(
        str(short),
        str(folder / "short.tsv"),
        subject=str(folder / "genes.fa"),
        task="blastn-short",
        outfmt="6 qseqid sseqid length",
        verbose=False,
    )
    assert result["rows"][0] == {"qseqid": "p", "sseqid": "cds1", "length": "22"}


@needs_blast
def test_blastdbcmd_range_and_batch(database):
    folder, records, _ = database
    one = blast.blastdbcmd_fetch_sequences(
        str(folder / "genes"),
        str(folder / "one.fa"),
        entries=["cds3"],
        sequence_range="11-40",
        verbose=False,
    )
    assert one["status"] == "completed"
    assert (folder / "one.fa").read_text().split("\n")[1] == records[3][10:40]
    (folder / "ids.txt").write_text("cds0\ncds4\n")
    blast.blastdbcmd_fetch_sequences(
        str(folder / "genes"),
        str(folder / "b.txt"),
        entry_batch=str(folder / "ids.txt"),
        outfmt="%a %l",
        verbose=False,
    )
    assert (folder / "b.txt").read_text().split("\n")[:2] == ["cds0 900", "cds4 900"]


@needs_blast
def test_makembindex_and_use_in_search(database):
    folder, _, _query = database
    if shutil.which("makembindex") is None:
        pytest.skip("makembindex is not installed")
    result = blast.makembindex_build_index(str(folder / "genes"), verbose=False)
    assert result["status"] == "completed" and result["outputs"][0].endswith(".00.idx")
    with pytest.raises(ConfigurationError):
        blast.makembindex_build_index(str(folder / "genes"), verbose=False)


@needs_blast
def test_windowmasker_pipeline_and_mask_conversion(tmp_path):
    if shutil.which("windowmasker") is None or shutil.which("convert2blastmask") is None:
        pytest.skip("windowmasker is not installed")
    rng = random.Random(4)
    genome = (
        "".join(rng.choice("ACGT") for _ in range(20000))
        + "AC" * 400
        + "".join(rng.choice("ACGT") for _ in range(20000))
    )
    fasta = tmp_path / "g.fa"
    fasta.write_text(
        ">chr1\n" + "\n".join(genome[i : i + 70] for i in range(0, len(genome), 70)) + "\n"
    )
    counts = blast.windowmasker_count_units(str(fasta), str(tmp_path / "c.txt"), verbose=False)
    assert counts["status"] == "completed"
    binary = blast.windowmasker_convert_counts(
        str(tmp_path / "c.txt"), str(tmp_path / "c.obin"), counts_format="obinary", verbose=False
    )
    assert binary["status"] == "completed"
    intervals = blast.windowmasker_mask_sequences(
        str(fasta), str(tmp_path / "c.obin"), str(tmp_path / "m.int"), verbose=False
    )
    soft = blast.windowmasker_mask_sequences(
        str(fasta),
        str(tmp_path / "c.txt"),
        str(tmp_path / "m.fa"),
        output_format="fasta",
        verbose=False,
    )
    assert intervals["status"] == soft["status"] == "completed"
    text = (tmp_path / "m.fa").read_text()
    assert any(char.islower() for char in text)
    mask = blast.convert2blastmask_from_fasta(
        str(tmp_path / "m.fa"),
        str(tmp_path / "mask.asn"),
        masking_algorithm="windowmasker",
        masking_options="-mk_counts",
        verbose=False,
    )
    assert mask["status"] == "completed"


def test_vdb_and_rps_validation(tmp_path):
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    with pytest.raises(ConfigurationError):
        blast.blast_vdb_search(str(query), [], str(tmp_path / "o"))
    with pytest.raises(ConfigurationError):
        blast.blast_vdb_search(
            str(query), ["SRR1"], str(tmp_path / "o"), program="blastn_vdb", genetic_code=1
        )
    with pytest.raises(ConfigurationError):
        blast.blast_vdb_cmd_fetch_reads(["SRR1"], str(tmp_path / "o"), view="sequences")
    with pytest.raises(ConfigurationError):
        blast.rpstblastn_search(str(query), str(tmp_path / "nodb"), str(tmp_path / "o"))


def test_vdb_programs_receive_the_expected_arguments(tmp_path):
    """A stand-in executable records its argument vector; no network is used."""
    folder = tmp_path / "bin"
    folder.mkdir()
    for name in ("blastn_vdb", "tblastn_vdb", "blast_vdb_cmd"):
        script = folder / name
        script.write_text(
            '#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a" >> "$0.log"; done\nwhile [ $# -gt 0 ]; do if [ "$1" = "-out" ]; then : > "$2"; fi; shift; done\n'
        )
        script.chmod(0o755)
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    result = blast.blast_vdb_search(
        str(query),
        ["SRR1", "SRR2"],
        str(tmp_path / "v.tsv"),
        program="tblastn_vdb",
        genetic_code=11,
        evalue=1e-3,
        blast_path=str(folder),
        verbose=False,
    )
    assert result["status"] == "completed"
    arguments = (folder / "tblastn_vdb.log").read_text().split("\n")
    assert arguments[arguments.index("-db") + 1] == "SRR1 SRR2"
    assert "-db_gencode" in arguments and "-evalue" in arguments
    fetched = blast.blast_vdb_cmd_fetch_reads(
        ["SRR1"], str(tmp_path / "f.txt"), entries=["all"], blast_path=str(folder), verbose=False
    )
    assert (
        fetched["status"] == "completed" and "-entry" in (folder / "blast_vdb_cmd.log").read_text()
    )
