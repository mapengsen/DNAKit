"""Tests for dnakit.databases.gget_gap_functions (gget downloads, cBioPortal, NCBI Virus and external alignments).

Network functions replace one gget upstream call with a recording fake, so no request is sent; the tests check
argument mapping, validation (the fake must not be called for invalid input) and the result document. The two
alignment programs run for real: MUSCLE 5 and DIAMOND are taken from gget's bundled Linux builds when gget is
installed, or from PATH, and the tests skip cleanly otherwise. One live check runs only with DNAKIT_LIVE_TESTS=1.
"""

from __future__ import annotations

import importlib
import inspect
import itertools
import os
import shutil
import stat
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from dnakit.databases import gget_gap_functions as gg
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    BackendVersionError,
    ConfigurationError,
)

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENES = DATA / "genes.fa"


def _gget_package() -> Any:
    return pytest.importorskip("gget", reason="gget is not installed")


def _bundled(name: str) -> str | None:
    """A gget-bundled Linux binary when gget is installed, else the binary on PATH."""
    try:
        package = importlib.import_module("gget")
        candidate = Path(package.__file__).parent / "bins" / "Linux" / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    except ImportError:
        pass
    return shutil.which(name)


def _fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return {key: "".join(parts) for key, parts in records.items()}


class Recorder:
    """A fake upstream call that records its arguments and returns a canned value."""

    def __init__(self, value: Any = None, *, effect: Callable[..., Any] | None = None) -> None:
        self.value = value
        self.effect = effect
        self.calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append((args, kwargs))
        if self.effect is not None:
            return self.effect(*args, **kwargs)
        return self.value


# --- module contract ---------------------------------------------------------------------------


def test_public_functions_have_keyword_only_verbose_and_docstrings() -> None:
    assert len(gg.__all__) == len(set(gg.__all__))
    for name in gg.__all__:
        function = getattr(gg, name)
        assert callable(function), name
        signature = inspect.signature(function)
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert function.__doc__, name
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, f"{name}.{parameter.name}"


# --- cBioPortal ---------------------------------------------------------------------------------


def test_cbio_plot_maps_arguments_and_reports_the_figure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    gget = _gget_package()

    def draw(studies: list[str], genes: list[str], **kwargs: Any) -> bool:
        Path(kwargs["figure_dir"], kwargs["figure_filename"]).write_bytes(b"\x89PNG")
        return True

    fake = Recorder(effect=draw)
    monkeypatch.setattr(gget, "cbio_plot", fake)
    record = gg.gget_cbio_plot(
        ["brca_tcga"],
        ["TP53", "BRCA1"],
        figure_dir=str(tmp_path / "figs"),
        figure_filename="heatmap.png",
        data_dir=str(tmp_path / "data"),
        verbose=False,
    )
    args, kwargs = fake.calls[0]
    assert args == (["brca_tcga"], ["TP53", "BRCA1"])
    assert (
        kwargs["stratification"] == "tissue"
        and kwargs["show"] is False
        and kwargs["confirm_download"] is False
    )
    assert record["figure"] == str(tmp_path / "figs" / "heatmap.png")
    assert record["service"].startswith("cBioPortal")


def test_cbio_plot_failure_and_bad_names_are_reported(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    gget = _gget_package()
    monkeypatch.setattr(gget, "cbio_plot", Recorder(value=False))
    with pytest.raises(BackendExecutionError):
        gg.gget_cbio_plot(
            ["brca_tcga"],
            ["TP53"],
            figure_dir=str(tmp_path / "f"),
            figure_filename="x.png",
            data_dir=str(tmp_path / "d"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gg.gget_cbio_plot(
            ["brca_tcga"],
            ["TP53"],
            figure_dir=str(tmp_path / "f"),
            figure_filename="../x.png",
            data_dir=str(tmp_path / "d"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gg.gget_cbio_plot(
            [],
            ["TP53"],
            figure_dir=str(tmp_path),
            figure_filename="x.png",
            data_dir=str(tmp_path),
            verbose=False,
        )


def test_cbio_download_writes_only_to_the_named_folder(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _gget_package()
    module = importlib.import_module("gget.gget_cbio")
    fake = Recorder(value=True)
    monkeypatch.setattr(module, "download_cbioportal_data", fake)
    record = gg.gget_cbio_download(["brca_tcga"], out_dir=str(tmp_path / "study"), verbose=False)
    assert fake.calls[0][0] == (["brca_tcga"],)
    assert fake.calls[0][1]["out_dir"] == str(tmp_path / "study")
    assert record["out_dir"] == str(tmp_path / "study")


def test_cbio_download_reports_incomplete_download(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _gget_package()
    module = importlib.import_module("gget.gget_cbio")
    monkeypatch.setattr(module, "download_cbioportal_data", Recorder(value=False))
    with pytest.raises(BackendExecutionError):
        gg.gget_cbio_download(["brca_tcga"], out_dir=str(tmp_path / "study"), verbose=False)


# --- Ensembl references -------------------------------------------------------------------------


def test_ref_maps_species_release_and_subset(monkeypatch: pytest.MonkeyPatch) -> None:
    gget = _gget_package()
    fake = Recorder(value={"gtf": {"url": "https://example.invalid/x.gtf.gz"}})
    monkeypatch.setattr(gget, "ref", fake)
    record = gg.gget_ref("homo_sapiens", which=["gtf", "dna"], release=112, ftp=True, verbose=False)
    args, kwargs = fake.calls[0]
    assert args == ("homo_sapiens",)
    assert (
        kwargs["which"] == ["gtf", "dna"] and kwargs["release"] == 112 and kwargs["save"] is False
    )
    assert record["links"]["gtf"]["url"].startswith("https://")


def test_ref_refuses_unknown_kinds_without_calling_upstream(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gget = _gget_package()
    fake = Recorder(value={})
    monkeypatch.setattr(gget, "ref", fake)
    with pytest.raises(ConfigurationError):
        gg.gget_ref("homo_sapiens", which=["gtf", "bam"], verbose=False)
    with pytest.raises(ConfigurationError):
        gg.gget_ref("Homo Sapiens", verbose=False)
    with pytest.raises(ConfigurationError):
        gg.gget_ref("homo_sapiens", release=0, verbose=False)
    assert fake.calls == []


def test_upstream_exception_is_wrapped_as_backend_error(monkeypatch: pytest.MonkeyPatch) -> None:
    gget = _gget_package()

    def broken(*args: Any, **kwargs: Any) -> Any:
        raise OSError("connection reset")

    monkeypatch.setattr(gget, "ref", broken)
    with pytest.raises(BackendExecutionError, match="connection reset"):
        gg.gget_ref("homo_sapiens", verbose=False)


# --- NCBI Virus and GenBank ----------------------------------------------------------------------


def test_virus_download_passes_filters_and_needs_an_output_folder(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    gget = _gget_package()
    fake = Recorder(value=str(tmp_path / "out"))
    monkeypatch.setattr(gget, "virus", fake)
    record = gg.gget_virus(
        "Zika virus",
        outfolder=str(tmp_path / "out"),
        min_seq_length=9000,
        max_seq_length=12000,
        annotated=True,
        min_release_date="2020-01-01",
        verbose=False,
    )
    args, kwargs = fake.calls[0]
    assert args == ("Zika virus",)
    assert (
        kwargs["min_seq_length"] == 9000
        and kwargs["annotated"] is True
        and kwargs["outfolder"] == str(tmp_path / "out")
    )
    assert kwargs["min_release_date"] == "2020-01-01" and kwargs["max_release_date"] is None
    assert record["virus"] == "Zika virus"
    with pytest.raises(ConfigurationError):
        gg.gget_virus(
            "Zika virus",
            outfolder=str(tmp_path / "out"),
            min_seq_length=10,
            max_seq_length=5,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gg.gget_virus("Zika; rm -rf", outfolder=str(tmp_path / "out2"), verbose=False)
    with pytest.raises(ConfigurationError):
        gg.gget_virus(
            "Zika virus",
            outfolder=str(tmp_path / "out"),
            min_release_date="2020/01/01",
            verbose=False,
        )


def test_virus_metadata_and_chunked_metadata_map_arguments(monkeypatch: pytest.MonkeyPatch) -> None:
    _gget_package()
    module = importlib.import_module("gget.gget_virus")
    rows = [{"accession": "MN908947.3", "length": 29903}]
    single = Recorder(value=rows)
    chunked = Recorder(value=rows)
    monkeypatch.setattr(module, "fetch_virus_metadata", single)
    monkeypatch.setattr(module, "fetch_virus_metadata_chunked", chunked)
    record = gg.gget_virus_metadata(
        "Severe acute respiratory syndrome coronavirus 2", complete_only=True, verbose=False
    )
    assert record["count"] == 1 and record["records"][0]["accession"] == "MN908947.3"
    assert single.calls[0][1]["complete_only"] is True and single.calls[0][1]["refseq_only"] is None
    chunk = gg.gget_virus_metadata_chunked(
        "Zika virus", max_release_date="2024-01-01", verbose=False
    )
    assert chunk["count"] == 1
    assert (
        chunked.calls[0][1]["max_release_date"] == "2024-01-01"
        and chunked.calls[0][1]["complete_only"] is False
    )


def test_sars_cov2_and_influenza_downloads_return_the_package_and_filters(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _gget_package()
    module = importlib.import_module("gget.gget_virus")
    zip_path = str(tmp_path / "sars.zip")
    monkeypatch.setattr(
        module, "download_sars_cov2_optimized", Recorder(value=(zip_path, ["host"], ["lineage"]))
    )
    monkeypatch.setattr(
        module, "download_alphainfluenza_optimized", Recorder(value=(zip_path, [], []))
    )
    record = gg.gget_sars_cov2_download(
        outdir=str(tmp_path / "sars"), lineage="BA.5", verbose=False
    )
    assert (
        record["zip_path"] == zip_path
        and record["applied_filters"] == ["host"]
        and record["missing_filters"] == ["lineage"]
    )
    flu = gg.gget_alphainfluenza_download(outdir=str(tmp_path / "flu"), verbose=False)
    assert flu["applied_filters"] == []


def test_sequences_by_accession_and_genbank_metadata_refuse_existing_outputs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _gget_package()
    module = importlib.import_module("gget.gget_virus")
    fasta = str(tmp_path / "seqs.fasta")
    monkeypatch.setattr(module, "download_sequences_by_accessions", Recorder(value=fasta))
    monkeypatch.setattr(module, "fetch_genbank_metadata", Recorder(value={"accessions": 2}))
    record = gg.gget_virus_sequences_by_accession(
        ["MN908947.3", "OK091006.1"], outdir=str(tmp_path / "seqs"), api_key="secret", verbose=False
    )
    assert record["fasta_path"] == fasta and record["accessions"] == 2
    assert "secret" not in str(record)
    xml, csv = tmp_path / "meta.xml", tmp_path / "meta.csv"
    summary = gg.gget_genbank_metadata(
        ["MN908947.3"], xml_path=str(xml), csv_path=str(csv), verbose=False
    )
    assert summary["summary"] == {"accessions": 2}
    xml.write_text("<x/>", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        gg.gget_genbank_metadata(
            ["MN908947.3"], xml_path=str(xml), csv_path=str(tmp_path / "other.csv"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        gg.gget_virus_sequences_by_accession(
            ["bad id with space"], outdir=str(tmp_path / "s2"), verbose=False
        )


def test_identifier_lists_must_be_lists_of_plain_identifiers() -> None:
    with pytest.raises(ConfigurationError):
        gg.gget_virus_sequences_by_accession("MN908947.3", outdir="/tmp/unused", verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        gg.gget_virus_sequences_by_accession([], outdir="/tmp/unused", verbose=False)


# --- external programs: MUSCLE 5 ----------------------------------------------------------------


@pytest.fixture()
def muscle_input(tmp_path: Path) -> Path:
    if not GENES.exists():
        pytest.skip("shared consistency genes.fa not present")
    records = _fasta_records(GENES)
    selected = list(records.items())[:3]
    path = tmp_path / "input.fa"
    path.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in selected), encoding="utf-8"
    )
    return path


def test_muscle_alignment_keeps_every_sequence(muscle_input: Path, tmp_path: Path) -> None:
    program = _bundled("muscle")
    if program is None:
        pytest.skip(
            "MUSCLE is not installed (gget bundles it on Linux; otherwise install MUSCLE 5)"
        )
    output = tmp_path / "aligned.afa"
    record = gg.gget_muscle_align(
        str(muscle_input), str(output), muscle_path=program, verbose=False
    )
    assert record["records"] == 3
    assert record["program"]["version"].split()[0] == "muscle"
    aligned = _fasta_records(output)
    original = _fasta_records(muscle_input)
    assert set(aligned) == set(original)
    for name, sequence in original.items():
        assert aligned[name].replace("-", "") == sequence.upper()
    assert record["aligned_length"] == len(next(iter(aligned.values())))


def test_muscle_matches_gget_muscle_output(muscle_input: Path, tmp_path: Path) -> None:
    gget = _gget_package()
    program = _bundled("muscle")
    if program is None:
        pytest.skip("MUSCLE is not installed")
    ours = tmp_path / "ours.afa"
    gg.gget_muscle_align(str(muscle_input), str(ours), muscle_path=program, verbose=False)
    theirs = tmp_path / "theirs.afa"
    gget.muscle(str(muscle_input), out=str(theirs), verbose=False)
    assert ours.read_bytes() == theirs.read_bytes()


def test_muscle_super5_runs_the_super5_algorithm(muscle_input: Path, tmp_path: Path) -> None:
    program = _bundled("muscle")
    if program is None:
        pytest.skip("MUSCLE is not installed")
    output = tmp_path / "super5.afa"
    record = gg.gget_muscle_align(
        str(muscle_input), str(output), super5=True, muscle_path=program, verbose=False
    )
    assert record["command"][1] == "-super5"
    assert record["records"] == 3


def test_missing_muscle_reports_an_install_hint(
    muscle_input: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(gg.shutil, "which", lambda name: None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        gg.gget_muscle_align(str(muscle_input), str(tmp_path / "x.afa"), verbose=False)
    assert "muscle" in str(excinfo.value).lower()


def _fake_program(directory: Path, name: str, body: str) -> str:
    path = directory / name
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return str(path)


def test_muscle_version_three_is_refused(muscle_input: Path, tmp_path: Path) -> None:
    old = _fake_program(tmp_path, "muscle", 'echo "muscle v3.8.1551"\n')
    with pytest.raises(BackendVersionError):
        gg.gget_muscle_align(
            str(muscle_input), str(tmp_path / "x.afa"), muscle_path=old, verbose=False
        )


def test_missing_explicit_program_path_is_unavailable(muscle_input: Path, tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        gg.gget_muscle_align(
            str(muscle_input),
            str(tmp_path / "x.afa"),
            muscle_path=str(tmp_path / "nope"),
            verbose=False,
        )


def test_existing_output_is_refused(muscle_input: Path, tmp_path: Path) -> None:
    output = tmp_path / "exists.afa"
    output.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        gg.gget_muscle_align(
            str(muscle_input), str(output), muscle_path=str(muscle_input), verbose=False
        )


# --- external programs: DIAMOND ------------------------------------------------------------------


@pytest.fixture()
def protein_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Protein query and reference: the reference holds three proteins and the query one of them with a change."""
    reference_seqs = {
        "ref_a": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAQ",
        "ref_b": "MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGIVHLTLVFFVHTEKLLPCIG",
        "ref_c": "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK",
    }
    reference = tmp_path / "reference.faa"
    reference.write_text(
        "".join(f">{k}\n{v}\n" for k, v in reference_seqs.items()), encoding="utf-8"
    )
    query_seq = reference_seqs["ref_c"][:40] + "W" + reference_seqs["ref_c"][41:]
    query = tmp_path / "query.faa"
    query.write_text(f">query_c\n{query_seq}\n", encoding="utf-8")
    return query, reference


def test_diamond_protein_search_finds_the_matching_reference(
    protein_pair: tuple[Path, Path], tmp_path: Path
) -> None:
    program = _bundled("diamond")
    if program is None:
        pytest.skip(
            "DIAMOND is not installed (gget bundles it on Linux; otherwise install DIAMOND)"
        )
    query, reference = protein_pair
    record = gg.gget_diamond_align(
        str(query), str(reference), str(tmp_path / "hits.tsv"), diamond_path=program, verbose=False
    )
    assert record["mode"] == "blastp"
    assert record["columns"][:3] == ["qseqid", "sseqid", "pident"]
    best = max(record["rows"], key=lambda row: row["bitscore"])
    assert best["qseqid"] == "query_c" and best["sseqid"] == "ref_c"
    assert 95.0 < best["pident"] < 100.0
    assert best["mismatch"] >= 1
    assert not list(tmp_path.glob("dnakit-diamond-*")), "the temporary database folder is removed"


def test_diamond_matches_gget_diamond_hits(protein_pair: tuple[Path, Path], tmp_path: Path) -> None:
    gget = _gget_package()
    program = _bundled("diamond")
    if program is None:
        pytest.skip("DIAMOND is not installed")
    query, reference = protein_pair
    ours = gg.gget_diamond_align(
        str(query), str(reference), str(tmp_path / "ours.tsv"), diamond_path=program, verbose=False
    )
    frame = gget.diamond(str(query), str(reference), out=str(tmp_path / "gget_out"), verbose=False)
    # gget renames the outfmt-6 columns; map them back to DIAMOND's names before comparing.
    renamed = {
        "query_accession": "qseqid",
        "subject_accession": "sseqid",
        "identity_percentage": "pident",
        "mismatches": "mismatch",
        "bit_score": "bitscore",
        "query_start": "qstart",
        "query_end": "qend",
        "subject_start": "sstart",
        "subject_end": "send",
    }
    theirs = [
        {renamed.get(key, key): value for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]
    keys = (
        "qseqid",
        "sseqid",
        "qstart",
        "qend",
        "sstart",
        "send",
        "pident",
        "mismatch",
        "bitscore",
    )
    ours_rows = sorted(
        tuple(float(row[key]) if key not in ("qseqid", "sseqid") else row[key] for key in keys)
        for row in ours["rows"]
    )
    theirs_rows = sorted(
        tuple(float(row[key]) if key not in ("qseqid", "sseqid") else str(row[key]) for key in keys)
        for row in theirs
    )
    assert ours_rows == theirs_rows
    assert ours_rows, "the query has at least one hit"


def test_diamond_translated_mode_uses_blastx(tmp_path: Path) -> None:
    program = _bundled("diamond")
    if program is None:
        pytest.skip("DIAMOND is not installed")
    if not GENES.exists():
        pytest.skip("shared consistency genes.fa not present")
    dna = next(iter(_fasta_records(GENES).values()))
    codons = {
        "".join(bases): amino
        for bases, amino in zip(
            itertools.product("TCAG", repeat=3),
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            strict=True,
        )
    }
    protein = "".join(codons.get(dna[i : i + 3], "X") for i in range(0, len(dna) - 2, 3)).replace(
        "*", "X"
    )
    reference = tmp_path / "protein.faa"
    reference.write_text(f">translated\n{protein}\n", encoding="utf-8")
    query = tmp_path / "dna.fa"
    query.write_text(f">dna\n{dna}\n", encoding="utf-8")
    record = gg.gget_diamond_align(
        str(query),
        str(reference),
        str(tmp_path / "blastx.tsv"),
        translated=True,
        diamond_path=program,
        verbose=False,
    )
    assert record["mode"] == "blastx"
    assert record["hits"] >= 1
    assert max(row["pident"] for row in record["rows"]) == 100.0


def test_diamond_rejects_bad_sensitivity_and_threads(
    protein_pair: tuple[Path, Path], tmp_path: Path
) -> None:
    query, reference = protein_pair
    with pytest.raises(ConfigurationError):
        gg.gget_diamond_align(
            str(query), str(reference), str(tmp_path / "x.tsv"), sensitivity="warp", verbose=False
        )
    with pytest.raises(ConfigurationError):
        gg.gget_diamond_align(
            str(query), str(reference), str(tmp_path / "x.tsv"), threads=0, verbose=False
        )


def test_missing_diamond_reports_an_install_hint(
    protein_pair: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query, reference = protein_pair
    monkeypatch.setattr(gg.shutil, "which", lambda name: None)
    with pytest.raises(BackendUnavailableError):
        gg.gget_diamond_align(str(query), str(reference), str(tmp_path / "x.tsv"), verbose=False)


# --- live (opt-in) -----------------------------------------------------------------------------------


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live check: set DNAKIT_LIVE_TESTS=1"
)
def test_live_ensembl_reference_links_for_human_gtf() -> None:
    _gget_package()
    record = gg.gget_ref("homo_sapiens", which="gtf", ftp=True, verbose=False)
    assert record["links"]
