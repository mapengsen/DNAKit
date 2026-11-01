"""Tests for dnakit.databases.genomepy_gap_functions (genomepy wrappers added for Fig. 2 coverage).

No test touches the network: genomepy provider functions are replaced by monkeypatch, and the UCSC
tools check is replaced where a conversion would otherwise need them. Set DNAKIT_LIVE_TESTS=1 to
run the one live catalogue query.
"""

from __future__ import annotations

import inspect
import os
import stat
from pathlib import Path

import pytest

pytest.importorskip("genomepy")

from dnakit.databases import genomepy_gap_functions as gap
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


@pytest.fixture()
def no_ucsc_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.annotation.utils as annotation

    def missing(tools: list | None = None) -> None:
        raise OSError("gtfToGenePred is not installed")

    monkeypatch.setattr(annotation, "check_ucsc_tools", missing)


@pytest.fixture()
def ucsc_tools_present(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.annotation.utils as annotation

    monkeypatch.setattr(annotation, "check_ucsc_tools", lambda tools=None: None)


def _fake_program(path: Path, name: str) -> Path:
    """A POSIX shell program that records its arguments and writes one index file."""
    path.mkdir(parents=True, exist_ok=True)
    program = path / name
    program.write_text('#!/bin/sh\necho "$@" > index.out\n', encoding="utf-8")
    program.chmod(program.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return program


def test_public_functions_are_verbose_keyword_only() -> None:
    assert len(gap.__all__) == 23
    for name in gap.__all__:
        parameter = inspect.signature(getattr(gap, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_list_available_genomes_forwards_provider_and_size(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.functions as functions

    seen: dict = {}

    def fake(provider=None, size=False):
        seen.update(provider=provider, size=size)
        yield ["hg38", "UCSC", "Human"]

    monkeypatch.setattr(functions, "list_available_genomes", fake)
    result = gap.genomepy_list_available_genomes(provider="ucsc", include_size=True)
    assert seen == {"provider": "ucsc", "size": True}
    assert result["rows"] == [["hg38", "UCSC", "Human"]] and result["row_count"] == 1


def test_search_genomes_uses_provider_search(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.providers as providers

    monkeypatch.setattr(providers, "search", lambda term, provider, exact, size: iter([[term]]))
    result = gap.genomepy_search_genomes("GRCh38", exact=True)
    assert result["rows"] == [["GRCh38"]]


def test_search_genomes_rejects_empty_term() -> None:
    with pytest.raises(ConfigurationError):
        gap.genomepy_search_genomes("  ")


def test_head_annotations_returns_printed_lines(
    monkeypatch: pytest.MonkeyPatch, ucsc_tools_present
) -> None:
    import genomepy.functions as functions

    def fake(name, provider=None, n=2):
        print("chr1\tgene A")
        print("chr1\tgene B")

    monkeypatch.setattr(functions, "head_annotations", fake)
    result = gap.genomepy_head_annotations("hg38", lines=2)
    assert result["lines"] == ["chr1\tgene A", "chr1\tgene B"]


def test_head_annotations_requires_ucsc_tools(
    monkeypatch: pytest.MonkeyPatch, no_ucsc_tools
) -> None:
    with pytest.raises(BackendUnavailableError, match="UCSC"):
        gap.genomepy_head_annotations("hg38")


def test_generate_annotation_needs_ucsc_tools(tmp_path: Path, no_ucsc_tools) -> None:
    template = tmp_path / "genes.gtf"
    template.write_text('chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "g";\n', encoding="utf-8")
    with pytest.raises(BackendUnavailableError) as info:
        gap.genomepy_generate_annotation(str(template), str(tmp_path / "out.genePred"))
    assert "ucsc" in info.value.hint.lower()


def test_generate_annotation_refuses_to_overwrite(tmp_path: Path, ucsc_tools_present) -> None:
    template = tmp_path / "genes.bed"
    template.write_text("chr1\t0\t10\tg\t0\t+\n", encoding="utf-8")
    target = tmp_path / "exists.genePred"
    target.write_text("old", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="overwrite"):
        gap.genomepy_generate_annotation(str(template), str(target))


def test_generate_annotation_writes_target_when_tools_present(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, ucsc_tools_present
) -> None:
    import genomepy.annotation.utils as annotation

    def fake_generate(template, target, overwrite=False):
        Path(target).write_text(Path(template).read_text(), encoding="utf-8")

    monkeypatch.setattr(annotation, "generate_annot", fake_generate)
    template = tmp_path / "genes.bed"
    template.write_text("chr1\t0\t10\tg\t0\t+\n", encoding="utf-8")
    result = gap.genomepy_generate_annotation(str(template), str(tmp_path / "out.genePred"))
    assert Path(result["outputs"]["annotation"]).read_text() == template.read_text()


def test_download_assembly_report_returns_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.providers.ncbi as ncbi
    import pandas as pd

    frame = pd.DataFrame({"Sequence-Name": ["1", "2"], "Sequence-Length": ["100", "200"]})
    monkeypatch.setattr(
        ncbi, "download_assembly_report", lambda acc, fname=None, quiet=False: frame
    )
    result = gap.genomepy_download_assembly_report("GCF_000001405.40")
    assert result["rows"][1] == {"Sequence-Name": "2", "Sequence-Length": "200"}


def test_download_assembly_report_writes_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import genomepy.providers.ncbi as ncbi

    def fake(acc, fname=None, quiet=False):
        Path(fname).write_text("Sequence-Name\tSequence-Length\n", encoding="utf-8")

    monkeypatch.setattr(ncbi, "download_assembly_report", fake)
    result = gap.genomepy_download_assembly_report(
        "GCA_000001405.1", output_path=str(tmp_path / "report.tsv")
    )
    assert Path(result["outputs"]["report"]).is_file()


def test_download_assembly_report_rejects_bad_accession_type() -> None:
    with pytest.raises(ConfigurationError):
        gap.genomepy_download_assembly_report(None)  # type: ignore[arg-type]


@pytest.mark.parametrize(("value", "expected"), [(3_100_000_000, 3_100_000_000), (-1, None)])
def test_genome_size_maps_missing_reports_to_none(
    monkeypatch: pytest.MonkeyPatch, value: int, expected: int | None
) -> None:
    import genomepy.providers.ncbi as ncbi

    monkeypatch.setattr(ncbi, "get_genome_size", lambda acc: value)
    assert gap.genomepy_genome_size("GCF_000001405.40")["genome_size"] == expected


def test_catalogue_wrappers_return_provider_dicts(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.providers.ensembl as ensembl
    import genomepy.providers.gencode as gencode
    import genomepy.providers.ncbi as ncbi
    import genomepy.providers.ucsc as ucsc

    catalogue = {"GRCh38": {"provider": "ncbi", "accession": "GCF_000001405.40"}}
    monkeypatch.setattr(ncbi, "get_genomes", lambda url: catalogue)
    monkeypatch.setattr(ensembl, "get_genomes", lambda url: catalogue)
    monkeypatch.setattr(ucsc, "get_genomes", lambda url: catalogue)
    monkeypatch.setattr(gencode, "get_genomes", lambda link: catalogue)
    for call in (
        gap.genomepy_ncbi_genomes,
        gap.genomepy_ensembl_genomes,
        gap.genomepy_ucsc_genomes,
        gap.genomepy_gencode_genomes,
    ):
        result = call("https://example.invalid/catalogue")
        assert result["genomes"] == catalogue and result["genome_count"] == 1


def test_gencode_releases_is_a_pure_listing_function() -> None:
    listing = ["release_38", "release_39", "release_21", "README"]
    assert gap.genomepy_gencode_releases(listing, "human")["releases"] == ["39", "38"]
    assert gap.genomepy_gencode_releases(listing, "mouse")["releases"] == ["M39", "M38"]


def test_gencode_releases_validates_species() -> None:
    with pytest.raises(ConfigurationError):
        gap.genomepy_gencode_releases(["release_38"], "zebrafish")


def test_gencode_to_ucsc_names_forwards_catalogue(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.providers.gencode as gencode

    monkeypatch.setattr(gencode, "get_gencode2ucsc", lambda genomes: {"chr1": "chr1"})
    result = gap.genomepy_gencode_to_ucsc_names({"human": {"ftp": "x"}})
    assert result["mapping"] == {"chr1": "chr1"}


def test_download_annotation_creates_named_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, ucsc_tools_present
) -> None:
    import genomepy.providers.base as base

    def fake(genomes_dir, annot_url, localname, n=None):
        directory = Path(genomes_dir) / localname
        directory.mkdir(parents=True)
        (directory / f"{localname}.annotation.gtf").write_text("chr1\tx\n", encoding="utf-8")

    monkeypatch.setattr(base, "download_annotation", fake)
    result = gap.genomepy_download_annotation(
        str(tmp_path), "https://example.invalid/a.gtf.gz", "hg38"
    )
    assert Path(result["outputs"]["directory"]).name == "hg38"
    assert any(value.endswith(".annotation.gtf") for value in result["outputs"].values())


def test_ucsc_download_annotation_uses_mysql_provider(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, ucsc_tools_present
) -> None:
    import genomepy.providers.ucsc as ucsc

    def fake(name, annot, genomes_dir, localname, n=None):
        (Path(genomes_dir) / localname).mkdir(parents=True)
        (Path(genomes_dir) / localname / "ucsc.genePred").write_text("x\n", encoding="utf-8")

    monkeypatch.setattr(ucsc, "download_annotation", fake)
    result = gap.genomepy_ucsc_download_annotation("hg38", "ncbiRefSeq", str(tmp_path), "hg38_ucsc")
    assert Path(result["outputs"]["directory"]).name == "hg38_ucsc"


def test_download_annotation_head_writes_first_lines(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import genomepy.providers.base as base

    def fake(annot_url, annot_file, n=5):
        Path(annot_file).write_text("".join(f"line{i}\n" for i in range(n)), encoding="utf-8")

    monkeypatch.setattr(base, "download_head", fake)
    result = gap.genomepy_download_annotation_head(
        "https://example.invalid/a.gtf.gz", str(tmp_path / "head.gtf"), n=3
    )
    assert Path(result["outputs"]["head"]).read_text().count("\n") == 3


def test_find_annotation_urls_forwards_url(monkeypatch: pytest.MonkeyPatch) -> None:
    import genomepy.providers.url as url_provider

    monkeypatch.setattr(
        url_provider,
        "search_url_for_annotations",
        lambda url: [f"{url}/genes.gtf.gz"],
    )
    result = gap.genomepy_find_annotation_urls("https://example.invalid/genome")
    assert result["urls"] == ["https://example.invalid/genome/genes.gtf.gz"]


def test_index_command_runs_only_the_aligner_program(tmp_path: Path) -> None:
    fake = _fake_program(tmp_path / "bin", "bwa")
    workdir = tmp_path / "work"
    workdir.mkdir()
    result = gap.genomepy_bwa_index_command(
        ["bwa", "index", "-p", "genes", "genes.fa"],
        working_directory=str(workdir),
        executable_path=str(fake),
    )
    created = Path(result["outputs"]["index_0"])
    assert created.name == "index.out"
    assert "index -p genes genes.fa" in created.read_text()
    assert result["backend"] == "bwa"


def test_index_command_refuses_other_executables(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="bowtie2-build"):
        gap.genomepy_bowtie2_index_command(["rm", "-rf", "x"], working_directory=str(tmp_path))


def test_index_command_refuses_shell_strings(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="no shell"):
        gap.genomepy_bowtie2_index_command(
            "bowtie2-build genes.fa genes",
            working_directory=str(tmp_path),  # type: ignore[arg-type]
        )


def test_index_command_needs_existing_working_directory(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="working_directory"):
        gap.genomepy_star_index_command(
            ["STAR", "--runMode", "genomeGenerate"], working_directory=str(tmp_path / "missing")
        )


def test_missing_index_program_gives_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(gap.shutil, "which", lambda name: None)
    with pytest.raises(BackendUnavailableError) as info:
        gap.genomepy_minimap2_index_command(
            ["minimap2", "-d", "genes.mmi", "genes.fa"], working_directory=str(tmp_path)
        )
    assert "minimap2" in info.value.hint


def test_run_index_cmd_checks_the_aligner_name(tmp_path: Path) -> None:
    fake = _fake_program(tmp_path / "bin", "bowtie2-build")
    workdir = tmp_path / "work"
    workdir.mkdir()
    result = gap.genomepy_run_index_cmd(
        "bowtie2",
        ["bowtie2-build", "--threads", "1", "genes.fa", "genes"],
        working_directory=str(workdir),
        executable_path=str(fake),
    )
    assert result["backend"] == "bowtie2"
    with pytest.raises(ConfigurationError, match="name must be one of"):
        gap.genomepy_run_index_cmd("bowtie3", ["bowtie3"], working_directory=str(workdir))


@pytest.mark.skipif(os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live NCBI request")
def test_live_genome_size_of_grch38() -> None:
    size = gap.genomepy_genome_size("GCF_000001405.40")["genome_size"]
    assert size is not None and 3_000_000_000 < size < 3_300_000_000
