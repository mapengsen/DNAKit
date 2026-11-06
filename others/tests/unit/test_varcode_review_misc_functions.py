"""Tests for dnakit.molecular_lab.varcode_review_misc_functions (Fig. 2 coverage gap, batch 0 review).

The reference-data download is tested with a fake pyensembl module that records calls; nothing is
downloaded and the network is never used. The gene-overlap wrapper runs on the installed varcode and
pyensembl packages with a small GTF written into tmp_path and the shared consistency genome. The tests
import the module directly, so they run before the area package is wired.
"""

from __future__ import annotations

import importlib
import inspect
import sys
import types
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    DownloadError,
)

misc = importlib.import_module("dnakit.molecular_lab.varcode_review_misc_functions")

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENOME = DATA / "genome.fa"


class _Cache:
    cache_directory_path = "/fake/pyensembl-cache"


class _Species:
    latin_name = "Homo sapiens"


class _FakeGenome:
    def __init__(
        self,
        name: str,
        *,
        present: bool,
        log: list[tuple[str, str]],
        download_error: Exception | None = None,
    ) -> None:
        self.name = name
        self.release = 110
        self.species = _Species()
        self.download_cache = _Cache()
        self._present = present
        self._log = log
        self._download_error = download_error

    def required_local_files_exist(self, empty_files_ok: bool = False) -> bool:
        return self._present

    def download(self, overwrite: bool = False, show_progress: bool = False) -> None:
        self._log.append(("download", self.name))
        if self._download_error is not None:
            raise self._download_error
        self._present = True

    def index(self, overwrite: bool = False, show_progress: bool = False) -> None:
        self._log.append(("index", self.name))


def _fake_pyensembl(monkeypatch: pytest.MonkeyPatch, genomes: dict[str, _FakeGenome]) -> None:
    fake = types.ModuleType("pyensembl")

    def genome_for_reference_name(reference_name: str, allow_older_downloaded_release: bool = True):
        if reference_name not in genomes:
            raise ValueError(f"no pyensembl release for {reference_name}")
        return genomes[reference_name]

    fake.genome_for_reference_name = genome_for_reference_name  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "pyensembl", fake)


def _gtf_row(contig: str, feature: str, start: int, end: int, attributes: str) -> str:
    return "\t".join([contig, "test", feature, str(start), str(end), ".", "+", ".", attributes])


def _write_gtf(path: Path) -> Path:
    gene = 'gene_id "G1"; gene_name "G1"; gene_biotype "protein_coding";'
    transcript = 'gene_id "G1"; transcript_id "T1"; transcript_biotype "protein_coding";'
    lines = [
        _gtf_row("chr1", "gene", 101, 400, gene),
        _gtf_row("chr1", "transcript", 101, 400, transcript),
        _gtf_row("chr1", "exon", 101, 400, f'{transcript} exon_number "1";'),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_public_functions_take_keyword_only_verbose() -> None:
    for function in (misc.varcode_install_reference_data, misc.varcode_variant_gene_overlaps):
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_missing_pyensembl_reports_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pyensembl", None)
    with pytest.raises(BackendUnavailableError) as caught:
        misc.varcode_install_reference_data(["GRCh38"])
    assert "bio-pyensembl" in str(caught.value.hint)


def test_install_downloads_and_indexes_missing_files(monkeypatch: pytest.MonkeyPatch) -> None:
    log: list[tuple[str, str]] = []
    _fake_pyensembl(monkeypatch, {"GRCh38": _FakeGenome("GRCh38", present=False, log=log)})
    result = misc.varcode_install_reference_data(["GRCh38"])
    assert log == [("download", "GRCh38"), ("index", "GRCh38")]
    row = result["references"][0]
    assert row["files_present_before"] is False
    assert row["downloaded"] is True and row["indexed"] is True
    assert row["release"] == 110 and row["species"] == "Homo sapiens"
    assert row["cache_directory"] == "/fake/pyensembl-cache"
    assert result["status"] == "completed" and result["backend"] == "varcode"


def test_install_reuses_present_files_without_network(monkeypatch: pytest.MonkeyPatch) -> None:
    log: list[tuple[str, str]] = []
    _fake_pyensembl(monkeypatch, {"GRCh38": _FakeGenome("GRCh38", present=True, log=log)})
    result = misc.varcode_install_reference_data(["GRCh38", "GRCh38"])
    assert log == []
    assert len(result["references"]) == 1
    row = result["references"][0]
    assert row["files_present_before"] is True
    assert row["downloaded"] is False and row["indexed"] is False


def test_install_download_failure_is_download_error(monkeypatch: pytest.MonkeyPatch) -> None:
    log: list[tuple[str, str]] = []
    genome = _FakeGenome(
        "GRCh38", present=False, log=log, download_error=OSError("connection refused")
    )
    _fake_pyensembl(monkeypatch, {"GRCh38": genome})
    with pytest.raises(DownloadError):
        misc.varcode_install_reference_data(["GRCh38"])
    assert ("index", "GRCh38") not in log


def test_install_unknown_reference_is_configuration_error(monkeypatch: pytest.MonkeyPatch) -> None:
    _fake_pyensembl(monkeypatch, {})
    with pytest.raises(ConfigurationError):
        misc.varcode_install_reference_data(["NoSuchAssembly"])


@pytest.mark.parametrize("names", [[], "GRCh38", ["bad name"], [""]])
def test_install_rejects_invalid_names(names: Any) -> None:
    with pytest.raises(ConfigurationError):
        misc.varcode_install_reference_data(names)


def test_gene_overlaps_validates_variants_before_building_reference(tmp_path: Path) -> None:
    gtf = _write_gtf(tmp_path / "annotation.gtf")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    bad_inputs: list[Any] = [
        [],
        "chr1:1:A>G",
        [{"contig": "chr1", "position": 0, "ref": "A", "alt": "G"}],
        [{"contig": "chr1", "position": 5, "ref": "X", "alt": "G"}],
        [{"contig": "chr1", "position": 5, "ref": "A"}],
        [{"contig": "chr1", "position": True, "ref": "A", "alt": "G"}],
    ]
    for variants in bad_inputs:
        with pytest.raises(ConfigurationError):
            misc.varcode_variant_gene_overlaps(str(gtf), str(fasta), variants)


def test_gene_overlaps_reports_genes_and_annotation_errors(tmp_path: Path) -> None:
    pytest.importorskip("varcode")
    pytest.importorskip("pyensembl")
    pysam = pytest.importorskip("pysam")
    if not GENOME.exists():
        pytest.skip("shared consistency genome.fa not present")
    gtf = _write_gtf(tmp_path / "annotation.gtf")
    with pysam.FastaFile(str(GENOME)) as fasta:
        base = fasta.fetch("chr1", 149, 150).upper()
    if base not in "ACGT":
        base = "A"
    variants = [
        {"contig": "chr1", "position": 150, "ref": base, "alt": "C" if base != "C" else "G"},
        {"contig": "chr1", "position": 1000, "ref": "A", "alt": "G"},
        {"contig": "chrZ", "position": 150, "ref": "A", "alt": "G"},
    ]
    result = misc.varcode_variant_gene_overlaps(
        str(gtf), str(GENOME), variants, cache_dir=str(tmp_path / "cache")
    )
    assert result["variant_count"] == 3
    inside, outside, unknown = result["variants"]
    assert [gene["gene_id"] for gene in inside["genes"]] == ["G1"]
    assert inside["genes"][0]["start"] == 101 and inside["genes"][0]["end"] == 400
    assert inside["genes"][0]["strand"] == "+"
    assert outside["genes"] == [] and outside["error"] is None
    assert unknown["error"] is not None and unknown["genes"] == []
    assert result["annotation_failures"] == 1
    assert result["reference"]["reused"] in (True, False)


def test_gene_overlaps_raise_on_error_raises(tmp_path: Path) -> None:
    pytest.importorskip("varcode")
    pytest.importorskip("pyensembl")
    pytest.importorskip("pysam")
    if not GENOME.exists():
        pytest.skip("shared consistency genome.fa not present")
    gtf = _write_gtf(tmp_path / "annotation.gtf")
    with pytest.raises(BackendExecutionError):
        misc.varcode_variant_gene_overlaps(
            str(gtf),
            str(GENOME),
            [{"contig": "chrZ", "position": 150, "ref": "A", "alt": "G"}],
            cache_dir=str(tmp_path / "cache"),
            raise_on_error=True,
        )
