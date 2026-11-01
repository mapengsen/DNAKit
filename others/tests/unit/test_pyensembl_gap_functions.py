"""Tests for dnakit.annotation.pyensembl_gap_functions (pyensembl listings and installed-cache lookups).

Listing functions replace one pyensembl upstream call with a recording fake, so no request is sent. The cache
functions run for real against a tiny Ensembl-style release written under tmp_path into a temporary
PYENSEMBL_CACHE_DIR; the release is indexed offline with pyensembl's own index step. Nothing is downloaded, and
the socket guard turns any network attempt into a test failure. Tests skip cleanly when pyensembl is missing.
"""

from __future__ import annotations

import datetime as dt
import gzip
import inspect
import socket
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from dnakit.annotation import pyensembl_gap_functions as pg
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

ACCESSION = "GCA_000001405.29"
GTF_TEXT = (
    '1\ttest\tgene\t100\t300\t.\t+\t.\tgene_id "ENSG00000000001"; gene_version "1"; '
    'gene_name "TESTG1"; gene_biotype "protein_coding";\n'
    '1\ttest\ttranscript\t100\t300\t.\t+\t.\tgene_id "ENSG00000000001"; gene_version "1"; '
    'transcript_id "ENST00000000001"; transcript_version "1"; transcript_name "TESTG1-201"; '
    'transcript_biotype "protein_coding"; gene_name "TESTG1"; gene_biotype "protein_coding";\n'
    '1\ttest\texon\t100\t300\t.\t+\t.\tgene_id "ENSG00000000001"; gene_version "1"; '
    'transcript_id "ENST00000000001"; transcript_version "1"; exon_number "1"; '
    'exon_id "ENSE00000000001"; exon_version "1"; gene_name "TESTG1"; '
    'gene_biotype "protein_coding"; transcript_name "TESTG1-201"; '
    'transcript_biotype "protein_coding";\n'
)


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


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


@pytest.fixture
def pyensembl() -> Any:
    return pytest.importorskip("pyensembl", reason="pyensembl is not installed")


@pytest.fixture
def cache_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pyensembl: Any) -> Path:
    """A temporary pyensembl cache; genome objects are cleared so no earlier test's path is reused."""
    root = tmp_path / "pyensembl-cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(root))
    pyensembl.ensembl_release.EnsemblRelease._genome_cache.clear()
    yield root
    pyensembl.ensembl_release.EnsemblRelease._genome_cache.clear()


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = text.encode("utf-8")
    if path.name.endswith(".gz"):
        with gzip.open(path, "wb") as handle:
            handle.write(data)
    else:
        path.write_bytes(data)


def _install_tiny_release(pyensembl: Any, release: int = 112) -> None:
    """Write a one-gene GTF and matching FASTA files into the cache, then index them offline."""
    genome = pyensembl.ensembl_release.EnsemblRelease(release, species="human")
    _write(Path(genome.download_cache.local_path(genome.gtf_url)), GTF_TEXT)
    cdna, ncrna = genome.transcript_fasta_urls
    _write(
        Path(genome.download_cache.local_path(cdna)),
        ">ENST00000000001.1 cdna\nACGTACGTACGTACGTACGT\n",
    )
    _write(
        Path(genome.download_cache.local_path(ncrna)), ">ENST00000000002.1 ncrna\nGGGGCCCCAAAA\n"
    )
    for url in genome.protein_fasta_urls:
        _write(Path(genome.download_cache.local_path(url)), ">ENSP00000000001.1 pep\nMKLLV\n")
    genome.index()
    assert genome.installed()


# --- module contract ---------------------------------------------------------------------------


def test_public_functions_have_keyword_only_verbose_and_docstrings() -> None:
    assert len(pg.__all__) == len(set(pg.__all__))
    for name in pg.__all__:
        function = getattr(pg, name)
        assert callable(function), name
        signature = inspect.signature(function)
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert signature.parameters["verbose"].default is True, name
        assert function.__doc__ and function.__doc__.strip(), name
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, f"{name}.{parameter.name}"


def test_aliases_are_the_same_upstream_objects(pyensembl: Any) -> None:
    """Rows that name a re-export or a shell.* entry point are aliases of one canonical function each."""
    import pyensembl.shell  # not imported by the package __init__

    dated = pyensembl.dated_releases
    assert pyensembl.available_annotation_dates is dated.available_annotation_dates
    assert pyensembl.available_dated_releases is dated.available_dated_releases
    assert pyensembl.shell.fetch_all_dated_releases is dated.fetch_all_dated_releases
    assert pyensembl.cached_release is pyensembl.ensembl_release.cached_release
    assert pyensembl.genome_for_reference_name is pyensembl.reference_name.genome_for_reference_name


def test_missing_pyensembl_reports_the_optional_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    # A None entry makes any import of pyensembl (or of a submodule already loaded) fail like a missing package.
    monkeypatch.setitem(sys.modules, "pyensembl", None)
    monkeypatch.setitem(sys.modules, "pyensembl.dated_releases", None)
    with pytest.raises(BackendUnavailableError) as raised:
        pg.list_ensembl_annotation_dates(ACCESSION)
    assert raised.value.hint is not None and "bio-pyensembl" in raised.value.hint


# --- Ensembl listings (network, faked) ----------------------------------------------------------


def test_annotation_dates_map_arguments_and_return_the_listing(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(["2025_01", "2026_04"])
    monkeypatch.setattr(pyensembl.dated_releases, "available_annotation_dates", fake)

    result = pg.list_ensembl_annotation_dates(ACCESSION, provider="ensembl", refresh=True)

    assert fake.calls == [((ACCESSION, "ensembl"), {"refresh": True})]
    assert result["status"] == "complete"
    assert result["backend"] == "pyensembl"
    assert result["request"] == {
        "assembly_accession": ACCESSION,
        "provider": "ensembl",
        "refresh": True,
    }
    assert result["data"] == {"dates": ["2025_01", "2026_04"], "count": 2}
    stamp = dt.datetime.fromisoformat(result["retrieved_utc"])
    assert stamp.utcoffset() == dt.timedelta(0)


def test_annotation_dates_refuse_bad_accessions_without_calling_upstream(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder([])
    monkeypatch.setattr(pyensembl.dated_releases, "available_annotation_dates", fake)
    for bad in ("GRCh38", "GCA_1.2", "", "   "):
        with pytest.raises(ConfigurationError):
            pg.list_ensembl_annotation_dates(bad)
    with pytest.raises(ConfigurationError):
        pg.list_ensembl_annotation_dates(ACCESSION, provider="Ensembl Genomes")
    with pytest.raises(ConfigurationError):
        pg.list_ensembl_annotation_dates(ACCESSION, refresh="yes")  # type: ignore[arg-type]
    assert fake.calls == []


def test_unreadable_listing_becomes_backend_error(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _unreachable(*args: Any, **kwargs: Any) -> Any:
        raise OSError("No annotation dates found at https://example.invalid/")

    monkeypatch.setattr(pyensembl.dated_releases, "available_annotation_dates", _unreachable)
    with pytest.raises(BackendExecutionError, match="No annotation dates found"):
        pg.list_ensembl_annotation_dates(ACCESSION)


def test_dated_releases_resolve_the_species_and_its_assembly(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(["2026_04"])
    monkeypatch.setattr(pyensembl.dated_releases, "available_dated_releases", fake)

    result = pg.list_ensembl_dated_releases("human", refresh=False)

    assert fake.calls == [(("human",), {"refresh": False})]
    data = result["data"]
    assert data["latin_name"] == "homo_sapiens"
    assert data["assembly_accession"] == ACCESSION
    assert data["provider"] == "ensembl"
    assert data["dates"] == ["2026_04"] and data["count"] == 1


def test_dated_releases_refuse_unknown_and_undated_species_before_listing(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(["2026_04"])
    monkeypatch.setattr(pyensembl.dated_releases, "available_dated_releases", fake)
    with pytest.raises(ConfigurationError, match="Species not found"):
        pg.list_ensembl_dated_releases("not_a_species")
    with pytest.raises(ConfigurationError, match="no dated Ensembl releases"):
        pg.list_ensembl_dated_releases("toxoplasma_gondii")
    assert fake.calls == []


def test_fetch_all_dated_releases_passes_expiry_and_summarises_species(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder({"mus_musculus": ["2025_10"], "homo_sapiens": ["2026_04"], "bos_taurus": None})
    monkeypatch.setattr(pyensembl.dated_releases, "fetch_all_dated_releases", fake)

    result = pg.fetch_all_ensembl_dated_releases(expire_after_days=2.5, timeout_seconds=5)

    (((), kwargs),) = fake.calls
    assert kwargs == {"expire_after": dt.timedelta(days=2.5), "timeout": 5}
    data = result["data"]
    assert data["dates_by_species"] == {
        "bos_taurus": None,
        "homo_sapiens": ["2026_04"],
        "mus_musculus": ["2025_10"],
    }
    assert data["species_count"] == 3
    assert data["species_with_dates"] == 2
    assert data["species_without_dates"] == ["bos_taurus"]


def test_fetch_all_dated_releases_without_any_listing_is_an_error(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(None)
    monkeypatch.setattr(pyensembl.dated_releases, "fetch_all_dated_releases", fake)
    with pytest.raises(BackendExecutionError, match="could not fetch any dated-release listing"):
        pg.fetch_all_ensembl_dated_releases()


def test_fetch_all_dated_releases_validates_numbers_before_listing(
    pyensembl: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder({})
    monkeypatch.setattr(pyensembl.dated_releases, "fetch_all_dated_releases", fake)
    with pytest.raises(ConfigurationError):
        pg.fetch_all_ensembl_dated_releases(expire_after_days=-1)
    with pytest.raises(ConfigurationError):
        pg.fetch_all_ensembl_dated_releases(timeout_seconds=0)
    with pytest.raises(ConfigurationError):
        pg.fetch_all_ensembl_dated_releases(expire_after_days=True)  # type: ignore[arg-type]
    assert fake.calls == []


# --- installed-cache lookups (real pyensembl, tiny fixture, no download) -----------------------


def test_cached_release_describes_an_installed_release_offline(
    pyensembl: Any, cache_root: Path
) -> None:
    _install_tiny_release(pyensembl, 112)

    result = pg.cached_ensembl_release(112, species="human", verbose=False)

    data = result["data"]
    assert result["status"] == "complete"
    assert result["operation"] == "cached_ensembl_release"
    assert result["request"] == {"release": 112, "species": "human"}
    assert data["species"] == "homo_sapiens"
    assert data["reference_name"] == "GRCh38"
    assert data["release"] == "112"
    assert data["installed"] is True
    assert data["annotation_status"] == "indexed"
    assert data["gene_count"] == 1
    assert data["contig_count"] == 1
    assert data["files"]["gtf"]["status"] == "available"
    assert Path(data["cache_directory"]).is_relative_to(cache_root)


def test_missing_cache_raises_with_install_hint_and_never_downloads(
    pyensembl: Any, cache_root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _forbidden(*args: Any, **kwargs: Any) -> None:
        raise AssertionError("a download was attempted")

    monkeypatch.setattr(pyensembl.genome.Genome, "download", _forbidden)
    monkeypatch.setattr(pyensembl.genome.Genome, "download_genome_fasta", _forbidden)

    with pytest.raises(BackendUnavailableError) as raised:
        pg.cached_ensembl_release(112)

    assert "not installed in the local pyensembl cache" in raised.value.message
    assert raised.value.hint is not None
    assert "pyensembl install --release 112 --species homo_sapiens" in raised.value.hint
    assert not cache_root.exists() or not any(cache_root.rglob("*.gtf.gz"))


def test_cached_release_refuses_invalid_release_values(pyensembl: Any, cache_root: Path) -> None:
    for bad in (0, -3, True, "", None):
        with pytest.raises(ConfigurationError):
            pg.cached_ensembl_release(bad)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        pg.cached_ensembl_release(112, species="not_a_species")


def test_reference_name_chooses_the_installed_release(pyensembl: Any, cache_root: Path) -> None:
    _install_tiny_release(pyensembl, 112)

    result = pg.ensembl_release_for_reference_name("GRCh38", verbose=False)

    assert result["data"]["release"] == "112"
    assert result["data"]["reference_name"] == "GRCh38"
    assert result["data"]["installed"] is True
    assert result["request"] == {
        "reference_name": "GRCh38",
        "allow_older_downloaded_release": True,
    }
    # Upstream returns the newest release (116) when older downloads are not allowed, although its docstring
    # says the newest downloaded one; that release is not installed here, so the lookup is refused.
    with pytest.raises(BackendUnavailableError) as raised:
        pg.ensembl_release_for_reference_name("grch38", allow_older_downloaded_release=False)
    assert raised.value.hint is not None
    assert "--release 116" in raised.value.hint


def test_reference_name_without_installed_release_raises_with_hint(
    pyensembl: Any, cache_root: Path
) -> None:
    with pytest.raises(BackendUnavailableError) as raised:
        pg.ensembl_release_for_reference_name("GRCh38")
    assert raised.value.hint is not None
    # With nothing installed, pyensembl names the newest release of the assembly; that is the one to install.
    assert "pyensembl install --release 116 --species homo_sapiens" in raised.value.hint


def test_unknown_reference_name_is_a_configuration_error(pyensembl: Any, cache_root: Path) -> None:
    with pytest.raises(ConfigurationError, match="not found"):
        pg.ensembl_release_for_reference_name("GRCh99")
    with pytest.raises(ConfigurationError):
        pg.ensembl_release_for_reference_name("")
