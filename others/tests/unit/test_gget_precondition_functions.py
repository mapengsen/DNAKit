"""Tests for dnakit.databases.gget_precondition_functions (CELLxGENE queries and COSMIC reference downloads).

No test uses the network. The CELLxGENE tests install a fake ``cellxgene_census`` module in ``sys.modules`` and
check the arguments gget passes to it. The COSMIC tests replace ``urllib.request.urlopen`` with a fake that serves
an in-memory archive, and check which requests carry the account, that the archive is extracted safely, and that
the password never reaches a result.
"""

from __future__ import annotations

import contextlib
import gzip
import inspect
import io
import json
import sys
import tarfile
import types
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("gget")
pd = pytest.importorskip("pandas")

from dnakit.databases import gget_precondition_functions as gp  # noqa: E402
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError  # noqa: E402


class _FakeAdata:
    """Stands in for the AnnData that cellxgene_census returns: two genes and three cells."""

    n_obs = 3
    n_vars = 2

    def __init__(self) -> None:
        self.obs = pd.DataFrame(
            {"tissue": ["lung", "lung", "blood"], "assay": ["10x", "10x", "smart-seq2"]}
        )
        self.var_names = pd.Index(["ACE2", "SLC5A1"])

    def write(self, path: str) -> None:
        Path(path).write_bytes(b"fake h5ad")


class _FakeReader:
    def concat(self) -> _FakeReader:
        return self

    def to_pandas(self) -> pd.DataFrame:
        return pd.DataFrame({"tissue": ["lung", "lung", "blood"], "cell_type": ["a", "b", "c"]})


def _fake_census(calls: list[tuple[Any, ...]]) -> types.ModuleType:
    module = types.ModuleType("cellxgene_census")

    @contextlib.contextmanager
    def open_soma(census_version: str = "stable"):  # type: ignore[no-untyped-def]
        calls.append(("open_soma", census_version))
        yield {
            "census_data": {
                "homo_sapiens": types.SimpleNamespace(
                    obs=types.SimpleNamespace(
                        read=lambda value_filter=None, column_names=None: _FakeReader()
                    )
                )
            }
        }

    def get_anndata(
        census, organism, var_value_filter=None, obs_value_filter=None, column_names=None
    ):  # type: ignore[no-untyped-def]
        calls.append(("get_anndata", organism, var_value_filter, obs_value_filter, column_names))
        return _FakeAdata()

    module.open_soma = open_soma  # type: ignore[attr-defined]
    module.get_anndata = get_anndata  # type: ignore[attr-defined]
    return module


def _tar_bytes(members: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return buffer.getvalue()


def _install_urlopen(
    monkeypatch: pytest.MonkeyPatch, archive: bytes
) -> list[tuple[str, dict[str, str]]]:
    calls: list[tuple[str, dict[str, str]]] = []

    def fake_urlopen(request: Any, timeout: float | None = None) -> io.BytesIO:
        calls.append((request.full_url, dict(request.header_items())))
        if request.full_url.startswith("https://cancer.sanger.ac.uk/"):
            return io.BytesIO(
                json.dumps({"url": "https://downloads.example.org/archive.tar"}).encode()
            )
        return io.BytesIO(archive)

    monkeypatch.setattr(gp.urllib.request, "urlopen", fake_urlopen)
    return calls


def test_public_functions_are_keyword_verbose() -> None:
    for name in (
        "gget_cellxgene_query",
        "gget_cosmic_download_reference",
        "gget_cosmic_select_reference",
    ):
        parameter = inspect.signature(getattr(gp, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


# --- CELLxGENE --------------------------------------------------------------------------------


def test_cellxgene_missing_package_raises_backend_unavailable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    with pytest.raises(BackendUnavailableError, match="cellxgene_census"):
        gp.gget_cellxgene_query(str(tmp_path / "out.h5ad"), genes=["ACE2"], verbose=False)


def test_cellxgene_rejects_bad_inputs_before_any_import(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"\.h5ad"):
        gp.gget_cellxgene_query(str(tmp_path / "out.csv"), genes=["ACE2"], verbose=False)
    with pytest.raises(ConfigurationError, match="species"):
        gp.gget_cellxgene_query(str(tmp_path / "a.h5ad"), species="gallus_gallus", verbose=False)
    with pytest.raises(ConfigurationError, match="genes"):
        gp.gget_cellxgene_query(str(tmp_path / "b.h5ad"), genes=["ACE2; DROP"], verbose=False)
    with pytest.raises(ConfigurationError, match="tissue"):
        gp.gget_cellxgene_query(str(tmp_path / "c.h5ad"), tissue=["lung'"], verbose=False)


def test_cellxgene_query_passes_gene_and_tissue_filters(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[Any, ...]] = []
    monkeypatch.setitem(sys.modules, "cellxgene_census", _fake_census(calls))
    out = tmp_path / "expression.h5ad"
    result = gp.gget_cellxgene_query(
        str(out),
        genes=["ACE2", "SLC5A1"],
        tissue=["lung"],
        census_version="2025-11-08",
        verbose=False,
    )
    assert out.read_bytes() == b"fake h5ad"
    assert ("open_soma", "2025-11-08") in calls
    get_call = next(call for call in calls if call[0] == "get_anndata")
    assert get_call[1] == "homo_sapiens"
    assert get_call[2] == "feature_name in ['ACE2', 'SLC5A1']"
    assert "tissue" in get_call[3]
    assert result["counts"] == {"cells": 3, "genes": 2}
    assert result["var_names"] == ["ACE2", "SLC5A1"]
    assert result["obs_columns"] == ["tissue", "assay"]
    assert result["filters"] == {"tissue": ["lung"]}


def test_cellxgene_metadata_only_writes_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[Any, ...]] = []
    monkeypatch.setitem(sys.modules, "cellxgene_census", _fake_census(calls))
    out = tmp_path / "metadata.csv"
    result = gp.gget_cellxgene_query(str(out), meta_only=True, tissue=["lung"], verbose=False)
    assert out.is_file()
    assert "tissue" in out.read_text(encoding="utf-8")
    assert result["counts"] == {"rows": 3}
    assert result["meta_only"] is True


# --- COSMIC -----------------------------------------------------------------------------------


def test_cosmic_select_downloads_extracts_and_unzips(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"Gene\tCount\nTP53\t3\n"
    archive = _tar_bytes({"Cosmic_MutantCensus_v99_GRCh37.tsv.gz": gzip.compress(payload)})
    calls = _install_urlopen(monkeypatch, archive)
    result = gp.gget_cosmic_select_reference(
        str(tmp_path / "ref"),
        "census",
        37,
        99,
        email="user@example.org",
        password="s3cret-pass",
        verbose=False,
    )
    path = Path(result["file_path"])
    assert path.read_bytes() == payload
    assert path.parent.name == "Cosmic_MutantCensus_Tsv_v99_GRCh37"
    assert result["downloaded"] is True
    assert "grch37/cosmic/v99/Cosmic_MutantCensus_Tsv_v99_GRCh37.tar" in calls[0][0]
    assert "Authorization" in calls[0][1]
    assert calls[1][0] == "https://downloads.example.org/archive.tar"
    assert "Authorization" not in calls[1][1]
    assert "s3cret-pass" not in json.dumps(result)


def test_cosmic_select_reuses_existing_file_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _install_urlopen(monkeypatch, b"")
    folder = tmp_path / "ref" / "Cosmic_MutantCensus_Tsv_v99_GRCh37"
    folder.mkdir(parents=True)
    existing = folder / "Cosmic_MutantCensus_v99_GRCh37.tsv"
    existing.write_text("kept\n", encoding="utf-8")
    result = gp.gget_cosmic_select_reference(
        str(tmp_path / "ref"),
        "census",
        37,
        99,
        email="user@example.org",
        password="pw",
        verbose=False,
    )
    assert result["downloaded"] is False
    assert result["file_path"] == str(existing)
    assert existing.read_text(encoding="utf-8") == "kept\n"
    assert calls == []


def test_cosmic_select_keeps_gz_when_not_unzipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = _tar_bytes({"Cosmic_GenomeScreensMutant_v98_GRCh38.tsv.gz": gzip.compress(b"x\n")})
    _install_urlopen(monkeypatch, archive)
    result = gp.gget_cosmic_select_reference(
        str(tmp_path / "ref"),
        "genome_screen",
        38,
        98,
        email="user@example.org",
        password="pw",
        unzip=False,
        verbose=False,
    )
    assert result["file_path"].endswith("Cosmic_GenomeScreensMutant_v98_GRCh38.tsv.gz")
    assert Path(result["file_path"]).is_file()


def test_cosmic_cancer_requires_grch37(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="GRCh37"):
        gp.gget_cosmic_select_reference(
            str(tmp_path), "cancer", 38, 99, email="user@example.org", password="pw", verbose=False
        )


def test_cosmic_rejects_unknown_project_and_bad_email(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="cosmic_project"):
        gp.gget_cosmic_select_reference(
            str(tmp_path), "transcriptome", 37, 99, email="u@e.org", password="pw", verbose=False
        )
    with pytest.raises(ConfigurationError, match="email"):
        gp.gget_cosmic_select_reference(
            str(tmp_path), "census", 37, 99, email="not-an-address", password="pw", verbose=False
        )
    with pytest.raises(ConfigurationError, match="password"):
        gp.gget_cosmic_select_reference(
            str(tmp_path), "census", 37, 99, email="u@e.org", password="", verbose=False
        )


def test_cosmic_download_refuses_other_hosts_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _install_urlopen(monkeypatch, b"")
    with pytest.raises(ConfigurationError, match=r"cancer\.sanger\.ac\.uk"):
        gp.gget_cosmic_download_reference(
            "https://example.org/api?path=a.tar",
            str(tmp_path),
            email="u@e.org",
            password="pw",
            verbose=False,
        )
    assert calls == []


def test_cosmic_download_rejects_unsafe_archive_members(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = _tar_bytes({"../escape.tsv": b"bad\n"})
    _install_urlopen(monkeypatch, archive)
    with pytest.raises(BackendExecutionError, match="Unsafe member"):
        gp.gget_cosmic_download_reference(
            "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch37/cosmic/v99/Bad_v99.tar&bucket=downloads",
            str(tmp_path / "out"),
            email="user@example.org",
            password="pw",
            verbose=False,
        )
    assert not (tmp_path / "out" / "escape.tsv").exists()
    assert not (tmp_path / "escape.tsv").exists()


def test_cosmic_download_reports_missing_download_address(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def failing(request: Any, timeout: float | None = None) -> io.BytesIO:
        return io.BytesIO(b"<html>login failed</html>")

    monkeypatch.setattr(gp.urllib.request, "urlopen", failing)
    with pytest.raises(BackendExecutionError, match="download address"):
        gp.gget_cosmic_download_reference(
            "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch37/cosmic/v99/X_v99.tar&bucket=downloads",
            str(tmp_path / "out"),
            email="user@example.org",
            password="pw",
            verbose=False,
        )
