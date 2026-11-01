"""Tests for allel_verify_functions.allel_vcf_to_zarr.

The VCF is written into tmp_path. The Zarr store written by DNAKit must hold the same arrays that
allel.read_vcf returns for the same file (the layout that allel.vcf_to_zarr intends), and the upstream
vcf_to_zarr call is checked to fail with the zarr 3 API so the reason stays recorded. No network.
"""

from __future__ import annotations

import importlib.util
import socket

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.population import allel_verify_functions as module

requires_allel = pytest.mark.skipif(
    importlib.util.find_spec("allel") is None or importlib.util.find_spec("zarr") is None,
    reason="scikit-allel and zarr are not installed",
)

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t10\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t25\t.\tCT\tC,T\t12.5\t.\t.\tGT\t1|0\t./.\t0/0
chr1\t80\trs3\tG\tA\t99\tPASS\t.\tGT\t1/1\t1/1\t0/1
"""


_ORIGINAL_GETADDRINFO = socket.getaddrinfo
_ORIGINAL_CONNECT = socket.socket.connect


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    # Name resolution and outbound connections are blocked; local socket pairs (used by asyncio in zarr) stay.
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)


@requires_allel
def test_store_matches_read_vcf(tmp_path):
    import allel
    import numpy as np
    import zarr

    vcf = tmp_path / "small.vcf"
    vcf.write_text(VCF)
    store = tmp_path / "small.zarr"
    result = module.allel_vcf_to_zarr(str(vcf), str(store), verbose=False)
    assert result["status"] == "completed"
    assert result["variants"] == 3
    assert result["samples"] == 3
    expected = allel.read_vcf(str(vcf))
    root = zarr.open_group(str(store), mode="r")
    assert set(expected) <= {"/".join(item["name"].split("/")) for item in result["arrays"]}
    for name, values in expected.items():
        stored = np.asarray(root[name][:])
        values = np.asarray(values)
        if values.dtype.kind == "O":
            values = values.astype(str)
        assert stored.shape == values.shape, name
        assert np.array_equal(stored, values), name


@requires_allel
def test_existing_store_is_never_overwritten(tmp_path):
    vcf = tmp_path / "small.vcf"
    vcf.write_text(VCF)
    store = tmp_path / "exists.zarr"
    store.mkdir()
    with pytest.raises(ConfigurationError):
        module.allel_vcf_to_zarr(str(vcf), str(store), verbose=False)


@requires_allel
def test_header_only_vcf_writes_no_store(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    store = tmp_path / "none.zarr"
    result = module.allel_vcf_to_zarr(str(vcf), str(store), verbose=False)
    assert result["variants"] == 0
    assert result["output_path"] is None
    assert not store.exists()


@requires_allel
def test_missing_vcf_is_a_configuration_error(tmp_path):
    with pytest.raises(ConfigurationError):
        module.allel_vcf_to_zarr(str(tmp_path / "missing.vcf"), str(tmp_path / "x.zarr"), verbose=False)


@requires_allel
def test_upstream_vcf_to_zarr_uses_removed_zarr2_api(tmp_path):
    import allel

    vcf = tmp_path / "small.vcf"
    vcf.write_text(VCF)
    with pytest.raises(AttributeError, match="create_dataset"):
        allel.vcf_to_zarr(str(vcf), str(tmp_path / "upstream.zarr"), overwrite=True)
