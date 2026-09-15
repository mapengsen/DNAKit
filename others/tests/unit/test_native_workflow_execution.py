"""Installed-style workflows, scientific invariants and ENCODE resource lifecycles."""

from __future__ import annotations

import asyncio
import hashlib
import json
import random
from pathlib import Path
from typing import Any

import httpx
import pytest

from dnakit import encode, workflows
from dnakit.exceptions import ConfigurationError


@pytest.mark.parametrize(
    "directory",
    [
        "databases/_providers",
        "encode/_backend",
        "workflows/common",
        "workflows/prs",
        "workflows/support",
    ],
)
def test_installed_backend_resources_match_recorded_checksums(directory: str) -> None:
    root = Path(workflows.__file__).parents[1] / directory
    manifest = json.loads((root / "UPSTREAM.json").read_text())
    for entry in manifest["files"]:
        resource = root / entry["path"]
        expected = entry.get("adapted_sha256", entry.get("sha256"))
        assert hashlib.sha256(resource.read_bytes()).hexdigest() == expected, (
            f"Bundled resource bytes changed: {directory}/{entry['path']}"
        )


def test_simulated_children_inherit_one_allele_from_each_parent_without_global_rng_change() -> None:
    locus = {
        "chromosome": "1",
        "position": 1,
        "ref": "A",
        "alt": "C",
        "dominance": "additive",
        "effect_size": 1,
    }
    father = {"id": "father", "sex": "Male", "loci": {"L1": {**locus, "alleles": ["A", "A"]}}}
    mother = {"id": "mother", "sex": "Female", "loci": {"L1": {**locus, "alleles": ["C", "C"]}}}
    traits = {"traits": {"dosage": {"loci": [dict(locus, id="L1", effect=1)]}}}
    state = random.getstate()
    original = json.dumps([father, mother, traits])
    result = workflows.simulate_inheritance(
        father,
        mother,
        traits,
        {"diseases": {}},
        mutation_rate=0,
        seed=13,
        offspring_count=8,
    )
    assert len(result["offspring"]) == 8
    for child in result["offspring"]:
        assert child["loci"]["L1"]["alleles"] == ["A", "C"]
        assert not child["mutations"]
    assert result == workflows.simulate_inheritance(
        father,
        mother,
        traits,
        {"diseases": {}},
        mutation_rate=0,
        seed=13,
        offspring_count=8,
    )
    assert random.getstate() == state
    assert json.dumps([father, mother, traits]) == original


def test_report_keeps_supplied_spanish_table_content_and_records_real_hashes(
    tmp_path: Path,
) -> None:
    pytest.importorskip("reportlab")
    pypdf = pytest.importorskip("pypdf")
    source = tmp_path / "report.md"
    source.write_text(
        "# Resultados\nTexto proporcionado.\n\n| Gen | Valor |\n|---|---|\n| BRCA1 | 12 |\n"
    )
    result = workflows.write_wes_report(source, tmp_path / "pdf", language="es")
    pdf = Path(result["output"])
    text = "\n".join(page.extract_text() for page in pypdf.PdfReader(pdf).pages)
    assert "Informe WES" in text and "BRCA1" in text and "12" in text
    assert result["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert result["output_sha256"] == hashlib.sha256(pdf.read_bytes()).hexdigest()
    with pytest.raises(ConfigurationError, match="already exists"):
        workflows.write_wes_report(source, tmp_path / "pdf")
    assert result["output_sha256"] == hashlib.sha256(pdf.read_bytes()).hexdigest()


def test_packaged_coordinate_script_uses_current_environment_and_leaves_bundle_unchanged(
    tmp_path: Path,
) -> None:
    from dnakit.skills import get_skill

    resource = workflows.resource_path("genomic-coordinates", "scripts/convert_coords.py")
    root = resource.parent.parent
    paths_before = {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}
    result = workflows.run(
        "genomic-coordinates",
        "scripts/convert_coords.py",
        arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
        workdir=tmp_path,
        verbose=False,
    )
    assert result["returncode"] == 0, result["output"]
    parsed = json.loads(result["output"])
    assert "1000" in json.dumps(parsed)
    assert result["upstream_revision"] == get_skill("genomic-coordinates")["upstream_revision"]
    assert paths_before == {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}


def test_support_and_skill_demo_resolve_shared_code_outside_source_directory(
    tmp_path: Path,
) -> None:
    result = workflows.run(
        "recombinator", "recombinator.py", arguments=["--demo"], workdir=tmp_path, verbose=False
    )
    assert result["returncode"] == 0, result["output"]
    assert "synthetic" in result["output"].lower()
    for name in (
        "clinical-variant-reporter",
        "clinical-trial-finder",
        "omics-target-evidence-mapper",
    ):
        result = workflows.run_support(name, arguments=["--help"], workdir=tmp_path, verbose=False)
        assert result["returncode"] == 0, result["output"]


@pytest.mark.parametrize("inside_event_loop", [False, True])
def test_encode_query_preserves_filters_and_closes_isolated_http_clients(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    inside_event_loop: bool,
) -> None:
    clients, requests = [], []
    original = httpx.AsyncClient

    def handle(request: Any) -> Any:
        requests.append(request)
        return httpx.Response(
            200,
            json={
                "total": 12,
                "@graph": [
                    {"accession": "ENCSR000AAA", "assay_title": "ATAC-seq", "status": "released"},
                ],
            },
        )

    def factory(**kwargs: Any) -> Any:
        client = original(**kwargs, transport=httpx.MockTransport(handle))
        clients.append(client)
        return client

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(httpx, "AsyncClient", factory)

    def query() -> Any:
        return json.loads(
            encode.search_experiments(
                assay_title="ATAC-seq",
                biosample_term_name="K562",
                limit=2,
                offset=4,
            )
        )

    if inside_event_loop:

        async def execute() -> Any:
            return query()

        result = asyncio.run(execute())
    else:
        result = query()
    assert result["results"][0]["accession"] == "ENCSR000AAA"
    assert result["has_more"] is True and result["next_offset"] == 6
    assert requests[0].url.params["biosample_ontology.term_name"] == "K562"
    assert requests[0].url.params["limit"] == "2"
    assert clients and all(client.is_closed for client in clients)
    assert not list(tmp_path.rglob("*.db")), "A public query should not create a tracker database"
