"""Scientific adapter contracts, aliases, output limits and workflow execution."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import pytest

from dnakit import databases
from dnakit.databases import _runtime
from dnakit.exceptions import BackendExecutionError, ConfigurationError


def provider(monkeypatch: pytest.MonkeyPatch, operation: Any, result: Any) -> Any:
    entry = databases.describe(operation)
    module = importlib.import_module("dnakit.databases._providers." + entry["module"])
    seen = []

    class Fake:
        def __init__(self, config: Any) -> None:
            self.config = config

        def run(self, arguments: Any) -> Any:
            seen.append(arguments)
            return result

    monkeypatch.setattr(module, entry["class"], Fake)
    return seen


def test_gene_alias_preserves_required_validation_and_result_shape(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen = provider(
        monkeypatch, "ensembl_lookup_gene", {"id": "ENSG00000141510", "display_name": "TP53"}
    )
    result = databases.ensembl_lookup_gene(id="TP53")
    assert isinstance(result, dict)
    assert result["display_name"] == "TP53"
    assert seen == [{"gene_id": "TP53", "species": "homo_sapiens"}]
    with pytest.raises(ConfigurationError):
        databases.ensembl_lookup_gene()
    with pytest.raises(ConfigurationError, match="Conflicting"):
        databases.ensembl_lookup_gene(id="BRCA1", gene_id="TP53")


def test_sequence_and_coordinate_aliases_are_translated(monkeypatch: pytest.MonkeyPatch) -> None:
    seen = provider(monkeypatch, "string_get_network", [{"preferredName_A": "TP53"}])
    assert isinstance(databases.string_get_network(protein_ids=["TP53", "BRCA1"]), list)
    assert seen[0]["identifiers"] == "TP53\rBRCA1"
    seen = provider(monkeypatch, "gnomad_get_sv_by_region", {"data": []})
    databases.gnomad_get_sv_by_region(chrom="17", start=10, end=20)
    assert seen[0]["stop"] == 20 and "end" not in seen[0]


def test_default_operation_is_applied(monkeypatch: pytest.MonkeyPatch) -> None:
    seen = provider(monkeypatch, "gtex_get_single_tissue_eqtls", {"status": "success", "data": []})
    databases.gtex_get_single_tissue_eqtls(gene_symbol="BRCA1")
    assert seen[0]["gencode_id"] == ["BRCA1"]
    assert seen[0]["operation"] == "get_single_tissue_eqtls"


def test_dynamic_dispatch_preserves_the_arguments_parameter(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen = provider(monkeypatch, "ensembl_lookup_gene", {"id": "ENSG1"})
    arguments = {"gene_id": "TP53"}
    result = databases.execute_tool(tool_name="ensembl_lookup_gene", arguments=arguments)
    assert result == {"id": "ENSG1"}
    assert seen == [{"gene_id": "TP53", "species": "homo_sapiens"}]
    assert arguments == {"gene_id": "TP53"}


def test_variant_region_retains_alt_allele(monkeypatch: pytest.MonkeyPatch) -> None:
    seen = provider(monkeypatch, "ensembl_vep_region", {"data": []})
    databases.ensembl_vep_region(region="17:10:10/T", species="homo_sapiens")
    assert seen[0]["region"] == "17:10-10" and seen[0]["allele"] == "T"


def test_page_limit_is_applied_and_truncation_is_visible(monkeypatch: pytest.MonkeyPatch) -> None:
    provider(
        monkeypatch,
        "open_targets_get_associated_drugs_by_target_ensembl_id",
        {"data": {"target": {"drugAndClinicalCandidates": {"rows": [1, 2, 3], "count": 30}}}},
    )
    result = databases.open_targets_get_associated_drugs_by_target_ensembl_id(
        ensemblId="ENSG1", size=2
    )
    assert isinstance(result, dict)
    collection = result["data"]["target"]["drugAndClinicalCandidates"]
    assert collection["rows"] == [1, 2] and collection["truncated"] is True
    assert collection["count"] == 30


@pytest.mark.parametrize(
    "response", [None, {"status": "error", "error": "failed"}, {"value": float("nan")}]
)
def test_invalid_results_are_not_success(monkeypatch: pytest.MonkeyPatch, response: Any) -> None:
    provider(monkeypatch, "ensembl_lookup_gene", response)
    with pytest.raises(BackendExecutionError):
        databases.ensembl_lookup_gene(gene_id="TP53")


def test_result_redaction_and_output_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("TEST_API_KEY", "test-secret-1234")
    provider(monkeypatch, "ensembl_lookup_gene", {"note": "test-secret-1234", "api_key": "hidden"})
    assert databases.ensembl_lookup_gene(gene_id="TP53") == {
        "note": "REDACTED",
        "api_key": "REDACTED",
    }
    with pytest.raises(BackendExecutionError, match="max_result_bytes"):
        databases.call("ensembl_lookup_gene", {"gene_id": "TP53"}, max_result_bytes=5)


def test_explicit_unknown_arguments_are_rejected_before_provider(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen = provider(monkeypatch, "ensembl_lookup_gene", {})
    with pytest.raises(ConfigurationError):
        databases.call("ensembl_lookup_gene", {"gene_id": "TP53", "typo": True})
    assert not seen


def test_scientific_mutations_are_marked_for_mcp() -> None:
    from dnakit.tools.catalog import build_tool_catalog

    specs = {s.name: s for s in build_tool_catalog(verbose=False)}
    assert specs["dnakit.workflows.run"].effect == "external_process"
    assert specs["dnakit.databases.call"].requires_confirmation
    assert specs["dnakit.encode.download_files"].requires_confirmation
    assert specs["dnakit.databases.cellx_gene_download_h5ad"].requires_confirmation
    assert not specs["dnakit.databases.search"].requires_confirmation


def test_enrichment_and_coordinate_names_exist_without_upstream_collection() -> None:
    assert callable(databases.string_get_functional_enrichment)
    assert callable(databases.pdb_get_structure)
    assert len(_runtime.definitions()) >= 445


def test_gene_peak_query_intersects_locus_in_bed_coordinates(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    provider(
        monkeypatch,
        "ensembl_lookup_gene",
        {
            "id": "ENSG1",
            "assembly_name": "GRCh38",
            "seq_region_name": "17",
            "start": 101,
            "end": 200,
        },
    )
    peak = tmp_path / "peaks.bed"
    peak.write_text(
        "chr17\t90\t100\texcluded-boundary\nchr17\t100\t101\tfirst-base\n"
        "chr17\t199\t200\tlast-base\nchr17\t200\t201\texcluded-right\n"
        "chr1\t100\t150\twrong-chrom\n"
    )
    result = databases.ch_ipatlas_get_peak_data(gene="BRCA1", peak_file=str(peak), flank_bp=0)
    assert isinstance(result, dict)
    assert result["matched"] == 2 and result["status"] == "success"
    assert [row["fields"][0] for row in result["data"]] == ["first-base", "last-base"]


def test_regulatory_instructions_are_not_reported_as_enrichment_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provider(
        monkeypatch,
        "ch_ipatlas_enrichment_analysis",
        {
            "status": "success",
            "data": {
                "message": "Submit gene list for enrichment analysis",
            },
        },
    )
    result = databases.ch_ipatlas_enrichment_analysis(gene="BRCA1")
    assert isinstance(result, dict)
    assert result["status"] == "requires_submission"
    assert result["analysis_completed"] is False
