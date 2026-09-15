"""Additional scientific query modes used by the bundled workflows."""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from typing import Any

from dnakit.exceptions import BackendExecutionError, ConfigurationError


def route(name: str, values: dict[str, Any]) -> Callable[[], object] | None:
    selectors: dict[str, tuple[str, Callable[[str, dict[str, Any]], object]]] = {
        "execute_tool": ("tool_name", _execute),
        "pdbe_get_uniprot_mappings": ("uniprot_id", _structures),
        "bio_rxiv_list_recent_preprints": ("query", _preprints),
        "med_rxiv_get_preprint": ("query", _preprints),
        "encode_search_experiments": ("biosample", _encode),
        "cellx_gene_get_cell_metadata": ("gene", _cell_metadata),
        "ch_ipatlas_get_peak_data": ("gene", _peaks),
    }
    selection = selectors.get(name)
    if selection is None or selection[0] not in values:
        return None
    _, function = selection
    return lambda: function(name, values)


def _execute(name: str, values: dict[str, Any]) -> object:
    from ._runtime import call, describe

    if describe(values["tool_name"])["name"] == "execute_tool":
        raise ConfigurationError("Recursive execute_tool calls are not supported.")
    return call(values["tool_name"], values.get("arguments"))


def _structures(name: str, values: dict[str, Any]) -> object:
    from ._runtime import call

    if values.get("pdb_id"):
        raise ConfigurationError(
            "Use pdb_id for PDB-to-UniProt mapping or uniprot_id for the reverse search."
        )
    response = call(
        "PDBeSIFTS_get_best_structures",
        {"uniprot_accession": values["uniprot_id"], "limit": values.get("limit", 20)},
    )
    data = response.get("data", response) if isinstance(response, dict) else response
    if isinstance(data, dict):
        data = data.get("structures", data.get(values["uniprot_id"], data))
    if not isinstance(data, list):
        raise BackendExecutionError("PDBe did not return the expected structure mapping records.")
    return data


def _preprints(name: str, values: dict[str, Any]) -> object:
    from ._runtime import call

    if values.get("doi"):
        raise ConfigurationError("Supply doi for an exact record or query for a preprint search.")
    # This mode explicitly queries indexed preprints; exact DOI and date-range modes
    # continue to use the source repository's own API.
    server = "medrxiv" if name.startswith("med_") else values.get("server", "biorxiv")
    query = f"({values['query']}) AND SRC:PPR AND PUBLISHER:{server}"
    for key in ("start_date", "end_date"):
        if key in values and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", values[key]):
            raise ConfigurationError(f"{key} must use YYYY-MM-DD.")
    if values.get("start_date") or values.get("end_date"):
        query += (
            f" AND FIRST_PDATE:[{values.get('start_date', '*')} TO {values.get('end_date', '*')}]"
        )
    result = call("europe_pmc_search_articles", {"query": query, "limit": values.get("limit", 10)})
    if isinstance(result, dict):
        result["dnakit_query_mode"] = {
            "provider": "Europe PMC",
            "server": server,
            "scope": "indexed preprints; not exhaustive repository coverage",
        }
    return result


def _encode(name: str, values: dict[str, Any]) -> object:
    from dnakit import encode

    values = dict(values)
    biosample = values.pop("biosample")
    if biosample != "all":
        values["biosample_term_name"] = biosample
    return json.loads(encode.search_experiments(**values))


def _cell_metadata(name: str, values: dict[str, Any]) -> object:
    import cellxgene_census

    organism = values.get("organism", "Homo sapiens")
    limit = values.get("max_cells", 100)
    obs_filter = values.get("obs_value_filter")
    tissue = values.get("tissue")
    if tissue and tissue != "all":
        condition = "tissue == " + json.dumps(tissue)
        obs_filter = f"({obs_filter}) and ({condition})" if obs_filter else condition
    with cellxgene_census.open_soma(
        census_version=values.get("census_version", "stable")
    ) as census:
        experiment = census["census_data"][organism.lower().replace(" ", "_")]
        # Read only the leading matching observation IDs, then fetch expression for
        # that explicitly bounded slice. It is not a representative population sample.
        ids: list[int] = []
        for table in experiment.obs.read(
            value_filter=obs_filter, column_names=["soma_joinid"]
        ).tables():
            ids.extend(table["soma_joinid"].slice(0, limit - len(ids)).to_pylist())
            if len(ids) >= limit:
                break
        if not ids:
            return {"status": "success", "data": [], "num_cells": 0, "gene": values["gene"]}
        adata = cellxgene_census.get_anndata(
            census,
            organism=organism,
            obs_coords=ids,
            var_value_filter="feature_name == " + json.dumps(values["gene"]),
            obs_column_names=values.get("column_names"),
        )
        if adata.n_vars == 0:
            raise ConfigurationError("Gene is absent from this Census release.")
        observations = adata.obs.reset_index().to_dict(orient="records")
        counts = adata.X.sum(axis=1)
        if hasattr(counts, "A1"):
            counts = counts.A1
        for item, count in zip(observations, counts, strict=True):
            item["queried_gene_count"] = float(count)
        # Normalize pandas/numpy scalars, missing values and categories to JSON.
        import pandas as pd

        records = json.loads(pd.DataFrame(observations).to_json(orient="records"))
        return {
            "status": "success",
            "data": records,
            "num_cells": len(records),
            "gene": values["gene"],
            "scope": (
                "leading matching cells in a bounded Census slice; no population-wide inference"
            ),
            "max_cells": limit,
        }


def _peaks(name: str, values: dict[str, Any]) -> object:
    from ._regulatory import gene_peaks

    return gene_peaks(values)


def pdb_get_structure(
    *, pdb_id: str, format: str = "cif", max_bytes: int = 10_000_000, verbose: bool = True
) -> dict[str, object]:
    """Retrieve a bounded RCSB PDB coordinate file as text, without writing it to disk.

    Service: https://www.rcsb.org/docs/programmatic-access/file-download-services
    """
    import hashlib

    import requests

    if not re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8}", pdb_id) or format not in {
        "cif",
        "pdb",
    }:
        raise ConfigurationError("Provide a PDB identifier and cif or pdb format.")
    if (
        isinstance(max_bytes, bool)
        or not isinstance(max_bytes, int)
        or not 1 <= max_bytes <= 100_000_000
    ):
        raise ConfigurationError("max_bytes must be between 1 and 100,000,000.")
    url = f"https://files.rcsb.org/download/{pdb_id}.{format}"
    chunks, size = [], 0
    with requests.get(url, timeout=(10, 60), stream=True) as response:
        response.raise_for_status()
        for chunk in response.iter_content(65536):
            size += len(chunk)
            if size > max_bytes:
                raise BackendExecutionError("PDB coordinate file exceeds max_bytes.")
            chunks.append(chunk)
    content = b"".join(chunks)
    text = content.decode("utf-8")
    if not content or (format == "cif" and not text.lstrip().startswith("data_")):
        raise BackendExecutionError("RCSB returned an invalid coordinate file.")
    return {
        "pdb_id": pdb_id,
        "format": format,
        "data": text,
        "url": url,
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def orphanet_get_disease_genes(*, orpha_code: str, verbose: bool = True) -> object:
    """Retrieve Orphanet's reported gene associations for an ORPHA disorder ID."""
    from ._runtime import call

    return call("Orphanet_get_genes", {"orpha_code": orpha_code})


def cadd_get_scores(*, variant: str, assembly: str = "hg38", verbose: bool = True) -> object:
    """Query indexed CADD annotations through MyVariant; preserve missing scores as missing."""
    from ._runtime import call

    result = call(
        "my_variant_query_variants", {"query": variant, "fields": "cadd", "assembly": assembly}
    )
    if isinstance(result, dict):
        result["dnakit_score_source"] = (
            "CADD annotations indexed by MyVariant; not a new model prediction"
        )
    return result


def string_get_functional_enrichment(
    *, identifiers: str, species: int = 9606, category: str | None = None, verbose: bool = True
) -> object:
    """Run STRING enrichment for carriage-return or newline-separated protein identifiers."""
    from ._runtime import call

    params = {"protein_ids": re.split(r"[\r\n]+", identifiers), "species": species}
    if category is not None:
        params["category"] = category
    return call("STRING_functional_enrichment", params)


EXTRA_FUNCTIONS: dict[str, Callable[..., object]] = {
    fn.__name__: fn
    for fn in (
        pdb_get_structure,
        orphanet_get_disease_genes,
        cadd_get_scores,
        string_get_functional_enrichment,
    )
}
