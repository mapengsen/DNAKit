"""Explicit parameter extensions for the scientific workflows bundled by DNAKit."""

from __future__ import annotations

import copy
import re
from typing import Any

from dnakit.exceptions import BackendExecutionError, ConfigurationError


def prepare(entry: dict[str, Any], values: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    extension = entry.get("extension", {})
    options = {}
    for alias, canonical in extension.get("aliases", {}).items():
        if alias not in values:
            continue
        value = values.pop(alias)
        if canonical in values and values[canonical] != value:
            raise ConfigurationError(f"Conflicting {alias} and {canonical}.")
        values[canonical] = value
    for name, conversion in extension.get("conversions", {}).items():
        if name not in values:
            continue
        value = values[name]
        if conversion == "join_cr" and isinstance(value, list):
            values[name] = "\r".join(value)
        elif conversion == "as_list" and isinstance(value, str):
            values[name] = [value]
        elif conversion == "ensembl_gene" and not value.startswith("ENSG"):
            from ._runtime import call

            result = call("ensembl_lookup_gene", {"gene_id": value, "species": "homo_sapiens"})
            if isinstance(result, dict) and isinstance(result.get("data"), dict):
                result = result["data"]
            if not isinstance(result, dict) or not re.fullmatch(
                r"ENSG[0-9]+(?:\.[0-9]+)?", str(result.get("id", ""))
            ):
                raise BackendExecutionError(
                    f"Gene symbol could not be resolved unambiguously: {value}"
                )
            values[name] = result["id"]
    handler = extension.get("handler")
    if handler == "page_size":
        if "size" in values:
            options["size"] = values.pop("size")
    elif handler == "p_value_filter":
        if "p_value_threshold" in values:
            options["p_value_threshold"] = values.pop("p_value_threshold")
    elif handler == "species_filter" and "species" in values:
        species = values.pop("species")
        taxon = {"human": "9606", "Homo sapiens": "9606", "mouse": "10090"}.get(species, species)
        if not str(taxon).isdigit():
            raise ConfigurationError("species must be a supported name or NCBI taxonomy ID.")
        values["query"] = f"({values['query']}) AND taxid:{taxon}"
    elif handler == "arxiv_category" and "category" in values:
        category = values.pop("category")
        if not re.fullmatch(r"[A-Za-z0-9.-]+", category):
            raise ConfigurationError("Invalid arXiv category.")
        values["query"] = f"({values['query']}) AND cat:{category}*"
    elif handler == "literature_source" and "source" in values:
        source = values.pop("source")
        if not re.fullmatch(r"[A-Za-z0-9]+", source):
            raise ConfigurationError("Invalid Europe PMC source.")
        values["query"] = f"({values['query']}) AND SRC:{source}"
    elif handler == "vep_region":
        region = values.get("region", "")
        if "/" in region and "allele" not in values:
            region, allele = region.rsplit("/", 1)
            parts = region.split(":")
            if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
                region = parts[0] + ":" + parts[1] + "-" + parts[2]
            values["region"], values["allele"] = region, allele
    elif handler == "region_coordinates":
        if any(key in values for key in ("chrom", "start", "end")):
            if not all(key in values for key in ("chrom", "start", "end")):
                raise ConfigurationError("chrom, start and end must be supplied together.")
            chrom, start, end = (values.pop(key) for key in ("chrom", "start", "end"))
            if end < start:
                raise ConfigurationError("end must not precede start.")
            region = f"{chrom}:{start}-{end}"
            if "region" in values and values["region"] != region:
                raise ConfigurationError("Conflicting region and coordinates.")
            values["region"] = region
        if values.get("species") == "human":
            values["species"] = "homo_sapiens"
    elif handler == "census_filters":
        import json

        gene, tissue = values.pop("gene", None), values.pop("tissue", None)
        for key, condition in (
            ("var_value_filter", f"feature_name == {json.dumps(gene)}" if gene else None),
            (
                "obs_value_filter",
                f"tissue == {json.dumps(tissue)}" if tissue and tissue != "all" else None,
            ),
        ):
            if condition:
                values[key] = f"({values[key]}) and ({condition})" if values.get(key) else condition
    if entry["name"] == "ch_ipatlas_enrichment_analysis" and values.get("cell_type_class") == "all":
        values.pop("cell_type_class")
    return values, options


def finish(result: object, options: dict[str, Any]) -> object:
    if not options:
        return result
    result = copy.deepcopy(result)
    if "size" in options:
        size = options["size"]
        found = False

        def visit(value: object) -> None:
            nonlocal found
            if isinstance(value, dict):
                if isinstance(value.get("rows"), list):
                    found = True
                    rows = value["rows"]
                    value["rows"] = rows[:size]
                    value["returned"] = len(value["rows"])
                    value["truncated"] = len(value["rows"]) < value.get("count", len(rows))
                    value["dnakit_limit_scope"] = (
                        "limits the provider page; additional pages are not fetched"
                    )
                else:
                    for item in value.values():
                        visit(item)

        visit(result)
        if not found:
            raise BackendExecutionError("Provider page has no rows collection to limit.")
    if "p_value_threshold" in options:
        if not isinstance(result, dict):
            raise BackendExecutionError("Provider result has no association records to filter.")
        rows = result.get("data") if isinstance(result, dict) else None
        if not isinstance(rows, list):
            raise BackendExecutionError("Provider result has no association records to filter.")
        threshold = options["p_value_threshold"]
        selected = []
        for row in rows:
            value = row.get("p_value") if isinstance(row, dict) else None
            try:
                if value is not None and 0 <= float(value) <= threshold:
                    selected.append(row)
            except (TypeError, ValueError):
                continue
        result["data"] = selected
        result["dnakit_filter"] = {
            "p_value_threshold": threshold,
            "input_page_count": len(rows),
            "returned": len(selected),
            "scope": "current provider page",
        }
    return result
