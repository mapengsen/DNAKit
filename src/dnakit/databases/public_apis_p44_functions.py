"""Fixed official Ensembl REST and NCBI Datasets v2 biological operations.

Endpoints and allowed GET arguments were checked against the official Ensembl
catalog and NCBI datasets.openapi.yaml on 2026-10-20. Responses retain upstream
page tokens; adapters never report a paginated first page as the complete set.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
from typing import Any
from urllib.parse import quote

from dnakit.databases._network_p44 import _params, _request
from dnakit.download.files import download_file
from dnakit.download.models import DownloadConfig, RemoteFile
from dnakit.exceptions import ConfigurationError
from dnakit.search._http import build_url, require_text


def _segment(value: str, name: str) -> str:
    text = require_text(value, name, max_length=4096)
    if text in {".", ".."}:
        raise ConfigurationError(f"{name} must be an identifier, not a path component.")
    return quote(text, safe="")


def _download(
    url: str, output_path: str, *, overwrite: bool, verbose: bool, timeout: float, max_bytes: int
) -> dict[str, Any]:
    from datetime import datetime, timezone
    from pathlib import Path

    item = download_file(
        RemoteFile(url, filename=Path(output_path).name),
        output_path,
        config=DownloadConfig(
            overwrite=overwrite,
            show_progress=verbose,
            timeout=timeout,
            max_file_bytes=max_bytes,
            max_total_bytes=max_bytes,
        ),
        verbose=verbose,
    )
    return {
        "status": "completed",
        "file": asdict(item),
        "provenance": {
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "service_version": "NCBI Datasets v2",
        },
    }


def ensembl_info_analysis(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List annotation analyses of a species. Official GET /info/analysis/:species; upstream
    defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/analysis/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_biotypes(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List gene and transcript biotypes of a species. Official GET /info/biotypes/:species;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/biotypes/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_biotypes_groups(
    group: str,
    object_type: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List biotype groups. Official GET /info/biotypes/groups/:group/:object_type; upstream
    defaults retained.
    """
    query = _params(parameters, {"object_type", "content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            (
                f"/info/biotypes/groups/{_segment(group, 'group')}"
                f"/{_segment(object_type, 'object_type')}"
            ),
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_biotypes_name(
    name: str,
    object_type: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Biotype properties by name. Official GET /info/biotypes/name/:name/:object_type;
    upstream defaults retained.
    """
    query = _params(parameters, {"object_type", "content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/biotypes/name/{_segment(name, 'name')}/{_segment(object_type, 'object_type')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_compara_methods(
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List comparative genomics methods. Official GET /info/compara/methods; upstream defaults
    retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url("https://rest.ensembl.org", "/info/compara/methods", query.items()),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_compara_species_sets(
    method: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List species sets analysed by a comparative method. Official GET
    /info/compara/species_sets/:method; upstream defaults retained.
    """
    query = _params(parameters, {"content-type", "filter"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/compara/species_sets/{_segment(method, 'method')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_external_dbs(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List external databases linked to a species. Official GET /info/external_dbs/:species;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/external_dbs/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_assembly(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return top-level assembly sequence metadata from GET info/assembly/:species."""
    query = _params(parameters, {"bands", "synonyms"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/assembly/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_genomes_accession(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome record by INSDC accession. Official GET /info/genomes/accession/:accession;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/genomes/accession/{_segment(accession, 'accession')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_genomes_assembly(
    assembly_id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome record by assembly ID. Official GET /info/genomes/assembly/:assembly_id; upstream
    defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/genomes/assembly/{_segment(assembly_id, 'assembly_id')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_genomes_division(
    division_name: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List genomes in a division. Official GET /info/genomes/division/:division_name; upstream
    defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/genomes/division/{_segment(division_name, 'division_name')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_genomes_taxonomy(
    taxon_name: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List genomes under a taxonomy node. Official GET /info/genomes/taxonomy/:taxon_name;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/genomes/taxonomy/{_segment(taxon_name, 'taxon_name')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_variation(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List variation sources of a species. Official GET /info/variation/:species; upstream
    defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/variation/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_variation_consequence_types(
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List variant consequence types. Official GET /info/variation/consequence_types; upstream
    defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url("https://rest.ensembl.org", "/info/variation/consequence_types", query.items()),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_variation_populations(
    species: str,
    population_name: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List individuals of a population. Official GET
    /info/variation/populations/:species/:population_name; upstream defaults retained.
    """
    query = _params(parameters, {"content-type", "filter"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            (
                f"/info/variation/populations/{_segment(species, 'species')}"
                f"/{_segment(population_name, 'population_name')}"
            ),
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_info_variation_populations_list(
    species: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List populations of a species. Official GET /info/variation/populations/:species;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type", "filter"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/info/variation/populations/{_segment(species, 'species')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_map_cds(
    id: str,
    region: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """CDS to genomic coordinate conversion. Official GET /map/cds/:id/:region; upstream
    defaults retained.
    """
    query = _params(parameters, {"coord_system", "content-type", "species"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/map/cds/{_segment(id, 'id')}/{_segment(region, 'region')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_map_translation(
    id: str,
    region: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Protein to genomic coordinate conversion. Official GET /map/translation/:id/:region;
    upstream defaults retained.
    """
    query = _params(parameters, {"coord_system", "content-type", "species"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/map/translation/{_segment(id, 'id')}/{_segment(region, 'region')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_overlap_translation(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Features of a translation (domains, variants). Official GET /overlap/translation/:id;
    upstream defaults retained.
    """
    query = _params(parameters, {"db_type", "type", "content-type", "feature"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/overlap/translation/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_phenotype_accession(
    species: str,
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phenotype annotations by ontology accession. Official GET
    /phenotype/accession/:species/:accession; upstream defaults retained.
    """
    query = _params(
        parameters,
        {
            "trait",
            "include_associated",
            "content-type",
            "include_pubmed_id",
            "source",
            "include_review_status",
        },
    )
    return _request(
        build_url(
            "https://rest.ensembl.org",
            (
                f"/phenotype/accession/{_segment(species, 'species')}"
                f"/{_segment(accession, 'accession')}"
            ),
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_phenotype_gene(
    species: str,
    gene: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phenotype annotations of a gene. Official GET /phenotype/gene/:species/:gene; upstream
    defaults retained.
    """
    query = _params(
        parameters,
        {
            "trait",
            "include_associated",
            "content-type",
            "include_pubmed_id",
            "source",
            "include_review_status",
        },
    )
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/phenotype/gene/{_segment(species, 'species')}/{_segment(gene, 'gene')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_phenotype_term(
    species: str,
    term: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phenotype annotations by ontology term. Official GET /phenotype/term/:species/:term;
    upstream defaults retained.
    """
    query = _params(
        parameters,
        {
            "trait",
            "include_associated",
            "content-type",
            "include_pubmed_id",
            "source",
            "include_review_status",
        },
    )
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/phenotype/term/{_segment(species, 'species')}/{_segment(term, 'term')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_vep_region(
    species: str,
    region: str,
    allele: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variant consequences by genomic region and allele. Official GET
    /vep/:species/region/:region/:allele/; upstream defaults retained.
    """
    query = _params(
        parameters,
        {
            "af_gnomade",
            "af",
            "af_1kg",
            "species",
            "numbers",
            "conservation",
            "af_gnomadg",
            "xref_refseq",
            "distance",
            "failed",
            "transcript_version",
            "content-type",
            "domains",
            "variant_class",
            "pick_allele",
            "everything",
            "symbol",
            "pick",
            "hgvs",
            "canonical",
            "protein",
            "minimal",
            "mane",
        },
    )
    return _request(
        build_url(
            "https://rest.ensembl.org",
            (
                f"/vep/{_segment(species, 'species')}/region/{_segment(region, 'region')}"
                f"/{_segment(allele, 'allele')}/"
            ),
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_variant_recoder(
    species: str,
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variant ID, HGVS and SPDI recoding. Official GET /variant_recoder/:species/:id; upstream
    defaults retained.
    """
    query = _params(parameters, {"vcf_string", "content-type", "fields", "var_synonyms"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/variant_recoder/{_segment(species, 'species')}/{_segment(id, 'id')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_variation_publication(
    species: str,
    publication: str,
    *,
    publication_type: str = "pmcid",
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variants in a publication by PMID or PMCID. Official GET
    /variation/:species/pmcid/:pmcid; upstream defaults retained.
    """
    query = _params(
        parameters,
        {
            "pops",
            "genotypes",
            "content-type",
            "population_genotypes",
            "phenotypes",
            "genotyping_chips",
            "pubmed_id",
        },
    )
    if publication_type not in {"pmid", "pmcid"}:
        raise ConfigurationError("publication_type must be pmid or pmcid.")
    return _request(
        build_url(
            "https://rest.ensembl.org",
            (
                f"/variation/{_segment(species, 'species')}/{publication_type}"
                f"/{_segment(publication, 'publication')}"
            ),
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_beacon_query(
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH Beacon allele presence query. Official GET /ga4gh/beacon/query; upstream defaults
    retained.
    """
    query = _params(
        parameters,
        {
            "content-type",
            "includeDatasetResponses",
            "referenceBases",
            "start",
            "alternateBases",
            "datasetIds",
            "assemblyId",
            "referenceName",
        },
    )
    return _request(
        build_url("https://rest.ensembl.org", "/ga4gh/beacon/query", query.items()),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_features_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH feature record by ID. Official GET /ga4gh/features/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/features/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_features_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH sequence features in a region. Official POST /ga4gh/features/search; returns the
    service page token.
    """
    body = _params(
        query,
        {
            "end",
            "featureTypes",
            "start",
            "pageToken",
            "parentId",
            "pageSize",
            "featureSetId",
            "referenceName",
        },
    )
    return _request(
        "https://rest.ensembl.org/ga4gh/features/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_callsets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH genotype callsets search. Official POST /ga4gh/callsets/search; returns the
    service page token.
    """
    body = _params(query, {"pageToken", "pageSize", "name", "variantSetId"})
    return _request(
        "https://rest.ensembl.org/ga4gh/callsets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_callsets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH callset by ID. Official GET /ga4gh/callsets/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/callsets/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_datasets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH datasets search. Official POST /ga4gh/datasets/search; returns the service page
    token.
    """
    body = _params(query, {"pageToken", "pageSize"})
    return _request(
        "https://rest.ensembl.org/ga4gh/datasets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_datasets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH dataset by ID. Official GET /ga4gh/datasets/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/datasets/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_featuresets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH feature sets search. Official POST /ga4gh/featuresets/search; returns the service
    page token.
    """
    body = _params(query, {"datasetId", "pageToken", "pageSize"})
    return _request(
        "https://rest.ensembl.org/ga4gh/featuresets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_featuresets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH feature set by ID. Official GET /ga4gh/featuresets/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/featuresets/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variants_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant record by ID. Official GET /ga4gh/variants/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/variants/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variantannotations_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant annotations in a region. Official POST /ga4gh/variantannotations/search;
    returns the service page token.
    """
    body = _params(
        query,
        {
            "end",
            "variantAnnotationSetId",
            "start",
            "pageToken",
            "referenceId",
            "pageSize",
            "effects",
            "referenceName",
        },
    )
    return _request(
        "https://rest.ensembl.org/ga4gh/variantannotations/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variants_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant calls in a region. Official POST /ga4gh/variants/search; returns the
    service page token.
    """
    body = _params(
        query,
        {"end", "pageToken", "callSetIds", "pageSize", "referenceName", "variantSetId", "start"},
    )
    return _request(
        "https://rest.ensembl.org/ga4gh/variants/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variantsets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant sets search. Official POST /ga4gh/variantsets/search; returns the service
    page token.
    """
    body = _params(query, {"datasetId", "pageToken", "pageSize"})
    return _request(
        "https://rest.ensembl.org/ga4gh/variantsets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variantsets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant set by ID. Official GET /ga4gh/variantsets/:id; upstream defaults retained."""
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/variantsets/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_references_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH reference sequences search. Official POST /ga4gh/references/search; returns the
    service page token.
    """
    body = _params(query, {"pageToken", "referenceSetId", "pageSize", "accessions", "md5checksums"})
    return _request(
        "https://rest.ensembl.org/ga4gh/references/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_references_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH reference sequence by ID. Official GET /ga4gh/references/:id; upstream defaults
    retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/references/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_referencesets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH reference sets search. Official POST /ga4gh/referencesets/search; returns the
    service page token.
    """
    body = _params(query, {"pageToken", "accessions", "pageSize", "assemblyId", "md5checksums"})
    return _request(
        "https://rest.ensembl.org/ga4gh/referencesets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_referencesets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH reference set by ID. Official GET /ga4gh/referencesets/:id; upstream defaults
    retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org", f"/ga4gh/referencesets/{_segment(id, 'id')}", query.items()
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variantannotationsets_search(
    query: Mapping[str, Any],
    *,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant annotation sets search. Official POST /ga4gh/variantannotationsets/search;
    returns the service page token.
    """
    body = _params(query, {"pageToken", "pageSize", "variantSetId"})
    return _request(
        "https://rest.ensembl.org/ga4gh/variantannotationsets/search",
        body=body,
        timeout=timeout,
        max_bytes=max_bytes,
    )


def ensembl_ga4gh_variantannotationsets_record(
    id: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """GA4GH variant annotation set by ID. Official GET /ga4gh/variantannotationsets/:id;
    upstream defaults retained.
    """
    query = _params(parameters, {"content-type"})
    return _request(
        build_url(
            "https://rest.ensembl.org",
            f"/ga4gh/variantannotationsets/{_segment(id, 'id')}",
            query.items(),
        ),
        timeout=timeout,
        max_bytes=max_bytes,
    )


def datasets_genome_report(
    identifier: str,
    *,
    by: str = "accession",
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome assembly dataset report by identifier (accession, WGS, assembly name). NCBI
    Datasets v2; one page, preserving next_page_token.
    """
    routes = {
        "accession": "/genome/accession/{accessions}/dataset_report",
        "assembly_name": "/genome/assembly_name/{assembly_names}/dataset_report",
        "wgs": "/genome/wgs/{wgs_accessions}/dataset_report",
    }
    if by not in routes:
        raise ConfigurationError(f"by must be one of {sorted(routes)}.")
    path = routes[by]
    import re

    for field in re.findall(r"{(\w+)}", path):
        value = identifier
        path = path.replace("{" + field + "}", _segment(value, field))
    query = _params(
        parameters,
        {
            "sort.direction",
            "filters.assembly_source",
            "filters.assembly_level",
            "tax_exact_match",
            "include_tabular_header",
            "page_token",
            "filters.has_annotation",
            "filters.is_ictv_exemplar",
            "returned_content",
            "page_size",
            "table_fields",
            "filters.search_text",
            "filters.assembly_version",
            "filters.is_type_material",
            "filters.exclude_atypical",
            "filters.first_release_date",
            "filters.exclude_multi_isolate",
            "filters.last_release_date",
            "filters.exclude_paired_reports",
            "sort.field",
            "filters.reference_only",
            "filters.is_metagenome_derived",
            "filters.type_material_category",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_download_summary(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome data package download summary (preview) by accession. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/genome/accession/{accessions}/download_summary"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(parameters, {"chromosomes", "include_annotation_type"})
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_sequence_accession_sequence_assemblies(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome assembly accession for a nucleotide sequence accession. NCBI Datasets v2; one
    page, preserving next_page_token.
    """
    path = "/genome/sequence_accession/{accession}/sequence_assemblies"
    path = path.replace("{accession}", _segment(accession, "accession"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_links(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome assembly external resource links by accession. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/genome/accession/{accessions}/links"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(parameters, {"assembly_link_types"})
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_taxon_checkm_histogram(
    species_taxon: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """CheckM quality histogram of genome assemblies by species taxon. NCBI Datasets v2; one
    page, preserving next_page_token.
    """
    path = "/genome/taxon/{species_taxon}/checkm_histogram"
    path = path.replace("{species_taxon}", _segment(species_taxon, "species_taxon"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_check(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome assembly accession validity check. NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/genome/accession/{accessions}/check"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_annotation_report(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome annotation report (gene records) by accession. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/genome/accession/{accession}/annotation_report"
    path = path.replace("{accession}", _segment(accession, "accession"))
    query = _params(
        parameters,
        {
            "sort.direction",
            "symbols",
            "page_size",
            "locations",
            "gene_types",
            "table_format",
            "sort.field",
            "annotation_ids",
            "include_tabular_header",
            "search_text",
            "page_token",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_annotation_summary(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome annotation report summary (facets) by accession. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/genome/accession/{accession}/annotation_summary"
    path = path.replace("{accession}", _segment(accession, "accession"))
    query = _params(parameters, {"sort.field", "sort.direction"})
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_annotation_report_download_summary(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome annotation data package download summary (preview) by accession. NCBI Datasets
    v2; one page, preserving next_page_token.
    """
    path = "/genome/accession/{accession}/annotation_report/download_summary"
    path = path.replace("{accession}", _segment(accession, "accession"))
    query = _params(
        parameters,
        {
            "sort.direction",
            "symbols",
            "locations",
            "gene_types",
            "include_annotation_type",
            "sort.field",
            "annotation_ids",
            "search_text",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_genome_accession_revision_history(
    accession: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome assembly revision history by accession. NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/genome/accession/{accession}/revision_history"
    path = path.replace("{accession}", _segment(accession, "accession"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_taxon_dataset_report(
    taxon: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene data report by taxon (search). NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/gene/taxon/{taxon}/dataset_report"
    path = path.replace("{taxon}", _segment(taxon, "taxon"))
    query = _params(
        parameters,
        {
            "returned_content",
            "page_size",
            "tax_search_subtree",
            "table_fields",
            "query",
            "table_format",
            "types",
            "include_tabular_header",
            "page_token",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_product_report(
    identifier: str,
    *,
    by: str = "id",
    taxon: str | None = None,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene product report by GeneID, RefSeq accession, symbol+taxon or locus tag. NCBI
    Datasets v2; one page, preserving next_page_token.
    """
    routes = {
        "id": "/gene/id/{gene_ids}/product_report",
        "accession": "/gene/accession/{accessions}/product_report",
        "locus_tag": "/gene/locus_tag/{locus_tags}/product_report",
        "symbol": "/gene/symbol/{symbols}/taxon/{taxon}/product_report",
    }
    if by not in routes:
        raise ConfigurationError(f"by must be one of {sorted(routes)}.")
    path = routes[by]
    import re

    for field in re.findall(r"{(\w+)}", path):
        value = taxon if field == "taxon" else identifier
        path = path.replace("{" + field + "}", _segment(value, field))
    query = _params(
        parameters,
        {
            "sort.direction",
            "tax_search_subtree",
            "sort.field",
            "include_tabular_header",
            "page_token",
            "accession_filter",
            "page_size",
            "table_fields",
            "query",
            "table_format",
            "types",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_taxon_product_report(
    taxon: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene product report by taxon (search). NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/gene/taxon/{taxon}/product_report"
    path = path.replace("{taxon}", _segment(taxon, "taxon"))
    query = _params(
        parameters,
        {
            "page_size",
            "tax_search_subtree",
            "table_fields",
            "query",
            "table_format",
            "types",
            "include_tabular_header",
            "page_token",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_id_download_summary(
    gene_ids: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene data package download summary (preview) by GeneID. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/gene/id/{gene_ids}/download_summary"
    path = path.replace("{gene_ids}", _segment(gene_ids, "gene_ids"))
    query = _params(
        parameters,
        {
            "returned_content",
            "include_annotation_type",
            "aux_report",
            "tabular_reports",
            "fasta_filter",
            "accession_filter",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_taxon_counts(
    taxon: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene counts by taxon. NCBI Datasets v2; one page, preserving next_page_token."""
    path = "/gene/taxon/{taxon}/counts"
    path = path.replace("{taxon}", _segment(taxon, "taxon"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_taxon_annotation_chromosome_summary(
    taxon: str,
    annotation_name: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene counts per chromosome by taxon and annotation name. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/gene/taxon/{taxon}/annotation/{annotation_name}/chromosome_summary"
    path = path.replace("{taxon}", _segment(taxon, "taxon"))
    path = path.replace("{annotation_name}", _segment(annotation_name, "annotation_name"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_id_expression_report(
    gene_ids: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene expression report by GeneID. NCBI Datasets v2; one page, preserving next_page_token."""
    path = "/gene/id/{gene_ids}/expression_report"
    path = path.replace("{gene_ids}", _segment(gene_ids, "gene_ids"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_gene_id_links(
    gene_ids: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene external resource links by GeneID. NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/gene/id/{gene_ids}/links"
    path = path.replace("{gene_ids}", _segment(gene_ids, "gene_ids"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_protein_accession_download(
    accessions: str,
    *,
    output_path: str,
    overwrite: bool = False,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prokaryote gene data package download by RefSeq protein accession. NCBI Datasets v2; one
    page, preserving next_page_token.
    """
    path = "/protein/accession/{accessions}/download"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(
        parameters, {"filename", "taxon", "include_annotation_type", "gene_flank_config.length"}
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _download(
        url, output_path, overwrite=overwrite, verbose=verbose, timeout=timeout, max_bytes=max_bytes
    )


def datasets_organelle_accession_download(
    accessions: str,
    *,
    output_path: str,
    overwrite: bool = False,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Organelle genome data package download by nucleotide accession. NCBI Datasets v2; one
    page, preserving next_page_token.
    """
    path = "/organelle/accession/{accessions}/download"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(parameters, {"filename", "include_annotation_type"})
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _download(
        url, output_path, overwrite=overwrite, verbose=verbose, timeout=timeout, max_bytes=max_bytes
    )


def datasets_organelle_accessions_dataset_report(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Organelle genome data report by nucleotide accession. NCBI Datasets v2; one page,
    preserving next_page_token.
    """
    path = "/organelle/accessions/{accessions}/dataset_report"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(
        parameters,
        {
            "returned_content",
            "sort.direction",
            "first_release_date",
            "last_release_date",
            "table_format",
            "sort.field",
            "organelle_types",
            "include_tabular_header",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_organelle_taxon_dataset_report(
    taxons: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Organelle genome data report by taxon (search). NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/organelle/taxon/{taxons}/dataset_report"
    path = path.replace("{taxons}", _segment(taxons, "taxons"))
    query = _params(
        parameters,
        {
            "returned_content",
            "sort.direction",
            "page_size",
            "first_release_date",
            "last_release_date",
            "table_format",
            "organelle_types",
            "sort.field",
            "tax_exact_match",
            "include_tabular_header",
            "page_token",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_virus_download_summary(
    identifier: str,
    *,
    by: str = "taxon",
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Virus genome data package download summary (preview) by taxon or SARS-CoV-2 protein.
    NCBI Datasets v2; one page, preserving next_page_token.
    """
    routes = {
        "taxon": "/virus/taxon/sars2/protein/{proteins}",
        "sars2_protein": "/virus/taxon/sars2/protein/{proteins}",
    }
    if by not in routes:
        raise ConfigurationError(f"by must be one of {sorted(routes)}.")
    path = routes[by]
    import re

    for field in re.findall(r"{(\w+)}", path):
        value = identifier
        path = path.replace("{" + field + "}", _segment(value, field))
    query = _params(
        parameters,
        {
            "geo_location",
            "usa_state",
            "include_sequence",
            "aux_report",
            "pangolin_classification",
            "host",
            "refseq_only",
            "updated_since",
            "complete_only",
            "annotated_only",
            "released_since",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_virus_taxon_annotation_report(
    taxon: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Virus annotation report by taxon (search). NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/virus/taxon/{taxon}/annotation_report"
    path = path.replace("{taxon}", _segment(taxon, "taxon"))
    query = _params(
        parameters,
        {
            "filter.annotated_only",
            "filter.updated_since",
            "filter.pangolin_classification",
            "page_size",
            "filter.refseq_only",
            "filter.released_since",
            "table_fields",
            "page_token",
            "filter.host",
            "filter.usa_state",
            "filter.geo_location",
            "filter.complete_only",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_virus_accession_annotation_report(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Virus annotation report by nucleotide accession. NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/virus/accession/{accessions}/annotation_report"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(
        parameters,
        {
            "filter.annotated_only",
            "filter.updated_since",
            "filter.pangolin_classification",
            "page_size",
            "filter.refseq_only",
            "filter.released_since",
            "table_fields",
            "page_token",
            "filter.host",
            "filter.usa_state",
            "filter.geo_location",
            "filter.complete_only",
        },
    )
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


def datasets_virus_accession_check(
    accessions: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Virus genome nucleotide accession validity check. NCBI Datasets v2; one page, preserving
    next_page_token.
    """
    path = "/virus/accession/{accessions}/check"
    path = path.replace("{accessions}", _segment(accessions, "accessions"))
    query = _params(parameters, set())
    url = build_url("https://api.ncbi.nlm.nih.gov/datasets/v2", path, query.items())
    return _request(url, timeout=timeout, max_bytes=max_bytes)


__all__ = [
    "datasets_gene_id_download_summary",
    "datasets_gene_id_expression_report",
    "datasets_gene_id_links",
    "datasets_gene_product_report",
    "datasets_gene_taxon_annotation_chromosome_summary",
    "datasets_gene_taxon_counts",
    "datasets_gene_taxon_dataset_report",
    "datasets_gene_taxon_product_report",
    "datasets_genome_accession_annotation_report",
    "datasets_genome_accession_annotation_report_download_summary",
    "datasets_genome_accession_annotation_summary",
    "datasets_genome_accession_check",
    "datasets_genome_accession_download_summary",
    "datasets_genome_accession_links",
    "datasets_genome_accession_revision_history",
    "datasets_genome_report",
    "datasets_genome_sequence_accession_sequence_assemblies",
    "datasets_genome_taxon_checkm_histogram",
    "datasets_organelle_accession_download",
    "datasets_organelle_accessions_dataset_report",
    "datasets_organelle_taxon_dataset_report",
    "datasets_protein_accession_download",
    "datasets_virus_accession_annotation_report",
    "datasets_virus_accession_check",
    "datasets_virus_download_summary",
    "datasets_virus_taxon_annotation_report",
    "ensembl_ga4gh_beacon_query",
    "ensembl_ga4gh_callsets_record",
    "ensembl_ga4gh_callsets_search",
    "ensembl_ga4gh_datasets_record",
    "ensembl_ga4gh_datasets_search",
    "ensembl_ga4gh_features_record",
    "ensembl_ga4gh_features_search",
    "ensembl_ga4gh_featuresets_record",
    "ensembl_ga4gh_featuresets_search",
    "ensembl_ga4gh_references_record",
    "ensembl_ga4gh_references_search",
    "ensembl_ga4gh_referencesets_record",
    "ensembl_ga4gh_referencesets_search",
    "ensembl_ga4gh_variantannotations_search",
    "ensembl_ga4gh_variantannotationsets_record",
    "ensembl_ga4gh_variantannotationsets_search",
    "ensembl_ga4gh_variants_record",
    "ensembl_ga4gh_variants_search",
    "ensembl_ga4gh_variantsets_record",
    "ensembl_ga4gh_variantsets_search",
    "ensembl_info_analysis",
    "ensembl_info_assembly",
    "ensembl_info_biotypes",
    "ensembl_info_biotypes_groups",
    "ensembl_info_biotypes_name",
    "ensembl_info_compara_methods",
    "ensembl_info_compara_species_sets",
    "ensembl_info_external_dbs",
    "ensembl_info_genomes_accession",
    "ensembl_info_genomes_assembly",
    "ensembl_info_genomes_division",
    "ensembl_info_genomes_taxonomy",
    "ensembl_info_variation",
    "ensembl_info_variation_consequence_types",
    "ensembl_info_variation_populations",
    "ensembl_info_variation_populations_list",
    "ensembl_map_cds",
    "ensembl_map_translation",
    "ensembl_overlap_translation",
    "ensembl_phenotype_accession",
    "ensembl_phenotype_gene",
    "ensembl_phenotype_term",
    "ensembl_variant_recoder",
    "ensembl_variation_publication",
    "ensembl_vep_region",
]
