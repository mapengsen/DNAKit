"""Ensembl REST adapters for coordinates, annotation, variants, and compara."""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from typing import Literal, cast
from urllib.parse import quote

from dnakit.exceptions import ConfigurationError, QueryError

from ._http import (
    build_url,
    limited_records,
    mapping_records,
    redact_url,
    request_json,
    require_text,
)
from ._shared import adapter_provenance, query_values, resolved_config
from .models import QueryProgress, QueryResult, SearchConfig

DEFAULT_ENSEMBL_REST_URL = "https://rest.ensembl.org"
_REGION = re.compile(r"^(?P<name>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)$")
_ENSEMBL_ID = re.compile(r"^ENS[A-Z0-9]*[GTP]\d+(?:\.\d+)?$", flags=re.IGNORECASE)
_FEATURES = frozenset(
    {
        "band",
        "cds",
        "constrained",
        "exon",
        "gene",
        "mane",
        "misc",
        "motif",
        "regulatory",
        "repeat",
        "simple",
        "somatic_structural_variation",
        "somatic_variation",
        "structural_variation",
        "transcript",
        "variation",
    }
)


def _species(value: str) -> str:
    return require_text(value, "species", max_length=256)


def _bounded_non_negative(value: int, name: str, *, maximum: int = 5_000_000) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= maximum:
        raise ConfigurationError(
            f"{name} must be an integer in [0, {maximum}].", code="INVALID_QUERY"
        )
    return value


def _parse_region(region: str, *, max_span: int) -> tuple[str, int, int]:
    value = require_text(region, "region", max_length=512)
    match = _REGION.fullmatch(value)
    if match is None:
        raise ConfigurationError(
            "region must use 'chromosome:start-end' with 0-based half-open coordinates.",
            code="INVALID_REGION_QUERY",
        )
    name = match.group("name")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end <= start:
        raise ConfigurationError(
            "region end must be greater than start.", code="INVALID_REGION_QUERY"
        )
    if end - start > max_span:
        raise ConfigurationError(
            "region exceeds this Ensembl endpoint's span limit.",
            code="REGION_QUERY_LIMIT",
            context={"span": end - start, "max_span": max_span},
        )
    return name, start, end


def _ensembl_region(region: str, *, strand: int, max_span: int) -> tuple[str, int, int, str]:
    if strand not in {-1, 1}:
        raise ConfigurationError("strand must be 1 or -1.", code="INVALID_REGION_QUERY")
    name, start, end = _parse_region(region, max_span=max_span)
    provider_region = f"{name}:{start + 1}..{end}:{strand}"
    return name, start, end, provider_region


def _response_records(
    payload: object, *, keys: Sequence[str] = ()
) -> tuple[dict[str, object], ...]:
    if isinstance(payload, list):
        return mapping_records(payload)
    if not isinstance(payload, Mapping):
        return ()
    for key in keys:
        records = mapping_records(payload, key=key)
        if records:
            return records
    return (dict(payload),)


def _with_zero_based(record: Mapping[str, object]) -> dict[str, object]:
    normalized = dict(record)
    start = record.get("start")
    end = record.get("end")
    if (
        isinstance(start, int)
        and not isinstance(start, bool)
        and isinstance(end, int)
        and not isinstance(end, bool)
        and start >= 1
        and end >= start
    ):
        normalized["start_0based"] = start - 1
        normalized["end_0based"] = end
    for key in ("mappings", "data", "Transcript", "Exon", "UTR"):
        nested = record.get(key)
        if isinstance(nested, list):
            normalized[key] = tuple(
                _with_zero_based(item) if isinstance(item, Mapping) else item for item in nested
            )
    for key in ("mapped", "original", "Translation"):
        nested_mapping = record.get(key)
        if isinstance(nested_mapping, Mapping):
            normalized[key] = _with_zero_based(nested_mapping)
    return normalized


def _result(
    query_type: str,
    url: str,
    payload: object,
    config: SearchConfig,
    *,
    record_keys: Sequence[str] = (),
    filters: Mapping[str, object] | None = None,
    metadata: Mapping[str, object] | None = None,
    normalize_coordinates: bool = False,
) -> QueryResult:
    records = _response_records(payload, keys=record_keys)
    if normalize_coordinates:
        records = tuple(_with_zero_based(record) for record in records)
    records = limited_records(records, config)
    return QueryResult(
        query_type,
        "Ensembl",
        redact_url(url),
        records,
        adapter_provenance(
            "Ensembl REST",
            citation_url="https://rest.ensembl.org/",
            filters=filters,
        ),
        total_count=len(records),
        metadata=metadata,
    )


def sequence(
    species: str,
    region: str | Sequence[str],
    *,
    strand: Literal[-1, 1] = 1,
    upstream: int = 0,
    downstream: int = 0,
    mask: Literal["hard", "soft"] | None = None,
    progress: Callable[[QueryProgress], None] | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Fetch one or more coordinate sequences using 0-based half-open regions."""

    resolved = resolved_config(config)
    organism = _species(species)
    regions = query_values(region, name="region", maximum=min(100, resolved.max_records))
    flank_5 = _bounded_non_negative(upstream, "upstream")
    flank_3 = _bounded_non_negative(downstream, "downstream")
    if mask not in {None, "hard", "soft"}:
        raise ConfigurationError("mask must be None, hard, or soft.", code="INVALID_QUERY")
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.", code="INVALID_QUERY")
    records: list[dict[str, object]] = []
    urls: list[str] = []
    for index, item in enumerate(regions, start=1):
        _, start, end, provider_region = _ensembl_region(item, strand=strand, max_span=10_000_000)
        url = build_url(
            api_base_url,
            f"/sequence/region/{quote(organism, safe='')}/{quote(provider_region, safe=':.-_')}",
            (("expand_5prime", flank_5), ("expand_3prime", flank_3), ("mask", mask)),
        )
        payload = request_json(url, resolved, provider="Ensembl")
        if not isinstance(payload, Mapping):
            raise QueryError(
                "Ensembl sequence returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
            )
        record = dict(payload)
        record.update(
            {
                "requested_region": item,
                "requested_start_0based": start,
                "requested_end_0based": end,
                "requested_strand": strand,
            }
        )
        records.append(record)
        urls.append(redact_url(url))
        if progress is not None:
            progress(QueryProgress(index, len(regions), item))
    records_tuple = limited_records(records, resolved)
    return QueryResult(
        "sequence",
        "Ensembl",
        urls[0] if len(urls) == 1 else f"{api_base_url.rstrip('/')}/sequence/region/[batch]",
        records_tuple,
        adapter_provenance(
            "Ensembl REST",
            citation_url="https://rest.ensembl.org/documentation/info/sequence_region",
            filters={
                "species": organism,
                "regions": regions,
                "strand": strand,
                "upstream": flank_5,
                "downstream": flank_3,
                "mask": mask,
            },
        ),
        total_count=len(records_tuple),
        metadata={
            "input_coordinate_system": "0-based half-open",
            "provider_coordinate_system": "1-based inclusive",
            "request_urls": tuple(urls),
        },
    )


def sequence_by_id(
    identifier: str,
    *,
    sequence_type: Literal["genomic", "cds", "cdna", "protein"] = "genomic",
    species: str | None = None,
    mask: Literal["hard", "soft"] | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Fetch genomic, cDNA, CDS, or protein sequence for an Ensembl stable ID."""

    resolved = resolved_config(config)
    stable_id = require_text(identifier, "identifier", max_length=256)
    if sequence_type not in {"genomic", "cds", "cdna", "protein"}:
        raise ConfigurationError("Unsupported sequence_type.", code="INVALID_QUERY")
    if mask not in {None, "hard", "soft"}:
        raise ConfigurationError("mask must be None, hard, or soft.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        f"/sequence/id/{quote(stable_id, safe='._-')}",
        (
            ("type", sequence_type),
            ("species", None if species is None else _species(species)),
            ("mask", mask),
        ),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "sequence",
        url,
        payload,
        resolved,
        filters={"identifier": stable_id, "sequence_type": sequence_type, "species": species},
        metadata={"sequence_type": sequence_type},
    )


def lookup(
    identifier: str,
    *,
    species: str | None = None,
    expand: bool = True,
    mane: bool = True,
    phenotypes: bool = False,
    utr: bool = True,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Look up one Ensembl gene/transcript/protein ID and connected features."""

    resolved = resolved_config(config)
    stable_id = require_text(identifier, "identifier", max_length=256)
    for name, value in (
        ("expand", expand),
        ("mane", mane),
        ("phenotypes", phenotypes),
        ("utr", utr),
    ):
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be boolean.", code="INVALID_QUERY")
    if not expand and (mane or utr):
        raise ConfigurationError("mane and utr require expand=True.", code="INVALID_QUERY")
    organism = None if species is None else _species(species)
    url = build_url(
        api_base_url,
        f"/lookup/id/{quote(stable_id, safe='._-')}",
        (
            ("species", organism),
            ("expand", expand),
            ("format", "full"),
            ("mane", mane),
            ("phenotypes", phenotypes),
            ("utr", utr),
        ),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "feature",
        url,
        payload,
        resolved,
        filters={
            "identifier": stable_id,
            "species": organism,
            "expand": expand,
            "mane": mane,
            "phenotypes": phenotypes,
            "utr": utr,
        },
        normalize_coordinates=True,
    )


def transcripts(
    query: str,
    *,
    species: str = "human",
    by: Literal["auto", "id", "symbol"] = "auto",
    expand: bool = True,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query transcripts, exons, translations, and canonical status by ID or symbol."""

    resolved = resolved_config(config)
    value = require_text(query, "query", max_length=256)
    organism = _species(species)
    if not isinstance(expand, bool):
        raise ConfigurationError("expand must be boolean.", code="INVALID_QUERY")
    if by == "auto":
        by = "id" if _ENSEMBL_ID.fullmatch(value) else "symbol"
    if by == "id":
        result = lookup(
            value,
            species=organism,
            expand=expand,
            mane=expand,
            utr=expand,
            config=resolved,
            api_base_url=api_base_url,
        )
    elif by == "symbol":
        url = build_url(
            api_base_url,
            f"/lookup/symbol/{quote(organism, safe='')}/{quote(value, safe='._-')}",
            (("expand", expand), ("format", "full")),
        )
        payload = request_json(url, resolved, provider="Ensembl")
        result = _result(
            "feature",
            url,
            payload,
            resolved,
            filters={"query": value, "species": organism, "by": by, "expand": expand},
            normalize_coordinates=True,
        )
    else:
        raise ConfigurationError(
            "transcripts by must be auto, id, or symbol.", code="INVALID_QUERY"
        )
    return QueryResult(
        "transcripts",
        result.provider,
        result.request_url,
        result.records,
        result.provenance,
        total_count=result.total_count,
        metadata={**dict(result.metadata), "query": value, "by": by},
    )


def annotation(
    species: str,
    region: str,
    *,
    features: Sequence[str] = ("gene", "transcript", "exon"),
    strand: Literal[-1, 1] = 1,
    biotype: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query genes, transcripts, repeats, variants, and regulation over a region."""

    resolved = resolved_config(config)
    organism = _species(species)
    _, start, end, provider_region = _ensembl_region(region, strand=strand, max_span=5_000_000)
    feature_values = tuple(require_text(item, "feature", max_length=64) for item in features)
    if not feature_values or any(item not in _FEATURES for item in feature_values):
        raise ConfigurationError(
            "features contains an unsupported Ensembl overlap feature.",
            code="INVALID_QUERY",
            context={"allowed": sorted(_FEATURES)},
        )
    url = build_url(
        api_base_url,
        f"/overlap/region/{quote(organism, safe='')}/{quote(provider_region, safe=':.-_')}",
        (("feature", feature_values), ("biotype", biotype)),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "annotation",
        url,
        payload,
        resolved,
        filters={"species": organism, "region": region, "features": feature_values},
        metadata={
            "input_coordinate_system": "0-based half-open",
            "provider_coordinate_system": "1-based inclusive",
            "requested_start_0based": start,
            "requested_end_0based": end,
        },
        normalize_coordinates=True,
    )


def nearby(
    species: str,
    region: str,
    *,
    feature: Literal["gene", "transcript", "repeat", "variation", "regulatory"] = "gene",
    distance: int = 100_000,
    limit: int = 20,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Find nearest Ensembl features within a bounded flanking window."""

    resolved = resolved_config(config)
    flank = _bounded_non_negative(distance, "distance", maximum=2_499_999)
    if (
        isinstance(limit, bool)
        or not isinstance(limit, int)
        or not 1 <= limit <= resolved.max_records
    ):
        raise ConfigurationError(
            "limit exceeds SearchConfig max_records.", code="QUERY_RECORD_LIMIT"
        )
    chromosome, start, end = _parse_region(region, max_span=5_000_000)
    expanded_start = max(0, start - flank)
    expanded_end = end + flank
    if expanded_end - expanded_start > 5_000_000:
        raise ConfigurationError(
            "region plus distance exceeds Ensembl's 5 Mb overlap limit.",
            code="REGION_QUERY_LIMIT",
        )
    overlap = annotation(
        species,
        f"{chromosome}:{expanded_start}-{expanded_end}",
        features=(feature,),
        config=resolved,
        api_base_url=api_base_url,
    )
    ranked: list[dict[str, object]] = []
    for record in overlap.records:
        item = cast(dict[str, object], dict(record))
        item_start = item.get("start_0based")
        item_end = item.get("end_0based")
        if not isinstance(item_start, int) or not isinstance(item_end, int):
            continue
        if item_end <= start:
            gap = start - item_end
        elif item_start >= end:
            gap = item_start - end
        else:
            gap = 0
        item["distance_bp"] = gap
        ranked.append(item)
    ranked.sort(
        key=lambda item: (
            item["distance_bp"] if isinstance(item.get("distance_bp"), int) else 2**63,
            str(item.get("id", "")),
        )
    )
    selected = tuple(ranked[:limit])
    return QueryResult(
        "nearby",
        "Ensembl",
        overlap.request_url,
        selected,
        overlap.provenance,
        total_count=len(ranked),
        metadata={
            "query_region": region,
            "distance": flank,
            "feature": feature,
            "input_coordinate_system": "0-based half-open",
        },
    )


def variant(
    query: str,
    *,
    species: str = "human",
    mode: Literal["auto", "id", "region", "hgvs"] = "auto",
    include_frequencies: bool = True,
    include_genotypes: bool = False,
    include_phenotypes: bool = True,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query variant records, population frequencies, or VEP consequences."""

    resolved = resolved_config(config)
    organism = _species(species)
    value = require_text(query, "query", max_length=2_048)
    if mode == "auto":
        if _REGION.fullmatch(value):
            mode = "region"
        elif re.search(r":[cgmnpr]\.\S+", value, flags=re.IGNORECASE):
            mode = "hgvs"
        else:
            mode = "id"
    if mode == "region":
        _, _, _, provider_region = _ensembl_region(value, strand=1, max_span=5_000_000)
        url = build_url(
            api_base_url,
            f"/overlap/region/{quote(organism, safe='')}/{quote(provider_region, safe=':.-_')}",
            (("feature", ("variation", "structural_variation")),),
        )
        keys: Sequence[str] = ()
    elif mode == "hgvs":
        url = build_url(
            api_base_url,
            f"/vep/{quote(organism, safe='')}/hgvs/{quote(value, safe='._:-')}",
            (("hgvs", True), ("mane", True), ("variant_class", True)),
        )
        keys = ()
    elif mode == "id":
        url = build_url(
            api_base_url,
            f"/variation/{quote(organism, safe='')}/{quote(value, safe='._:-')}",
            (
                ("pops", include_frequencies),
                ("population_genotypes", include_frequencies),
                ("genotypes", include_genotypes),
                ("phenotypes", include_phenotypes),
            ),
        )
        keys = ()
    else:
        raise ConfigurationError("Unsupported variant query mode.", code="INVALID_QUERY")
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "variant",
        url,
        payload,
        resolved,
        record_keys=keys,
        filters={"query": value, "species": organism, "mode": mode},
        metadata={"mode": mode},
        normalize_coordinates=True,
    )


def regulation(
    species: str,
    region: str,
    *,
    include_motifs: bool = True,
    include_constrained: bool = False,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query Ensembl regulatory features, motifs, and optional constrained elements."""

    features = ["regulatory"]
    if include_motifs:
        features.append("motif")
    if include_constrained:
        features.append("constrained")
    result = annotation(
        species,
        region,
        features=tuple(features),
        config=config,
        api_base_url=api_base_url,
    )
    return QueryResult(
        "regulation",
        result.provider,
        result.request_url,
        result.records,
        result.provenance,
        total_count=result.total_count,
        metadata=dict(result.metadata),
    )


def homology(
    identifier: str,
    *,
    species: str = "human",
    by: Literal["auto", "id", "symbol"] = "auto",
    homology_type: Literal["orthologues", "paralogues", "projections", "all"] = "all",
    target_species: str | None = None,
    sequence_type: Literal["none", "cdna", "protein"] = "none",
    aligned: bool = False,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query orthologs, paralogs, projections, and gene-tree alignments."""

    resolved = resolved_config(config)
    value = require_text(identifier, "identifier", max_length=256)
    organism = _species(species)
    if by == "auto":
        by = "id" if _ENSEMBL_ID.fullmatch(value) else "symbol"
    if by == "id":
        path = f"/homology/id/{quote(organism, safe='')}/{quote(value, safe='._-')}"
    elif by == "symbol":
        path = f"/homology/symbol/{quote(organism, safe='')}/{quote(value, safe='._-')}"
    else:
        raise ConfigurationError("homology by must be auto, id, or symbol.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        path,
        (
            ("type", homology_type),
            ("target_species", target_species),
            ("sequence", sequence_type),
            ("aligned", aligned),
            ("cigar_line", aligned),
        ),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "homology",
        url,
        payload,
        resolved,
        record_keys=("data",),
        filters={
            "identifier": value,
            "species": organism,
            "by": by,
            "type": homology_type,
            "target_species": target_species,
        },
    )


def id_convert(
    identifier: str,
    *,
    species: str = "human",
    by: Literal["auto", "id", "symbol"] = "auto",
    all_levels: bool = True,
    external_db: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Convert Ensembl IDs/symbols through Ensembl external cross-references."""

    resolved = resolved_config(config)
    value = require_text(identifier, "identifier", max_length=256)
    if by == "auto":
        by = "id" if value.upper().startswith("ENS") else "symbol"
    if by == "id":
        path = f"/xrefs/id/{quote(value, safe='._-')}"
        params: tuple[tuple[str, object], ...] = (
            ("all_levels", all_levels),
            ("external_db", external_db),
        )
    elif by == "symbol":
        path = f"/xrefs/symbol/{quote(_species(species), safe='')}/{quote(value, safe='._-')}"
        params = (("external_db", external_db),)
    else:
        raise ConfigurationError("id_convert by must be auto, id, or symbol.", code="INVALID_QUERY")
    url = build_url(api_base_url, path, params)
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "id_conversion",
        url,
        payload,
        resolved,
        filters={"identifier": value, "species": species, "by": by, "external_db": external_db},
    )


def map_coordinates(
    species: str,
    region: str,
    *,
    source_assembly: str,
    target_assembly: str,
    strand: Literal[-1, 1] = 1,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Map a 0-based half-open region between two assemblies."""

    resolved = resolved_config(config)
    organism = _species(species)
    _, start, end, provider_region = _ensembl_region(region, strand=strand, max_span=5_000_000)
    source = require_text(source_assembly, "source_assembly", max_length=128)
    target = require_text(target_assembly, "target_assembly", max_length=128)
    url = build_url(
        api_base_url,
        (
            f"/map/{quote(organism, safe='')}/{quote(source, safe='._-')}/"
            f"{quote(provider_region, safe=':.-_')}/{quote(target, safe='._-')}"
        ),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "coordinate_mapping",
        url,
        payload,
        resolved,
        record_keys=("mappings",),
        filters={
            "species": organism,
            "region": region,
            "source_assembly": source,
            "target_assembly": target,
        },
        metadata={
            "input_coordinate_system": "0-based half-open",
            "provider_coordinate_system": "1-based inclusive",
            "requested_start_0based": start,
            "requested_end_0based": end,
        },
        normalize_coordinates=True,
    )


def comparative_alignment(
    species: str,
    region: str,
    *,
    strand: Literal[-1, 1] = 1,
    method: str = "EPO",
    species_set_group: str | None = "mammals",
    display_species: Sequence[str] = (),
    aligned: bool = True,
    mask: Literal["hard", "soft"] | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Retrieve bounded Ensembl Compara genomic alignment blocks."""

    resolved = resolved_config(config)
    organism = _species(species)
    _, start, end, provider_region = _ensembl_region(region, strand=strand, max_span=10_000_000)
    method_value = require_text(method, "method", max_length=64)
    displays = tuple(_species(item) for item in display_species)
    if mask not in {None, "hard", "soft"}:
        raise ConfigurationError("mask must be None, hard, or soft.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        f"/alignment/region/{quote(organism, safe='')}/{quote(provider_region, safe=':.-_')}",
        (
            ("method", method_value),
            ("species_set_group", species_set_group),
            ("display_species_set", displays),
            ("aligned", aligned),
            ("mask", mask),
        ),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "comparative_alignment",
        url,
        payload,
        resolved,
        filters={
            "species": organism,
            "region": region,
            "method": method_value,
            "species_set_group": species_set_group,
            "display_species": displays,
        },
        metadata={
            "input_coordinate_system": "0-based half-open",
            "provider_coordinate_system": "1-based inclusive",
            "requested_start_0based": start,
            "requested_end_0based": end,
        },
        normalize_coordinates=True,
    )


def genome_info(
    species: str,
    *,
    expand_sequences: bool = False,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Return Ensembl assembly, genebuild, release, and sequence metadata."""

    resolved = resolved_config(config)
    organism = _species(species)
    url = build_url(
        api_base_url,
        f"/info/genomes/{quote(organism, safe='._-')}",
        (("expand", expand_sequences),),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "database_version",
        url,
        payload,
        resolved,
        filters={"species": organism, "expand_sequences": expand_sequences},
    )


def chromosome(
    species: str,
    name: str,
    *,
    include_bands: bool = False,
    include_synonyms: bool = True,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> QueryResult:
    """Query chromosome, organelle, plasmid, contig, or scaffold metadata."""

    resolved = resolved_config(config)
    organism = _species(species)
    sequence_name = require_text(name, "name", max_length=256)
    url = build_url(
        api_base_url,
        f"/info/assembly/{quote(organism, safe='')}/{quote(sequence_name, safe='._-')}",
        (("bands", include_bands), ("synonyms", include_synonyms)),
    )
    payload = request_json(url, resolved, provider="Ensembl")
    return _result(
        "chromosome",
        url,
        payload,
        resolved,
        filters={"species": organism, "name": sequence_name},
    )


__all__ = [
    "DEFAULT_ENSEMBL_REST_URL",
    "annotation",
    "chromosome",
    "comparative_alignment",
    "genome_info",
    "homology",
    "id_convert",
    "lookup",
    "map_coordinates",
    "nearby",
    "regulation",
    "sequence",
    "sequence_by_id",
    "transcripts",
    "variant",
]
