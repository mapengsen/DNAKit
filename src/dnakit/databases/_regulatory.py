"""Gene-locus peak queries; never substitute an antigen search for locus overlap."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from dnakit.exceptions import BackendExecutionError, ConfigurationError


def gene_peaks(values: dict[str, Any]) -> dict[str, Any]:
    """Resolve a locus, select an assembled track, and optionally intersect supplied BED data.

    Protocol: https://chip-atlas.org/agents (POST /download).
    A download link is a retrieval plan, not a collection of observed binding peaks.
    """
    import requests

    from ._runtime import call

    if values.get("experiment_id"):
        raise ConfigurationError("Use experiment_id for a track or gene for a locus query.")
    genome = values.get("genome", "hg38")
    if genome != "hg38":
        raise ConfigurationError(
            "Gene-symbol resolution currently requires hg38; use an explicit experiment track "
            "for other assemblies. Coordinates are never silently lifted between assemblies."
        )
    gene = call("ensembl_lookup_gene", {"gene_id": values["gene"], "species": "homo_sapiens"})
    if isinstance(gene, dict) and isinstance(gene.get("data"), dict):
        gene = gene["data"]
    if not isinstance(gene, dict) or gene.get("assembly_name") != "GRCh38":
        raise BackendExecutionError("Ensembl must resolve a unique GRCh38 gene locus.")
    flank = values.get("flank_bp", 5000)
    chrom = "chr" + str(gene["seq_region_name"]).removeprefix("chr")
    start, end = max(0, int(gene["start"]) - 1 - flank), int(gene["end"]) + flank
    region = {"chrom": chrom, "start": start, "end": end, "coordinate_system": "BED0"}
    context = {"gene": values["gene"], "assembly": genome, "region": region}
    if values.get("peak_file"):
        source = Path(values["peak_file"]).expanduser().resolve(strict=True)
        if source.stat().st_size > values.get("max_file_bytes", 100_000_000):
            raise ConfigurationError("Peak BED exceeds max_file_bytes.")
        peaks: list[dict[str, Any]] = []
        matches = 0
        with source.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 3:
                    raise ConfigurationError(f"Invalid BED row {line_number}.")
                left, right = int(fields[1]), int(fields[2])
                if left < 0 or right <= left:
                    raise ConfigurationError(f"Invalid BED interval at row {line_number}.")
                if fields[0] == chrom and left < end and right > start:
                    matches += 1
                    if len(peaks) < values.get("limit", 1000):
                        peaks.append(
                            {"chrom": chrom, "start": left, "end": right, "fields": fields[3:]}
                        )
        return {
            "status": "success",
            "data": peaks,
            **context,
            "matched": matches,
            "truncated": matches > len(peaks),
            "source": str(source),
            "scope": (
                "overlap with the supplied BED; assembly and TF-track provenance "
                "are caller-supplied"
            ),
        }
    experiment = values.get("experiment_type", "TF")
    condition = {
        "genome": genome,
        "agClass": "TFs and others" if experiment == "TF" else experiment,
        "agSubClass": "All antigens",
        "clClass": "All cell types",
        "clSubClass": "All cell types",
        "qval": str(int(values.get("threshold", "05"))),
    }
    response = requests.post(
        "https://chip-atlas.org/download", json={"condition": condition}, timeout=30
    )
    response.raise_for_status()
    url = response.json().get("url")
    if not isinstance(url, str) or urlsplit(url).hostname != "chip-atlas.dbcls.jp":
        raise BackendExecutionError("ChIP-Atlas did not provide an assembled peak track URL.")
    return {
        "status": "requires_data",
        "data": {"url": url, "track_condition": condition},
        **context,
        "analysis_completed": False,
        "next_step": (
            "Download the selected BED, verify its assembly and track, "
            "then call again with peak_file."
        ),
        "scope": "retrieval plan; no peaks have been fetched or motif disruption evaluated",
    }
