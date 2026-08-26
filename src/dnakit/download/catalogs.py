"""Curated public ClinVar and GEO download locations."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

from dnakit.exceptions import ConfigurationError

from .files import dataset
from .models import DatasetDownloadResult, DownloadConfig, DownloadProgress, RemoteFile

_GEO_ACCESSION = re.compile(r"^(GSE|GDS|GPL|GSM)(\d+)$", flags=re.IGNORECASE)


def variants(
    output_dir: str | Path,
    *,
    source: str = "clinvar",
    format: str = "vcf",
    assembly: str = "GRCh38",
    include_index: bool = True,
    include_checksum_files: bool = True,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
) -> DatasetDownloadResult:
    """Download current public ClinVar VCF, XML, or variant-summary TSV."""

    source_value = source.casefold()
    if source_value not in {"clinvar", "dbsnp"}:
        raise ConfigurationError(
            "variants source must be clinvar or dbsnp; use dataset() for other archives."
        )
    kind = format.casefold()
    resources: list[RemoteFile]
    if source_value == "dbsnp":
        if kind != "vcf":
            raise ConfigurationError("dbSNP downloads currently use format='vcf'.")
        accessions = {
            "GRCh37": "GCF_000001405.25",
            "GRCh38": "GCF_000001405.40",
        }
        accession = accessions.get(assembly)
        if accession is None:
            raise ConfigurationError("dbSNP VCF assembly must be GRCh37 or GRCh38.")
        base = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF"
        resources = [RemoteFile(f"{base}/{accession}.gz", filename=f"dbsnp_{assembly}.vcf.gz")]
        if include_index:
            resources.append(
                RemoteFile(f"{base}/{accession}.gz.tbi", filename=f"dbsnp_{assembly}.vcf.gz.tbi")
            )
        if include_checksum_files:
            resources.append(
                RemoteFile(f"{base}/{accession}.gz.md5", filename=f"dbsnp_{assembly}.vcf.gz.md5")
            )
            if include_index:
                resources.append(
                    RemoteFile(
                        f"{base}/{accession}.gz.tbi.md5",
                        filename=f"dbsnp_{assembly}.vcf.gz.tbi.md5",
                    )
                )
    elif kind == "vcf":
        if assembly not in {"GRCh37", "GRCh38"}:
            raise ConfigurationError("ClinVar VCF assembly must be GRCh37 or GRCh38.")
        base = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{assembly}"
        resources = [RemoteFile(f"{base}/clinvar.vcf.gz", filename=f"clinvar_{assembly}.vcf.gz")]
        if include_index:
            resources.append(
                RemoteFile(f"{base}/clinvar.vcf.gz.tbi", filename=f"clinvar_{assembly}.vcf.gz.tbi")
            )
    elif kind == "xml":
        resources = [
            RemoteFile(
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz",
                filename="ClinVarVCVRelease_00-latest.xml.gz",
            )
        ]
    elif kind in {"tsv", "summary"}:
        resources = [
            RemoteFile(
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                filename="clinvar_variant_summary.txt.gz",
            )
        ]
    else:
        raise ConfigurationError("ClinVar format must be vcf, xml, or tsv.")
    return dataset(
        resources,
        output_dir,
        kind=f"{source_value}_{kind}",
        source="NCBI ClinVar" if source_value == "clinvar" else "NCBI dbSNP",
        metadata={
            "assembly": assembly if kind == "vcf" else None,
            "format": kind,
            "provider_checksum_files_included": include_checksum_files
            if source_value == "dbsnp"
            else False,
        },
        config=config,
        progress=progress,
    )


def _geo_range(accession: str) -> tuple[str, str, str]:
    match = _GEO_ACCESSION.fullmatch(accession.strip())
    if match is None:
        raise ConfigurationError("GEO accession must start with GSE, GDS, GPL, or GSM.")
    prefix = match.group(1).upper()
    digits = match.group(2)
    range_name = f"{prefix}{digits[:-3]}nnn" if len(digits) > 3 else f"{prefix}nnn"
    category = {"GSE": "series", "GDS": "datasets", "GPL": "platforms", "GSM": "samples"}[prefix]
    return prefix + digits, category, range_name


def expression(
    accession: str,
    output_dir: str | Path,
    *,
    format: str = "soft",
    platform: str | None = None,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
) -> DatasetDownloadResult:
    """Download a standard GEO SOFT, MINiML, or common Series RAW archive."""

    value, category, range_name = _geo_range(accession)
    kind = format.casefold()
    root = f"https://ftp.ncbi.nlm.nih.gov/geo/{category}/{range_name}/{value}"
    if kind == "soft":
        url = f"{root}/soft/{value}_family.soft.gz"
        filename = f"{value}_family.soft.gz"
    elif kind == "miniml":
        url = f"{root}/miniml/{value}_family.xml.tgz"
        filename = f"{value}_family.xml.tgz"
    elif kind in {"raw", "supplementary"}:
        if not value.startswith("GSE"):
            raise ConfigurationError("The common RAW archive convention applies to GSE accessions.")
        url = f"{root}/suppl/{value}_RAW.tar"
        filename = f"{value}_RAW.tar"
    elif kind in {"matrix", "series_matrix"}:
        if not value.startswith("GSE"):
            raise ConfigurationError("Series Matrix downloads require a GSE accession.")
        platform_value: str | None = None
        if platform is not None:
            platform_value, _, _ = _geo_range(platform)
            if not platform_value.startswith("GPL"):
                raise ConfigurationError("platform must be a GPL accession.")
        stem = value if platform_value is None else f"{value}-{platform_value}"
        filename = f"{stem}_series_matrix.txt.gz"
        url = f"{root}/matrix/{filename}"
    else:
        raise ConfigurationError("GEO format must be soft, miniml, raw, or matrix.")
    return dataset(
        (RemoteFile(url, filename=filename),),
        output_dir,
        kind=f"geo_{kind}",
        source="NCBI GEO",
        metadata={"accession": value, "format": kind, "platform": platform},
        config=config,
        progress=progress,
    )


__all__ = ["expression", "variants"]
