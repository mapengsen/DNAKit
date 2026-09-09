"""Resolve, validate and pin matching reference FASTA/GFF3 resources."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import tempfile
import zipfile
from collections.abc import Mapping
from pathlib import Path, PurePosixPath
from typing import Any, TextIO

from dnakit.core import Interval, Strand
from dnakit.core._json import to_json_compatible
from dnakit.download import DownloadConfig, genome_package
from dnakit.download.files import _write_manifest
from dnakit.exceptions import ConfigurationError, DownloadError, InputFormatError
from dnakit.io import read_gff3
from dnakit.references import resolve_genome_assembly, supported_genome_aliases
from dnakit.search import SearchConfig, assembly, taxonomy

from ._gene_mapping import ReferenceGene
from .gene_models import GeneAnnotationConfig

_ACCESSION = re.compile(r"GC[AF]_[0-9]+\.[0-9]+", re.IGNORECASE)
_GENE_TYPES = {
    "gene",
    "pseudogene",
    "C_gene_segment",
    "D_gene_segment",
    "J_gene_segment",
    "V_gene_segment",
}
_SCHEMA = "dnakit-gene-reference-v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1_048_576):
            digest.update(block)
    return digest.hexdigest()


def _text_file(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def fasta_lengths(path: Path, config: GeneAnnotationConfig) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    total = 0
    with _text_file(path) as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise InputFormatError("Reference FASTA has missing or duplicate sequence IDs.")
                name = fields[0]
                lengths[name] = 0
                if len(lengths) > 1_000_000:
                    raise InputFormatError("Reference FASTA has too many sequences.")
            else:
                if name is None or re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line) is None:
                    raise InputFormatError("Reference FASTA must contain ungapped IUPAC DNA.")
                lengths[name] += len(line)
                total += len(line)
                if total > config.max_reference_bases:
                    raise ConfigurationError("Reference exceeds max_reference_bases.")
    if not lengths or any(size == 0 for size in lengths.values()):
        raise InputFormatError("Reference FASTA contains no bases or an empty sequence.")
    return lengths


def read_genes(
    path: Path, lengths: dict[str, int], config: GeneAnnotationConfig
) -> tuple[tuple[ReferenceGene, ...], set[str], dict[str, str]]:
    """Stream gene-level rows through DNAKit's strict GFF3 codec.

    Transcript/CDS rows are not materialized. A contig is considered represented
    only if the GFF has a sequence-region directive, region, or gene on it.
    """

    selected = io.StringIO()
    represented: set[str] = set()
    circular: set[str] = set()
    metadata: dict[str, str] = {}
    count, gene_count = 0, 0
    with _text_file(path) as stream:
        for line_number, raw in enumerate(stream, 1):
            if len(raw) > 1_000_000:
                raise InputFormatError("GFF3 line exceeds one million characters.")
            line = raw.rstrip("\r\n")
            if line == "##FASTA":
                break
            if line.startswith("##gff-version"):
                if line_number != 1:
                    raise InputFormatError("GFF3 version directive must be the first line.")
                selected.write(line + "\n")
            elif line.startswith("##sequence-region "):
                fields = line.split()
                if len(fields) != 4:
                    raise InputFormatError("Invalid GFF3 sequence-region directive.")
                represented.add(fields[1])
                try:
                    if not 1 <= int(fields[2]) <= int(fields[3]) <= lengths[fields[1]]:
                        raise ValueError("out of bounds")
                except (ValueError, KeyError) as exc:
                    raise InputFormatError("FASTA and GFF3 sequence-region disagree.") from exc
            elif line.startswith("#!"):
                fields = line[2:].split(maxsplit=1)
                if len(fields) == 2:
                    metadata[fields[0]] = fields[1]
            elif line and not line.startswith("#"):
                count += 1
                if count > config.max_gff_records:
                    raise InputFormatError("GFF3 exceeds max_gff_records.")
                fields = line.split("\t")
                if len(fields) != 9:
                    raise InputFormatError(f"GFF3 line {line_number} needs nine columns.")
                if fields[2] == "region":
                    represented.add(fields[0])
                    if "Is_circular=true" in fields[8].split(";"):
                        circular.add(fields[0])
                if fields[2] in _GENE_TYPES:
                    gene_count += 1
                    if gene_count > 1_000_000:
                        raise InputFormatError("GFF3 exceeds one million gene records.")
                    selected.write(line + "\n")
    selected.seek(0)
    document = read_gff3(selected)
    if not document.entries:
        raise InputFormatError(
            "No gene-level features found in GFF3; a CDS/transcript-only file is insufficient.",
            code="GENE_ANNOTATION_UNAVAILABLE",
        )
    genes: dict[tuple[str, str], ReferenceGene] = {}
    for entry in document.entries:
        feature = entry.feature
        if entry.sequence_id not in lengths or not isinstance(feature.location, Interval):
            raise InputFormatError("GFF3 sequence IDs must exactly match the reference FASTA.")
        size = lengths[entry.sequence_id]
        part = feature.location
        parts: tuple[Interval, ...]
        if part.end > size:
            if entry.sequence_id not in circular or part.start >= size or len(part) > size:
                raise InputFormatError("GFF3 gene is outside the reference sequence bounds.")
            parts = (Interval(part.start, size), Interval(0, part.end - size))
        else:
            parts = (part,)
        represented.add(entry.sequence_id)
        qualifiers = feature.qualifiers

        def attribute(name: str, values: Mapping[str, object] = qualifiers) -> str | None:
            value = values.get(name)
            return value if isinstance(value, str) and value.strip() else None

        locus_tag = attribute("locus_tag")
        name = attribute("gene_name") or attribute("gene")
        if name is None and feature.label not in {None, locus_tag, feature.id}:
            name = feature.label
        dbxref = qualifiers.get("Dbxref", ())
        if isinstance(dbxref, str):
            dbxref = (dbxref,)
        gene_id = attribute("gene_id")
        if gene_id is None and isinstance(dbxref, tuple):
            gene_id = next(
                (
                    value
                    for value in dbxref
                    if isinstance(value, str) and value.startswith("GeneID:")
                ),
                None,
            )
        feature_id = feature.id or locus_tag or f"{entry.sequence_id}:{part.start + 1}-{part.end}"
        gene = ReferenceGene(
            entry.sequence_id,
            feature_id,
            gene_id or feature_id,
            name,
            locus_tag,
            attribute("gene_biotype") or attribute("biotype") or attribute("gene_type"),
            {Strand.FORWARD: "+", Strand.REVERSE: "-"}.get(feature.strand, "."),
            parts,
        )
        key = (entry.sequence_id, feature_id)
        if key in genes:
            previous = genes[key]
            if (previous.gene_id, previous.gene_name, previous.strand, previous.biotype) != (
                gene.gene_id,
                gene.gene_name,
                gene.strand,
                gene.biotype,
            ):
                raise InputFormatError("Inconsistent multi-part gene ID in GFF3.")
            gene = ReferenceGene(
                gene.sequence_id,
                feature_id,
                gene.gene_id,
                gene.gene_name,
                gene.locus_tag,
                gene.biotype,
                gene.strand,
                tuple(dict.fromkeys((*previous.parts, *parts))),
            )
        genes[key] = gene
    if represented - lengths.keys():
        raise InputFormatError("GFF3 contains contigs absent from FASTA; use matching files.")
    return tuple(genes.values()), represented, metadata


def _mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _normalize_report(value: Any) -> Any:
    """NCBI REST uses snake_case; downloaded protobuf JSONL uses camelCase."""

    if isinstance(value, dict):
        return {
            re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", key).lower(): _normalize_report(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_normalize_report(item) for item in value]
    return value


def _taxon(query: str, search_config: SearchConfig) -> Mapping[str, Any]:
    result = taxonomy(query, config=search_config)
    if len(result.records) != 1 or (result.total_count or 0) > 1:
        raise ConfigurationError("Species is ambiguous; use a scientific name or TaxID.")
    record = _mapping(result.records[0].get("taxonomy"))
    if not isinstance(record.get("tax_id"), int):
        raise ConfigurationError("NCBI did not resolve a valid TaxID.")
    if not _mapping(_mapping(record.get("classification")).get("species")):
        raise ConfigurationError("Provide a species or strain, not a higher taxonomic rank.")
    return record


def _validate_species(tax_id: int, organism: Mapping[str, Any], config: SearchConfig) -> None:
    actual = organism.get("tax_id")
    if isinstance(actual, str) and actual.isdigit():
        actual = int(actual)
    if isinstance(actual, bool) or not isinstance(actual, int) or actual < 1:
        raise ConfigurationError("Assembly metadata is missing a valid organism TaxID.")
    if actual == tax_id:
        return
    record = _taxon(str(actual), config)
    if record.get("tax_id") != tax_id and tax_id not in record.get("parents", ()):
        raise ConfigurationError(
            "The requested assembly does not belong to the supplied species/strain.",
            code="ASSEMBLY_SPECIES_MISMATCH",
            context={"requested_tax_id": tax_id, "assembly_tax_id": actual},
        )


def _resolve_report(
    tax_id: int, requested: str | None, search_config: SearchConfig
) -> tuple[str, Mapping[str, Any]]:
    if requested is None:
        # A species' reference may be registered to a descendant strain TaxID.
        result = assembly(
            str(tax_id),
            by="taxon",
            reference_only=True,
            assembly_source="refseq",
            current_only=True,
            page_size=1000,
            config=search_config,
        )
        candidates = sorted(
            {
                str(record["accession"])
                for record in result.records
                if _ACCESSION.fullmatch(str(record.get("accession", "")))
            }
        )
        if not candidates:
            raise ConfigurationError(
                "No current RefSeq reference assembly found; provide an accession or local files.",
                code="GENOME_ASSEMBLY_NOT_FOUND",
            )
        if len(candidates) != 1 or result.next_page_token or (result.total_count or 0) > 1:
            raise ConfigurationError(
                "Multiple reference assemblies found; specify a strain or an assembly accession.",
                code="AMBIGUOUS_GENOME_ASSEMBLY",
                context={"candidates": candidates},
            )
        accession = candidates[0]
    elif _ACCESSION.fullmatch(requested):
        accession = requested.upper()
    elif requested.casefold() in supported_genome_aliases():
        accession = resolve_genome_assembly(requested).accession
    else:
        # Assembly names are species-scoped. Do not select a newer accession in
        # response to an explicit old version or silently take a first hit.
        matches: set[str] = set()
        token = None
        seen_tokens: set[str] = set()
        for _ in range(100):
            result = assembly(
                str(tax_id),
                by="taxon",
                current_only=False,
                page_size=1000,
                page_token=token,
                config=search_config,
            )
            for record in result.records:
                name = _mapping(record.get("assembly_info")).get("assembly_name")
                if isinstance(name, str) and name.casefold() == requested.casefold():
                    matches.add(str(record["accession"]))
            token = result.next_page_token
            if token is None:
                break
            if token in seen_tokens:
                raise ConfigurationError("NCBI repeated an assembly pagination token.")
            seen_tokens.add(token)
        else:
            raise ConfigurationError("Assembly search is too broad; use a versioned accession.")
        # GenBank and RefSeq counterparts may share a name. The choice of
        # annotation source is material, so ask for the exact accession here.
        if len(matches) != 1:
            raise ConfigurationError(
                "Assembly name did not resolve uniquely; provide a versioned GCA_/GCF_ accession.",
                code="AMBIGUOUS_GENOME_ASSEMBLY" if matches else "GENOME_ASSEMBLY_NOT_FOUND",
                context={"candidates": sorted(matches)},
            )
        accession = next(iter(matches))
    result = assembly(accession, by="accession", current_only=False, config=search_config)
    exact = [record for record in result.records if record.get("accession") == accession]
    if len(exact) != 1:
        raise ConfigurationError("NCBI did not return the exact requested assembly version.")
    report = exact[0]
    _validate_species(tax_id, _mapping(report.get("organism")), search_config)
    if not _mapping(report.get("annotation_info")):
        raise ConfigurationError(
            "This assembly has no available NCBI annotation; provide matching local FASTA/GFF3.",
            code="GENE_ANNOTATION_UNAVAILABLE",
        )
    return accession, report


def _unpack_reference(
    package: Path, destination: Path, accession: str, max_bytes: int
) -> Mapping[str, Any]:
    """Extract only the selected pair, verifying the package's MD5 checksums."""

    with zipfile.ZipFile(package) as archive:
        names = archive.namelist()
        if len(names) > 1000 or len(set(names)) != len(names):
            raise DownloadError("Genome package has too many entries or duplicate names.")
        prefix = f"ncbi_dataset/data/{accession}/"
        fasta = [
            name for name in names if name.startswith(prefix) and name.endswith("_genomic.fna")
        ]
        gff = prefix + "genomic.gff"
        report_path = "ncbi_dataset/data/assembly_data_report.jsonl"
        if len(fasta) != 1 or gff not in names:
            raise DownloadError(
                "NCBI package must contain one reference FASTA and matching genomic.gff.",
                code="GENE_ANNOTATION_UNAVAILABLE",
            )
        for name in ("md5sum.txt", report_path):
            if name not in names or archive.getinfo(name).file_size > 20_000_000:
                raise DownloadError("Genome package metadata is missing or oversized.")
        checksums = {}
        for line in archive.read("md5sum.txt").decode("utf-8").splitlines():
            checksum, name = line.split(maxsplit=1)
            checksums[name.lstrip("*").removeprefix("./")] = checksum.lower()
        reports = [
            _normalize_report(json.loads(line))
            for line in archive.read(report_path).decode("utf-8").splitlines()
            if line.strip()
        ]
        exact = [report for report in reports if report.get("accession") == accession]
        if len(exact) != 1:
            raise DownloadError("Genome package report does not match the requested accession.")
        total = 0
        for source, filename in ((fasta[0], "reference.fna"), (gff, "annotation.gff3")):
            parts = PurePosixPath(source).parts
            if ".." in parts or source.startswith("/") or source not in checksums:
                raise DownloadError("Genome package member/checksum is invalid.")
            total += archive.getinfo(source).file_size
            if total > max_bytes:
                raise DownloadError("Uncompressed genome pair exceeds max_download_bytes.")
            digest = hashlib.md5(usedforsecurity=False)
            with archive.open(source) as stream, (destination / filename).open("xb") as output:
                size = 0
                while block := stream.read(1_048_576):
                    size += len(block)
                    if size > archive.getinfo(source).file_size:
                        raise DownloadError("Genome package member exceeds its declared size.")
                    digest.update(block)
                    output.write(block)
            if digest.hexdigest() != checksums[source]:
                raise DownloadError("Genome package failed its MD5 checksum verification.")
    return _mapping(exact[0])


def _json_object(path: Path) -> dict[str, Any]:
    if path.stat().st_size > 20_000_000:
        raise ConfigurationError("Reference cache manifest is oversized.")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ConfigurationError("Reference cache manifest must be an object.")
    return value


def _cached_reference(root: Path, accession: str) -> tuple[Path, Path, dict[str, Any]]:
    if _ACCESSION.fullmatch(accession) is None:
        raise ConfigurationError("Invalid accession in reference cache.")
    directory = root / "assemblies" / accession
    metadata = _json_object(directory / "reference.json")
    if metadata.get("schema") != _SCHEMA or metadata.get("assembly") != accession:
        raise ConfigurationError("Incompatible reference cache manifest.")
    fasta, gff = directory / "reference.fna", directory / "annotation.gff3"
    if sha256_file(fasta) != metadata.get("fasta_sha256") or sha256_file(gff) != metadata.get(
        "gff3_sha256"
    ):
        raise ConfigurationError(
            "Cached reference checksum mismatch; use a new cache directory.",
            code="REFERENCE_CACHE_INVALID",
        )
    return fasta, gff, metadata


def prepare_reference(
    species: str, requested: str | None, config: GeneAnnotationConfig
) -> tuple[Path, Path, dict[str, Any]]:
    """Persist the exact selected pair; subsequent identical requests are offline."""

    root = Path(config.cache_dir).expanduser().resolve()
    key = hashlib.sha256(
        json.dumps(
            [species.casefold(), None if requested is None else requested.casefold()]
        ).encode()
    ).hexdigest()
    pin = root / "requests" / f"{key}.json"
    if pin.is_file():
        value = _json_object(pin)
        cached_request = value.get("requested")
        if value.get("species") != species.casefold() or (
            None if cached_request is None else str(cached_request).casefold()
        ) != (None if requested is None else requested.casefold()):
            raise ConfigurationError("Reference request cache mismatch.")
        fasta, gff, metadata = _cached_reference(root, str(value.get("assembly")))
        return fasta, gff, {**metadata, "selection": "cached", "requested_assembly": requested}
    if not config.allow_download:
        raise ConfigurationError(
            "No pinned reference for this request and downloads are disabled.",
            code="REFERENCE_NOT_CACHED",
            hint="Provide local reference_fasta and annotation_gff, or enable downloads.",
        )
    search_config = SearchConfig(api_key=config.api_key)
    taxon = _taxon(species, search_config)
    tax_id = int(taxon["tax_id"])
    accession, report = _resolve_report(tax_id, requested, search_config)
    parent = root / "assemblies"
    parent.mkdir(parents=True, exist_ok=True)
    directory = parent / accession
    if not directory.exists():
        lock = parent / f"{accession}.lock"
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError as exc:
            raise ConfigurationError(
                "This assembly is being prepared; retry after it finishes."
            ) from exc
        os.close(descriptor)
        try:
            with tempfile.TemporaryDirectory(prefix=".gene-reference-", dir=parent) as raw:
                work = Path(raw)
                data = work / "pair"
                data.mkdir()
                package = genome_package(
                    accession,
                    work / "download",
                    include=("GENOME_FASTA", "GENOME_GFF", "SEQUENCE_REPORT"),
                    config=DownloadConfig(
                        api_key=config.api_key,
                        max_file_bytes=config.max_download_bytes,
                        max_total_bytes=config.max_download_bytes,
                    ),
                )
                report = _unpack_reference(
                    Path(package.files[0].path), data, accession, config.max_download_bytes
                )
                _validate_species(tax_id, _mapping(report.get("organism")), search_config)
                lengths = fasta_lengths(data / "reference.fna", config)
                _, _, headers = read_genes(data / "annotation.gff3", lengths, config)
                declared = _ACCESSION.search(headers.get("genome-build-accession", ""))
                if declared is not None and declared[0].upper() != accession:
                    raise DownloadError("GFF3 declares a different assembly than its data package.")
                metadata = {
                    "schema": _SCHEMA,
                    "assembly": accession,
                    "assembly_name": _mapping(report.get("assembly_info")).get("assembly_name"),
                    "organism": to_json_compatible(_mapping(report.get("organism"))),
                    "annotation": to_json_compatible(_mapping(report.get("annotation_info"))),
                    "annotation_headers": headers,
                    "source": "NCBI Datasets",
                    "fasta_sha256": sha256_file(data / "reference.fna"),
                    "gff3_sha256": sha256_file(data / "annotation.gff3"),
                    "package_sha256": package.files[0].sha256,
                    "species_validated": True,
                }
                _write_manifest(metadata, data / "reference.json", overwrite=False)
                os.rename(data, directory)
        finally:
            lock.unlink()
    fasta, gff, metadata = _cached_reference(root, accession)
    _validate_species(tax_id, _mapping(metadata.get("organism")), search_config)
    pin.parent.mkdir(parents=True, exist_ok=True)
    payload = {"species": species.casefold(), "requested": requested, "assembly": accession}
    try:
        _write_manifest(payload, pin, overwrite=False)
    except FileExistsError:
        if _json_object(pin) != payload:
            raise ConfigurationError("A concurrent request pinned a different assembly.") from None
    return (
        fasta,
        gff,
        {
            **metadata,
            "selection": "automatic" if requested is None else "explicit",
            "requested_assembly": requested,
        },
    )
