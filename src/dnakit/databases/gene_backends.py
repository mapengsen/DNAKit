"""Optional gene-level database clients: Ensembl, MyGene, MyVariant, BioMart and local GTF/FASTA.

The functions in this module are explicit adapters over optional packages:

- ``gget`` searches Ensembl genes, fetches gene metadata and sequences, and applies
  cDNA mutations to local sequence windows.
- ``pyensembl`` indexes a local Ensembl-style GTF and queries genes, loci and gene
  sequences from a local genome FASTA. It never downloads data.
- ``genomepy`` summarizes a local FASTA, extracts spliced and BED-defined sequences,
  and installs genome assemblies.
- ``biothings_client`` formats genomic HGVS identifiers from VCF-style records.
- ``mygene`` and ``myvariant`` query MyGene.info and MyVariant.info.
- ``pybiomart`` queries BioMart datasets; ``ncbi-datasets-pylib`` queries NCBI gene records.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint. Functions that
call live services are network-bound; their tests replace the client call and check the
parsing and validation.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import math
import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_GGET = "bio-gget"
_PYENSEMBL = "bio-pyensembl"
_GENOMEPY = "bio-genomepy"
_BIOTHINGS = "bio-biothings"
_MYGENE = "bio-mygene"
_MYVARIANT = "bio-myvariant"
_PYBIOMART = "bio-pybiomart"
_NCBI_DATASETS = "bio-ncbi-datasets"

_PYENSEMBL_REFERENCE = "dnakit_custom"
_BIOMART_HOST = "http://www.ensembl.org"
_MYGENE_FIELDS = "symbol,name,entrezgene,ensembl.gene,taxid"
_MYVARIANT_FIELDS = "dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,dbnsfp.genename"
_MAX_IDENTIFIERS = 200
_MAX_BATCH = 1000
_MAX_MUTATIONS = 1000
_MAX_SEQUENCE_LENGTH = 1_000_000
_MAX_REGIONS = 1_000_000
_MAX_RECORDS = 1_000_000

_ID_PATTERN = re.compile(r"[A-Za-z0-9._:/+_-]{1,64}")
_SPECIES_PATTERN = re.compile(r"[a-z][a-z0-9_]{1,127}")
_TAXA_PATTERN = re.compile(r"[A-Za-z0-9_,]{1,128}")
_FIELDS_PATTERN = re.compile(r"[A-Za-z0-9_.,*]{1,512}")
_SCOPES_PATTERN = re.compile(r"[A-Za-z0-9_.,]{1,128}")
_DNA_PATTERN = re.compile(r"[ACGTNRYSWKMBDHV]+")
_BASE_PATTERN = re.compile(r"[ACGTN]+")
_HGVS_CHROM_PATTERN = re.compile(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", re.IGNORECASE)
_HGVS_ID_PATTERN = re.compile(r"(?:chr)?[0-9A-Za-z]{1,6}:g\.[^\s]{1,256}")
_ASSEMBLY_PATTERN = re.compile(r"[A-Za-z0-9_.-]{1,128}")
_PROVIDER_PATTERN = re.compile(r"[a-z][a-z0-9_]{1,63}")
_BIOMART_NAME_PATTERN = re.compile(r"[A-Za-z0-9_]{1,128}")
_SYMBOL_PATTERN = re.compile(r"[A-Za-z0-9_.@-]{1,64}")
_TAXON_PATTERN = re.compile(r"[A-Za-z0-9 _.-]{1,64}")


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value or len(value) > limit:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return value.strip()


def _pattern(value: object, name: str, pattern: re.Pattern[str], *, limit: int = 256) -> str:
    text = _text(value, name, limit=limit)
    if pattern.fullmatch(text) is None:
        raise ConfigurationError(f"{name} has an unsupported format: {value!r}.")
    return text


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _sequence_of(values: object, name: str, *, limit: int) -> list[object]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} accepts at most {limit} items; got {len(values)}.")
    return list(values)


def _identifiers(values: object, name: str) -> list[str]:
    items = _sequence_of(values, name, limit=_MAX_IDENTIFIERS)
    return [_pattern(item, name, _ID_PATTERN) for item in items]


def _words(value: object) -> list[str]:
    if isinstance(value, str):
        return [_text(value, "searchwords", limit=128)]
    items = _sequence_of(value, "searchwords", limit=50)
    return [_text(item, "searchwords", limit=128) for item in items]


def _fields(value: object) -> str:
    if isinstance(value, (list, tuple)):
        value = ",".join(_text(item, "fields", limit=128) for item in value)
    return _pattern(value, "fields", _FIELDS_PATTERN, limit=512)


def _dna(value: object, name: str, *, max_length: int = _MAX_SEQUENCE_LENGTH) -> str:
    """IUPAC DNA string; used for sequence windows and genome intervals."""

    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a DNA string.")
    text = value.strip().upper()
    if not text or len(text) > max_length or _DNA_PATTERN.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} must be 1 to {max_length} IUPAC DNA symbols (ACGTNRYSWKMBDHV)."
        )
    return text


def _bases(value: object, name: str, *, max_length: int) -> str:
    """Uppercase A, C, G, T or N only; used for HGVS identifiers and mutation inputs."""

    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string of A, C, G, T or N.")
    text = value.strip().upper()
    if not text or len(text) > max_length or _BASE_PATTERN.fullmatch(text) is None:
        raise ConfigurationError(f"{name} must be 1 to {max_length} symbols from A, C, G, T, N.")
    return text


def _vcf_position(text: str, record: int) -> int:
    if not text.isdigit():
        raise ConfigurationError(f"VCF record {record} has a non-integer POS: {text!r}.")
    return _int(int(text), "POS", 1, 2_147_483_647)


def _plain(value: object) -> Any:
    """Convert numpy, pandas and float values to JSON-compatible Python objects."""

    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return str(value)


def _provenance(distribution: str, **notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"{distribution} {version}", **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": provenance}


def _upstream(package: str, operation: str, call: Callable[[], Any]) -> Any:
    """Run one upstream call; third-party clients raise many exception types."""

    try:
        return call()
    except Exception as exc:
        raise BackendExecutionError(f"{package} {operation} failed: {exc}") from exc


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _frame_rows(frame: Any) -> tuple[list[str], list[dict[str, Any]]]:
    """Return JSON-ready rows from a pandas frame; missing values become None."""

    columns = [str(column) for column in frame.columns]
    clean = frame.astype(object).where(frame.notna(), None)
    rows = [_plain(row) for row in clean.to_dict(orient="records")]
    return columns, rows


def _new_output(value: object, name: str) -> Path:
    path = Path(_text(value, name, limit=4096)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


# --- gget -----------------------------------------------------------------------------


def search_ensembl_genes(
    searchwords: str | Sequence[str],
    *,
    species: str = "homo_sapiens",
    release: int | None = None,
    andor: str = "or",
    max_results: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search Ensembl genes by free-text words through gget (network)."""

    gget = require("gget", _GGET)
    words = _words(searchwords)
    taxon = _pattern(species, "species", _SPECIES_PATTERN, limit=128)
    chosen_release = None if release is None else _int(release, "release", 1, 10_000)
    mode = _choice(andor, "andor", ("or", "and"))
    limit = _int(max_results, "max_results", 1, 1000)
    results = _upstream(
        "gget",
        "search",
        lambda: gget.search(
            words,
            species=taxon,
            release=chosen_release,
            andor=mode,
            limit=limit,
            json=True,
            verbose=False,
        ),
    )
    records = list(results or [])
    return _document(
        "search_ensembl_genes",
        _provenance("gget", source="Ensembl REST API via gget"),
        searchwords=words,
        species=taxon,
        release=chosen_release,
        count=len(records),
        results=_plain(records),
    )


def fetch_ensembl_gene_info(
    ensembl_ids: Sequence[str],
    *,
    external_sources: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch Ensembl gene metadata for each ID through gget (network)."""

    gget = require("gget", _GGET)
    ids = _identifiers(ensembl_ids, "ensembl_ids")
    external = _flag(external_sources, "external_sources")
    info = _upstream(
        "gget",
        "info",
        lambda: gget.info(ids, ncbi=external, uniprot=external, json=True, verbose=False),
    )
    mapping = info or {}
    records = {identifier: _plain(mapping.get(identifier)) for identifier in ids}
    return _document(
        "fetch_ensembl_gene_info",
        _provenance("gget", source="Ensembl REST API via gget", external_sources=external),
        requested=len(ids),
        found=sum(1 for value in records.values() if value is not None),
        records=records,
    )


def fetch_ensembl_sequences(
    ensembl_ids: Sequence[str],
    *,
    translate: bool = False,
    isoforms: bool = False,
    max_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch gene or transcript sequences from gget as parsed FASTA records (network)."""

    gget = require("gget", _GGET)
    ids = _identifiers(ensembl_ids, "ensembl_ids")
    protein = _flag(translate, "translate")
    isoform_flag = _flag(isoforms, "isoforms")
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    raw = _upstream(
        "gget",
        "seq",
        lambda: gget.seq(ids, translate=protein, isoforms=isoform_flag, verbose=False),
    )
    items = list(raw or [])
    if len(items) % 2:
        raise BackendExecutionError("gget seq returned a FASTA header without a sequence.")
    if len(items) // 2 > limit:
        raise ConfigurationError(f"gget returned more than max_records={limit} sequences.")
    records = []
    for header, sequence in zip(items[0::2], items[1::2], strict=True):
        if not isinstance(header, str) or not header.startswith(">"):
            raise BackendExecutionError("gget seq returned an unexpected FASTA layout.")
        identifier, _, description = header[1:].partition(" ")
        records.append({"id": identifier, "description": description, "sequence": sequence})
    return _document(
        "fetch_ensembl_sequences",
        _provenance("gget", source="Ensembl nucleotide or UniProt protein sequences via gget"),
        requested=len(ids),
        translate=protein,
        isoforms=isoform_flag,
        alphabet="protein" if protein else "nucleotide",
        count=len(records),
        records=records,
    )


def apply_cdna_mutations(
    sequences: Mapping[str, str],
    mutations: Mapping[str, str],
    *,
    flank: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply one c. notation mutation per sequence with gget.mutate (local, no network)."""

    gget = require("gget", _GGET)
    if not isinstance(sequences, Mapping) or not sequences:
        raise ConfigurationError("sequences must map sequence IDs to DNA strings.")
    if not isinstance(mutations, Mapping) or not mutations:
        raise ConfigurationError("mutations must map sequence IDs to c. notation strings.")
    if len(sequences) > _MAX_MUTATIONS or len(mutations) > _MAX_MUTATIONS:
        raise ConfigurationError(f"accepts at most {_MAX_MUTATIONS} sequences and mutations.")
    window = _int(flank, "flank", 1, 1000)
    clean: dict[str, str] = {}
    for key, value in sequences.items():
        clean[_text(key, "sequence ID", limit=64)] = _bases(value, "sequence", max_length=100_000)
    checked: list[tuple[str, str]] = []
    for key, mutation in mutations.items():
        seq_id = _text(key, "mutation sequence ID", limit=64)
        if seq_id not in clean:
            raise ConfigurationError(f"mutation refers to an unknown sequence ID: {seq_id!r}.")
        notation = _text(mutation, "mutation", limit=64)
        if not notation.startswith("c.") or any(char.isspace() for char in notation):
            raise ConfigurationError("mutation must use c. notation such as 'c.2C>T'.")
        checked.append((seq_id, notation))
    rows = []
    for seq_id, notation in checked:
        output = _upstream(
            "gget",
            "mutate",
            lambda seq=clean[seq_id], mut=notation: gget.mutate(
                [seq], [mut], k=window, merge_identical=False, verbose=False
            ),
        )
        mutant = output[0] if output else None
        rows.append(
            {
                "seq_id": seq_id,
                "mutation": notation,
                "status": "applied" if mutant is not None else "dropped",
                "mutant_window": mutant,
            }
        )
    return _document(
        "apply_cdna_mutations",
        _provenance("gget", source="gget.mutate", flank=window),
        requested=len(rows),
        applied=sum(1 for row in rows if row["status"] == "applied"),
        mutations=rows,
        window_convention="flank bases each side of the 1-based c. position, clipped at ends",
    )


# --- pyensembl ------------------------------------------------------------------------


def _gtf_genome(gtf: object, cache_dir: object, genome_fasta: object = None) -> Any:
    pyensembl = require("pyensembl", _PYENSEMBL)
    path = input_path(gtf)  # type: ignore[arg-type]
    cache = None
    if cache_dir is not None:
        cache = str(Path(_text(str(cache_dir), "cache_dir", limit=4096)).expanduser().resolve())
    fasta = None if genome_fasta is None else str(input_path(genome_fasta))  # type: ignore[arg-type]
    # The annotation name carries a content digest, so an edited GTF gets a fresh index.
    annotation = f"gtf_{_file_digest(path)[:16]}"
    genome = _upstream(
        "pyensembl",
        "Genome",
        lambda: pyensembl.Genome(
            reference_name=_PYENSEMBL_REFERENCE,
            annotation_name=annotation,
            gtf_path_or_url=str(path),
            genome_fasta_path_or_url=fasta,
            cache_directory_path=cache,
        ),
    )
    _upstream("pyensembl", "index", genome.index)
    return genome


def _gene_record(genome: Any, gene: Any) -> dict[str, Any]:
    return {
        "gene_id": gene.gene_id,
        "gene_name": gene.gene_name,
        "contig": gene.contig,
        "start": int(gene.start),
        "end": int(gene.end),
        "strand": gene.strand,
        "biotype": gene.biotype,
        "transcript_ids": list(genome.transcript_ids_of_gene_id(gene.gene_id)),
    }


def index_gene_models(
    gtf: str | Path, *, cache_dir: str | Path | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Index a local Ensembl-style GTF with pyensembl and count its features."""

    genome = _gtf_genome(gtf, cache_dir)
    return _document(
        "index_gene_models",
        _provenance("pyensembl", coordinate_system="1-based inclusive", reference="custom GTF"),
        source=str(input_path(gtf)),
        cache_directory=str(genome.download_cache.cache_directory_path),
        contigs=len(genome.contigs()),
        genes=len(genome.gene_ids()),
        transcripts=len(genome.transcript_ids()),
        exons=len(genome.exon_ids()),
    )


def query_genes_by_name(
    gtf: str | Path,
    gene_names: Sequence[str],
    *,
    cache_dir: str | Path | None = None,
    max_results: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up genes by exact gene name in a local GTF through pyensembl."""

    names = [
        _text(item, "gene_names", limit=128)
        for item in _sequence_of(gene_names, "gene_names", limit=_MAX_IDENTIFIERS)
    ]
    limit = _int(max_results, "max_results", 1, _MAX_RECORDS)
    genome = _gtf_genome(gtf, cache_dir)
    results: dict[str, list[dict[str, Any]]] = {}
    total = 0
    for name in names:
        try:
            gene_ids = genome.gene_ids_of_gene_name(name)
        except ValueError:  # pyensembl reports a name with no matches as ValueError
            gene_ids = []
        total += len(gene_ids)
        if total > limit:
            raise ConfigurationError(f"query matches more than max_results={limit} genes.")
        results[name] = [_gene_record(genome, genome.gene_by_id(gid)) for gid in gene_ids]
    return _document(
        "query_genes_by_name",
        _provenance("pyensembl", coordinate_system="1-based inclusive"),
        source=str(input_path(gtf)),
        matches=total,
        results=results,
    )


def query_genes_at_locus(
    gtf: str | Path,
    contig: str,
    start: int,
    end: int | None = None,
    *,
    strand: str | None = None,
    cache_dir: str | Path | None = None,
    max_results: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """List genes overlapping a 1-based inclusive locus in a local GTF through pyensembl."""

    chrom = _text(contig, "contig", limit=64)
    first = _int(start, "start", 1, 2_147_483_647)
    last = first if end is None else _int(end, "end", first, 2_147_483_647)
    chosen_strand = None if strand is None else _choice(strand, "strand", ("+", "-"))
    limit = _int(max_results, "max_results", 1, _MAX_RECORDS)
    genome = _gtf_genome(gtf, cache_dir)
    if chrom not in genome.contigs():
        raise ConfigurationError(f"contig {chrom!r} is not in the annotation.")
    genes = _upstream(
        "pyensembl",
        "genes_at_locus",
        lambda: genome.genes_at_locus(chrom, first, last, strand=chosen_strand),
    )
    if len(genes) > limit:
        raise ConfigurationError(f"locus overlaps more than max_results={limit} genes.")
    return _document(
        "query_genes_at_locus",
        _provenance("pyensembl", coordinate_system="1-based inclusive"),
        contig=chrom,
        start=first,
        end=last,
        strand=chosen_strand,
        count=len(genes),
        genes=[_gene_record(genome, gene) for gene in genes],
    )


def fetch_gene_sequence(
    gtf: str | Path,
    genome_fasta: str | Path,
    gene_id: str,
    *,
    cache_dir: str | Path | None = None,
    max_length: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a gene's genomic locus from a local genome FASTA, oriented by gene strand."""

    identifier = _pattern(gene_id, "gene_id", _ID_PATTERN, limit=64)
    limit = _int(max_length, "max_length", 1, 10_000_000)
    genome = _gtf_genome(gtf, cache_dir, genome_fasta)
    try:
        gene = genome.gene_by_id(identifier)
    except ValueError as exc:
        raise ConfigurationError(f"gene_id {identifier!r} is not in the annotation.") from exc
    length = int(gene.end) - int(gene.start) + 1
    if length > limit:
        raise ConfigurationError(f"gene locus length {length} exceeds max_length={limit}.")
    sequence = _upstream(
        "pyensembl",
        "sequence",
        lambda: genome.sequence(gene.contig, gene.start, gene.end, strand=gene.strand),
    )
    return _document(
        "fetch_gene_sequence",
        _provenance("pyensembl", coordinate_system="1-based inclusive", orientation="gene strand"),
        gene_id=gene.gene_id,
        gene_name=gene.gene_name,
        contig=gene.contig,
        start=int(gene.start),
        end=int(gene.end),
        strand=gene.strand,
        length=length,
        sequence=str(sequence),
    )


# --- genomepy -------------------------------------------------------------------------


def read_genome_summary(
    fasta: str | Path, *, max_contigs: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """Summarize a local genome FASTA with genomepy; writes index and sidecar files beside it."""

    genomepy = require("genomepy", _GENOMEPY)
    path = input_path(fasta)
    limit = _int(max_contigs, "max_contigs", 1, 10_000_000)
    genome = _upstream("genomepy", "Genome", lambda: genomepy.Genome(str(path)))
    sizes = dict(genome.sizes)
    if len(sizes) > limit:
        raise ConfigurationError(f"FASTA has more than max_contigs={limit} records.")
    gaps = dict(genome.gaps or {})
    contigs = [
        {"name": name, "length": int(length), "n_bases": int(gaps.get(name, 0))}
        for name, length in sizes.items()
    ]
    tax_id = None if genome.tax_id in (None, "na") else str(genome.tax_id)
    accession = None
    if genome.assembly_accession not in (None, "na"):
        accession = str(genome.assembly_accession)
    return _document(
        "read_genome_summary",
        _provenance(
            "genomepy", coordinate_system="1-based inclusive", side_effects="sidecar files"
        ),
        source=str(path),
        name=genome.name,
        tax_id=tax_id,
        assembly_accession=accession,
        contig_count=len(contigs),
        total_length=sum(item["length"] for item in contigs),
        contigs=contigs,
    )


def fetch_spliced_sequence(
    fasta: str | Path,
    contig: str,
    intervals: Sequence[Sequence[int]],
    *,
    reverse_complement: bool = False,
    max_length: int = _MAX_SEQUENCE_LENGTH,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join 1-based inclusive intervals from one contig, as for spliced exons."""

    genomepy = require("genomepy", _GENOMEPY)
    path = input_path(fasta)
    chrom = _text(contig, "contig", limit=128)
    flag = _flag(reverse_complement, "reverse_complement")
    limit = _int(max_length, "max_length", 1, 10_000_000)
    pairs: list[tuple[int, int]] = []
    for item in _sequence_of(intervals, "intervals", limit=10_000):
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)) or len(item) != 2:
            raise ConfigurationError("each interval must be a [start, end] pair.")
        start = _int(item[0], "interval start", 1, 2_147_483_647)
        end = _int(item[1], "interval end", start, 2_147_483_647)
        pairs.append((start, end))
    if sum(end - start + 1 for start, end in pairs) > limit:
        raise ConfigurationError(f"spliced sequence exceeds max_length={limit}.")
    genome = _upstream("genomepy", "Genome", lambda: genomepy.Genome(str(path)))
    if chrom not in genome.sizes:
        raise ConfigurationError(f"contig {chrom!r} is not in the FASTA.")
    for start, end in pairs:
        if end > genome.sizes[chrom]:
            raise ConfigurationError(f"interval {start}-{end} exceeds the contig length.")
    spliced = _upstream(
        "genomepy",
        "get_spliced_seq",
        lambda: genome.get_spliced_seq(chrom, pairs, rc=flag),
    )
    return _document(
        "fetch_spliced_sequence",
        _provenance("genomepy", coordinate_system="1-based inclusive intervals"),
        source=str(path),
        contig=chrom,
        intervals=[[start, end] for start, end in pairs],
        reverse_complement=flag,
        length=len(spliced.seq),
        sequence=str(spliced.seq),
    )


def write_bed_sequences(
    fasta: str | Path,
    bed: str | Path,
    output: str | Path,
    *,
    stranded: bool = False,
    extend_up: int = 0,
    extend_down: int = 0,
    max_regions: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write FASTA records for BED regions with genomepy; never overwrites the output."""

    genomepy = require("genomepy", _GENOMEPY)
    fasta_path = input_path(fasta)
    bed_path = input_path(bed)
    out_path = _new_output(str(output), "output")
    flag = _flag(stranded, "stranded")
    up = _int(extend_up, "extend_up", 0, 1_000_000)
    down = _int(extend_down, "extend_down", 0, 1_000_000)
    limit = _int(max_regions, "max_regions", 1, _MAX_REGIONS)
    regions = 0
    with bed_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                regions += 1
                if regions > limit:
                    raise ConfigurationError(f"BED has more than max_regions={limit} regions.")
    genome = _upstream("genomepy", "Genome", lambda: genomepy.Genome(str(fasta_path)))
    _upstream(
        "genomepy",
        "track2fasta",
        lambda: genome.track2fasta(
            str(bed_path),
            fastafile=str(out_path),
            stranded=flag,
            extend_up=up,
            extend_down=down,
        ),
    )
    with out_path.open("r", encoding="utf-8") as handle:
        records = sum(1 for line in handle if line.startswith(">"))
    return _document(
        "write_bed_sequences",
        _provenance("genomepy", coordinate_system="BED 0-based half-open"),
        source=str(fasta_path),
        bed=str(bed_path),
        output=str(out_path),
        regions=regions,
        records=records,
        skipped=regions - records,
        stranded=flag,
        extend_up=up,
        extend_down=down,
    )


def install_genome_assembly(
    name: str,
    *,
    provider: str | None = None,
    genomes_dir: str | Path | None = None,
    annotation: bool = False,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and install a genome assembly with genomepy (network, large files)."""

    genomepy = require("genomepy", _GENOMEPY)
    assembly = _pattern(name, "name", _ASSEMBLY_PATTERN, limit=128)
    chosen_provider = None
    if provider is not None:
        chosen_provider = _pattern(provider, "provider", _PROVIDER_PATTERN, limit=64)
    directory = None
    if genomes_dir is not None:
        directory = str(Path(_text(str(genomes_dir), "genomes_dir", limit=4096)).expanduser())
    flag = _flag(annotation, "annotation")
    workers = _int(threads, "threads", 1, 64)
    genome = _upstream(
        "genomepy",
        "install_genome",
        lambda: genomepy.install_genome(
            assembly,
            provider=chosen_provider,
            genomes_dir=directory,
            annotation=flag,
            threads=workers,
        ),
    )
    return _document(
        "install_genome_assembly",
        _provenance("genomepy", source=chosen_provider or "genomepy default provider"),
        name=assembly,
        provider=chosen_provider,
        fasta=str(genome.genome_file),
        annotation_gtf=genome.annotation_gtf_file,
    )


# --- biothings_client -----------------------------------------------------------------


def format_hgvs_id(
    chrom: str, pos: int, ref: str, alt: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Format a genomic HGVS identifier from one VCF-style record (local, no network)."""

    biothings = require("biothings_client", _BIOTHINGS)
    chromosome = _text(chrom, "chrom", limit=16)
    if _HGVS_CHROM_PATTERN.fullmatch(chromosome) is None:
        raise ConfigurationError(f"chrom must be a human chromosome name; got {chrom!r}.")
    position = _int(pos, "pos", 1, 2_147_483_647)
    reference = _bases(ref, "ref", max_length=1000)
    alternate = _bases(alt, "alt", max_length=1000)
    if reference == alternate:
        raise ConfigurationError("ref and alt must differ.")
    client = biothings.get_client("variant")
    try:
        hgvs = client.format_hgvs(chromosome, position, reference, alternate)
    except ValueError as exc:
        raise ConfigurationError(f"cannot format an HGVS id: {exc}") from exc
    return _document(
        "format_hgvs_id",
        _provenance(
            "biothings_client", coordinate_system="VCF 1-based POS", assembly="not encoded"
        ),
        hgvs_id=hgvs,
        chrom=chromosome,
        pos=position,
        ref=reference,
        alt=alternate,
    )


def hgvs_ids_from_vcf(
    path: str | Path, *, max_records: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """Format genomic HGVS ids for every sequence ALT in a VCF through biothings_client."""

    biothings = require("biothings_client", _BIOTHINGS)
    vcf = input_path(path)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    client = biothings.get_client("variant")
    identifiers: list[str] = []
    skipped = {"unsupported_contig": 0, "non_sequence_allele": 0, "same_as_ref": 0}
    records = 0
    with vcf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            records += 1
            if records > limit:
                raise ConfigurationError(f"VCF has more than max_records={limit} records.")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ConfigurationError(f"VCF record {records} has fewer than five columns.")
            position = _vcf_position(fields[1], records)
            chrom = fields[0]
            if _HGVS_CHROM_PATTERN.fullmatch(chrom) is None:
                skipped["unsupported_contig"] += 1
                continue
            ref = fields[3].upper()
            for alt_text in fields[4].split(","):
                alt = alt_text.upper()
                if _BASE_PATTERN.fullmatch(alt) is None or _BASE_PATTERN.fullmatch(ref) is None:
                    skipped["non_sequence_allele"] += 1
                    continue
                if alt == ref:
                    skipped["same_as_ref"] += 1
                    continue
                identifiers.append(client.format_hgvs(chrom, position, ref, alt))
    return _document(
        "hgvs_ids_from_vcf",
        _provenance(
            "biothings_client", coordinate_system="VCF 1-based POS", assembly="not encoded"
        ),
        source=str(vcf),
        records=records,
        hgvs_ids=identifiers,
        skipped=skipped,
    )


# --- mygene ---------------------------------------------------------------------------


def query_mygene_genes(
    query: str,
    *,
    species: str = "human",
    fields: str | Sequence[str] = _MYGENE_FIELDS,
    size: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search MyGene.info gene records with a query string (network)."""

    mygene = require("mygene", _MYGENE)
    text = _text(query, "query", limit=512)
    taxon = _pattern(species, "species", _TAXA_PATTERN, limit=128)
    field_text = _fields(fields)
    count = _int(size, "size", 1, 1000)
    client = mygene.MyGeneInfo()
    result = _upstream(
        "mygene",
        "query",
        lambda: client.query(text, species=taxon, fields=field_text, size=count),
    )
    hits = list((result or {}).get("hits", []))
    return _document(
        "query_mygene_genes",
        _provenance("mygene", source="MyGene.info"),
        query=text,
        species=taxon,
        total=int((result or {}).get("total", len(hits))),
        hits=_plain(hits),
    )


def fetch_mygene_gene(
    gene_id: str,
    *,
    species: str = "human",
    fields: str | Sequence[str] = _MYGENE_FIELDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch one MyGene.info gene record by Entrez or Ensembl ID (network)."""

    mygene = require("mygene", _MYGENE)
    identifier = _pattern(gene_id, "gene_id", _ID_PATTERN, limit=64)
    taxon = _pattern(species, "species", _TAXA_PATTERN, limit=128)
    field_text = _fields(fields)
    client = mygene.MyGeneInfo()
    gene = _upstream(
        "mygene",
        "getgene",
        lambda: client.getgene(identifier, fields=field_text, species=taxon),
    )
    return _document(
        "fetch_mygene_gene",
        _provenance("mygene", source="MyGene.info"),
        gene_id=identifier,
        found=gene is not None,
        gene=_plain(gene),
    )


def map_gene_identifiers(
    identifiers: Sequence[str],
    *,
    scopes: str = "entrezgene",
    species: str = "human",
    fields: str | Sequence[str] = "symbol,name,taxid",
    verbose: bool = True,
) -> dict[str, Any]:
    """Map gene identifiers of one scope to MyGene.info records in batch (network)."""

    mygene = require("mygene", _MYGENE)
    terms = [
        _text(item, "identifiers", limit=64)
        for item in _sequence_of(identifiers, "identifiers", limit=_MAX_BATCH)
    ]
    scope = _pattern(scopes, "scopes", _SCOPES_PATTERN, limit=128)
    taxon = _pattern(species, "species", _TAXA_PATTERN, limit=128)
    field_text = _fields(fields)
    client = mygene.MyGeneInfo()
    rows = _upstream(
        "mygene",
        "querymany",
        lambda: client.querymany(terms, scopes=scope, fields=field_text, species=taxon),
    )
    return _document(
        "map_gene_identifiers",
        _provenance("mygene", source="MyGene.info batch query"),
        scopes=scope,
        requested=len(terms),
        results=_plain(list(rows or [])),
    )


# --- myvariant ------------------------------------------------------------------------


def fetch_variant_annotations(
    hgvs_ids: Sequence[str],
    *,
    assembly: str = "hg38",
    fields: str | Sequence[str] = _MYVARIANT_FIELDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch MyVariant.info annotations for genomic HGVS ids (network)."""

    myvariant = require("myvariant", _MYVARIANT)
    ids = [
        _pattern(item, "hgvs_ids", _HGVS_ID_PATTERN, limit=300)
        for item in _sequence_of(hgvs_ids, "hgvs_ids", limit=_MAX_BATCH)
    ]
    build = _choice(assembly, "assembly", ("hg19", "hg38"))
    field_text = _fields(fields)
    client = myvariant.MyVariantInfo()
    rows = _upstream(
        "myvariant",
        "getvariants",
        lambda: client.getvariants(ids, fields=field_text, assembly=build),
    )
    return _document(
        "fetch_variant_annotations",
        _provenance("myvariant", source="MyVariant.info", assembly=build),
        requested=len(ids),
        results=_plain(list(rows or [])),
    )


def query_variant_annotations(
    query: str,
    *,
    assembly: str = "hg38",
    fields: str | Sequence[str] = _MYVARIANT_FIELDS,
    size: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search MyVariant.info variant records with a query string (network)."""

    myvariant = require("myvariant", _MYVARIANT)
    text = _text(query, "query", limit=512)
    build = _choice(assembly, "assembly", ("hg19", "hg38"))
    field_text = _fields(fields)
    count = _int(size, "size", 1, 1000)
    client = myvariant.MyVariantInfo()
    result = _upstream(
        "myvariant",
        "query",
        lambda: client.query(text, fields=field_text, assembly=build, size=count),
    )
    hits = list((result or {}).get("hits", []))
    return _document(
        "query_variant_annotations",
        _provenance("myvariant", source="MyVariant.info", assembly=build),
        query=text,
        total=int((result or {}).get("total", len(hits))),
        hits=_plain(hits),
    )


# --- pybiomart ------------------------------------------------------------------------


def query_biomart_dataset(
    dataset: str,
    attributes: Sequence[str],
    *,
    filters: Mapping[str, object] | None = None,
    host: str = _BIOMART_HOST,
    only_unique: bool = True,
    use_attr_names: bool = False,
    max_rows: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query a BioMart dataset with explicit attributes and filters (network)."""

    pybiomart = require("pybiomart", _PYBIOMART)
    name = _pattern(dataset, "dataset", _BIOMART_NAME_PATTERN, limit=128)
    names = [
        _pattern(item, "attributes", _BIOMART_NAME_PATTERN, limit=128)
        for item in _sequence_of(attributes, "attributes", limit=200)
    ]
    filter_map: dict[str, object] | None = None
    if filters is not None:
        if not isinstance(filters, Mapping) or not filters:
            raise ConfigurationError("filters must be a non-empty mapping or None.")
        filter_map = {
            _pattern(key, "filter name", _BIOMART_NAME_PATTERN, limit=128): value
            for key, value in filters.items()
        }
    unique = _flag(only_unique, "only_unique")
    attribute_names = _flag(use_attr_names, "use_attr_names")
    limit = _int(max_rows, "max_rows", 1, 1_000_000)
    host_text = _text(host, "host", limit=256)
    handle = pybiomart.Dataset(name=name, host=host_text)
    frame = _upstream(
        "pybiomart",
        "query",
        lambda: handle.query(
            attributes=names,
            filters=filter_map,
            only_unique=unique,
            use_attr_names=attribute_names,
        ),
    )
    if len(frame) > limit:
        raise ConfigurationError(f"BioMart returned more than max_rows={limit} rows.")
    columns, rows = _frame_rows(frame)
    return _document(
        "query_biomart_dataset",
        _provenance("pybiomart", source=host_text),
        dataset=name,
        host=host_text,
        columns=columns,
        row_count=len(rows),
        rows=rows,
    )


def list_biomart_attributes(
    dataset: str,
    *,
    host: str = _BIOMART_HOST,
    max_rows: int = 5000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the attributes a BioMart dataset exposes (network)."""

    pybiomart = require("pybiomart", _PYBIOMART)
    name = _pattern(dataset, "dataset", _BIOMART_NAME_PATTERN, limit=128)
    limit = _int(max_rows, "max_rows", 1, 1_000_000)
    host_text = _text(host, "host", limit=256)
    handle = pybiomart.Dataset(name=name, host=host_text)
    frame = _upstream("pybiomart", "list_attributes", handle.list_attributes)
    if len(frame) > limit:
        raise ConfigurationError(f"dataset lists more than max_rows={limit} attributes.")
    columns, rows = _frame_rows(frame)
    return _document(
        "list_biomart_attributes",
        _provenance("pybiomart", source=host_text),
        dataset=name,
        columns=columns,
        count=len(rows),
        attributes=rows,
    )


def list_biomart_datasets(
    mart: str = "ENSEMBL_MART_ENSEMBL",
    *,
    host: str = _BIOMART_HOST,
    max_rows: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the datasets in one BioMart mart (network)."""

    pybiomart = require("pybiomart", _PYBIOMART)
    mart_name = _pattern(mart, "mart", _BIOMART_NAME_PATTERN, limit=128)
    limit = _int(max_rows, "max_rows", 1, 1_000_000)
    host_text = _text(host, "host", limit=256)
    server = pybiomart.Server(host=host_text)
    marts = _upstream("pybiomart", "marts", lambda: server.marts)
    if mart_name not in marts:
        raise ConfigurationError(f"mart {mart_name!r} is not listed by {host_text}.")
    frame = _upstream("pybiomart", "list_datasets", marts[mart_name].list_datasets)
    if len(frame) > limit:
        raise ConfigurationError(f"mart lists more than max_rows={limit} datasets.")
    columns, rows = _frame_rows(frame)
    return _document(
        "list_biomart_datasets",
        _provenance("pybiomart", source=host_text),
        mart=mart_name,
        columns=columns,
        count=len(rows),
        datasets=rows,
    )


# --- ncbi-datasets-pylib --------------------------------------------------------------


def _ncbi_gene_api() -> Any:
    ncbi_datasets = require("ncbi.datasets", _NCBI_DATASETS)
    return _upstream("ncbi-datasets-pylib", "GeneApi", ncbi_datasets.GeneApi)


def _ncbi_payload(response: Any) -> dict[str, Any]:
    """Convert an openapi model response to plain JSON-compatible data."""

    if response is None:
        return {}
    to_dict = getattr(response, "to_dict", None)
    return _plain(to_dict() if callable(to_dict) else response)


def fetch_ncbi_gene_by_id(
    gene_ids: Sequence[int], *, max_genes: int = 100, verbose: bool = True
) -> dict[str, Any]:
    """Fetch NCBI gene metadata for Gene IDs through ncbi-datasets-pylib (network)."""

    ids = [
        _int(item, "gene_ids", 1, 2**63 - 1)
        for item in _sequence_of(gene_ids, "gene_ids", limit=_MAX_IDENTIFIERS)
    ]
    limit = _int(max_genes, "max_genes", 1, 1000)
    api = _ncbi_gene_api()
    response = _upstream(
        "ncbi-datasets-pylib",
        "gene_metadata_by_id",
        lambda: api.gene_metadata_by_id(ids, limit=str(limit)),
    )
    payload = _ncbi_payload(response)
    genes = list(payload.get("genes") or [])
    return _document(
        "fetch_ncbi_gene_by_id",
        _provenance("ncbi-datasets-pylib", source="NCBI Datasets gene API"),
        requested=len(ids),
        total_count=payload.get("total_count"),
        genes=genes,
    )


def fetch_ncbi_gene_by_symbol(
    symbols: Sequence[str],
    taxon: str,
    *,
    max_genes: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch NCBI gene metadata for gene symbols within one taxon (network)."""

    names = [
        _pattern(item, "symbols", _SYMBOL_PATTERN, limit=64)
        for item in _sequence_of(symbols, "symbols", limit=_MAX_IDENTIFIERS)
    ]
    organism = _pattern(taxon, "taxon", _TAXON_PATTERN, limit=64)
    limit = _int(max_genes, "max_genes", 1, 1000)
    api = _ncbi_gene_api()
    response = _upstream(
        "ncbi-datasets-pylib",
        "gene_metadata_by_tax_and_symbol",
        lambda: api.gene_metadata_by_tax_and_symbol(names, organism, limit=str(limit)),
    )
    payload = _ncbi_payload(response)
    genes = list(payload.get("genes") or [])
    return _document(
        "fetch_ncbi_gene_by_symbol",
        _provenance("ncbi-datasets-pylib", source="NCBI Datasets gene API"),
        symbols=names,
        taxon=organism,
        total_count=payload.get("total_count"),
        genes=genes,
    )


def fetch_ncbi_gene_orthologs(
    gene_id: int,
    *,
    taxon_filter: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch the NCBI ortholog set for one Gene ID, optionally filtered by taxon (network)."""

    identifier = _int(gene_id, "gene_id", 1, 2**63 - 1)
    taxa = None
    if taxon_filter is not None:
        taxa = [
            _pattern(item, "taxon_filter", _TAXON_PATTERN, limit=64)
            for item in _sequence_of(taxon_filter, "taxon_filter", limit=100)
        ]
    extra: dict[str, Any] = {} if taxa is None else {"taxon_filter": taxa}
    api = _ncbi_gene_api()
    response = _upstream(
        "ncbi-datasets-pylib",
        "gene_orthologs_by_id",
        lambda: api.gene_orthologs_by_id(identifier, **extra),
    )
    payload = _ncbi_payload(response)
    return _document(
        "fetch_ncbi_gene_orthologs",
        _provenance("ncbi-datasets-pylib", source="NCBI Datasets gene API"),
        gene_id=identifier,
        taxon_filter=taxa,
        ortholog_set_id=payload.get("ortholog_set_id"),
        method=payload.get("method"),
        genes=list(payload.get("genes") or []),
    )


__all__ = [
    "apply_cdna_mutations",
    "fetch_ensembl_gene_info",
    "fetch_ensembl_sequences",
    "fetch_gene_sequence",
    "fetch_mygene_gene",
    "fetch_ncbi_gene_by_id",
    "fetch_ncbi_gene_by_symbol",
    "fetch_ncbi_gene_orthologs",
    "fetch_spliced_sequence",
    "fetch_variant_annotations",
    "format_hgvs_id",
    "hgvs_ids_from_vcf",
    "index_gene_models",
    "install_genome_assembly",
    "list_biomart_attributes",
    "list_biomart_datasets",
    "map_gene_identifiers",
    "query_biomart_dataset",
    "query_genes_at_locus",
    "query_genes_by_name",
    "query_mygene_genes",
    "query_variant_annotations",
    "read_genome_summary",
    "search_ensembl_genes",
    "write_bed_sequences",
]

_register_public_verbose(globals(), tuple(__all__))
