"""varcode operations on a local Ensembl-style reference (a GTF plus a genome FASTA), without downloads.

varcode predicts variant effects with pyensembl transcript models. ``build_local_reference`` indexes a
user-supplied GTF and genome FASTA as a pyensembl reference: it derives the spliced cDNA and the translated
CDS (protein) FASTA from the GTF and the genome, then indexes them with pyensembl. Nothing is downloaded.
The other functions take the same ``gtf`` and ``genome_fasta`` paths, reuse the indexed reference, and call
varcode with it. Variants are passed as mappings with ``contig``, ``start`` (1-based), ``ref`` and ``alt``.

Functions that take a transcript take its ``transcript_id`` string. Results are JSON-serialisable dictionaries
or lists of them; varcode effect objects are summarised as class name and description, the transcript and the
gene where they apply.

``varcode_infer_genome`` and ``varcode_infer_genome_for_reference_name`` resolve Ensembl releases through
pyensembl; the names are looked up in locally installed releases and may need the network when a release is
not installed. The other functions never use the network.

Packages are imported only when a function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.metadata
import itertools
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_VARCODE_EXTRA = "bio-varcode"
_PYENSEMBL_EXTRA = "bio-pyensembl"
_READS_EXTRA = "bio-reads"
_REFERENCE_NAME = "dnakit_local"
_BUILDER_VERSION = "dnakit-local-reference-1"
_CODON_TCAG = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = {
    "".join(bases): aa
    for bases, aa in zip(itertools.product("TCAG", repeat=3), _CODON_TCAG, strict=True)
}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ATTRIBUTE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s+"([^"]*)"')
_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:+|-]{0,127}")
_BASES = re.compile(r"[ACGTNacgtn]*")
_STRANDS = ("+", "-")
_ANNOTATORS = ("fast", "protein_diff", "transcript_model", "realized")
_FAILURE_CLASSES = frozenset({"Failure", "Unresolved"})
_SV_TYPES = ("DEL", "DUP", "INV", "INS", "CNV")
_MAX_TEXT = 10_000
_MAX_POSITION = 2_000_000_000
_MAX_VARIANTS = 1_000
_MAX_RANDOM = 10_000
_MAX_CANDIDATES = 1_000
_DEFAULT_CACHE = Path.home() / ".cache" / "dnakit" / "varcode_reference"
_MEMO: dict[str, Any] = {}


# --- validation and conversion ------------------------------------------------------------


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > _MAX_TEXT:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _identifier(value: object, name: str) -> str:
    text = _text(value, name)
    if _IDENTIFIER.fullmatch(text) is None:
        raise ConfigurationError(f"{name} must be a plain identifier (letters, digits and ._:+|-).")
    return text


def _int(value: object, name: str, *, minimum: int = 0, maximum: int = _MAX_POSITION) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _bases(value: object, name: str, *, allow_empty: bool) -> str:
    if not isinstance(value, str) or _BASES.fullmatch(value) is None:
        raise ConfigurationError(f"{name} must contain only A, C, G, T or N.")
    if not value and not allow_empty:
        raise ConfigurationError(f"{name} must not be empty.")
    return value.upper()


def _mapping(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping with contig, start, ref and alt.")
    return value


def _variants(values: object, name: str) -> list[Mapping[str, Any]]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of variant mappings.")
    if len(values) > _MAX_VARIANTS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_VARIANTS} variants.")
    return [_mapping(item, f"{name}[{index}]") for index, item in enumerate(values)]


def _path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if not target.is_file():
        raise ConfigurationError(f"{name} does not exist: {target}")
    return target


def _cache(value: object) -> Path:
    if value is None:
        return _DEFAULT_CACHE
    return Path(_text(str(value), "cache_dir")).expanduser().resolve()


def _plain(value: Any) -> Any:
    """Make results JSON-serialisable: scalars, mappings and sequences only."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_plain(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _plain(getattr(value, field.name)) for field in dataclasses.fields(value)
        }
    return str(value)


def _varcode_version() -> str:
    try:
        return importlib.metadata.version("varcode")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _document(operation: str, **fields: Any) -> dict[str, Any]:
    return {
        "operation": operation,
        "status": "completed",
        "backend": "varcode",
        "package": f"varcode {_varcode_version()}",
        "coordinates": "1-based inclusive",
        **{key: _plain(value) for key, value in fields.items()},
    }


def _vc(*parts: str) -> Any:
    """Import one varcode submodule (for example ``("effects",)``) on demand."""
    return require(".".join(["varcode", *parts]), _VARCODE_EXTRA)


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run one varcode/pyensembl call; reject upstream value errors as configuration errors."""
    try:
        return function(*args, **kwargs)
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc


# --- local reference -----------------------------------------------------------------------


def _revcomp(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _translate(coding: str) -> str:
    protein = "".join(
        _CODONS.get(coding[index : index + 3].upper(), "X")
        for index in range(0, len(coding) - 2, 3)
    )
    return protein.rstrip("*")


def _parse_gtf(path: Path) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """Group exon and CDS lines by transcript and return the normalised GTF lines."""
    transcripts: dict[str, dict[str, Any]] = {}
    lines: list[str] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip() or raw.startswith("#"):
            lines.append(raw)
            continue
        fields = raw.split("\t")
        if len(fields) != 9:
            raise ConfigurationError(f"GTF line {number} must have 9 tab-separated fields.")
        contig, feature, strand = fields[0], fields[2], fields[6]
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError as exc:
            raise ConfigurationError(f"GTF line {number} has non-integer coordinates.") from exc
        if start < 1 or end < start or strand not in _STRANDS:
            raise ConfigurationError(f"GTF line {number} has invalid coordinates or strand.")
        attributes = dict(_ATTRIBUTE.findall(fields[8]))
        transcript_id = attributes.get("transcript_id")
        if feature in ("exon", "CDS") and transcript_id is None:
            raise ConfigurationError(f"GTF line {number} ({feature}) lacks transcript_id.")
        if feature == "CDS" and "protein_id" not in attributes:
            attributes["protein_id"] = transcript_id or ""
            fields[8] = fields[8].rstrip() + f' protein_id "{transcript_id}";'
        if feature == "exon" and "exon_id" not in attributes:
            # pyensembl indexes exons by exon_id; GTFs without one get a stable id from transcript and order.
            ordinal = sum(1 for key in transcripts.get(transcript_id, {}).get("exons", ())) + 1
            fields[8] = fields[8].rstrip() + f' exon_id "{transcript_id}.exon{ordinal}";'
        if transcript_id is not None and feature in ("exon", "CDS", "transcript"):
            entry = transcripts.setdefault(
                transcript_id,
                {"contig": contig, "strand": strand, "exons": [], "cds": [], "protein_id": None},
            )
            if entry["contig"] != contig or entry["strand"] != strand:
                raise ConfigurationError(
                    f"Transcript {transcript_id} has segments on different contigs or strands."
                )
            if feature == "exon":
                entry["exons"].append((start, end))
            elif feature == "CDS":
                entry["cds"].append((start, end))
                entry["protein_id"] = attributes["protein_id"]
        lines.append("\t".join(fields))
    return transcripts, lines


def _derived_fastas(
    transcripts: Mapping[str, Mapping[str, Any]],
    genome_fasta: Path,
    cdna_path: Path,
    protein_path: Path,
) -> tuple[int, int]:
    pysam = require("pysam", _READS_EXTRA)
    fasta = pysam.FastaFile(str(genome_fasta))
    cdna_count = protein_count = 0
    try:
        with (
            cdna_path.open("w", encoding="utf-8") as cdna,
            protein_path.open("w", encoding="utf-8") as proteins,
        ):
            for transcript_id, entry in transcripts.items():
                if entry["contig"] not in fasta.references:
                    raise ConfigurationError(f"Contig {entry['contig']} is not in {genome_fasta}.")
                if not entry["exons"]:
                    continue
                exons = sorted(entry["exons"])
                sequence = "".join(
                    fasta.fetch(entry["contig"], start - 1, end) for start, end in exons
                )
                if entry["strand"] == "-":
                    sequence = _revcomp(sequence)
                cdna.write(f">{transcript_id}\n{sequence}\n")
                cdna_count += 1
                if entry["cds"]:
                    coding = "".join(
                        fasta.fetch(entry["contig"], start - 1, end)
                        for start, end in sorted(entry["cds"])
                    )
                    if entry["strand"] == "-":
                        coding = _revcomp(coding)
                    proteins.write(f">{entry['protein_id']}\n{_translate(coding)}\n")
                    protein_count += 1
    finally:
        fasta.close()
    return cdna_count, protein_count


def _digest(gtf: Path, genome_fasta: Path) -> str:
    stat = genome_fasta.stat()
    hasher = hashlib.sha256()
    hasher.update(_BUILDER_VERSION.encode())
    hasher.update(gtf.read_bytes())
    hasher.update(f"{genome_fasta}|{stat.st_size}|{stat.st_mtime_ns}".encode())
    return hasher.hexdigest()[:20]


def _genome(gtf: object, genome_fasta: object, cache_dir: object) -> tuple[Any, dict[str, Any]]:
    """Index (or reuse) the local pyensembl reference for a GTF and genome FASTA."""
    gtf_path = _path(gtf, "gtf")
    fasta_path = _path(genome_fasta, "genome_fasta")
    digest = _digest(gtf_path, fasta_path)
    root = _cache(cache_dir) / digest
    cached = _MEMO.get(digest)
    if cached is not None:
        genome, indexed_in = cached
        return genome, {"digest": digest, "cache_directory": str(indexed_in), "reused": True}
    pyensembl = require("pyensembl", _PYENSEMBL_EXTRA)
    root.mkdir(parents=True, exist_ok=True)
    transcripts, lines = _parse_gtf(gtf_path)
    normalised = root / "annotation.gtf"
    normalised.write_text("\n".join(lines) + "\n", encoding="utf-8")
    cdna_path, protein_path = root / "cdna.fa", root / "protein.fa"
    cdna_count, protein_count = _derived_fastas(transcripts, fasta_path, cdna_path, protein_path)
    if cdna_count == 0:
        raise ConfigurationError("The GTF has no transcript with exon lines.")
    genome = pyensembl.Genome(
        reference_name=_REFERENCE_NAME,
        annotation_name=f"local_{digest}",
        gtf_path_or_url=str(normalised),
        transcript_fasta_paths_or_urls=[str(cdna_path)],
        protein_fasta_paths_or_urls=[str(protein_path)] if protein_count else None,
        genome_fasta_path_or_url=str(fasta_path),
        cache_directory_path=str(root / "pyensembl"),
    )
    try:
        genome.index()
    except (ValueError, KeyError, OSError) as exc:
        raise BackendExecutionError(
            f"pyensembl could not index the local reference: {exc}"
        ) from exc
    _MEMO[digest] = (genome, root)
    info = {
        "digest": digest,
        "cache_directory": str(root),
        "reused": False,
        "transcripts_with_cdna": cdna_count,
    }
    return genome, info


def _variant(genome: Any, record: object, name: str = "variant") -> Any:
    item = _mapping(record, name)
    missing = sorted({"contig", "start", "ref", "alt"} - set(item))
    if missing:
        raise ConfigurationError(f"{name} is missing fields: {', '.join(missing)}.")
    contig = _identifier(item["contig"], f"{name}.contig")
    start = _int(item["start"], f"{name}.start", minimum=1)
    ref = _bases(item["ref"], f"{name}.ref", allow_empty=True)
    alt = _bases(item["alt"], f"{name}.alt", allow_empty=True)
    if not ref and not alt:
        raise ConfigurationError(f"{name} needs a non-empty ref or alt.")
    varcode = require("varcode", _VARCODE_EXTRA)
    return _call(varcode.Variant, contig, start, ref, alt, genome=genome)


def _variant_record(variant: Any) -> dict[str, Any]:
    return {
        "contig": str(variant.contig),
        "start": int(variant.start),
        "ref": str(variant.ref),
        "alt": str(variant.alt),
    }


def _annotator(name: object) -> Any:
    if name is None:
        return None
    if name not in _ANNOTATORS:
        raise ConfigurationError(f"annotator must be one of: {', '.join(_ANNOTATORS)}.")
    registry = require("varcode", _VARCODE_EXTRA).annotators.registry
    return registry.get_annotator(str(name))


def _transcript(genome: Any, transcript_id: object) -> Any:
    identifier = _identifier(transcript_id, "transcript_id")
    try:
        return genome.transcript_by_id(identifier)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Transcript {identifier} is not in the local reference.") from exc


def _describe(effect: Any) -> str:
    description = getattr(effect, "short_description", None)
    if callable(description):
        description = description()
    return str(description) if description is not None else str(effect)


def _effect_record(effect: Any) -> dict[str, Any]:
    record: dict[str, Any] = {"class": type(effect).__name__, "description": _describe(effect)}
    variant = getattr(effect, "variant", None)
    if variant is not None and hasattr(variant, "contig"):
        record["variant"] = _variant_record(variant)
    transcript = getattr(effect, "transcript", None)
    if transcript is not None:
        record["transcript_id"] = str(getattr(transcript, "transcript_id", ""))
        record["gene_name"] = getattr(transcript, "gene_name", None)
    return record


def _candidate_record(candidate: Any) -> dict[str, Any]:
    if hasattr(candidate, "effect"):
        return {
            "source": str(getattr(candidate, "source", "")),
            "evidence": _plain(getattr(candidate, "evidence", {}) or {}),
            "effect": _effect_record(candidate.effect),
        }
    return _effect_record(candidate)


def _mutant_record(mutant: Any) -> dict[str, Any]:
    reference = getattr(mutant, "reference_transcript", None)
    return {
        "transcript_id": str(getattr(reference, "transcript_id", "")),
        "annotator": str(getattr(mutant, "annotator_name", "")),
        "cdna_sequence": getattr(mutant, "cdna_sequence", None),
        "mutant_protein_sequence": getattr(mutant, "mutant_protein_sequence", None),
        "edits": [_plain(edit) for edit in getattr(mutant, "edits", ())],
    }


def _layout(genome: Any, transcript: Any, flank: object) -> Any:
    """A GenomicLayout around one transcript, with ``flank`` bases on each side.

    The sequence provider reads the local reference's genome FASTA (1-based inclusive ``reference_range``).
    """
    window = _int(flank, "flank", minimum=0, maximum=100_000)
    varcode = require("varcode", _VARCODE_EXTRA)

    def provider(contig: str, start: int, end: int) -> str:
        return str(varcode.reference_range(genome, contig, start, end))

    return _call(
        _vc("genomic_layout").GenomicLayout.from_transcript,
        transcript,
        flank=window,
        sequence_provider=provider,
    )


def _sequence_of(values: object, name: str) -> list[Any]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if len(values) > _MAX_CANDIDATES:
        raise ConfigurationError(f"{name} may hold at most {_MAX_CANDIDATES} items.")
    return list(values)


# --- public: local reference and reference sequence ------------------------------------------


def build_local_reference(
    gtf: str,
    genome_fasta: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index a local GTF and genome FASTA as a pyensembl reference for varcode, without downloading.

    The spliced cDNA and the translated CDS are derived from the GTF exon and CDS lines and the genome
    FASTA. Exons must be on contigs of the FASTA. The index is cached under ``cache_dir`` (default
    ``~/.cache/dnakit/varcode_reference``) and keyed by the content of the GTF and the FASTA identity. An index
    already loaded in this process is reused, and its original cache folder is reported.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    return _document(
        "build_local_reference",
        reference_name=genome.reference_name,
        annotation_name=genome.annotation_name,
        cache_directory=info["cache_directory"],
        reused=info["reused"],
        contigs=len(genome.contigs()),
        genes=len(genome.gene_ids()),
        transcripts=len(genome.transcript_ids()),
        exons=len(genome.exon_ids()),
    )


def varcode_reference_base(
    gtf: str,
    genome_fasta: str,
    contig: str,
    position: int,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the forward-strand reference base at a 1-based position (varcode.reference_base)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    name = _identifier(contig, "contig")
    base = _call(varcode.reference_base, genome, name, _int(position, "position", minimum=1))
    return _document(
        "reference_base",
        contig=name,
        position=position,
        base=base,
        cache_directory=info["cache_directory"],
    )


def varcode_reference_range(
    gtf: str,
    genome_fasta: str,
    contig: str,
    start: int,
    end: int,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the forward-strand reference sequence over an inclusive 1-based interval (varcode.reference_range)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    first = _int(start, "start", minimum=1)
    last = _int(end, "end", minimum=first)
    name = _identifier(contig, "contig")
    sequence = _call(varcode.reference_range, genome, name, first, last)
    return _document(
        "reference_range",
        contig=name,
        start=first,
        end=last,
        sequence=sequence,
        cache_directory=info["cache_directory"],
    )


def varcode_left_align_indels(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shift indels to their leftmost equivalent position using the local reference (varcode.transforms)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    collection = varcode.VariantCollection(
        [_variant(genome, item) for item in _variants(variants, "variants")]
    )
    aligned = _call(_vc("transforms").left_align_indels, collection)
    rows = [_variant_record(variant) for variant in aligned]
    return _document(
        "left_align_indels", variants=rows, count=len(rows), cache_directory=info["cache_directory"]
    )


def varcode_random_variants(
    gtf: str,
    genome_fasta: str,
    count: int,
    *,
    deletions: bool = True,
    insertions: bool = True,
    random_seed: int | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw random variants that overlap complete coding transcripts of the local reference (varcode.util)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    total = _int(count, "count", minimum=1, maximum=_MAX_RANDOM)
    if not isinstance(deletions, bool) or not isinstance(insertions, bool):
        raise ConfigurationError("deletions and insertions must be booleans.")
    seed = (
        None
        if random_seed is None
        else _int(random_seed, "random_seed", minimum=0, maximum=2**32 - 1)
    )
    drawn = _call(
        _vc("util").random_variants,
        total,
        deletions=deletions,
        insertions=insertions,
        random_seed=seed,
        ensembl=genome,
    )
    rows = [_variant_record(variant) for variant in drawn]
    return _document(
        "random_variants",
        count=len(rows),
        random_seed=seed,
        variants=rows,
        cache_directory=info["cache_directory"],
    )


def varcode_infer_genome(reference: str | int, *, verbose: bool = True) -> dict[str, Any]:
    """Resolve an Ensembl release by integer or reference name through pyensembl (varcode.infer_genome).

    An integer is an Ensembl release number. A string is a reference name such as ``GRCh38`` or ``GRCh38:93``;
    UCSC names such as ``hg19`` are converted to the Ensembl equivalent. The release object is created lazily.
    """
    require("varcode", _VARCODE_EXTRA)
    if isinstance(reference, bool) or not isinstance(reference, (str, int)):
        raise ConfigurationError("reference must be an Ensembl release number or a reference name.")
    if isinstance(reference, str):
        _text(reference, "reference")
    genome, converted = _call(_vc("reference").infer_genome, reference)
    return _document(
        "infer_genome",
        requested=reference,
        reference_name=str(getattr(genome, "reference_name", "")),
        release=getattr(genome, "release", None),
        species=str(getattr(getattr(genome, "species", None), "latin_name", "")),
        converted_from_ucsc=bool(converted),
    )


def varcode_infer_genome_for_reference_name(
    reference_name: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Resolve a reference name, optionally with a release, to an Ensembl release (varcode.reference)."""
    require("varcode", _VARCODE_EXTRA)
    name = _text(reference_name, "reference_name")
    genome, converted = _call(_vc("reference").infer_genome_for_reference_name, name)
    return _document(
        "infer_genome_for_reference_name",
        requested=name,
        reference_name=str(getattr(genome, "reference_name", "")),
        release=getattr(genome, "release", None),
        species=str(getattr(getattr(genome, "species", None), "latin_name", "")),
        converted_from_ucsc=bool(converted),
    )


# --- public: effect prediction on the local reference ----------------------------------------


def varcode_predict_variant_effects(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
    annotator: str | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the effect of one variant on every transcript it overlaps (varcode.effects.predict_variant_effects)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    effects = _call(
        _vc("effects").predict_variant_effects,
        item,
        raise_on_error=bool(raise_on_error),
        annotator=_annotator(annotator),
    )
    rows = [_effect_record(effect) for effect in effects]
    return _document(
        "predict_variant_effects",
        variant=_variant_record(item),
        effects=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def varcode_predict_variant_effect_on_transcript(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    *,
    annotator: str | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the effect (such as a frameshift) of one variant on one transcript."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    effect = _call(
        _vc("effects").predict_variant_effect_on_transcript,
        item,
        _transcript(genome, transcript_id),
        annotator=_annotator(annotator),
    )
    return _document(
        "predict_variant_effect_on_transcript",
        variant=_variant_record(item),
        transcript_id=_identifier(transcript_id, "transcript_id"),
        effect=_effect_record(effect),
        cache_directory=info["cache_directory"],
    )


def varcode_predict_variant_effect_on_transcript_or_failure(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict a variant's effect on one transcript; report a failure instead of raising (varcode)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    effect = _vc("effects").predict_variant_effect_on_transcript_or_failure(
        item, _transcript(genome, transcript_id)
    )
    failed = isinstance(effect, Exception) or type(effect).__name__ in _FAILURE_CLASSES
    return _document(
        "predict_variant_effect_on_transcript_or_failure",
        variant=_variant_record(item),
        transcript_id=_identifier(transcript_id, "transcript_id"),
        status="failed" if failed else "completed",
        effect=None if failed else _effect_record(effect),
        failure=str(effect) if failed else None,
        cache_directory=info["cache_directory"],
    )


def varcode_predict_realized_effect(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    transcript_id: str,
    *,
    germline_variants: Sequence[Mapping[str, Any]] = (),
    max_hypotheses: int = 64,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the realised top effect of somatic variants on one transcript, with alternatives (varcode).

    Somatic variants are applied together on the transcript; ``germline_variants`` are phased against them
    when varcode can. ``candidates`` lists the alternative outcomes that were kept (``max_hypotheses``).
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    somatic = [_variant(genome, item, "variants[]") for item in _variants(variants, "variants")]
    germline = [
        _variant(genome, item, "germline_variants[]")
        for item in _variants(germline_variants, "germline_variants")
    ]
    limit = _int(max_hypotheses, "max_hypotheses", minimum=1, maximum=10_000)
    effect = _call(
        _vc("realized_effects").predict_realized_effect,
        somatic,
        _transcript(genome, transcript_id),
        germline_variants=germline,
        max_hypotheses=limit,
    )
    candidates = [
        _candidate_record(item)
        for item in list(getattr(effect, "candidates", ()))[:_MAX_CANDIDATES]
    ]
    return _document(
        "predict_realized_effect",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        top_effect=_effect_record(effect),
        candidates=candidates,
        cache_directory=info["cache_directory"],
    )


def varcode_predict_germline_aware_effect(
    gtf: str,
    genome_fasta: str,
    somatic_variant: Mapping[str, Any],
    transcript_id: str,
    germline_variants: Sequence[Mapping[str, Any]],
    *,
    annotator: str | None = None,
    max_hypotheses: int | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict a somatic variant's effect against the germline-applied transcript (varcode.germline).

    The germline variants form a complete ``GermlineContext``; phase is inferred only when varcode can resolve
    it, and otherwise the effect is reported with the hypotheses varcode enumerated.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    somatic = _variant(genome, somatic_variant, "somatic_variant")
    germline = [
        _variant(genome, item, "germline_variants[]")
        for item in _variants(germline_variants, "germline_variants")
    ]
    context = _vc("germline").GermlineContext.from_variants(germline)
    selected = _annotator(annotator) or _vc("annotators", "registry").get_default_annotator()
    limit = (
        None
        if max_hypotheses is None
        else _int(max_hypotheses, "max_hypotheses", minimum=1, maximum=10_000)
    )
    effect = _call(
        _vc("germline").predict_germline_aware_effect,
        somatic,
        _transcript(genome, transcript_id),
        context,
        selected,
        max_hypotheses=limit,
    )
    return _document(
        "predict_germline_aware_effect",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        somatic_variant=_variant_record(somatic),
        germline_count=len(germline),
        effect=_effect_record(effect),
        cache_directory=info["cache_directory"],
    )


def varcode_default_germline_window(
    gtf: str,
    genome_fasta: str,
    somatic_variant: Mapping[str, Any],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the genomic window searched for germline variants that can phase with a somatic variant."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    somatic = _variant(genome, somatic_variant, "somatic_variant")
    contig, start, end = _call(
        _vc("germline").default_germline_window, somatic, _transcript(genome, transcript_id)
    )
    return _document(
        "default_germline_window",
        contig=str(contig),
        start=int(start),
        end=int(end),
        cache_directory=info["cache_directory"],
    )


def varcode_apply_germline_to_transcript(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    germline_variants: Sequence[Mapping[str, Any]],
    *,
    somatic_variant: Mapping[str, Any] | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply germline variants to a transcript and return the resulting mutant transcript (varcode.germline)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    germline = [
        _variant(genome, item, "germline_variants[]")
        for item in _variants(germline_variants, "germline_variants")
    ]
    context = _vc("germline").GermlineContext.from_variants(germline)
    somatic = (
        None if somatic_variant is None else _variant(genome, somatic_variant, "somatic_variant")
    )
    mutant = _call(
        _vc("germline").apply_germline_to_transcript,
        _transcript(genome, transcript_id),
        context,
        somatic,
    )
    return _document(
        "apply_germline_to_transcript",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        mutant_transcript=_mutant_record(mutant),
        cache_directory=info["cache_directory"],
    )


def varcode_apply_variant_to_transcript(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply one variant to a transcript and return the mutant cDNA and protein (varcode.mutant_transcript)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    mutant = _call(
        _vc("mutant_transcript").apply_variant_to_transcript,
        item,
        _transcript(genome, transcript_id),
    )
    return _document(
        "apply_variant_to_transcript",
        variant=_variant_record(item),
        transcript_id=_identifier(transcript_id, "transcript_id"),
        mutant_transcript=_mutant_record(mutant),
        cache_directory=info["cache_directory"],
    )


def varcode_apply_variants_to_transcript(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply several variants together to one transcript and return the resulting mutant (varcode.mutant_transcript)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    items = [_variant(genome, item, "variants[]") for item in _variants(variants, "variants")]
    mutant = _call(
        _vc("mutant_transcript").apply_variants_to_transcript,
        items,
        _transcript(genome, transcript_id),
    )
    return _document(
        "apply_variants_to_transcript",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        variants=[_variant_record(item) for item in items],
        mutant_transcript=_mutant_record(mutant),
        cache_directory=info["cache_directory"],
    )


def varcode_classify_from_protein_diff(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    ref_protein: str,
    mut_protein: str,
    *,
    length_delta: int = 0,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify a variant's protein change by diffing a reference and a mutant protein (varcode.annotators)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    ref = _text(ref_protein, "ref_protein")
    mut = _text(mut_protein, "mut_protein")
    delta = _int(length_delta, "length_delta", minimum=-_MAX_POSITION, maximum=_MAX_POSITION)
    effect = _call(
        _vc("annotators", "protein_diff").classify_from_protein_diff,
        item,
        _transcript(genome, transcript_id),
        ref,
        mut,
        length_delta=delta,
    )
    return _document(
        "classify_from_protein_diff",
        variant=_variant_record(item),
        transcript_id=_identifier(transcript_id, "transcript_id"),
        effect=_effect_record(effect),
        cache_directory=info["cache_directory"],
    )


def varcode_codon_table_for_transcript(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the codon table for a transcript from its contig (varcode.annotators.protein_diff)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    transcript = _transcript(genome, transcript_id)
    table = _call(_vc("annotators", "protein_diff").codon_table_for_transcript, transcript)
    return _document(
        "codon_table_for_transcript",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        contig=str(transcript.contig),
        codon_table={
            "name": str(getattr(table, "name", "")),
            "ncbi_id": getattr(table, "ncbi_id", None),
            "start_codons": sorted(getattr(table, "start_codons", ()) or ()),
            "stop_codons": sorted(getattr(table, "stop_codons", ()) or ()),
        },
        cache_directory=info["cache_directory"],
    )


def varcode_reference_selenocysteine(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the cDNA offsets of annotated selenocysteine codons in a transcript (varcode.effects.selenocysteine)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    require("varcode", _VARCODE_EXTRA)
    offsets = _call(
        _vc("effects", "selenocysteine").reference_selenocysteine,
        _transcript(genome, transcript_id),
    )
    return _document(
        "reference_selenocysteine",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        cdna_offsets=[int(offset) for offset in offsets],
        cache_directory=info["cache_directory"],
    )


def varcode_candidates_from_effects(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    *,
    source: str = "varcode",
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict a variant's effects and wrap each as an effect candidate (varcode.candidates_from_effects).

    Effect objects cannot be passed in as JSON, so this predicts the variant's effects on the local reference
    first and then wraps them with varcode's ``candidates_from_effects`` under ``source``.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    text = _text(source, "source")
    effects = tuple(_vc("effects").predict_variant_effects(item))
    candidates = _call(varcode.candidates_from_effects, effects, source=text)
    rows = [_candidate_record(candidate) for candidate in list(candidates)[:_MAX_CANDIDATES]]
    return _document(
        "candidates_from_effects",
        variant=_variant_record(item),
        source=text,
        candidates=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def varcode_make_rna_outcome(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    *,
    source: str = "rna",
    transcript_model_id: str | None = None,
    read_count: int | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Wrap a variant's effect on one transcript as an RNA-sourced candidate (varcode.make_rna_outcome).

    The caller supplies the RNA source label and optional model id and read count. No RNA reads are read here;
    the effect itself is predicted on the local reference.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    effect = _vc("effects").predict_variant_effect_on_transcript(
        item, _transcript(genome, transcript_id)
    )
    count = (
        None
        if read_count is None
        else _int(read_count, "read_count", minimum=0, maximum=_MAX_POSITION)
    )
    candidate = _call(
        varcode.make_rna_outcome,
        effect,
        source=_text(source, "source"),
        transcript_model_id=None
        if transcript_model_id is None
        else _identifier(transcript_model_id, "transcript_model_id"),
        read_count=count,
    )
    return _document(
        "make_rna_outcome",
        variant=_variant_record(item),
        transcript_id=_identifier(transcript_id, "transcript_id"),
        candidate=_candidate_record(candidate),
        cache_directory=info["cache_directory"],
    )


def varcode_predict_structural_variant_effect(
    gtf: str,
    genome_fasta: str,
    structural_variant: Mapping[str, Any],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify a structural variant (deletion, duplication, inversion, insertion or CNV) on one transcript.

    ``structural_variant`` has ``contig``, ``start``, ``end`` (1-based inclusive) and ``sv_type`` (DEL, DUP, INV,
    INS or CNV). It is built as a varcode ``StructuralVariant`` on the local reference.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    item = _mapping(structural_variant, "structural_variant")
    missing = sorted({"contig", "start", "end", "sv_type"} - set(item))
    if missing:
        raise ConfigurationError(f"structural_variant is missing fields: {', '.join(missing)}.")
    contig = _identifier(item["contig"], "structural_variant.contig")
    start = _int(item["start"], "structural_variant.start", minimum=1)
    end = _int(item["end"], "structural_variant.end", minimum=start)
    sv_type = _text(item["sv_type"], "structural_variant.sv_type").upper()
    if sv_type not in _SV_TYPES:
        raise ConfigurationError(
            f"structural_variant.sv_type must be one of: {', '.join(_SV_TYPES)}."
        )
    sv = _call(
        _vc("structural_variant").StructuralVariant, contig, start, sv_type, end=end, genome=genome
    )
    effect = _call(
        _vc("effects", "structural").predict_structural_variant_effect,
        sv,
        _transcript(genome, transcript_id),
    )
    return _document(
        "predict_structural_variant_effect",
        structural_variant={"contig": contig, "start": start, "end": end, "sv_type": sv_type},
        transcript_id=_identifier(transcript_id, "transcript_id"),
        effect=_effect_record(effect),
        cache_directory=info["cache_directory"],
    )


def varcode_enumerate_from_structural_variant(
    gtf: str,
    genome_fasta: str,
    structural_variant: Mapping[str, Any],
    *,
    flank: int = 500,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate cryptic-exon candidates in the flanks of a structural variant breakpoint (varcode.cryptic_exons).

    ``structural_variant`` is a mapping as in ``varcode_predict_structural_variant_effect``. Candidates are scored
    with varcode's default splice-site scorer and reported with their genomic coordinates.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    item = _mapping(structural_variant, "structural_variant")
    missing = sorted({"contig", "start", "end", "sv_type"} - set(item))
    if missing:
        raise ConfigurationError(f"structural_variant is missing fields: {', '.join(missing)}.")
    contig = _identifier(item["contig"], "structural_variant.contig")
    start = _int(item["start"], "structural_variant.start", minimum=1)
    end = _int(item["end"], "structural_variant.end", minimum=start)
    sv_type = _text(item["sv_type"], "structural_variant.sv_type").upper()
    if sv_type not in _SV_TYPES:
        raise ConfigurationError(
            f"structural_variant.sv_type must be one of: {', '.join(_SV_TYPES)}."
        )
    sv = _call(
        _vc("structural_variant").StructuralVariant, contig, start, sv_type, end=end, genome=genome
    )
    window = _int(flank, "flank", minimum=0, maximum=100_000)
    candidates = _call(_vc("cryptic_exons").enumerate_from_structural_variant, sv, flank=window)
    rows = [_plain(candidate) for candidate in list(candidates)[:_MAX_CANDIDATES]]
    return _document(
        "enumerate_from_structural_variant",
        structural_variant={"contig": contig, "start": start, "end": end, "sv_type": sv_type},
        flank=window,
        candidates=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def varcode_apply_phase_resolver_to_effects(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    phased_vcf: str,
    sample: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict a variant's effects and resolve them with the phase in a VCF sample (varcode.phasing).

    ``phased_vcf`` must carry the phase of ``sample`` in its GT (for example ``0|1``) and PS fields. The resolver
    is ``VCFPhaseResolver``; effects are reported after varcode applies it.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    varcode = require("varcode", _VARCODE_EXTRA)
    item = _variant(genome, variant)
    vcf = _path(phased_vcf, "phased_vcf")
    name = _text(sample, "sample")
    collection = _call(varcode.load_vcf, str(vcf), genome=genome)
    resolver = _call(_vc("phasing").VCFPhaseResolver, collection, name)
    effects = _call(_vc("effects").predict_variant_effects, item)
    resolved = _call(varcode.apply_phase_resolver_to_effects, effects, resolver)
    rows = [_effect_record(effect) for effect in list(resolved)[:_MAX_CANDIDATES]]
    return _document(
        "apply_phase_resolver_to_effects",
        variant=_variant_record(item),
        sample=name,
        effects=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def varcode_build_exon_runs(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    flank: int = 50,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project a transcript's sense-strand exons onto its genomic layout (varcode.transcript_layout.build_exon_runs).

    The layout is the transcript with ``flank`` bases on each side (``GenomicLayout.from_transcript``).
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    transcript = _transcript(genome, transcript_id)
    layout = _layout(genome, transcript, flank)
    runs = _call(_vc("transcript_layout").build_exon_runs, transcript, layout)
    return _document(
        "build_exon_runs",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        flank=int(flank),
        runs=[_plain(run) for run in list(runs)[:_MAX_CANDIDATES]],
        cache_directory=info["cache_directory"],
    )


def varcode_enumerate_layout_splice_plans(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    flank: int = 50,
    max_plans: int = 64,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate graph-valid splice plans on a transcript layout (varcode.splice_graph.enumerate_layout_splice_plans)."""
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    transcript = _transcript(genome, transcript_id)
    layout = _layout(genome, transcript, flank)
    limit = _int(max_plans, "max_plans", minimum=1, maximum=10_000)
    plans = list(
        _call(
            _vc("splice_graph").enumerate_layout_splice_plans, transcript, layout, max_plans=limit
        )
    )
    return _document(
        "enumerate_layout_splice_plans",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        flank=int(flank),
        plans=[_plain(plan) for plan in plans],
        count=len(plans),
        cache_directory=info["cache_directory"],
    )


def varcode_realize_exon_path(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    kept_run_keys: Sequence[str] | None = None,
    flank: int = 50,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Splice selected exon runs of a transcript and translate from its mapped CDS start (varcode).

    ``kept_run_keys`` selects runs by the ``key`` values reported by ``varcode_build_exon_runs``; ``None`` keeps all.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    transcript = _transcript(genome, transcript_id)
    layout = _layout(genome, transcript, flank)
    keys = None
    if kept_run_keys is not None:
        keys = [
            _text(key, "kept_run_keys[]") for key in _sequence_of(kept_run_keys, "kept_run_keys")
        ]
    product = _call(_vc("transcript_layout").realize_exon_path, transcript, layout, keys)
    return _document(
        "realize_exon_path",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        product=_plain(product),
        cache_directory=info["cache_directory"],
    )


def varcode_realize_splice_plan(
    gtf: str,
    genome_fasta: str,
    transcript_id: str,
    *,
    plan_index: int = 0,
    flank: int = 50,
    scan_flank: int = 50,
    max_plans: int = 64,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Realise one enumerated splice plan (exon skip, retention or cryptic choices) on the transcript layout.

    ``plan_index`` selects a plan from ``varcode_enumerate_layout_splice_plans`` (the same ``max_plans``). The
    disrupted splice sites are computed with varcode's ``disrupted_splice_sites``.
    """
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    transcript = _transcript(genome, transcript_id)
    layout = _layout(genome, transcript, flank)
    graph = _vc("splice_graph")
    index = _int(plan_index, "plan_index", minimum=0, maximum=max_plans)
    limit = _int(max_plans, "max_plans", minimum=1, maximum=10_000)
    plans = list(_call(graph.enumerate_layout_splice_plans, transcript, layout, max_plans=limit))
    if index >= len(plans):
        raise ConfigurationError(
            f"plan_index {index} is out of range; the layout has {len(plans)} plans."
        )
    runs = _call(_vc("transcript_layout").build_exon_runs, transcript, layout)
    statuses = _call(graph.disrupted_splice_sites, transcript, layout, runs)
    product, unresolved = _call(
        _vc("transcript_layout").realize_splice_plan,
        transcript,
        layout,
        plans[index],
        statuses,
        scan_flank=_int(scan_flank, "scan_flank", minimum=0, maximum=10_000),
    )
    return _document(
        "realize_splice_plan",
        transcript_id=_identifier(transcript_id, "transcript_id"),
        plan_index=index,
        plan_count=len(plans),
        product=_plain(product),
        unresolved_mechanisms=_plain(unresolved),
        cache_directory=info["cache_directory"],
    )


__all__ = [
    "build_local_reference",
    "varcode_apply_germline_to_transcript",
    "varcode_apply_variant_to_transcript",
    "varcode_apply_variants_to_transcript",
    "varcode_candidates_from_effects",
    "varcode_classify_from_protein_diff",
    "varcode_codon_table_for_transcript",
    "varcode_default_germline_window",
    "varcode_infer_genome",
    "varcode_infer_genome_for_reference_name",
    "varcode_left_align_indels",
    "varcode_make_rna_outcome",
    "varcode_predict_germline_aware_effect",
    "varcode_predict_realized_effect",
    "varcode_predict_variant_effect_on_transcript",
    "varcode_predict_variant_effect_on_transcript_or_failure",
    "varcode_predict_variant_effects",
    "varcode_apply_phase_resolver_to_effects",
    "varcode_build_exon_runs",
    "varcode_enumerate_from_structural_variant",
    "varcode_enumerate_layout_splice_plans",
    "varcode_predict_structural_variant_effect",
    "varcode_random_variants",
    "varcode_realize_exon_path",
    "varcode_realize_splice_plan",
    "varcode_reference_base",
    "varcode_reference_range",
    "varcode_reference_selenocysteine",
]
_register_public_verbose(globals())
