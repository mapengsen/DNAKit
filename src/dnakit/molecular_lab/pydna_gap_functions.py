"""pydna wrappers: codon adaptation, SnapGene history, gel images, gel size, BLAST and NEB Tm.

Packages are imported only when a function runs. Two wrappers use the network through pydna and
are documented as network wrappers: ``pydna_blast_best_nt_hit`` (NCBI BLAST, via pydna's
``sequence_picker.genbank_accession``) and ``pydna_neb_melting_temperature`` (the NEB Tm API, via
pydna's ``tm.tm_neb``). Both reach the network only through the pydna module attribute, so tests
replace it. Gel images are written as new PNG files. Output files are never overwritten.

``pydna_codon_adaptation_index`` computes the index from explicit codon weights. The organism table
route needs the cai2 package, which pydna imports before it computes the index.
"""

from __future__ import annotations

import importlib
import importlib.util
import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_GEL_HINT = 'Install the gel extra of pydna: python -m pip install "pydna[gel]" (Pillow and SciPy).'
_CAI2_HINT = (
    'Install cai2 through the express extra of pydna: python -m pip install "pydna[express]".'
)
_REQUESTS_HINT = (
    'Install the download extra of pydna: python -m pip install "pydna[download]" (requests).'
)

_EXTRA = "bio-assembly"
_DNA = re.compile(r"^[ACGT]+$")
_PRODCODE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
_UNSCORED_CODONS = frozenset({"ATG", "TGG", "TAA", "TAG", "TGA"})
_LADDERS = ("GeneRuler_1kb_plus", "HI_LO_DNA_MARKER")
_MAX_CODING_BP = 300_000
_MAX_BLAST_BP = 20_000
_MAX_PRIMER = 200
_MAX_LANES = 20
_MAX_BANDS = 50
_MAX_FRAGMENT = 10_000_000


# Validation and conversion helpers.


def _optional(module: str, hint: str) -> Any:
    """Import a module that a pydna feature needs beyond the bio-assembly extra."""

    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(f"{module} could not be imported: {exc}", hint=hint) from exc


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not minimum <= float(value) <= maximum:
        raise ConfigurationError(f"{name} must be in [{minimum}, {maximum}].")
    return float(value)


def _choice(value: object, name: str, options: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {list(options)}; got {value!r}.")
    return value


def _dna(value: object, name: str, *, maximum: int) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} has {len(text)} bases; the limit is {maximum}.")
    if not _DNA.match(text):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _local_file(value: object, name: str, suffix: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a local file path.")
    if "://" in value:
        raise ConfigurationError(f"{name} must be a local path; remote URLs are not read.")
    path = input_path(value)
    if path.suffix.lower() != suffix:
        raise ConfigurationError(f"{name} must be a file ending in {suffix}; got {path.name!r}.")
    return path


def _new_png(value: object) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("output_path must be a path string.")
    target = Path(value).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    if target.exists():
        raise ConfigurationError(
            f"Output file already exists and will not be overwritten: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _provenance(convention: str) -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    try:
        package = f"pydna {version('pydna')}"
    except PackageNotFoundError:
        package = "pydna"
    return {"package": package, "convention": convention}


def _plain(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


# Codon adaptation.


def _codon_weights(weights: object) -> dict[str, float]:
    if not isinstance(weights, dict) or not weights:
        raise ConfigurationError("weights must be a non-empty mapping from codon to weight.")
    table: dict[str, float] = {}
    for codon, value in weights.items():
        if not isinstance(codon, str) or len(codon) != 3 or not _DNA.match(codon.upper()):
            raise ConfigurationError(f"weights key {codon!r} is not a DNA codon.")
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ConfigurationError(f"weight for {codon} must be a positive finite number.")
        table[codon.upper()] = float(value)
    return table


def pydna_codon_adaptation_index(
    seq: str,
    *,
    weights: dict[str, float] | None = None,
    organism: str = "sce",
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the codon adaptation index (CAI) of a coding sequence.

    With ``weights`` (codon to relative adaptiveness, each in (0, 1]), the index is the geometric
    mean of the weights of the scored codons, as defined by Sharp and Li (1987). ATG, TGG and stop
    codons are not scored. Missing codons raise ConfigurationError. Without ``weights``, pydna's
    organism table is used through pydna.utils.cai, which needs the cai2 package.
    """

    text = _dna(seq, "seq", maximum=_MAX_CODING_BP)
    if len(text) < 3 or len(text) % 3:
        raise ConfigurationError("seq must be a whole number of codons (at least one).")
    codons = [text[index : index + 3] for index in range(0, len(text), 3)]
    if weights is None:
        if importlib.util.find_spec("cai2") is None:
            raise BackendUnavailableError(
                "cai2 is not installed, so pydna's organism tables cannot be used.",
                hint=_CAI2_HINT + " Or pass explicit codon weights to compute the index.",
            )
        utils = require("pydna.utils", _EXTRA)
        table = require("pydna.codon", _EXTRA).weights
        key = _choice(organism, "organism", tuple(table))
        value = round(float(utils.cai(text, organism=key)), 3)
        return {
            "status": "completed",
            "operation": "codon_adaptation_index",
            "method": f"pydna.utils.cai table for {key} (cai2)",
            "cai": value,
            "codon_count": len(codons),
            "provenance": _provenance("CAI from pydna's organism table through cai2"),
        }
    table = _codon_weights(weights)
    scored = [codon for codon in codons if codon not in _UNSCORED_CODONS]
    if not scored:
        raise ConfigurationError(
            "seq has no codon left to score after ATG, TGG and stop codons are removed."
        )
    missing = sorted({codon for codon in scored if codon not in table})
    if missing:
        raise ConfigurationError(
            f"weights lack {len(missing)} codon(s) used in seq, for example {missing[:5]}."
        )
    mean_log = sum(math.log(table[codon]) for codon in scored) / len(scored)
    return {
        "status": "completed",
        "operation": "codon_adaptation_index",
        "method": (
            "explicit weights; geometric mean (Sharp and Li, 1987); ATG, TGG and stops not scored"
        ),
        "cai": round(math.exp(mean_log), 3),
        "codon_count": len(codons),
        "scored_codon_count": len(scored),
        "provenance": _provenance(
            "CAI computed from the explicit weights; rounded to 3 decimals as pydna does"
        ),
    }


# SnapGene history files.


def _feature_summary(feature: Any) -> dict[str, Any]:
    location = feature.location
    return {
        "type": str(feature.type),
        "start": int(location.start),
        "end": int(location.end),
        "strand": location.strand,
        "label": str(feature.qualifiers.get("label", [""])[0]) if feature.qualifiers else "",
    }


def pydna_parse_snapgene_history(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a local SnapGene .dna file with pydna, resolving its cloning history when it has one.

    The result gives the root sequence (name, length, topology, bases), its features, and whether
    the file carries a cloning history. A file without a history returns the root sequence only.
    """

    source = _local_file(path, "path", ".dna")
    module = require("pydna.snapgene_history_parser", _EXTRA)
    sgffp = require("sgffp", _EXTRA)
    record = module.parse_snapgene_history(str(source))
    features = [_feature_summary(feature) for feature in record.features]
    history = getattr(record, "source", None)
    return {
        "status": "completed",
        "operation": "parse_snapgene_history",
        "name": str(record.name),
        "length": len(record),
        "circular": bool(record.circular),
        "sequence": str(record.seq).upper(),
        "feature_count": len(features),
        "features": features,
        "has_history": bool(sgffp.SgffReader.from_file(str(source)).has_history),
        "source_type": None if history is None else type(history).__name__,
        "provenance": _provenance("SnapGene .dna parsed by pydna.snapgene_history_parser"),
    }


# Gel images and the gel size model.


def _ladder(name: str) -> Any:
    ladders = require("pydna.ladders", _EXTRA)
    return getattr(ladders, _choice(name, "ladder", _LADDERS))


def _size_model(name: str) -> Any:
    gel = require("pydna.gel", _EXTRA)
    return gel.interpolator(mwstd=_ladder(name))


def pydna_predict_gel_relative_front(
    lengths: list[int],
    *,
    ladder: str = "GeneRuler_1kb_plus",
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the relative migration front (Rf) of DNA fragments on an agarose gel.

    pydna's gel size model is a natural cubic spline through the Rf values of a built-in ladder,
    evaluated at each fragment length. Lengths outside the ladder range have no prediction (None).
    No figure is drawn; use pydna_draw_gel for an image.
    """

    if not isinstance(lengths, list) or not lengths or len(lengths) > 500:
        raise ConfigurationError("lengths must be a list of 1 to 500 fragment lengths in bp.")
    sizes = [_integer(item, "length", 1, _MAX_FRAGMENT) for item in lengths]
    ladder_name = _choice(ladder, "ladder", _LADDERS)
    _optional("scipy.interpolate", _GEL_HINT)
    model = _size_model(ladder_name)
    bands = _ladder(ladder_name)[::-1]  # the same knot order pydna.gel.interpolator uses
    predicted: list[float | None] = []
    for size in sizes:
        value = float(model(size))
        predicted.append(value if math.isfinite(value) else None)
    return {
        "status": "completed",
        "operation": "predict_gel_relative_front",
        "ladder": ladder_name,
        "ladder_sizes": [len(band) for band in bands],
        "ladder_relative_fronts": [float(band.rf) for band in bands],
        "lengths": sizes,
        "relative_front": predicted,
        "provenance": _provenance(
            "natural cubic spline through the ladder's Rf values (pydna.gel.interpolator)"
        ),
    }


def pydna_draw_gel(
    output_path: str,
    *,
    lanes: list[list[int]] | None = None,
    ladder: str = "GeneRuler_1kb_plus",
    gel_length: int = 600,
    margin: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw an agarose gel image from fragment lengths and write it as a new PNG.

    Each lane is a list of fragment lengths in bp. Bands are placed with pydna's gel size model for
    the chosen ladder, so every length must lie within the ladder's range. With ``lanes`` left as
    None, the ladder itself is drawn as one lane (pydna's default).
    """

    target = _new_png(output_path)
    name = _choice(ladder, "ladder", _LADDERS)
    _optional("PIL.Image", _GEL_HINT)
    _optional("scipy.interpolate", _GEL_HINT)
    height = _integer(gel_length, "gel_length", 100, 5000)
    offset = _integer(margin, "margin", 0, height - 10)
    model = _size_model(name)
    low, high = float(model.x[0]), float(model.x[-1])
    gel = require("pydna.gel", _EXTRA)
    fakeseq = require("pydna.fakeseq", _EXTRA).FakeSeq
    if lanes is None:
        lane_sizes = [[len(fragment) for fragment in _ladder(name)]]
    else:
        if not isinstance(lanes, list) or not 1 <= len(lanes) <= _MAX_LANES:
            raise ConfigurationError(f"lanes must be a list of 1 to {_MAX_LANES} lanes.")
        lane_sizes = []
        for lane in lanes:
            if not isinstance(lane, list) or not 1 <= len(lane) <= _MAX_BANDS:
                raise ConfigurationError(f"each lane must list 1 to {_MAX_BANDS} fragment lengths.")
            sizes = [_integer(item, "band length", 1, _MAX_FRAGMENT) for item in lane]
            outside = [size for size in sizes if not low <= size <= high]
            if outside:
                raise ConfigurationError(
                    f"band lengths {outside[:5]} are outside the {name} range "
                    f"{int(low)}-{int(high)} bp."
                )
            lane_sizes.append(sizes)
    samples = [[fakeseq(size) for size in lane] for lane in lane_sizes]
    image = gel.gel(samples=samples, gel_length=height, margin=offset, interpolator=model)
    image.save(target, format="PNG")
    return {
        "path": str(target),
        "data": {
            "ladder": name,
            "lanes": lane_sizes,
            "lane_count": len(lane_sizes),
            "band_count": sum(len(lane) for lane in lane_sizes),
            "gel_length": height,
            "margin": offset,
            "image_width": int(image.width),
            "image_height": int(image.height),
            "provenance": _provenance(
                "pydna.gel.gel drawn with Pillow; bands placed by the ladder's Rf spline"
            ),
        },
    }


# Network wrappers.


def pydna_blast_best_nt_hit(
    sequence: str, *, name: str = "query", verbose: bool = True
) -> dict[str, Any]:
    """Find the best NCBI nt hit for a DNA sequence with BLAST and return the matched region.

    Network wrapper. pydna's ``sequence_picker.genbank_accession`` submits a megablast search of the
    sequence against NCBI's nt database and returns the best hit region as a record. It looks up no
    accession, despite its name. Requires network access and the NCBI BLAST service.
    """

    text = _dna(sequence, "sequence", maximum=_MAX_BLAST_BP)
    label = _text_name(name)
    picker = require("pydna.sequence_picker", _EXTRA)
    record = picker.genbank_accession(text)
    return {
        "service": "NCBI BLAST (blastn, nt, megablast)",
        "source": "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
        "result": {
            "query_name": label,
            "sequence": str(record.seq).upper(),
            "length": len(record),
            "circular": bool(record.circular),
            "description": str(record.description),
        },
    }


def _text_name(value: object) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 128:
        raise ConfigurationError("name must be a non-empty string of at most 128 characters.")
    return value


def pydna_neb_melting_temperature(
    primer: str,
    *,
    conc: float = 0.5,
    prodcode: str = "q5-0",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the melting temperature (degrees Celsius) of a primer from the NEB Tm calculator.

    Network wrapper. pydna's ``tm.tm_neb`` posts the primer to https://tmapi.neb.com/tm. ``conc``
    is the primer concentration in micromolar and ``prodcode`` is an NEB product code, such as
    ``q5-0`` for Q5 polymerase.
    """

    text = _dna(primer, "primer", maximum=_MAX_PRIMER)
    concentration = _number(conc, "conc", minimum=0.001, maximum=100.0)
    if not isinstance(prodcode, str) or not _PRODCODE.match(prodcode):
        raise ConfigurationError("prodcode must be 1 to 64 letters, digits, '-' or '_'.")
    _optional("requests", _REQUESTS_HINT)
    tm_module = require("pydna.tm", _EXTRA)
    value = tm_module.tm_neb(text, conc=concentration, prodcode=prodcode)
    return {
        "service": "NEB Tm calculator",
        "source": "https://tmapi.neb.com/tm",
        "result": {
            "primer": text,
            "conc_uM": concentration,
            "prodcode": prodcode,
            "tm_celsius": float(value),
        },
    }


__all__ = [
    "pydna_blast_best_nt_hit",
    "pydna_codon_adaptation_index",
    "pydna_draw_gel",
    "pydna_neb_melting_temperature",
    "pydna_parse_snapgene_history",
    "pydna_predict_gel_relative_front",
]

_register_public_verbose(globals())
