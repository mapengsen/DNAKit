"""DnaCauldron Golden Gate and Type-2S functions exposed through ``dnakit.molecular_lab``.

Each function imports DnaCauldron (and Biopython) only when it runs. A missing
package raises ``BackendUnavailableError`` with the ``dnakit[bio-dnacauldron]``
install hint. Coordinates are 0-based with ``start`` inclusive and ``end``
exclusive. Every function returns a JSON-safe dictionary with ``status``,
``operation``, the result fields and ``provenance``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnacauldron"
_MAX_DNA = 200_000
_MAX_PARTS = 200
_MAX_FEATURES = 1000
_TYPE2S_ENZYMES = ("BsmBI", "BsaI", "BbsI", "AarI", "SapI")


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _dna_list(value: object, name: str) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of DNA strings.")
    if not 1 <= len(value) <= _MAX_PARTS:
        raise ConfigurationError(f"{name} must contain between 1 and {_MAX_PARTS} sequences.")
    texts = [_dna(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if sum(len(item) for item in texts) > _MAX_DNA:
        raise ConfigurationError(f"{name} exceeds {_MAX_DNA} bases in total.")
    return texts


def _enzyme(value: object, name: str) -> str:
    if value not in _TYPE2S_ENZYMES:
        allowed = ", ".join(_TYPE2S_ENZYMES)
        raise ConfigurationError(f"{name} must be one of {allowed}.")
    return str(value)


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnacauldron")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnacauldron", "version": installed, **notes}


def _linear_records(texts: list[str], dnacauldron: Any) -> list[Any]:
    records = []
    for index, text in enumerate(texts, start=1):
        record = dnacauldron.sequence_to_biopython_record(
            text, id=f"part_{index}", name=f"part_{index}"
        )
        record.annotations["topology"] = "linear"
        records.append(record)
    return records


def autoselect_golden_gate_enzyme(
    parts: Sequence[str],
    *,
    enzymes: Sequence[str] = _TYPE2S_ENZYMES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pick the Type-2S enzyme that a set of Golden Gate part sequences was designed for.

    The enzyme is the candidate that cuts each part into the expected number of
    pieces with DnaCauldron's site counting. ``enzymes`` restricts the candidates
    to the given subset of ``BsmBI``, ``BsaI``, ``BbsI``, ``AarI`` and ``SapI``.
    """

    dnacauldron = require("dnacauldron", _EXTRA)
    from dnacauldron.biotools.autoselect_enzyme import autoselect_enzyme

    texts = _dna_list(parts, "parts")
    if isinstance(enzymes, (str, bytes)) or not isinstance(enzymes, Sequence) or not enzymes:
        raise ConfigurationError("enzymes must be a non-empty list of enzyme names.")
    candidates = tuple(_enzyme(item, "enzymes") for item in enzymes)
    if len(set(candidates)) != len(candidates):
        raise ConfigurationError("enzymes must not repeat a name.")
    records = _linear_records(texts, dnacauldron)
    try:
        chosen = autoselect_enzyme(records, enzymes=candidates)
    except (AssertionError, IndexError, KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DnaCauldron could not select an enzyme: {exc}") from exc
    if not chosen:
        raise ConfigurationError("DnaCauldron found no enzyme for these parts.")
    return {
        "status": "completed",
        "operation": "autoselect_golden_gate_enzyme",
        "part_count": len(texts),
        "candidate_enzymes": list(candidates),
        "enzyme": str(chosen),
        "provenance": _provenance(method="dnacauldron.biotools.autoselect_enzyme"),
    }


def list_type2s_digestion_overhangs(
    sequences: Sequence[str],
    *,
    enzyme: str = "auto",
    parts_only: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the sticky overhangs that Type-2S digestion creates in DNA sequences.

    ``enzyme`` is ``"auto"`` or one of the Type-2S enzymes ``BsmBI``, ``BsaI``,
    ``BbsI``, ``AarI`` and ``SapI``. With ``parts_only`` true, only overhangs on
    part fragments are kept. The result is sorted and unique; DnaCauldron's
    ``"None"`` placeholder for missing overhangs is removed and counted.
    """

    dnacauldron = require("dnacauldron", _EXTRA)
    from dnacauldron.utils import list_digestion_overhangs

    texts = _dna_list(sequences, "sequences")
    if enzyme != "auto":
        enzyme = _enzyme(enzyme, "enzyme")
    if not isinstance(parts_only, bool):
        raise ConfigurationError("parts_only must be a boolean.")
    records = _linear_records(texts, dnacauldron)
    try:
        raw = [
            str(item)
            for item in list_digestion_overhangs(records, enzyme=enzyme, parts_only=parts_only)
        ]
    except (AssertionError, IndexError, KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DnaCauldron could not list overhangs: {exc}") from exc
    overhangs = sorted({item for item in raw if item != "None"})
    return {
        "status": "completed",
        "operation": "list_type2s_digestion_overhangs",
        "enzyme": enzyme,
        "parts_only": parts_only,
        "sequence_count": len(texts),
        "overhangs": overhangs,
        "overhang_count": len(overhangs),
        "none_placeholder_count": sum(1 for item in raw if item == "None"),
        "provenance": _provenance(method="dnacauldron.utils.list_digestion_overhangs"),
    }


def crop_dna_record_with_features(
    sequence: str,
    *,
    start: int,
    end: int,
    features: Sequence[Mapping[str, Any]] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """Crop a DNA record and keep every feature that overlaps the cropped segment.

    ``features`` is a list of mappings with ``type``, ``start``, ``end`` and an
    optional ``label``. DnaCauldron's saddling crop decides which features are
    kept and cuts the sequence. DnaCauldron 2.0.12 leaves partly included
    features in absolute coordinates, so this wrapper recomputes their
    coordinates by clipping to the window and shifting by ``start``.
    """

    require("Bio.SeqFeature", _EXTRA)
    dnacauldron = require("dnacauldron", _EXTRA)
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from dnacauldron.biotools.record_operations import crop_record_with_saddling_features

    text = _dna(sequence, "sequence")
    left = _integer(start, "start", minimum=0, maximum=len(text) - 1)
    right = _integer(end, "end", minimum=left + 1, maximum=len(text))
    if isinstance(features, (str, bytes)) or not isinstance(features, Sequence):
        raise ConfigurationError("features must be a list of feature mappings.")
    if len(features) > _MAX_FEATURES:
        raise ConfigurationError(f"features accepts at most {_MAX_FEATURES} entries.")
    record = dnacauldron.sequence_to_biopython_record(text, id="record", name="record")
    record.annotations["topology"] = "linear"
    spans: list[tuple[int, int, str, str]] = []
    for index, item in enumerate(features):
        if not isinstance(item, Mapping) or not isinstance(item.get("type"), str):
            raise ConfigurationError("Each feature needs a string 'type' key.")
        feature_type = item["type"].strip()
        if not feature_type or len(feature_type) > 64:
            raise ConfigurationError("Feature type must be 1 to 64 characters.")
        f_start = _integer(
            item.get("start"), f"features[{index}].start", minimum=0, maximum=len(text) - 1
        )
        f_end = _integer(
            item.get("end"), f"features[{index}].end", minimum=f_start + 1, maximum=len(text)
        )
        label = "" if item.get("label") is None else str(item["label"])
        qualifiers: dict[str, list[str]] = {"dnakit_index": [str(index)]}
        if label:
            qualifiers["label"] = [label]
        record.features.append(
            SeqFeature(FeatureLocation(f_start, f_end), type=feature_type, qualifiers=qualifiers)
        )
        spans.append((f_start, f_end, feature_type, label))
    try:
        cropped = crop_record_with_saddling_features(record, left, right)
    except (AssertionError, IndexError, KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DnaCauldron could not crop the record: {exc}") from exc
    rows = []
    for feature in cropped.features:
        f_start, f_end, feature_type, label = spans[int(feature.qualifiers["dnakit_index"][0])]
        rows.append(
            {
                "type": feature_type,
                "start": max(f_start, left) - left,
                "end": min(f_end, right) - left,
                "label": label,
            }
        )
    rows.sort(key=lambda row: (row["start"], row["end"], row["type"], row["label"]))
    cropped_sequence = str(cropped.seq)
    return {
        "status": "completed",
        "operation": "crop_dna_record_with_features",
        "source_length": len(text),
        "start": left,
        "end": right,
        "sequence": cropped_sequence,
        "sequence_length": len(cropped_sequence),
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(
            method="dnacauldron.biotools.crop_record_with_saddling_features",
            coordinates="clipped and shifted by dnakit; upstream leaves partial features absolute",
        ),
    }


__all__ = [
    "autoselect_golden_gate_enzyme",
    "crop_dna_record_with_features",
    "list_type2s_digestion_overhangs",
]

_register_public_verbose(globals(), tuple(__all__))
