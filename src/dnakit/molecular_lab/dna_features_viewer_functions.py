"""DNA Features Viewer annotation and graphical-translation helpers (optional).

dna_features_viewer and Biopython are imported only inside each function. The
complement, reverse-complement and translate helpers of the package duplicate
:mod:`dnakit.ops` and are not wrapped here.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require, sequence_text
from dnakit.exceptions import ConfigurationError

_MAX_SEQUENCE = 1_000_000
_MAX_TEXT = 512
_RESERVED_TAGS = frozenset({"seqrecord", "location", "feature_type", "margin"})


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version("dna-features-viewer")
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": "dna_features_viewer", "version": version, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _interval(start: object, end: object, length: int) -> tuple[int, int]:
    begin = _integer(start, "start", minimum=0, maximum=length - 1)
    stop = length if end is None else _integer(end, "end", minimum=1, maximum=length)
    if begin >= stop:
        raise ConfigurationError("start must be smaller than end.")
    return begin, stop


def _strand(value: object, allowed: tuple[int, ...]) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value not in allowed:
        raise ConfigurationError(f"strand must be one of {list(allowed)}.")
    return value


def _tags(qualifiers: Mapping[str, str] | None) -> dict[str, str]:
    if qualifiers is None:
        return {}
    if not isinstance(qualifiers, Mapping):
        raise ConfigurationError("qualifiers must be a mapping of strings to strings.")
    tags: dict[str, str] = {}
    for key, value in sorted(qualifiers.items(), key=lambda item: str(item[0])):
        if not isinstance(key, str) or not key.strip() or key in _RESERVED_TAGS:
            raise ConfigurationError(f"qualifier name is empty or reserved: {key!r}")
        if not isinstance(value, str) or len(value) > _MAX_TEXT or "\x00" in value:
            raise ConfigurationError(
                f"qualifier {key!r} must be a string of at most {_MAX_TEXT} characters."
            )
        tags[key] = value
    return tags


def _feature_row(feature: Any) -> dict[str, Any]:
    location = feature.location
    qualifiers = {}
    for key, value in sorted(feature.qualifiers.items()):
        if isinstance(value, (list, tuple)):
            value = ", ".join(str(item) for item in value)
        qualifiers[str(key)] = str(value)
    return {
        "type": str(feature.type),
        "start": int(location.start),
        "end": int(location.end),
        "strand": None if location.strand is None else int(location.strand),
        "qualifiers": qualifiers,
    }


def annotate_sequence_feature(
    sequence: str,
    *,
    start: int = 0,
    end: int | None = None,
    strand: int = 1,
    feature_type: str = "misc_feature",
    qualifiers: Mapping[str, str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate a DNA string with one feature using DNA Features Viewer.

    ``start`` is 0-based and inclusive; ``end`` is exclusive and defaults to the
    sequence length. ``strand`` is 1, -1 or 0. ``qualifiers`` become the
    feature's string qualifiers.
    """

    dfv = require("dna_features_viewer", "bio-dna-features-viewer")
    bio_seq = require("Bio.Seq", "bio-dna-features-viewer")
    bio_record = require("Bio.SeqRecord", "bio-dna-features-viewer")
    text = sequence_text(sequence)
    if len(text) > _MAX_SEQUENCE:
        raise ConfigurationError(f"sequence exceeds {_MAX_SEQUENCE} bases.")
    begin, stop = _interval(start, end, len(text))
    feature_strand = _strand(strand, (-1, 0, 1))
    if not isinstance(feature_type, str) or not feature_type.strip() or len(feature_type) > 64:
        raise ConfigurationError(
            "feature_type must be a non-empty string of at most 64 characters."
        )
    tags = _tags(qualifiers)
    record = bio_record.SeqRecord(bio_seq.Seq(text), id="sequence")
    dfv.annotate_biopython_record(
        record,
        location=(begin, stop, feature_strand),
        feature_type=feature_type,
        margin=0,
        **tags,
    )
    return {
        "status": "completed",
        "operation": "annotate_sequence_feature",
        "sequence_length": len(text),
        "feature": _feature_row(record.features[-1]),
        "provenance": _provenance(annotator="annotate_biopython_record"),
    }


def graphical_translation(
    sequence: str,
    *,
    start: int = 0,
    end: int | None = None,
    strand: int = 1,
    long_form: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Translate a subsegment of a DNA string in DNA Features Viewer's graphical order.

    For ``strand=-1`` the amino acids come from the reverse complement, and the
    returned string is mirrored so that it reads left to right under the forward
    sequence. ``long_form=True`` returns three-letter residue names.
    """

    biotools = require("dna_features_viewer.biotools", "bio-dna-features-viewer")
    text = sequence_text(sequence)
    begin, stop = _interval(start, end, len(text))
    feature_strand = _strand(strand, (-1, 1))
    if not isinstance(long_form, bool):
        raise ConfigurationError("long_form must be a boolean.")
    value = biotools.extract_graphical_translation(
        text, (begin, stop, feature_strand), long_form=long_form
    )
    translation: str | list[str] = list(value) if long_form else str(value)
    return {
        "status": "completed",
        "operation": "graphical_translation",
        "sequence_length": len(text),
        "start": begin,
        "end": stop,
        "strand": feature_strand,
        "long_form": long_form,
        "translation": translation,
        "residue_count": len(translation),
        "provenance": _provenance(translator="dna_features_viewer.biotools"),
    }


__all__ = ["annotate_sequence_feature", "graphical_translation"]

_register_public_verbose(globals(), tuple(__all__))
