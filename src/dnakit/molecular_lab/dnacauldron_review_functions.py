"""Review wrappers for DnaCauldron Type-2S (Golden Gate) assembly helpers (G18n).

Each function imports DnaCauldron (and Biopython) only when it runs. A missing package
raises ``BackendUnavailableError`` with the ``dnakit[bio-dnacauldron]`` install hint.
Parts and vectors are DNA strings of A, C, G and T; they are upper-cased. Every
function returns a JSON-safe dictionary with ``status``, ``operation``, the result
fields and ``provenance``.

Two upstream quirks are handled explicitly. DnaCauldron's backbone detection
(``record_contains_backbone``) always uses BsmBI, whatever enzyme is passed, so the
backbone functions here accept only BsmBI. Upstream asserts or raises bare errors on
inputs with the wrong site layout; those are reported as ``ConfigurationError``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnacauldron"
_TYPE2S = ("BbsI", "AarI", "BsaI", "BsmBI", "SapI")
_MAX_PART_DNA = 20_000
_MAX_PARTS = 200
_MAX_TOTAL_DNA = 200_000
_MAX_PRODUCTS = 200
_MAX_SUBSTITUTIONS = 8


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _choice(value: object, name: str, allowed: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value not in allowed:
        options = ", ".join(allowed)
        raise ConfigurationError(f"{name} must be one of: {options}.")
    return value


def _dna(value: object, name: str, *, maximum: int = _MAX_PART_DNA) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 80:
        raise ConfigurationError(f"{name} must be a non-empty name of at most 80 characters.")
    return value.strip()


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnacauldron")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnacauldron", "version": installed, **notes}


def _run(action: Any, message: str) -> Any:
    mix_error = require("dnacauldron.AssemblyMix", _EXTRA).AssemblyMixError
    try:
        return action()
    except (
        AssertionError,
        AttributeError,
        IndexError,
        RuntimeError,
        TypeError,
        ValueError,
        mix_error,
    ) as exc:
        detail = str(exc) or type(exc).__name__
        raise ConfigurationError(f"{message}: {detail}") from exc


def _record(name: str, text: str, topology: str) -> Any:
    seq_class = require("Bio.Seq", _EXTRA).Seq
    record_class = require("Bio.SeqRecord", _EXTRA).SeqRecord
    record = record_class(seq_class(text), id=name, name=name, annotations={"molecule_type": "DNA"})
    require("dnacauldron.biotools", _EXTRA).set_record_topology(record, topology)
    return record


def _named(value: object, name: str, topology: str) -> list[Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of name to DNA string.")
    if not 1 <= len(value) <= _MAX_PARTS:
        raise ConfigurationError(f"{name} must contain between 1 and {_MAX_PARTS} items.")
    records = []
    total = 0
    for key, text in value.items():
        label = _name(key, f"{name} key")
        dna = _dna(text, f"{name}[{label!r}]")
        total += len(dna)
        records.append(_record(label, dna, topology))
    if total > _MAX_TOTAL_DNA:
        raise ConfigurationError(f"{name} exceeds {_MAX_TOTAL_DNA} bases in total.")
    return sorted(records, key=lambda item: item.id)


def _product_rows(records: list[Any], *, listed: int) -> tuple[list[dict[str, Any]], int]:
    rows = sorted(
        (
            {"sequence": str(item.seq).upper(), "length": len(item), "circular": True}
            for item in records
        ),
        key=lambda row: (row["length"], row["sequence"]),
    )
    return rows[:listed], len(rows)


def _type2s(value: object, name: str) -> str:
    return _choice(value, name, _TYPE2S)


def simulate_type2s_assembly(
    parts: Mapping[str, str],
    *,
    enzyme: str = "auto",
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate a Type-2S (Golden Gate) assembly of linear parts into circular products.

    ``parts`` maps part names to linear DNA strings whose ends carry the overhangs.
    ``enzyme`` is one of BbsI, AarI, BsaI, BsmBI or SapI, or ``"auto"`` to let
    DnaCauldron pick the enzyme that cuts each part about twice. Products are the
    circular assemblies DnaCauldron finds, sorted by length and sequence.
    """

    mix_module = require("dnacauldron.AssemblyMix", _EXTRA)
    biotools = require("dnacauldron.biotools", _EXTRA)
    records = _named(parts, "parts", "linear")
    if enzyme == "auto":
        chosen = _run(lambda: biotools.autoselect_enzyme(records), "DnaCauldron could not pick")
    else:
        chosen = _type2s(enzyme, "enzyme")
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    products = _run(
        lambda: list(
            mix_module.generate_type2s_restriction_mix(
                parts=records, enzyme=chosen
            ).compute_circular_assemblies()
        ),
        "DnaCauldron could not assemble the parts",
    )
    rows, total = _product_rows(products, listed=listed)
    return {
        "status": "completed",
        "operation": "simulate_type2s_assembly",
        "enzyme": chosen,
        "part_names": [item.id for item in records],
        "circular_product_count": total,
        "products": rows,
        "provenance": _provenance(method="dnacauldron.AssemblyMix.generate_type2s_restriction_mix"),
    }


def clone_into_donor_vector(
    donor_sequence: str,
    insert_sequence: str,
    *,
    enzyme: str = "BsmBI",
    verbose: bool = True,
) -> dict[str, Any]:
    """Clone an insert into a donor vector by swapping its Type-2S stuffer.

    The donor is a circular vector with one enzyme-site pair around the stuffer it
    holds. The insert is a circular plasmid whose site pair flanks the insert, so its
    digest leaves exactly one site-free fragment. The insert must carry the same two
    overhangs as the donor's stuffer; otherwise upstream fails and the call raises
    ``ConfigurationError``. The product is the donor with its stuffer replaced by the
    insert.
    """

    swap = require("dnacauldron.utils.insert_parts_on_backbones", _EXTRA).swap_donor_vector_part
    chosen = _type2s(enzyme, "enzyme")
    donor = _record(
        "donor", _dna(donor_sequence, "donor_sequence", maximum=_MAX_TOTAL_DNA), "circular"
    )
    insert = _record(
        "insert", _dna(insert_sequence, "insert_sequence", maximum=_MAX_TOTAL_DNA), "circular"
    )
    product = _run(
        lambda: swap(donor_vector=donor, insert=insert, enzyme=chosen),
        "DnaCauldron could not clone the insert into the donor",
    )
    sequence = str(product.seq).upper()
    return {
        "status": "completed",
        "operation": "clone_into_donor_vector",
        "enzyme": chosen,
        "product_sequence": sequence,
        "product_length": len(sequence),
        "provenance": _provenance(
            method="dnacauldron.utils.insert_parts_on_backbones.swap_donor_vector_part"
        ),
    }


def assign_parts_to_backbones(
    parts: Mapping[str, str],
    backbones: Mapping[str, str],
    *,
    min_backbone_length: int = 500,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign each circular Golden Gate part to the backbone with matching overhangs.

    ``parts`` and ``backbones`` map names to circular DNA strings carrying BsmBI sites.
    A backbone must hold one stuffer, and its overhangs are matched to the part's
    overhangs (in either orientation). Each part is reported as ``assigned`` (with the
    product sequence), ``no_backbone`` or ``error``. Results are sorted by part name.
    """

    ipb = require("dnacauldron.utils.insert_parts_on_backbones", _EXTRA).insert_parts_on_backbones
    part_records = _named(parts, "parts", "circular")
    backbone_records = _named(backbones, "backbones", "circular")
    minimum = _integer(min_backbone_length, "min_backbone_length", minimum=1, maximum=100_000)
    choices = _run(
        lambda: ipb(
            part_records,
            backbone_records,
            enzyme="BsmBI",
            min_backbone_length=minimum,
            process_parts_with_backbone=False,
            default_backbone_choice=None,
        ),
        "DnaCauldron could not assign the parts",
    )
    rows = []
    for choice in choices:
        backbone = choice.backbone_record
        backbone_name = backbone.id if hasattr(backbone, "id") else None
        final = choice.final_record
        error = None if choice.error is None else (str(choice.error) or "backbone error")
        if error is not None:
            state = "error"
        elif final is not None:
            state = "assigned"
        else:
            state = "no_backbone"
        rows.append(
            {
                "part": choice.record.id,
                "status_per_part": state,
                "already_on_backbone": bool(choice.already_on_backbone),
                "backbone": backbone_name,
                "product_sequence": None if final is None else str(final.seq).upper(),
                "product_length": None if final is None else len(final),
                "error": error,
            }
        )
    rows.sort(key=lambda row: row["part"])
    return {
        "status": "completed",
        "operation": "assign_parts_to_backbones",
        "enzyme": "BsmBI",
        "min_backbone_length": minimum,
        "part_count": len(rows),
        "assigned_count": sum(1 for row in rows if row["status_per_part"] == "assigned"),
        "assignments": rows,
        "provenance": _provenance(
            method="dnacauldron.utils.insert_parts_on_backbones.insert_parts_on_backbones"
        ),
    }


def substitute_type2s_overhangs(
    circular_sequence: str,
    substitutions: Mapping[str, str],
    *,
    enzyme: str = "BsmBI",
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace Type-2S overhangs on a circular vector, which moves a part in an assembly standard.

    ``substitutions`` maps each existing 4-base overhang to its new 4-base overhang.
    Only overhangs that occur as sticky ends of the vector's digest are applied, and
    they are reported in ``applied_overhangs``. The circular product is returned as
    ``product_sequence``.
    """

    utils = require("dnacauldron.utils.utils", _EXTRA)
    generate = require("dnacauldron.AssemblyMix", _EXTRA).generate_type2s_restriction_mix
    chosen = _type2s(enzyme, "enzyme")
    if not isinstance(substitutions, Mapping) or not 1 <= len(substitutions) <= _MAX_SUBSTITUTIONS:
        raise ConfigurationError(
            f"substitutions must map 1 to {_MAX_SUBSTITUTIONS} overhangs to new overhangs."
        )
    table = {}
    for old, new in substitutions.items():
        old_text = _dna(old, "substitution key", maximum=4)
        new_text = _dna(new, "substitution value", maximum=4)
        if len(old_text) != 4 or len(new_text) != 4:
            raise ConfigurationError("overhangs must be exactly 4 bases long.")
        table[old_text] = new_text
    vector = _record(
        "vector", _dna(circular_sequence, "circular_sequence", maximum=_MAX_TOTAL_DNA), "circular"
    )
    product = _run(
        lambda: utils.substitute_overhangs(vector, table, enzyme=chosen),
        "DnaCauldron could not substitute the overhangs",
    )
    mix = _run(lambda: generate(parts=[vector], enzyme=chosen), "DnaCauldron could not digest")
    ends = {
        str(end).upper()
        for fragment in mix.fragments
        if not fragment.is_reversed
        for end in (fragment.seq.left_end, fragment.seq.right_end)
    }
    sequence = str(product.seq).upper()
    return {
        "status": "completed",
        "operation": "substitute_type2s_overhangs",
        "enzyme": chosen,
        "requested": sorted(table),
        "applied_overhangs": sorted(key for key in table if key in ends),
        "product_sequence": sequence,
        "product_length": len(sequence),
        "provenance": _provenance(method="dnacauldron.utils.utils.substitute_overhangs"),
    }


def record_has_backbone(
    circular_sequence: str,
    *,
    min_backbone_length: int = 500,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tell whether a circular BsmBI construct still carries a backbone besides its insert.

    DnaCauldron cuts the vector with BsmBI and looks for one site-free fragment (the
    insert). The record has a backbone when the rest of the circle is longer than
    ``min_backbone_length``. The result reports ``has_backbone``, ``insert_length`` and
    ``backbone_length``. The enzyme is fixed to BsmBI, as upstream fixes it.
    """

    module = require("dnacauldron.utils.insert_parts_on_backbones", _EXTRA)
    generate = require("dnacauldron.AssemblyMix", _EXTRA).generate_type2s_restriction_mix
    minimum = _integer(min_backbone_length, "min_backbone_length", minimum=1, maximum=100_000)
    record = _record(
        "construct",
        _dna(circular_sequence, "circular_sequence", maximum=_MAX_TOTAL_DNA),
        "circular",
    )
    has = _run(
        lambda: module.record_contains_backbone(record, min_backbone_length=minimum),
        "DnaCauldron could not digest the construct",
    )
    mix = _run(lambda: generate(parts=[record], enzyme="BsmBI"), "DnaCauldron could not digest")
    inserts = [fragment for fragment in mix.filtered_fragments if not fragment.is_reversed]
    if not inserts:
        raise ConfigurationError("no site-free insert fragment was found in the construct.")
    insert_length = len(inserts[0])
    return {
        "status": "completed",
        "operation": "record_has_backbone",
        "enzyme": "BsmBI",
        "min_backbone_length": minimum,
        "construct_length": len(record),
        "insert_length": insert_length,
        "backbone_length": len(record) - insert_length,
        "has_backbone": bool(has),
        "provenance": _provenance(
            method="dnacauldron.utils.insert_parts_on_backbones.record_contains_backbone"
        ),
    }


__all__ = [
    "assign_parts_to_backbones",
    "clone_into_donor_vector",
    "record_has_backbone",
    "simulate_type2s_assembly",
    "substitute_type2s_overhangs",
]

_register_public_verbose(globals(), tuple(__all__))
