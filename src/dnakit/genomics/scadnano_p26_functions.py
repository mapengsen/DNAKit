"""scadnano DNA nanostructure design operations as DNAKit functions (pack P26).

Every function works on a design given as a scadnano JSON mapping (the content of a ``.sc``
file). Functions that edit a design never change their argument: they load a private copy with
``scadnano.Design.from_scadnano_json_map``, apply the scadnano operation and return the new design
as a JSON mapping under the key ``"design"``. Strands are addressed by their 0-based position in
``design["strands"]`` and domains by their 0-based position in the strand; helix indices and offsets
are scadnano's own. Functions that write a file take an ``output_path`` that must not exist yet.

scadnano (https://github.com/UC-Davis-molecular-computing/scadnano-python-package) is imported only
when a function runs; install it with ``pip install "dnakit[dna-nanotechnology]"``. Ordinary tests
never use the network.
"""

from __future__ import annotations

import copy
import importlib.metadata
import json
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "dna-nanotechnology"
_MAX_TEXT = 50_000_000
_MAX_ITEMS = 200_000
_MAX_SEQUENCE = 1_000_000
_BACKEND_ERRORS = (
    ValueError,
    TypeError,
    KeyError,
    IndexError,
    AttributeError,
    AssertionError,
    ArithmeticError,
    NotImplementedError,
)
_PLATE_TYPES = (96, 384)
_MOD_LOCATIONS = ("5p", "3p", "internal")
_M13_VARIANTS = ("p7249", "p7560", "p8064", "p8634")
_BUILDER_STEPS = {
    "to": ("offset",),
    "move": ("delta",),
    "cross": ("helix", "offset", "move"),
    "loopout": ("helix", "length", "offset", "move"),
    "extension_3p": ("num_bases", "display_length", "display_angle"),
    "extension_5p": ("num_bases", "display_length", "display_angle"),
    "with_deletions": ("deletions",),
    "with_insertions": ("insertions",),
    "with_sequence": ("sequence",),
    "with_domain_sequence": ("sequence",),
    "with_name": ("name",),
    "with_label": ("label",),
    "with_domain_name": ("name",),
    "as_scaffold": (),
    "as_circular": (),
    "with_vendor_fields": ("scale", "purification", "plate", "well"),
}


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _sc() -> Any:
    return require("scadnano", _EXTRA)


def _provenance() -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version("scadnano")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "scadnano", "package_version": installed}


@contextmanager
def _guard(action: str) -> Iterator[None]:
    """Turn scadnano's own validation errors into ``ConfigurationError``."""

    try:
        yield
    except ConfigurationError:
        raise
    except _BACKEND_ERRORS as exc:
        raise ConfigurationError(f"scadnano could not {action}: {exc}") from exc


def _int(value: object, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, allowed: Sequence[Any], name: str) -> Any:
    if value not in allowed or isinstance(value, bool):
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _text(value: object, name: str, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value):
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def _sequence_text(value: object, name: str) -> str:
    text = _text(value, name)
    if len(text) > _MAX_SEQUENCE:
        raise ConfigurationError(f"{name} may be at most {_MAX_SEQUENCE} characters.")
    return text


def _new_output(path: object) -> Path:
    if not isinstance(path, (str, Path)) or not str(path):
        raise ConfigurationError("output_path must be a path.")
    target = Path(path)
    if target.exists():
        raise ConfigurationError(f"output_path already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"the directory of output_path does not exist: {target.parent}")
    return target


def _write(path: object, text: str) -> str | None:
    if path is None:
        return None
    target = _new_output(path)
    target.write_text(text, encoding="utf-8")
    return str(target)


def _design_json(design: object) -> dict[str, Any]:
    if not isinstance(design, Mapping):
        raise ConfigurationError("design must be a scadnano JSON mapping.")
    try:
        text = json.dumps(design, allow_nan=False)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError("design must contain only JSON values.") from exc
    if len(text) > _MAX_TEXT:
        raise ConfigurationError("design JSON is larger than 50 MB.")
    return dict(json.loads(text))


def _model(design: object) -> Any:
    """Return a private ``scadnano.Design`` built from a JSON mapping."""

    data = _design_json(design)
    sc = _sc()
    with _guard("load the design"):
        return sc.Design.from_scadnano_json_map(data)


def _dump(model: Any) -> dict[str, Any]:
    with _guard("serialise the design"):
        return dict(json.loads(model.to_json()))


def _strand(model: Any, index: object) -> Any:
    position = _int(index, "strand_index", minimum=0)
    if position >= len(model.strands):
        raise ConfigurationError(
            f"strand_index {position} is out of range; the design has {len(model.strands)} strands."
        )
    return model.strands[position]


def _domain(strand: Any, index: object) -> Any:
    position = _int(index, "domain_index", minimum=0)
    if position >= len(strand.domains):
        raise ConfigurationError(
            f"domain_index {position} is out of range; "
            f"the strand has {len(strand.domains)} domains."
        )
    return strand.domains[position]


def _bound_domain(strand: Any, index: object) -> Any:
    sc = _sc()
    domain = _domain(strand, index)
    if not isinstance(domain, sc.Domain):
        raise ConfigurationError("domain_index must name a helix-bound domain, not a loopout.")
    return domain


def _helix(model: Any, index: object) -> Any:
    position = _int(index, "helix", minimum=0)
    if position not in model.helices:
        raise ConfigurationError(f"helix {position} does not exist in the design.")
    return model.helices[position]


def _domain_kind(domain: Any) -> str:
    sc = _sc()
    if isinstance(domain, sc.Loopout):
        return "loopout"
    if isinstance(domain, sc.Extension):
        return "extension"
    return "domain"


def _domain_info(domain: Any) -> dict[str, Any]:
    info = dict(json.loads(json.dumps(domain.to_json_serializable(suppress_indent=False))))
    info["kind"] = _domain_kind(domain)
    return info


def _strand_info(strand: Any, index: int) -> dict[str, Any]:
    return {
        "strand_index": index,
        "name": strand.name,
        "scaffold": bool(strand.is_scaffold),
        "circular": bool(strand.circular),
        "domains": len(strand.domains),
        "dna_length": strand.dna_length(),
        "dna_sequence": strand.dna_sequence,
    }


def _summary(model: Any) -> dict[str, Any]:
    return {
        "grid": str(model.grid.value),
        "helices": len(model.helices),
        "strands": len(model.strands),
        "scaffolds": sum(1 for strand in model.strands if strand.is_scaffold),
        "strand_names": [strand.name for strand in model.strands][:200],
        "strand_lengths": [strand.dna_length() for strand in model.strands][:_MAX_ITEMS],
    }


def _edit_result(model: Any, **extra: Any) -> dict[str, Any]:
    return {"design": _dump(model), **_summary(model), **extra, "provenance": _provenance()}


def _domain_from_json(sc: Any, mapping: object) -> Any:
    if not isinstance(mapping, Mapping):
        raise ConfigurationError("a domain must be a scadnano JSON mapping.")
    data = dict(mapping)
    with _guard("read the domain"):
        if "loopout" in data:
            return sc.Loopout.from_json(data)
        if "extension_num_bases" in data:
            return sc.Extension.from_json(data)
        return sc.Domain.from_json(data)


def _plate_type(sc: Any, value: object) -> Any:
    _choice(value, _PLATE_TYPES, "plate_type")
    return sc.PlateType.wells96 if value == 96 else sc.PlateType.wells384


# ---------------------------------------------------------------------------------------------
# design input and output
# ---------------------------------------------------------------------------------------------


def load_scadnano_design(
    design: Mapping[str, Any] | None = None,
    *,
    path: str | Path | None = None,
    json_text: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load and validate a scadnano design from a ``.sc`` file, a JSON string or a JSON mapping.

    Covers ``Design.from_scadnano_file`` (``path``), ``Design.from_scadnano_json_str``
    (``json_text``) and ``Design.from_scadnano_json_map`` (``design``); give exactly one of them.
    Returns the normalised design as a JSON mapping plus grid, helix and strand counts.
    """

    given = [value is not None for value in (design, path, json_text)]
    if sum(given) != 1:
        raise ConfigurationError("give exactly one of design, path or json_text.")
    sc = _sc()
    with _guard("load the design"):
        if path is not None:
            source = Path(path)
            if not source.is_file():
                raise ConfigurationError(f"design file does not exist: {source}")
            if source.stat().st_size > _MAX_TEXT:
                raise ConfigurationError("design file is larger than 50 MB.")
            model = sc.Design.from_scadnano_file(str(source))
            origin = "file"
        elif json_text is not None:
            if len(_text(json_text, "json_text")) > _MAX_TEXT:
                raise ConfigurationError("json_text is larger than 50 MB.")
            model = sc.Design.from_scadnano_json_str(json_text)
            origin = "json_string"
        else:
            model = _model(design)
            origin = "json_map"
    return _edit_result(model, source_kind=origin)


def export_scadnano_design(
    design: Mapping[str, Any],
    *,
    output_path: str | Path | None = None,
    suppress_indent: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Serialise a design to scadnano ``.sc`` JSON text and optionally write it to a new file.

    Covers ``Design.to_json`` and ``Design.write_scadnano_file``; ``suppress_indent=True`` is
    scadnano's compact default layout.
    """

    model = _model(design)
    with _guard("export the design"):
        text = model.to_json(_flag(suppress_indent, "suppress_indent"))
    return {
        "json_text": text,
        "output_path": _write(output_path, text),
        **_summary(model),
        "provenance": _provenance(),
    }


def import_cadnano_v2_design(
    source: str | Path | Mapping[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Convert a cadnano v2 design (a ``.json`` file or its parsed mapping) to a scadnano design.

    Covers ``Design.from_cadnano_v2``.
    """

    sc = _sc()
    with _guard("import the cadnano v2 design"):
        if isinstance(source, Mapping):
            model = sc.Design.from_cadnano_v2(json_dict=copy.deepcopy(dict(source)))
        else:
            file = Path(source)
            if not file.is_file():
                raise ConfigurationError(f"cadnano v2 file does not exist: {file}")
            if file.stat().st_size > _MAX_TEXT:
                raise ConfigurationError("cadnano v2 file is larger than 50 MB.")
            model = sc.Design.from_cadnano_v2(filename=str(file))
    return _edit_result(model)


def export_cadnano_v2_design(
    design: Mapping[str, Any],
    *,
    name: str = "",
    whitespace: bool = True,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a design as cadnano v2 JSON (``to_cadnano_v2_json``, ``write_cadnano_v2_file``).

    Returns the JSON text and the same content as a mapping; optionally writes a new file.
    """

    model = _model(design)
    with _guard("export the cadnano v2 design"):
        text = model.to_cadnano_v2_json(
            name=_text(name, "name", empty=True), whitespace=_flag(whitespace, "whitespace")
        )
    return {
        "json_text": text,
        "cadnano_v2": json.loads(text),
        "output_path": _write(output_path, text),
        "provenance": _provenance(),
    }


def export_oxview_design(
    design: Mapping[str, Any],
    *,
    use_strand_colors: bool = True,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a design as oxView JSON.

    Covers ``to_oxview_format``, ``to_oxview_json`` and ``write_oxview_file``.

    Every strand needs a complete DNA sequence only if you want sequences in oxView; scadnano
    fills unassigned bases itself. Returns the JSON text, strand and nucleotide counts, and
    optionally writes a new file.
    """

    model = _model(design)
    with _guard("export the oxView design"):
        text = model.to_oxview_format(
            warn_duplicate_strand_names=False,
            use_strand_colors=_flag(use_strand_colors, "use_strand_colors"),
        )
    parsed = json.loads(text)
    systems = parsed.get("systems", [])
    strands = [s for system in systems for s in system.get("strands", [])]
    nucleotides = sum(len(s.get("monomers", [])) for s in strands)
    return {
        "json_text": text,
        "strands": len(strands),
        "nucleotides": nucleotides,
        "output_path": _write(output_path, text),
        "provenance": _provenance(),
    }


def export_idt_order(
    design: Mapping[str, Any],
    *,
    order_format: str = "bulk_text",
    output_path: str | Path | None = None,
    delimiter: str = ",",
    domain_delimiter: str = "",
    only_strands_with_vendor_fields: bool = False,
    export_scaffold: bool = False,
    export_non_modified_strand_version: bool = False,
    plate_type: int = 96,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export the strands of a design as an IDT order.

    ``order_format="bulk_text"`` covers ``Design.to_idt_bulk_input_format`` and
    ``write_idt_bulk_input_file`` (name, sequence with vendor modification codes, scale,
    purification, one line per strand). ``order_format="plate_excel"`` covers
    ``write_idt_plate_excel_file`` and needs ``output_path`` (a new ``.xlsx`` file); it uses
    ``plate_type`` 96 or 384. The scaffold is left out unless ``export_scaffold=True``.
    """

    _choice(order_format, ("bulk_text", "plate_excel"), "order_format")
    sc = _sc()
    model = _model(design)
    options = {
        "only_strands_with_vendor_fields": _flag(
            only_strands_with_vendor_fields, "only_strands_with_vendor_fields"
        ),
        "export_scaffold": _flag(export_scaffold, "export_scaffold"),
        "export_non_modified_strand_version": _flag(
            export_non_modified_strand_version, "export_non_modified_strand_version"
        ),
    }
    if order_format == "bulk_text":
        with _guard("export the IDT bulk input"):
            text = model.to_idt_bulk_input_format(
                delimiter=_text(delimiter, "delimiter"),
                domain_delimiter=_text(domain_delimiter, "domain_delimiter", empty=True),
                warn_duplicate_name=False,
                **options,
            )
        return {
            "order_format": order_format,
            "text": text,
            "lines": len(text.splitlines()),
            "output_path": _write(output_path, text),
            "provenance": _provenance(),
        }
    if output_path is None:
        raise ConfigurationError("order_format='plate_excel' needs output_path.")
    target = _new_output(output_path)
    if target.suffix.lower() != ".xlsx":
        raise ConfigurationError("output_path must end in .xlsx for the plate workbook.")
    require("openpyxl", _EXTRA)
    with _guard("export the IDT plate workbook"):
        model.write_idt_plate_excel_file(
            directory=str(target.parent),
            filename=target.name,
            warn_duplicate_name=False,
            warn_using_default_plates=False,
            plate_type=_plate_type(sc, plate_type),
            **options,
        )
    if not target.is_file():
        raise ConfigurationError("scadnano did not write the plate workbook.")
    return {
        "order_format": order_format,
        "output_path": str(target),
        "bytes": target.stat().st_size,
        "plate_type": plate_type,
        "provenance": _provenance(),
    }


def build_plate_maps(
    design: Mapping[str, Any],
    *,
    plate_type: int = 96,
    strand_indices: Sequence[int] | None = None,
    table_format: str = "pipe",
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign strands to 96 or 384 well synthesis plates (``Design.plate_maps``).

    Only strands whose vendor fields name a plate and a well are placed (scadnano's rule); one
    plate map is returned per plate name. Each holds the well to strand-name map and the printable
    table (``PlateMap.to_table``; ``table_format`` is a tabulate format such as ``pipe``).
    """

    sc = _sc()
    model = _model(design)
    chosen = None
    if strand_indices is not None:
        chosen = [_strand(model, i) for i in strand_indices]
    with _guard("build plate maps"):
        maps = model.plate_maps(
            warn_duplicate_strand_names=False,
            plate_type=_plate_type(sc, plate_type),
            strands=chosen,
        )
        plates = [
            {
                "plate_name": plate.plate_name,
                "wells": {
                    well: strand.name or strand.default_export_name()
                    for well, strand in plate.well_to_strand.items()
                },
                "table": plate.to_table(
                    tablefmt=_text(table_format, "table_format"),
                    warn_unsupported_title_format=False,
                ),
            }
            for plate in maps
        ]
    return {"plate_type": plate_type, "plates": plates, "provenance": _provenance()}


# ---------------------------------------------------------------------------------------------
# DNA sequences
# ---------------------------------------------------------------------------------------------


def get_strand_sequence(
    design: Mapping[str, Any],
    strand_index: int,
    *,
    domain_index: int | None = None,
    delimiter: str = "",
    vendor_format: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the DNA sequence of a strand, or of one of its domains, loopouts or extensions.

    ``vendor_format=False`` gives the plain sequence with ``delimiter`` between domains
    (``Strand.dna_sequence_delimited``). ``vendor_format=True`` gives the sequence a synthesis
    vendor expects, with modification codes inserted (``Strand.vendor_dna_sequence`` with
    ``delimiter`` as ``domain_delimiter``; with ``domain_index`` the ``vendor_dna_sequence`` of
    that domain, loopout or extension). Unassigned bases are ``?``.
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    text = _text(delimiter, "delimiter", empty=True)
    with _guard("read the sequence"):
        if domain_index is not None:
            domain = _domain(strand, domain_index)
            sequence = domain.vendor_dna_sequence() if vendor_format else domain.dna_sequence
        elif _flag(vendor_format, "vendor_format"):
            sequence = strand.vendor_dna_sequence(domain_delimiter=text)
        else:
            sequence = strand.dna_sequence_delimited(text) if strand.dna_sequence else None
    return {
        "strand_index": strand_index,
        "domain_index": domain_index,
        "sequence": sequence,
        "vendor_format": vendor_format,
        "provenance": _provenance(),
    }


def get_domain_subsequence(
    design: Mapping[str, Any],
    strand_index: int,
    domain_index: int,
    offset_left: int,
    offset_right: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """DNA subsequence of a domain between two helix offsets (``Domain.dna_sequence_in``).

    Both offsets are helix offsets and BOTH ARE INCLUSIVE (scadnano's convention for this call);
    the result is written 5' to 3' on the strand and accounts for deletions and insertions.
    """

    model = _model(design)
    domain = _bound_domain(_strand(model, strand_index), domain_index)
    left = _int(offset_left, "offset_left")
    right = _int(offset_right, "offset_right")
    with _guard("read the domain subsequence"):
        sequence = domain.dna_sequence_in(left, right)
    return {
        "sequence": sequence,
        "offset_left": left,
        "offset_right": right,
        "provenance": _provenance(),
    }


def assign_dna_sequence(
    design: Mapping[str, Any],
    sequence: str,
    *,
    strand_index: int,
    domain_index: int | None = None,
    assign_complement: bool = True,
    check_length: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign DNA to a strand, or to one domain of a strand (``Design.assign_dna``).

    Covers ``Design.assign_dna``, ``Strand.set_dna_sequence`` and the builder's ``with_sequence``
    and ``with_domain_sequence``. With ``assign_complement=True`` the Watson-Crick complement is
    written to the strands bound to the assigned region. ``?`` is the unknown-base wildcard.
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    domain = None if domain_index is None else _domain(strand, domain_index)
    with _guard("assign the sequence"):
        model.assign_dna(
            strand,
            _sequence_text(sequence, "sequence"),
            assign_complement=_flag(assign_complement, "assign_complement"),
            domain=domain,
            check_length=_flag(check_length, "check_length"),
        )
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def assign_complement_sequence(
    design: Mapping[str, Any],
    strand_index: int,
    other_strand_index: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Give a strand the Watson-Crick complement of the strand it is bound to.

    Covers ``Strand.assign_dna_complement_from``: the DNA of ``other_strand_index`` is copied as
    its complement onto the bound regions of ``strand_index``; unbound bases stay unassigned.
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    other = _strand(model, other_strand_index)
    with _guard("assign the complement"):
        strand.assign_dna_complement_from(other)
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def assign_helix_domain_sequences(
    design: Mapping[str, Any],
    strand_index: int,
    *,
    helix: int,
    forward: bool,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign a strand the DNA of the domains on one helix (``Design.assign_dna_sequence``).

    ``forward`` is the direction of the strand being assigned. As in scadnano, the sequences of
    the domains on ``helix`` running in the OPPOSITE direction are concatenated in helix order
    and written to ``strand_index`` as they are (no complement is taken, no complement is
    written to other strands).
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    _helix(model, helix)
    with _guard("assign the sequence from the helix"):
        model.assign_dna_sequence(strand, _flag(forward, "forward"), helix)
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def check_complementarity(
    first: str,
    second: str,
    *,
    mode: str = "bases",
    allow_wildcard: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test Watson-Crick complementarity of two bases or of two sequences.

    ``mode="bases"`` covers ``scadnano.bases_complementary`` (one base each);
    ``mode="sequences"`` covers ``scadnano.reverse_complementary``, which compares the two
    sequences base by base from the first position (so ``second`` is the complement of ``first``
    written in the same order, scadnano's definition). ``allow_wildcard=True`` lets ``?`` match
    anything.
    """

    _choice(mode, ("bases", "sequences"), "mode")
    sc = _sc()
    a = _sequence_text(first, "first")
    b = _sequence_text(second, "second")
    wild = _flag(allow_wildcard, "allow_wildcard")
    with _guard("compare the sequences"):
        if mode == "bases":
            result = sc.bases_complementary(a, b, allow_wildcard=wild)
        else:
            result = sc.reverse_complementary(a, b, allow_wildcard=wild)
    return {"complementary": bool(result), "mode": mode, "provenance": _provenance()}


def get_m13_sequence(
    *, rotation: int = 5587, variant: str = "p7249", verbose: bool = True
) -> dict[str, Any]:
    """Return the M13mp18 scaffold sequence, rotated (``scadnano.m13``).

    ``variant`` is ``p7249`` (default, 7249 nt), ``p7560``, ``p8064`` or ``p8634``; ``rotation`` is
    scadnano's default 5587, the cadnano convention for the 5' end.
    """

    _choice(variant, _M13_VARIANTS, "variant")
    sc = _sc()
    spin = _int(rotation, "rotation", minimum=0)
    with _guard("read the M13 sequence"):
        sequence = sc.m13(rotation=spin, variant=sc.M13Variant[variant])
    return {
        "sequence": sequence,
        "length": len(sequence),
        "variant": variant,
        "rotation": spin,
        "provenance": _provenance(),
    }


def assign_m13_scaffold(
    design: Mapping[str, Any],
    *,
    rotation: int = 5587,
    variant: str = "p7249",
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign M13mp18 DNA to the scaffold strand and complements to the staples
    (``Design.assign_m13_to_scaffold``).
    """

    _choice(variant, _M13_VARIANTS, "variant")
    sc = _sc()
    model = _model(design)
    spin = _int(rotation, "rotation", minimum=0)
    with _guard("assign M13 to the scaffold"):
        model.assign_m13_to_scaffold(rotation=spin, variant=sc.M13Variant[variant])
    return _edit_result(model, variant=variant, rotation=spin)


def get_strand_without_modifications(
    design: Mapping[str, Any], strand_index: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Return a copy of a strand with all 5', 3' and internal modifications removed
    (``Strand.no_modifications_version``); the design itself is not changed.
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    with _guard("copy the strand"):
        plain = strand.no_modifications_version()
        info = dict(json.loads(json.dumps(plain.to_json_serializable(suppress_indent=False))))
        sequence = plain.dna_sequence
    return {"strand": info, "dna_sequence": sequence, "provenance": _provenance()}


def inline_deletions_insertions(
    design: Mapping[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Fold deletions and insertions into domain lengths (``Design.inline_deletions_insertions``).

    The design afterwards has no deletion or insertion entries; domain starts and ends shift so
    that the total number of bases stays the same.
    """

    model = _model(design)
    before = [strand.dna_length() for strand in model.strands]
    with _guard("inline deletions and insertions"):
        model.inline_deletions_insertions()
    return _edit_result(model, strand_lengths_before=before)


# ---------------------------------------------------------------------------------------------
# modifications and base pairs
# ---------------------------------------------------------------------------------------------


def set_strand_modification(
    design: Mapping[str, Any],
    strand_index: int,
    *,
    location: str,
    display_text: str,
    vendor_code: str,
    connector_length: int = 4,
    index: int | None = None,
    allowed_bases: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attach a 5', 3' or internal DNA modification to a strand.

    ``location`` is ``5p`` (``Strand.set_modification_5p``), ``3p`` (``set_modification_3p``) or
    ``internal`` (``set_modification_internal``; ``index`` is the 0-based base after which, or on
    which, the modification sits, and ``allowed_bases`` restricts the bases it may attach to).
    ``vendor_code`` is the string the vendor uses, for example ``/5Biosg/``.
    """

    _choice(location, _MOD_LOCATIONS, "location")
    sc = _sc()
    model = _model(design)
    strand = _strand(model, strand_index)
    text = _text(display_text, "display_text")
    code = _text(vendor_code, "vendor_code")
    length = _int(connector_length, "connector_length", minimum=0)
    with _guard("set the modification"):
        if location == "5p":
            strand.set_modification_5p(sc.Modification5Prime(text, code, length))
        elif location == "3p":
            strand.set_modification_3p(sc.Modification3Prime(text, code, length))
        else:
            if index is None:
                raise ConfigurationError("an internal modification needs index.")
            bases = None if allowed_bases is None else frozenset(allowed_bases)
            modification = sc.ModificationInternal(text, code, length, allowed_bases=bases)
            strand.set_modification_internal(
                _int(index, "index", minimum=0), modification, warn_on_no_dna=False
            )
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def remove_strand_modification(
    design: Mapping[str, Any],
    strand_index: int,
    *,
    location: str,
    index: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the 5', 3' or an internal modification from a strand
    (``remove_modification_5p``, ``remove_modification_3p``, ``remove_modification_internal``).
    """

    _choice(location, _MOD_LOCATIONS, "location")
    model = _model(design)
    strand = _strand(model, strand_index)
    with _guard("remove the modification"):
        if location == "5p":
            strand.remove_modification_5p()
        elif location == "3p":
            strand.remove_modification_3p()
        else:
            if index is None:
                raise ConfigurationError("an internal modification needs index.")
            strand.remove_modification_internal(_int(index, "index", minimum=0))
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def list_design_modifications(
    design: Mapping[str, Any], *, location: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List the distinct DNA modifications used in a design (``Design.modifications``).

    ``location`` limits the list to ``5p``, ``3p`` or ``internal`` modifications.
    """

    sc = _sc()
    model = _model(design)
    kind = None
    if location is not None:
        _choice(location, _MOD_LOCATIONS, "location")
        kind = {
            "5p": sc.ModificationType.five_prime,
            "3p": sc.ModificationType.three_prime,
            "internal": sc.ModificationType.internal,
        }[location]
    with _guard("list the modifications"):
        found = model.modifications(kind) if kind is not None else model.modifications()
        rows = [
            dict(json.loads(json.dumps(m.to_json_serializable(suppress_indent=False))))
            for m in found
        ]
    return {"modifications": rows, "count": len(rows), "provenance": _provenance()}


def compute_base_pairs(
    design: Mapping[str, Any], *, allow_mismatches: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Compute which bases of the design are paired (``Design.base_pairs``).

    Returns a map from helix index to the list of paired offsets (`forward` and reverse strand at
    the same offset), counting offsets in which both strands carry a base. With
    ``allow_mismatches=False`` only Watson-Crick pairs count and an error is raised for strands
    without a complete sequence, as in scadnano.
    """

    model = _model(design)
    with _guard("compute base pairs"):
        pairs = model.base_pairs(allow_mismatches=_flag(allow_mismatches, "allow_mismatches"))
    table = {str(helix): sorted(int(o) for o in offsets) for helix, offsets in pairs.items()}
    return {
        "base_pairs": table,
        "total": sum(len(v) for v in table.values()),
        "provenance": _provenance(),
    }


# ---------------------------------------------------------------------------------------------
# helices, strands and domains
# ---------------------------------------------------------------------------------------------


def _locator(model: Any) -> dict[int, tuple[int, int]]:
    return {
        id(domain): (si, di)
        for si, strand in enumerate(model.strands)
        for di, domain in enumerate(strand.domains)
    }


def _located(model: Any, domains: Sequence[Any]) -> list[dict[str, Any]]:
    where = _locator(model)
    rows = []
    for domain in domains:
        si, di = where.get(id(domain), (None, None))
        row = {"strand_index": si, "domain_index": di}
        row.update(_domain_info(domain))
        rows.append(row)
    return rows


def add_design_helix(
    design: Mapping[str, Any],
    *,
    idx: int,
    max_offset: int,
    min_offset: int = 0,
    grid_position: Sequence[int] | None = None,
    roll: float = 0.0,
    group: str = "default_group",
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a helix with index ``idx`` to a design (``Design.add_helix``).

    ``max_offset`` is the exclusive upper offset, so the helix holds ``max_offset - min_offset``
    base pairs. ``grid_position`` is (h, v) on the design's grid; ``roll`` is in degrees.
    """

    sc = _sc()
    model = _model(design)
    position = None
    if grid_position is not None:
        if len(grid_position) not in (2, 3) or any(
            isinstance(v, bool) or not isinstance(v, int) for v in grid_position
        ):
            raise ConfigurationError("grid_position must be two or three integers.")
        position = tuple(int(v) for v in grid_position)
    if isinstance(roll, bool) or not isinstance(roll, (int, float)):
        raise ConfigurationError("roll must be a number.")
    with _guard("add the helix"):
        helix = sc.Helix(
            grid_position=position,
            max_offset=_int(max_offset, "max_offset"),
            min_offset=_int(min_offset, "min_offset"),
            roll=float(roll),
            group=_text(group, "group"),
        )
        position_idx = _int(idx, "idx", minimum=0)
        model.add_helix(position_idx, helix)
        # scadnano's add_helix leaves helix.idx at -1, which makes the exported JSON unreadable.
        model.helices[position_idx].idx = position_idx
    return _edit_result(model)


def add_design_strand(
    design: Mapping[str, Any], strand: Mapping[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Add a strand, given as its scadnano JSON mapping, to a design (``Design.add_strand``).

    The mapping has the same shape as an entry of ``design["strands"]``: ``domains`` (each a
    helix domain, ``{"loopout": n}`` or ``{"extension_num_bases": n, ...}``), and optionally
    ``sequence``, ``color``, ``is_scaffold``, ``circular``, ``name``, ``vendor`` and
    modification entries. The strand is appended to the end of the strand list.
    """

    if not isinstance(strand, Mapping) or not isinstance(strand.get("domains"), Sequence):
        raise ConfigurationError("strand must be a scadnano strand mapping with a domains list.")
    data = _design_json(design)
    strands = list(data.get("strands", []))
    strands.append(json.loads(json.dumps(strand)))
    data["strands"] = strands
    model = _model(data)
    return _edit_result(model, added_strand_index=len(model.strands) - 1)


def build_strand(
    design: Mapping[str, Any],
    *,
    helix: int,
    offset: int,
    steps: Sequence[Mapping[str, Any]],
    verbose: bool = True,
) -> dict[str, Any]:
    """Build one strand domain by domain and add it to the design (``Design.draw_strand``).

    This is scadnano's fluent builder (``StrandBuilder``) with the first domain starting at
    ``helix`` and ``offset``. Each step is a mapping with ``"op"`` and its arguments:
    ``to`` (offset), ``move`` (delta), ``cross`` (helix, offset or move: crossover to another
    helix, covering ``StrandBuilder.cross``), ``loopout`` (helix, length, offset or move),
    ``extension_3p`` / ``extension_5p`` (num_bases, display_length, display_angle),
    ``with_deletions`` (deletions), ``with_insertions`` (list of [offset, length]),
    ``with_sequence`` / ``with_domain_sequence`` (sequence), ``with_name``, ``with_label``,
    ``with_domain_name``, ``as_scaffold``, ``as_circular`` and ``with_vendor_fields``
    (scale, purification, plate, well). Example: ``[{"op": "to", "offset": 16},
    {"op": "cross", "helix": 1}, {"op": "move", "delta": -16}]``.
    """

    model = _model(design)
    _helix(model, helix)
    if not isinstance(steps, Sequence) or not steps:
        raise ConfigurationError("steps must be a non-empty list of mappings.")
    if len(steps) > 10_000:
        raise ConfigurationError("steps may hold at most 10000 entries.")
    before = len(model.strands)
    with _guard("build the strand"):
        builder = model.draw_strand(helix, _int(offset, "offset"))
        for position, step in enumerate(steps):
            if not isinstance(step, Mapping) or step.get("op") not in _BUILDER_STEPS:
                raise ConfigurationError(
                    f"step {position} needs an 'op' from {sorted(_BUILDER_STEPS)}."
                )
            op = step["op"]
            extra = set(step) - {"op"} - set(_BUILDER_STEPS[op])
            if extra:
                raise ConfigurationError(
                    f"step {position} ({op}) has unknown keys {sorted(extra)}."
                )
            arguments = {k: v for k, v in step.items() if k != "op"}
            if op == "with_insertions":
                arguments["insertions"] = [tuple(pair) for pair in arguments["insertions"]]
            builder = getattr(builder, op)(**arguments)
    if len(model.strands) != before + 1:
        raise ConfigurationError("the steps did not produce a strand; add at least one domain.")
    index = len(model.strands) - 1
    return _edit_result(
        model, added_strand_index=index, strand=_strand_info(model.strands[index], index)
    )


def remove_design_strand(
    design: Mapping[str, Any], strand_index: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Remove a strand from a design (``Design.remove_strand``); later strands move down one."""

    model = _model(design)
    strand = _strand(model, strand_index)
    with _guard("remove the strand"):
        model.remove_strand(strand)
    return _edit_result(model, removed_strand_index=strand_index)


def find_strand_by_name(
    design: Mapping[str, Any], name: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Find the strand with a given name (``Design.strand_with_name``).

    ``found`` is False if no strand has that name.
    """

    model = _model(design)
    wanted = _text(name, "name")
    for index, strand in enumerate(model.strands):
        if strand.name == wanted:
            return {"found": True, **_strand_info(strand, index), "provenance": _provenance()}
    return {"found": False, "name": wanted, "provenance": _provenance()}


def list_strands_at_helix_end(
    design: Mapping[str, Any], helix: int, *, end: str = "5p", verbose: bool = True
) -> dict[str, Any]:
    """List the strands whose 5' end (``end="5p"``, ``Design.strands_starting_on_helix``) or 3' end
    (``end="3p"``, ``Design.strands_ending_on_helix``) lies on a helix.
    """

    _choice(end, ("5p", "3p"), "end")
    model = _model(design)
    _helix(model, helix)
    with _guard("list strands"):
        strands = (
            model.strands_starting_on_helix(helix)
            if end == "5p"
            else model.strands_ending_on_helix(helix)
        )
    where = {id(s): i for i, s in enumerate(model.strands)}
    rows = [_strand_info(s, where[id(s)]) for s in strands]
    return {"helix": helix, "end": end, "strands": rows, "provenance": _provenance()}


def find_domains_at_offset(
    design: Mapping[str, Any],
    helix: int,
    offset: int,
    *,
    forward: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the domains that cover a helix offset (``Design.domains_at``).

    With ``forward`` True or False only the domain on that strand direction is returned
    (``Design.domain_at``; ``None`` is returned if there is none).
    """

    model = _model(design)
    _helix(model, helix)
    position = _int(offset, "offset")
    with _guard("find domains"):
        if forward is None:
            found = list(model.domains_at(helix, position))
        else:
            one = model.domain_at(helix, position, _flag(forward, "forward"))
            found = [] if one is None else [one]
    return {"domains": _located(model, found), "provenance": _provenance()}


def insert_strand_domain(
    design: Mapping[str, Any],
    strand_index: int,
    domain: Mapping[str, Any],
    *,
    order: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert a domain, loopout or extension into a strand.

    ``domain`` is a scadnano domain mapping, for example ``{"helix": 1, "forward": false,
    "start": 0, "end": 8}``. ``order`` is the position in the strand (``Strand.insert_domain``,
    ``Design.insert_domain``); ``None`` appends it at the 3' end (``Design.append_domain``).
    """

    sc = _sc()
    model = _model(design)
    strand = _strand(model, strand_index)
    new = _domain_from_json(sc, domain)
    with _guard("insert the domain"):
        if order is None:
            model.append_domain(strand, new)
        else:
            model.insert_domain(strand, _int(order, "order", minimum=0), new)
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def remove_strand_domain(
    design: Mapping[str, Any], strand_index: int, domain_index: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Remove a domain from a strand (``Design.remove_domain``, ``Strand.remove_domain``)."""

    model = _model(design)
    strand = _strand(model, strand_index)
    domain = _domain(strand, domain_index)
    with _guard("remove the domain"):
        model.remove_domain(strand, domain)
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def set_strand_domains(
    design: Mapping[str, Any],
    strand_index: int,
    domains: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the whole domain list of a strand (``Strand.set_domains``).

    ``domains`` are scadnano domain mappings.
    """

    sc = _sc()
    model = _model(design)
    strand = _strand(model, strand_index)
    if not isinstance(domains, Sequence) or not domains:
        raise ConfigurationError("domains must be a non-empty list of domain mappings.")
    items = [_domain_from_json(sc, d) for d in domains]
    with _guard("set the domains"):
        strand.set_domains(items)
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def rotate_strand_domains(
    design: Mapping[str, Any],
    strand_index: int,
    rotation: int,
    *,
    forward: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rotate the domain order of a strand by ``rotation`` places (``Strand.rotate_domains``)."""

    model = _model(design)
    strand = _strand(model, strand_index)
    with _guard("rotate the domains"):
        strand.rotate_domains(_int(rotation, "rotation"), forward=_flag(forward, "forward"))
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


def reverse_strands(
    design: Mapping[str, Any], *, strand_index: int | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Reverse the polarity of one strand (``Strand.reverse``) or of every strand
    (``Design.reverse_all``, ``strand_index=None``): domain order and directions are flipped.
    """

    model = _model(design)
    with _guard("reverse the strands"):
        if strand_index is None:
            model.reverse_all()
        else:
            _strand(model, strand_index).reverse()
    return _edit_result(model)


def set_strand_properties(
    design: Mapping[str, Any],
    strand_index: int,
    *,
    scaffold: bool | None = None,
    circular: bool | None = None,
    vendor_fields: Mapping[str, str | None] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Change strand properties.

    ``scaffold`` marks the strand as scaffold or not (``Strand.set_scaffold``). ``circular`` makes
    it circular or linear (``Strand.set_circular``, ``Strand.set_linear``). ``vendor_fields`` is a
    mapping with ``scale``, ``purification``, ``plate`` and ``well`` (``VendorFields``). Leave an
    argument as ``None`` to keep it.
    """

    sc = _sc()
    model = _model(design)
    strand = _strand(model, strand_index)
    with _guard("set the strand properties"):
        if scaffold is not None:
            strand.set_scaffold(_flag(scaffold, "scaffold"))
        if circular is not None:
            if _flag(circular, "circular"):
                strand.set_circular()
            else:
                strand.set_linear()
        if vendor_fields is not None:
            allowed = {"scale", "purification", "plate", "well"}
            if not isinstance(vendor_fields, Mapping) or set(vendor_fields) - allowed:
                raise ConfigurationError(f"vendor_fields keys must be among {sorted(allowed)}.")
            strand.vendor_fields = sc.VendorFields(**dict(vendor_fields))
    return _edit_result(model, assigned=_strand_info(strand, strand_index))


# ---------------------------------------------------------------------------------------------
# overlaps, crossovers and indices
# ---------------------------------------------------------------------------------------------


def check_domain_overlap(
    design: Mapping[str, Any],
    strand_index: int,
    other_strand_index: int,
    *,
    domain_index: int | None = None,
    other_domain_index: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether two domains, or two whole strands, overlap on a helix.

    With both domain indices: ``overlaps`` is ``Domain.overlaps`` (same helix, overlapping offsets
    and OPPOSITE directions, so the domains can pair); ``overlap_interval`` is
    ``Domain.compute_overlap`` (start inclusive, end exclusive; ``None`` when the offsets do not
    overlap on the same helix); ``overlaps_illegally`` is ``Domain.overlaps_illegally`` (overlap
    in the SAME direction, which scadnano forbids). Without them: ``Strand.overlaps``.
    """

    model = _model(design)
    first = _strand(model, strand_index)
    second = _strand(model, other_strand_index)
    with _guard("test the overlap"):
        if (domain_index is None) != (other_domain_index is None):
            raise ConfigurationError("give both domain indices or neither.")
        if domain_index is None:
            return {"overlaps": bool(first.overlaps(second)), "provenance": _provenance()}
        a = _bound_domain(first, domain_index)
        b = _bound_domain(second, other_domain_index)
        overlaps = bool(a.overlaps(b))
        low, high = a.compute_overlap(b)
        interval = [int(low), int(high)] if low >= 0 else None
        illegal = bool(a.overlaps_illegally(b))
    return {
        "overlaps": overlaps,
        "overlap_interval": interval,
        "overlaps_illegally": illegal,
        "provenance": _provenance(),
    }


def find_overlapping_domains(
    design: Mapping[str, Any], strand_index: int, domain_index: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Find the domains whose offset range overlaps the range of one domain.

    Covers ``Domain.find_overlapping_ranges``, which compares offset ranges only (cadnano's
    ``_findOverlappingRanges``); it does not look at the helix or the direction, so filter the
    returned rows by ``helix`` and ``forward`` if needed. Loopouts are ignored and the domain
    itself is not listed.
    """

    model = _model(design)
    domain = _bound_domain(_strand(model, strand_index), domain_index)
    others = [
        d
        for s in model.strands
        for d in s.domains
        if isinstance(d, _sc().Domain) and d is not domain
    ]
    with _guard("find overlapping domains"):
        found = list(domain.find_overlapping_ranges(others))
    return {"domains": _located(model, found), "provenance": _provenance()}


def find_helix_domain_overlap_pairs(
    design: Mapping[str, Any], *, helix: int | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List the pairs of overlapping domains of opposite direction on a helix
    (``scadnano.find_overlapping_domains_on_helix``).

    Each pair is a forward and a reverse domain whose offset ranges overlap on the helix, which
    is how scadnano finds the places where bases can pair. Two domains of the same direction that
    overlap are illegal in scadnano; that test is ``check_domain_overlap`` (``overlaps_illegally``)
    and loading such a design fails. ``helix=None`` checks every helix.
    """

    sc = _sc()
    model = _model(design)
    indices = sorted(model.helices) if helix is None else [helix]
    where = _locator(model)
    pairs = []
    with _guard("check overlaps"):
        for idx in indices:
            for a, b in sc.find_overlapping_domains_on_helix(_helix(model, idx)):
                pairs.append(
                    {
                        "helix": idx,
                        "first": dict(
                            zip(("strand_index", "domain_index"), where[id(a)], strict=True)
                        ),
                        "second": dict(
                            zip(("strand_index", "domain_index"), where[id(b)], strict=True)
                        ),
                    }
                )
    return {"pairs": pairs, "count": len(pairs), "provenance": _provenance()}


def domain_has_crossover(
    design: Mapping[str, Any],
    strand_index: int,
    domain_index: int,
    offset: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether a domain has a crossover at a helix offset (``Domain.has_crossover_at``)."""

    model = _model(design)
    domain = _bound_domain(_strand(model, strand_index), domain_index)
    with _guard("test the crossover"):
        result = domain.has_crossover_at(_int(offset, "offset"))
    return {"has_crossover": bool(result), "provenance": _provenance()}


def map_domain_offset_to_dna_index(
    design: Mapping[str, Any],
    strand_index: int,
    domain_index: int,
    *,
    offset: int | None = None,
    offset_closer_to_5p: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map positions of a domain to indices in the strand's DNA sequence.

    Always returns ``start_index``: the index in the strand DNA where the domain begins
    (``Domain.get_seq_start_idx``, ``Loopout.get_seq_start_idx``,
    ``Strand.dna_index_start_domain``).
    With ``offset`` it also returns ``dna_index``, the strand index of that helix offset
    (``Domain.domain_offset_to_strand_dna_idx``); ``offset_closer_to_5p`` picks the first or the
    last base when the offset holds an insertion.
    """

    model = _model(design)
    strand = _strand(model, strand_index)
    domain = _domain(strand, domain_index)
    with _guard("map the offset"):
        result: dict[str, Any] = {"start_index": int(strand.dna_index_start_domain(domain))}
        if offset is not None:
            bound = _bound_domain(strand, domain_index)
            result["dna_index"] = int(
                bound.domain_offset_to_strand_dna_idx(
                    _int(offset, "offset"), _flag(offset_closer_to_5p, "offset_closer_to_5p")
                )
            )
    return {**result, "provenance": _provenance()}


def list_helix_crossovers(
    design: Mapping[str, Any],
    helix: int,
    *,
    allow_intrahelix: bool = True,
    allow_intergroup: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the crossovers incident to a helix as (other helix, offset, forward) triples
    (``Helix.crossover_addresses``).
    """

    model = _model(design)
    one = _helix(model, helix)
    with _guard("list crossovers"):
        rows = one.crossover_addresses(
            model.helices,
            allow_intrahelix=_flag(allow_intrahelix, "allow_intrahelix"),
            allow_intergroup=_flag(allow_intergroup, "allow_intergroup"),
        )
    crossovers = [{"helix": int(h), "offset": int(o), "forward": bool(f)} for h, o, f in rows]
    return {"crossovers": crossovers, "count": len(crossovers), "provenance": _provenance()}


def split_helix_strands(
    design: Mapping[str, Any], helix: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Split the domains on a helix into staple and scaffold domains by cadnano's parity rule
    (``Helix.get_strand_sets``): on even helices the forward domains are scaffold.
    """

    model = _model(design)
    with _guard("split the helix strands"):
        staples, scaffolds = _helix(model, helix).get_strand_sets()
    return {
        "staple_domains": _located(model, staples),
        "scaffold_domains": _located(model, scaffolds),
        "provenance": _provenance(),
    }


def list_crossover_points(design: Mapping[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Return cadnano's crossover-point tables of the design grid (``Design.get_crossover_points``).

    Four lists per grid (square or honeycomb): the scaffold low and high offsets and the staple
    low and high offsets, in each lattice step.
    """

    model = _model(design)
    with _guard("get the crossover points"):
        scaf_low, scaf_high, stap_low, stap_high = model.get_crossover_points()
    return {
        "grid": str(model.grid.value),
        "scaffold_low": scaf_low,
        "scaffold_high": scaf_high,
        "staple_low": stap_low,
        "staple_high": stap_high,
        "provenance": _provenance(),
    }


def list_potential_crossovers(
    design: Mapping[str, Any], helix: int, *, index: int | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List the offsets where a crossover from a helix to a neighbour is possible
    (``Design.potential_crossover_list``): neighbour helix, offset, ``scaffold`` or ``staple``,
    and whether it is the low site of the crossover pair. ``index`` limits the search window.
    """

    model = _model(design)
    one = _helix(model, helix)
    with _guard("list potential crossovers"):
        rows = model.potential_crossover_list(
            one, None if index is None else _int(index, "index", minimum=0)
        )
    sites = [
        {"neighbor_helix": int(h.idx), "offset": int(i), "strand_type": t, "is_low": bool(low)}
        for h, i, t, low in (rows or [])
    ]
    return {"sites": sites, "count": len(sites), "provenance": _provenance()}


def get_helix_neighbors(
    design: Mapping[str, Any], helix: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Return the neighbouring helices of a helix on the lattice, ordered by parity
    (``Design.get_helix_neighbors``).
    """

    model = _model(design)
    one = _helix(model, helix)
    with _guard("get the neighbours"):
        neighbours = model.get_helix_neighbors(one)
    return {
        "neighbors": [None if h is None else int(h.idx) for h in neighbours],
        "provenance": _provenance(),
    }


# ---------------------------------------------------------------------------------------------
# editing operations on the design
# ---------------------------------------------------------------------------------------------


def add_crossover(
    design: Mapping[str, Any],
    *,
    helix: int,
    helix2: int,
    offset: int,
    forward: bool,
    offset2: int | None = None,
    forward2: bool | None = None,
    full: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a crossover between two helices.

    ``full=False`` adds a half crossover that joins two strand ends (``Design.add_half_crossover``).
    ``full=True`` adds a full crossover, two half crossovers one on each strand direction
    (``Design.add_full_crossover``). ``offset2`` and ``forward2`` default to ``offset`` and the
    opposite direction.
    """

    model = _model(design)
    _helix(model, helix)
    _helix(model, helix2)
    with _guard("add the crossover"):
        method = model.add_full_crossover if _flag(full, "full") else model.add_half_crossover
        method(
            helix,
            helix2,
            _int(offset, "offset"),
            _flag(forward, "forward"),
            None if offset2 is None else _int(offset2, "offset2"),
            None if forward2 is None else _flag(forward2, "forward2"),
        )
    return _edit_result(model)


def add_design_deletion(
    design: Mapping[str, Any], helix: int, offset: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Add a single-base deletion at a helix offset (``Design.add_deletion``)."""

    model = _model(design)
    _helix(model, helix)
    with _guard("add the deletion"):
        model.add_deletion(helix, _int(offset, "offset"))
    return _edit_result(model)


def add_design_insertion(
    design: Mapping[str, Any], helix: int, offset: int, length: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Add an insertion of ``length`` extra bases at a helix offset (``Design.add_insertion``).

    Strands that already carry a DNA sequence keep the old sequence, whose length no longer fits;
    assign the sequences after adding insertions.
    """

    model = _model(design)
    _helix(model, helix)
    with _guard("add the insertion"):
        model.add_insertion(helix, _int(offset, "offset"), _int(length, "length", minimum=1))
    return _edit_result(model)


def add_design_nick(
    design: Mapping[str, Any],
    helix: int,
    offset: int,
    forward: bool,
    *,
    new_color: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nick the domain at a helix offset, splitting its strand (``Design.add_nick``).

    The nick lies between ``offset - 1`` and ``offset``; the new strand gets a new colour when
    ``new_color`` is True.
    """

    model = _model(design)
    _helix(model, helix)
    with _guard("add the nick"):
        model.add_nick(
            helix,
            _int(offset, "offset"),
            _flag(forward, "forward"),
            new_color=_flag(new_color, "new_color"),
        )
    return _edit_result(model)


def ligate_design_nick(
    design: Mapping[str, Any], helix: int, offset: int, forward: bool, *, verbose: bool = True
) -> dict[str, Any]:
    """Join two strands that meet in a nick at a helix offset (``Design.ligate``)."""

    model = _model(design)
    _helix(model, helix)
    with _guard("ligate the nick"):
        model.ligate(helix, _int(offset, "offset"), _flag(forward, "forward"))
    return _edit_result(model)


def move_strands(
    design: Mapping[str, Any], delta: int, *, axis: str = "offset", verbose: bool = True
) -> dict[str, Any]:
    """Shift every strand of the design.

    ``axis="offset"`` moves all domains by ``delta`` offsets along their helices
    (``Design.move_strand_offsets``); ``axis="helix"`` moves them by ``delta`` helix indices
    (``Design.move_strands_on_helices``).
    """

    _choice(axis, ("offset", "helix"), "axis")
    model = _model(design)
    step = _int(delta, "delta")
    with _guard("move the strands"):
        if axis == "offset":
            model.move_strand_offsets(step)
        else:
            model.move_strands_on_helices(step)
    return _edit_result(model)


def autostaple_design(design: Mapping[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Generate staples over the scaffold with cadnano's autostaple (``Design.autostaple``).

    Existing staples are removed first; the designs must not hold loopouts or extensions.
    """

    model = _model(design)
    with _guard("autostaple the design"):
        model.autostaple()
    return _edit_result(model)


def relax_helix_rolls(
    design: Mapping[str, Any], *, compute_only: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Relax the roll angle of every helix so that crossovers face their neighbours.

    ``compute_only=False`` applies ``Design.relax_helix_rolls`` and returns the new design.
    ``compute_only=True`` only reports ``Helix.compute_relaxed_roll_delta`` for each helix and
    leaves the design unchanged.
    """

    model = _model(design)
    with _guard("relax the helix rolls"):
        if _flag(compute_only, "compute_only"):
            deltas = {
                str(idx): float(
                    h.compute_relaxed_roll_delta(model.helices, model.grid, model.geometry)
                )
                for idx, h in model.helices.items()
            }
            return {
                "design": _dump(model),
                "roll_deltas": deltas,
                "applied": False,
                "provenance": _provenance(),
            }
        model.relax_helix_rolls()
        rolls = {str(idx): float(h.roll) for idx, h in model.helices.items()}
    return _edit_result(model, rolls=rolls, applied=True)


def get_backbone_angle(
    design: Mapping[str, Any], helix: int, offset: int, forward: bool, *, verbose: bool = True
) -> dict[str, Any]:
    """Angle in degrees of the backbone at a helix offset (``Helix.backbone_angle_at_offset``).

    Uses the helix roll and the design geometry (``bases_per_turn``, ``minor_groove_angle``).
    """

    model = _model(design)
    one = _helix(model, helix)
    with _guard("compute the backbone angle"):
        angle = one.backbone_angle_at_offset(
            _int(offset, "offset"), _flag(forward, "forward"), model.geometry
        )
    return {"angle_degrees": float(angle), "provenance": _provenance()}


# ---------------------------------------------------------------------------------------------
# origami rectangle
# ---------------------------------------------------------------------------------------------

_NICK_PATTERNS = ("staggered", "staggered_opposite", "even", "odd")


def create_origami_rectangle(
    *,
    num_helices: int,
    num_cols: int,
    assign_seq: bool = True,
    seam_left_column: int = -1,
    nick_pattern: str = "staggered",
    twist_correction_deletion_spacing: int = 0,
    twist_correction_start_col: int = 1,
    twist_correction_deletion_offset: int = -1,
    num_flanking_columns: int = 1,
    num_flanking_helices: int = 0,
    custom_scaffold: str | None = None,
    edge_staples: bool = True,
    scaffold_nick_offset: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a DNA origami rectangle design (``origami_rectangle.create``).

    The scaffold runs through ``num_helices`` helices in ``num_cols`` columns of 16 bases with a
    two-column seam; staples are crossed over and nicked by ``nick_pattern`` (``staggered``,
    ``staggered_opposite``, ``even`` or ``odd``). With ``assign_seq=True`` M13mp18 (or
    ``custom_scaffold``) is assigned and staple sequences are derived. Twist-correction deletions
    are added every ``twist_correction_deletion_spacing`` columns when it is positive.
    """

    _sc()
    rectangle = require("scadnano.origami_rectangle", _EXTRA)
    _choice(nick_pattern, _NICK_PATTERNS, "nick_pattern")
    with _guard("create the origami rectangle"):
        model = rectangle.create(
            num_helices=_int(num_helices, "num_helices", minimum=1),
            num_cols=_int(num_cols, "num_cols", minimum=1),
            assign_seq=_flag(assign_seq, "assign_seq"),
            seam_left_column=_int(seam_left_column, "seam_left_column"),
            nick_pattern=rectangle.NickPattern[nick_pattern],
            twist_correction_deletion_spacing=_int(
                twist_correction_deletion_spacing, "twist_correction_deletion_spacing", minimum=0
            ),
            twist_correction_start_col=_int(
                twist_correction_start_col, "twist_correction_start_col"
            ),
            twist_correction_deletion_offset=_int(
                twist_correction_deletion_offset, "twist_correction_deletion_offset"
            ),
            num_flanking_columns=_int(num_flanking_columns, "num_flanking_columns", minimum=0),
            num_flanking_helices=_int(num_flanking_helices, "num_flanking_helices", minimum=0),
            custom_scaffold=None
            if custom_scaffold is None
            else _sequence_text(custom_scaffold, "custom_scaffold"),
            edge_staples=_flag(edge_staples, "edge_staples"),
            scaffold_nick_offset=_int(scaffold_nick_offset, "scaffold_nick_offset"),
        )
    return _edit_result(model)


def add_twist_correction_deletions(
    design: Mapping[str, Any],
    *,
    offset_start: int,
    deletion_spacing: int,
    deletion_start_col: int,
    deletion_offset: int,
    num_helices: int,
    num_cols: int,
    num_flanking_helices: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add twist-correction deletions to an origami rectangle
    (``origami_rectangle.add_twist_correction_deletions``): one deletion per helix every
    ``deletion_spacing`` columns from column ``deletion_start_col``. ``deletion_offset`` is the
    position inside the 16-base column; negative picks the offset furthest from the edges.
    """

    rectangle = require("scadnano.origami_rectangle", _EXTRA)
    model = _model(design)
    with _guard("add the twist-correction deletions"):
        rectangle.add_twist_correction_deletions(
            model,
            _int(offset_start, "offset_start"),
            _int(deletion_spacing, "deletion_spacing", minimum=1),
            _int(deletion_start_col, "deletion_start_col"),
            _int(deletion_offset, "deletion_offset"),
            _int(num_helices, "num_helices", minimum=1),
            _int(num_cols, "num_cols", minimum=1),
            _int(num_flanking_helices, "num_flanking_helices", minimum=0),
        )
    return _edit_result(model)


def add_deletion_in_range(
    design: Mapping[str, Any],
    *,
    helix: int,
    start: int,
    end: int,
    deletion_offset: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add one deletion on a helix at a free offset inside ``[start, end)``
    (``origami_rectangle.add_deletion_in_range``).

    Only offsets covered by exactly two domains, not at a domain end and not already a deletion or
    insertion are free. ``deletion_offset`` is relative to ``start``; negative picks the free
    offset furthest from the interval edges.
    """

    rectangle = require("scadnano.origami_rectangle", _EXTRA)
    model = _model(design)
    _helix(model, helix)
    with _guard("add the deletion"):
        rectangle.add_deletion_in_range(
            model,
            helix,
            _int(start, "start"),
            _int(end, "end"),
            _int(deletion_offset, "deletion_offset"),
        )
    return _edit_result(model)


__all__ = [
    "add_crossover",
    "add_deletion_in_range",
    "add_design_deletion",
    "add_design_helix",
    "add_design_insertion",
    "add_design_nick",
    "add_design_strand",
    "add_twist_correction_deletions",
    "assign_complement_sequence",
    "assign_dna_sequence",
    "assign_helix_domain_sequences",
    "assign_m13_scaffold",
    "autostaple_design",
    "build_plate_maps",
    "build_strand",
    "check_complementarity",
    "check_domain_overlap",
    "compute_base_pairs",
    "create_origami_rectangle",
    "domain_has_crossover",
    "export_cadnano_v2_design",
    "export_idt_order",
    "export_oxview_design",
    "export_scadnano_design",
    "find_domains_at_offset",
    "find_helix_domain_overlap_pairs",
    "find_overlapping_domains",
    "find_strand_by_name",
    "get_backbone_angle",
    "get_domain_subsequence",
    "get_helix_neighbors",
    "get_m13_sequence",
    "get_strand_sequence",
    "get_strand_without_modifications",
    "import_cadnano_v2_design",
    "inline_deletions_insertions",
    "insert_strand_domain",
    "ligate_design_nick",
    "list_crossover_points",
    "list_design_modifications",
    "list_helix_crossovers",
    "list_potential_crossovers",
    "list_strands_at_helix_end",
    "load_scadnano_design",
    "map_domain_offset_to_dna_index",
    "move_strands",
    "relax_helix_rolls",
    "remove_design_strand",
    "remove_strand_domain",
    "remove_strand_modification",
    "reverse_strands",
    "rotate_strand_domains",
    "set_strand_domains",
    "set_strand_modification",
    "set_strand_properties",
    "split_helix_strands",
]

_register_public_verbose(globals(), tuple(__all__))
