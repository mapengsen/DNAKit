"""Optional-backend DNA design, record, assembly and visualization adapters.

The functions in this module are explicit adapters over optional packages:

- ``evaluate_dnachisel_specs()``, ``optimize_dnachisel_record()`` and
  ``reverse_translate_protein()`` use DnaChisel.
- ``read_dnacauldron_parts()`` and ``write_dnacauldron_constructs()`` use DnaCauldron.
- ``design_pcr_primers()``, ``design_assembly_primers()`` and
  ``find_recombination_sites()`` use pydna.
- ``build_sbol3_document()``, ``read_sbol3_sequences()`` and
  ``validate_sbol3_document()`` use pySBOL3.
- ``read_snapgene_sequence()``, ``read_snapgene_features()`` and
  ``convert_snapgene_to_genbank()`` use snapgene_reader, which only reads files.
- ``build_feature_map()``, ``render_feature_map_svg()`` and
  ``build_genbank_feature_map()`` use DNA Features Viewer.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint. The
DnaChisel optimizer, DnaCauldron simulation, pySBOL3 inventory and DNA Features
Viewer summary already have adapters in :mod:`dnakit.integrations`; the pydna
Gibson product adapter is in :mod:`dnakit.bioinformatics.cloning`. This module
does not repeat them.

Coordinates are 0-based, with ``start`` inclusive and ``end`` exclusive, unless a
function says otherwise. Every function returns a JSON-safe dictionary.
"""

from __future__ import annotations

import importlib.metadata
import math
import re
import struct
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any
from xml.parsers.expat import ExpatError

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_MAX_DNA = 2_000_000
_MAX_FILE_BYTES = 50_000_000
_MAX_SPECS_LIMIT = 200
_MAX_FEATURES_LIMIT = 100_000
_MAX_PROTEIN = 10_000
_MAX_TEMPLATE = 100_000
_MAX_FRAGMENTS = 20
_MAX_COMPONENTS = 200
_MAX_PARTS = 200
_PROTEIN_PATTERN = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]+\*?")
_COLOR_PATTERN = re.compile(r"#[0-9A-Fa-f]{6}")
_DNACHISEL_SPECS = frozenset(
    {
        "AvoidChanges",
        "AvoidHairpins",
        "AvoidPattern",
        "AvoidStopCodons",
        "CodonOptimize",
        "EnforceChanges",
        "EnforceGCContent",
        "EnforceMeltingTemperature",
        "EnforceTranslation",
        "MaximizeCAI",
    }
)
_GENBANK_SUFFIXES = (".gb", ".gbk", ".genbank")
_SBOL_FORMATS = ("turtle", "ntriples", "rdfxml", "jsonld")


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a finite number in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str, *, maximum: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum or "\x00" in value:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    return value.strip()


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a DNA string.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    return sequence_text(value)


def _canonical_dna(value: object, name: str, *, maximum: int) -> str:
    text = _dna(value, name, maximum=maximum)
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _source_file(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("source must be a local file path.")
    path = input_path(value)
    if path.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"Input file exceeds {_MAX_FILE_BYTES} bytes: {path.name}")
    return path


def _output_file(value: object, name: str, suffixes: tuple[str, ...]) -> Path:
    text = _text(str(value) if isinstance(value, Path) else value, name, maximum=4096)
    path = Path(text).expanduser().resolve()
    if suffixes and path.suffix.lower() not in suffixes:
        raise ConfigurationError(f"{name} must end with {' or '.join(suffixes)}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new file: {path}")
    return path


def _provenance(distribution: str, **notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": distribution, "version": installed, **notes}


def _float_or_none(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _locations(value: object) -> list[list[int]]:
    rows: list[list[int]] = []
    for location in value or ():
        start = getattr(location, "start", None)
        end = getattr(location, "end", None)
        if start is None or end is None:
            continue
        rows.append([int(start), int(end)])
    return rows


def _strand(value: object) -> int | None:
    if value is None:
        return None
    return int(value)


# DnaChisel


def _dnachisel_specs(value: object, field: str, dnachisel: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{field} must be a list of specification mappings.")
    if len(value) > _MAX_SPECS_LIMIT:
        raise ConfigurationError(f"{field} accepts at most {_MAX_SPECS_LIMIT} specifications.")
    built: list[Any] = []
    for item in value:
        if not isinstance(item, Mapping) or not isinstance(item.get("name"), str):
            raise ConfigurationError(f"Each {field} entry needs a string 'name' key.")
        name = item["name"]
        if name not in _DNACHISEL_SPECS:
            allowed = ", ".join(sorted(_DNACHISEL_SPECS))
            raise ConfigurationError(
                f"Unsupported DnaChisel specification {name!r}; allowed: {allowed}."
            )
        kwargs: dict[str, Any] = {}
        for key, argument in item.items():
            if key == "name":
                continue
            if key == "location" and isinstance(argument, list):
                argument = tuple(argument)
            kwargs[str(key)] = argument
        try:
            built.append(getattr(dnachisel, name)(**kwargs))
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(f"Invalid DnaChisel specification {name!r}: {exc}") from exc
    return built


def _dnachisel_problem(
    dnachisel: Any, sequence: Any, constraints: list, objectives: list, circular: bool
) -> Any:
    class_name = "CircularDnaOptimizationProblem" if circular else "DnaOptimizationProblem"
    problem_class = getattr(dnachisel, class_name, None)
    if problem_class is None:
        raise ConfigurationError(f"The installed DnaChisel package lacks {class_name}.")
    return problem_class(
        sequence=sequence, constraints=constraints, objectives=objectives, logger=None
    )


def _evaluation_rows(evaluations: Iterable[Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in evaluations:
        specification = item.specification
        rows.append(
            {
                "name": type(specification).__name__,
                "label": str(specification),
                "passes": bool(item.passes),
                "score": _float_or_none(item.score),
                "locations": _locations(item.locations),
                "message": str(item.message),
            }
        )
    return rows


def evaluate_dnachisel_specs(
    sequence: str,
    *,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score DnaChisel constraints and objectives without changing the sequence.

    Each specification is a mapping with a ``name`` from the allow-list and
    keyword arguments, for example ``{"name": "AvoidPattern", "pattern": "GGTCTC"}``.
    """

    dnachisel = require("dnachisel", "bio-dnachisel")
    text = _dna(sequence, "sequence")
    circular = _flag(circular, "circular")
    built_constraints = _dnachisel_specs(constraints, "constraints", dnachisel)
    built_objectives = _dnachisel_specs(objectives, "objectives", dnachisel)
    if not built_constraints and not built_objectives:
        raise ConfigurationError("Provide at least one constraint or objective.")
    try:
        problem = _dnachisel_problem(dnachisel, text, built_constraints, built_objectives, circular)
        constraint_rows = _evaluation_rows(problem.constraints_evaluations())
        objective_rows = _evaluation_rows(problem.objectives_evaluations())
        all_pass = bool(problem.all_constraints_pass())
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaChisel evaluation failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "evaluate_dnachisel_specs",
        "sequence_length": len(text),
        "circular": circular,
        "all_constraints_pass": all_pass,
        "failed_constraints": [row["label"] for row in constraint_rows if not row["passes"]],
        "constraints": constraint_rows,
        "objectives": objective_rows,
        "provenance": _provenance("dnachisel", optimization="none; the sequence is only scored"),
    }


def optimize_dnachisel_record(
    source: str | Path,
    *,
    output_path: str | Path,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Optimize a local GenBank, FASTA or SnapGene record and write annotated GenBank.

    The output file is new and is never overwritten. DnaChisel adds constraint
    and objective annotations to the GenBank features it writes.
    """

    dnachisel = require("dnachisel", "bio-dnachisel")
    path = _source_file(source)
    target = _output_file(output_path, "output_path", _GENBANK_SUFFIXES)
    circular = _flag(circular, "circular")
    built_constraints = _dnachisel_specs(constraints, "constraints", dnachisel)
    built_objectives = _dnachisel_specs(objectives, "objectives", dnachisel)
    if not built_constraints and not built_objectives:
        raise ConfigurationError("Provide at least one constraint or objective.")
    record_id = re.sub(r"[^A-Za-z0-9_]", "_", path.stem)[:16] or "record"
    try:
        record = dnachisel.load_record(
            str(path), linear=not circular, name=record_id, file_format="auto"
        )
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaChisel could not read {path.name}: {exc}") from exc
    if len(record.seq) > _MAX_DNA:
        raise ConfigurationError(f"Record exceeds {_MAX_DNA} bases.")
    try:
        problem = _dnachisel_problem(
            dnachisel, record, built_constraints, built_objectives, circular
        )
        problem.resolve_constraints()
        if built_objectives:
            problem.optimize()
        all_pass = bool(problem.all_constraints_pass())
        edits = int(problem.number_of_edits())
        problem.to_record(
            str(target),
            with_constraints=True,
            with_objectives=True,
            record_id=record_id,
        )
    except (dnachisel.NoSolutionError, OSError, TypeError, ValueError, RuntimeError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"DnaChisel optimization failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "optimize_dnachisel_record",
        "source": str(path),
        "output_path": str(target),
        "output_format": "genbank",
        "sequence_length": len(problem.sequence),
        "edits": edits,
        "all_constraints_pass": all_pass,
        "provenance": _provenance(
            "dnachisel",
            optimization="resolve constraints; optimize objectives when given",
            output="GenBank with DnaChisel constraint and objective annotations",
        ),
    }


def reverse_translate_protein(
    protein: str,
    *,
    table: str = "Standard",
    max_residues: int = _MAX_PROTEIN,
    verbose: bool = True,
) -> dict[str, Any]:
    """Back-translate a protein with DnaChisel, choosing the first codon in table order.

    The result is deterministic. It does not use codon usage; use
    ``optimize_dnachisel_record()`` with ``CodonOptimize`` for host-specific codons.
    A trailing ``*`` is kept as a stop codon.
    """

    dnachisel = require("dnachisel", "bio-dnachisel")
    limit = _integer(max_residues, "max_residues", minimum=1, maximum=_MAX_PROTEIN)
    text = _text(protein, "protein", maximum=limit + 1).upper()
    if not _PROTEIN_PATTERN.fullmatch(text) or len(text.rstrip("*")) > limit:
        raise ConfigurationError(
            "protein must contain only the 20 standard residues, optionally ending in '*'."
        )
    table_name = _text(table, "table", maximum=64)
    table_names = tuple(getattr(getattr(dnachisel, "biotools", None), "CODON_TABLE_NAMES", ()))
    if table_names and table_name not in table_names:
        raise ConfigurationError(f"Unknown codon table {table_name!r}.")
    try:
        dna = str(dnachisel.reverse_translate(text, table=table_name))
    except (KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DnaChisel could not back-translate the protein: {exc}") from exc
    return {
        "status": "completed",
        "operation": "reverse_translate_protein",
        "protein_length": len(text),
        "sequence": dna,
        "sequence_length": len(dna),
        "table": table_name,
        "randomized": False,
        "provenance": _provenance("dnachisel", codon_choice="first codon in table order"),
    }


# DnaCauldron


def read_dnacauldron_parts(
    sources: Sequence[str | Path],
    *,
    topology: str = "default_to_linear",
    max_parts: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load local GenBank, FASTA or SnapGene part files through DnaCauldron.

    Each part reports its ID, length, sequence and the overhangs DnaCauldron
    finds in its annotations. Part IDs must be unique.
    """

    dnacauldron = require("dnacauldron", "bio-dnacauldron")
    if isinstance(sources, (str, bytes, Path)) or not isinstance(sources, Sequence):
        raise ConfigurationError("sources must be a list of local file paths.")
    limit = _integer(max_parts, "max_parts", minimum=1, maximum=_MAX_PARTS)
    if not sources or len(sources) > limit:
        raise ConfigurationError(f"sources must contain between 1 and {limit} files.")
    if topology not in {"default_to_linear", "default_to_circular", "linear", "circular"}:
        raise ConfigurationError(
            "topology must be default_to_linear, default_to_circular, linear or circular."
        )
    parts: list[dict[str, Any]] = []
    seen: set[str] = set()
    total = 0
    for source in sources:
        path = _source_file(source)
        try:
            record = dnacauldron.load_record(str(path), topology=topology)
            overhangs = list(dnacauldron.list_overhangs_from_record_annotations(record))
        except (OSError, TypeError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(f"DnaCauldron could not read {path.name}: {exc}") from exc
        sequence = str(record.seq).upper()
        total += len(sequence)
        if total > _MAX_DNA:
            raise ConfigurationError(f"Parts exceed {_MAX_DNA} bases in total.")
        part_id = str(record.id)
        if part_id in seen:
            raise ConfigurationError(f"Duplicate part ID {part_id!r}.")
        seen.add(part_id)
        parts.append(
            {
                "id": part_id,
                "file": path.name,
                "length": len(sequence),
                "topology": str(record.annotations.get("topology", "linear")),
                "overhangs": [str(item) for item in overhangs],
                "sequence": sequence,
            }
        )
    return {
        "status": "completed",
        "operation": "read_dnacauldron_parts",
        "part_count": len(parts),
        "parts": parts,
        "provenance": _provenance("dnacauldron", topology=topology),
    }


def write_dnacauldron_constructs(
    parts: Mapping[str, str],
    *,
    output_dir: str | Path,
    method: str = "gibson",
    enzyme: str = "auto",
    max_constructs: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate a DnaCauldron assembly and write each construct as an annotated GenBank file.

    ``method`` is ``"gibson"`` or ``"type2s"``. The output directory is new and
    is created by this call. Unlike ``simulate_dnacauldron_assembly()``, the
    constructs keep their feature annotations in the written files.
    """

    dnacauldron = require("dnacauldron", "bio-dnacauldron")
    if not isinstance(parts, Mapping) or not 2 <= len(parts) <= _MAX_PARTS:
        raise ConfigurationError(f"parts must contain between 2 and {_MAX_PARTS} named sequences.")
    normalized: dict[str, str] = {}
    total = 0
    for name, value in parts.items():
        part_name = _text(name, "part name", maximum=64)
        sequence = _dna(value, f"parts[{part_name!r}]", maximum=_MAX_DNA)
        total += len(sequence)
        if total > _MAX_DNA:
            raise ConfigurationError(f"Assembly input exceeds {_MAX_DNA} bases in total.")
        if part_name in normalized:
            raise ConfigurationError(f"Duplicate part name {part_name!r}.")
        normalized[part_name] = sequence
    if method not in {"gibson", "type2s"}:
        raise ConfigurationError("method must be 'gibson' or 'type2s'.")
    enzyme_name = _text(enzyme, "enzyme", maximum=64)
    limit = _integer(max_constructs, "max_constructs", minimum=1, maximum=100)
    directory = (
        Path(
            _text(
                str(output_dir) if isinstance(output_dir, Path) else output_dir,
                "output_dir",
                maximum=4096,
            )
        )
        .expanduser()
        .resolve()
    )
    if not directory.parent.is_dir():
        raise ConfigurationError(f"output_dir parent does not exist: {directory.parent}")
    if directory.exists():
        raise ConfigurationError(f"output_dir already exists; choose a new directory: {directory}")

    repository = dnacauldron.SequenceRepository()
    for name, sequence in normalized.items():
        repository.add_record((name, sequence))
    common = {
        "parts": list(normalized),
        "max_constructs": limit,
        "expected_constructs": "any_number",
        "expect_no_unused_parts": False,
    }
    if method == "gibson":
        assembly = dnacauldron.GibsonAssembly(**common)
    else:
        assembly = dnacauldron.Type2sRestrictionAssembly(enzyme=enzyme_name, **common)
    try:
        simulation = assembly.simulate(sequence_repository=repository)
        records = list(simulation.construct_records)[:limit]
        warnings = [str(item) for item in simulation.warnings]
        errors = [str(item) for item in simulation.errors]
    except (AttributeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaCauldron assembly simulation failed: {exc}") from exc
    directory.mkdir(parents=True, exist_ok=False)
    files: list[dict[str, Any]] = []
    try:
        for index, record in enumerate(records, start=1):
            safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", str(record.id))[:60] or "construct"
            target = directory / f"{index:03d}_{safe_id}.gb"
            dnacauldron.write_record(record, str(target), fmt="genbank")
            files.append(
                {
                    "id": str(record.id),
                    "file": target.name,
                    "length": len(record.seq),
                    "topology": str(record.annotations.get("topology", "")),
                }
            )
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaCauldron could not write constructs: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_dnacauldron_constructs",
        "method": method,
        "parts": list(normalized),
        "output_dir": str(directory),
        "constructs": files,
        "warnings": warnings,
        "errors": errors,
        "provenance": _provenance(
            "dnacauldron",
            output="GenBank files; prediction only, no wet-lab result is inferred",
        ),
    }


# pydna


def _primer_design(
    pydna_design: Any, dseqrecord: Any, text: str, target_tm: float, limit: int
) -> Any:
    try:
        return pydna_design.primer_design(dseqrecord(text), target_tm=target_tm, limit=limit)
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not design primers: {exc}") from exc


def design_pcr_primers(
    template: str,
    *,
    target_tm: float = 55.0,
    min_primer_length: int = 13,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design a forward and reverse primer pair that spans a linear template with pydna.

    The primers start at the template ends, so the product is the whole template.
    ``target_tm`` is passed to pydna's primer design.
    """

    pydna_design = require("pydna.design", "bio-assembly")
    dseqrecord = require("pydna.dseqrecord", "bio-assembly").Dseqrecord
    text = _dna(template, "template", maximum=_MAX_TEMPLATE)
    tm = _number(target_tm, "target_tm", minimum=40.0, maximum=75.0)
    limit = _integer(min_primer_length, "min_primer_length", minimum=10, maximum=40)
    if len(text) < 2 * limit:
        raise ConfigurationError("template is too short for two primers of min_primer_length.")
    amplicon = _primer_design(pydna_design, dseqrecord, text, tm, limit)
    forward = str(amplicon.forward_primer.seq)
    reverse = str(amplicon.reverse_primer.seq)
    return {
        "status": "completed",
        "operation": "design_pcr_primers",
        "template_length": len(text),
        "forward_primer": forward,
        "reverse_primer": reverse,
        "forward_length": len(forward),
        "reverse_length": len(reverse),
        "product_length": len(amplicon),
        "provenance": _provenance("pydna", design="pydna.design.primer_design", target_tm=tm),
    }


def design_assembly_primers(
    templates: Sequence[str],
    *,
    overlap: int = 35,
    max_link: int = 40,
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add overlap tails to primers so that adjacent PCR fragments can be assembled.

    Each template is first given a primer pair with pydna, then pydna's
    ``assembly_fragments()`` adds tails for Gibson assembly or fusion PCR. The
    fragments are joined in the order given.
    """

    pydna_design = require("pydna.design", "bio-assembly")
    dseqrecord = require("pydna.dseqrecord", "bio-assembly").Dseqrecord
    if isinstance(templates, (str, bytes)) or not isinstance(templates, Sequence):
        raise ConfigurationError("templates must be a list of DNA strings.")
    if not 2 <= len(templates) <= _MAX_FRAGMENTS:
        raise ConfigurationError(
            f"templates must contain between 2 and {_MAX_FRAGMENTS} fragments."
        )
    overlap_bases = _integer(overlap, "overlap", minimum=10, maximum=100)
    link = _integer(max_link, "max_link", minimum=overlap_bases, maximum=200)
    flag = _flag(circular, "circular")
    texts = [
        _dna(item, f"templates[{index}]", maximum=_MAX_TEMPLATE)
        for index, item in enumerate(templates)
    ]
    if sum(len(item) for item in texts) > _MAX_DNA:
        raise ConfigurationError(f"Templates exceed {_MAX_DNA} bases in total.")
    amplicons = [_primer_design(pydna_design, dseqrecord, text, 55.0, 13) for text in texts]
    try:
        fragments = pydna_design.assembly_fragments(
            amplicons, overlap=overlap_bases, maxlink=link, circular=flag
        )
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not add assembly tails: {exc}") from exc
    rows = [
        {
            "index": index,
            "forward_primer": str(item.forward_primer.seq),
            "reverse_primer": str(item.reverse_primer.seq),
            "product_length": len(item),
        }
        for index, item in enumerate(fragments)
    ]
    return {
        "status": "completed",
        "operation": "design_assembly_primers",
        "overlap": overlap_bases,
        "max_link": link,
        "circular": flag,
        "fragments": rows,
        "provenance": _provenance(
            "pydna",
            design="pydna.design.assembly_fragments after pydna.design.primer_design",
        ),
    }


def find_recombination_sites(
    sequence: str,
    *,
    system: str = "loxp",
    greedy: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find loxP or Gateway att sites in a DNA sequence with pydna.

    ``system`` is ``"loxp"`` or ``"gateway"``. Positions are 0-based half-open
    intervals on the forward strand. ``greedy`` is passed to pydna's Gateway
    site finder and is ignored for loxP.
    """

    dseqrecord = require("pydna.dseqrecord", "bio-assembly").Dseqrecord
    text = _dna(sequence, "sequence")
    if system not in {"loxp", "gateway"}:
        raise ConfigurationError("system must be 'loxp' or 'gateway'.")
    flag = _flag(greedy, "greedy")
    record = dseqrecord(text)
    try:
        if system == "loxp":
            found = require("pydna.cre_lox", "bio-assembly").find_loxP_sites(record)
        else:
            found = require("pydna.gateway", "bio-assembly").find_gateway_sites(record, greedy=flag)
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna site search failed: {exc}") from exc
    sites = [
        {
            "name": str(name),
            "start": int(location.start),
            "end": int(location.end),
            "strand": _strand(location.strand),
        }
        for name, locations in found.items()
        for location in locations
    ]
    sites.sort(key=lambda row: (row["start"], row["end"], row["name"]))
    return {
        "status": "completed",
        "operation": "find_recombination_sites",
        "system": system,
        "sequence_length": len(text),
        "site_count": len(sites),
        "sites": sites,
        "provenance": _provenance("pydna", greedy=flag if system == "gateway" else None),
    }


# pySBOL3


def _sbol_format(value: object, sbol3: Any) -> tuple[str, str]:
    if value not in _SBOL_FORMATS:
        raise ConfigurationError(f"file_format must be one of {', '.join(_SBOL_FORMATS)}.")
    mapping = {
        "turtle": sbol3.TURTLE,
        "ntriples": sbol3.NTRIPLES,
        "rdfxml": sbol3.RDF_XML,
        "jsonld": sbol3.JSONLD,
    }
    return str(value), mapping[str(value)]


def _namespace(value: object) -> str:
    text = _text(value, "namespace", maximum=512)
    if not text.startswith(("http://", "https://")) or text[-1] not in "/#":
        raise ConfigurationError("namespace must be an http(s) URI ending in '/' or '#'.")
    return text


def _issue_row(issue: Any) -> dict[str, Any]:
    return {
        "object": None if getattr(issue, "object_id", None) is None else str(issue.object_id),
        "rule": None if getattr(issue, "rule_id", None) is None else str(issue.rule_id),
        "message": str(getattr(issue, "message", issue)),
    }


def build_sbol3_document(
    components: Mapping[str, str],
    *,
    namespace: str,
    file_format: str = "turtle",
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an SBOL3 document with one DNA Component and Sequence per named sequence.

    Each display ID becomes ``<namespace><display_id>``; its Sequence is
    ``<namespace><display_id>_sequence``. Without ``output_path`` the serialized
    document is returned in the ``document`` field. With ``output_path`` the new
    file is written and the field is None.
    """

    sbol3 = require("sbol3", "bio-pysbol3")
    ns = _namespace(namespace)
    format_name, format_constant = _sbol_format(file_format, sbol3)
    if not isinstance(components, Mapping) or not 1 <= len(components) <= _MAX_COMPONENTS:
        raise ConfigurationError(
            f"components must contain between 1 and {_MAX_COMPONENTS} entries."
        )
    target = None
    if output_path is not None:
        target = _output_file(
            output_path, "output_path", (".ttl", ".nt", ".xml", ".rdf", ".json", ".jsonld")
        )
    document = sbol3.Document()
    ids: list[str] = []
    total = 0
    for display_id, value in components.items():
        if not isinstance(display_id, str) or not sbol3.is_valid_display_id(display_id):
            raise ConfigurationError(f"Invalid SBOL display ID {display_id!r}.")
        text = _dna(value, f"components[{display_id!r}]", maximum=_MAX_DNA)
        total += len(text)
        if total > _MAX_DNA:
            raise ConfigurationError(f"Components exceed {_MAX_DNA} bases in total.")
        sequence = sbol3.Sequence(
            ns + display_id + "_sequence",
            elements=text,
            encoding=sbol3.IUPAC_DNA_ENCODING,
            namespace=ns,
        )
        component = sbol3.Component(
            ns + display_id,
            sbol3.SBO_DNA,
            sequences=[sequence],
            namespace=ns,
            name=display_id,
        )
        document.add(sequence)
        document.add(component)
        ids.append(display_id)
    report = document.validate()
    try:
        if target is None:
            text_output = document.write_string(format_constant)
        else:
            document.write(str(target), format_constant)
            text_output = None
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"pySBOL3 could not serialize the document: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_sbol3_document",
        "namespace": ns,
        "file_format": format_name,
        "component_count": len(ids),
        "components": ids,
        "sequence_bases": total,
        "validation": {
            "error_count": len(report.errors),
            "warning_count": len(report.warnings),
            "errors": [_issue_row(item) for item in report.errors],
            "warnings": [_issue_row(item) for item in report.warnings],
        },
        "output_path": None if target is None else str(target),
        "document": text_output,
        "provenance": _provenance(
            "sbol3", standard="SBOL3", identity_rule="namespace + display ID"
        ),
    }


def _read_sbol3(sbol3: Any, path: Path) -> Any:
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, SyntaxError, TypeError, ValueError, RuntimeError, sbol3.SBOLError) as exc:
        raise ConfigurationError(f"pySBOL3 could not read {path.name}: {exc}") from exc
    return document


def read_sbol3_sequences(
    source: str | Path,
    *,
    max_objects: int = 50_000,
    max_sequences: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every top-level SBOL3 Sequence from a local document, with its elements."""

    sbol3 = require("sbol3", "bio-pysbol3")
    path = _source_file(source)
    object_limit = _integer(max_objects, "max_objects", minimum=1, maximum=500_000)
    sequence_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=10_000)
    document = _read_sbol3(sbol3, path)
    objects = 0
    rows: list[dict[str, Any]] = []
    for item in document:
        objects += 1
        if objects > object_limit:
            raise ConfigurationError(f"SBOL3 document exceeds max_objects={object_limit}.")
        if not isinstance(item, sbol3.Sequence):
            continue
        if len(rows) >= sequence_limit:
            raise ConfigurationError(f"SBOL3 document exceeds max_sequences={sequence_limit}.")
        elements = item.elements or ""
        rows.append(
            {
                "identity": str(item.identity),
                "display_id": None if item.display_id is None else str(item.display_id),
                "encoding": None if item.encoding is None else str(item.encoding),
                "length": len(elements),
                "elements": elements,
            }
        )
    return {
        "status": "completed",
        "operation": "read_sbol3_sequences",
        "source": str(path),
        "object_count": objects,
        "sequence_count": len(rows),
        "sequences": rows,
        "provenance": _provenance("sbol3", standard="SBOL3", scope="top-level Sequence objects"),
    }


def validate_sbol3_document(
    source: str | Path,
    *,
    max_issues: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a local SBOL3 document with pySBOL3 and list its errors and warnings."""

    sbol3 = require("sbol3", "bio-pysbol3")
    path = _source_file(source)
    limit = _integer(max_issues, "max_issues", minimum=1, maximum=10_000)
    document = _read_sbol3(sbol3, path)
    report = document.validate()
    return {
        "status": "completed",
        "operation": "validate_sbol3_document",
        "source": str(path),
        "valid": not report.errors,
        "error_count": len(report.errors),
        "warning_count": len(report.warnings),
        "errors": [_issue_row(item) for item in report.errors[:limit]],
        "warnings": [_issue_row(item) for item in report.warnings[:limit]],
        "provenance": _provenance("sbol3", validator="pySBOL3 Document.validate"),
    }


# snapgene_reader


def _snapgene_dict(snapgene: Any, path: Path) -> dict[str, Any]:
    try:
        data = snapgene.snapgene_file_to_dict(filepath=str(path))
    except (
        ExpatError,
        IndexError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        struct.error,
    ) as exc:
        raise ConfigurationError(f"snapgene_reader could not parse {path.name}: {exc}") from exc
    return data


def read_snapgene_sequence(
    source: str | Path, *, max_bases: int = _MAX_DNA, verbose: bool = True
) -> dict[str, Any]:
    """Read the DNA sequence and topology from a SnapGene .dna file.

    snapgene_reader only reads files; this function does not write .dna files.
    """

    snapgene = require("snapgene_reader", "bio-snapgene")
    path = _source_file(source)
    limit = _integer(max_bases, "max_bases", minimum=1, maximum=_MAX_DNA)
    data = _snapgene_dict(snapgene, path)
    if "dna" not in data or "seq" not in data:
        raise ConfigurationError(f"{path.name} has no DNA sequence block.")
    sequence = str(data["seq"])
    if len(sequence) > limit:
        raise ConfigurationError(f"SnapGene sequence exceeds max_bases={limit}.")
    dna = data["dna"]
    return {
        "status": "completed",
        "operation": "read_snapgene_sequence",
        "source": str(path),
        "sequence": sequence,
        "length": len(sequence),
        "topology": str(dna["topology"]),
        "strandedness": str(dna["strandedness"]),
        "dam_methylated": bool(dna["damMethylated"]),
        "dcm_methylated": bool(dna["dcmMethylated"]),
        "eco_ki_methylated": bool(dna["ecoKIMethylated"]),
        "export_version": int(data["exportVersion"]),
        "import_version": int(data["importVersion"]),
        "provenance": _provenance("snapgene_reader", format="SnapGene .dna, read only"),
    }


def read_snapgene_features(
    source: str | Path, *, max_features: int = 10_000, verbose: bool = True
) -> dict[str, Any]:
    """List the features of a SnapGene .dna file with their 0-based coordinates.

    ``strand`` is ``"+"``, ``"-"``, ``"."`` (none) or ``"="`` (both), as stored by SnapGene.
    """

    snapgene = require("snapgene_reader", "bio-snapgene")
    path = _source_file(source)
    limit = _integer(max_features, "max_features", minimum=1, maximum=_MAX_FEATURES_LIMIT)
    data = _snapgene_dict(snapgene, path)
    features = list(data.get("features", []))
    if len(features) > limit:
        raise ConfigurationError(f"SnapGene file exceeds max_features={limit}.")
    rows: list[dict[str, Any]] = []
    for feature in features:
        qualifiers = feature.get("qualifiers", {})
        label = qualifiers.get("label") if isinstance(qualifiers, Mapping) else None
        rows.append(
            {
                "name": str(feature["name"]),
                "type": str(feature["type"]),
                "start": int(feature["start"]),
                "end": int(feature["end"]),
                "strand": str(feature["strand"]),
                "color": str(feature["color"]),
                "label": label if isinstance(label, str) else str(feature["name"]),
            }
        )
    return {
        "status": "completed",
        "operation": "read_snapgene_features",
        "source": str(path),
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(
            "snapgene_reader", coordinates="0-based, start inclusive, end exclusive"
        ),
    }


def convert_snapgene_to_genbank(
    source: str | Path, output_path: str | Path, *, verbose: bool = True
) -> dict[str, Any]:
    """Convert a SnapGene .dna file to a new GenBank file with snapgene_reader.

    The GenBank text is produced by snapgene_reader's own converter. Its LOCUS
    date is fixed by that converter.
    """

    snapgene = require("snapgene_reader", "bio-snapgene")
    path = _source_file(source)
    target = _output_file(output_path, "output_path", _GENBANK_SUFFIXES)
    try:
        with path.open("rb") as reader, target.open("w", encoding="utf-8") as writer:
            snapgene.snapgene_file_to_gbk(reader, writer)
    except (ExpatError, IndexError, KeyError, OSError, TypeError, ValueError, struct.error) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"snapgene_reader could not convert {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "convert_snapgene_to_genbank",
        "source": str(path),
        "output_path": str(target),
        "bytes": target.stat().st_size,
        "provenance": _provenance("snapgene_reader", output="GenBank text from snapgene_reader"),
    }


# DNA Features Viewer


def _feature_rows(features: object, sequence_length: int, limit: int) -> list[dict[str, Any]]:
    if isinstance(features, (str, bytes)) or not isinstance(features, Sequence):
        raise ConfigurationError("features must be a list of mappings with start and end.")
    if len(features) > limit:
        raise ConfigurationError(f"features exceed max_features={limit}.")
    rows: list[dict[str, Any]] = []
    for item in features:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Each feature must be a mapping with start and end.")
        start = _integer(item.get("start"), "feature start", minimum=0, maximum=sequence_length - 1)
        end = _integer(item.get("end"), "feature end", minimum=start + 1, maximum=sequence_length)
        strand = item.get("strand", 0)
        if isinstance(strand, bool) or strand not in (1, -1, 0):
            raise ConfigurationError("feature strand must be 1, -1 or 0.")
        label = item.get("label", "")
        if label is None:
            label = ""
        if not isinstance(label, str) or len(label) > 200 or "\x00" in label:
            raise ConfigurationError("feature label must be a string of at most 200 characters.")
        color = item.get("color", "#000080")
        if not isinstance(color, str) or not _COLOR_PATTERN.fullmatch(color):
            raise ConfigurationError("feature color must be a #RRGGBB string.")
        rows.append(
            {"label": label, "start": start, "end": end, "strand": int(strand), "color": color}
        )
    return rows


def _graphic_features(dfv: Any, rows: list[dict[str, Any]]) -> list[Any]:
    return [
        dfv.GraphicFeature(
            start=row["start"],
            end=row["end"],
            strand=row["strand"] or None,
            label=row["label"] or None,
            color=row["color"],
        )
        for row in rows
    ]


def _feature_levels(graphic: list[Any]) -> list[float | None]:
    levels_module = require(
        "dna_features_viewer.compute_features_levels", "bio-dna-features-viewer"
    )
    levels = levels_module.compute_features_levels(graphic)
    return [_float_or_none(levels.get(item)) for item in graphic]


def build_feature_map(
    sequence_length: int,
    features: Sequence[Mapping[str, Any]],
    *,
    max_features: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return feature-map data (coordinates, colors and lanes) for a linear sequence.

    ``lane`` is the vertical level that DNA Features Viewer assigns so that
    overlapping features do not collide. Lanes start at 0.0.
    """

    dfv = require("dna_features_viewer", "bio-dna-features-viewer")
    length = _integer(sequence_length, "sequence_length", minimum=1, maximum=_MAX_DNA)
    limit = _integer(max_features, "max_features", minimum=1, maximum=_MAX_FEATURES_LIMIT)
    rows = _feature_rows(features, length, limit)
    graphic = _graphic_features(dfv, rows)
    lanes = _feature_levels(graphic)
    result = [{**row, "lane": lane} for row, lane in zip(rows, lanes, strict=True)]
    known = [lane for lane in lanes if lane is not None]
    return {
        "status": "completed",
        "operation": "build_feature_map",
        "sequence_length": length,
        "feature_count": len(result),
        "lane_count": int(max(known)) + 1 if known else 0,
        "features": result,
        "provenance": _provenance(
            "dna-features-viewer",
            coordinates="0-based, start inclusive, end exclusive",
            lanes="compute_features_levels from DNA Features Viewer",
        ),
    }


def render_feature_map_svg(
    sequence_length: int,
    features: Sequence[Mapping[str, Any]],
    *,
    output_path: str | Path,
    figure_width: float = 8.0,
    max_features: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a linear feature map with DNA Features Viewer and write it as a new SVG file.

    The Matplotlib figure is rendered with an Agg canvas so no display is needed.
    The SVG date metadata is omitted, so the same input gives the same bytes.
    """

    dfv = require("dna_features_viewer", "bio-dna-features-viewer")
    matplotlib = require("matplotlib", "bio-plots")
    figure_module = require("matplotlib.figure", "bio-plots")
    canvas_module = require("matplotlib.backends.backend_agg", "bio-plots")
    length = _integer(sequence_length, "sequence_length", minimum=1, maximum=_MAX_DNA)
    limit = _integer(max_features, "max_features", minimum=1, maximum=_MAX_FEATURES_LIMIT)
    width = _number(figure_width, "figure_width", minimum=2.0, maximum=20.0)
    target = _output_file(output_path, "output_path", (".svg",))
    rows = _feature_rows(features, length, limit)
    record = dfv.GraphicRecord(sequence_length=length, features=_graphic_features(dfv, rows))
    figure = figure_module.Figure(figsize=(width, 2.0))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot()
    try:
        record.plot(ax=axes)
        # A fixed hash salt makes SVG IDs repeat; rc_context leaves global rcParams unchanged.
        with matplotlib.rc_context({"svg.hashsalt": "dnakit-feature-map"}):
            figure.savefig(str(target), format="svg", metadata={"Date": None})
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"DNA Features Viewer could not render the map: {exc}") from exc
    return {
        "status": "completed",
        "operation": "render_feature_map_svg",
        "output_path": str(target),
        "format": "svg",
        "bytes": target.stat().st_size,
        "sequence_length": length,
        "feature_count": len(rows),
        "provenance": _provenance(
            "dna-features-viewer", renderer="matplotlib Agg canvas, SVG output"
        ),
    }


def build_genbank_feature_map(
    source: str | Path, *, max_features: int = 10_000, verbose: bool = True
) -> dict[str, Any]:
    """Return full feature-map data (coordinates, colors and lanes) for a local GenBank record.

    Unlike ``summarize_dna_features_viewer()``, which previews only the first
    features, every feature is returned here, with its lane.
    """

    dfv = require("dna_features_viewer", "bio-dna-features-viewer")
    path = _source_file(source)
    limit = _integer(max_features, "max_features", minimum=1, maximum=_MAX_FEATURES_LIMIT)
    translator_class = getattr(dfv, "BiopythonTranslator", None)
    if not callable(translator_class):
        raise ConfigurationError(
            "The installed DNA Features Viewer package lacks BiopythonTranslator."
        )
    try:
        record = translator_class().translate_record(str(path))
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DNA Features Viewer could not read {path.name}: {exc}") from exc
    graphic = list(record.features)
    if len(graphic) > limit:
        raise ConfigurationError(f"GenBank record exceeds max_features={limit}.")
    rows: list[dict[str, Any]] = []
    for feature in graphic:
        if feature.start is None or feature.end is None:
            raise ConfigurationError("DNA Features Viewer returned a feature without bounds.")
        rows.append(
            {
                "label": "" if feature.label is None else str(feature.label),
                "start": int(feature.start),
                "end": int(feature.end),
                "strand": _strand(feature.strand) or 0,
                "color": str(feature.color),
            }
        )
    lanes = _feature_levels(graphic)
    known = [lane for lane in lanes if lane is not None]
    return {
        "status": "completed",
        "operation": "build_genbank_feature_map",
        "source": str(path),
        "sequence_length": int(record.sequence_length),
        "feature_count": len(rows),
        "lane_count": int(max(known)) + 1 if known else 0,
        "features": [{**row, "lane": lane} for row, lane in zip(rows, lanes, strict=True)],
        "provenance": _provenance(
            "dna-features-viewer",
            coordinates="0-based, start inclusive, end exclusive",
            reader="BiopythonTranslator",
        ),
    }


# pyGenomeViz


def _genbank_label(qualifiers: Mapping[str, Any]) -> str:
    for key in ("gene", "locus_tag", "product", "label"):
        value = qualifiers.get(key)
        if isinstance(value, list) and value:
            value = value[0]
        if isinstance(value, str) and value.strip():
            return value.strip()[:80]
    return ""


def read_genbank_genome_tracks(
    source: str | Path,
    *,
    feature_type: str = "CDS",
    max_records: int = 100,
    max_features: int = 50_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read GenBank records as genome tracks: one size and feature list per record.

    ``feature_type`` is one GenBank feature key such as ``"CDS"`` or ``"gene"``.
    Sizes come from pyGenomeViz's ``Genbank`` parser. The feature filter is
    DNAKit's own loop over the parsed records, because pyGenomeViz 0.4.4's
    ``get_seqid2features()`` reads ``SeqFeature.strand``, which Biopython 1.88
    no longer has. Multi-part locations use their overall span. Each feature is
    labelled from its ``gene``, ``locus_tag``, ``product`` or ``label`` qualifier.
    Coordinates are 0-based, end-exclusive.
    """

    parser = require("pygenomeviz.parser.genbank", "bio-pygenomeviz")
    path = _source_file(source)
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=1000)
    feature_limit = _integer(max_features, "max_features", minimum=1, maximum=500_000)
    key = _text(feature_type, "feature_type", maximum=32)
    try:
        parsed = parser.Genbank(str(path))
        sizes = parsed.get_seqid2size()
        records = list(parsed.records)
    except (OSError, TypeError, ValueError, RuntimeError, IndexError, KeyError) as exc:
        raise ConfigurationError(f"pyGenomeViz could not read {path.name}: {exc}") from exc
    if not sizes:
        raise ConfigurationError(f"{path.name} contains no GenBank records.")
    if len(sizes) > record_limit:
        raise ConfigurationError(f"{path.name} exceeds max_records={record_limit}.")
    tracks: list[dict[str, Any]] = []
    total = 0
    for seqid, length in sizes.items():
        rows: list[dict[str, Any]] = []
        for record in (item for item in records if str(item.id) == seqid):
            for feature in record.features:
                if feature.type != key or feature.location is None:
                    continue
                total += 1
                if total > feature_limit:
                    raise ConfigurationError(f"{path.name} exceeds max_features={feature_limit}.")
                rows.append(
                    {
                        "type": str(feature.type),
                        "start": int(feature.location.start),
                        "end": int(feature.location.end),
                        "strand": _strand(feature.location.strand) or 0,
                        "label": _genbank_label(feature.qualifiers),
                    }
                )
        tracks.append({"name": str(seqid), "length": int(length), "features": rows})
    return {
        "status": "completed",
        "operation": "read_genbank_genome_tracks",
        "source": str(path),
        "feature_type": key,
        "track_count": len(tracks),
        "feature_count": total,
        "tracks": tracks,
        "provenance": _provenance(
            "pygenomeviz",
            parser="pygenomeviz.parser.genbank.Genbank",
            coordinates="0-based, start inclusive, end exclusive",
        ),
    }


def render_genome_tracks_svg(
    tracks: Sequence[Mapping[str, Any]],
    *,
    output_path: str | Path,
    links: Sequence[Mapping[str, Any]] = (),
    fig_width: float = 15.0,
    max_tracks: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw feature tracks, with optional links between adjacent tracks, as a new SVG file.

    Each track is ``{"name", "length", "features"}``; each feature is
    ``{"start", "end", "strand", "label"}`` with strand 1, -1 or 0 (drawn as 1).
    Coordinates are 0-based, end-exclusive. Each link is
    ``{"track_a", "start_a", "end_a", "track_b", "start_b", "end_b"}`` and must
    join neighbouring tracks in the order given, because pyGenomeViz draws links
    only between adjacent tracks. pyGenomeViz writes a creation timestamp and
    random group IDs into the SVG, so the bytes differ between runs.
    """

    pgv = require("pygenomeviz", "bio-pygenomeviz")
    target = _output_file(output_path, "output_path", (".svg",))
    width = _number(fig_width, "fig_width", minimum=4.0, maximum=60.0)
    limit = _integer(max_tracks, "max_tracks", minimum=1, maximum=50)
    if isinstance(tracks, (str, bytes)) or not isinstance(tracks, Sequence) or not tracks:
        raise ConfigurationError("tracks must be a non-empty list of track mappings.")
    if len(tracks) > limit:
        raise ConfigurationError(f"tracks exceed max_tracks={limit}.")
    names: list[str] = []
    lengths: dict[str, int] = {}
    track_rows: list[list[dict[str, Any]]] = []
    feature_total = 0
    for track in tracks:
        if not isinstance(track, Mapping):
            raise ConfigurationError("Each track must be a mapping with name, length and features.")
        name = _text(track.get("name"), "track name", maximum=64)
        if name in lengths:
            raise ConfigurationError(f"Duplicate track name {name!r}.")
        length = _integer(
            track.get("length"), f"length of track {name!r}", minimum=1, maximum=_MAX_DNA
        )
        names.append(name)
        lengths[name] = length
        feature_rows = _feature_rows(track.get("features", []), length, 10_000)
        feature_total += len(feature_rows)
        track_rows.append(feature_rows)
    if feature_total > 10_000:
        raise ConfigurationError("tracks contain more than 10000 features in total.")
    if isinstance(links, (str, bytes)) or not isinstance(links, Sequence):
        raise ConfigurationError("links must be a list of link mappings.")
    if len(links) > 10_000:
        raise ConfigurationError("links accept at most 10000 entries.")
    index = {name: position for position, name in enumerate(names)}
    link_rows: list[tuple[str, int, int, str, int, int]] = []
    for link in links:
        if not isinstance(link, Mapping):
            raise ConfigurationError("Each link must be a mapping.")
        first = _text(link.get("track_a"), "track_a", maximum=64)
        second = _text(link.get("track_b"), "track_b", maximum=64)
        if first not in index or second not in index:
            raise ConfigurationError("Link tracks must be among the given track names.")
        if abs(index[first] - index[second]) != 1:
            raise ConfigurationError("Links must join adjacent tracks in the order given.")
        start_a = _integer(link.get("start_a"), "start_a", minimum=0, maximum=lengths[first] - 1)
        end_a = _integer(link.get("end_a"), "end_a", minimum=start_a + 1, maximum=lengths[first])
        start_b = _integer(link.get("start_b"), "start_b", minimum=0, maximum=lengths[second] - 1)
        end_b = _integer(link.get("end_b"), "end_b", minimum=start_b + 1, maximum=lengths[second])
        link_rows.append((first, start_a, end_a, second, start_b, end_b))
    genome = pgv.GenomeViz(fig_width=width)
    try:
        for name, rows in zip(names, track_rows, strict=True):
            feature_track = genome.add_feature_track(name, lengths[name])
            for row in rows:
                # pyGenomeViz draws strand 0 as forward.
                feature_track.add_feature(
                    row["start"], row["end"], row["strand"] or 1, label=row["label"]
                )
        for first, start_a, end_a, second, start_b, end_b in link_rows:
            genome.add_link((first, start_a, end_a), (second, start_b, end_b))
        genome.savefig(str(target))
    except (OSError, TypeError, ValueError, RuntimeError, KeyError, IndexError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"pyGenomeViz could not draw the tracks: {exc}") from exc
    return {
        "status": "completed",
        "operation": "render_genome_tracks_svg",
        "output_path": str(target),
        "format": "svg",
        "bytes": target.stat().st_size,
        "track_count": len(tracks),
        "feature_count": feature_total,
        "link_count": len(link_rows),
        "provenance": _provenance("pygenomeviz", renderer="GenomeViz.savefig", output="SVG"),
    }


# primer3-py (GPL-2.0). It is imported only inside these functions; no primer3 code is copied.


def _thermo_conditions(
    mv_conc: object, dv_conc: object, dntp_conc: object, dna_conc: object
) -> tuple[float, float, float, float]:
    return (
        _number(mv_conc, "mv_conc", minimum=0.0, maximum=1000.0),
        _number(dv_conc, "dv_conc", minimum=0.0, maximum=100.0),
        _number(dntp_conc, "dntp_conc", minimum=0.0, maximum=100.0),
        _number(dna_conc, "dna_conc", minimum=0.0, maximum=10_000.0),
    )


def _thermo_row(operation: str, result: Any, length: int, conditions: dict[str, float]) -> dict:
    found = bool(result.structure_found)

    def value(item: object) -> float | None:
        return _float_or_none(item) if found else None

    return {
        "status": "completed",
        "operation": operation,
        "sequence_length": length,
        "conditions": conditions,
        "structure_found": found,
        "tm_c": value(result.tm),
        "dg_cal_per_mol": value(result.dg),
        "dh_cal_per_mol": value(result.dh),
        "ds_cal_per_mol_k": value(result.ds),
    }


def calculate_primer3_tm(
    sequence: str,
    *,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Melting temperature of one primer with primer3-py (SantaLucia nearest-neighbour model).

    Concentrations are mM for monovalent (``mv_conc``), divalent (``dv_conc``)
    and dNTP (``dntp_conc``), and nM for the primer (``dna_conc``).
    """

    primer3 = require("primer3", "bio-primer3")
    text = _canonical_dna(sequence, "sequence", maximum=200)
    mv, dv, dntp, dna = _thermo_conditions(mv_conc, dv_conc, dntp_conc, dna_conc)
    try:
        tm = primer3.calc_tm(text, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 could not compute Tm: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_primer3_tm",
        "sequence_length": len(text),
        "tm_c": _float_or_none(tm),
        "conditions": {"mv_mM": mv, "dv_mM": dv, "dntp_mM": dntp, "dna_nM": dna},
        "provenance": _provenance("primer3-py", tm_method="santalucia", license="GPL-2.0"),
    }


def calculate_primer3_hairpin(
    sequence: str,
    *,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hairpin stability of one primer with primer3-py.

    The thermodynamic fields are None when primer3 finds no hairpin.
    Energies are in cal/mol and entropy in cal/(mol*K).
    """

    primer3 = require("primer3", "bio-primer3")
    text = _canonical_dna(sequence, "sequence", maximum=200)
    mv, dv, dntp, dna = _thermo_conditions(mv_conc, dv_conc, dntp_conc, dna_conc)
    temperature = _number(temp_c, "temp_c", minimum=0.0, maximum=100.0)
    loop = _integer(max_loop, "max_loop", minimum=3, maximum=100)
    try:
        result = primer3.calc_hairpin(
            text,
            mv_conc=mv,
            dv_conc=dv,
            dntp_conc=dntp,
            dna_conc=dna,
            temp_c=temperature,
            max_loop=loop,
        )
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 could not compute the hairpin: {exc}") from exc
    row = _thermo_row(
        "calculate_primer3_hairpin",
        result,
        len(text),
        {"mv_mM": mv, "dv_mM": dv, "dntp_mM": dntp, "dna_nM": dna, "temp_c": temperature},
    )
    row["max_loop"] = loop
    row["provenance"] = _provenance("primer3-py", license="GPL-2.0")
    return row


def calculate_primer3_dimer(
    sequence: str,
    *,
    partner: str | None = None,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Self-dimer of one primer, or hetero-dimer with ``partner``, using primer3-py.

    The thermodynamic fields are None when primer3 finds no dimer structure.
    """

    primer3 = require("primer3", "bio-primer3")
    text = _canonical_dna(sequence, "sequence", maximum=200)
    other = None if partner is None else _canonical_dna(partner, "partner", maximum=200)
    mv, dv, dntp, dna = _thermo_conditions(mv_conc, dv_conc, dntp_conc, dna_conc)
    temperature = _number(temp_c, "temp_c", minimum=0.0, maximum=100.0)
    loop = _integer(max_loop, "max_loop", minimum=3, maximum=100)
    options = {
        "mv_conc": mv,
        "dv_conc": dv,
        "dntp_conc": dntp,
        "dna_conc": dna,
        "temp_c": temperature,
        "max_loop": loop,
    }
    try:
        if other is None:
            result = primer3.calc_homodimer(text, **options)
        else:
            result = primer3.calc_heterodimer(text, other, **options)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 could not compute the dimer: {exc}") from exc
    row = _thermo_row(
        "calculate_primer3_dimer",
        result,
        len(text),
        {"mv_mM": mv, "dv_mM": dv, "dntp_mM": dntp, "dna_nM": dna, "temp_c": temperature},
    )
    row["dimer_type"] = "homodimer" if other is None else "heterodimer"
    row["partner_length"] = None if other is None else len(other)
    row["max_loop"] = loop
    row["provenance"] = _provenance("primer3-py", license="GPL-2.0")
    return row


__all__ = [
    "build_feature_map",
    "build_genbank_feature_map",
    "build_sbol3_document",
    "calculate_primer3_dimer",
    "calculate_primer3_hairpin",
    "calculate_primer3_tm",
    "convert_snapgene_to_genbank",
    "design_assembly_primers",
    "design_pcr_primers",
    "evaluate_dnachisel_specs",
    "find_recombination_sites",
    "optimize_dnachisel_record",
    "read_dnacauldron_parts",
    "read_genbank_genome_tracks",
    "read_sbol3_sequences",
    "read_snapgene_features",
    "read_snapgene_sequence",
    "render_feature_map_svg",
    "render_genome_tracks_svg",
    "reverse_translate_protein",
    "validate_sbol3_document",
    "write_dnacauldron_constructs",
]

_register_public_verbose(globals(), tuple(__all__))
