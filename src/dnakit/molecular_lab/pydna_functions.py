"""pydna primer, PCR, overlap and assembly functions exposed through ``dnakit.molecular_lab``.

Each function imports pydna (and Biopython) only when it runs. A missing package
raises ``BackendUnavailableError`` with the ``dnakit[bio-assembly]`` install hint.
Positions are 0-based; overlap matches are ``(x_start, y_start, length)``. Every
function returns a JSON-safe dictionary with ``status``, ``operation``, the result
fields and ``provenance``. ``fetch_genbank_record()`` needs network access and is
covered by mocked tests only.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-assembly"
_MAX_PAIR = 20_000
_MAX_ASSEMBLY = 200_000
_MAX_FRAGMENTS = 20
_MAX_PRIMERS = 1000
_MAX_PRODUCTS = 200
_MAX_MATCHES = 1000
_MAX_TEXT_BYTES = 1_000_000
_ACCESSION = re.compile(r"[A-Za-z]{1,4}_?[0-9]{5,9}(\.[0-9]{1,3})?")
_OVERLAP_METHODS = {
    "gibson": "gibson_assembly",
    "in_fusion": "in_fusion_assembly",
    "fusion_pcr": "fusion_pcr_assembly",
    "in_vivo": "in_vivo_assembly",
}
_ENZYME_METHODS = {
    "golden_gate": "golden_gate_assembly",
    "restriction_ligation": "restriction_ligation_assembly",
}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _dna(value: object, name: str, *, maximum: int = _MAX_PAIR) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _fragments(value: object, *, maximum: int = _MAX_ASSEMBLY) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("fragments must be a list of DNA strings.")
    if not 2 <= len(value) <= _MAX_FRAGMENTS:
        raise ConfigurationError(f"fragments must contain between 2 and {_MAX_FRAGMENTS} items.")
    texts = [_dna(item, f"fragments[{index}]") for index, item in enumerate(value)]
    if sum(len(item) for item in texts) > maximum:
        raise ConfigurationError(f"fragments exceed {maximum} bases in total.")
    return texts


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("pydna")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "pydna", "version": installed, **notes}


def _matches(raw: Sequence[Sequence[int]]) -> tuple[list[dict[str, int]], int]:
    rows = sorted(
        ({"x_start": int(a), "y_start": int(b), "length": int(c)} for a, b, c in raw),
        key=lambda row: (row["x_start"], row["y_start"], row["length"]),
    )
    return rows[:_MAX_MATCHES], len(rows)


def simulate_pcr_amplicon(
    template: str,
    forward_primer: str,
    reverse_primer: str,
    *,
    min_annealing: int = 14,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate PCR of a linear DNA template with a primer pair using pydna.

    Each primer must anneal over at least ``min_annealing`` bases. The product
    runs from the 5' end of the forward primer to the 5' end of the reverse
    primer. A template without a product raises ``ConfigurationError``.
    """

    amplify = require("pydna.amplify", _EXTRA)
    primer_class = require("pydna.primer", _EXTRA).Primer
    dseq_class = require("pydna.dseqrecord", _EXTRA).Dseqrecord
    text = _dna(template, "template", maximum=_MAX_PAIR)
    forward = _dna(forward_primer, "forward_primer", maximum=200)
    reverse = _dna(reverse_primer, "reverse_primer", maximum=200)
    limit = _integer(min_annealing, "min_annealing", minimum=8, maximum=60)
    try:
        amplicon = amplify.pcr(
            primer_class(forward), primer_class(reverse), dseq_class(text), limit=limit
        )
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not amplify the template: {exc}") from exc
    product = str(amplicon.seq).upper()
    return {
        "status": "completed",
        "operation": "simulate_pcr_amplicon",
        "template_length": len(text),
        "forward_primer": forward,
        "reverse_primer": reverse,
        "min_annealing": limit,
        "product_sequence": product,
        "product_length": len(product),
        "provenance": _provenance(method="pydna.amplify.pcr"),
    }


def assemble_pydna_fragments(
    fragments: Sequence[str],
    *,
    method: str = "gibson",
    min_overlap: int = 25,
    enzymes: Sequence[str] = (),
    circular_only: bool = False,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate assembly products of DNA fragments with a pydna assembly algorithm.

    Overlap methods are ``gibson``, ``in_fusion``, ``fusion_pcr`` and ``in_vivo``;
    they use ``min_overlap``. Enzyme methods are ``golden_gate`` and
    ``restriction_ligation``; they need ``enzymes``, names from ``Bio.Restriction``
    such as ``"BsaI"``. Products are sorted by length and then sequence. At most
    ``max_products`` are listed; ``product_count`` gives the total.
    """

    assembly = require("pydna.assembly2", _EXTRA)
    dseq_class = require("pydna.dseqrecord", _EXTRA).Dseqrecord
    texts = _fragments(fragments)
    if method not in _OVERLAP_METHODS and method not in _ENZYME_METHODS:
        allowed = ", ".join(sorted([*_OVERLAP_METHODS, *_ENZYME_METHODS]))
        raise ConfigurationError(f"method must be one of: {allowed}.")
    circular = _flag(circular_only, "circular_only")
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    records = [dseq_class(text) for text in texts]
    try:
        if method in _OVERLAP_METHODS:
            if enzymes:
                raise ConfigurationError(
                    "enzymes applies only to golden_gate and restriction_ligation."
                )
            overlap = _integer(min_overlap, "min_overlap", minimum=10, maximum=200)
            function = getattr(assembly, _OVERLAP_METHODS[method])
            products = function(records, limit=overlap, circular_only=circular)
            chosen_enzymes: list[str] = []
            overlap_value: int | None = overlap
        else:
            if not enzymes or isinstance(enzymes, str):
                raise ConfigurationError(
                    "enzymes must list at least one enzyme name for this method."
                )
            chosen_enzymes = _enzyme_names(enzymes)
            classes = [_enzyme_class(name) for name in chosen_enzymes]
            function = getattr(assembly, _ENZYME_METHODS[method])
            products = function(records, enzymes=classes, circular_only=circular)
            overlap_value = None
    except ConfigurationError:
        raise
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not assemble the fragments: {exc}") from exc
    rows = sorted(
        (
            {
                "sequence": str(product.seq).upper(),
                "length": len(product),
                "circular": bool(product.circular),
            }
            for product in products
        ),
        key=lambda row: (row["length"], row["sequence"], row["circular"]),
    )
    return {
        "status": "completed",
        "operation": "assemble_pydna_fragments",
        "method": method,
        "fragment_count": len(texts),
        "min_overlap": overlap_value,
        "enzymes": chosen_enzymes,
        "circular_only": circular,
        "product_count": len(rows),
        "products": rows[:listed],
        "provenance": _provenance(
            method=f"pydna.assembly2.{_OVERLAP_METHODS.get(method) or _ENZYME_METHODS[method]}"
        ),
    }


def _enzyme_names(value: Sequence[str]) -> list[str]:
    if isinstance(value, bytes) or not isinstance(value, Sequence):
        raise ConfigurationError("enzymes must be a list of enzyme names.")
    if len(value) > 10:
        raise ConfigurationError("enzymes accepts at most 10 names.")
    names = [str(item).strip() for item in value]
    if any(not name or not name.isalnum() for name in names):
        raise ConfigurationError("Each enzyme name must be alphanumeric.")
    if len(set(names)) != len(names):
        raise ConfigurationError("enzymes must not repeat a name.")
    return names


def _enzyme_class(name: str) -> Any:
    restriction = require("Bio.Restriction", _EXTRA)
    cut_class = require("Bio.Restriction.Restriction", _EXTRA).AbstractCut
    candidate = getattr(restriction, name, None)
    if not (isinstance(candidate, type) and issubclass(candidate, cut_class)):
        raise ConfigurationError(f"Unknown restriction enzyme {name!r}.")
    return candidate


def find_common_substrings(
    sequence_x: str,
    sequence_y: str,
    *,
    min_length: int = 25,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find every common substring of at least ``min_length`` bases between two sequences.

    Each match is ``(x_start, y_start, length)`` with 0-based starts. Matches are
    sorted and capped at 1000 rows; ``match_count`` gives the total.
    """

    module = require("pydna.common_sub_strings", _EXTRA)
    first = _dna(sequence_x, "sequence_x")
    second = _dna(sequence_y, "sequence_y")
    limit = _integer(min_length, "min_length", minimum=1, maximum=_MAX_PAIR)
    try:
        raw = module.common_sub_strings(first, second, limit)
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not compare the sequences: {exc}") from exc
    rows, total = _matches(raw)
    return {
        "status": "completed",
        "operation": "find_common_substrings",
        "min_length": limit,
        "match_count": total,
        "matches": rows,
        "provenance": _provenance(method="pydna.common_sub_strings.common_sub_strings"),
    }


def find_terminal_overlaps(
    sequence_x: str,
    sequence_y: str,
    *,
    min_length: int = 15,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find overlaps where the end of ``sequence_x`` matches the start of ``sequence_y``.

    This is the terminal-overlap test used for Gibson-style joins. Each match is
    ``(x_start, y_start, length)``; the match ends at the end of ``sequence_x`` and
    starts at the beginning of ``sequence_y`` in pydna's overlap search.
    """

    module = require("pydna.common_sub_strings", _EXTRA)
    first = _dna(sequence_x, "sequence_x")
    second = _dna(sequence_y, "sequence_y")
    limit = _integer(min_length, "min_length", minimum=1, maximum=_MAX_PAIR)
    try:
        raw = module.terminal_overlap(first, second, limit)
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not compare the sequences: {exc}") from exc
    rows, total = _matches(raw)
    return {
        "status": "completed",
        "operation": "find_terminal_overlaps",
        "min_length": limit,
        "match_count": total,
        "matches": rows,
        "provenance": _provenance(method="pydna.common_sub_strings.terminal_overlap"),
    }


def primer_pair_hybridization(
    forward_primer: str,
    reverse_primer: str,
    *,
    minimal_annealing: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the duplex products that two primers can form with each other.

    The primers are 5'-to-3' oligonucleotides. ``minimal_annealing`` is the
    shortest shared complementary stretch. Primers that anneal with mismatches
    raise ``ConfigurationError``, as pydna does.
    """

    module = require("pydna.oligonucleotide_hybridization", _EXTRA)
    primer_class = require("pydna.primer", _EXTRA).Primer
    forward = _dna(forward_primer, "forward_primer", maximum=200)
    reverse = _dna(reverse_primer, "reverse_primer", maximum=200)
    minimum = _integer(minimal_annealing, "minimal_annealing", minimum=4, maximum=60)
    try:
        overhangs = sorted(
            int(item)
            for item in module.oligonucleotide_hybridization_overhangs(forward, reverse, minimum)
        )
        products = module.oligonucleotide_hybridization(
            primer_class(forward), primer_class(reverse), minimum
        )
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not hybridize the primers: {exc}") from exc
    rows = sorted(
        ({"sequence": str(item.seq).upper(), "length": len(item)} for item in products),
        key=lambda row: (row["length"], row["sequence"]),
    )
    return {
        "status": "completed",
        "operation": "primer_pair_hybridization",
        "minimal_annealing": minimum,
        "overhangs": overhangs,
        "product_count": len(rows),
        "products": rows,
        "provenance": _provenance(method="pydna.oligonucleotide_hybridization"),
    }


def primer_melting_temperature(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Compute the nearest-neighbour melting temperature of a primer with pydna's default table.

    The value is rounded to three decimals. Degenerate or non-ACGT input is
    rejected, as pydna's strict mode does.
    """

    tm_module = require("pydna.tm", _EXTRA)
    text = _dna(sequence, "sequence", maximum=200)
    if len(text) < 6:
        raise ConfigurationError("sequence must be at least 6 bases long.")
    try:
        value = float(tm_module.tm_default(text, check=True, strict=True))
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not compute Tm: {exc}") from exc
    return {
        "status": "completed",
        "operation": "primer_melting_temperature",
        "sequence": text,
        "length": len(text),
        "tm_celsius": round(value, 3),
        "provenance": _provenance(method="pydna.tm.tm_default", rounding="3 decimals"),
    }


def parse_primer_sequences(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read primer sequences from FASTA-like text or a local file.

    ``source`` is either a path to an existing file or the text itself. Names come
    from the FASTA header when present; the order of the input is kept.
    """

    parser = require("pydna.parsers", _EXTRA)
    if isinstance(source, Path) or (
        isinstance(source, str) and "\n" not in source and Path(source).expanduser().is_file()
    ):
        path = Path(source).expanduser().resolve()
        if path.stat().st_size > _MAX_TEXT_BYTES:
            raise ConfigurationError(f"Primer file exceeds {_MAX_TEXT_BYTES} bytes.")
        text = path.read_text(encoding="utf-8")
        origin = path.name
    elif isinstance(source, str) and source.strip():
        if len(source.encode("utf-8")) > _MAX_TEXT_BYTES:
            raise ConfigurationError(f"Primer text exceeds {_MAX_TEXT_BYTES} bytes.")
        text = source
        origin = "text"
    else:
        raise ConfigurationError("source must be a file path or non-empty primer text.")
    try:
        primers = parser.parse_primers(text)
    except (IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not read primers: {exc}") from exc
    if len(primers) > _MAX_PRIMERS:
        raise ConfigurationError(f"Input holds more than {_MAX_PRIMERS} primers.")
    rows = []
    for index, item in enumerate(primers):
        sequence = str(item.seq).upper()
        if set(sequence) - set("ACGT"):
            raise ConfigurationError(f"Primer {index + 1} must contain only A, C, G and T.")
        rows.append({"name": str(item.name), "sequence": sequence, "length": len(sequence)})
    return {
        "status": "completed",
        "operation": "parse_primer_sequences",
        "source": origin,
        "primer_count": len(rows),
        "primers": rows,
        "provenance": _provenance(method="pydna.parsers.parse_primers"),
    }


def fetch_genbank_record(accession: str, *, email: str, verbose: bool = True) -> dict[str, Any]:
    """Download a GenBank nucleotide record by accession through pydna's NCBI client.

    The accession looks like ``CS570233.1`` or ``NC_000913.3``. ``email`` is the
    contact address that NCBI requires; it is passed to the client and not stored
    in the result. The call needs network access.
    """

    module = require("pydna.genbank", _EXTRA)
    if not isinstance(accession, str) or not _ACCESSION.fullmatch(accession.strip()):
        raise ConfigurationError("accession must look like CS570233.1 or NC_000913.3.")
    if not isinstance(email, str) or "@" not in email or len(email) > 254 or " " in email:
        raise ConfigurationError("email must be a contact address for NCBI.")
    code = accession.strip()
    try:
        record = module.genbank(code, email=email)
    except (OSError, RuntimeError, ValueError) as exc:
        raise ConfigurationError(f"pydna could not fetch {code}: {exc}") from exc
    sequence = str(record.seq).upper()
    return {
        "status": "completed",
        "operation": "fetch_genbank_record",
        "accession": code,
        "name": str(getattr(record, "name", "")),
        "description": str(getattr(record, "description", "")),
        "topology": "circular" if bool(getattr(record, "circular", False)) else "linear",
        "length": len(sequence),
        "sequence": sequence,
        "provenance": _provenance(method="pydna.genbank.genbank", network="NCBI Entrez"),
    }


__all__ = [
    "assemble_pydna_fragments",
    "fetch_genbank_record",
    "find_common_substrings",
    "find_terminal_overlaps",
    "parse_primer_sequences",
    "primer_melting_temperature",
    "primer_pair_hybridization",
    "simulate_pcr_amplicon",
]

_register_public_verbose(globals(), tuple(__all__))
