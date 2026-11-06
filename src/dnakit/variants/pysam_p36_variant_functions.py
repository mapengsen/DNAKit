"""pysam variant and annotation-record operations (pack P36): writing VCF/BCF, header editing, record
translation, call comparison, IUPAC genotype codes and GTF/GFF3 attribute strings.

Each function calls the pysam class or method that carries the operation (``VariantFile``,
``VariantHeader``, ``VariantRecord``, ``VCF.compare_calls``, ``Pileup.encodeGenotype``, ``GTFProxy`` and
``GFF3Proxy``) with the caller's arguments. VCF positions in the inputs and outputs are 1-based, as in the
file format. Output files are written only to new paths; the suffix selects ``.vcf``, ``.vcf.gz`` or
``.bcf``. pysam is imported when a function runs; a missing package raises ``BackendUnavailableError``
with the ``bio-reads`` hint.

Not offered on purpose: ``GTFProxy.from_dict`` (it raises ``KeyError: attribute_dict2string`` in pysam
0.23.3; the working ``dict2attribute_string`` is used instead) and the legacy ``samtools pileup -c`` text
reader ``pysam.Pileup`` (that format no longer exists).
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-reads"
_MAX_ROWS = 1_000_000
_ERRORS = (ValueError, KeyError, OSError, TypeError, AttributeError, IndexError, LookupError)
_MODES = {".vcf": "w", ".gz": "wz", ".bcf": "wb"}


def _pysam() -> Any:
    return require("pysam", _EXTRA)


def _provenance(operation: str) -> dict[str, str]:
    return {
        "package": "pysam",
        "version": importlib.metadata.version("pysam"),
        "operation": operation,
        "coordinates": "1-based (VCF text)",
    }


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _rows(value: object, name: str, *, allow_empty: bool = True) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if (not allow_empty and not value) or len(value) > _MAX_ROWS:
        raise ConfigurationError(
            f"{name} must hold between {0 if allow_empty else 1} and {_MAX_ROWS} entries."
        )
    return list(value)


def _mapping(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping.")
    return value


def _new_output(value: object, name: str) -> tuple[Path, str]:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}")
    if path.name.endswith(".vcf.gz"):
        return path, "wz"
    if path.suffix == ".vcf":
        return path, "w"
    if path.suffix == ".bcf":
        return path, "wb"
    raise ConfigurationError(f"{name} must end with .vcf, .vcf.gz or .bcf.")


def _header_from_spec(
    pysam: Any,
    contigs: Sequence[Mapping[str, Any]],
    info_fields: Sequence[Mapping[str, Any]],
    format_fields: Sequence[Mapping[str, Any]],
    filters: Sequence[Mapping[str, Any]],
    meta: Sequence[Mapping[str, Any]],
    header_lines: Sequence[str],
    samples: Sequence[str],
) -> Any:
    header = pysam.VariantHeader()
    for item in meta:
        header.add_meta(
            _text(item.get("key"), "meta key"), _text(str(item.get("value")), "meta value")
        )
    for item in contigs:
        contig_id = _text(item.get("id"), "contig id")
        if item.get("length") is None:
            header.contigs.add(contig_id)
        else:
            header.contigs.add(contig_id, length=int(item["length"]))
    for kind, entries in (("INFO", info_fields), ("FORMAT", format_fields)):
        for item in entries:
            header.add_meta(
                kind,
                items=[
                    ("ID", _text(item.get("id"), f"{kind} id")),
                    ("Number", item.get("number", 1)),
                    ("Type", item.get("type", "String")),
                    ("Description", item.get("description", "")),
                ],
            )
    for item in filters:
        header.add_meta(
            "FILTER",
            items=[
                ("ID", _text(item.get("id"), "filter id")),
                ("Description", item.get("description", "")),
            ],
        )
    for line in header_lines:
        header.add_line(_text(line, "header line"))
    for name in samples:
        header.add_sample(_text(name, "sample name"))
    return header


def _as_value(value: Any) -> Any:
    return tuple(_as_value(v) for v in value) if isinstance(value, (list, tuple)) else value


def write_vcf_records(
    output_path: str | Path,
    records: Sequence[Mapping[str, Any]],
    *,
    contigs: Sequence[Mapping[str, Any]] = (),
    info_fields: Sequence[Mapping[str, Any]] = (),
    format_fields: Sequence[Mapping[str, Any]] = (),
    filters: Sequence[Mapping[str, Any]] = (),
    meta: Sequence[Mapping[str, Any]] = (),
    header_lines: Sequence[str] = (),
    samples: Sequence[str] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a VCF header and write VCF, bgzipped VCF or BCF records (``VariantFile.write``).

    The header is made with ``VariantHeader.add_meta`` / ``add_line`` / ``add_sample``: ``contigs`` are
    ``{"id", "length"}``; ``info_fields`` and ``format_fields`` are ``{"id", "number", "type",
    "description"}``; ``filters`` are ``{"id", "description"}``; ``meta`` are ``{"key", "value"}`` lines such
    as ``fileDate``; ``header_lines`` are raw ``##...`` lines; ``samples`` are the sample names.
    Each record is ``{"contig", "pos", "ref", "alts", "id", "qual", "filter", "info", "samples"}`` with a
    1-based ``pos``; ``samples`` maps a sample name to its FORMAT values (``GT`` as a list such as
    ``[0, 1]``; ``"phased": true`` marks a phased genotype). The output suffix selects the format.
    """

    target, mode = _new_output(output_path, "output_path")
    rows = _rows(records, "records")
    pysam = _pysam()
    try:
        header = _header_from_spec(
            pysam,
            [_mapping(c, "contig") for c in _rows(contigs, "contigs")],
            [_mapping(c, "info field") for c in _rows(info_fields, "info_fields")],
            [_mapping(c, "format field") for c in _rows(format_fields, "format_fields")],
            [_mapping(c, "filter") for c in _rows(filters, "filters")],
            [_mapping(c, "meta") for c in _rows(meta, "meta")],
            [str(c) for c in _rows(header_lines, "header_lines")],
            [str(c) for c in _rows(samples, "samples")],
        )
        with pysam.VariantFile(str(target), mode, header=header) as out:
            for index, item in enumerate(rows):
                row = _mapping(item, f"records[{index}]")
                pos = int(row["pos"])
                ref = _text(row.get("ref"), "ref")
                alts = tuple(str(a) for a in row.get("alts", ()))
                record = out.new_record(
                    contig=_text(row.get("contig"), "contig"),
                    start=pos - 1,
                    alleles=(ref, *alts),
                    id=row.get("id"),
                    qual=row.get("qual"),
                    filter=row.get("filter") or None,
                    info={k: _as_value(v) for k, v in dict(row.get("info", {})).items()},
                )
                for sample, values in dict(row.get("samples", {})).items():
                    for key, value in dict(values).items():
                        if key == "phased":
                            record.samples[sample].phased = bool(value)
                        else:
                            record.samples[sample][key] = _as_value(value)
                out.write(record)
    except _ERRORS as exc:
        target.unlink(missing_ok=True)
        raise BackendExecutionError(f"pysam could not write the VCF: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_vcf_records",
        "output": str(target),
        "format": {"w": "VCF", "wz": "VCF (bgzip)", "wb": "BCF"}[mode],
        "record_count": len(rows),
        "sample_count": len(samples),
        "header_text": str(header),
        "provenance": _provenance("VariantFile.write"),
    }


def merge_vcf_headers(
    first_vcf: str | Path,
    second_vcf: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the header of a second VCF into a copy of the first header and write it (``VariantHeader.merge``).

    The output is a header-only VCF. Lines present in the second header but not in the first are added;
    pysam does not merge sample columns, so the samples of the first header are kept.
    """

    pysam = _pysam()
    first, second = input_path(first_vcf), input_path(second_vcf)
    target, mode = _new_output(output_path, "output_path")
    try:
        with pysam.VariantFile(str(first)) as one, pysam.VariantFile(str(second)) as two:
            merged = one.header.copy()
            merged.merge(two.header)
            before = {str(r) for r in one.header.records}
            added = [str(r).strip() for r in merged.records if str(r) not in before]
        with pysam.VariantFile(str(target), mode, header=merged):
            pass
    except _ERRORS as exc:
        target.unlink(missing_ok=True)
        raise BackendExecutionError(f"pysam header merge failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "merge_vcf_headers",
        "output": str(target),
        "added_header_lines": added,
        "samples": list(merged.samples),
        "header_text": str(merged),
        "provenance": _provenance("VariantHeader.merge"),
    }


def translate_vcf_records(
    source_vcf: str | Path,
    header_vcf: str | Path,
    output_path: str | Path,
    *,
    max_records: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rewrite the records of one VCF under the header of another VCF (``VariantRecord.translate``).

    Contigs, INFO, FORMAT and FILTER definitions are looked up in the target header; a record that uses a
    definition missing there fails with an error that names it. Sample columns of the target header that
    the source lacks stay missing.
    """

    pysam = _pysam()
    source, target_header = input_path(source_vcf), input_path(header_vcf)
    target, mode = _new_output(output_path, "output_path")
    limit = int(max_records)
    if not 1 <= limit <= _MAX_ROWS:
        raise ConfigurationError(f"max_records must be in [1, {_MAX_ROWS}].")
    count = 0
    try:
        with pysam.VariantFile(str(source)) as src, pysam.VariantFile(str(target_header)) as dst:
            header = dst.header
            with pysam.VariantFile(str(target), mode, header=header) as out:
                for record in src:
                    if count >= limit:
                        raise ConfigurationError(
                            f"source holds more than max_records={limit} records."
                        )
                    copy = record.copy()
                    copy.translate(header)
                    out.write(copy)
                    count += 1
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    except _ERRORS as exc:
        target.unlink(missing_ok=True)
        raise BackendExecutionError(f"pysam translate failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "translate_vcf_records",
        "output": str(target),
        "record_count": count,
        "provenance": _provenance("VariantRecord.translate"),
    }


def compare_variant_calls(
    pos1: int,
    ref1: str,
    alt1: str,
    pos2: int,
    ref2: str,
    alt2: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decide whether two variant calls are the same event (``pysam.VCF.compare_calls``).

    The positions must be equal; then the shared rightmost bases of REF and ALT are trimmed from each call
    and only the remaining ALT strings are compared (pysam does not compare the trimmed REF strings).
    ``ATT>AT`` and ``AT>A`` at the same position are the same deletion of one T; ``A>C`` and ``A>G`` are not.
    """

    pysam = _pysam()
    values = (
        int(pos1),
        _text(ref1, "ref1"),
        _text(alt1, "alt1"),
        int(pos2),
        _text(ref2, "ref2"),
        _text(alt2, "alt2"),
    )
    try:
        same = bool(pysam.VCF().compare_calls(*values))
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam compare_calls failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "compare_variant_calls",
        "call1": {"pos": values[0], "ref": values[1], "alt": values[2]},
        "call2": {"pos": values[3], "ref": values[4], "alt": values[5]},
        "equal": same,
        "provenance": _provenance("VCF.compare_calls"),
    }


def convert_iupac_genotype_codes(
    values: Sequence[str],
    *,
    direction: str = "encode",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert diploid genotypes to one-letter IUPAC codes and back (``Pileup.encodeGenotype`` / ``decodeGenotype``).

    ``direction="encode"`` maps two-base genotypes such as ``AG`` to ``r`` (lower case when the first base
    sorts before the second, upper case otherwise; homozygous ``GG`` gives ``G``). ``direction="decode"``
    maps codes such as ``m`` or ``M`` back to the two-base genotype. Unknown values are reported in
    ``errors`` and give ``None``.
    """

    if direction not in ("encode", "decode"):
        raise ConfigurationError("direction must be encode or decode.")
    pysam = _pysam()
    from importlib import import_module

    pileup = import_module("pysam.Pileup")
    del pysam
    rows = [_text(v, "values item") for v in _rows(values, "values", allow_empty=False)]
    function = pileup.encodeGenotype if direction == "encode" else pileup.decodeGenotype
    out: dict[str, str | None] = {}
    errors: list[str] = []
    for value in rows:
        try:
            out[value] = function(value)
        except KeyError:
            out[value] = None
            errors.append(value)
    return {
        "status": "completed",
        "operation": "convert_iupac_genotype_codes",
        "direction": direction,
        "results": out,
        "errors": errors,
        "provenance": _provenance(f"Pileup.{direction}Genotype"),
    }


def convert_gff_attributes(
    *,
    attribute_string: str | None = None,
    attributes: Mapping[str, Any] | None = None,
    flavor: str = "gtf",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the attribute column of a GTF or GFF3 line between its text form and a dictionary.

    Give ``attribute_string`` (for example ``gene_id "A"; gene_name "B";`` for ``flavor="gtf"`` or
    ``ID=a;Name=b`` for ``flavor="gff3"``) to get the parsed dictionary (``GTFProxy.attribute_string2dict`` /
    ``GFF3Proxy.attribute_string2dict``), or ``attributes`` to get the formatted string
    (``dict2attribute_string``). GFF3 values stay percent-encoded as in the file; the string form does not
    add escapes.
    """

    if flavor not in ("gtf", "gff3"):
        raise ConfigurationError("flavor must be gtf or gff3.")
    if (attribute_string is None) == (attributes is None):
        raise ConfigurationError("Give exactly one of attribute_string and attributes.")
    pysam = _pysam()
    line = b"chr1\tsrc\tfeature\t1\t2\t.\t+\t.\t."
    try:
        parser = pysam.asGTF() if flavor == "gtf" else pysam.asGFF3()
        proxy = parser(line, len(line))
        if attribute_string is not None:
            text = _text(attribute_string, "attribute_string")
            parsed = {str(k): v for k, v in proxy.attribute_string2dict(text).items()}
            result: dict[str, Any] = {"attributes": parsed}
            call = "attribute_string2dict"
        else:
            mapping = _mapping(attributes, "attributes")
            if not mapping:
                raise ConfigurationError("attributes must not be empty.")
            result = {"attribute_string": proxy.dict2attribute_string(dict(mapping))}
            call = "dict2attribute_string"
    except ConfigurationError:
        raise
    except _ERRORS as exc:
        raise BackendExecutionError(f"pysam attribute conversion failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "convert_gff_attributes",
        "flavor": flavor,
        **result,
        "provenance": _provenance(f"{'GTFProxy' if flavor == 'gtf' else 'GFF3Proxy'}.{call}"),
    }


__all__ = [
    "compare_variant_calls",
    "convert_gff_attributes",
    "convert_iupac_genotype_codes",
    "merge_vcf_headers",
    "translate_vcf_records",
    "write_vcf_records",
]
_register_public_verbose(globals())
