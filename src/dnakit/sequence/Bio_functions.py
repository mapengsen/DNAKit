"""Biopython-backed DNA and gene-level functions for ``dnakit.sequence``.

Each public function imports Biopython only when it is called, so the module
imports without the optional ``bio-formats`` extra. Inputs are explicit
strings, mappings or local files. Network services, external binaries and
structure or protein-only utilities are not wrapped here.

Coordinates in the results are 1-based and inclusive unless a field name says
otherwise. Results are plain JSON-serializable dictionaries.
"""

from __future__ import annotations

import importlib.metadata
import math
import random
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_UNAMBIGUOUS_DNA = frozenset("ACGT")
_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY*X-")
_MAX_SEQUENCE = 1_000_000
_MAX_ROWS = 500
_MAX_RECORDS = 100_000
_MAX_REPLICATES = 1_000
_MAX_GENETIC_CODE = 33
_DNA_MODELS = ("NG86", "LWL85")
_QUALITY_DIRECTIONS = ("solexa_to_phred", "phred_to_solexa")
_TREE_METHODS = ("nj", "upgma")
_DISTANCE_MODELS = ("identity", "blastn", "trans")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _sequence(value: object, name: str, *, alphabet: frozenset[str], maximum: int) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds the limit of {maximum} symbols.")
    invalid = sorted(set(text) - alphabet)
    if invalid:
        raise ConfigurationError(f"{name} contains unsupported symbols {invalid}.")
    return text


def _named_sequences(
    value: object,
    name: str,
    *,
    alphabet: frozenset[str],
    maximum: int,
    minimum_rows: int = 2,
    max_rows: int = _MAX_ROWS,
    equal_length: bool = True,
) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping from identifier to sequence.")
    if not minimum_rows <= len(value) <= max_rows:
        raise ConfigurationError(
            f"{name} must contain between {minimum_rows} and {max_rows} sequences."
        )
    items: dict[str, str] = {}
    for key in sorted(value, key=str):
        identifier = str(key).strip()
        if not identifier:
            raise ConfigurationError(f"{name} identifiers must be non-empty.")
        if identifier in items:
            raise ConfigurationError(f"{name} identifiers must be unique.")
        items[identifier] = _sequence(
            value[key], f"{name}[{identifier}]", alphabet=alphabet, maximum=maximum
        )
    lengths = {len(text) for text in items.values()}
    if equal_length and len(lengths) != 1:
        raise ConfigurationError(f"all rows in {name} must have equal length.")
    return items


def _codon_table(genetic_code: object) -> tuple[int, Any]:
    table_id = _integer(genetic_code, "genetic_code", minimum=1, maximum=_MAX_GENETIC_CODE)
    codon_table = require("Bio.Data.CodonTable", _EXTRA)
    table = codon_table.unambiguous_dna_by_id.get(table_id)
    if table is None:
        raise ConfigurationError(f"genetic_code {table_id} is not available in Biopython.")
    return table_id, table


def _seq_records(items: Mapping[str, str], bio_seq: Any, bio_record: Any) -> list[Any]:
    return [bio_record(bio_seq(text), id=identifier) for identifier, text in items.items()]


def _numeric_alignment(sequences: Mapping[str, str]) -> Any:
    """Build an unambiguous Biopython Alignment with ungapped identical coordinates."""

    np = require("numpy", "bio-analysis")
    align = require("Bio.Align", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    length = len(next(iter(sequences.values())))
    records = _seq_records(sequences, seq_module.Seq, record_module.SeqRecord)
    coordinates = np.array([[0, length] for _ in records])
    return align.Alignment(records, coordinates)


def _matrix_from_lower(values: Sequence[Sequence[float]], size: int) -> list[list[float | None]]:
    full: list[list[float | None]] = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row):
            value = values[row][column]
            cell = None if value is None or value < 0 else float(value)
            full[row][column] = cell
            full[column][row] = cell
    return full


def gc_content_by_codon_position(
    sequence: str,
    *,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return total and codon-position G+C percentages with Biopython's GC123."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    text = _sequence(sequence, "sequence", alphabet=_IUPAC_DNA, maximum=limit)
    seq_utils = require("Bio.SeqUtils", _EXTRA)
    total, first, second, third = seq_utils.GC123(text)
    return {
        "status": "completed",
        "operation": "gc_content_by_codon_position",
        "sequence_length": len(text),
        "gc_percent": {
            "total": float(total),
            "first": float(first),
            "second": float(second),
            "third": float(third),
        },
        "provenance": _provenance(
            parser="Bio.SeqUtils.GC123",
            positions="codon positions 1, 2 and 3 count from the first base",
            ambiguity="IUPAC symbols other than G and C count only toward the length",
        ),
    }


def gc_skew_windows(
    sequence: str,
    *,
    window: int = 100,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return (G-C)/(G+C) for consecutive non-overlapping windows."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    text = _sequence(sequence, "sequence", alphabet=_IUPAC_DNA, maximum=limit)
    size = _integer(window, "window", minimum=1, maximum=len(text))
    seq_utils = require("Bio.SeqUtils", _EXTRA)
    skews = seq_utils.GC_skew(text, window=size)
    windows = [
        {
            "start": index * size + 1,
            "end": min((index + 1) * size, len(text)),
            "skew": float(skew),
        }
        for index, skew in enumerate(skews)
    ]
    return {
        "status": "completed",
        "operation": "gc_skew_windows",
        "sequence_length": len(text),
        "window": size,
        "window_count": len(windows),
        "windows": windows,
        "provenance": _provenance(
            parser="Bio.SeqUtils.GC_skew",
            windowing="consecutive non-overlapping windows; the last window may be shorter",
            empty_skew="windows without G or C report skew 0.0",
        ),
    }


def local_composition_complexity(
    sequence: str,
    *,
    window: int = 20,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the local composition complexity (LCC) of every sliding window."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    text = _sequence(sequence, "sequence", alphabet=_UNAMBIGUOUS_DNA, maximum=limit)
    size = _integer(window, "window", minimum=1, maximum=len(text))
    lcc = require("Bio.SeqUtils.lcc", _EXTRA)
    values = [float(value) for value in lcc.lcc_mult(text, size)]
    return {
        "status": "completed",
        "operation": "local_composition_complexity",
        "sequence_length": len(text),
        "window": size,
        "value_count": len(values),
        "values": values,
        "minimum": min(values),
        "maximum": max(values),
        "mean": math.fsum(values) / len(values),
        "provenance": _provenance(
            parser="Bio.SeqUtils.lcc.lcc_mult",
            indexing="values[i] describes the window that starts at 1-based position i + 1",
        ),
    }


def ambiguous_codon_amino_acids(
    codon: str,
    *,
    genetic_code: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return every amino acid encoded by one IUPAC-ambiguous codon."""

    text = _sequence(codon, "codon", alphabet=_IUPAC_DNA, maximum=3)
    if len(text) != 3:
        raise ConfigurationError("codon must contain exactly three nucleotide symbols.")
    table_id, table = _codon_table(genetic_code)
    codon_table = require("Bio.Data.CodonTable", _EXTRA)
    iupac = require("Bio.Data.IUPACData", _EXTRA)
    stop_only = False
    try:
        amino_acids = sorted(
            codon_table.list_possible_proteins(
                text, table.forward_table, iupac.ambiguous_dna_values
            )
        )
    except KeyError:
        stop_only = True
        amino_acids = []
    except codon_table.TranslationError as exc:
        raise ConfigurationError(str(exc)) from exc
    return {
        "status": "completed",
        "operation": "ambiguous_codon_amino_acids",
        "codon": text,
        "genetic_code": table_id,
        "amino_acids": amino_acids,
        "amino_acid_count": len(amino_acids),
        "stop_only": stop_only,
        "provenance": _provenance(
            parser="Bio.Data.CodonTable.list_possible_proteins",
            table=f"NCBI genetic code {table_id}",
        ),
    }


def expand_ambiguous_codons(
    codons: Sequence[str],
    *,
    max_codons: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add IUPAC-ambiguous codons that are fully covered by an input codon list."""

    limit = _integer(max_codons, "max_codons", minimum=1, maximum=10_000)
    if isinstance(codons, str) or not isinstance(codons, Sequence) or not codons:
        raise ConfigurationError("codons must be a non-empty sequence of three-letter strings.")
    if len(codons) > limit:
        raise ConfigurationError(f"codons exceeds the limit of {limit} entries.")
    cleaned = sorted(
        {_sequence(codon, "codon", alphabet=_UNAMBIGUOUS_DNA, maximum=3) for codon in codons}
    )
    if any(len(codon) != 3 for codon in cleaned):
        raise ConfigurationError("every codon must contain exactly three nucleotide symbols.")
    codon_table = require("Bio.Data.CodonTable", _EXTRA)
    iupac = require("Bio.Data.IUPACData", _EXTRA)
    expanded = sorted(set(codon_table.list_ambiguous_codons(cleaned, iupac.ambiguous_dna_values)))
    added = [codon for codon in expanded if codon not in cleaned]
    return {
        "status": "completed",
        "operation": "expand_ambiguous_codons",
        "input_codons": cleaned,
        "added_codons": added,
        "expanded_codons": expanded,
        "provenance": _provenance(
            parser="Bio.Data.CodonTable.list_ambiguous_codons",
            rule="an ambiguous codon is added only when every concrete codon it covers is present",
        ),
    }


def sequence_checksum(
    sequence: str,
    *,
    method: str = "seguid",
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a SEGUID, CRC32, CRC64 or GCG checksum of one sequence."""

    selected = _choice(method, ("seguid", "crc32", "crc64", "gcg"), "method")
    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    text = _sequence(sequence, "sequence", alphabet=_IUPAC_DNA, maximum=limit)
    checksum = require("Bio.SeqUtils.CheckSum", _EXTRA)
    value = getattr(checksum, selected)(text)
    return {
        "status": "completed",
        "operation": "sequence_checksum",
        "method": selected,
        "sequence_length": len(text),
        "value": value,
        "provenance": _provenance(
            parser=f"Bio.SeqUtils.CheckSum.{selected}",
            input="whitespace removed and letters uppercased before the checksum",
        ),
    }


def mcdonald_kreitman_test(
    sequences: Mapping[str, str],
    species: Mapping[str, str],
    *,
    genetic_code: int = 1,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the McDonald-Kreitman neutrality test on coding sequences from two or more species."""

    limit = _integer(max_length, "max_length", minimum=3, maximum=_MAX_SEQUENCE)
    items = _named_sequences(
        sequences,
        "sequences",
        alphabet=_UNAMBIGUOUS_DNA,
        maximum=limit,
        max_rows=_MAX_ROWS,
    )
    length = len(next(iter(items.values())))
    if length % 3:
        raise ConfigurationError("coding sequences must have a length that is a multiple of 3.")
    if not isinstance(species, Mapping) or set(map(str, species)) != set(items):
        raise ConfigurationError("species must map every sequence identifier to a species label.")
    labels = {identifier: str(species[identifier]).strip() for identifier in items}
    if any(not label for label in labels.values()):
        raise ConfigurationError("species labels must be non-empty.")
    if len(set(labels.values())) < 2:
        raise ConfigurationError("the McDonald-Kreitman test needs at least two species.")
    table_id, table = _codon_table(genetic_code)
    analysis = require("Bio.Align.analysis", _EXTRA)
    alignment = _numeric_alignment(items)
    try:
        p_value = analysis.mktest(
            alignment, species=[labels[key] for key in items], codon_table=table
        )
    except ZeroDivisionError as exc:
        raise ConfigurationError(
            "the McDonald-Kreitman 2x2 table has an empty cell; "
            "the test is undefined for these counts."
        ) from exc
    return {
        "status": "completed",
        "operation": "mcdonald_kreitman_test",
        "genetic_code": table_id,
        "sequence_count": len(items),
        "codon_count": length // 3,
        "species": sorted(set(labels.values())),
        "p_value": float(p_value),
        "provenance": _provenance(
            parser="Bio.Align.analysis.mktest",
            counting="changes are counted per codon; fixed and polymorphic sites are separated",
        ),
    }


def pairwise_dn_ds_matrix(
    sequences: Mapping[str, str],
    *,
    method: str = "NG86",
    genetic_code: int = 1,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return pairwise dN and dS matrices for codon-aligned coding sequences."""

    selected = _choice(method, _DNA_MODELS, "method")
    limit = _integer(max_length, "max_length", minimum=3, maximum=_MAX_SEQUENCE)
    items = _named_sequences(
        sequences,
        "sequences",
        alphabet=_UNAMBIGUOUS_DNA,
        maximum=limit,
        max_rows=_MAX_ROWS,
    )
    if len(next(iter(items.values()))) % 3:
        raise ConfigurationError("coding sequences must have a length that is a multiple of 3.")
    table_id, table = _codon_table(genetic_code)
    analysis = require("Bio.Align.analysis", _EXTRA)
    alignment = _numeric_alignment(items)
    dn_matrix, ds_matrix = analysis.calculate_dn_ds_matrix(
        alignment, method=selected, codon_table=table
    )
    ids = list(items)
    size = len(ids)
    dn_full = _matrix_from_lower(dn_matrix.matrix, size)
    ds_full = _matrix_from_lower(ds_matrix.matrix, size)
    saturated = sum(
        1
        for row in range(size)
        for column in range(row)
        if dn_full[row][column] is None or ds_full[row][column] is None
    )
    return {
        "status": "completed",
        "operation": "pairwise_dn_ds_matrix",
        "method": selected,
        "genetic_code": table_id,
        "ids": ids,
        "dn": dn_full,
        "ds": ds_full,
        "saturated_pair_count": saturated,
        "provenance": _provenance(
            parser=f"Bio.Align.analysis.calculate_dn_ds_matrix ({selected})",
            saturation="null marks a -1 value, which Biopython reports when changes saturate",
        ),
    }


def build_codon_alignment(
    protein_alignment: Mapping[str, str],
    nucleotide_sequences: Mapping[str, str],
    *,
    genetic_code: int = 1,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Thread ungapped DNA codons through a protein alignment."""

    limit = _integer(max_length, "max_length", minimum=3, maximum=_MAX_SEQUENCE)
    proteins = _named_sequences(
        protein_alignment,
        "protein_alignment",
        alphabet=_AMINO_ACIDS,
        maximum=limit,
        max_rows=_MAX_ROWS,
    )
    nucleotides = _named_sequences(
        nucleotide_sequences,
        "nucleotide_sequences",
        alphabet=_UNAMBIGUOUS_DNA,
        maximum=limit,
        max_rows=_MAX_ROWS,
        equal_length=False,
    )
    if set(proteins) != set(nucleotides):
        raise ConfigurationError(
            "protein_alignment and nucleotide_sequences need the same identifiers."
        )
    table_id, table = _codon_table(genetic_code)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    codonalign = require("Bio.codonalign", _EXTRA)
    protein_rows = [
        record_module.SeqRecord(seq_module.Seq(text), id=identifier)
        for identifier, text in proteins.items()
    ]
    nucleotide_rows = [
        record_module.SeqRecord(seq_module.Seq(nucleotides[identifier]), id=identifier)
        for identifier in proteins
    ]
    try:
        multiple = require("Bio.Align", _EXTRA).MultipleSeqAlignment(protein_rows)
        result = codonalign.build(multiple, nucleotide_rows, codon_table=table)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"Biopython could not build the codon alignment: {exc}") from exc
    rows = {str(record.id): str(record.seq) for record in result}
    ordered = {identifier: rows[identifier] for identifier in proteins}
    return {
        "status": "completed",
        "operation": "build_codon_alignment",
        "genetic_code": table_id,
        "ids": list(proteins),
        "alignment_length": result.get_alignment_length(),
        "sequences": ordered,
        "provenance": _provenance(
            parser="Bio.codonalign.build",
            gaps="protein gaps become three-character codon gaps; stop codons are not added",
        ),
    }


def parse_primersearch_output(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read EMBOSS primersearch output into amplicon sequences grouped by primer name."""

    source = input_path(path)
    emboss = require("Bio.Emboss.PrimerSearch", _EXTRA)
    with source.open(encoding="utf-8") as handle:
        record = emboss.read(handle)
    primers: dict[str, list[dict[str, Any]]] = {}
    for name in sorted(record.amplifiers):
        amplicons = []
        for amplifier in record.amplifiers[name]:
            lines = amplifier.hit_info.splitlines()
            sequence = lines[0].strip() if lines else ""
            details = [line.strip() for line in lines[1:] if line.strip()]
            amplicons.append(
                {"sequence": sequence, "length": int(amplifier.length), "details": details}
            )
        primers[name] = amplicons
    return {
        "status": "completed",
        "operation": "parse_primersearch_output",
        "source": str(source),
        "primer_count": len(primers),
        "amplicon_count": sum(len(items) for items in primers.values()),
        "primers": primers,
        "provenance": _provenance(
            parser="Bio.Emboss.PrimerSearch.read",
            scope="parses a saved primersearch report; does not run EMBOSS",
        ),
    }


def read_ace_assembly(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Parse an ACE assembly file into contigs, their reads and read placements."""

    source = input_path(path)
    ace = require("Bio.Sequencing.Ace", _EXTRA)
    with source.open(encoding="utf-8") as handle:
        try:
            record = ace.read(handle)
        except ValueError as exc:
            raise ConfigurationError(f"{source.name} is not a valid ACE file: {exc}") from exc
    contigs = []
    for contig in record.contigs:
        contigs.append(
            {
                "name": str(contig.name),
                "length": len(contig.sequence),
                "sequence": str(contig.sequence),
                "reads": [
                    {
                        "name": str(item.rd.name),
                        "sequence": str(item.rd.sequence),
                        "quality_clipping": [
                            int(item.qa.qual_clipping_start),
                            int(item.qa.qual_clipping_end),
                        ],
                        "alignment_clipping": [
                            int(item.qa.align_clipping_start),
                            int(item.qa.align_clipping_end),
                        ],
                    }
                    for item in contig.reads
                ],
                "placements": [
                    {
                        "read": str(item.name),
                        "strand": str(item.coru),
                        "padded_start": int(item.padded_start),
                    }
                    for item in contig.af
                ],
            }
        )
    return {
        "status": "completed",
        "operation": "read_ace_assembly",
        "source": str(source),
        "contig_count": len(contigs),
        "read_count": sum(len(contig["reads"]) for contig in contigs),
        "contigs": contigs,
        "provenance": _provenance(
            parser="Bio.Sequencing.Ace.read",
            coordinates="clipping and padded positions are 1-based as written in the ACE file",
        ),
    }


def read_phd_sequence(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read one PHRED base-call (PHD) record with its per-base qualities and peaks."""

    source = input_path(path)
    phd = require("Bio.Sequencing.Phd", _EXTRA)
    try:
        record = phd.read(str(source))
    except ValueError as exc:
        raise ConfigurationError(f"{source.name} is not a valid PHD record: {exc}") from exc
    bases = [str(site[0]) for site in record.sites]
    qualities = [int(site[1]) for site in record.sites]
    peaks = [
        int(site[2]) if len(site) > 2 and site[2] is not None else None for site in record.sites
    ]
    comments = {key: value for key, value in sorted(record.comments.items()) if value is not None}
    return {
        "status": "completed",
        "operation": "read_phd_sequence",
        "source": str(source),
        "file_name": str(record.file_name),
        "sequence": "".join(bases),
        "length": len(bases),
        "qualities": qualities,
        "peak_locations": peaks,
        "comments": comments,
        "provenance": _provenance(
            parser="Bio.Sequencing.Phd.read",
            positions="peak locations are 1-based positions in the trace, when the file gives them",
        ),
    }


def read_genepop_genotypes(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a GenePop genotype file into loci, populations and two-allele genotypes."""

    source = input_path(path)
    genepop = require("Bio.PopGen.GenePop", _EXTRA)
    with source.open(encoding="utf-8") as handle:
        try:
            record = genepop.read(handle)
        except (ValueError, IndexError) as exc:
            raise ConfigurationError(f"{source.name} is not a valid GenePop file: {exc}") from exc
    populations = []
    for index, population in enumerate(record.populations, start=1):
        individuals = []
        for name, genotypes in population:
            individuals.append(
                {
                    "name": str(name).strip(),
                    "genotypes": [[allele_a, allele_b] for allele_a, allele_b in genotypes],
                }
            )
        populations.append(
            {"index": index, "individual_count": len(individuals), "individuals": individuals}
        )
    return {
        "status": "completed",
        "operation": "read_genepop_genotypes",
        "source": str(source),
        "comment": str(record.comment_line).strip(),
        "marker_length": int(record.marker_len),
        "loci": [str(locus) for locus in record.loci_list],
        "population_count": len(populations),
        "individual_count": sum(item["individual_count"] for item in populations),
        "populations": populations,
        "provenance": _provenance(
            parser="Bio.PopGen.GenePop.read",
            missing_data="a genotype with 0000 is returned as [null, null]",
        ),
    }


def read_gaf_annotations(
    path: str | Path,
    *,
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse GAF 1.0 or 2.x gene-annotation records with GO term and aspect counts."""

    source = input_path(path)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    goa = require("Bio.UniProt.GOA", _EXTRA)
    records: list[dict[str, Any]] = []
    go_counts: dict[str, int] = {}
    aspect_counts: dict[str, int] = {}
    with source.open(encoding="utf-8") as handle:
        try:
            for item in goa.gafiterator(handle):
                if len(records) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
                go_id = str(item.get("GO_ID", ""))
                aspect = str(item.get("Aspect", ""))
                go_counts[go_id] = go_counts.get(go_id, 0) + 1
                aspect_counts[aspect] = aspect_counts.get(aspect, 0) + 1
                records.append(
                    {
                        "db": str(item.get("DB", "")),
                        "db_object_id": str(item.get("DB_Object_ID", "")),
                        "db_object_symbol": str(item.get("DB_Object_Symbol", "")),
                        "db_object_name": str(item.get("DB_Object_Name", "")),
                        "db_object_type": str(item.get("DB_Object_Type", "")),
                        "qualifiers": [str(value) for value in item.get("Qualifier", []) if value],
                        "go_id": go_id,
                        "references": [
                            str(value) for value in item.get("DB:Reference", []) if value
                        ],
                        "evidence": str(item.get("Evidence", "")),
                        "with": [str(value) for value in item.get("With", []) if value],
                        "aspect": aspect,
                        "synonyms": [str(value) for value in item.get("Synonym", []) if value],
                        "taxon_ids": [str(value) for value in item.get("Taxon_ID", []) if value],
                        "date": str(item.get("Date", "")),
                        "assigned_by": str(item.get("Assigned_By", "")),
                    }
                )
        except ValueError as exc:
            raise ConfigurationError(f"{source.name} is not a valid GAF file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_gaf_annotations",
        "source": str(source),
        "record_count": len(records),
        "go_term_counts": dict(sorted(go_counts.items())),
        "aspect_counts": dict(sorted(aspect_counts.items())),
        "records": records,
        "provenance": _provenance(
            parser="Bio.UniProt.GOA.gafiterator",
            scope="annotation records are read from a local file only",
        ),
    }


def parse_kegg_gene_entries(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Parse KEGG GENE flat-file entries: names, orthology, organism, position and links."""

    source = input_path(path)
    kegg_gene = require("Bio.KEGG.Gene", _EXTRA)
    entries: list[dict[str, Any]] = []
    with source.open(encoding="utf-8") as handle:
        try:
            for item in kegg_gene.parse(handle):
                organism = getattr(item, "organism", None)
                entries.append(
                    {
                        "entry": str(item.entry),
                        "names": [str(name) for name in item.name],
                        "definition": str(item.definition),
                        "orthology": [
                            {"id": str(key), "name": str(name)} for key, name in item.orthology
                        ],
                        "organism": (
                            {"code": str(organism[0]), "name": str(organism[1])}
                            if organism
                            else None
                        ),
                        "position": str(item.position) if item.position else None,
                        "dblinks": [
                            {
                                "database": str(database),
                                "identifiers": [str(value) for value in values],
                            }
                            for database, values in item.dblinks
                        ],
                    }
                )
        except (ValueError, IndexError, TypeError) as exc:
            raise ConfigurationError(f"{source.name} is not a valid KEGG GENE file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_kegg_gene_entries",
        "source": str(source),
        "entry_count": len(entries),
        "entries": entries,
        "provenance": _provenance(
            parser="Bio.KEGG.Gene.parse",
            scope="parses a saved KEGG flat file; no KEGG REST request is made",
        ),
    }


def load_substitution_matrix(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Load a bundled Biopython substitution matrix, including DNA matrices such as NUC.4.4."""

    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("name must be a non-empty matrix name.")
    substitution = require("Bio.Align.substitution_matrices", _EXTRA)
    available = sorted(substitution.load())
    selected = name.strip()
    if selected not in available:
        raise ConfigurationError(
            f"unknown substitution matrix {selected!r}; choose from {available}."
        )
    matrix = substitution.load(selected)
    alphabet = str(matrix.alphabet)
    values = [[float(matrix[(row, column)]) for column in alphabet] for row in alphabet]
    symmetric = all(
        values[i][j] == values[j][i] for i in range(len(alphabet)) for j in range(len(alphabet))
    )
    return {
        "status": "completed",
        "operation": "load_substitution_matrix",
        "name": selected,
        "alphabet": alphabet,
        "size": len(alphabet),
        "symmetric": symmetric,
        "matrix": values,
        "provenance": _provenance(
            parser="Bio.Align.substitution_matrices.load",
            order="rows and columns follow the alphabet order",
        ),
    }


def convert_quality_score(
    value: float,
    *,
    direction: str = "solexa_to_phred",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert one Solexa and PHRED quality score using Biopython's FASTQ formulas."""

    selected = _choice(direction, _QUALITY_DIRECTIONS, "direction")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError("value must be a finite number.")
    if selected == "phred_to_solexa" and value < 0:
        raise ConfigurationError("a PHRED quality score cannot be negative.")
    quality_io = require("Bio.SeqIO.QualityIO", _EXTRA)
    if selected == "solexa_to_phred":
        converted = quality_io.phred_quality_from_solexa(float(value))
    else:
        converted = quality_io.solexa_quality_from_phred(float(value))
    return {
        "status": "completed",
        "operation": "convert_quality_score",
        "direction": selected,
        "input": float(value),
        "output": float(converted),
        "provenance": _provenance(
            parser="Bio.SeqIO.QualityIO.phred_quality_from_solexa"
            if selected == "solexa_to_phred"
            else "Bio.SeqIO.QualityIO.solexa_quality_from_phred",
        ),
    }


def read_psl_alignments(
    path: str | Path,
    *,
    max_alignments: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse PSL pairwise genomic alignments into target and query coordinate blocks."""

    source = input_path(path)
    limit = _integer(max_alignments, "max_alignments", minimum=1, maximum=_MAX_RECORDS)
    align = require("Bio.Align", _EXTRA)
    alignments: list[dict[str, Any]] = []
    try:
        for item in align.parse(str(source), "psl"):
            if len(alignments) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_alignments={limit}.")
            target_coordinates = [int(value) for value in item.coordinates[0]]
            query_coordinates = [int(value) for value in item.coordinates[1]]
            alignments.append(
                {
                    "target_id": str(item.target.id),
                    "query_id": str(item.query.id),
                    "target_start": target_coordinates[0],
                    "target_end": target_coordinates[-1],
                    "query_start": query_coordinates[0],
                    "query_end": query_coordinates[-1],
                    "target_coordinates": target_coordinates,
                    "query_coordinates": query_coordinates,
                }
            )
    except ValueError as exc:
        raise ConfigurationError(f"{source.name} is not a valid PSL file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_psl_alignments",
        "source": str(source),
        "alignment_count": len(alignments),
        "alignments": alignments,
        "provenance": _provenance(
            parser="Bio.Align.parse(format='psl')",
            coordinates="boundary positions of aligned blocks; equal steps are ungapped",
        ),
    }


def bootstrap_consensus_tree(
    sequences: Mapping[str, str],
    *,
    replicates: int = 100,
    method: str = "nj",
    distance_model: str = "identity",
    seed: int = 0,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a majority-rule bootstrap consensus tree from an aligned set of DNA sequences."""

    tree_method = _choice(method, _TREE_METHODS, "method")
    model = _choice(distance_model, _DISTANCE_MODELS, "distance_model")
    count = _integer(replicates, "replicates", minimum=1, maximum=_MAX_REPLICATES)
    seed_value = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    items = _named_sequences(
        sequences,
        "sequences",
        alphabet=_UNAMBIGUOUS_DNA,
        maximum=limit,
        minimum_rows=3,
        max_rows=200,
    )
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    align = require("Bio.Align", _EXTRA)
    consensus = require("Bio.Phylo.Consensus", _EXTRA)
    construction = require("Bio.Phylo.TreeConstruction", _EXTRA)
    multiple = align.MultipleSeqAlignment(
        [
            record_module.SeqRecord(seq_module.Seq(text), id=identifier)
            for identifier, text in items.items()
        ]
    )
    constructor = construction.DistanceTreeConstructor(
        construction.DistanceCalculator(model), tree_method
    )
    state = random.getstate()
    random.seed(seed_value)
    try:
        tree = consensus.bootstrap_consensus(
            multiple, count, constructor, consensus.majority_consensus
        )
    finally:
        random.setstate(state)
    newick = tree.format("newick").strip()
    return {
        "status": "completed",
        "operation": "bootstrap_consensus_tree",
        "ids": list(items),
        "method": tree_method,
        "distance_model": model,
        "replicates": count,
        "seed": seed_value,
        "newick": newick,
        "provenance": _provenance(
            parser="Bio.Phylo.Consensus.bootstrap_consensus",
            randomness="the global random state is seeded for the call and restored afterwards",
            consensus="majority rule with cutoff 0.5",
        ),
    }


def parse_primer3_output(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read one Primer3 (boulder-style) design report into primer pairs and internal oligos."""

    source = input_path(path)
    primer3 = require("Bio.Emboss.Primer3", _EXTRA)
    with source.open(encoding="utf-8") as handle:
        try:
            record = primer3.read(handle)
        except ValueError as exc:
            raise ConfigurationError(
                f"{source.name} is not a single Primer3 record: {exc}"
            ) from exc
    pairs = []
    for primer in record.primers:
        entry: dict[str, Any] = {
            "product_size": int(primer.size),
            "forward": {
                "sequence": str(primer.forward_seq),
                "start": int(primer.forward_start),
                "length": int(primer.forward_length),
                "tm": float(primer.forward_tm),
                "gc_percent": float(primer.forward_gc),
            },
            "reverse": {
                "sequence": str(primer.reverse_seq),
                "start": int(primer.reverse_start),
                "length": int(primer.reverse_length),
                "tm": float(primer.reverse_tm),
                "gc_percent": float(primer.reverse_gc),
            },
            "internal": None,
        }
        if int(primer.internal_length):
            entry["internal"] = {
                "sequence": str(primer.internal_seq),
                "start": int(primer.internal_start),
                "length": int(primer.internal_length),
                "tm": float(primer.internal_tm),
                "gc_percent": float(primer.internal_gc),
            }
        pairs.append(entry)
    return {
        "status": "completed",
        "operation": "parse_primer3_output",
        "source": str(source),
        "pair_count": len(pairs),
        "primer_pairs": pairs,
        "comments": str(record.comments).strip(),
        "provenance": _provenance(
            parser="Bio.Emboss.Primer3.read",
            scope="parses saved Primer3 output; does not run Primer3",
        ),
    }


__all__ = [
    "ambiguous_codon_amino_acids",
    "bootstrap_consensus_tree",
    "build_codon_alignment",
    "convert_quality_score",
    "expand_ambiguous_codons",
    "gc_content_by_codon_position",
    "gc_skew_windows",
    "load_substitution_matrix",
    "local_composition_complexity",
    "mcdonald_kreitman_test",
    "pairwise_dn_ds_matrix",
    "parse_kegg_gene_entries",
    "parse_primer3_output",
    "parse_primersearch_output",
    "read_ace_assembly",
    "read_gaf_annotations",
    "read_genepop_genotypes",
    "read_phd_sequence",
    "read_psl_alignments",
    "sequence_checksum",
]

_register_public_verbose(globals(), tuple(__all__))
