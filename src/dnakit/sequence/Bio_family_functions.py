"""Biopython functions added by family group H3 (sequence area).

Each function covers a Biopython reader or writer that the existing wrappers cannot reach through
a parameter: multi-record forms of a single-record reader (several trees, alignment blocks, Primer3
reports or PHRED records in one file), the pairwise formats of ``Bio.Align``, the BLAST XML2
format, GenBank header fields, custom substitution matrices and the GOA gene-product files.
Biopython is imported only when a function runs. Every function reads a local file or explicit
strings and writes nothing. No function changes global state: the random state used by the
bootstrap function is saved and restored.
"""

from __future__ import annotations

import importlib.metadata
import io
import random
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_UNAMBIGUOUS = frozenset("ACGT")
_NUCLEOTIDE_MATRIX_SYMBOLS = frozenset("ACGTURYSWKMBDHVN-")
_MAX_SEQUENCE = 1_000_000
_MAX_RECORDS = 100_000
_MAX_TREES = 10_000
_MAX_REPLICATES = 1_000
_MAX_FILE_BYTES = 256 * 1024 * 1024
_ALIGNMENT_FILE_FORMATS = ("bed", "clustal", "fasta", "maf", "phylip", "psl", "sam")
_ALIGNMENT_BLOCK_FORMATS = ("maf", "fasta-m10")
_TREE_FORMATS = ("newick", "nexus", "nexml", "phyloxml", "cdao")
_TREE_METHODS = ("nj", "upgma")
_DISTANCE_MODELS = ("identity", "blastn", "trans")
_GOA_FORMATS = ("gpa", "gpi")
_UNDEFINED_LENGTH = 2**40


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


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _source(path: object) -> Path:
    if not isinstance(path, (str, Path)):
        raise ConfigurationError("path must be a file path.")
    source = input_path(path)
    if source.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"{source.name} is larger than {_MAX_FILE_BYTES} bytes.")
    return source


def _plain(value: object) -> Any:
    """Keep strings, numbers and booleans; turn containers into plain lists and dictionaries."""

    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in sorted(value.items(), key=str)}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = sorted(value, key=str) if isinstance(value, (set, frozenset)) else value
        return [_plain(item) for item in items]
    return str(value)


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _nucleotides(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty sequence string.")
    text = value.strip().upper()
    if len(text) > _MAX_SEQUENCE:
        raise ConfigurationError(f"{name} may contain at most {_MAX_SEQUENCE} symbols.")
    if set(text) - _UNAMBIGUOUS:
        raise ConfigurationError(f"{name} may contain only A, C, G and T.")
    return text


def _aligned_rows(alignment: object, count: int) -> list[str] | None:
    """Return the aligned rows of a pairwise alignment, or None when a sequence is undefined."""

    try:
        return [str(alignment[index]) for index in range(count)]  # type: ignore[index]
    except Exception:
        return None


def _sequence_length(sequence: object) -> int | None:
    try:
        length = len(sequence)  # type: ignore[arg-type]
    except Exception:
        return None
    return length if 0 < length < _UNDEFINED_LENGTH else None


def read_alignment_file(
    path: str | Path,
    *,
    file_format: str = "psl",
    exactly_one: bool = False,
    max_alignments: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read pairwise alignments from a BED, Clustal, FASTA, MAF, PHYLIP, PSL or SAM file.

    Calls ``Bio.Align.parse``. With ``exactly_one=True`` the call is ``Bio.Align.read``: it
    raises ``ConfigurationError`` unless the file holds exactly one alignment. Each alignment
    reports its sequence names, the aligned columns, the score when the format stores one, the
    block coordinates (one row per sequence) and the gapped rows when the file carries the
    sequences. PSL, SAM and BED files store no sequence, so ``aligned_rows`` is null for them.
    """

    source = _source(path)
    fmt = _choice(file_format, _ALIGNMENT_FILE_FORMATS, "file_format")
    single = _flag(exactly_one, "exactly_one")
    limit = _integer(max_alignments, "max_alignments", minimum=1, maximum=_MAX_RECORDS)
    align = require("Bio.Align", _EXTRA)
    alignments: list[dict[str, Any]] = []
    try:
        if single:
            found: list[Any] = [align.read(str(source), fmt)]
        else:
            found = []
            for item in align.parse(str(source), fmt):
                if len(found) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_alignments={limit}.")
                found.append(item)
        for index, item in enumerate(found):
            rows = len(item.sequences)
            alignments.append(
                {
                    "index": index,
                    "sequence_ids": [str(getattr(seq, "id", "")) for seq in item.sequences],
                    "sequence_lengths": [_sequence_length(seq) for seq in item.sequences],
                    "aligned_columns": int(item.length),
                    "score": _optional_float(getattr(item, "score", None)),
                    "coordinates": [[int(v) for v in row] for row in item.coordinates.tolist()],
                    "aligned_rows": _aligned_rows(item, rows),
                }
            )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(
            f"{source.name} is not a valid {fmt} alignment file: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "read_alignment_file",
        "source": str(source),
        "file_format": fmt,
        "exactly_one": single,
        "alignment_count": len(alignments),
        "alignments": alignments,
        "provenance": _provenance(
            parser="Bio.Align.read" if single else "Bio.Align.parse",
            coordinates="block boundaries, 0-based and half-open, one row per sequence",
        ),
    }


def read_alignment_blocks(
    path: str | Path,
    *,
    file_format: str = "maf",
    max_blocks: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every alignment block of a MAF or FASTA-M10 file with ``Bio.AlignIO.parse``.

    A MAF file holds one block per aligned region and a FASTA-M10 file one alignment per hit,
    so the single-alignment reader (``read_alignment``) rejects them. Each block lists its rows
    (identifier, gapped sequence and the per-row annotations such as the MAF ``start``,
    ``size``, ``strand`` and ``srcSize``) and the block annotations such as the MAF ``score``.
    """

    source = _source(path)
    fmt = _choice(file_format, _ALIGNMENT_BLOCK_FORMATS, "file_format")
    limit = _integer(max_blocks, "max_blocks", minimum=1, maximum=_MAX_RECORDS)
    align_io = require("Bio.AlignIO", _EXTRA)
    blocks: list[dict[str, Any]] = []
    try:
        with source.open(encoding="utf-8") as handle:
            for block in align_io.parse(handle, fmt):
                if len(blocks) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_blocks={limit}.")
                blocks.append(
                    {
                        "index": len(blocks),
                        "row_count": len(block),
                        "alignment_length": int(block.get_alignment_length()),
                        "annotations": _plain(block.annotations),
                        "rows": [
                            {
                                "id": str(record.id),
                                "name": str(record.name),
                                "sequence": str(record.seq),
                                "annotations": _plain(record.annotations),
                            }
                            for record in block
                        ],
                    }
                )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_alignment_blocks",
        "source": str(source),
        "file_format": fmt,
        "block_count": len(blocks),
        "blocks": blocks,
        "provenance": _provenance(
            parser="Bio.AlignIO.parse",
            coordinates="MAF start values are 0-based, as written in the file",
        ),
    }


def read_substitution_matrix_file(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a custom nucleotide substitution matrix from a text file.

    Calls ``Bio.Align.substitution_matrices.read``. The file has optional ``#`` comment lines, a
    header row of symbols and one row per symbol. The alphabet must contain only nucleotide
    symbols (A, C, G, T, U and the IUPAC codes) so the result can score DNA alignments. Bundled
    matrices are loaded by name with ``load_substitution_matrix``.
    """

    source = _source(path)
    substitution = require("Bio.Align.substitution_matrices", _EXTRA)
    try:
        with source.open(encoding="utf-8") as handle:
            matrix = substitution.read(handle)
    except Exception as exc:
        raise ConfigurationError(
            f"{source.name} is not a valid substitution matrix: {exc}"
        ) from exc
    alphabet = str(matrix.alphabet)
    if not alphabet or set(alphabet.upper()) - _NUCLEOTIDE_MATRIX_SYMBOLS:
        raise ConfigurationError(
            f"{source.name} has a non-nucleotide alphabet {alphabet!r}; "
            "only A, C, G, T, U and IUPAC symbols are accepted."
        )
    values = [[float(matrix[(row, column)]) for column in alphabet] for row in alphabet]
    header = [str(line).strip() for line in getattr(matrix, "header", None) or []]
    return {
        "status": "completed",
        "operation": "read_substitution_matrix_file",
        "source": str(source),
        "alphabet": alphabet,
        "size": len(alphabet),
        "symmetric": all(
            values[i][j] == values[j][i] for i in range(len(alphabet)) for j in range(i)
        ),
        "matrix": values,
        "comments": [line for line in header if line],
        "provenance": _provenance(
            parser="Bio.Align.substitution_matrices.read",
            order="rows and columns follow the alphabet order",
        ),
    }


def _hsp_row(hsp: Any) -> dict[str, Any]:
    annotations = hsp.annotations
    coordinates = [[int(v) for v in row] for row in hsp.coordinates.tolist()]
    return {
        "score": _optional_float(hsp.score),
        "bit_score": _optional_float(annotations.get("bit score")),
        "evalue": _optional_float(annotations.get("evalue")),
        "identity": annotations.get("identity"),
        "positive": annotations.get("positive"),
        "gaps": annotations.get("gaps"),
        "aligned_columns": int(hsp.length),
        "target_coordinates": coordinates[0],
        "query_coordinates": coordinates[1],
    }


def read_blast_records(
    path: str | Path,
    *,
    exactly_one: bool = False,
    max_records: int = 1_000,
    max_hits: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read BLAST XML or XML2 output with ``Bio.Blast.parse`` (one record per query).

    With ``exactly_one=True`` the call is ``Bio.Blast.read``: it raises ``ConfigurationError``
    unless the file holds exactly one query. Each record lists the query, the program, version
    and database, and every hit with its high-scoring pairs (score, bit score, E-value, identity
    counts and the target and query block coordinates, 0-based and half-open). The SearchIO
    parser behind ``parse_search_results`` cannot read the XML2 layout.
    """

    source = _source(path)
    single = _flag(exactly_one, "exactly_one")
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    hit_limit = _integer(max_hits, "max_hits", minimum=1, maximum=_MAX_RECORDS)
    blast = require("Bio.Blast", _EXTRA)
    rows: list[dict[str, Any]] = []
    try:
        if single:
            found: list[Any] = [blast.read(str(source))]
        else:
            found = []
            with blast.parse(str(source)) as records:
                for record in records:
                    if len(found) >= record_limit:
                        raise ConfigurationError(
                            f"{source.name} exceeds max_records={record_limit}."
                        )
                    for key in ("source", "program", "version", "reference", "db"):
                        if hasattr(records, key):
                            setattr(record, key, getattr(records, key))
                    found.append(record)
        for record in found:
            if len(record) > hit_limit:
                raise ConfigurationError(f"a record exceeds max_hits={hit_limit}.")
            rows.append(
                {
                    "query_id": str(record.query.id),
                    "query_description": str(record.query.description),
                    "query_length": len(record.query.seq) if record.query.seq is not None else None,
                    "program": getattr(record, "program", None),
                    "version": getattr(record, "version", None),
                    "database": getattr(record, "db", None),
                    "hit_count": len(record),
                    "hits": [
                        {
                            "target_id": str(hit.target.id),
                            "target_name": str(hit.target.name),
                            "target_description": str(hit.target.description),
                            "target_length": _sequence_length(hit.target.seq),
                            "hsps": [_hsp_row(hsp) for hsp in hit],
                        }
                        for hit in record
                    ],
                }
            )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not valid BLAST XML output: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_blast_records",
        "source": str(source),
        "exactly_one": single,
        "record_count": len(rows),
        "records": rows,
        "provenance": _provenance(
            parser="Bio.Blast.read" if single else "Bio.Blast.parse",
            coordinates="HSP block boundaries are 0-based and half-open",
        ),
    }


def _primer_pair(primer: Any) -> dict[str, Any]:
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
    return entry


def parse_primer3_records(
    path: str | Path,
    *,
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every record of an EMBOSS eprimer3 or Primer3 report (one per target sequence).

    Calls ``Bio.Emboss.Primer3.parse``. A report for several target sequences holds one record
    per target; ``parse_primer3_output`` accepts only a single-record report. Each record gives
    its header comments and the primer pairs with an optional internal oligo.
    """

    source = _source(path)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    primer3 = require("Bio.Emboss.Primer3", _EXTRA)
    text = source.read_text(encoding="utf-8")
    if not text.strip():
        # Bio.Emboss.Primer3.parse loops forever on an empty stream.
        raise ConfigurationError(f"{source.name} is empty.")
    records: list[dict[str, Any]] = []
    try:
        for record in primer3.parse(io.StringIO(text)):
            if len(records) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
            pairs = [_primer_pair(primer) for primer in record.primers]
            records.append(
                {
                    "index": len(records),
                    "comments": str(record.comments).strip(),
                    "pair_count": len(pairs),
                    "primer_pairs": pairs,
                }
            )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid Primer3 report: {exc}") from exc
    if not records:
        raise ConfigurationError(f"{source.name} holds no Primer3 record.")
    return {
        "status": "completed",
        "operation": "parse_primer3_records",
        "source": str(source),
        "record_count": len(records),
        "pair_count": sum(item["pair_count"] for item in records),
        "records": records,
        "provenance": _provenance(
            parser="Bio.Emboss.Primer3.parse",
            scope="parses saved Primer3 output; does not run Primer3",
        ),
    }


def read_genbank_record_fields(
    path: str | Path,
    *,
    exactly_one: bool = False,
    max_records: int = 1_000,
    max_bases: int = 50_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the header fields and feature table of GenBank records with ``Bio.GenBank.parse``.

    The low-level GenBank parser keeps fields that ``read_sequence_records`` drops: the LOCUS
    line (name, length, molecule type, topology text, division and date), accessions, version,
    keywords, source and organism, taxonomy, the number of references, the comment, the base
    counts and every feature with its location text and qualifiers. With ``exactly_one=True``
    the call is ``Bio.GenBank.read`` and a file with several records raises
    ``ConfigurationError``.
    """

    source = _source(path)
    single = _flag(exactly_one, "exactly_one")
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    base_limit = _integer(max_bases, "max_bases", minimum=1, maximum=1_000_000_000)
    genbank = require("Bio.GenBank", _EXTRA)
    rows: list[dict[str, Any]] = []
    total = 0
    try:
        with source.open(encoding="utf-8") as handle:
            found = [genbank.read(handle)] if single else genbank.parse(handle)
            for record in found:
                if len(rows) >= record_limit:
                    raise ConfigurationError(f"{source.name} exceeds max_records={record_limit}.")
                sequence = str(record.sequence).upper()
                total += len(sequence)
                if total > base_limit:
                    raise ConfigurationError(f"{source.name} exceeds max_bases={base_limit}.")
                rows.append(
                    {
                        "locus": str(record.locus),
                        "size": int(record.size) if str(record.size).isdigit() else None,
                        "residue_type": " ".join(str(record.residue_type).split()),
                        "data_file_division": str(record.data_file_division),
                        "date": str(record.date),
                        "accessions": [str(item) for item in record.accession],
                        "version": str(record.version),
                        "keywords": [str(item) for item in record.keywords if str(item)],
                        "source": str(record.source),
                        "organism": str(record.organism),
                        "taxonomy": [str(item) for item in record.taxonomy],
                        "reference_count": len(record.references),
                        "comment": str(record.comment),
                        "base_counts": str(record.base_counts).strip(),
                        "sequence_length": len(sequence),
                        "sequence": sequence,
                        "features": [
                            {
                                "type": str(feature.key),
                                "location": str(feature.location),
                                "qualifiers": [
                                    {"key": str(item.key).strip("=/"), "value": str(item.value)}
                                    for item in feature.qualifiers
                                ],
                            }
                            for feature in record.features
                        ],
                    }
                )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid GenBank file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_genbank_record_fields",
        "source": str(source),
        "exactly_one": single,
        "record_count": len(rows),
        "total_bases": total,
        "records": rows,
        "provenance": _provenance(
            parser="Bio.GenBank.read" if single else "Bio.GenBank.parse",
            coordinates="feature locations are the 1-based, inclusive text written in the file",
        ),
    }


def read_phylogeny_collection(
    path: str | Path,
    *,
    file_format: str = "newick",
    max_trees: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every tree of a Newick, NEXUS, NeXML, phyloXML or CDAO file.

    Calls ``Bio.Phylo.parse``, which reaches the codec parsers ``NexusIO.parse``,
    ``NeXMLIO.parse``, ``PhyloXMLIO.parse`` and ``CDAOIO.parse``. ``read_tree`` reads one tree
    and raises on a file with several, so replicate sets such as bootstrap or posterior trees
    need this function. Each tree reports its name, rooted flag, tip names (sorted), total
    branch length and a Newick string.
    """

    source = _source(path)
    fmt = _choice(file_format, _TREE_FORMATS, "file_format")
    limit = _integer(max_trees, "max_trees", minimum=1, maximum=_MAX_TREES)
    phylo = require("Bio.Phylo", _EXTRA)
    trees: list[dict[str, Any]] = []
    try:
        for tree in phylo.parse(str(source), fmt):
            if len(trees) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_trees={limit}.")
            tips = sorted(str(clade.name) for clade in tree.get_terminals())
            trees.append(
                {
                    "index": len(trees),
                    "name": None if not getattr(tree, "name", None) else str(tree.name),
                    "rooted": bool(tree.rooted),
                    "tip_count": len(tips),
                    "tips": tips,
                    "total_branch_length": float(tree.total_branch_length()),
                    "newick": str(tree.format("newick")).strip(),
                }
            )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} tree file: {exc}") from exc
    if not trees:
        raise ConfigurationError(f"{source.name} holds no {fmt} tree.")
    return {
        "status": "completed",
        "operation": "read_phylogeny_collection",
        "source": str(source),
        "file_format": fmt,
        "tree_count": len(trees),
        "trees": trees,
        "provenance": _provenance(parser="Bio.Phylo.parse"),
    }


def bootstrap_replicate_trees(
    sequences: Mapping[str, str],
    *,
    replicates: int = 10,
    method: str = "nj",
    distance_model: str = "identity",
    seed: int = 0,
    max_length: int = _MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the bootstrap replicate trees of an aligned set of DNA sequences.

    Calls ``Bio.Phylo.Consensus.bootstrap_trees``: alignment columns are resampled with
    replacement and one tree is built per replicate. ``bootstrap_consensus_tree`` returns only
    the majority-rule consensus; this function returns the replicate trees themselves, which
    ``bootstrap_branch_support`` can count against a target tree. The global random state is
    seeded for the call and restored afterwards.
    """

    tree_method = _choice(method, _TREE_METHODS, "method")
    model = _choice(distance_model, _DISTANCE_MODELS, "distance_model")
    count = _integer(replicates, "replicates", minimum=1, maximum=_MAX_REPLICATES)
    seed_value = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE)
    if not isinstance(sequences, Mapping) or not 3 <= len(sequences) <= 200:
        raise ConfigurationError("sequences must map 3 to 200 identifiers to aligned DNA strings.")
    items: dict[str, str] = {}
    for key in sorted(sequences, key=str):
        identifier = str(key).strip()
        if not identifier or identifier in items:
            raise ConfigurationError("sequence identifiers must be unique and non-empty.")
        text = _nucleotides(sequences[key], f"sequences[{identifier}]")
        if len(text) > limit:
            raise ConfigurationError(f"sequences[{identifier}] is longer than {limit}.")
        items[identifier] = text
    if len({len(text) for text in items.values()}) != 1:
        raise ConfigurationError("all sequences must have the same length.")
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
        trees = [
            str(tree.format("newick")).strip()
            for tree in consensus.bootstrap_trees(multiple, count, constructor)
        ]
    finally:
        random.setstate(state)
    return {
        "status": "completed",
        "operation": "bootstrap_replicate_trees",
        "ids": list(items),
        "method": tree_method,
        "distance_model": model,
        "replicates": count,
        "seed": seed_value,
        "trees": trees,
        "provenance": _provenance(
            parser="Bio.Phylo.Consensus.bootstrap_trees",
            randomness="the global random state is seeded for the call and restored afterwards",
        ),
    }


def read_phd_records(
    path: str | Path, *, max_records: int = 10_000, verbose: bool = True
) -> dict[str, Any]:
    """Read every PHRED base-call record of a PHD file with ``Bio.Sequencing.Phd.parse``.

    ``read_phd_sequence`` reads a file with exactly one record. Each record here gives the file
    name, the base calls, the per-base PHRED qualities, the peak locations (null when the file
    gives none) and the comment fields.
    """

    source = _source(path)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    phd = require("Bio.Sequencing.Phd", _EXTRA)
    records: list[dict[str, Any]] = []
    try:
        for record in phd.parse(str(source)):
            if len(records) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
            sites = record.sites
            records.append(
                {
                    "index": len(records),
                    "file_name": str(record.file_name),
                    "sequence": "".join(str(site[0]) for site in sites),
                    "length": len(sites),
                    "qualities": [int(site[1]) for site in sites],
                    "peak_locations": [
                        int(site[2]) if len(site) > 2 and site[2] is not None else None
                        for site in sites
                    ],
                    "comments": {
                        key: value
                        for key, value in sorted(record.comments.items())
                        if value is not None
                    },
                }
            )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid PHD file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_phd_records",
        "source": str(source),
        "record_count": len(records),
        "records": records,
        "provenance": _provenance(
            parser="Bio.Sequencing.Phd.parse",
            positions="peak locations are 1-based positions in the trace, when the file gives them",
        ),
    }


def read_goa_gene_product_file(
    path: str | Path,
    *,
    file_format: str = "gpa",
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GOA gene-product association (GPA) or gene-product information (GPI) file.

    ``file_format="gpa"`` calls ``Bio.UniProt.GOA.gpa_iterator`` (versions 1.0 and 1.1:
    gene product, GO term, evidence and reference per line); ``"gpi"`` calls ``gpi_iterator``
    (versions 1.0, 1.1 and 1.2: gene product symbol, name, synonyms, type, taxon and cross
    references). The version comes from the ``!gpa-version`` or ``!gpi-version`` header line,
    which is required. ``read_gaf_annotations`` reads the GAF layout.
    """

    source = _source(path)
    fmt = _choice(file_format, _GOA_FORMATS, "file_format")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    goa = require("Bio.UniProt.GOA", _EXTRA)
    records: list[dict[str, Any]] = []
    try:
        with source.open(encoding="utf-8") as handle:
            header = handle.readline().strip()
            handle.seek(0)
            iterator = goa.gpa_iterator(handle) if fmt == "gpa" else goa.gpi_iterator(handle)
            for item in iterator:
                if len(records) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
                records.append(
                    {
                        str(key): [part for part in value if part]
                        if isinstance(value, list)
                        else value
                        for key, value in item.items()
                    }
                )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} file: {exc}") from exc
    version = header.split(":", 1)[1].strip() if ":" in header else None
    return {
        "status": "completed",
        "operation": "read_goa_gene_product_file",
        "source": str(source),
        "file_format": fmt,
        "version": version,
        "record_count": len(records),
        "records": records,
        "provenance": _provenance(
            parser="Bio.UniProt.GOA.gpa_iterator"
            if fmt == "gpa"
            else "Bio.UniProt.GOA.gpi_iterator",
            scope="annotation records are read from a local file only",
        ),
    }


__all__ = [
    "bootstrap_replicate_trees",
    "parse_primer3_records",
    "read_alignment_blocks",
    "read_alignment_file",
    "read_blast_records",
    "read_genbank_record_fields",
    "read_goa_gene_product_file",
    "read_phd_records",
    "read_phylogeny_collection",
    "read_substitution_matrix_file",
]

_register_public_verbose(globals(), tuple(__all__))
