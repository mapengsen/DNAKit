"""Biopython functions added by review group G18o (sequence area).

Each function imports Biopython only when it runs, so this module imports without the
``bio-formats`` extra. Inputs are DNA or RNA strings, Newick text, aligned or unaligned records
and local files. Outputs are plain dictionaries that can be serialized to JSON. Writers never
overwrite an existing file: the output path must be new.
"""

from __future__ import annotations

import importlib.metadata
import io
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_DNA_SYMBOLS = frozenset("ACGTRYSWKMBDHVN")
_RNA_SYMBOLS = frozenset("ACGURYSWKMBDHVN")
_UNAMBIGUOUS = frozenset("ACGT")
_MAX_SEQUENCE = 1_000_000
_MAX_RECORDS = 100_000
_MAX_TREES = 1_000
_MAX_TREE_TEXT = 5_000_000
_MAX_MOTIF_INSTANCES = 10_000
_MAX_MOTIFS = 1_000
_ALIGNMENT_FORMATS = ("bed", "clustal", "fasta", "maf", "phylip", "psl", "sam")
_ALIGNMENT_MODES = ("global", "local")
_SEARCH_INPUT_FORMATS = (
    "blast-tab",
    "blast-xml",
    "blat-psl",
    "hmmer3-tab",
    "hmmscan3-domtab",
    "hmmsearch3-domtab",
    "phmmer3-domtab",
)
_SEARCH_OUTPUT_FORMATS = _SEARCH_INPUT_FORMATS
_MOTIF_FORMATS = ("clusterbuster", "jaspar", "pfm", "transfac")
_TM_VALUESETS = tuple(range(1, 8))
_TM_SALT_METHODS = tuple(range(0, 8))
_CODON = re.compile(r"^[ACGT]{3}$")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _sequence(value: object, name: str, alphabet: frozenset[str], *, gaps: bool = False) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty sequence string.")
    text = value.strip().upper()
    if len(text) > _MAX_SEQUENCE:
        raise ConfigurationError(f"{name} may contain at most {_MAX_SEQUENCE} symbols.")
    allowed = alphabet | {"-"} if gaps else alphabet
    if set(text) - allowed:
        raise ConfigurationError(f"{name} contains symbols outside the allowed alphabet.")
    return text


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float = 0.0) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    number = float(value)
    if number != number or number in (float("inf"), float("-inf")) or number < minimum:
        raise ConfigurationError(f"{name} must be a finite number >= {minimum}.")
    return number


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _new_output(value: str | Path, name: str) -> Path:
    """Return a path that does not exist yet and whose directory exists."""

    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new file: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"{name} parent directory does not exist: {path.parent}")
    return path


def _text(value: object, name: str, *, limit: int = _MAX_TREE_TEXT) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} is longer than {limit} characters.")
    return value


def _read_newick(newick: object, name: str) -> Any:
    phylo = require("Bio.Phylo", _EXTRA)
    text = _text(newick, name)
    try:
        return phylo.read(io.StringIO(text), "newick")
    except Exception as exc:
        raise ConfigurationError(f"{name} is not a valid Newick tree: {exc}") from exc


def _newick_text(tree: Any) -> str:
    return str(tree.format("newick")).strip()


def write_pairwise_alignment_file(
    query: str,
    target: str,
    output_path: str | Path,
    *,
    file_format: str = "fasta",
    mode: str = "global",
    query_id: str = "query",
    target_id: str = "target",
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two DNA sequences with Biopython and write the best alignment to a new file.

    The aligner uses fixed scores (match 2, mismatch -1, gap open -2, gap extend -1), so the
    result is reproducible. ``file_format`` is one of the text formats that Bio.Align writes for
    pairwise alignments: fasta, clustal, maf, phylip, psl, sam or bed.
    """

    fmt = _choice(file_format, _ALIGNMENT_FORMATS, "file_format")
    alignment_mode = _choice(mode, _ALIGNMENT_MODES, "mode")
    query_text = _sequence(query, "query", _DNA_SYMBOLS)
    target_text = _sequence(target, "target", _DNA_SYMBOLS)
    for label, identifier in (("query_id", query_id), ("target_id", target_id)):
        if (
            not isinstance(identifier, str)
            or not identifier
            or any(c.isspace() for c in identifier)
        ):
            raise ConfigurationError(f"{label} must be a non-empty string without spaces.")
    destination = _new_output(output_path, "output_path")
    align = require("Bio.Align", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    aligner = align.PairwiseAligner(
        mode=alignment_mode,
        match_score=2.0,
        mismatch_score=-1.0,
        open_gap_score=-2.0,
        extend_gap_score=-1.0,
    )
    left = record_module.SeqRecord(seq_module.Seq(query_text), id=query_id, description="")
    right = record_module.SeqRecord(seq_module.Seq(target_text), id=target_id, description="")
    best = aligner.align(left, right)[0]
    with destination.open("x", encoding="utf-8") as handle:
        written = align.write([best], handle, fmt)
    return {
        "status": "completed",
        "operation": "write_pairwise_alignment_file",
        "output_path": str(destination),
        "file_format": fmt,
        "mode": alignment_mode,
        "score": float(best.score),
        "aligned_query": str(best[0]),
        "aligned_target": str(best[1]),
        "alignment_count": int(written),
        "provenance": _provenance(
            method="Bio.Align.PairwiseAligner with match 2, mismatch -1, gap open -2, extend -1",
            writer="Bio.Align.write",
        ),
    }


def convert_search_results(
    source: str | Path,
    source_format: str,
    output_path: str | Path,
    output_format: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a local sequence-search report (BLAST, BLAT or HMMER) to another format.

    The output file must not exist. Formats follow Bio.SearchIO names, for example
    ``blast-tab`` or ``blast-xml``. The returned count is the number of query results.
    """

    src = input_path(source)
    in_format = _choice(source_format, _SEARCH_INPUT_FORMATS, "source_format")
    out_format = _choice(output_format, _SEARCH_OUTPUT_FORMATS, "output_format")
    destination = _new_output(output_path, "output_path")
    searchio = require("Bio.SearchIO", _EXTRA)
    try:
        counted = searchio.convert(str(src), in_format, str(destination), out_format)
    except (AttributeError, ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"Bio.SearchIO could not convert {src.name}: {exc}") from exc
    counts = tuple(counted) if isinstance(counted, tuple) else (counted, None, None, None)
    return {
        "status": "completed",
        "operation": "convert_search_results",
        "source": str(src),
        "source_format": in_format,
        "output_path": str(destination),
        "output_format": out_format,
        "query_count": int(counts[0]),
        "hit_count": None if counts[1] is None else int(counts[1]),
        "hsp_count": None if counts[2] is None else int(counts[2]),
        "fragment_count": None if counts[3] is None else int(counts[3]),
        "provenance": _provenance(method="Bio.SearchIO.convert"),
    }


def codon_back_table(
    genetic_code: int = 1,
    *,
    stop_codon: str = "TAA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return one codon for each amino acid of an NCBI genetic code, with stop as ``*``.

    Biopython keeps a single codon per amino acid: the last one in sorted codon order. Stop is
    keyed as ``*`` and mapped to ``stop_codon``, which must be a stop codon of the table.
    """

    codon_table = require("Bio.Data.CodonTable", _EXTRA)
    table_id = _integer(genetic_code, "genetic_code", minimum=1, maximum=33)
    try:
        table = codon_table.unambiguous_dna_by_id[table_id]
    except KeyError as exc:
        raise ConfigurationError(
            f"genetic_code {table_id} is not an unambiguous DNA table."
        ) from exc
    stop = _text(stop_codon, "stop_codon", limit=3).upper()
    if not _CODON.match(stop) or stop not in table.stop_codons:
        raise ConfigurationError(f"stop_codon must be a stop codon of table {table_id}.")
    back = codon_table.make_back_table(table.forward_table, stop)
    mapping = {amino: back[amino] for amino in sorted(key for key in back if key is not None)}
    mapping["*"] = stop
    return {
        "status": "completed",
        "operation": "codon_back_table",
        "genetic_code": table_id,
        "table_name": table.names[0],
        "codon_count": len(table.forward_table),
        "back_table": mapping,
        "provenance": _provenance(method="Bio.Data.CodonTable.make_back_table"),
    }


def adam_consensus_tree(newick_trees: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Return the Adam consensus of several Newick trees.

    Each tree must parse as Newick. Biopython's Adam consensus is returned unchanged, so its
    topology is the upstream result for the given trees.
    """

    if not isinstance(newick_trees, Sequence) or isinstance(newick_trees, str):
        raise ConfigurationError("newick_trees must be a sequence of Newick strings.")
    if not 1 <= len(newick_trees) <= _MAX_TREES:
        raise ConfigurationError(f"newick_trees must contain 1 to {_MAX_TREES} trees.")
    trees = [_read_newick(text, f"newick_trees[{i}]") for i, text in enumerate(newick_trees)]
    consensus = require("Bio.Phylo.Consensus", _EXTRA)
    result = consensus.adam_consensus(trees)
    return {
        "status": "completed",
        "operation": "adam_consensus_tree",
        "tree_count": len(trees),
        "newick": _newick_text(result),
        "provenance": _provenance(method="Bio.Phylo.Consensus.adam_consensus"),
    }


def bootstrap_branch_support(
    target_newick: str,
    replicate_newicks: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate a target tree with the percentage of bootstrap trees that contain each clade.

    Support is written as the confidence of each internal clade, on a 0 to 100 scale, in the
    returned Newick string. Replicate trees must carry the same tip names as the target.
    """

    if not isinstance(replicate_newicks, Sequence) or isinstance(replicate_newicks, str):
        raise ConfigurationError("replicate_newicks must be a sequence of Newick strings.")
    if not 1 <= len(replicate_newicks) <= _MAX_TREES:
        raise ConfigurationError(f"replicate_newicks must contain 1 to {_MAX_TREES} trees.")
    target = _read_newick(target_newick, "target_newick")
    replicates = [
        _read_newick(text, f"replicate_newicks[{i}]") for i, text in enumerate(replicate_newicks)
    ]
    consensus = require("Bio.Phylo.Consensus", _EXTRA)
    annotated = consensus.get_support(target, replicates)
    return {
        "status": "completed",
        "operation": "bootstrap_branch_support",
        "replicate_count": len(replicates),
        "newick": _newick_text(annotated),
        "scale": "percent",
        "provenance": _provenance(method="Bio.Phylo.Consensus.get_support"),
    }


def tree_to_graph_edges(newick: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the nodes and parent-to-child edges of a Newick tree as plain lists.

    Node ids are assigned in pre-order from 0 at the root, so the output is stable. This is the
    same graph that Bio.Phylo.to_networkx and Bio.Phylo.to_igraph build, without the optional
    graph libraries.
    """

    tree = _read_newick(newick, "newick")
    clades = list(tree.find_clades(order="preorder"))
    index = {id(clade): position for position, clade in enumerate(clades)}
    nodes = [
        {
            "id": index[id(clade)],
            "name": clade.name,
            "branch_length": None if clade.branch_length is None else float(clade.branch_length),
            "is_tip": not clade.clades,
        }
        for clade in clades
    ]
    edges: list[dict[str, Any]] = []
    for clade in clades:
        for child in clade.clades:
            edges.append(
                {
                    "source": index[id(clade)],
                    "target": index[id(child)],
                    "branch_length": None
                    if child.branch_length is None
                    else float(child.branch_length),
                }
            )
    return {
        "status": "completed",
        "operation": "tree_to_graph_edges",
        "node_count": len(nodes),
        "edge_count": len(edges),
        "tip_count": sum(1 for node in nodes if node["is_tip"]),
        "nodes": nodes,
        "edges": edges,
        "provenance": _provenance(method="Bio.Phylo tree structure in pre-order"),
    }


def complement_rna_sequence(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the RNA complement of an RNA string (A-U, C-G) with IUPAC codes.

    The output keeps the U alphabet. Upper-case input is required after normalization.
    """

    seq_module = require("Bio.Seq", _EXTRA)
    text = _sequence(sequence, "sequence", _RNA_SYMBOLS)
    result = str(seq_module.complement_rna(text))
    return {
        "status": "completed",
        "operation": "complement_rna_sequence",
        "input": text,
        "sequence": result,
        "length": len(result),
        "provenance": _provenance(method="Bio.Seq.complement_rna"),
    }


def reverse_complement_rna_sequence(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the reverse complement of an RNA string, keeping U and IUPAC codes."""

    seq_module = require("Bio.Seq", _EXTRA)
    text = _sequence(sequence, "sequence", _RNA_SYMBOLS)
    result = str(seq_module.reverse_complement_rna(text))
    return {
        "status": "completed",
        "operation": "reverse_complement_rna_sequence",
        "input": text,
        "sequence": result,
        "length": len(result),
        "provenance": _provenance(method="Bio.Seq.reverse_complement_rna"),
    }


def read_fasta_qual_pairs(
    fasta_path: str | Path,
    qual_path: str | Path,
    *,
    max_records: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read matched FASTA and QUAL files as sequences with PHRED qualities.

    Records are paired by position, as Biopython does. Qualities are integers. The result is
    limited to ``max_records`` records.
    """

    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    fasta = input_path(fasta_path)
    qual = input_path(qual_path)
    quality_io = require("Bio.SeqIO.QualityIO", _EXTRA)
    records: list[dict[str, Any]] = []
    try:
        for record in quality_io.PairedFastaQualIterator(str(fasta), str(qual)):
            if len(records) >= limit:
                raise ConfigurationError(f"The FASTA/QUAL pair holds more than {limit} records.")
            records.append(
                {
                    "id": record.id,
                    "description": record.description,
                    "sequence": str(record.seq),
                    "qualities": [
                        int(value) for value in record.letter_annotations["phred_quality"]
                    ],
                }
            )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"FASTA and QUAL files do not pair: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_fasta_qual_pairs",
        "fasta_path": str(fasta),
        "qual_path": str(qual),
        "record_count": len(records),
        "records": records,
        "provenance": _provenance(method="Bio.SeqIO.QualityIO.PairedFastaQualIterator"),
    }


def write_qual_file(
    records: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write PHRED qualities of DNA records to a new QUAL file.

    Each record needs an ``id``, a DNA ``sequence`` and one integer quality per base in
    ``qualities`` (0 to 93). The sequence is used only to check the quality count.
    """

    if not isinstance(records, Sequence) or isinstance(records, str) or not records:
        raise ConfigurationError("records must be a non-empty sequence of mappings.")
    if len(records) > _MAX_RECORDS:
        raise ConfigurationError(f"records may contain at most {_MAX_RECORDS} entries.")
    destination = _new_output(output_path, "output_path")
    require("Bio.SeqIO.QualityIO", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    built = []
    for index, item in enumerate(records):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"records[{index}] must be a mapping.")
        identifier = item.get("id")
        if (
            not isinstance(identifier, str)
            or not identifier
            or any(c.isspace() for c in identifier)
        ):
            raise ConfigurationError(
                f"records[{index}].id must be a non-empty string without spaces."
            )
        text = _sequence(item.get("sequence"), f"records[{index}].sequence", _DNA_SYMBOLS)
        qualities = item.get("qualities")
        if not isinstance(qualities, Sequence) or isinstance(qualities, str):
            raise ConfigurationError(f"records[{index}].qualities must be a list of integers.")
        if len(qualities) != len(text):
            raise ConfigurationError(
                f"records[{index}] has {len(qualities)} qualities for {len(text)} bases."
            )
        scores = [
            _integer(value, f"records[{index}].qualities", minimum=0, maximum=93)
            for value in qualities
        ]
        record = record_module.SeqRecord(seq_module.Seq(text), id=identifier, description="")
        record.letter_annotations["phred_quality"] = scores
        built.append(record)
    with destination.open("x", encoding="utf-8") as handle:
        for record in built:
            handle.write(format(record, "qual"))
    return {
        "status": "completed",
        "operation": "write_qual_file",
        "output_path": str(destination),
        "record_count": len(built),
        "provenance": _provenance(method="Bio.SeqIO format 'qual' (QualityIO writer)"),
    }


def melting_temperature_gc(
    sequence: str,
    *,
    valueset: int = 7,
    na: float = 50.0,
    k: float = 0.0,
    tris: float = 0.0,
    mg: float = 0.0,
    dntps: float = 0.0,
    saltcorr: int = 0,
    mismatch: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the melting temperature of an unambiguous DNA primer from its GC content.

    This is the empirical GC-based formula of Biopython (``Tm_GC``), not the nearest-neighbour
    model of ``melting_temperature``. ``valueset`` selects one of seven published coefficient
    sets. Salt concentrations are in mM and ``saltcorr`` (0 to 7) selects a salt correction.
    """

    sequence_text = _sequence(sequence, "sequence", _UNAMBIGUOUS)
    chosen_valueset = _integer(
        valueset, "valueset", minimum=_TM_VALUESETS[0], maximum=_TM_VALUESETS[-1]
    )
    salt_method = _integer(
        saltcorr, "saltcorr", minimum=_TM_SALT_METHODS[0], maximum=_TM_SALT_METHODS[-1]
    )
    salts = {
        "Na": _number(na, "na"),
        "K": _number(k, "k"),
        "Tris": _number(tris, "tris"),
        "Mg": _number(mg, "mg"),
        "dNTPs": _number(dntps, "dntps"),
    }
    flag = _flag(mismatch, "mismatch")
    seq_utils = require("Bio.SeqUtils.MeltingTemp", _EXTRA)
    tm = seq_utils.Tm_GC(
        sequence_text,
        check=True,
        strict=True,
        valueset=chosen_valueset,
        Na=salts["Na"],
        K=salts["K"],
        Tris=salts["Tris"],
        Mg=salts["Mg"],
        dNTPs=salts["dNTPs"],
        saltcorr=salt_method,
        mismatch=flag,
    )
    return {
        "status": "completed",
        "operation": "melting_temperature_gc",
        "sequence": sequence_text,
        "length": len(sequence_text),
        "tm_celsius": float(tm),
        "valueset": chosen_valueset,
        "saltcorr": salt_method,
        "salt_mM": salts,
        "mismatch": flag,
        "provenance": _provenance(method="Bio.SeqUtils.MeltingTemp.Tm_GC (empirical GC formula)"),
    }


_GAF_FIELDS = (
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence",
    "With",
    "Aspect",
    "DB_Object_Name",
    "Synonym",
    "DB_Object_Type",
    "Taxon_ID",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID",
)


def write_gaf_file(
    records: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write gene-ontology annotations to a new GAF 2.0 file.

    Each record is a mapping with GAF 2.0 field names. A missing field is written as an empty
    field. A list value is written with ``|`` between its items, as in GAF. Tabs and newlines
    are rejected because they would break the table.
    """

    if not isinstance(records, Sequence) or isinstance(records, str) or not records:
        raise ConfigurationError("records must be a non-empty sequence of GAF mappings.")
    if len(records) > _MAX_RECORDS:
        raise ConfigurationError(f"records may contain at most {_MAX_RECORDS} entries.")
    destination = _new_output(output_path, "output_path")
    goa = require("Bio.UniProt.GOA", _EXTRA)
    lines: list[dict[str, Any]] = []
    for index, item in enumerate(records):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"records[{index}] must be a mapping.")
        unknown = sorted(set(item) - set(_GAF_FIELDS))
        if unknown:
            raise ConfigurationError(f"records[{index}] has unknown GAF fields: {unknown}.")
        row: dict[str, Any] = {}
        for field in _GAF_FIELDS:
            value = item.get(field, "")
            if value == []:
                value = ""
            values = value if isinstance(value, list) else [value]
            if any(
                not isinstance(part, str) or "\t" in part or "\n" in part or "|" in part
                for part in values
            ):
                raise ConfigurationError(
                    f"records[{index}].{field} must be text without tabs, newlines or '|'."
                )
            row[field] = value
        lines.append(row)
    with destination.open("x", encoding="utf-8") as handle:
        for row in lines:
            goa.writerec(row, handle, fields=_GAF_FIELDS)
    return {
        "status": "completed",
        "operation": "write_gaf_file",
        "output_path": str(destination),
        "record_count": len(lines),
        "gaf_version": "2.0",
        "provenance": _provenance(method="Bio.UniProt.GOA.writerec"),
    }


def write_motifs(
    motifs: Mapping[str, Sequence[str]],
    file_format: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format DNA motifs (equal-length sites per motif) as JASPAR, PFM, TRANSFAC or ClusterBuster.

    ``motifs`` maps each motif name to its instances, for example ``{"M1": ["ACGT", "ACGA"]}``.
    Names are sorted in the output. The ``pfm`` format holds one motif, so it rejects more.
    The result contains the formatted text and each motif's length and consensus.
    """

    fmt = _choice(file_format, _MOTIF_FORMATS, "file_format")
    if not isinstance(motifs, Mapping) or not motifs:
        raise ConfigurationError("motifs must be a non-empty mapping of name to instances.")
    if len(motifs) > _MAX_MOTIFS:
        raise ConfigurationError(f"motifs may contain at most {_MAX_MOTIFS} entries.")
    if fmt == "pfm" and len(motifs) != 1:
        raise ConfigurationError("The pfm format holds one motif; use jaspar for several.")
    motif_module = require("Bio.motifs", _EXTRA)
    built = []
    summaries = []
    for name in sorted(motifs):
        if not isinstance(name, str) or not name or any(c.isspace() for c in name):
            raise ConfigurationError("motif names must be non-empty strings without spaces.")
        instances = motifs[name]
        if not isinstance(instances, Sequence) or isinstance(instances, str) or not instances:
            raise ConfigurationError(f"motif {name!r} needs a non-empty list of instances.")
        if len(instances) > _MAX_MOTIF_INSTANCES:
            raise ConfigurationError(
                f"motif {name!r} may hold at most {_MAX_MOTIF_INSTANCES} instances."
            )
        sites = [_sequence(site, f"motif {name!r} instance", _UNAMBIGUOUS) for site in instances]
        if len({len(site) for site in sites}) != 1:
            raise ConfigurationError(f"motif {name!r} instances must share one length.")
        motif = motif_module.create(sites)
        motif.name = name
        motif.matrix_id = name
        built.append(motif)
        summaries.append({"name": name, "length": len(sites[0]), "consensus": str(motif.consensus)})
    text = motif_module.write(built, fmt)
    return {
        "status": "completed",
        "operation": "write_motifs",
        "file_format": fmt,
        "motif_count": len(built),
        "motifs": summaries,
        "text": text,
        "provenance": _provenance(method="Bio.motifs.write"),
    }


def read_mast_report(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a local MEME Suite MAST XML report: motifs, search settings and sequence diagrams.

    The diagram of a sequence lists its motif hits in order, for example ``3-[+1]-5``: gaps are
    numbers, ``[+1]`` is a forward hit of motif 1 and ``[-1]`` a reverse-complement hit.
    """

    source = input_path(path)
    mast = require("Bio.motifs.mast", _EXTRA)
    with source.open("rb") as handle:
        try:
            record = mast.read(handle)
        except (AttributeError, TypeError, ValueError) as exc:
            raise ConfigurationError(f"{source.name} is not a MAST XML report: {exc}") from exc
    motifs = [
        {
            "name": motif.name,
            "id": motif.id,
            "alt_id": motif.alt_id,
            "length": int(motif.length),
        }
        for motif in record
    ]
    sequences = [{"name": name, "diagram": record.diagrams[name]} for name in record.sequences]
    return {
        "status": "completed",
        "operation": "read_mast_report",
        "path": str(source),
        "version": record.version,
        "database": record.database,
        "alphabet": record.alphabet,
        "strand_handling": record.strand_handling,
        "motif_count": len(motifs),
        "motifs": motifs,
        "sequence_count": len(sequences),
        "sequences": sequences,
        "provenance": _provenance(method="Bio.motifs.mast.read"),
    }


__all__ = [
    "adam_consensus_tree",
    "bootstrap_branch_support",
    "codon_back_table",
    "complement_rna_sequence",
    "convert_search_results",
    "melting_temperature_gc",
    "read_fasta_qual_pairs",
    "read_mast_report",
    "reverse_complement_rna_sequence",
    "tree_to_graph_edges",
    "write_gaf_file",
    "write_motifs",
    "write_pairwise_alignment_file",
    "write_qual_file",
]

_register_public_verbose(globals(), tuple(__all__))
