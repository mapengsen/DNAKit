"""Biopython operations that DNAKit did not provide before wave 1a (package ``Bio``).

Two kinds of function live here.

* Pure-Python string operations that reproduce Biopython's string semantics exactly: letter case is
  kept, IUPAC codes are complemented, and characters outside the alphabet are passed through
  unchanged, as ``Bio.Seq`` does. The DNAKit counterparts upper-case their output, so they are not
  identical to the Biopython originals.
* Thin functions over Biopython (``Bio.Seq``, ``Bio.SeqUtils``, ``Bio.SeqIO``, ``Bio.SearchIO``,
  ``Bio.motifs``, ``Bio.Phylo``, ``Bio.PDB``). They are imported only when called, through the
  ``bio-formats`` (or ``bio-plots``) extra, so this module imports without Biopython. They keep the
  upstream defaults where the DNAKit counterpart has different defaults, a length limit, or a
  different output shape.

The two external programs, NACCESS and PSEA, are run through the bounded runner in
``dnakit.bioinformatics._common`` with fixed arguments and no shell. When the program is not
installed the function raises ``BackendUnavailableError`` with an install hint.

Plot functions write a PNG file with the matplotlib object-oriented API (no pyplot, no global
backend change) and return the file path together with the plotted data.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from io import StringIO
from itertools import accumulate
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, input_path, require, run
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_BIO_EXTRA = "bio-formats"
_PLOT_EXTRA = "bio-plots"
_MAX_SEQUENCE = 10_000_000
_MAX_TREE_TEXT = 5_000_000
_MAX_RECORDS = 100_000
_FORMAT = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,40}$")
_SEQ_TYPES = ("DNA", "RNA", "protein")
_NN_TABLES = ("DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4", "RNA_NN1", "RNA_NN2", "RNA_NN3")
_FASTA_LAYOUTS = ("wrapped", "two_line")
_NACCESS_INSTALL = (
    "Install NACCESS (http://www.bioinf.manchester.ac.uk/naccess/), put the naccess executable on "
    "PATH or pass naccess_path."
)
_PSEA_INSTALL = (
    "Install PSEA (https://github.com/mmbio/psea or the Bioconda package), put the psea executable "
    "on PATH or pass psea_path."
)
_EXTERNAL_TIMEOUT_SECONDS = 600.0

# Biopython's ``complement`` tables: IUPAC nucleotide codes in both cases; U is complemented as T.
_DNA_COMPLEMENT = str.maketrans(
    "ACGTUMRWSYKVHDBXNacgtumrwsykvhdbxn",
    "TGCAAKYWSRMBDHVXNtgcaakywsrmbdhvxn",
)
# Biopython's ``complement_rna`` table: T is complemented as U (to A), as in ``Bio.Seq``.
_RNA_COMPLEMENT = str.maketrans(
    "ACGTUMRWSYKVHDBXNacgtumrwsykvhdbxn",
    "UGCAAKYWSRMBDHVXNugcaakywsrmbdhvxn",
)
# Bio.Seq.IUPACData.ambiguous_dna_values, copied exactly: the pattern built by nt_search depends on
# these strings, so the order of letters is kept.
_IUPAC_DNA_VALUES: Mapping[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
}


def _sequence(value: object, name: str) -> str:
    """Return a bounded text sequence or raise a DNAKit configuration error."""

    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a text sequence.", code="BIO_SEQUENCE_TYPE")
    if len(value) > _MAX_SEQUENCE:
        raise ConfigurationError(
            f"{name} is longer than {_MAX_SEQUENCE} characters.", code="BIO_SEQUENCE_TOO_LONG"
        )
    return value


def _format_name(value: object, name: str) -> str:
    if not isinstance(value, str) or not _FORMAT.match(value):
        raise ConfigurationError(f"{name} must be a Biopython format name.", code="BIO_FORMAT")
    return value


def _path_text(value: str | Path) -> str:
    return str(input_path(value))


def reverse_complement_keep_case(sequence: str, *, verbose: bool = True) -> str:
    """Return the reverse complement of a DNA string as ``Bio.Seq.reverse_complement`` does.

    Letter case is kept, IUPAC codes are complemented (Y to R, M to K, ...), U is treated as T and
    any other character is left unchanged.
    """

    return _sequence(sequence, "sequence").translate(_DNA_COMPLEMENT)[::-1]


def complement_keep_case(sequence: str, *, verbose: bool = True) -> str:
    """Return the complement of a DNA string (no reversal) as ``Bio.Seq.complement`` does.

    Letter case is kept, IUPAC codes are complemented, U is treated as T and other characters are
    left unchanged.
    """

    return _sequence(sequence, "sequence").translate(_DNA_COMPLEMENT)


def complement_rna_keep_case(sequence: str, *, verbose: bool = True) -> str:
    """Return the RNA complement (A to U, T and U to A) as ``Bio.SeqUtils.complement_rna`` does.

    Letter case is kept and characters outside the IUPAC alphabet are left unchanged.
    """

    return _sequence(sequence, "sequence").translate(_RNA_COMPLEMENT)


def transcribe_keep_case(sequence: str, *, verbose: bool = True) -> str:
    """Transcribe a coding-strand DNA string into RNA by replacing T with U, as ``Bio.Seq.transcribe``.

    Only T and t are replaced; ambiguity codes and other characters are unchanged.
    """

    return _sequence(sequence, "sequence").replace("T", "U").replace("t", "u")


def back_transcribe_keep_case(sequence: str, *, verbose: bool = True) -> str:
    """Back-transcribe an RNA string into DNA by replacing U with T, as ``Bio.Seq.back_transcribe``.

    Only U and u are replaced; other characters are unchanged.
    """

    return _sequence(sequence, "sequence").replace("U", "T").replace("u", "t")


def search_nucleotide_iupac(sequence: str, subsequence: str, *, verbose: bool = True) -> list[str | int]:
    """Find every (overlapping) occurrence of a DNA pattern that may contain IUPAC codes.

    Reproduces ``Bio.SeqUtils.nt_search``: the result list starts with the regular-expression pattern
    built from ``subsequence`` (ambiguity codes become character classes, such as ``[AG]`` for R) and
    continues with the 0-based start positions of all matches on the forward strand. The search is
    case-sensitive, as in Biopython. An empty ``subsequence`` raises ``ValueError`` here; Biopython
    never returns for it.
    """

    text = _sequence(sequence, "sequence")
    query = _sequence(subsequence, "subsequence")
    if not query:
        raise ValueError("subsequence must not be empty.")
    parts: list[str] = []
    for symbol in query:
        if symbol not in _IUPAC_DNA_VALUES:
            raise ConfigurationError(
                f"{symbol!r} is not an IUPAC DNA code.", code="BIO_IUPAC_SYMBOL"
            )
        value = _IUPAC_DNA_VALUES[symbol]
        parts.append(value if len(value) == 1 else f"[{value}]")
    pattern = "".join(parts)
    result: list[str | int] = [pattern]
    position = -1
    while True:
        position += 1
        match = re.search(pattern, text[position:])
        if not match:
            break
        position += int(match.start(0))
        result.append(position)
    return result


def translate_nucleotide_bio(
    sequence: str,
    *,
    table: int | str = 1,
    stop_symbol: str = "*",
    to_stop: bool = False,
    cds: bool = False,
    gap: str | None = None,
    verbose: bool = True,
) -> str:
    """Translate a nucleotide string into protein letters with Biopython's ``Bio.Seq.translate``.

    ``table`` is an NCBI genetic-code id (default 1, the standard code) or a Biopython table name.
    Degenerate codons that code for one amino acid are translated (for example GCN gives A), other
    degenerate codons give X, and partial trailing codons are dropped with a Biopython warning.
    ``to_stop`` stops at the first stop codon, ``cds`` checks start and stop codons as for a coding
    sequence, and ``gap`` names a gap character passed through as-is. DNAKit's
    ``translate_dna_with_genetic_code`` has none of these options and raises on degenerate codons.
    """

    text = _sequence(sequence, "sequence")
    if isinstance(table, bool) or not isinstance(table, (int, str)):
        raise ConfigurationError("table must be an NCBI genetic-code id or a table name.")
    seq_module = require("Bio.Seq", _BIO_EXTRA)
    return str(
        seq_module.translate(
            text, table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap
        )
    )


def six_frame_translation_report(sequence: str, *, genetic_code: int = 1, verbose: bool = True) -> str:
    """Return the six-frame translation as the text report of ``Bio.SeqUtils.six_frame_translations``.

    The report has a header with the A, T, G and C counts, the length and the GC content, then one
    block per 60 nucleotides showing the three forward frames, the sequence, its reverse complement
    and the three reverse frames. It is produced by Biopython because the DNAKit ``translate_six_frames``
    returns a dictionary of six strings and cannot reproduce this layout.
    """

    text = _sequence(sequence, "sequence")
    if isinstance(genetic_code, bool) or not isinstance(genetic_code, int):
        raise ConfigurationError("genetic_code must be an integer NCBI table id.")
    seq_utils = require("Bio.SeqUtils", _BIO_EXTRA)
    return str(seq_utils.six_frame_translations(text, genetic_code))


def molecular_weight_biopython(
    sequence: str,
    *,
    seq_type: str = "DNA",
    double_stranded: bool = False,
    circular: bool = False,
    monoisotopic: bool = False,
    verbose: bool = True,
) -> float:
    """Return the molecular weight in daltons as ``Bio.SeqUtils.molecular_weight`` does.

    Biopython uses its own residue masses (a dNMP residue for DNA) and subtracts water for linear
    chains. ``double_stranded`` doubles the value, ``circular`` removes the terminal water and
    ``monoisotopic`` selects monoisotopic instead of average masses. DNAKit's ``molecular_weight``
    uses other residue masses and does not offer ``circular`` or ``monoisotopic``.
    """

    text = _sequence(sequence, "sequence")
    if seq_type not in _SEQ_TYPES:
        raise ConfigurationError(f"seq_type must be one of {_SEQ_TYPES}.")
    seq_utils = require("Bio.SeqUtils", _BIO_EXTRA)
    return float(
        seq_utils.molecular_weight(
            text,
            seq_type=seq_type,
            double_stranded=double_stranded,
            circular=circular,
            monoisotopic=monoisotopic,
        )
    )


def melting_temperature_biopython(
    sequence: str,
    *,
    method: str = "nearest_neighbor",
    check: bool = True,
    strict: bool = True,
    shift: int = 0,
    nn_table: str = "DNA_NN3",
    saltcorr: int = 0,
    dnac1: float = 25.0,
    dnac2: float = 25.0,
    na: float = 50.0,
    k: float = 0.0,
    tris: float = 0.0,
    mg: float = 0.0,
    dntps: float = 0.0,
    verbose: bool = True,
) -> float:
    """Return the melting temperature in Celsius with Biopython's ``Bio.SeqUtils.MeltingTemp``.

    ``method="wallace"`` calls ``Tm_Wallace`` (2 degrees per A or T, 4 per G or C; no length limit).
    ``method="nearest_neighbor"`` calls ``Tm_NN`` with the chosen nearest-neighbour table (default
    ``DNA_NN3``), strand concentrations ``dnac1``/``dnac2`` (nM), monovalent salt ``na``/``k``/``tris``
    (mM), divalent ``mg`` and ``dntps`` (mM) and salt correction ``saltcorr`` (0 means none, as in
    Biopython). DNAKit's ``melting_temperature`` uses the SantaLucia 1998 table and its own salt model.
    """

    text = _sequence(sequence, "sequence")
    if method not in ("wallace", "nearest_neighbor"):
        raise ConfigurationError("method must be 'wallace' or 'nearest_neighbor'.")
    if nn_table not in _NN_TABLES:
        raise ConfigurationError(f"nn_table must be one of {_NN_TABLES}.")
    mt = require("Bio.SeqUtils.MeltingTemp", _BIO_EXTRA)
    if method == "wallace":
        return float(mt.Tm_Wallace(text, check=check, strict=strict))
    return float(
        mt.Tm_NN(
            text,
            check=check,
            strict=strict,
            shift=shift,
            nn_table=getattr(mt, nn_table),
            saltcorr=saltcorr,
            dnac1=dnac1,
            dnac2=dnac2,
            Na=na,
            K=k,
            Tris=tris,
            Mg=mg,
            dNTPs=dntps,
        )
    )


def gc_fraction_biopython(sequence: str, *, ambiguous: str = "remove", verbose: bool = True) -> float:
    """Return the G+C fraction (0 to 1) as ``Bio.SeqUtils.gc_fraction`` does.

    ``ambiguous`` is ``"remove"`` (ambiguous symbols are not counted in the denominator), ``"ignore"``
    (they count in the denominator as non-GC) or ``"weighted"`` (partly-GC codes such as S count as
    fractions). DNAKit's ``gc_at_content`` offers only ``error`` and ``ignore`` policies.
    """

    text = _sequence(sequence, "sequence")
    if ambiguous not in ("remove", "ignore", "weighted"):
        raise ConfigurationError("ambiguous must be 'remove', 'ignore' or 'weighted'.")
    seq_utils = require("Bio.SeqUtils", _BIO_EXTRA)
    return float(seq_utils.gc_fraction(text, ambiguous=ambiguous))


def codon_dn_ds_biopython(
    sequence_a: str,
    sequence_b: str,
    *,
    method: str = "NG86",
    codon_table: int = 1,
    k: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return unrounded dN and dS between two codon-aligned coding sequences (Biopython).

    ``method`` is ``"NG86"`` (Nei and Gojobori 1986, as ``Bio.codonalign.codonseq.cal_dn_ds``) or
    ``"LWL85"`` (Li, Wu and Luo 1985, with the transition/transversion ratio ``k``). ``codon_table``
    is an NCBI id of an unambiguous genetic code. The sequences must have the same length, a
    multiple of three, and no stop codons. DNAKit's ``calculate_dn_ds`` supports NG86 only and
    reports values rounded to four decimals.
    """

    first = _sequence(sequence_a, "sequence_a")
    second = _sequence(sequence_b, "sequence_b")
    if method not in ("NG86", "LWL85"):
        raise ConfigurationError("method must be 'NG86' or 'LWL85'.")
    if isinstance(codon_table, bool) or not isinstance(codon_table, int):
        raise ConfigurationError("codon_table must be an NCBI genetic-code id.")
    if isinstance(k, bool) or not isinstance(k, (int, float)) or k <= 0:
        raise ConfigurationError("k must be a positive number.")
    codonseq = require("Bio.codonalign.codonseq", _BIO_EXTRA)
    codon_tables = require("Bio.Data.CodonTable", _BIO_EXTRA)
    if codon_table not in codon_tables.unambiguous_dna_by_id:
        raise ConfigurationError(f"codon_table {codon_table} is not a known unambiguous table.")
    try:
        dn, ds = codonseq.cal_dn_ds(
            codonseq.CodonSeq(first),
            codonseq.CodonSeq(second),
            method=method,
            codon_table=codon_tables.unambiguous_dna_by_id[codon_table],
            k=float(k),
        )
    except (ValueError, ZeroDivisionError, codon_tables.TranslationError) as exc:
        raise ConfigurationError(f"dN/dS could not be computed: {exc}") from exc
    dn = float(dn)
    ds = float(ds)
    return {
        "dN": dn,
        "dS": ds,
        "dN_dS": (dn / ds) if ds > 0 else None,
        "method": method,
        "codon_table": codon_table,
        "codon_count": len(first) // 3,
    }


def parse_fasta_title_sequence(
    path: str | Path, *, layout: str = "wrapped", verbose: bool = True
) -> list[tuple[str, str]]:
    """Read FASTA records as ``(title, sequence)`` pairs with Biopython's simple parsers.

    ``layout="wrapped"`` uses ``SimpleFastaParser`` (sequence lines may be wrapped; the title is the
    header text after ``>``). ``layout="two_line"`` uses ``FastaTwoLineParser`` (each record is
    exactly one header line and one sequence line). Sequence text is returned as written.
    """

    if layout not in _FASTA_LAYOUTS:
        raise ConfigurationError(f"layout must be one of {_FASTA_LAYOUTS}.")
    fasta = require("Bio.SeqIO.FastaIO", _BIO_EXTRA)
    parser = fasta.SimpleFastaParser if layout == "wrapped" else fasta.FastaTwoLineParser
    records: list[tuple[str, str]] = []
    with input_path(path).open(encoding="utf-8") as handle:
        for title, sequence in parser(handle):
            records.append((title, sequence))
            if len(records) > _MAX_RECORDS:
                raise ConfigurationError("FASTA file has too many records.", code="BIO_TOO_MANY")
    return records


def iterate_fastq_tuples(path: str | Path, *, verbose: bool = True) -> list[tuple[str, str, str]]:
    """Read FASTQ records as ``(title, sequence, quality)`` string triples.

    Uses Biopython's ``FastqGeneralIterator``, so no SeqRecord objects are built and quality strings
    are returned as Phred+33 text exactly as in the file.
    """

    quality_io = require("Bio.SeqIO.QualityIO", _BIO_EXTRA)
    records: list[tuple[str, str, str]] = []
    with input_path(path).open(encoding="utf-8") as handle:
        for title, sequence, quality in quality_io.FastqGeneralIterator(handle):
            records.append((title, sequence, quality))
            if len(records) > _MAX_RECORDS:
                raise ConfigurationError("FASTQ file has too many records.", code="BIO_TOO_MANY")
    return records


def _record_summary(record: Any) -> dict[str, Any]:
    sequence = str(record.seq)
    return {
        "id": record.id,
        "description": record.description,
        "length": len(sequence),
        "sequence": sequence,
    }


def read_single_sequence_record(
    path: str | Path, format: str = "fasta", *, verbose: bool = True
) -> dict[str, Any]:
    """Read a file that holds exactly one sequence record, as ``Bio.SeqIO.read`` does.

    Returns the id, description, length and sequence. Biopython raises ``ValueError`` when the file
    has no record or more than one; this function raises ``ConfigurationError`` in the same cases.
    Any Biopython sequence format name is accepted.
    """

    _format_name(format, "format")
    seq_io = require("Bio.SeqIO", _BIO_EXTRA)
    try:
        with input_path(path).open(encoding="utf-8") as handle:
            record = seq_io.read(handle, format)
    except ValueError as exc:
        raise ConfigurationError(f"{Path(path).name} does not hold exactly one record: {exc}") from exc
    return _record_summary(record)


def sequence_records_by_id(
    path: str | Path, format: str = "fasta", *, verbose: bool = True
) -> dict[str, dict[str, Any]]:
    """Read sequence records into a dictionary keyed by record id, as ``Bio.SeqIO.to_dict`` does.

    Duplicate ids raise ``ConfigurationError`` (Biopython raises ``ValueError``).
    """

    _format_name(format, "format")
    seq_io = require("Bio.SeqIO", _BIO_EXTRA)
    try:
        with input_path(path).open(encoding="utf-8") as handle:
            mapping = seq_io.to_dict(seq_io.parse(handle, format))
    except ValueError as exc:
        raise ConfigurationError(f"{Path(path).name} could not be indexed by id: {exc}") from exc
    return {key: _record_summary(record) for key, record in mapping.items()}


def _hsp_summary(hsp: Any) -> dict[str, Any]:
    return {
        "query_range": [int(hsp.query_range[0]), int(hsp.query_range[1])],
        "hit_range": [int(hsp.hit_range[0]), int(hsp.hit_range[1])],
        "evalue": float(hsp.evalue),
        "bitscore": float(hsp.bitscore),
        "identities": int(getattr(hsp, "ident_num", 0) or 0),
    }


def _query_summary(query: Any) -> dict[str, Any]:
    hits = []
    for hit in query.hits:
        hits.append(
            {
                "id": hit.id,
                "description": hit.description,
                "hsps": [_hsp_summary(hsp) for hsp in hit.hsps],
            }
        )
    return {
        "id": query.id,
        "description": query.description,
        "length": int(query.seq_len),
        "hit_count": len(hits),
        "hits": hits,
    }


def read_single_search_result(
    path: str | Path, format: str = "blast-xml", *, verbose: bool = True
) -> dict[str, Any]:
    """Read a search-output file with exactly one query, as ``Bio.SearchIO.read`` does.

    The result is a JSON-serialisable summary of the query: its id, length and hits, each with its
    HSPs (query and hit ranges, E-value, bit score, identities). Files with zero or several queries
    raise ``ConfigurationError``, as ``SearchIO.read`` raises ``ValueError``.
    """

    _format_name(format, "format")
    searchio = require("Bio.SearchIO", _BIO_EXTRA)
    try:
        with input_path(path).open("rb") as handle:
            query = searchio.read(handle, format)
    except ValueError as exc:
        raise ConfigurationError(f"{Path(path).name} does not hold exactly one query: {exc}") from exc
    return _query_summary(query)


def search_results_by_id(
    path: str | Path, format: str = "blast-xml", *, verbose: bool = True
) -> dict[str, dict[str, Any]]:
    """Read every query of a search-output file into a dictionary keyed by query id.

    Each value has the summary of :func:`read_single_search_result`. Duplicate query ids raise
    ``ConfigurationError``, as ``SearchIO.to_dict`` does.
    """

    _format_name(format, "format")
    searchio = require("Bio.SearchIO", _BIO_EXTRA)
    try:
        with input_path(path).open("rb") as handle:
            mapping = searchio.to_dict(searchio.parse(handle, format))
    except ValueError as exc:
        raise ConfigurationError(f"{Path(path).name} could not be indexed by id: {exc}") from exc
    return {key: _query_summary(query) for key, query in mapping.items()}


def read_single_motif(path: str | Path, format: str = "jaspar", *, verbose: bool = True) -> dict[str, Any]:
    """Read a motif file that holds exactly one motif, as ``Bio.motifs.read`` does.

    Returns the motif name, length, consensus and the per-position base counts. Biopython raises
    ``ValueError`` for zero or several motifs; this function raises ``ConfigurationError``.
    """

    _format_name(format, "format")
    motifs = require("Bio.motifs", _BIO_EXTRA)
    try:
        with input_path(path).open(encoding="utf-8") as handle:
            motif = motifs.read(handle, format)
    except ValueError as exc:
        raise ConfigurationError(f"{Path(path).name} does not hold exactly one motif: {exc}") from exc
    counts = {letter: [int(value) for value in motif.counts[letter]] for letter in "ACGT"}
    return {
        "name": motif.name,
        "length": int(motif.length),
        "consensus": str(motif.consensus),
        "counts": counts,
    }


def draw_ascii_tree(newick: str, *, column_width: int = 80, verbose: bool = True) -> dict[str, Any]:
    """Draw a Newick tree as ASCII art with ``Bio.Phylo.draw_ascii`` and return the text and the tips.

    ``column_width`` is the width the drawing may use (Biopython default 80). The result holds the
    drawing text, its line count and the tip names in tree order.
    """

    text = _sequence(newick, "newick")
    if len(text) > _MAX_TREE_TEXT:
        raise ConfigurationError("Newick text is too large.", code="BIO_TREE_TOO_LARGE")
    if isinstance(column_width, bool) or not isinstance(column_width, int) or column_width < 1:
        raise ConfigurationError("column_width must be a positive integer.")
    phylo = require("Bio.Phylo", _BIO_EXTRA)
    tree = phylo.read(StringIO(text), "newick")
    buffer = StringIO()
    phylo.draw_ascii(tree, file=buffer, column_width=column_width)
    drawing = buffer.getvalue()
    return {
        "text": drawing,
        "line_count": len(drawing.splitlines()),
        "tip_names": [str(tip.name) for tip in tree.get_terminals()],
        "column_width": column_width,
    }


def draw_tree_figure(
    newick: str,
    image_path: str | Path,
    *,
    show_confidence: bool = True,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a Newick tree with ``Bio.Phylo.draw`` into a PNG file and return the path and the tips.

    The figure is drawn with the matplotlib object API on a private canvas, so the global pyplot
    state is not changed. The result holds the PNG path, the tip names and the clade count.
    """

    text = _sequence(newick, "newick")
    if len(text) > _MAX_TREE_TEXT:
        raise ConfigurationError("Newick text is too large.", code="BIO_TREE_TOO_LARGE")
    phylo = require("Bio.Phylo", _BIO_EXTRA)
    figure_module = require("matplotlib.figure", _PLOT_EXTRA)
    canvas_module = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    tree = phylo.read(StringIO(text), "newick")
    figure = figure_module.Figure(figsize=(6, 4))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(1, 1, 1)
    phylo.draw(tree, do_show=False, axes=axes, show_confidence=show_confidence)
    target = Path(image_path)
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, format="png", dpi=dpi)
    return {
        "image_path": str(target.resolve()),
        "format": "png",
        "tip_names": [str(tip.name) for tip in tree.get_terminals()],
        "clade_count": sum(1 for _ in tree.find_clades()),
    }


def gc_skew_figure(
    sequence: str,
    image_path: str | Path,
    *,
    window: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the GC skew per window and the accumulated skew, and save a line plot as PNG.

    The skew values come from ``Bio.SeqUtils.GC_skew(sequence, window)``; the accumulated skew is
    their running sum, as in ``Bio.SeqUtils.xGC_skew``. Biopython draws the circular map on a Tk
    canvas (it needs a display); this function writes a rectangular line plot instead, and returns
    the numbers so that they can be checked without the image.
    """

    text = _sequence(sequence, "sequence")
    if isinstance(window, bool) or not isinstance(window, int) or window < 1:
        raise ConfigurationError("window must be a positive integer.")
    seq_utils = require("Bio.SeqUtils", _BIO_EXTRA)
    figure_module = require("matplotlib.figure", _PLOT_EXTRA)
    canvas_module = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    skew = [float(value) for value in seq_utils.GC_skew(text, window)]
    accumulated = list(accumulate(skew))
    target = Path(image_path)
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    target.parent.mkdir(parents=True, exist_ok=True)
    figure = figure_module.Figure(figsize=(7, 3.5))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(1, 1, 1)
    positions = list(range(len(skew)))
    axes.plot(positions, skew, color="tab:blue", label="GC skew")
    axes.plot(positions, accumulated, color="tab:purple", label="Accumulated GC skew")
    axes.set_xlabel(f"window index (window = {window} nt)")
    axes.set_ylabel("skew")
    axes.legend()
    figure.tight_layout()
    figure.savefig(target, format="png", dpi=150)
    return {
        "image_path": str(target.resolve()),
        "format": "png",
        "window": window,
        "window_count": len(skew),
        "gc_skew": skew,
        "accumulated_gc_skew": accumulated,
    }


def _external_executable(name: str, requested: str | None, hint: str) -> str:
    program = requested or name
    resolved = shutil.which(program)
    if resolved is None:
        raise BackendUnavailableError(f"External program {program!r} was not found.", hint=hint)
    return str(Path(resolved).resolve())


def run_naccess_accessibility(
    pdb_path: str | Path,
    *,
    probe_size: float | None = None,
    z_slice: float | None = None,
    naccess_path: str | None = None,
    work_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the NACCESS program on a PDB file and return the per-residue solvent accessibility.

    The file is copied into a working directory under a ``.pdb`` name (NACCESS requires it) and run
    as ``naccess <file> [-p probe_size] [-z z_slice]``. The ``.rsa`` (relative) and ``.asa`` (absolute)
    output files are read with Biopython's ``process_rsa_data``. NACCESS is not bundled: if it is not
    installed a ``BackendUnavailableError`` with an install hint is raised. The output is written under
    ``work_dir`` (a new temporary directory when not given).
    """

    source = input_path(pdb_path)
    if not source.is_file():
        raise ConfigurationError(f"PDB file not found: {source.name}")
    program = _external_executable("naccess", naccess_path, _NACCESS_INSTALL)
    if probe_size is not None and not 0 < float(probe_size) <= 10:
        raise ConfigurationError("probe_size must be in (0, 10] angstrom.")
    if z_slice is not None and not 0 < float(z_slice) <= 10:
        raise ConfigurationError("z_slice must be in (0, 10] angstrom.")
    directory = Path(work_dir) if work_dir is not None else Path(tempfile.mkdtemp(prefix="naccess_"))
    directory.mkdir(parents=True, exist_ok=True)
    model = directory / "model.pdb"
    model.write_bytes(source.read_bytes())
    arguments = [model.name]
    if probe_size is not None:
        arguments += ["-p", str(float(probe_size))]
    if z_slice is not None:
        arguments += ["-z", str(float(z_slice))]
    rsa_path = directory / "model.rsa"
    asa_path = directory / "model.asa"
    log = run(
        program,
        arguments,
        directory,
        ExecutionConfig(timeout_seconds=_EXTERNAL_TIMEOUT_SECONDS),
        outputs=(rsa_path, asa_path),
    )
    naccess = require("Bio.PDB.NACCESS", _BIO_EXTRA)
    rsa_lines = rsa_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    residues = []
    for key, values in naccess.process_rsa_data(rsa_lines).items():
        residues.append({"residue": str(key), "values": [float(value) for value in values]})
    return {
        "status": "completed",
        "operation": "run_naccess_accessibility",
        "backend": "naccess",
        "program": program,
        "rsa_path": str(rsa_path),
        "asa_path": str(asa_path),
        "residue_count": len(residues),
        "residues": residues,
        "log_tail": log[-4000:],
    }


def _parse_psea_secondary_structure(path: Path) -> str:
    """Return the secondary-structure string from the ``>p-sea`` block of a PSEA ``.sea`` file."""

    structure = ""
    started = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[0:6] == ">p-sea":
                started = True
                continue
            if not started:
                continue
            if line[0] == "\n":
                break
            structure += line[0:-1]
    return structure


def run_psea_secondary_structure(
    pdb_path: str | Path,
    *,
    psea_path: str | None = None,
    work_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the PSEA program on a PDB file and return its secondary-structure string.

    PSEA is run as ``psea <file>`` in the working directory and writes ``<name>.sea``; the
    secondary-structure block of that file is returned, one character per residue (H helix, E strand,
    C coil), as ``Bio.PDB.PSEA.psea`` does. A missing program raises ``BackendUnavailableError``.
    """

    source = input_path(pdb_path)
    if not source.is_file():
        raise ConfigurationError(f"PDB file not found: {source.name}")
    program = _external_executable("psea", psea_path, _PSEA_INSTALL)
    directory = Path(work_dir) if work_dir is not None else Path(tempfile.mkdtemp(prefix="psea_"))
    directory.mkdir(parents=True, exist_ok=True)
    stem = source.stem
    copied = directory / f"{stem}.pdb"
    copied.write_bytes(source.read_bytes())
    sea_path = directory / f"{stem}.sea"
    log = run(
        program,
        [copied.name],
        directory,
        ExecutionConfig(timeout_seconds=_EXTERNAL_TIMEOUT_SECONDS),
        outputs=(sea_path,),
    )
    structure = _parse_psea_secondary_structure(sea_path)
    if not structure:
        raise BackendExecutionError("psea produced no secondary-structure block.", context={"log": log})
    return {
        "status": "completed",
        "operation": "run_psea_secondary_structure",
        "backend": "psea",
        "program": program,
        "sea_path": str(sea_path),
        "length": len(structure),
        "secondary_structure": structure,
        "log_tail": log[-4000:],
    }


__all__ = [
    "back_transcribe_keep_case",
    "complement_keep_case",
    "codon_dn_ds_biopython",
    "complement_rna_keep_case",
    "draw_ascii_tree",
    "draw_tree_figure",
    "gc_fraction_biopython",
    "gc_skew_figure",
    "iterate_fastq_tuples",
    "melting_temperature_biopython",
    "molecular_weight_biopython",
    "parse_fasta_title_sequence",
    "read_single_motif",
    "read_single_search_result",
    "read_single_sequence_record",
    "reverse_complement_keep_case",
    "run_naccess_accessibility",
    "run_psea_secondary_structure",
    "search_nucleotide_iupac",
    "search_results_by_id",
    "sequence_records_by_id",
    "six_frame_translation_report",
    "transcribe_keep_case",
    "translate_nucleotide_bio",
]

_register_public_verbose(globals(), tuple(__all__))
