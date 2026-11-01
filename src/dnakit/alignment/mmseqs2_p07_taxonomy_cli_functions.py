"""MMseqs2 taxonomy assignment and multi-hit search.

Command-line facades over MMseqs2 (integration pack P07). Each public function runs one fixed
MMseqs2 module with an explicit argument list (no shell) through
``dnakit.backends.external.execute_bounded_command``. The user installs MMseqs2 (conda or a release
binary) and puts ``mmseqs`` on PATH or passes ``mmseqs_path``; the tested version is the build with
commit 76de5d8ae114750adc59288516d947b6a72a4b05. MMseqs2 databases are addressed by path prefix;
options that are not given are not passed, so MMseqs2's own defaults apply. Results are
JSON-serialisable dictionaries with status, version, command, outputs, entry or line counts and a
bounded log tail.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment import _mmseqs2_p07_core as core

_TAXONOMY_SUFFIXES = (
    "_mapping",
    "_taxonomy",
    "_nodes.dmp",
    "_names.dmp",
    "_merged.dmp",
    "_delnodes.dmp",
)


def mmseqs_taxonomy(
    query: str | Sequence[str],
    target_db: str,
    output: str,
    *,
    input_kind: str = "fasta",
    lca_mode: int | None = None,
    majority: float | None = None,
    vote_mode: int | None = None,
    lca_ranks: str | None = None,
    tax_lineage: int | None = None,
    blacklist: str | None = None,
    taxon_list: str | None = None,
    orf_filter: bool | None = None,
    tax_output_mode: int | None = None,
    sensitivity: float | None = None,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    max_seqs: int | None = None,
    kmer_size: int | None = None,
    search_type: int | None = None,
    min_length: int | None = None,
    translation_table: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Taxonomic classification of sequences (``mmseqs easy-taxonomy``, ``mmseqs taxonomy``).

    The classification searches the target (protein, or nucleotide contigs translated on the fly),
    and assigns each query the lowest common ancestor of its top hits. ``target_db`` must carry
    taxonomy (see ``mmseqs_createtaxdb``). ``input_kind="fasta"`` runs ``easy-taxonomy``: FASTA
    query file(s) and a new file prefix as ``output``; ``<output>_lca.tsv``, ``<output>_report``,
    ``<output>_tophit_aln`` and ``<output>_tophit_report`` are written.
    ``input_kind="database"`` runs
    ``taxonomy``: a query database and the prefix of a new taxonomy result database.

    Parameters:
    - query, target_db, output, input_kind: see above.
    - lca_mode: 1 single-search LCA, 2 or 3 approximate 2bLCA, 4 top hit (default 3).
    - majority, vote_mode, lca_ranks, tax_lineage, blacklist, taxon_list: LCA and output options
      (majority default 0.5; vote_mode 0 uniform, 1 minus log E-value, 2 score; tax_lineage 0 none,
      1 names, 2 taxon ids).
    - orf_filter: pre-filter query ORFs with a non-selective search (nucleotide queries,
      default on).
    - tax_output_mode: database mode only; 0 LCA, 1 alignment, 2 both.
    - sensitivity, e_value, min_seq_id, coverage, cov_mode, max_seqs, kmer_size, search_type:
      search options.
    - min_length, translation_table: ORF extraction for nucleotide queries.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.

    The fasta mode adds ``records`` (up to 1000 rows of ``_lca.tsv`` as query, taxon_id, rank,
    name).
    """
    kind = core.choose_kind(input_kind)
    options: list[str] = []
    options += core.opt_choice("--lca-mode", lca_mode, (1, 2, 3, 4), "lca_mode")
    options += core.opt_float("--majority", majority, "majority", 0.0, 1.0)
    options += core.opt_choice("--vote-mode", vote_mode, (0, 1, 2), "vote_mode")
    options += core.opt_text("--lca-ranks", lca_ranks, "lca_ranks", r"^[A-Za-z_ ]+(,[A-Za-z_ ]+)*$")
    options += core.opt_choice("--tax-lineage", tax_lineage, (0, 1, 2), "tax_lineage")
    options += core.opt_text("--blacklist", blacklist, "blacklist", r"^[0-9A-Za-z_:, .]+$")
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    options += core.opt_bool("--orf-filter", orf_filter, "orf_filter")
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_int("--min-length", min_length, "min_length", 1, None)
    options += core.opt_int("--translation-table", translation_table, "translation_table", 1, 31)
    target_path = core.input_database(target_db, "target_db")
    if kind == "database":
        options += core.opt_choice(
            "--tax-output-mode", tax_output_mode, (0, 1, 2), "tax_output_mode"
        )
        query_path = core.input_database(query if isinstance(query, str) else "", "query")
        output_path = core.output_database(output, "output")
        return core.run(
            "mmseqs_taxonomy",
            "taxonomy",
            [query_path, target_path, output_path],
            options,
            databases={"taxonomy_db": output_path},
            needs_tmp=True,
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    core.reject_with(
        tax_output_mode is not None, "tax_output_mode applies to input_kind='database' only."
    )
    inputs = core.fasta_files(query, "query")
    prefix, family = core.output_family(
        output,
        {
            "lca": "_lca.tsv",
            "report": "_report",
            "top_hit_alignment": "_tophit_aln",
            "top_hit_report": "_tophit_report",
        },
        "output",
    )
    result = core.run(
        "mmseqs_taxonomy",
        "easy-taxonomy",
        [*inputs, target_path, prefix],
        options,
        files=family,
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["records"] = core.table(family["lca"], ("query", "taxon_id", "rank", "name"))
    return result


def mmseqs_createtaxdb(
    sequence_db: str,
    *,
    ncbi_tax_dump: str | None = None,
    tax_mapping_file: str | None = None,
    tax_mapping_mode: int | None = None,
    tax_db_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate a sequence database with taxonomy labels (``mmseqs createtaxdb``).

    The taxonomy files are written next to the database and the database is then usable for
    taxonomic classification: ``<sequence_db>_mapping`` (sequence to taxon id) and, in the
    default binary
    mode, ``<sequence_db>_taxonomy`` (or ``_nodes.dmp``, ``_names.dmp``, ``_merged.dmp`` with
    tax_db_mode 0).
    Existing taxonomy files are never overwritten. The taxonomy comes from an NCBI-style
    taxonomy dump
    (nodes.dmp, names.dmp, merged.dmp) and a mapping from sequence identifiers to taxon ids.

    Parameters:
    - sequence_db: prefix of the sequence database (modified in place by adding the taxonomy files).
    - ncbi_tax_dump: directory with the NCBI taxonomy dump.
    - tax_mapping_file: tab-separated file, sequence identifier then taxon id.
    - tax_mapping_mode: 0 map through the .lookup file, 1 through the .source file (default 0).
    - tax_db_mode: 0 write .dmp flat files, 1 binary dump (default 1).
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.
    """
    database = core.input_database(sequence_db, "sequence_db")
    mapping = Path(f"{database}_mapping")
    existing = [
        str(Path(f"{database}{suffix}"))
        for suffix in _TAXONOMY_SUFFIXES
        if Path(f"{database}{suffix}").exists()
    ]
    core.reject_with(bool(existing), f"The database already carries taxonomy files: {existing}")
    options: list[str] = []
    options += core.opt_directory("--ncbi-tax-dump", ncbi_tax_dump, "ncbi_tax_dump")
    options += core.opt_file("--tax-mapping-file", tax_mapping_file, "tax_mapping_file")
    options += core.opt_choice("--tax-mapping-mode", tax_mapping_mode, (0, 1), "tax_mapping_mode")
    options += core.opt_choice("--tax-db-mode", tax_db_mode, (0, 1), "tax_db_mode")
    result = core.run(
        "mmseqs_createtaxdb",
        "createtaxdb",
        [database],
        options,
        extra_outputs={"mapping": mapping},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["taxonomy_files"] = [
        str(Path(f"{database}{suffix}"))
        for suffix in _TAXONOMY_SUFFIXES
        if Path(f"{database}{suffix}").exists()
    ]
    return result


def mmseqs_aggregate_taxonomy(
    sequence_db: str,
    set_to_sequence_db: str,
    taxonomy_per_sequence_db: str,
    output: str,
    *,
    alignment_per_sequence_db: str | None = None,
    swap_headers: bool = False,
    majority: float | None = None,
    vote_mode: int | None = None,
    lca_ranks: str | None = None,
    tax_lineage: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Aggregate taxon labels of a set of sequences into one (``mmseqs aggregatetax``,
    ``aggregatetaxweights``).

    A set is, for example, the ORFs of one contig. Without ``alignment_per_sequence_db`` the labels
    are combined by their lowest common ancestor (``aggregatetax``). With it, the labels are voted
    with weights from the alignment E-value or score and must reach the ``majority`` fraction
    (``aggregatetaxweights``); the alignments come from
    ``mmseqs_taxonomy(..., input_kind="database", tax_output_mode=2)``.

    Parameters:
    - sequence_db: taxonomy-annotated target sequence database.
    - set_to_sequence_db: database that maps every set key to its member keys. For ORFs of
      contigs it
      is the transposed header database of the ORF database; pass ``<orf_db>_h`` with
      ``swap_headers=True`` and the facade transposes it (``mmseqs swapdb``) first.
    - taxonomy_per_sequence_db: taxonomy result of the sequences.
    - output: prefix of the new taxonomy result database per set.
    - alignment_per_sequence_db: alignment result of the sequences; selects the weighted variant.
    - swap_headers: transpose ``set_to_sequence_db`` before use (see above); default False.
    - majority: weighted variant only; minimal fraction of agreement among the labelled sequences of
      a set (MMseqs2 default 0.5).
    - vote_mode: weighted variant only; 0 uniform, 1 minus log E-value, 2 score (default 1).
    - lca_ranks, tax_lineage: extra output columns as in mmseqs_majoritylca.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.
    """
    target = core.input_database(sequence_db, "sequence_db")
    set_map = core.input_database(set_to_sequence_db, "set_to_sequence_db")
    per_sequence = core.input_database(taxonomy_per_sequence_db, "taxonomy_per_sequence_db")
    output_path = core.output_database(output, "output")
    options: list[str] = []
    options += core.opt_text("--lca-ranks", lca_ranks, "lca_ranks", r"^[A-Za-z_ ]+(,[A-Za-z_ ]+)*$")
    options += core.opt_choice("--tax-lineage", tax_lineage, (0, 1, 2), "tax_lineage")
    if alignment_per_sequence_db is None:
        core.reject_with(
            majority is not None or vote_mode is not None,
            "majority and vote_mode need alignment_per_sequence_db (aggregatetaxweights).",
        )
        command = "aggregatetax"
        alignment_path = None
    else:
        command = "aggregatetaxweights"
        alignment_path = core.input_database(alignment_per_sequence_db, "alignment_per_sequence_db")
        options += core.opt_float("--majority", majority, "majority", 0.0, 1.0)
        options += core.opt_choice("--vote-mode", vote_mode, (0, 1, 2), "vote_mode")
    scratch = (
        tempfile.mkdtemp(prefix=".mmseqs_swap_", dir=output_path.parent) if swap_headers else None
    )
    try:
        if scratch is not None:
            swapped = Path(scratch) / "set_to_sequence"
            core.run(
                "mmseqs_aggregate_taxonomy",
                "swapdb",
                [set_map, swapped],
                [],
                databases={"swapped": swapped},
                threads=threads,
                mmseqs_path=mmseqs_path,
                verbose=False,
            )
            # swapdb writes a generic type; the aggregation modules expect a taxonomy-map type.
            Path(f"{swapped}.dbtype").write_bytes(b"\x05\x00\x00\x00")
            set_map = swapped
        positionals = [target, set_map, per_sequence]
        if alignment_path is not None:
            positionals.append(alignment_path)
        positionals.append(output_path)
        return core.run(
            "mmseqs_aggregate_taxonomy",
            command,
            positionals,
            options,
            databases={"taxonomy_db": output_path},
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    finally:
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)


def mmseqs_addtaxonomy(
    target_db: str,
    result_db: str,
    annotated_db: str,
    *,
    tax_lineage: int | None = None,
    lca_ranks: str | None = None,
    pick_id_from: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add taxonomic labels to a result database (`mmseqs addtaxonomy`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - target_db: taxonomy-annotated target sequence database.
    - result_db: result (alignment) database.
    - annotated_db: new result database with taxon columns.
    - tax_lineage: Lineage column: 0 none, 1 all lineage names, 2 all lineage taxon ids
      (default 0).
    - lca_ranks: Add a column with these ranks, comma separated, e.g. 'phylum,genus'.
    - pick_id_from: Take the taxon id from 1 the query or 2 the target (default 2).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    target_db_path = core.input_database(target_db, "target_db")
    result_db_path = core.input_database(result_db, "result_db")
    annotated_db_path = core.output_database(annotated_db, "annotated_db")
    positionals: list[Path] = []
    positionals.append(target_db_path)
    positionals.append(result_db_path)
    positionals.append(annotated_db_path)
    options: list[str] = []
    options += core.opt_choice("--tax-lineage", tax_lineage, (0, 1, 2), "tax_lineage")
    options += core.opt_text("--lca-ranks", lca_ranks, "lca_ranks", "^[A-Za-z_ ]+(,[A-Za-z_ ]+)*$")
    options += core.opt_choice("--pick-id-from", pick_id_from, (1, 2), "pick_id_from")
    result = core.run(
        "mmseqs_addtaxonomy",
        "addtaxonomy",
        positionals,
        options,
        databases={"annotated_db": annotated_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_taxonomyreport(
    sequence_db: str,
    result_db: str,
    output: str,
    *,
    report_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a Kraken-style or Krona taxonomy report from a taxonomy result (`mmseqs
    taxonomyreport`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: taxonomy-annotated sequence database.
    - result_db: taxonomy result (or any result/sequence) database.
    - output: new report file.
    - report_mode: Report format: 0 Kraken, 1 Krona (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    The result adds ``rows`` (up to 1000 report rows) for the Kraken format.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    result_db_path = core.input_database(result_db, "result_db")
    output_path = core.output_file(output, "output")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(result_db_path)
    positionals.append(output_path)
    options: list[str] = []
    options += core.opt_choice("--report-mode", report_mode, (0, 1), "report_mode")
    result = core.run(
        "mmseqs_taxonomyreport",
        "taxonomyreport",
        positionals,
        options,
        databases={},
        files={"output": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    if report_mode in (None, 0):
        result["rows"] = core.table(output_path)
    return result


def mmseqs_filtertaxdb(
    target_db: str,
    taxonomy_db: str,
    filtered_db: str,
    *,
    taxon_list: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the entries of a taxonomy result whose taxon matches a taxon expression (`mmseqs
    filtertaxdb`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - target_db: taxonomy-annotated target database.
    - taxonomy_db: taxonomy result database.
    - filtered_db: new taxonomy result database.
    - taxon_list: Taxon ids to keep, e.g. '9606||810' (or), '2759&&!4751' (and not); only
      digits, ',', '&', '|', '!' and brackets.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    target_db_path = core.input_database(target_db, "target_db")
    taxonomy_db_path = core.input_database(taxonomy_db, "taxonomy_db")
    filtered_db_path = core.output_database(filtered_db, "filtered_db")
    positionals: list[Path] = []
    positionals.append(target_db_path)
    positionals.append(taxonomy_db_path)
    positionals.append(filtered_db_path)
    options: list[str] = []
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    result = core.run(
        "mmseqs_filtertaxdb",
        "filtertaxdb",
        positionals,
        options,
        databases={"filtered_db": filtered_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_filtertaxseqdb(
    sequence_db: str,
    filtered_db: str,
    *,
    taxon_list: str | None = None,
    subdb_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the sequences of a taxonomy-annotated database that belong to a taxon expression
    (`mmseqs filtertaxseqdb`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: taxonomy-annotated sequence database.
    - filtered_db: new sequence database.
    - taxon_list: Taxon ids to keep, e.g. '9606||810' (or), '2759&&!4751' (and not); only
      digits, ',', '&', '|', '!' and brackets.
    - subdb_mode: 0 copy the data, 1 soft-link the data and write a new index (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    filtered_db_path = core.output_database(filtered_db, "filtered_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(filtered_db_path)
    options: list[str] = []
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    options += core.opt_choice("--subdb-mode", subdb_mode, (0, 1), "subdb_mode")
    result = core.run(
        "mmseqs_filtertaxseqdb",
        "filtertaxseqdb",
        positionals,
        options,
        databases={"filtered_db": filtered_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_lcaalign(
    query_db: str,
    target_db: str,
    result_db: str,
    alignment_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    add_backtrace: bool | None = None,
    alignment_mode: int | None = None,
    alignment_output_mode: int | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gapped alignment tailored to the lowest-common-ancestor computation (`mmseqs lcaalign`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - result_db: prefilter result database.
    - alignment_db: new alignment database.
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - min_aln_len: Minimum alignment length (default 0).
    - seq_id_mode: Identity is computed over 0 the alignment length, 1 the shorter or 2 the
      longer sequence (default 0).
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - add_backtrace: Store the backtrace string so that convertalis can write alignments.
    - alignment_mode: How the alignment is computed: 0 automatic, 1 score and end position
      only, 2 also start position and coverage, 3 also sequence identity, 4 ungapped only.
    - alignment_output_mode: What the alignment reports: 0 automatic, 1 score and end
      position only, 2 also start position and coverage, 3 also sequence identity, 4
      ungapped only.
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    result_db_path = core.input_database(result_db, "result_db")
    alignment_db_path = core.output_database(alignment_db, "alignment_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(result_db_path)
    positionals.append(alignment_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    options += core.opt_choice(
        "--alignment-output-mode", alignment_output_mode, (0, 1, 2, 3, 4), "alignment_output_mode"
    )
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_lcaalign",
        "lcaalign",
        positionals,
        options,
        databases={"alignment_db": alignment_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_majoritylca(
    target_db: str,
    result_db: str,
    taxonomy_db: str,
    *,
    majority: float | None = None,
    vote_mode: int | None = None,
    lca_ranks: str | None = None,
    blacklist: str | None = None,
    tax_lineage: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lowest common ancestor of the hits of every query by majority voting (`mmseqs
    majoritylca`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - target_db: taxonomy-annotated target sequence database.
    - result_db: alignment or prefilter result database.
    - taxonomy_db: new taxonomy result database.
    - majority: Minimal fraction of agreement among the taxonomically assigned sequences of
      a set (default 0.5).
    - vote_mode: Weights for the majority vote: 0 uniform, 1 minus log E-value, 2 score
      (default 1).
    - lca_ranks: Add a column with these ranks, comma separated, e.g. 'phylum,genus'.
    - blacklist: Comma-separated taxa ignored in the LCA computation (default 12908,28384).
    - tax_lineage: Lineage column: 0 none, 1 all lineage names, 2 all lineage taxon ids
      (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    target_db_path = core.input_database(target_db, "target_db")
    result_db_path = core.input_database(result_db, "result_db")
    taxonomy_db_path = core.output_database(taxonomy_db, "taxonomy_db")
    positionals: list[Path] = []
    positionals.append(target_db_path)
    positionals.append(result_db_path)
    positionals.append(taxonomy_db_path)
    options: list[str] = []
    options += core.opt_float("--majority", majority, "majority", 0.0, 1.0)
    options += core.opt_choice("--vote-mode", vote_mode, (0, 1, 2), "vote_mode")
    options += core.opt_text("--lca-ranks", lca_ranks, "lca_ranks", "^[A-Za-z_ ]+(,[A-Za-z_ ]+)*$")
    options += core.opt_text("--blacklist", blacklist, "blacklist", "^[0-9A-Za-z_:, .]+$")
    options += core.opt_choice("--tax-lineage", tax_lineage, (0, 1, 2), "tax_lineage")
    result = core.run(
        "mmseqs_majoritylca",
        "majoritylca",
        positionals,
        options,
        databases={"taxonomy_db": taxonomy_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_multihitdb(
    fasta_files: str | Sequence[str],
    set_db: str,
    *,
    dbtype: int | None = None,
    shuffle: bool | None = None,
    createdb_mode: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    max_gaps: int | None = None,
    contig_start_mode: int | None = None,
    contig_end_mode: int | None = None,
    orf_start_mode: int | None = None,
    forward_frames: str | None = None,
    reverse_frames: str | None = None,
    translation_table: int | None = None,
    use_all_table_starts: bool | None = None,
    id_offset: int | None = None,
    add_orf_stop: bool | None = None,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a grouped sequence database (sets of ORFs per contig) for multi-hit search (`mmseqs
    multihitdb`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - fasta_files: nucleotide FASTA file(s); every file becomes one set (for example a
      genome), its ORFs are the members.
    - set_db: prefix of the new set database.
    - dbtype: Database type: 0 automatic, 1 amino acid, 2 nucleotide (default 0).
    - shuffle: Shuffle the input order inside the database (default on).
    - createdb_mode: 0 copy the data, 1 soft-link the FASTA and write a new index (single-
      line FASTA/FASTQ only; default 0).
    - min_length: Minimum number of codons of an open reading frame (default 30).
    - max_length: Maximum number of codons of an open reading frame (default 32734).
    - max_gaps: Maximum number of codons with gaps or unknown residues before an ORF is
      rejected.
    - contig_start_mode: A contig start may be 0 incomplete, 1 complete or 2 both (default
      2).
    - contig_end_mode: A contig end may be 0 incomplete, 1 complete or 2 both (default 2).
    - orf_start_mode: ORF fragments run 0 from start to stop, 1 from any codon to stop, 2
      from the last start codon to stop (default 1).
    - forward_frames: Comma-separated forward-strand frames to extract, for example '1,2,3'
      (default 1,2,3).
    - reverse_frames: Comma-separated reverse-strand frames to extract, for example '1,3'
      (default 1,2,3).
    - translation_table: NCBI genetic code used for translation (default 1, standard).
    - use_all_table_starts: Accept every start codon of the genetic code, not only ATG.
    - id_offset: Offset of the numeric ids in the index file (default 0).
    - add_orf_stop: Add a stop symbol at complete ORF starts and ends.
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    fasta_files_path = core.fasta_files(fasta_files, "fasta_files")
    set_db_path = core.output_database(set_db, "set_db")
    positionals: list[Path] = []
    positionals.extend(fasta_files_path)
    positionals.append(set_db_path)
    options: list[str] = []
    options += core.opt_choice("--dbtype", dbtype, (0, 1, 2), "dbtype")
    options += core.opt_bool("--shuffle", shuffle, "shuffle")
    options += core.opt_choice("--createdb-mode", createdb_mode, (0, 1), "createdb_mode")
    options += core.opt_int("--min-length", min_length, "min_length", 1, None)
    options += core.opt_int("--max-length", max_length, "max_length", 1, None)
    options += core.opt_int("--max-gaps", max_gaps, "max_gaps", 0, None)
    options += core.opt_choice(
        "--contig-start-mode", contig_start_mode, (0, 1, 2), "contig_start_mode"
    )
    options += core.opt_choice("--contig-end-mode", contig_end_mode, (0, 1, 2), "contig_end_mode")
    options += core.opt_choice("--orf-start-mode", orf_start_mode, (0, 1, 2), "orf_start_mode")
    options += core.frames("--forward-frames", forward_frames, "forward_frames")
    options += core.frames("--reverse-frames", reverse_frames, "reverse_frames")
    options += core.opt_int("--translation-table", translation_table, "translation_table", 1, 31)
    options += core.opt_bool("--use-all-table-starts", use_all_table_starts, "use_all_table_starts")
    options += core.opt_int("--id-offset", id_offset, "id_offset", 0, None)
    options += core.opt_bool("--add-orf-stop", add_orf_stop, "add_orf_stop")
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    result = core.run(
        "mmseqs_multihitdb",
        "multihitdb",
        positionals,
        options,
        databases={"set_db": set_db_path},
        files={},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_multihitsearch(
    query_set_db: str,
    target_set_db: str,
    result_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    sensitivity: float | None = None,
    max_seqs: int | None = None,
    kmer_size: int | None = None,
    search_type: int | None = None,
    simple_best_hit: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a grouped query set database against a grouped target set database (`mmseqs
    multihitsearch`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_set_db: query set database from multihitdb.
    - target_set_db: target set database from multihitdb.
    - result_db: new result database.
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - min_aln_len: Minimum alignment length (default 0).
    - seq_id_mode: Identity is computed over 0 the alignment length, 1 the shorter or 2 the
      longer sequence (default 0).
    - sensitivity: Search sensitivity, 1.0 fast to 7.5 sensitive (default 4.0 for the
      prefilter, workflow-specific otherwise).
    - max_seqs: Maximum number of prefilter results per query (affects sensitivity; default
      300).
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - search_type: Search type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4
      translated nucleotide alignment (default 0).
    - simple_best_hit: Update the p-value with the single best hit (on) or with the best and
      second best hit (off); default on.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_set_db_path = core.input_database(query_set_db, "query_set_db")
    target_set_db_path = core.input_database(target_set_db, "target_set_db")
    result_db_path = core.output_database(result_db, "result_db")
    positionals: list[Path] = []
    positionals.append(query_set_db_path)
    positionals.append(target_set_db_path)
    positionals.append(result_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_bool("--simple-best-hit", simple_best_hit, "simple_best_hit")
    result = core.run(
        "mmseqs_multihitsearch",
        "multihitsearch",
        positionals,
        options,
        databases={"result_db": result_db_path},
        files={},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_besthitperset(
    query_set_db: str,
    target_set_db: str,
    result_db: str,
    best_hit_db: str,
    *,
    simple_best_hit: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """For every sequence set keep the best hit of the sequence-level result and update its
    p-value (`mmseqs besthitperset`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_set_db: query set database from mmseqs_multihitdb.
    - target_set_db: target set database from mmseqs_multihitdb.
    - result_db: sequence-level search or alignment result of the two set databases.
    - best_hit_db: new result database with one best hit per set pair.
    - simple_best_hit: Update the p-value with the single best hit (on) or with the best and
      second best hit (off); default on.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_set_db_path = core.input_database(query_set_db, "query_set_db")
    target_set_db_path = core.input_database(target_set_db, "target_set_db")
    result_db_path = core.input_database(result_db, "result_db")
    best_hit_db_path = core.output_database(best_hit_db, "best_hit_db")
    positionals: list[Path] = []
    positionals.append(query_set_db_path)
    positionals.append(target_set_db_path)
    positionals.append(result_db_path)
    positionals.append(best_hit_db_path)
    options: list[str] = []
    options += core.opt_bool("--simple-best-hit", simple_best_hit, "simple_best_hit")
    result = core.run(
        "mmseqs_besthitperset",
        "besthitperset",
        positionals,
        options,
        databases={"best_hit_db": best_hit_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_combinepvalperset(
    query_set_db: str,
    target_set_db: str,
    result_db: str,
    pvalue_db: str,
    *,
    alpha: float | None = None,
    aggregation_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the hit p-values of every pair of sequence sets (`mmseqs combinepvalperset`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_set_db: query set database from mmseqs_multihitdb.
    - target_set_db: target set database from mmseqs_multihitdb.
    - result_db: set-level result, for example the output of mmseqs_multihitsearch.
    - pvalue_db: new p-value database.
    - alpha: Alpha for combining p-values during aggregation (default 1).
    - aggregation_mode: Combined p-value from 0 multi-hit, 1 minimum of all p-values, 2
      product of p-values, 3 truncated product (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_set_db_path = core.input_database(query_set_db, "query_set_db")
    target_set_db_path = core.input_database(target_set_db, "target_set_db")
    result_db_path = core.input_database(result_db, "result_db")
    pvalue_db_path = core.output_database(pvalue_db, "pvalue_db")
    positionals: list[Path] = []
    positionals.append(query_set_db_path)
    positionals.append(target_set_db_path)
    positionals.append(result_db_path)
    positionals.append(pvalue_db_path)
    options: list[str] = []
    options += core.opt_float("--alpha", alpha, "alpha", 0.0, None)
    options += core.opt_choice(
        "--aggregation-mode", aggregation_mode, (0, 1, 2, 3), "aggregation_mode"
    )
    result = core.run(
        "mmseqs_combinepvalperset",
        "combinepvalperset",
        positionals,
        options,
        databases={"pvalue_db": pvalue_db_path},
        files={},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_mergeresultsbyset(
    set_to_member_db: str,
    result_db: str,
    merged_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the results of the ORFs of every set (contig) back onto the set (`mmseqs
    mergeresultsbyset`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - set_to_member_db: the <set_db>_set_to_member database that mmseqs_multihitdb wrote.
    - result_db: sequence-level result to merge, for example from mmseqs_besthitperset.
    - merged_db: new result database with one entry per set.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    set_to_member_db_path = core.input_database(set_to_member_db, "set_to_member_db")
    result_db_path = core.input_database(result_db, "result_db")
    merged_db_path = core.output_database(merged_db, "merged_db")
    positionals: list[Path] = []
    positionals.append(set_to_member_db_path)
    positionals.append(result_db_path)
    positionals.append(merged_db_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_mergeresultsbyset",
        "mergeresultsbyset",
        positionals,
        options,
        databases={"merged_db": merged_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


__all__ = [
    "mmseqs_addtaxonomy",
    "mmseqs_aggregate_taxonomy",
    "mmseqs_besthitperset",
    "mmseqs_combinepvalperset",
    "mmseqs_createtaxdb",
    "mmseqs_filtertaxdb",
    "mmseqs_filtertaxseqdb",
    "mmseqs_lcaalign",
    "mmseqs_majoritylca",
    "mmseqs_mergeresultsbyset",
    "mmseqs_multihitdb",
    "mmseqs_multihitsearch",
    "mmseqs_taxonomy",
    "mmseqs_taxonomyreport",
]

_register_public_verbose(globals(), tuple(__all__))
