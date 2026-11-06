"""MMseqs2 search, profile and database-utility facades (integration pack P30).

Command-line facades over MMseqs2, complementing the P07 facades. Each public function runs one
fixed MMseqs2 module with an explicit argument list (no shell) through
``dnakit.backends.external.execute_bounded_command`` and the private runner of the P07 facades
(``_mmseqs2_p07_core``). The user installs MMseqs2 (conda or a release binary) and puts ``mmseqs``
on PATH or passes ``mmseqs_path``; the tested version is MMseqs2 release 18-8cc5c. MMseqs2
databases are addressed by path prefix; options that are not given are not passed, so MMseqs2's
own defaults apply. Results are JSON-serialisable dictionaries with status, version, command,
outputs, entry or line counts and a bounded log tail. Outputs are never overwritten.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment import _mmseqs2_p07_core as core
from dnakit.bioinformatics import _gap_cli

_DATABASES = (
    "UniRef100",
    "UniRef90",
    "UniRef50",
    "UniProtKB",
    "UniProtKB/TrEMBL",
    "UniProtKB/Swiss-Prot",
    "NR",
    "NT",
    "GTDB",
    "PDB",
    "PDB70",
    "Pfam-A.full",
    "Pfam-A.seed",
    "Pfam-B",
    "CDD",
    "eggNOG",
    "VOGDB",
    "dbCAN2",
    "SILVA",
    "Resfinder",
    "Kalamari",
)


def _alignment_options(
    *,
    e_value: float | None,
    min_seq_id: float | None,
    coverage: float | None,
    cov_mode: int | None,
    min_aln_len: int | None,
    add_backtrace: bool | None,
    realign: bool | None,
    gap_open: int | None,
    gap_extend: int | None,
    max_accept: int | None,
    max_rejected: int | None,
    zdrop: int | None,
) -> list[str]:
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_bool("--realign", realign, "realign")
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    options += core.opt_int("--max-accept", max_accept, "max_accept", 1, None)
    options += core.opt_int("--max-rejected", max_rejected, "max_rejected", 0, None)
    options += core.opt_int("--zdrop", zdrop, "zdrop", 0, None)
    return options


def mmseqs_search(
    query_db: str,
    target_db: str,
    alignment_db: str,
    *,
    sensitivity: float | None = None,
    kmer_length: int | None = None,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    max_seqs: int | None = None,
    search_type: int | None = None,
    num_iterations: int | None = None,
    start_sens: float | None = None,
    sens_steps: int | None = None,
    mask: bool | None = None,
    strand: int | None = None,
    add_backtrace: bool | None = None,
    realign: bool | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    max_accept: int | None = None,
    max_rejected: int | None = None,
    zdrop: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a query database against a target database (`mmseqs search`).

    Databases are given as path prefixes (create them with ``mmseqs_createdb``). The result is an
    alignment database that ``mmseqs_convertalis`` turns into a table. An option left at None is
    not passed, so the MMseqs2 default applies.

    Parameters:
    - query_db, target_db: prefixes of the query and target sequence (or profile) databases.
    - alignment_db: prefix of the new alignment database.
    - sensitivity: -s, 1.0 faster, 4.0 fast, 7.5 sensitive (default 5.7).
    - kmer_length: -k, k-mer length (0 chooses automatically).
    - e_value: -e, list matches below this E-value (default 0.001).
    - min_seq_id: --min-seq-id, minimum sequence identity 0.0-1.0.
    - coverage, cov_mode: -c and --cov-mode (0 query and target, 1 target, 2 query, 3-5 length
      criteria).
    - min_aln_len: minimum alignment length.
    - max_seqs: --max-seqs, results per query passing the prefilter (default 300).
    - search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
      nucleotide alignment.
    - num_iterations: iterative profile search rounds (default 1); start_sens and sens_steps set
      the sensitivity ladder.
    - mask: low-complexity masking with tantan in the prefilter (default on).
    - strand: 0 reverse, 1 forward, 2 both (nucleotide searches only; default 1).
    - add_backtrace: -a, store the alignment backtrace so that convertalis can write CIGARs.
    - realign: shorter, more conservative alignments.
    - gap_open, gap_extend, max_accept, max_rejected, zdrop: alignment scoring and stopping
      parameters (zdrop applies to nucleotide alignment).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_path = core.input_database(query_db, "query_db")
    target_path = core.input_database(target_db, "target_db")
    output_path = core.output_database(alignment_db, "alignment_db")
    options: list[str] = []
    options += core.opt_float("-s", sensitivity, "sensitivity", 1.0, 8.0)
    options += core.opt_int("-k", kmer_length, "kmer_length", 0, 32)
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_int("--num-iterations", num_iterations, "num_iterations", 1, 20)
    options += core.opt_float("--start-sens", start_sens, "start_sens", 1.0, 8.0)
    options += core.opt_int("--sens-steps", sens_steps, "sens_steps", 1, 20)
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_choice("--strand", strand, (0, 1, 2), "strand")
    options += _alignment_options(
        e_value=e_value,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
        min_aln_len=min_aln_len,
        add_backtrace=add_backtrace,
        realign=realign,
        gap_open=gap_open,
        gap_extend=gap_extend,
        max_accept=max_accept,
        max_rejected=max_rejected,
        zdrop=zdrop,
    )
    return core.run(
        "mmseqs_search",
        "search",
        [query_path, target_path, output_path],
        options,
        databases={"alignment_db": output_path},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_align(
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
    alignment_mode: int | None = None,
    alignment_output_mode: int | None = None,
    add_backtrace: bool | None = None,
    realign: bool | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    max_accept: int | None = None,
    max_rejected: int | None = None,
    zdrop: int | None = None,
    wrapped_scoring: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Re-align the hits of a prefilter or search result with Smith-Waterman (`mmseqs align`).

    Databases are given as path prefixes. ``result_db`` holds the candidate hits (for example
    the output of ``mmseqs_prefilter``); the new ``alignment_db`` holds scored alignments that
    ``mmseqs_convertalis`` can export. An option left at None is not passed.

    Parameters:
    - query_db, target_db, result_db: existing databases; alignment_db: new alignment database.
    - e_value, min_seq_id, coverage, cov_mode, min_aln_len: hit filters as in ``mmseqs_search``.
    - alignment_mode: 0 automatic, 1 score and end position, 2 also start and coverage, 3 also
      sequence identity.
    - alignment_output_mode: 0 automatic, 1-3 as above, 4 ungapped only, 5 score only.
    - add_backtrace, realign, gap_open, gap_extend, max_accept, max_rejected, zdrop: as in
      ``mmseqs_search``.
    - wrapped_scoring: double the nucleotide query to allow diagonals that wrap around its ends
      (circular sequences).
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    query_path = core.input_database(query_db, "query_db")
    target_path = core.input_database(target_db, "target_db")
    result_path = core.input_database(result_db, "result_db")
    output_path = core.output_database(alignment_db, "alignment_db")
    options: list[str] = []
    options += core.opt_choice("--alignment-mode", alignment_mode, (0, 1, 2, 3), "alignment_mode")
    options += core.opt_choice(
        "--alignment-output-mode",
        alignment_output_mode,
        (0, 1, 2, 3, 4, 5),
        "alignment_output_mode",
    )
    options += core.opt_bool("--wrapped-scoring", wrapped_scoring, "wrapped_scoring")
    options += _alignment_options(
        e_value=e_value,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
        min_aln_len=min_aln_len,
        add_backtrace=add_backtrace,
        realign=realign,
        gap_open=gap_open,
        gap_extend=gap_extend,
        max_accept=max_accept,
        max_rejected=max_rejected,
        zdrop=zdrop,
    )
    return core.run(
        "mmseqs_align",
        "align",
        [query_path, target_path, result_path, output_path],
        options,
        databases={"alignment_db": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_createindex(
    sequence_db: str,
    *,
    kmer_length: int | None = None,
    sensitivity: float | None = None,
    mask: bool | None = None,
    search_type: int | None = None,
    split: int | None = None,
    spaced_kmer_mode: bool | None = None,
    spaced_kmer_pattern: str | None = None,
    max_seqs: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Precompute the k-mer index of a sequence database (`mmseqs createindex`).

    The index is written next to the database as ``<sequence_db>.idx*`` and speeds up later
    searches against it. The call refuses a database that already has an index.

    Parameters:
    - sequence_db: prefix of an existing sequence database.
    - kmer_length: -k, k-mer length (0 chooses automatically).
    - sensitivity: -s, sensitivity the index is built for (default 7.5).
    - mask: mask low-complexity regions with tantan (default on).
    - search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
      nucleotide alignment.
    - split: number of index chunks (0 chooses automatically).
    - spaced_kmer_mode, spaced_kmer_pattern: spaced k-mers (default on) and a custom 0/1 pattern.
    - max_seqs: maximum prefilter results per query that the index is sized for.
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    database = core.input_database(sequence_db, "sequence_db")
    index_marker = Path(f"{database}.idx.dbtype")
    index_files = Path(f"{database}.idx")
    core.reject_with(
        index_marker.exists() or index_files.exists(),
        f"{database} already has an index; remove the .idx files first.",
    )
    options: list[str] = []
    options += core.opt_int("-k", kmer_length, "kmer_length", 0, 32)
    options += core.opt_float("-s", sensitivity, "sensitivity", 1.0, 8.0)
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_int("--split", split, "split", 0, None)
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_text(
        "--spaced-kmer-pattern", spaced_kmer_pattern, "spaced_kmer_pattern", r"[01]+"
    )
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    result = core.run(
        "mmseqs_createindex",
        "createindex",
        [database],
        options,
        databases={},
        extra_outputs={"index": Path(f"{database}.idx")},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        anchor=database.parent,
        verbose=verbose,
    )
    result["outputs"]["index_dbtype"] = str(index_marker)
    return result


def mmseqs_createsubdb(
    subset: str,
    input_db: str,
    output_db: str,
    *,
    subdb_mode: int | None = None,
    id_mode: int | None = None,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract a subset of the entries of a database into a new database (`mmseqs createsubdb`).

    Parameters:
    - subset: text file with one database key (or FASTA identifier, see ``id_mode``) per line,
      or the prefix of another database whose keys are used.
    - input_db: prefix of the existing database.
    - output_db: prefix of the new database.
    - subdb_mode: 0 copy the data, 1 soft-link the data and write a new index (default 0).
    - id_mode: 0 select by database key, 1 select by FASTA identifier through the .lookup file
      (default 0).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    subset_path = Path(subset).expanduser().resolve()
    if not subset_path.is_file():
        subset_path = core.input_database(subset, "subset")
    database = core.input_database(input_db, "input_db")
    output_path = core.output_database(output_db, "output_db")
    options: list[str] = []
    options += core.opt_choice("--subdb-mode", subdb_mode, (0, 1), "subdb_mode")
    options += core.opt_choice("--id-mode", id_mode, (0, 1), "id_mode")
    return core.run(
        "mmseqs_createsubdb",
        "createsubdb",
        [subset_path, database, output_path],
        options,
        databases={"output_db": output_path},
        needs_tmp=False,
        threads=None,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_databases(
    name: str,
    sequence_db: str,
    *,
    remove_tmp_files: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and prepare a named reference database (`mmseqs databases`).

    This contacts the database provider over the network and can download many gigabytes; it
    is never run by the ordinary tests.

    Parameters:
    - name: one of UniRef100, UniRef90, UniRef50, UniProtKB, UniProtKB/TrEMBL,
      UniProtKB/Swiss-Prot, NR, NT, GTDB, PDB, PDB70, Pfam-A.full, Pfam-A.seed, Pfam-B, CDD,
      eggNOG, VOGDB, dbCAN2, SILVA, Resfinder or Kalamari (NT, SILVA, Resfinder and Kalamari
      hold nucleotide sequences).
    - sequence_db: prefix of the new database.
    - remove_tmp_files: delete the downloaded temporary files afterwards.
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    chosen = _gap_cli.choice(name, _DATABASES, "name")
    output_path = core.output_database(sequence_db, "sequence_db")
    options = core.opt_bool("--remove-tmp-files", remove_tmp_files, "remove_tmp_files")
    result = core.run(
        "mmseqs_databases",
        "databases",
        [chosen, output_path],
        options,
        databases={"sequence_db": output_path},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["database_name"] = chosen
    return result


def mmseqs_convertmsa(
    msa_file: str,
    msa_db: str,
    *,
    identifier_field: int | None = None,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Stockholm multiple alignments into an MMseqs2 MSA database (`mmseqs convertmsa`).

    Parameters:
    - msa_file: Stockholm file (``.sto``, optionally gzip-compressed) with one or many
      alignments.
    - msa_db: prefix of the new MSA database.
    - identifier_field: Stockholm comment used as the MSA identifier: 0 ID, 1 AC (default 0;
      when the comment is missing the name of the first sequence is used).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    source = core.input_file(msa_file, "msa_file")
    output_path = core.output_database(msa_db, "msa_db")
    options = core.opt_choice("--identifier-field", identifier_field, (0, 1), "identifier_field")
    return core.run(
        "mmseqs_convertmsa",
        "convertmsa",
        [source, output_path],
        options,
        databases={"msa_db": output_path},
        needs_tmp=False,
        threads=None,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_msa2profile(
    msa_db: str,
    profile_db: str,
    *,
    msa_type: int | None = None,
    match_mode: int | None = None,
    match_ratio: float | None = None,
    skip_query: bool | None = None,
    filter_msa: bool | None = None,
    max_seq_id: float | None = None,
    coverage: float | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Turn an MSA database into a profile database (`mmseqs msa2profile`).

    Parameters:
    - msa_db: prefix of an MSA database (from ``mmseqs_convertmsa`` or ``mmseqs_result2msa``).
    - profile_db: prefix of the new profile database.
    - msa_type: 0 cA3M, 1 A3M, 2 FASTA (default 2).
    - match_mode: 0 keep columns with a residue in the first sequence, 1 keep columns with a
      residue in ``match_ratio`` of all sequences (default 0).
    - match_ratio: column occupancy for match_mode 1 (default 0.5).
    - skip_query: leave the first (query) sequence out of the profile.
    - filter_msa: filter the MSA before building the profile (default on).
    - max_seq_id: maximum pairwise identity kept by the filter (default 0.9).
    - coverage: minimum fraction of the query covered by member sequences (default 0.0).
    - gap_open, gap_extend: gap costs used for the profile.
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    source = core.input_database(msa_db, "msa_db")
    output_path = core.output_database(profile_db, "profile_db")
    options: list[str] = []
    options += core.opt_choice("--msa-type", msa_type, (0, 1, 2), "msa_type")
    options += core.opt_choice("--match-mode", match_mode, (0, 1), "match_mode")
    options += core.opt_float("--match-ratio", match_ratio, "match_ratio", 0.0, 1.0)
    options += core.opt_bool("--skip-query", skip_query, "skip_query")
    options += core.opt_bool("--filter-msa", filter_msa, "filter_msa")
    options += core.opt_float("--max-seq-id", max_seq_id, "max_seq_id", 0.0, 1.0)
    options += core.opt_float("--cov", coverage, "coverage", 0.0, 1.0)
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    return core.run(
        "mmseqs_msa2profile",
        "msa2profile",
        [source, output_path],
        options,
        databases={"profile_db": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_result2profile(
    query_db: str,
    target_db: str,
    result_db: str,
    profile_db: str,
    *,
    e_profile: float | None = None,
    mask_profile: bool | None = None,
    filter_msa: bool | None = None,
    max_seq_id: float | None = None,
    coverage: float | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a profile database from search results (`mmseqs result2profile`).

    For each query the hits of ``result_db`` (an alignment database from ``mmseqs_search``
    run with backtraces, ``add_backtrace=True``) are stacked into an MSA and converted to a
    profile.

    Parameters:
    - query_db, target_db, result_db: existing databases; profile_db: new profile database.
    - e_profile: include hits with an E-value below this threshold (default 0.001).
    - mask_profile: mask the query with tantan before building the profile (default on).
    - filter_msa, max_seq_id, coverage: MSA filtering as in ``mmseqs_msa2profile``.
    - gap_open, gap_extend: gap costs used for the profile.
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    query_path = core.input_database(query_db, "query_db")
    target_path = core.input_database(target_db, "target_db")
    result_path = core.input_database(result_db, "result_db")
    output_path = core.output_database(profile_db, "profile_db")
    options: list[str] = []
    options += core.opt_float("--e-profile", e_profile, "e_profile", 0.0, None)
    options += core.opt_bool("--mask-profile", mask_profile, "mask_profile")
    options += core.opt_bool("--filter-msa", filter_msa, "filter_msa")
    options += core.opt_float("--max-seq-id", max_seq_id, "max_seq_id", 0.0, 1.0)
    options += core.opt_float("--cov", coverage, "coverage", 0.0, 1.0)
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    return core.run(
        "mmseqs_result2profile",
        "result2profile",
        [query_path, target_path, result_path, output_path],
        options,
        databases={"profile_db": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_profile2consensus(
    profile_db: str,
    sequence_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the consensus sequence of every profile into a sequence database
    (`mmseqs profile2consensus`).

    Parameters:
    - profile_db: prefix of a profile database.
    - sequence_db: prefix of the new sequence database (one consensus per profile).
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    source = core.input_database(profile_db, "profile_db")
    output_path = core.output_database(sequence_db, "sequence_db")
    return core.run(
        "mmseqs_profile2consensus",
        "profile2consensus",
        [source, output_path],
        [],
        databases={"sequence_db": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_profile2repseq(
    profile_db: str,
    sequence_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the representative member sequence of every profile (`mmseqs profile2repseq`).

    The representative is the profile's own query sequence as stored in the profile database.

    Parameters:
    - profile_db: prefix of a profile database.
    - sequence_db: prefix of the new sequence database.
    - threads, mmseqs_path: as in ``mmseqs_search``.
    """
    source = core.input_database(profile_db, "profile_db")
    output_path = core.output_database(sequence_db, "sequence_db")
    return core.run(
        "mmseqs_profile2repseq",
        "profile2repseq",
        [source, output_path],
        [],
        databases={"sequence_db": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_profile2pssm(
    profile_db: str,
    pssm_file: str,
    *,
    db_output: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export profiles as position-specific scoring matrices (`mmseqs profile2pssm`).

    Parameters:
    - profile_db: prefix of a profile database.
    - pssm_file: new text file with one PSSM per profile (a database prefix when
      ``db_output`` is true).
    - db_output: return a result database instead of a text file (default off).
    - threads, mmseqs_path: as in ``mmseqs_search``.

    For the text output the result adds ``pssm_preview`` (the first 20 lines).
    """
    source = core.input_database(profile_db, "profile_db")
    options = core.opt_bool("--db-output", db_output, "db_output")
    if db_output:
        output_db = core.output_database(pssm_file, "pssm_file")
        return core.run(
            "mmseqs_profile2pssm",
            "profile2pssm",
            [source, output_db],
            options,
            databases={"pssm_db": output_db},
            needs_tmp=False,
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    output_path = core.output_file(pssm_file, "pssm_file")
    result = core.run(
        "mmseqs_profile2pssm",
        "profile2pssm",
        [source, output_path],
        options,
        files={"pssm": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["pssm_preview"] = [row for row in core.table(output_path)[:20]]
    return result


def mmseqs_unpackdb(
    database: str,
    output_dir: str,
    *,
    unpack_name_mode: int | None = None,
    unpack_suffix: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unpack every entry of a database into its own file (`mmseqs unpackdb`).

    Typical use is an MSA database, which is unpacked into one alignment file per entry.

    Parameters:
    - database: prefix of an existing database.
    - output_dir: new folder that receives the files.
    - unpack_name_mode: 0 name files by database key, 1 by accession through the .lookup file
      (default 1).
    - unpack_suffix: file suffix for the unpacked files, for example ``.a3m``; add ``.gz`` to
      write compressed files.
    - threads, mmseqs_path: as in ``mmseqs_search``.

    The result adds ``files`` (the unpacked file names).
    """
    source = core.input_database(database, "database")
    target = Path(output_dir).expanduser().resolve()
    core.reject_with(target.exists(), f"output_dir must not exist yet: {target}")
    core.reject_with(
        not target.parent.is_dir(), f"The folder for output_dir does not exist: {target.parent}"
    )
    options: list[str] = []
    options += core.opt_choice("--unpack-name-mode", unpack_name_mode, (0, 1), "unpack_name_mode")
    options += core.opt_text("--unpack-suffix", unpack_suffix, "unpack_suffix", r"[A-Za-z0-9._-]+")
    target.mkdir()
    try:
        result = core.run(
            "mmseqs_unpackdb",
            "unpackdb",
            [source, target],
            options,
            databases={},
            needs_tmp=False,
            threads=threads,
            mmseqs_path=mmseqs_path,
            anchor=target,
            verbose=verbose,
        )
    except Exception:
        import shutil

        shutil.rmtree(target, ignore_errors=True)
        raise
    names = sorted(path.name for path in target.iterdir())
    result["outputs"]["output_dir"] = str(target)
    result["files"] = names
    result["file_count"] = len(names)
    return result


__all__ = [
    "mmseqs_align",
    "mmseqs_convertmsa",
    "mmseqs_createindex",
    "mmseqs_createsubdb",
    "mmseqs_databases",
    "mmseqs_msa2profile",
    "mmseqs_profile2consensus",
    "mmseqs_profile2pssm",
    "mmseqs_profile2repseq",
    "mmseqs_result2profile",
    "mmseqs_search",
    "mmseqs_unpackdb",
]

_register_public_verbose(globals(), tuple(__all__))
