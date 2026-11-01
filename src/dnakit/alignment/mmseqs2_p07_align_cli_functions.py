"""MMseqs2 prefilters, alignment modules and ORF utilities.

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

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment import _mmseqs2_p07_core as core


def mmseqs_prefilter(
    query_db: str,
    target_db: str,
    prefilter_db: str,
    *,
    sensitivity: float | None = None,
    kmer_size: int | None = None,
    max_seqs: int | None = None,
    min_ungapped_score: int | None = None,
    comp_bias_corr: bool | None = None,
    add_self_matches: bool | None = None,
    spaced_kmer_mode: bool | None = None,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    taxon_list: str | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """k-mer prefilter with double consecutive diagonal hits (`mmseqs prefilter`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - prefilter_db: new prefilter result database.
    - sensitivity: Search sensitivity, 1.0 fast to 7.5 sensitive (default 4.0 for the
      prefilter, workflow-specific otherwise).
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - max_seqs: Maximum number of prefilter results per query (affects sensitivity; default
      300).
    - min_ungapped_score: Accept only matches whose ungapped alignment score is at least
      this (default 15).
    - comp_bias_corr: Correct for locally biased composition (default on).
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - spaced_kmer_mode: Use spaced k-mers instead of consecutive k-mers.
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - taxon_list: Taxon ids to keep, e.g. '9606||810' (or), '2759&&!4751' (and not); only
      digits, ',', '&', '|', '!' and brackets.
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    prefilter_db_path = core.output_database(prefilter_db, "prefilter_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(prefilter_db_path)
    options: list[str] = []
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int(
        "--min-ungapped-score", min_ungapped_score, "min_ungapped_score", 0, None
    )
    options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_prefilter",
        "prefilter",
        positionals,
        options,
        databases={"prefilter_db": prefilter_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_ungappedprefilter(
    query_db: str,
    target_db: str,
    prefilter_db: str,
    *,
    max_seqs: int | None = None,
    min_ungapped_score: int | None = None,
    comp_bias_corr: bool | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    e_value: float | None = None,
    taxon_list: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prefilter that scores the optimal ungapped diagonal of every query-target pair (`mmseqs
    ungappedprefilter`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - prefilter_db: new prefilter result database.
    - max_seqs: Maximum number of prefilter results per query (affects sensitivity; default
      300).
    - min_ungapped_score: Accept only matches whose ungapped alignment score is at least
      this (default 15).
    - comp_bias_corr: Correct for locally biased composition (default on).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - taxon_list: Taxon ids to keep, e.g. '9606||810' (or), '2759&&!4751' (and not); only
      digits, ',', '&', '|', '!' and brackets.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    prefilter_db_path = core.output_database(prefilter_db, "prefilter_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(prefilter_db_path)
    options: list[str] = []
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int(
        "--min-ungapped-score", min_ungapped_score, "min_ungapped_score", 0, None
    )
    options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    result = core.run(
        "mmseqs_ungappedprefilter",
        "ungappedprefilter",
        positionals,
        options,
        databases={"prefilter_db": prefilter_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_gappedprefilter(
    query_db: str,
    target_db: str,
    prefilter_db: str,
    *,
    max_seqs: int | None = None,
    min_ungapped_score: int | None = None,
    comp_bias_corr: bool | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    e_value: float | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    taxon_list: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prefilter that scores the optimal gapped Smith-Waterman alignment of every query-target
    pair, slow (`mmseqs gappedprefilter`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - prefilter_db: new prefilter result database.
    - max_seqs: Maximum number of prefilter results per query (affects sensitivity; default
      300).
    - min_ungapped_score: Accept only matches whose ungapped alignment score is at least
      this (default 15).
    - comp_bias_corr: Correct for locally biased composition (default on).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - taxon_list: Taxon ids to keep, e.g. '9606||810' (or), '2759&&!4751' (and not); only
      digits, ',', '&', '|', '!' and brackets.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    prefilter_db_path = core.output_database(prefilter_db, "prefilter_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(prefilter_db_path)
    options: list[str] = []
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int(
        "--min-ungapped-score", min_ungapped_score, "min_ungapped_score", 0, None
    )
    options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    options += core.taxon_list("--taxon-list", taxon_list, "taxon_list")
    result = core.run(
        "mmseqs_gappedprefilter",
        "gappedprefilter",
        positionals,
        options,
        databases={"prefilter_db": prefilter_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_kmermatcher(
    sequence_db: str,
    prefilter_db: str,
    *,
    kmer_size: int | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    kmer_per_seq: int | None = None,
    spaced_kmer_mode: bool | None = None,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    adjust_kmer_len: bool | None = None,
    hash_shift: int | None = None,
    ignore_multi_kmer: bool | None = None,
    linclust_version: int | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find bottom-m-hashed k-mer matches within one sequence database (`mmseqs kmermatcher`);
    the building block of linear-time clustering.

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - prefilter_db: new prefilter result database.
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - kmer_per_seq: Number of k-mers selected per sequence (linear-time workflows).
    - spaced_kmer_mode: Use spaced k-mers instead of consecutive k-mers.
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - adjust_kmer_len: Adjust the k-mer length by specificity (nucleotides only).
    - hash_shift: Shift of the k-mer hash initialisation (default 67).
    - ignore_multi_kmer: Skip k-mers that occur more than once in a sequence.
    - linclust_version: Linclust algorithm version 1 or 2 (default 2).
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    prefilter_db_path = core.output_database(prefilter_db, "prefilter_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(prefilter_db_path)
    options: list[str] = []
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--kmer-per-seq", kmer_per_seq, "kmer_per_seq", 0, None)
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    options += core.opt_bool("--adjust-kmer-len", adjust_kmer_len, "adjust_kmer_len")
    options += core.opt_int("--hash-shift", hash_shift, "hash_shift", 0, None)
    options += core.opt_bool("--ignore-multi-kmer", ignore_multi_kmer, "ignore_multi_kmer")
    options += core.opt_choice("--linclust-version", linclust_version, (1, 2), "linclust_version")
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_kmermatcher",
        "kmermatcher",
        positionals,
        options,
        databases={"prefilter_db": prefilter_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_kmersearch(
    query_db: str,
    kmer_index_db: str,
    prefilter_db: str,
    *,
    coverage: float | None = None,
    cov_mode: int | None = None,
    kmer_per_seq: int | None = None,
    pick_n_sim_kmer: int | None = None,
    result_direction: int | None = None,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find bottom-m-hashed k-mer matches between a query database and a k-mer index of a target
    database (`mmseqs kmersearch`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - kmer_index_db: k-mer index of the target from mmseqs_createlinindex (the .linidx
      database).
    - prefilter_db: new prefilter result database.
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - kmer_per_seq: Number of k-mers selected per sequence (linear-time workflows).
    - pick_n_sim_kmer: Number of similar k-mers added to the search (default 1).
    - result_direction: Result is 0 query-centric or 1 target-centric (default 1).
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    kmer_index_db_path = core.input_database(kmer_index_db, "kmer_index_db")
    prefilter_db_path = core.output_database(prefilter_db, "prefilter_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(kmer_index_db_path)
    positionals.append(prefilter_db_path)
    options: list[str] = []
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--kmer-per-seq", kmer_per_seq, "kmer_per_seq", 0, None)
    options += core.opt_int("--pick-n-sim-kmer", pick_n_sim_kmer, "pick_n_sim_kmer", 0, None)
    options += core.opt_choice("--result-direction", result_direction, (0, 1), "result_direction")
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_kmersearch",
        "kmersearch",
        positionals,
        options,
        databases={"prefilter_db": prefilter_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_alignall(
    sequence_db: str,
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
    add_self_matches: bool | None = None,
    add_backtrace: bool | None = None,
    alignment_mode: int | None = None,
    comp_bias_corr: bool | None = None,
    score_bias: float | None = None,
    zdrop: int | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gapped local all-versus-all alignment of the members inside every result entry (`mmseqs
    alignall`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - result_db: result (for example cluster) database.
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
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - add_backtrace: Store the backtrace string so that convertalis can write alignments.
    - alignment_mode: How the alignment is computed: 0 automatic, 1 score and end position
      only, 2 also start position and coverage, 3 also sequence identity, 4 ungapped only.
    - comp_bias_corr: Correct for locally biased composition (default on).
    - score_bias: Score bias added when computing the Smith-Waterman alignment, in bits
      (default 0).
    - zdrop: Largest allowed score drop before a nucleotide alignment is truncated (default
      40).
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    result_db_path = core.input_database(result_db, "result_db")
    alignment_db_path = core.output_database(alignment_db, "alignment_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
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
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_float("--score-bias", score_bias, "score_bias", None, None)
    options += core.opt_int("--zdrop", zdrop, "zdrop", 0, None)
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_alignall",
        "alignall",
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


def mmseqs_transitivealign(
    sequence_db: str,
    alignment_db: str,
    transitive_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    add_self_matches: bool | None = None,
    add_backtrace: bool | None = None,
    alignment_mode: int | None = None,
    alignment_output_mode: int | None = None,
    alt_ali: int | None = None,
    max_rejected: int | None = None,
    max_accept: int | None = None,
    score_bias: float | None = None,
    realign: bool | None = None,
    zdrop: int | None = None,
    wrapped_scoring: bool | None = None,
    comp_bias_corr: bool | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transfer alignments through a shared sequence by transitivity (`mmseqs transitivealign`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - alignment_db: alignment database to extend.
    - transitive_db: new alignment database.
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
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - add_backtrace: Store the backtrace string so that convertalis can write alignments.
    - alignment_mode: How the alignment is computed: 0 automatic, 1 score and end position
      only, 2 also start position and coverage, 3 also sequence identity, 4 ungapped only.
    - alignment_output_mode: What the alignment reports: 0 automatic, 1 score and end
      position only, 2 also start position and coverage, 3 also sequence identity, 4
      ungapped only.
    - alt_ali: Show up to this many alternative alignments (default 0).
    - max_rejected: Stop aligning a query after this many rejected alignments.
    - max_accept: Stop aligning a query after this many accepted alignments.
    - score_bias: Score bias added when computing the Smith-Waterman alignment, in bits
      (default 0).
    - realign: Compute more conservative, shorter alignments (scores and E-values
      unchanged).
    - zdrop: Largest allowed score drop before a nucleotide alignment is truncated (default
      40).
    - wrapped_scoring: Double the nucleotide query during scoring so that diagonals may wrap
      around the ends (circular sequences).
    - comp_bias_corr: Correct for locally biased composition (default on).
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    transitive_db_path = core.output_database(transitive_db, "transitive_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(alignment_db_path)
    positionals.append(transitive_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    options += core.opt_choice(
        "--alignment-output-mode", alignment_output_mode, (0, 1, 2, 3, 4), "alignment_output_mode"
    )
    options += core.opt_int("--alt-ali", alt_ali, "alt_ali", 0, None)
    options += core.opt_int("--max-rejected", max_rejected, "max_rejected", 0, None)
    options += core.opt_int("--max-accept", max_accept, "max_accept", 0, None)
    options += core.opt_float("--score-bias", score_bias, "score_bias", None, None)
    options += core.opt_bool("--realign", realign, "realign")
    options += core.opt_int("--zdrop", zdrop, "zdrop", 0, None)
    options += core.opt_bool("--wrapped-scoring", wrapped_scoring, "wrapped_scoring")
    options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_transitivealign",
        "transitivealign",
        positionals,
        options,
        databases={"transitive_db": transitive_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_rescorediagonal(
    query_db: str,
    target_db: str,
    prefilter_db: str,
    result_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    rescore_mode: int | None = None,
    add_self_matches: bool | None = None,
    add_backtrace: bool | None = None,
    filter_hits: bool | None = None,
    sort_results: int | None = None,
    wrapped_scoring: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the sequence identity of every prefilter diagonal by Hamming distance or local or
    global alignment (`mmseqs rescorediagonal`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - prefilter_db: prefilter result database.
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
    - rescore_mode: Rescore diagonals with 0 Hamming distance, 1 local alignment (score
      only), 2 local alignment, 3 global alignment, 4 longest alignment fulfilling the
      window criterion (default 0).
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - add_backtrace: Store the backtrace string so that convertalis can write alignments.
    - filter_hits: Filter hits by sequence identity and coverage.
    - sort_results: 0 do not sort, 1 sort by E-value (alignments) or identity (Hamming)
      (default 0).
    - wrapped_scoring: Double the nucleotide query during scoring so that diagonals may wrap
      around the ends (circular sequences).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    prefilter_db_path = core.input_database(prefilter_db, "prefilter_db")
    result_db_path = core.output_database(result_db, "result_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(prefilter_db_path)
    positionals.append(result_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_choice("--rescore-mode", rescore_mode, (0, 1, 2, 3, 4), "rescore_mode")
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_bool("--filter-hits", filter_hits, "filter_hits")
    options += core.opt_choice("--sort-results", sort_results, (0, 1), "sort_results")
    options += core.opt_bool("--wrapped-scoring", wrapped_scoring, "wrapped_scoring")
    result = core.run(
        "mmseqs_rescorediagonal",
        "rescorediagonal",
        positionals,
        options,
        databases={"result_db": result_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_fwbw(
    query_db: str,
    target_db: str,
    alignment_db: str,
    result_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    mact: float | None = None,
    fwbw_gapopen: float | None = None,
    fwbw_gapextend: float | None = None,
    temperature: float | None = None,
    blocklen: int | None = None,
    fwbw_backtrace_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Forward-backward (maximum accuracy) alignment of the pairs of an alignment database
    (`mmseqs fwbw`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - alignment_db: alignment database to realign.
    - result_db: new alignment database.
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - min_aln_len: Minimum alignment length (default 0).
    - seq_id_mode: Identity is computed over 0 the alignment length, 1 the shorter or 2 the
      longer sequence (default 0).
    - mact: Maximum-accuracy threshold of forward-backward alignment (default 0.035).
    - fwbw_gapopen: Gap open penalty of forward-backward alignment (default 10).
    - fwbw_gapextend: Gap extension penalty of forward-backward alignment (default 2).
    - temperature: Temperature of the forward-backward algorithm (default 1).
    - blocklen: Block length of the forward-backward algorithm (default 16).
    - fwbw_backtrace_mode: Backtrace: 0 none, 1 local (default 1).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    result_db_path = core.output_database(result_db, "result_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(alignment_db_path)
    positionals.append(result_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_float("--mact", mact, "mact", 0.0, None)
    options += core.opt_float("--fwbw-gapopen", fwbw_gapopen, "fwbw_gapopen", 0.0, None)
    options += core.opt_float("--fwbw-gapextend", fwbw_gapextend, "fwbw_gapextend", 0.0, None)
    options += core.opt_float("--temperature", temperature, "temperature", 0.0, None)
    options += core.opt_int("--blocklen", blocklen, "blocklen", 1, None)
    options += core.opt_choice(
        "--fwbw-backtrace-mode", fwbw_backtrace_mode, (0, 1), "fwbw_backtrace_mode"
    )
    result = core.run(
        "mmseqs_fwbw",
        "fwbw",
        positionals,
        options,
        databases={"result_db": result_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_alignbykmer(
    query_db: str,
    target_db: str,
    result_db: str,
    alignment_db: str,
    *,
    kmer_size: int | None = None,
    spaced_kmer_mode: bool | None = None,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    add_self_matches: bool | None = None,
    filter_hits: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Heuristic gapped local alignment seeded by k-mer matches (`mmseqs alignbykmer`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - result_db: prefilter result database.
    - alignment_db: new alignment database.
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - spaced_kmer_mode: Use spaced k-mers instead of consecutive k-mers.
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - min_aln_len: Minimum alignment length (default 0).
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - filter_hits: Filter hits by sequence identity and coverage.
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
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_bool("--filter-hits", filter_hits, "filter_hits")
    result = core.run(
        "mmseqs_alignbykmer",
        "alignbykmer",
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


def mmseqs_extractorfs(
    sequence_db: str,
    orf_db: str,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    max_gaps: int | None = None,
    contig_start_mode: int | None = None,
    contig_end_mode: int | None = None,
    orf_start_mode: int | None = None,
    forward_frames: str | None = None,
    reverse_frames: str | None = None,
    translation_table: int | None = None,
    translate: bool | None = None,
    use_all_table_starts: bool | None = None,
    id_offset: int | None = None,
    create_lookup: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Six-frame extraction of open reading frames from nucleotide sequences (`mmseqs
    extractorfs`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: nucleotide sequence database.
    - orf_db: new ORF database (DNA, or protein with translate=True).
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
    - translate: Translate the ORFs to amino acids (default off: DNA ORFs).
    - use_all_table_starts: Accept every start codon of the genetic code, not only ATG.
    - id_offset: Offset of the numeric ids in the index file (default 0).
    - create_lookup: Create a lookup file for the produced database (can be large; default
      off).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    orf_db_path = core.output_database(orf_db, "orf_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(orf_db_path)
    options: list[str] = []
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
    options += core.opt_bool("--translate", translate, "translate")
    options += core.opt_bool("--use-all-table-starts", use_all_table_starts, "use_all_table_starts")
    options += core.opt_int("--id-offset", id_offset, "id_offset", 0, None)
    options += core.opt_bool("--create-lookup", create_lookup, "create_lookup")
    result = core.run(
        "mmseqs_extractorfs",
        "extractorfs",
        positionals,
        options,
        databases={"orf_db": orf_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_extractframes(
    sequence_db: str,
    frame_db: str,
    *,
    forward_frames: str | None = None,
    reverse_frames: str | None = None,
    translation_table: int | None = None,
    translate: bool | None = None,
    create_lookup: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the reading frames of nucleotide sequences, translated or not (`mmseqs
    extractframes`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: nucleotide sequence database.
    - frame_db: new database of the frames.
    - forward_frames: Comma-separated forward-strand frames to extract, for example '1,2,3'
      (default 1,2,3).
    - reverse_frames: Comma-separated reverse-strand frames to extract, for example '1,3'
      (default 1,2,3).
    - translation_table: NCBI genetic code used for translation (default 1, standard).
    - translate: Translate the ORFs to amino acids (default off: DNA ORFs).
    - create_lookup: Create a lookup file for the produced database (can be large; default
      off).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    frame_db_path = core.output_database(frame_db, "frame_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(frame_db_path)
    options: list[str] = []
    options += core.frames("--forward-frames", forward_frames, "forward_frames")
    options += core.frames("--reverse-frames", reverse_frames, "reverse_frames")
    options += core.opt_int("--translation-table", translation_table, "translation_table", 1, 31)
    options += core.opt_bool("--translate", translate, "translate")
    options += core.opt_bool("--create-lookup", create_lookup, "create_lookup")
    result = core.run(
        "mmseqs_extractframes",
        "extractframes",
        positionals,
        options,
        databases={"frame_db": frame_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_orftocontig(
    contig_db: str,
    orf_db: str,
    alignment_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the ORF locations on their contigs in alignment format (`mmseqs orftocontig`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - contig_db: nucleotide sequence database of the contigs.
    - orf_db: database of the extracted ORFs (extractorfs output); its header database
      <orf_db>_h is read.
    - alignment_db: new alignment database.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    contig_db_path = core.input_database(contig_db, "contig_db")
    orf_db_path = core.input_database(orf_db, "orf_db")
    alignment_db_path = core.output_database(alignment_db, "alignment_db")
    positionals: list[Path] = []
    positionals.append(contig_db_path)
    positionals.append(orf_db_path)
    positionals.append(alignment_db_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_orftocontig",
        "orftocontig",
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


def mmseqs_recoverlongestorf(
    orf_db: str,
    result_db: str,
    output: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recover the longest ORF per query contig of a result for taxonomic assignment (`mmseqs
    recoverlongestorf`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - orf_db: ORF database.
    - result_db: result database of the ORFs.
    - output: new TSV file.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    The result adds ``rows`` (up to 1000 tab-separated rows as lists of text).
    """
    orf_db_path = core.input_database(orf_db, "orf_db")
    result_db_path = core.input_database(result_db, "result_db")
    output_path = core.output_file(output, "output")
    positionals: list[Path] = []
    positionals.append(orf_db_path)
    positionals.append(result_db_path)
    positionals.append(output_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_recoverlongestorf",
        "recoverlongestorf",
        positionals,
        options,
        databases={},
        files={"output": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["rows"] = core.table(output_path)
    return result


def mmseqs_splitsequence(
    sequence_db: str,
    split_db: str,
    *,
    sequence_overlap: int | None = None,
    max_seq_len: int | None = None,
    sequence_split_mode: int | None = None,
    headers_split_mode: int | None = None,
    create_lookup: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split long sequences into overlapping pieces (`mmseqs splitsequence`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - split_db: new sequence database of the pieces.
    - sequence_overlap: Overlap between consecutive pieces (default 300).
    - max_seq_len: Maximum sequence length (default 65535).
    - sequence_split_mode: 0 copy the data, 1 soft-link the data and write a new index
      (default 1).
    - headers_split_mode: Headers of the pieces: 0 with the split position, 1 the original
      header (default 0).
    - create_lookup: Create a lookup file for the produced database (can be large; default
      off).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    split_db_path = core.output_database(split_db, "split_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(split_db_path)
    options: list[str] = []
    options += core.opt_int("--sequence-overlap", sequence_overlap, "sequence_overlap", 0, None)
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    options += core.opt_choice(
        "--sequence-split-mode", sequence_split_mode, (0, 1), "sequence_split_mode"
    )
    options += core.opt_choice(
        "--headers-split-mode", headers_split_mode, (0, 1), "headers_split_mode"
    )
    options += core.opt_bool("--create-lookup", create_lookup, "create_lookup")
    result = core.run(
        "mmseqs_splitsequence",
        "splitsequence",
        positionals,
        options,
        databases={"split_db": split_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_masksequence(
    sequence_db: str,
    masked_db: str,
    *,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    mask_n_repeat: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Soft-mask low-complexity regions of a sequence database with tantan (`mmseqs
    masksequence`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - masked_db: new sequence database with the masked residues in lower case.
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - mask_n_repeat: Mask letters that repeat more than this many times in a row (default 0,
      off).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    masked_db_path = core.output_database(masked_db, "masked_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(masked_db_path)
    options: list[str] = []
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    options += core.opt_int("--mask-n-repeat", mask_n_repeat, "mask_n_repeat", 0, None)
    result = core.run(
        "mmseqs_masksequence",
        "masksequence",
        positionals,
        options,
        databases={"masked_db": masked_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_offsetalignment(
    query_db: str,
    query_orf_db: str,
    target_db: str,
    target_orf_db: str,
    alignment_db: str,
    offset_db: str,
    *,
    search_type: int | None = None,
    chain_alignments: bool | None = None,
    merge_query: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shift the alignments of ORFs or split sequences by their offset on the original sequence
    (`mmseqs offsetalignment`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: original query sequence database.
    - query_orf_db: query ORF or piece database.
    - target_db: original target sequence database.
    - target_orf_db: target ORF or piece database.
    - alignment_db: alignment database of the ORFs.
    - offset_db: new alignment database in original coordinates.
    - search_type: Search type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4
      translated nucleotide alignment (default 0).
    - chain_alignments: Chain overlapping alignments.
    - merge_query: Combine the ORFs or pieces of a query into one entry (default on).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    query_orf_db_path = core.input_database(query_orf_db, "query_orf_db")
    target_db_path = core.input_database(target_db, "target_db")
    target_orf_db_path = core.input_database(target_orf_db, "target_orf_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    offset_db_path = core.output_database(offset_db, "offset_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(query_orf_db_path)
    positionals.append(target_db_path)
    positionals.append(target_orf_db_path)
    positionals.append(alignment_db_path)
    positionals.append(offset_db_path)
    options: list[str] = []
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_bool("--chain-alignments", chain_alignments, "chain_alignments")
    options += core.opt_bool("--merge-query", merge_query, "merge_query")
    result = core.run(
        "mmseqs_offsetalignment",
        "offsetalignment",
        positionals,
        options,
        databases={"offset_db": offset_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_proteinaln2nucl(
    nucleotide_query_db: str,
    nucleotide_target_db: str,
    protein_query_db: str,
    protein_target_db: str,
    alignment_db: str,
    nucleotide_alignment_db: str,
    *,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map protein alignments back to nucleotide alignments (`mmseqs proteinaln2nucl`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - nucleotide_query_db: nucleotide query database.
    - nucleotide_target_db: nucleotide target database.
    - protein_query_db: protein query database (translated).
    - protein_target_db: protein target database (translated).
    - alignment_db: protein alignment database.
    - nucleotide_alignment_db: new nucleotide alignment database.
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    nucleotide_query_db_path = core.input_database(nucleotide_query_db, "nucleotide_query_db")
    nucleotide_target_db_path = core.input_database(nucleotide_target_db, "nucleotide_target_db")
    protein_query_db_path = core.input_database(protein_query_db, "protein_query_db")
    protein_target_db_path = core.input_database(protein_target_db, "protein_target_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    nucleotide_alignment_db_path = core.output_database(
        nucleotide_alignment_db, "nucleotide_alignment_db"
    )
    positionals: list[Path] = []
    positionals.append(nucleotide_query_db_path)
    positionals.append(nucleotide_target_db_path)
    positionals.append(protein_query_db_path)
    positionals.append(protein_target_db_path)
    positionals.append(alignment_db_path)
    positionals.append(nucleotide_alignment_db_path)
    options: list[str] = []
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    result = core.run(
        "mmseqs_proteinaln2nucl",
        "proteinaln2nucl",
        positionals,
        options,
        databases={"nucleotide_alignment_db": nucleotide_alignment_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


__all__ = [
    "mmseqs_alignall",
    "mmseqs_alignbykmer",
    "mmseqs_extractframes",
    "mmseqs_extractorfs",
    "mmseqs_fwbw",
    "mmseqs_gappedprefilter",
    "mmseqs_kmermatcher",
    "mmseqs_kmersearch",
    "mmseqs_masksequence",
    "mmseqs_offsetalignment",
    "mmseqs_orftocontig",
    "mmseqs_prefilter",
    "mmseqs_proteinaln2nucl",
    "mmseqs_recoverlongestorf",
    "mmseqs_rescorediagonal",
    "mmseqs_splitsequence",
    "mmseqs_transitivealign",
    "mmseqs_ungappedprefilter",
]

_register_public_verbose(globals(), tuple(__all__))
