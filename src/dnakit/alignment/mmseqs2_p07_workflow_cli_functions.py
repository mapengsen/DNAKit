"""MMseqs2 database, search, clustering and export workflows.

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

import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment import _mmseqs2_p07_core as core

_DEFAULT_COLUMNS = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"


def mmseqs_linsearch(
    query: str,
    target: str,
    output: str,
    *,
    input_kind: str = "fasta",
    search_type: int | None = None,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    kmer_per_seq: int | None = None,
    alignment_mode: int | None = None,
    format_mode: int | None = None,
    format_output: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Linear-time, less sensitive sequence search (``mmseqs easy-linsearch``, ``mmseqs
    linsearch``).

    The search matches bottom-m-hashed k-mers and is built for large, highly similar sets. It covers
    two entry points of MMseqs2 through ``input_kind``: ``"fasta"`` runs ``easy-linsearch`` (FASTA
    query and FASTA or database target, a new alignment table as ``output``) and ``"database"`` runs
    ``linsearch`` (query and target are database prefixes, ``output`` is the prefix of a new
    alignment
    database). In database mode the k-mer index ``<target>.linidx`` must exist; create it with
    ``mmseqs_createlinindex``. Nucleotide input needs ``search_type=3``; on protein input the linear
    search finds only near-identical matches.

    Parameters:
    - query, target, output: see above.
    - input_kind: "fasta" (easy-linsearch, default) or "database" (linsearch).
    - search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
      nucleotide alignment.
    - e_value, min_seq_id, coverage, cov_mode, min_aln_len, seq_id_mode: match filters as in MMseqs2
      (defaults 0.001, 0, 0, 0, 0, 0).
    - kmer_per_seq: number of k-mers selected per sequence (MMseqs2 default 21).
    - alignment_mode: 0 automatic, 1 score and end position only, 2 also start position and
      coverage,
      3 also sequence identity, 4 ungapped only.
    - format_mode, format_output: table format and columns (fasta mode only), as for convertalis.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.

    The result adds ``records`` (up to 1000 rows) in fasta mode with format_mode 0 or 4.
    """
    kind = core.choose_kind(input_kind)
    options: list[str] = []
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_int("--kmer-per-seq", kmer_per_seq, "kmer_per_seq", 0, None)
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    if kind == "database":
        core.reject_with(
            format_mode is not None or format_output is not None,
            "format_mode and format_output apply to input_kind='fasta'; use mmseqs_convertalis.",
        )
        query_path = core.input_database(query, "query")
        target_path = core.input_database(target, "target")
        output_path = core.output_database(output, "output")
        core.reject_with(
            not Path(f"{target_path}.linidx.dbtype").is_file(),
            f"The k-mer index {target_path}.linidx is missing; run mmseqs_createlinindex first.",
        )
        return core.run(
            "mmseqs_linsearch",
            "linsearch",
            [query_path, target_path, output_path],
            options,
            databases={"alignment_db": output_path},
            needs_tmp=True,
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    options += core.opt_choice("--format-mode", format_mode, (0, 1, 2, 3, 4), "format_mode")
    options += core.columns("--format-output", format_output, "format_output")
    query_path = core.input_file(query, "query")
    target_path = core.fasta_or_database(target, "target")
    output_path = core.output_file(output, "output")
    result = core.run(
        "mmseqs_linsearch",
        "easy-linsearch",
        [query_path, target_path, output_path],
        options,
        files={"output": output_path},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    if format_mode in (None, 0):
        result["records"] = core.table(output_path, (format_output or _DEFAULT_COLUMNS).split(","))
    elif format_mode == 4:
        rows = core.table(output_path)
        result["records"] = (
            [dict(zip(rows[0], row, strict=False)) for row in rows[1:]] if rows else []
        )
    return result


def mmseqs_linclust(
    sequences: str | Sequence[str],
    output: str,
    *,
    input_kind: str = "fasta",
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    e_value: float | None = None,
    cluster_mode: int | None = None,
    similarity_type: int | None = None,
    kmer_size: int | None = None,
    kmer_per_seq: int | None = None,
    alignment_mode: int | None = None,
    spaced_kmer_mode: bool | None = None,
    adjust_kmer_len: bool | None = None,
    ignore_multi_kmer: bool | None = None,
    hash_shift: int | None = None,
    linclust_version: int | None = None,
    mask: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Linear-time, less sensitive clustering (``mmseqs easy-linclust``, ``mmseqs linclust``).

    ``input_kind="fasta"`` runs ``easy-linclust``: ``sequences`` are FASTA file(s) and
    ``output`` is a
    new file prefix; ``<output>_cluster.tsv``, ``<output>_rep_seq.fasta`` and
    ``<output>_all_seqs.fasta``
    are written. ``input_kind="database"`` runs ``linclust``: ``sequences`` is a sequence database
    prefix and ``output`` the prefix of a new cluster database. MMseqs2's own defaults apply (for
    easy-linclust identity 0 and coverage 0.8 with cov_mode 0).

    Parameters:
    - sequences, output, input_kind: see above.
    - min_seq_id, coverage, cov_mode, e_value: match filters (coverage default 0.8).
    - cluster_mode: 0 set cover (greedy), 1 connected component, 2 or 3 greedy by length.
    - similarity_type: 1 alignment score, 2 sequence identity (default 2).
    - kmer_size, kmer_per_seq, spaced_kmer_mode, adjust_kmer_len, ignore_multi_kmer, hash_shift,
      linclust_version, mask: k-mer stage options (kmer_size 0 = automatic).
    - alignment_mode: 0 automatic, 1 score and end position, 2 also start position and
      coverage, 3 also identity.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.

    The fasta mode adds ``clusters`` (representative to member ids), ``cluster_count`` and
    ``sequence_count``.
    """
    kind = core.choose_kind(input_kind)
    options: list[str] = []
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_choice("--cluster-mode", cluster_mode, (0, 1, 2, 3), "cluster_mode")
    options += core.opt_choice("--similarity-type", similarity_type, (1, 2), "similarity_type")
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_int("--kmer-per-seq", kmer_per_seq, "kmer_per_seq", 0, None)
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_bool("--adjust-kmer-len", adjust_kmer_len, "adjust_kmer_len")
    options += core.opt_bool("--ignore-multi-kmer", ignore_multi_kmer, "ignore_multi_kmer")
    options += core.opt_int("--hash-shift", hash_shift, "hash_shift", 0, None)
    options += core.opt_choice("--linclust-version", linclust_version, (1, 2), "linclust_version")
    options += core.opt_bool("--mask", mask, "mask")
    if kind == "database":
        sequences_path = core.input_database(
            sequences if isinstance(sequences, str) else "", "sequences"
        )
        output_path = core.output_database(output, "output")
        return core.run(
            "mmseqs_linclust",
            "linclust",
            [sequences_path, output_path],
            options,
            databases={"cluster_db": output_path},
            needs_tmp=True,
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    inputs = core.fasta_files(sequences, "sequences")
    prefix, family = core.output_family(
        output,
        {
            "clusters": "_cluster.tsv",
            "representatives": "_rep_seq.fasta",
            "all_sequences": "_all_seqs.fasta",
        },
        "output",
    )
    result = core.run(
        "mmseqs_linclust",
        "easy-linclust",
        [*inputs, prefix],
        options,
        files=family,
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    clusters: dict[str, list[str]] = {}
    for row in core.table(family["clusters"]):
        if len(row) >= 2:
            clusters.setdefault(row[0], []).append(row[1])
    result["clusters"] = clusters
    result["cluster_count"] = len(clusters)
    result["sequence_count"] = sum(len(members) for members in clusters.values())
    return result


def mmseqs_rbh(
    query: str,
    target: str,
    output: str,
    *,
    input_kind: str = "fasta",
    sensitivity: float | None = None,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    max_seqs: int | None = None,
    kmer_size: int | None = None,
    search_type: int | None = None,
    alignment_mode: int | None = None,
    format_mode: int | None = None,
    format_output: str | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reciprocal best-hit search of two sequence sets (``mmseqs easy-rbh``, ``mmseqs rbh``).

    A pair is reported when each sequence is the other's best hit in a search in both directions.
    ``input_kind="fasta"`` runs ``easy-rbh``: FASTA (or database) inputs and a new alignment
    table as
    ``output``. ``input_kind="database"`` runs ``rbh``: database prefixes and a new alignment
    database;
    ``rbh`` filters the merged two-way result with ``result2rbh`` itself.

    Parameters:
    - query, target, output, input_kind: see above.
    - sensitivity: search sensitivity (MMseqs2 default 5.7 for this workflow).
    - e_value, min_seq_id, coverage, cov_mode, min_aln_len, max_seqs, kmer_size: search and
      match filters.
    - search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
      nucleotide alignment.
    - alignment_mode: 0 automatic, 1 score and end position, 2 also start position and
      coverage, 3 also identity.
    - format_mode, format_output: table format and columns (fasta mode only).
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.

    The fasta mode adds ``records`` (up to 1000 rows) for format_mode 0 and 4 and ``pair_count``.
    """
    kind = core.choose_kind(input_kind)
    options: list[str] = []
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    if kind == "database":
        core.reject_with(
            format_mode is not None or format_output is not None,
            "format_mode and format_output apply to input_kind='fasta'; use mmseqs_convertalis.",
        )
        query_path = core.input_database(query, "query")
        target_path = core.input_database(target, "target")
        output_path = core.output_database(output, "output")
        return core.run(
            "mmseqs_rbh",
            "rbh",
            [query_path, target_path, output_path],
            options,
            databases={"alignment_db": output_path},
            needs_tmp=True,
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    options += core.opt_choice("--format-mode", format_mode, (0, 1, 2, 3, 4), "format_mode")
    options += core.columns("--format-output", format_output, "format_output")
    query_path = core.input_file(query, "query")
    target_path = core.fasta_or_database(target, "target")
    output_path = core.output_file(output, "output")
    result = core.run(
        "mmseqs_rbh",
        "easy-rbh",
        [query_path, target_path, output_path],
        options,
        files={"output": output_path},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["pair_count"] = result["line_counts"]["output"] - (1 if format_mode == 4 else 0)
    if format_mode in (None, 0):
        result["records"] = core.table(output_path, (format_output or _DEFAULT_COLUMNS).split(","))
    elif format_mode == 4:
        rows = core.table(output_path)
        result["records"] = (
            [dict(zip(rows[0], row, strict=False)) for row in rows[1:]] if rows else []
        )
    return result


def mmseqs_createlinindex(
    sequence_db: str,
    *,
    kmer_size: int | None = None,
    kmer_per_seq: int | None = None,
    min_seq_id: float | None = None,
    search_type: int | None = None,
    spaced_kmer_mode: bool | None = None,
    adjust_kmer_len: bool | None = None,
    ignore_multi_kmer: bool | None = None,
    hash_shift: int | None = None,
    mask: bool | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the bottom-m-hash k-mer index of a sequence database (``mmseqs createlinindex``).

    The index ``<sequence_db>.linidx`` is what ``mmseqs_kmersearch`` and the database mode of
    ``mmseqs_linsearch`` search against. It is written next to the sequence database and must
    not exist yet.

    Parameters:
    - sequence_db: prefix of the sequence database to index.
    - kmer_size: k-mer length (0 = automatic). kmer_per_seq: k-mers selected per sequence
      (default 21).
    - min_seq_id, search_type, spaced_kmer_mode, adjust_kmer_len, ignore_multi_kmer,
      hash_shift, mask,
      max_seq_len: as for mmseqs_kmermatcher; search_type 3 for nucleotides.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.
    """
    database = core.input_database(sequence_db, "sequence_db")
    index = core.output_database(f"{database}.linidx", "the k-mer index")
    options: list[str] = []
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_int("--kmer-per-seq", kmer_per_seq, "kmer_per_seq", 0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_bool("--spaced-kmer-mode", spaced_kmer_mode, "spaced_kmer_mode")
    options += core.opt_bool("--adjust-kmer-len", adjust_kmer_len, "adjust_kmer_len")
    options += core.opt_bool("--ignore-multi-kmer", ignore_multi_kmer, "ignore_multi_kmer")
    options += core.opt_int("--hash-shift", hash_shift, "hash_shift", 0, None)
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    return core.run(
        "mmseqs_createlinindex",
        "createlinindex",
        [database],
        options,
        databases={"kmer_index_db": index},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_result2msa(
    query_db: str,
    target_db: str,
    result_db: str,
    msa_db: str,
    *,
    dna: bool = False,
    msa_format_mode: int | None = None,
    filter_msa: bool | None = None,
    filter_min_enable: int | None = None,
    max_seq_id: float | None = None,
    qid: str | None = None,
    qsc: float | None = None,
    cov: float | None = None,
    diff: int | None = None,
    allow_deletion: bool | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    comp_bias_corr: bool | None = None,
    skip_query: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multiple sequence alignment per query (``mmseqs result2msa``, ``mmseqs result2dnamsa``).

    ``dna=False`` runs ``result2msa`` (profile-style MSA with insertions, any alphabet, filters and
    output formats). ``dna=True`` runs ``result2dnamsa`` (nucleotide MSA without insertions in the
    query); it only accepts ``skip_query`` besides the databases. The result is an MSA database; use
    ``mmseqs_convert2fasta`` or ``mmseqs_createtsv`` on it, or msa_format_mode 4 or 5 for flat
    formats.

    Parameters:
    - query_db, target_db, result_db, msa_db: databases (msa_db is new).
    - dna: use result2dnamsa.
    - msa_format_mode: 0 binary cA3M, 1 cA3M with consensus, 2 aligned FASTA DB, 3 aligned
      FASTA with
      header summary, 4 Stockholm flat file, 5 A3M, 6 A3M with alignment info (default 2).
    - filter_msa, filter_min_enable, max_seq_id, qid, qsc, cov, diff, allow_deletion: MSA filters
      (defaults off, 0, 0.9, none, -20, 0, 1000... as in MMseqs2).
    - gap_open, gap_extend, comp_bias_corr: realignment options.
    - skip_query: leave the query out of the MSA.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.
    """
    query_path = core.input_database(query_db, "query_db")
    target_path = core.input_database(target_db, "target_db")
    result_path = core.input_database(result_db, "result_db")
    output_path = core.output_database(msa_db, "msa_db")
    options: list[str] = []
    if dna:
        core.reject_with(
            any(
                value is not None
                for value in (
                    msa_format_mode,
                    filter_msa,
                    filter_min_enable,
                    max_seq_id,
                    qid,
                    qsc,
                    cov,
                    diff,
                    allow_deletion,
                    gap_open,
                    gap_extend,
                    comp_bias_corr,
                )
            ),
            "dna=True (result2dnamsa) accepts only skip_query.",
        )
        command = "result2dnamsa"
    else:
        command = "result2msa"
        options += core.opt_choice(
            "--msa-format-mode", msa_format_mode, (0, 1, 2, 3, 4, 5, 6), "msa_format_mode"
        )
        options += core.opt_bool("--filter-msa", filter_msa, "filter_msa")
        options += core.opt_int(
            "--filter-min-enable", filter_min_enable, "filter_min_enable", 0, None
        )
        options += core.opt_float("--max-seq-id", max_seq_id, "max_seq_id", 0.0, 1.0)
        options += core.opt_text("--qid", qid, "qid", r"^[0-9.]+(,[0-9.]+)*$")
        options += core.opt_float("--qsc", qsc, "qsc", -50.0, 100.0)
        options += core.opt_float("--cov", cov, "cov", 0.0, 1.0)
        options += core.opt_int("--diff", diff, "diff", 0, None)
        options += core.opt_bool("--allow-deletion", allow_deletion, "allow_deletion")
        options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
        options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
        options += core.opt_bool("--comp-bias-corr", comp_bias_corr, "comp_bias_corr")
    options += core.opt_bool("--skip-query", skip_query, "skip_query")
    return core.run(
        "mmseqs_result2msa",
        command,
        [query_path, target_path, result_path, output_path],
        options,
        databases={"msa_db": output_path},
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )


def mmseqs_maskbygff(
    gff_file: str,
    sequence_db: str,
    masked_db: str,
    *,
    gff_type: str | None = None,
    id_offset: int | None = None,
    names_as_ids: bool = True,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask the regions of selected GFF3 features in a sequence database (``mmseqs maskbygff``).

    The masked residues are written as lower case or X as MMseqs2 does. MMseqs2 matches the first
    GFF3 column against the numeric database keys (the first column of
    ``<sequence_db>.lookup``), not
    against sequence names. With ``names_as_ids=True`` (default) the facade therefore writes a
    temporary copy of the GFF3 file in which sequence names are replaced by their database keys;
    lines whose sequence name is not in the lookup are kept unchanged. The original file is never
    modified.

    Parameters:
    - gff_file: GFF3 file. sequence_db: sequence database of the genome. masked_db: new database.
    - gff_type: comma-separated feature types to select (default all).
    - id_offset: offset of the numeric ids in the index file (default 0).
    - names_as_ids: translate sequence names to database keys (see above); False passes the
      file as is.
    - threads: CPU cores (default 1). mmseqs_path: path of the mmseqs executable.
    """
    gff_path = core.input_file(gff_file, "gff_file")
    database = core.input_database(sequence_db, "sequence_db")
    output_path = core.output_database(masked_db, "masked_db")
    options: list[str] = []
    options += core.opt_text("--gff-type", gff_type, "gff_type", r"^[A-Za-z0-9_.:,-]+$")
    options += core.opt_int("--id-offset", id_offset, "id_offset", 0, None)
    with tempfile.TemporaryDirectory(prefix=".mmseqs_gff_", dir=output_path.parent) as scratch:
        used = gff_path
        if names_as_ids:
            lookup = Path(f"{database}.lookup")
            core.reject_with(
                not lookup.is_file(), f"{lookup} is needed to translate sequence names."
            )
            keys: dict[str, str] = {}
            with lookup.open(encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 2:
                        keys[fields[1]] = fields[0]
            used = Path(scratch) / "features.gff3"
            with (
                gff_path.open(encoding="utf-8") as source,
                used.open("w", encoding="utf-8") as target,
            ):
                for line in source:
                    fields = line.rstrip("\n").split("\t")
                    if not line.startswith("#") and len(fields) > 1 and fields[0] in keys:
                        fields[0] = keys[fields[0]]
                    target.write("\t".join(fields) + "\n")
        result = core.run(
            "mmseqs_maskbygff",
            "maskbygff",
            [used, database, output_path],
            options,
            databases={"masked_db": output_path},
            threads=threads,
            mmseqs_path=mmseqs_path,
            verbose=verbose,
        )
    return result


def mmseqs_createdb(
    fasta_files: str | Sequence[str],
    database: str,
    *,
    dbtype: int | None = None,
    shuffle: bool | None = None,
    createdb_mode: int | None = None,
    write_lookup: bool | None = None,
    mask: bool | None = None,
    mask_prob: float | None = None,
    mask_lower_case: bool | None = None,
    id_offset: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert FASTA/FASTQ files into an MMseqs2 sequence database (`mmseqs createdb`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - fasta_files: FASTA or FASTQ file(s), plain, gzip or bzip2.
    - database: prefix of the new sequence database.
    - dbtype: Database type: 0 automatic, 1 amino acid, 2 nucleotide (default 0).
    - shuffle: Shuffle the input order inside the database (default on).
    - createdb_mode: 0 copy the data, 1 soft-link the FASTA and write a new index (single-
      line FASTA/FASTQ only; default 0).
    - write_lookup: Write the .lookup file that maps internal ids to FASTA ids (default on).
    - mask: Mask low-complexity regions with tantan before the k-mer stage.
    - mask_prob: Probability above which a position is masked (default 0.9).
    - mask_lower_case: Exclude lower-case residues from the k-mer search.
    - id_offset: Offset of the numeric ids in the index file (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    fasta_files_path = core.fasta_files(fasta_files, "fasta_files")
    database_path = core.output_database(database, "database")
    positionals: list[Path] = []
    positionals.extend(fasta_files_path)
    positionals.append(database_path)
    options: list[str] = []
    options += core.opt_choice("--dbtype", dbtype, (0, 1, 2), "dbtype")
    options += core.opt_bool("--shuffle", shuffle, "shuffle")
    options += core.opt_choice("--createdb-mode", createdb_mode, (0, 1), "createdb_mode")
    options += core.opt_bool("--write-lookup", write_lookup, "write_lookup")
    options += core.opt_bool("--mask", mask, "mask")
    options += core.opt_float("--mask-prob", mask_prob, "mask_prob", 0.0, 1.0)
    options += core.opt_bool("--mask-lower-case", mask_lower_case, "mask_lower_case")
    options += core.opt_int("--id-offset", id_offset, "id_offset", 0, None)
    result = core.run(
        "mmseqs_createdb",
        "createdb",
        positionals,
        options,
        databases={"database": database_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_map(
    query_db: str,
    target_db: str,
    alignment_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    min_aln_len: int | None = None,
    seq_id_mode: int | None = None,
    add_backtrace: bool | None = None,
    sensitivity: float | None = None,
    kmer_size: int | None = None,
    max_seqs: int | None = None,
    search_type: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map nearly identical sequences of a query database onto a target database (`mmseqs map`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
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
    - add_backtrace: Store the backtrace string so that convertalis can write alignments.
    - sensitivity: Search sensitivity, 1.0 fast to 7.5 sensitive (default 4.0 for the
      prefilter, workflow-specific otherwise).
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - max_seqs: Maximum number of prefilter results per query (affects sensitivity; default
      300).
    - search_type: Search type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4
      translated nucleotide alignment (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    alignment_db_path = core.output_database(alignment_db, "alignment_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(alignment_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_int("--min-aln-len", min_aln_len, "min_aln_len", 0, None)
    options += core.opt_choice("--seq-id-mode", seq_id_mode, (0, 1, 2), "seq_id_mode")
    options += core.opt_bool("-a", add_backtrace, "add_backtrace")
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_int("--max-seqs", max_seqs, "max_seqs", 1, None)
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    result = core.run(
        "mmseqs_map",
        "map",
        positionals,
        options,
        databases={"alignment_db": alignment_db_path},
        files={},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_clusterupdate(
    old_sequence_db: str,
    new_sequence_db: str,
    old_cluster_db: str,
    mapped_sequence_db: str,
    new_cluster_db: str,
    *,
    e_value: float | None = None,
    min_seq_id: float | None = None,
    coverage: float | None = None,
    cov_mode: int | None = None,
    cluster_mode: int | None = None,
    sensitivity: float | None = None,
    kmer_size: int | None = None,
    alignment_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Update an earlier clustering with new sequences (`mmseqs clusterupdate`); MMseqs2 deletes
    the old cluster database.

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - old_sequence_db: sequence database that was clustered before.
    - new_sequence_db: sequence database with the old and the new sequences.
    - old_cluster_db: previous cluster database.
    - mapped_sequence_db: new sequence database whose keys are mapped to the old ones.
    - new_cluster_db: new cluster database.
    - e_value: Maximum E-value of a reported match (MMseqs2 default 0.001).
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - coverage: Minimum fraction of aligned residues, defined by cov_mode (default 0).
    - cov_mode: Coverage definition for coverage: 0 query and target, 1 target, 2 query, 3
      target length at least x% of query length, 4 query length at least x% of target
      length, 5 shorter sequence at least x% of the other (MMseqs2 default 0).
    - cluster_mode: Clustering: 0 set cover (greedy), 1 connected component (BLASTclust), 2
      or 3 greedy by sequence length (CD-HIT like) (default 0).
    - sensitivity: Search sensitivity, 1.0 fast to 7.5 sensitive (default 4.0 for the
      prefilter, workflow-specific otherwise).
    - kmer_size: k-mer length; 0 lets MMseqs2 choose the optimum (default 0).
    - alignment_mode: How the alignment is computed: 0 automatic, 1 score and end position
      only, 2 also start position and coverage, 3 also sequence identity, 4 ungapped only.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    MMseqs2 deletes ``old_cluster_db`` while it works; the result reports this as
    ``old_cluster_db_removed``. Copy the database first if you still need it.
    """
    old_sequence_db_path = core.input_database(old_sequence_db, "old_sequence_db")
    new_sequence_db_path = core.input_database(new_sequence_db, "new_sequence_db")
    old_cluster_db_path = core.input_database(old_cluster_db, "old_cluster_db")
    mapped_sequence_db_path = core.output_database(mapped_sequence_db, "mapped_sequence_db")
    new_cluster_db_path = core.output_database(new_cluster_db, "new_cluster_db")
    positionals: list[Path] = []
    positionals.append(old_sequence_db_path)
    positionals.append(new_sequence_db_path)
    positionals.append(old_cluster_db_path)
    positionals.append(mapped_sequence_db_path)
    positionals.append(new_cluster_db_path)
    options: list[str] = []
    options += core.opt_float("-e", e_value, "e_value", 0.0, None)
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_float("-c", coverage, "coverage", 0.0, 1.0)
    options += core.opt_choice("--cov-mode", cov_mode, (0, 1, 2, 3, 4, 5), "cov_mode")
    options += core.opt_choice("--cluster-mode", cluster_mode, (0, 1, 2, 3), "cluster_mode")
    options += core.opt_float("-s", sensitivity, "sensitivity", 0.0, None)
    options += core.opt_int("-k", kmer_size, "kmer_size", 0, None)
    options += core.opt_choice(
        "--alignment-mode", alignment_mode, (0, 1, 2, 3, 4), "alignment_mode"
    )
    result = core.run(
        "mmseqs_clusterupdate",
        "clusterupdate",
        positionals,
        options,
        databases={
            "mapped_sequence_db": mapped_sequence_db_path,
            "new_cluster_db": new_cluster_db_path,
        },
        files={},
        needs_tmp=True,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    result["old_cluster_db_removed"] = not Path(f"{old_cluster_db_path}.dbtype").exists()
    return result


def mmseqs_convertalis(
    query_db: str,
    target_db: str,
    alignment_db: str,
    output: str,
    *,
    format_mode: int | None = None,
    format_output: str | None = None,
    search_type: int | None = None,
    translation_table: int | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export an alignment database to a BLAST-tab, SAM or custom-column table (`mmseqs
    convertalis`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - alignment_db: alignment database.
    - output: new output table.
    - format_mode: Output format: 0 BLAST-tab, 1 SAM, 2 BLAST-tab with query and target
      length, 3 pretty HTML, 4 BLAST-tab with a header line (default 0).
    - format_output: Comma-separated output columns for formats 0 and 4, e.g.
      'query,target,fident,alnlen,evalue,bits' (MMseqs2 offers query, target, evalue,
      gapopen, pident, fident, nident, qstart, qend, qlen, tstart, tend, tlen, alnlen, raw,
      bits, cigar, qseq, tseq, qheader, theader, qaln, taln, qframe, tframe, mismatch, qcov,
      tcov, qset, qsetid, tset, tsetid, taxid, taxname, taxlineage, qorfstart, qorfend,
      torfstart, torfend, ppos).
    - search_type: Search type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4
      translated nucleotide alignment (default 0).
    - translation_table: NCBI genetic code used for translation (default 1, standard).
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    The result adds ``records`` (up to 1000 rows as dictionaries) for format_mode 0 and 4.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    output_path = core.output_file(output, "output")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(alignment_db_path)
    positionals.append(output_path)
    options: list[str] = []
    options += core.opt_choice("--format-mode", format_mode, (0, 1, 2, 3, 4), "format_mode")
    options += core.columns("--format-output", format_output, "format_output")
    options += core.opt_choice("--search-type", search_type, (0, 1, 2, 3, 4), "search_type")
    options += core.opt_int("--translation-table", translation_table, "translation_table", 1, 31)
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    result = core.run(
        "mmseqs_convertalis",
        "convertalis",
        positionals,
        options,
        databases={},
        files={"output": output_path},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    if format_mode in (None, 0):
        names = (format_output or _DEFAULT_COLUMNS).split(",")
        result["records"] = core.table(output_path, names)
    elif format_mode == 4:
        rows = core.table(output_path)
        result["records"] = (
            [dict(zip(rows[0], row, strict=False)) for row in rows[1:]] if rows else []
        )
    return result


def mmseqs_createtsv(
    query_db: str,
    result_db: str,
    output: str,
    *,
    target_db: str | None = None,
    first_seq_as_repr: bool | None = None,
    target_column: int | None = None,
    full_header: bool | None = None,
    idx_seq_src: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a result database (clusters, alignments, taxonomy) to a tab-separated file (`mmseqs
    createtsv`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database; omit it for results that have no target (clusters
      of one database).
    - result_db: result database.
    - output: new TSV file.
    - first_seq_as_repr: Use the first sequence of each cluster as its representative.
    - target_column: Column holding the target id (default 1; 0 when the results have no
      target id).
    - full_header: Replace database ids by the full FASTA header.
    - idx_seq_src: Take sequence names from 0 automatic, 1 split or translated sequences, 2
      input sequences (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    The result adds ``rows`` (up to 1000 tab-separated rows as lists of text).
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = None if target_db is None else core.input_database(target_db, "target_db")
    result_db_path = core.input_database(result_db, "result_db")
    output_path = core.output_file(output, "output")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    if target_db_path is not None:
        positionals.append(target_db_path)
    positionals.append(result_db_path)
    positionals.append(output_path)
    options: list[str] = []
    options += core.opt_bool("--first-seq-as-repr", first_seq_as_repr, "first_seq_as_repr")
    options += core.opt_int("--target-column", target_column, "target_column", 0, None)
    options += core.opt_bool("--full-header", full_header, "full_header")
    options += core.opt_choice("--idx-seq-src", idx_seq_src, (0, 1, 2), "idx_seq_src")
    result = core.run(
        "mmseqs_createtsv",
        "createtsv",
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


def mmseqs_convert2fasta(
    sequence_db: str,
    output: str,
    *,
    use_header_file: bool | None = None,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a sequence database as a FASTA file (`mmseqs convert2fasta`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - output: new FASTA file.
    - use_header_file: Map entry keys through the header database instead of the sequence
      database.
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.

    The result adds ``sequence_count``.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    output_path = core.output_file(output, "output")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(output_path)
    options: list[str] = []
    options += core.opt_bool("--use-header-file", use_header_file, "use_header_file")
    result = core.run(
        "mmseqs_convert2fasta",
        "convert2fasta",
        positionals,
        options,
        databases={},
        files={"output": output_path},
        needs_tmp=False,
        threads=None,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    with output_path.open(encoding="utf-8", errors="replace") as handle:
        result["sequence_count"] = sum(1 for line in handle if line.startswith(">"))
    return result


def mmseqs_clust(
    sequence_db: str,
    result_db: str,
    cluster_db: str,
    *,
    cluster_mode: int | None = None,
    max_iterations: int | None = None,
    similarity_type: int | None = None,
    set_mode: bool | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster an alignment or prefilter result graph with set-cover, connected-component or
    greedy-by-length clustering (`mmseqs clust`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - result_db: alignment or prefilter result database.
    - cluster_db: new cluster database.
    - cluster_mode: Clustering: 0 set cover (greedy), 1 connected component (BLASTclust), 2
      or 3 greedy by sequence length (CD-HIT like) (default 0).
    - max_iterations: Maximum depth of the breadth-first search of connected-component
      clustering (default 1000).
    - similarity_type: Score used for clustering: 1 alignment score, 2 sequence identity
      (default 2).
    - set_mode: Cluster by set instead of by entry.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    result_db_path = core.input_database(result_db, "result_db")
    cluster_db_path = core.output_database(cluster_db, "cluster_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(result_db_path)
    positionals.append(cluster_db_path)
    options: list[str] = []
    options += core.opt_choice("--cluster-mode", cluster_mode, (0, 1, 2, 3), "cluster_mode")
    options += core.opt_int("--max-iterations", max_iterations, "max_iterations", 1, None)
    options += core.opt_choice("--similarity-type", similarity_type, (1, 2), "similarity_type")
    options += core.opt_bool("--set-mode", set_mode, "set_mode")
    result = core.run(
        "mmseqs_clust",
        "clust",
        positionals,
        options,
        databases={"cluster_db": cluster_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_clusthash(
    sequence_db: str,
    alignment_db: str,
    *,
    min_seq_id: float | None = None,
    max_seq_len: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster equal-length sequences by hashing them (`mmseqs clusthash`); the result is an
    alignment database of identical or near-identical pairs.

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - alignment_db: new alignment database.
    - min_seq_id: Minimum sequence identity of a reported match, 0 to 1 (default 0).
    - max_seq_len: Maximum sequence length (default 65535).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    alignment_db_path = core.output_database(alignment_db, "alignment_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(alignment_db_path)
    options: list[str] = []
    options += core.opt_float("--min-seq-id", min_seq_id, "min_seq_id", 0.0, 1.0)
    options += core.opt_int("--max-seq-len", max_seq_len, "max_seq_len", 1, None)
    result = core.run(
        "mmseqs_clusthash",
        "clusthash",
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


def mmseqs_pickconsensusrep(
    sequence_db: str,
    cluster_db: str,
    new_cluster_db: str,
    *,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Choose a new representative sequence per cluster based on the cluster consensus (`mmseqs
    pickconsensusrep`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - cluster_db: cluster database.
    - new_cluster_db: new cluster database.
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    cluster_db_path = core.input_database(cluster_db, "cluster_db")
    new_cluster_db_path = core.output_database(new_cluster_db, "new_cluster_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(cluster_db_path)
    positionals.append(new_cluster_db_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_pickconsensusrep",
        "pickconsensusrep",
        positionals,
        options,
        databases={"new_cluster_db": new_cluster_db_path},
        files={},
        needs_tmp=True,
        threads=None,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_result2repseq(
    sequence_db: str,
    result_db: str,
    representative_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the representative sequence of every cluster of a result database (`mmseqs
    result2repseq`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - sequence_db: sequence database.
    - result_db: cluster or result database.
    - representative_db: new sequence database of the representatives.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    result_db_path = core.input_database(result_db, "result_db")
    representative_db_path = core.output_database(representative_db, "representative_db")
    positionals: list[Path] = []
    positionals.append(sequence_db_path)
    positionals.append(result_db_path)
    positionals.append(representative_db_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_result2repseq",
        "result2repseq",
        positionals,
        options,
        databases={"representative_db": representative_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_summarizealis(
    alignment_db: str,
    summary_db: str,
    *,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the alignments of every query into one row with unique coverage, coverage and
    mean sequence identity (`mmseqs summarizealis`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - alignment_db: alignment database.
    - summary_db: new summary database.
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    summary_db_path = core.output_database(summary_db, "summary_db")
    positionals: list[Path] = []
    positionals.append(alignment_db_path)
    positionals.append(summary_db_path)
    options: list[str] = []
    result = core.run(
        "mmseqs_summarizealis",
        "summarizealis",
        positionals,
        options,
        databases={"summary_db": summary_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_filterresult(
    query_db: str,
    target_db: str,
    result_db: str,
    filtered_db: str,
    *,
    max_seq_id: float | None = None,
    qid: str | None = None,
    qsc: float | None = None,
    cov: float | None = None,
    diff: int | None = None,
    filter_min_enable: int | None = None,
    allow_deletion: bool | None = None,
    add_self_matches: bool | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter the pairwise alignment results of every query by redundancy, identity to the query,
    score and coverage (`mmseqs filterresult`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - result_db: alignment result database.
    - filtered_db: new filtered result database.
    - max_seq_id: Reduce redundancy: maximum pairwise sequence identity (default 0.9).
    - qid: Reduce diversity: minimum sequence identity with the query, comma-separated for
      several thresholds (0 to 1).
    - qsc: Reduce diversity: minimum score per aligned residue with the query (default -20).
    - cov: Minimum fraction of query residues covered by a kept sequence (default 0).
    - diff: Keep the most diverse sequences, at least this many per block of 50 columns.
    - filter_min_enable: Filter only MSAs with more than this many sequences; 0 always
      filters.
    - allow_deletion: Allow deletions in the MSA.
    - add_self_matches: Add each query as a match of itself (used for clustering).
    - gap_open: Gap open cost (MMseqs2 default 11 for amino acids, 5 for nucleotides).
    - gap_extend: Gap extension cost (default 1 for amino acids, 2 for nucleotides).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    result_db_path = core.input_database(result_db, "result_db")
    filtered_db_path = core.output_database(filtered_db, "filtered_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(result_db_path)
    positionals.append(filtered_db_path)
    options: list[str] = []
    options += core.opt_float("--max-seq-id", max_seq_id, "max_seq_id", 0.0, 1.0)
    options += core.opt_text("--qid", qid, "qid", "^[0-9.]+(,[0-9.]+)*$")
    options += core.opt_float("--qsc", qsc, "qsc", -50.0, 100.0)
    options += core.opt_float("--cov", cov, "cov", 0.0, 1.0)
    options += core.opt_int("--diff", diff, "diff", 0, None)
    options += core.opt_int("--filter-min-enable", filter_min_enable, "filter_min_enable", 0, None)
    options += core.opt_bool("--allow-deletion", allow_deletion, "allow_deletion")
    options += core.opt_bool("--add-self-matches", add_self_matches, "add_self_matches")
    options += core.opt_int("--gap-open", gap_open, "gap_open", 0, None)
    options += core.opt_int("--gap-extend", gap_extend, "gap_extend", 0, None)
    result = core.run(
        "mmseqs_filterresult",
        "filterresult",
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


def mmseqs_extractalignedregion(
    query_db: str,
    target_db: str,
    alignment_db: str,
    region_db: str,
    *,
    extract_mode: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut the aligned region out of the query or the target sequences (`mmseqs
    extractalignedregion`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - query_db: query sequence database.
    - target_db: target sequence database.
    - alignment_db: alignment database.
    - region_db: new sequence database of the aligned regions.
    - extract_mode: Extract the aligned region from 1 the query or 2 the target (default 2).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    query_db_path = core.input_database(query_db, "query_db")
    target_db_path = core.input_database(target_db, "target_db")
    alignment_db_path = core.input_database(alignment_db, "alignment_db")
    region_db_path = core.output_database(region_db, "region_db")
    positionals: list[Path] = []
    positionals.append(query_db_path)
    positionals.append(target_db_path)
    positionals.append(alignment_db_path)
    positionals.append(region_db_path)
    options: list[str] = []
    options += core.opt_choice("--extract-mode", extract_mode, (1, 2), "extract_mode")
    result = core.run(
        "mmseqs_extractalignedregion",
        "extractalignedregion",
        positionals,
        options,
        databases={"region_db": region_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


def mmseqs_gff2db(
    gff_files: str | Sequence[str],
    sequence_db: str,
    feature_db: str,
    *,
    gff_type: str | None = None,
    id_offset: int | None = None,
    threads: int = 1,
    mmseqs_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the regions of GFF3 features from a sequence database into a new database (`mmseqs
    gff2db`).

    Databases are given as path prefixes. An option left at None is not passed, so the
    MMseqs2 default applies. Outputs must not exist yet.

    Parameters:
    - gff_files: GFF3 file(s).
    - sequence_db: nucleotide sequence database of the genome.
    - feature_db: new sequence database of the feature sequences.
    - gff_type: Comma-separated GFF3 feature types to select (default: all).
    - id_offset: Offset of the numeric ids in the index file (default 0).
    - threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
    - mmseqs_path: path of the mmseqs executable; default is the one on PATH.
    """
    gff_files_path = core.fasta_files(gff_files, "gff_files")
    sequence_db_path = core.input_database(sequence_db, "sequence_db")
    feature_db_path = core.output_database(feature_db, "feature_db")
    positionals: list[Path] = []
    positionals.extend(gff_files_path)
    positionals.append(sequence_db_path)
    positionals.append(feature_db_path)
    options: list[str] = []
    options += core.opt_text("--gff-type", gff_type, "gff_type", "^[A-Za-z0-9_.:,-]+$")
    options += core.opt_int("--id-offset", id_offset, "id_offset", 0, None)
    result = core.run(
        "mmseqs_gff2db",
        "gff2db",
        positionals,
        options,
        databases={"feature_db": feature_db_path},
        files={},
        needs_tmp=False,
        threads=threads,
        mmseqs_path=mmseqs_path,
        verbose=verbose,
    )
    return result


__all__ = [
    "mmseqs_clust",
    "mmseqs_clusterupdate",
    "mmseqs_clusthash",
    "mmseqs_convert2fasta",
    "mmseqs_convertalis",
    "mmseqs_createdb",
    "mmseqs_createlinindex",
    "mmseqs_createtsv",
    "mmseqs_extractalignedregion",
    "mmseqs_filterresult",
    "mmseqs_gff2db",
    "mmseqs_linclust",
    "mmseqs_linsearch",
    "mmseqs_map",
    "mmseqs_maskbygff",
    "mmseqs_pickconsensusrep",
    "mmseqs_rbh",
    "mmseqs_result2msa",
    "mmseqs_result2repseq",
    "mmseqs_summarizealis",
]

_register_public_verbose(globals(), tuple(__all__))
