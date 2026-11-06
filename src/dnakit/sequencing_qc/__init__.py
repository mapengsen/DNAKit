"""Local candidate adapters for the sequencing qc capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit._candidate_runtime import run_candidate_operation as run_candidate_operation

from . import _impl_metrics as _m
from . import _impl_report as _r
from ._external import run_cutadapt as _run_cutadapt_backend
from ._external import run_fastp as _run_fastp_backend
from ._external import run_fastqc as _run_fastqc_backend
from ._external import run_multiqc as _run_multiqc_backend
from ._impl_reads import run_op as _run
from .long_reads import (
    filter_long_reads,
    filter_nanopore_reads,
    summarize_long_reads,
    trim_porechop,
)
from .short_reads import run_afterqc, trim_trimmomatic


def assign_sample_barcodes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Assign each read to a sample by matching its barcode against a sample sheet with a mismatch tolerance.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences). A single read set only: a {'read1': ..., 'read2': ...} pair mapping is not accepted (INVALID_INPUT). For paired data pass the one read that carries the sample barcode (usually the index read or read 1).
    options: samples (required; mapping sample -> barcode, dual indexes as 'i7+i5'), barcode_location ('read_start' default,
    'read_end' or 'header' = index at the end of the Illumina header comment), barcode_length (default: the sample-sheet
    barcode length), max_mismatches (default 1), min_distance_gap (default 1: the best sample must be closer than the
    second best by at least this many mismatches, otherwise the read is 'ambiguous'), max_listed (per-read rows returned,
    default 100000), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Result: per-read assignment rows, reads per sample, assigned/unassigned/ambiguous counts. Hamming distance, N counts as a mismatch.
    """
    return _run("sequencing_qc", "assign_sample_barcodes", _m.op_assign_barcodes, input, options)


def basecall_nanopore_reads(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Basecall Oxford Nanopore raw signal (FAST5/POD5) into bases with a neural-network basecaller.

    Not implemented: Basecalling converts raw FAST5/POD5 current signal to bases with a neural-network model; no basecaller or model is bundled. Always returns status "not_implemented".
    """
    return candidate_not_implemented(
        "sequencing_qc",
        "basecall_nanopore_reads",
        input,
        options=options,
        requires="Oxford Nanopore basecaller (dorado/guppy) and trained model weights",
        message="Basecalling converts raw FAST5/POD5 current signal to bases with a neural-network model; no basecaller or model is bundled.",
    )


def calculate_adapter_content(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the cumulative adapter content along reads (FastQC-style k-mer rule).

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: adapters (mapping name -> sequence or list; default FastQC's Illumina universal, small-RNA 3'/5', Nextera and
    SOLiD adapters), kmer_size (default 12), max_mismatches (Hamming, default 0), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    A read hits an adapter when the first min(kmer_size, len(adapter)) adapter bases occur in the read; the cumulative
    percentage at position p counts reads whose first hit starts at or before p. Adapters truncated by the read end are not seen.
    """
    return _run("sequencing_qc", "calculate_adapter_content", _m.op_adapter_content, input, options)


def calculate_adapter_content_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate adapter, template-switch-oligo and poly-A content of 10x Chromium read 2 (cDNA).

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Same k-mer rule and options as calculate_adapter_content (adapters, kmer_size, max_mismatches); the default list adds
    the 10x template switch oligo (AAGCAGTGGTATCAACGCAGAGTACATGGG) and a 12-mer poly-A.
    """
    return _run(
        "sequencing_qc",
        "calculate_adapter_content_tenx",
        _m.op_adapter_content_tenx,
        input,
        options,
    )


def calculate_barcode_quality(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the sequencing quality and validity of read barcodes (sample or cell barcodes).

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences). A single read set only: a {'read1': ..., 'read2': ...} pair mapping is not accepted (INVALID_INPUT). For paired data pass the one read that carries the barcode.
    options: barcode_location ('read_start' default, 'read_end' or 'header'), barcode_length (default 8; not used for 'header'),
    whitelist (list or newline text of valid barcodes), max_mismatches (default 1, used with a whitelist), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Result: distinct barcodes, N-containing barcodes, mean quality and Q30 fraction per position (when qualities exist),
    barcode entropy in bits and, with a whitelist, exact / corrected / ambiguous / invalid fractions.
    """
    return _run("sequencing_qc", "calculate_barcode_quality", _m.op_barcode_quality, input, options)


def calculate_barcode_quality_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate quality and whitelist validity of 10x cell barcodes and UMIs from read 1.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: whitelist (valid cell barcodes), max_mismatches (default 1). Result: cell-barcode report (see
    calculate_barcode_quality) plus UMI distinct count, N count, mean quality, Q30 fraction and homopolymer UMIs.
    """
    return _run(
        "sequencing_qc",
        "calculate_barcode_quality_tenx",
        _m.op_barcode_quality_tenx,
        input,
        options,
    )


def calculate_basecalling_qc(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate basecalling quality metrics: no-call rate, Q30 fraction, expected errors per read and Illumina pass-filter fraction.

    input (Phred qualities required): FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: min_mean_quality (default 20; reads with lower mean Q are counted), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    The pass-filter fraction comes from the CASAVA header flag ('1:Y:0:...' = filtered, 'N' = passing) and is null when absent.
    Instrument-side metrics (phasing, intensity) need vendor files and are not computed.
    """
    return _run("sequencing_qc", "calculate_basecalling_qc", _m.op_basecalling_qc, input, options)


def calculate_basecalling_qc_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate Cell Ranger-style basecalling metrics (Q30 in barcode, UMI and RNA read) for 10x reads.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000). Quality strings are required.
    Extra options: min_mean_quality (default 20). Each region gets the metrics of calculate_basecalling_qc.
    """
    return _run(
        "sequencing_qc", "calculate_basecalling_qc_tenx", _m.op_basecalling_qc_tenx, input, options
    )


def calculate_contamination(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Screen reads for contamination by exact k-mer containment against user-supplied reference sequences.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: references (required; mapping name -> sequence, or FASTA text), kmer_size (default 21), min_kmer_fraction (default
    0.5: share of a read's canonical k-mers present in a reference for a hit), expected_reference (optional name; reads that
    hit only other references are then reported as contamination_fraction), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    This is a k-mer screen, not an alignment, and only detects contaminants whose sequence the caller supplies.
    """
    return _run("sequencing_qc", "calculate_contamination", _m.op_contamination, input, options)


def calculate_contamination_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Screen the cDNA reads (read 2) of 10x data for contamination by k-mer containment against supplied references.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options as calculate_contamination: references (required), kmer_size, min_kmer_fraction, expected_reference.
    """
    return _run(
        "sequencing_qc", "calculate_contamination_tenx", _m.op_contamination_tenx, input, options
    )


def calculate_duplicate_rate(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the exact-sequence duplicate rate of reads: (reads - distinct sequences) / reads.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: prefix_length (compare only the first N bases; default 0 = whole read), reverse_complement_aware (default False),
    phred_offset (33 or 64, default 33), max_reads (default 1000000).
    For pairs the key is the concatenation of both mates. Positions are unknown, so this is sequence duplication, not mapping duplication.
    """
    return _run("sequencing_qc", "calculate_duplicate_rate", _m.op_duplicate_rate, input, options)


def calculate_duplication(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate the FastQC-style sequence duplication level histogram of reads.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: prefix_length (default 0), reverse_complement_aware (default False), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Levels: 1-9 copies, >10, >50, >100, >500, >1k, >5k, >10k+; reported as percent of total reads and of distinct sequences.
    No FastQC-style extrapolation to unseen reads is applied.
    """
    return _run("sequencing_qc", "calculate_duplication", _m.op_duplication, input, options)


def calculate_duplication_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the duplication of 10x reads, keyed on cell barcode + UMI + the first bases of read 2.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: cdna_prefix_length (default 25). Heuristic proxy for sequencing saturation: Cell Ranger uses
    (barcode, UMI, gene) from an alignment, which is unavailable here.
    """
    return _run(
        "sequencing_qc", "calculate_duplication_tenx", _m.op_duplication_tenx, input, options
    )


def calculate_error_profile(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the per-position expected error profile of reads from their Phred scores.

    input (Phred qualities required): FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Expected error rate at a position = mean of 10^(-Q/10); the result also gives cumulative expected errors, expected
    errors per read and the N fraction per position. This is the error the quality scores claim, not an observed mismatch
    rate (use estimate_read_error_model with a reference for that).
    """
    return _run("sequencing_qc", "calculate_error_profile", _m.op_error_profile, input, options)


def calculate_error_profile_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate Phred-derived expected error profiles for the barcode, UMI and RNA regions of 10x reads.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000). Quality strings are required. Each region gets the result of calculate_error_profile.
    """
    return _run(
        "sequencing_qc", "calculate_error_profile_tenx", _m.op_error_profile_tenx, input, options
    )


def calculate_gc_bias(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate the GC content distribution of reads and its deviation from a fitted normal distribution (FastQC convention).

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: reference_gc (optional expected GC fraction in [0, 1]; the mean shift is reported), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Result: per-read GC% histogram, mean/sd/mode, theoretical normal counts, deviation percentage and GC% by position.
    This is a read-composition profile, not Picard's coverage-versus-GC bias, which needs alignments.
    """
    return _run("sequencing_qc", "calculate_gc_bias", _m.op_gc_bias, input, options)


def calculate_gc_bias_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate the GC content distribution of 10x cDNA reads (read 2).

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: reference_gc (optional GC fraction in [0, 1]). Same result as calculate_gc_bias.
    """
    return _run("sequencing_qc", "calculate_gc_bias_tenx", _m.op_gc_bias_tenx, input, options)


def calculate_index_hopping(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the index-hopping rate of dual-indexed reads from the sample indexes in their headers.

    input: one read set (FASTQ text or path, or a read list; not a {'read1': ..., 'read2': ...} pair, which raises INVALID_INPUT; for paired data pass one mate, both mates carry the same indexes) whose Illumina header comments end in 'i7+i5' (for example '1:N:0:ACGTAC+TTGGCA').
    options: samples (required; mapping sample -> 'i7+i5'), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    A read is hopped when its i7 and i5 both occur in the sample sheet but not as a pair; rate = hopped / (hopped + expected
    pair reads). Indexes must match exactly. Unknown-index reads are reported separately.
    """
    return _run("sequencing_qc", "calculate_index_hopping", _m.op_index_hopping, input, options)


def calculate_index_hopping_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the index-hopping rate of dual-indexed 10x libraries from the sample indexes in read headers.

    input: one read set (FASTQ text, path or read list; a {'read1': ..., 'read2': ...} pair is not accepted and raises INVALID_INPUT, pass one mate) whose header comment ends in 'i7+i5'; options as calculate_index_hopping (samples required).
    Applies to 10x dual-index kits; single-index kits (one index per sample) carry no hopping signal and raise an error.
    """
    return _run(
        "sequencing_qc", "calculate_index_hopping_tenx", _m.op_index_hopping, input, options
    )


def calculate_insert_size(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate insert (fragment) size statistics from read-pair overlap, SAM text, a 10x fragments table or a list of sizes.

    input: {'read1': ..., 'read2': ...} (size = overlap of read 1 with the reverse complement of read 2), SAM text/path
    (|TLEN| of mapped first-in-pair reads), a fragments/BED-like TSV (chrom, start, end[, barcode, count]; size = end - start,
    weighted by count) or a list of numbers.
    options: min_overlap (default 10), max_mismatch_rate (default 0.1), max_pairs (pairs sized by overlap, default 200000),
    proper_pairs_only (SAM, default True), max_insert_size (default 100000), phred_offset, max_reads. Pairs without a detectable overlap are excluded and counted.
    """
    return _run("sequencing_qc", "calculate_insert_size", _m.op_insert_size, input, options)


def calculate_insert_size_distribution(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the insert (fragment) size histogram with cumulative fractions.

    input and options as calculate_insert_size, plus bin_width (default 1, widened automatically for wide ranges).
    Result: the summary statistics of calculate_insert_size and a histogram of {start, end, count, fraction, cumulative_fraction}.
    """
    return _run(
        "sequencing_qc",
        "calculate_insert_size_distribution",
        _m.op_insert_size_distribution,
        input,
        options,
    )


def calculate_insert_size_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the fragment size distribution of 10x data from a fragments.tsv table, SAM text or a list of sizes.

    input: 10x scATAC fragments.tsv text/path (chrom, start, end, barcode, count; size = end - start weighted by count), SAM
    text, a list of numbers or a read pair. Options and result as calculate_insert_size. 10x scRNA read 1 holds the barcode
    and UMI, so scRNA libraries have no read-pair insert size.
    """
    return _run("sequencing_qc", "calculate_insert_size_tenx", _m.op_insert_size, input, options)


def calculate_phasing_qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate Illumina phasing and pre-phasing percentages from the sequencer's InterOp run metrics.

    Not implemented: Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads. Always returns status "not_implemented".
    """
    return candidate_not_implemented(
        "sequencing_qc",
        "calculate_phasing_qc",
        input,
        options=options,
        requires="Illumina InterOp run-folder files (phasing/prephasing metrics)",
        message="Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads.",
    )


def calculate_phasing_qc_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate Illumina phasing and pre-phasing percentages of a 10x library run from the sequencer's InterOp run metrics.

    Not implemented: Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads. Always returns status "not_implemented".
    """
    return candidate_not_implemented(
        "sequencing_qc",
        "calculate_phasing_qc_tenx",
        input,
        options=options,
        requires="Illumina InterOp run-folder files (phasing/prephasing metrics)",
        message="Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads.",
    )


def calculate_qscore_summary(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate a global Q-score summary: mean/median Q, fractions of bases at Q10-Q40 and read-level Q thresholds.

    input (Phred qualities required): FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: read_q_method ('probability' default = -10*log10 of the mean error probability, as MinKNOW/NanoPlot; or 'arithmetic'),
    phred_offset (33 or 64, default 33), max_reads (default 1000000).
    """
    return _run("sequencing_qc", "calculate_qscore_summary", _m.op_qscore_summary, input, options)


def calculate_quality(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate per-position base quality statistics (mean, median, quartiles, 10th/90th percentile) and the per-read mean quality distribution.

    input (Phred qualities required): FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Percentiles use linear interpolation between order statistics and are rounded to 6 decimals. Per-read mean quality is the arithmetic mean of Q (FastQC).
    """
    return _run("sequencing_qc", "calculate_quality", _m.op_quality, input, options)


def calculate_quality_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate per-position quality statistics for the barcode, UMI and RNA regions of 10x reads.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000). Quality strings are required. Each region gets the result of calculate_quality.
    """
    return _run("sequencing_qc", "calculate_quality_tenx", _m.op_quality_tenx, input, options)


def calculate_read_length(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate the read length distribution: min/max/mean/median/sd, N50, N90, L50 and a histogram.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: bin_width (default 1, automatic for wide ranges), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    N50 is the length L such that reads of length >= L contain at least half of all bases.
    """
    return _run("sequencing_qc", "calculate_read_length", _m.op_read_length, input, options)


def calculate_read_length_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Calculate the read length distribution of read 1 and read 2 of 10x data.

    input: a pair {'read1': ..., 'read2': ...}. options: bin_width, phred_offset, max_reads. Result as calculate_read_length per mate.
    """
    return _run(
        "sequencing_qc", "calculate_read_length_tenx", _m.op_read_length_tenx, input, options
    )


def calculate_yield(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate sequencing yield: reads, bases, mean length, N bases, GC fraction and Q20/Q30 bases.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: genome_size (optional; adds expected_coverage = bases / genome_size), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    """
    return _run("sequencing_qc", "calculate_yield", _m.op_yield, input, options)


def calculate_yield_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate the yield of 10x reads per region (barcode, UMI, RNA read) plus barcode counts.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: whitelist (valid cell barcodes), max_mismatches (default 1). Result adds distinct barcodes, median reads
    per barcode and, with a whitelist, the valid-barcode fraction.
    """
    return _run("sequencing_qc", "calculate_yield_tenx", _m.op_yield_tenx, input, options)


def collapse_umi_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Collapse reads that share a UMI (and read start) into one consensus read per molecule.

    input: one read set whose names end in '<umi_separator><UMI>' (umi_tools convention), e.g. 'read1_ACGTACGT'. A {'read1': ..., 'read2': ...} pair is not accepted (INVALID_INPUT); for paired data collapse the read that carries the UMI name.
    options: umi_separator ('_'), group_by ('umi_and_start' default or 'umi'), start_length (default 10), correct_errors
    (default False; directional UMI network first), threshold (Hamming, default 1), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Without alignments the read start approximates the mapping position. The consensus is a quality-weighted vote per
    position over positions covered by at least half of the group. Consensus quality = clip(sum Q agreeing - sum Q disagreeing, 0, 93).
    """
    return _run("sequencing_qc", "collapse_umi_reads", _m.op_collapse_umi_reads, input, options)


def correct_long_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Correct sequencing errors in long reads by overlap-based consensus (canu/ratatosk style).

    Not implemented: Long-read self-correction needs all-vs-all overlap detection and consensus calling that DNAKit does not implement natively. Always returns status "not_implemented".
    """
    return candidate_not_implemented(
        "sequencing_qc",
        "correct_long_reads",
        input,
        options=options,
        requires="long-read error-correction program (canu, ratatosk or similar)",
        message="Long-read self-correction needs all-vs-all overlap detection and consensus calling that DNAKit does not implement natively.",
    )


def correct_umi_errors(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Correct UMI sequencing errors by clustering UMIs on a Hamming-distance network (UMI-tools methods).

    input: a list of UMI strings, a mapping UMI -> read count, or reads whose names end in '_<UMI>' (at most 5000 distinct UMIs).
    options: method ('directional' default, 'adjacency' or 'cluster'), threshold (Hamming, default 1), umi_separator ('_'), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    directional: edge a->b if count(a) >= 2*count(b)-1; groups grown by breadth-first search from the most abundant UMI.
    The most abundant UMI of each group is its representative.
    """
    return _run("sequencing_qc", "correct_umi_errors", _m.op_correct_umi_errors, input, options)


def demultiplex_barcoded_reads(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Demultiplex reads into per-sample read sets by barcode, trimming the barcode when it is part of the read.

    input and options as assign_sample_barcodes: one read set only, a {'read1': ..., 'read2': ...} pair is not accepted (INVALID_INPUT); for paired data demultiplex the read that carries the barcode (samples required, barcode_location, barcode_length, max_mismatches,
    min_distance_gap), plus trim_barcode (default True for 'read_start'/'read_end'), output_directory (optional existing
    directory: one <sample>.fastq per sample is written there, never overwriting) and max_listed (reads returned per sample, default 1000).
    Reads without an unambiguous match go to 'unassigned'.
    """
    return _run("sequencing_qc", "demultiplex_barcoded_reads", _m.op_demultiplex, input, options)


def detect_cross_sample_contamination(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Flag sample pairs that share an unexpectedly large part of their k-mers (possible cross-sample contamination).

    input: mapping sample name -> reads (FASTQ text/path or read list); at least two samples.
    options: kmer_size (default 21), threshold (containment in [0, 1], default 0.1), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    containment[a][b] = |K(a) intersect K(b)| / |K(a)| over canonical k-mers. Heuristic: only meaningful when the samples are
    expected to be unrelated (different organisms/individuals); related samples share k-mers legitimately.
    """
    return _run(
        "sequencing_qc",
        "detect_cross_sample_contamination",
        _m.op_cross_sample_contamination,
        input,
        options,
    )


def detect_platform(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect the sequencing platform from the read-name format (rule-based).

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Recognises Illumina (instrument:run:flowcell:lane:tile:x:y), Oxford Nanopore (runid=/ch=/UUID), PacBio (movie/zmw/ccs),
    MGI/BGI (L#C#R#) and Ion Torrent names; SRA-renamed reads carry no platform information. The platform is the majority
    layout when it covers at least half of the reads. The instrument model is not inferred.
    """
    return _run("sequencing_qc", "detect_platform", _m.op_detect_platform, input, options)


def detect_platform_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect the sequencing platform of 10x reads from read-2 names and infer the Chromium barcode+UMI layout from read-1 length.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Only the read-1 length is used (24 / 26 / 28 maps to 14+10 / 16+10 / 16+12 bp layouts and the compatible Chromium chemistries); read-1 sequences are not compared with a chemistry or whitelist. This is a
    length heuristic; Cell Ranger's whitelist-based chemistry autodetection is not performed.
    """
    return _run("sequencing_qc", "detect_platform_tenx", _m.op_detect_platform_tenx, input, options)


def estimate_library_complexity(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Estimate library complexity (distinct molecules) with the Lander-Waterman / Picard model.

    input: reads (exact-sequence duplicates define the distinct set; pairs use both mates) or a mapping
    {'total_reads': N, 'unique_reads': C} (also accepted in options).
    options: prefix_length, reverse_complement_aware, phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Solves C / X = 1 - exp(-N / X) for the library size X and predicts the unique reads at 0.5x, 1x, 2x, 5x and 10x depth.
    X is unbounded (null) when every read is distinct.
    """
    return _run(
        "sequencing_qc", "estimate_library_complexity", _m.op_library_complexity, input, options
    )


def estimate_read_error_model(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Estimate a read error model: observed mismatch rates against a reference, or the Phred-derived expected error profile.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences). A single read set only: a {'read1': ..., 'read2': ...} pair mapping is not accepted (INVALID_INPUT); for paired data run the function once per mate.
    options: reference (sequence or FASTA text up to 5 Mb; without it the model is the quality-derived profile), seed_length
    (default 11), max_mismatch_fraction (default 0.2), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    With a reference each read is placed by k-mer seed voting (both strands) and scored by gap-free Hamming distance. Result:
    mismatch rate overall and per read position, substitution counts and empirical error rate per reported Q. Indels are not modelled.
    """
    return _run(
        "sequencing_qc", "estimate_read_error_model", _m.op_read_error_model, input, options
    )


def extract_umi_sequences(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Extract UMI sequences from reads and append them to the read names (umi_tools convention).

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences). A single read set only: a {'read1': ..., 'read2': ...} pair mapping is not accepted (INVALID_INPUT). For paired data pass the one read that carries the UMI (usually read 1).
    options: umi_length (default 8) with umi_location ('5prime' default or '3prime'), or pattern (5' only; N = UMI base,
    C = cell-barcode base, X = discarded base, e.g. 'CCCCCCNNNNXX'), umi_separator ('_'), append_to_name (default True), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Result: records with id, umi, cell_barcode, trimmed sequence and quality.
    """
    return _run("sequencing_qc", "extract_umi_sequences", _m.op_extract_umi, input, options)


def generate_qc_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Generate a FastQC-style QC report: yield, length, quality, GC, adapter, duplication and platform sections with pass/warn/fail flags.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: adapters, output_path (optional .txt/.md file; refuses to overwrite unless overwrite=True), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Flags use the FastQC thresholds listed in the result (flag_definitions). The result also contains report_text.
    """
    return _run("sequencing_qc", "generate_qc_report", _r.op_report, input, options)


def generate_qc_report_tenx(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Generate a QC report for 10x Chromium reads: yield and basecalling per region, barcode/UMI validity, duplication proxy and FastQC-style flags for the RNA read.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: whitelist, adapters, output_path (optional .txt/.md; no overwrite unless overwrite=True).
    """
    return _run("sequencing_qc", "generate_qc_report_tenx", _r.op_report_tenx, input, options)


def plot_per_base_quality(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a per-base quality box plot (FastQC style) as an SVG string.

    input (Phred qualities required): FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: mate ('read1' default or 'read2' for pairs), max_columns (default 150; longer reads are grouped), output_path
    (optional .svg; no overwrite unless overwrite=True), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    Boxes show lower to upper quartile, the red line the median, whiskers the 10th-90th percentile and the blue line the mean.
    """
    return _run(
        "sequencing_qc", "plot_per_base_quality", _r.op_plot_per_base_quality, input, options
    )


def plot_qc_summary(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a bar chart of key QC percentages (Q20, Q30, GC, N, duplication, adapter) as an SVG string.

    input: FASTQ text or file path (.gz allowed), a list of reads ((id, sequence, quality) tuples, {'sequence', 'quality'} mappings, FASTQ record objects or bare sequences), or a pair {'read1': ..., 'read2': ...}.
    options: adapters, output_path (optional .svg; no overwrite unless overwrite=True), phred_offset (33 or 64, default 33), max_reads (default 1000000).
    """
    return _run("sequencing_qc", "plot_qc_summary", _r.op_plot_qc_summary, input, options)


def plot_qc_summary_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a bar chart of 10x QC percentages (Q30 in barcode/UMI/RNA read, adapter, duplication proxy, GC) as an SVG string.

    input: a pair {'read1': <barcode+UMI reads>, 'read2': <cDNA reads>} (FASTQ text, path or read list). options: chemistry ('v1' 14+10, 'v2' 16+10, 'v3' 16+12 bp barcode+UMI; default 'v3'), barcode_length, umi_length, phred_offset (33), max_reads (1000000).
    Extra options: whitelist (adds the valid-barcode bar), adapters, output_path (optional .svg).
    """
    return _run("sequencing_qc", "plot_qc_summary_tenx", _r.op_plot_qc_summary_tenx, input, options)


def run_cutadapt(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_cutadapt_backend(input, options=options, verbose=verbose)


def run_fastp(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_fastp_backend(input, options=options, verbose=verbose)


def run_fastqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_fastqc_backend(input, options=options, verbose=verbose)


def run_multiqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_multiqc_backend(input, options=options, verbose=verbose)


__all__ = [
    "assign_sample_barcodes",
    "basecall_nanopore_reads",
    "calculate_adapter_content",
    "calculate_adapter_content_tenx",
    "calculate_barcode_quality",
    "calculate_barcode_quality_tenx",
    "calculate_basecalling_qc",
    "calculate_basecalling_qc_tenx",
    "calculate_contamination",
    "calculate_contamination_tenx",
    "calculate_duplicate_rate",
    "calculate_duplication",
    "calculate_duplication_tenx",
    "calculate_error_profile",
    "calculate_error_profile_tenx",
    "calculate_gc_bias",
    "calculate_gc_bias_tenx",
    "calculate_index_hopping",
    "calculate_index_hopping_tenx",
    "calculate_insert_size",
    "calculate_insert_size_distribution",
    "calculate_insert_size_tenx",
    "calculate_phasing_qc",
    "calculate_phasing_qc_tenx",
    "calculate_qscore_summary",
    "calculate_quality",
    "calculate_quality_tenx",
    "calculate_read_length",
    "calculate_read_length_tenx",
    "calculate_yield",
    "calculate_yield_tenx",
    "collapse_umi_reads",
    "correct_long_reads",
    "correct_umi_errors",
    "demultiplex_barcoded_reads",
    "detect_cross_sample_contamination",
    "detect_platform",
    "detect_platform_tenx",
    "estimate_library_complexity",
    "estimate_read_error_model",
    "extract_umi_sequences",
    "filter_long_reads",
    "filter_nanopore_reads",
    "generate_qc_report",
    "generate_qc_report_tenx",
    "plot_per_base_quality",
    "plot_qc_summary",
    "plot_qc_summary_tenx",
    "run_afterqc",
    "run_cutadapt",
    "run_fastp",
    "run_fastqc",
    "run_multiqc",
    "summarize_long_reads",
    "trim_porechop",
    "trim_trimmomatic",
]

__candidate_adapter_operations__ = set(
    [
        "assign_sample_barcodes",
        "basecall_nanopore_reads",
        "calculate_adapter_content",
        "calculate_adapter_content_tenx",
        "calculate_barcode_quality",
        "calculate_barcode_quality_tenx",
        "calculate_basecalling_qc",
        "calculate_basecalling_qc_tenx",
        "calculate_contamination",
        "calculate_contamination_tenx",
        "calculate_duplicate_rate",
        "calculate_duplication",
        "calculate_duplication_tenx",
        "calculate_error_profile",
        "calculate_error_profile_tenx",
        "calculate_gc_bias",
        "calculate_gc_bias_tenx",
        "calculate_index_hopping",
        "calculate_index_hopping_tenx",
        "calculate_insert_size",
        "calculate_insert_size_distribution",
        "calculate_insert_size_tenx",
        "calculate_phasing_qc",
        "calculate_phasing_qc_tenx",
        "calculate_qscore_summary",
        "calculate_quality",
        "calculate_quality_tenx",
        "calculate_read_length",
        "calculate_read_length_tenx",
        "calculate_yield",
        "calculate_yield_tenx",
        "collapse_umi_reads",
        "correct_long_reads",
        "correct_umi_errors",
        "demultiplex_barcoded_reads",
        "detect_cross_sample_contamination",
        "detect_platform",
        "detect_platform_tenx",
        "estimate_library_complexity",
        "estimate_read_error_model",
        "extract_umi_sequences",
        "generate_qc_report",
        "generate_qc_report_tenx",
        "plot_per_base_quality",
        "plot_qc_summary",
        "plot_qc_summary_tenx",
        "run_cutadapt",
        "run_fastp",
        "run_fastqc",
        "run_multiqc",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

# Optional-backend exports (see docs/api/features)
from .backends import (
    calculate_expected_errors,
    cluster_umi_counts,
    extract_read_barcodes,
    find_adapter_match,
    find_quality_trim_bounds,
    map_barcodes_to_whitelist,
    read_slow5_signal,
    summarize_slow5_signals,
    trim_fastq_adapters,
    write_slow5_signals,
)

__all__ = [
    *__all__,
    "calculate_expected_errors",
    "cluster_umi_counts",
    "extract_read_barcodes",
    "find_adapter_match",
    "find_quality_trim_bounds",
    "map_barcodes_to_whitelist",
    "read_slow5_signal",
    "summarize_slow5_signals",
    "trim_fastq_adapters",
    "write_slow5_signals",
]

# Full-pass exports (see docs/api/features)
from .cutadapt_functions import (
    adapter_kmer_positions,
    adapter_kmer_table,
    edit_neighbours,
    expand_adapter_braces,
    hamming_neighbours,
    nextseq_trim_bounds,
    parse_adapter_spec,
    poly_a_tail_start,
    read_adapter_fasta,
    sequence_kmer_chunks,
)
from .umi_tools_functions import (
    average_umi_distance,
    barcode_indel_matches,
    call_cell_whitelist,
    cigar_splice_offset,
    count_umis_per_gene,
    detect_bam_tags,
    extract_barcode_positions,
    find_umi_neighbours,
    mask_low_quality_umi,
    parse_umi_read_name,
    read_barcode_whitelist,
    tag_read_name_with_umi,
    umi_connected_component,
)

__all__ = [
    *__all__,
    "adapter_kmer_positions",
    "adapter_kmer_table",
    "edit_neighbours",
    "expand_adapter_braces",
    "hamming_neighbours",
    "nextseq_trim_bounds",
    "parse_adapter_spec",
    "poly_a_tail_start",
    "read_adapter_fasta",
    "sequence_kmer_chunks",
    "average_umi_distance",
    "barcode_indel_matches",
    "call_cell_whitelist",
    "cigar_splice_offset",
    "count_umis_per_gene",
    "detect_bam_tags",
    "extract_barcode_positions",
    "find_umi_neighbours",
    "mask_low_quality_umi",
    "parse_umi_read_name",
    "read_barcode_whitelist",
    "tag_read_name_with_umi",
    "umi_connected_component",
]

# Review-pass exports (see docs/api/features)
from .cutadapt_review_functions import (
    adapter_untrimmed_interval,
)
from .umi_tools_review_functions import (
    barcode_from_bam_tag,
    metacontig_contig_map,
    parse_prefixed_umi_name,
    resolve_above_knee_barcodes,
    soft_clip_adjusted_position,
    tag_reads_by_metacontig,
)

__all__ = [
    *__all__,
    "adapter_untrimmed_interval",
    "barcode_from_bam_tag",
    "metacontig_contig_map",
    "parse_prefixed_umi_name",
    "resolve_above_knee_barcodes",
    "soft_clip_adjusted_position",
    "tag_reads_by_metacontig",
]


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .umi_tools_gap_wrappers import (
    umi_count_per_gene,
    umi_count_table,
    umi_dedup_bam,
    umi_extract_to_names,
    umi_group_bam,
    umi_prepare_bam_for_rsem,
    umi_whitelist_cells,
)

__all__ = [
    *__all__,
    "umi_count_per_gene",
    "umi_count_table",
    "umi_dedup_bam",
    "umi_extract_to_names",
    "umi_group_bam",
    "umi_prepare_bam_for_rsem",
    "umi_whitelist_cells",
]

# Public integration of the reviewed P00–P45 operation packs.
from .afterqc_p15_cli_functions import (
    afterqc_debubble,
    afterqc_move_barcodes_to_read_names,
    afterqc_overlap_analysis,
)
from .fastp_p14_cli_functions import (
    fastp_process_reads,
)
from .gatk_p00_metrics_cli_functions import (
    gatk_calculate_fingerprint_metrics,
    gatk_calculate_read_group_checksum,
    gatk_check_duplicate_marking,
    gatk_check_fingerprint,
    gatk_check_pileup,
    gatk_check_terminator_block,
    gatk_collect_alignment_summary_metrics,
    gatk_collect_base_distribution_by_cycle,
    gatk_collect_hiseqx_pf_fail_metrics,
    gatk_collect_hs_metrics,
    gatk_collect_illumina_basecalling_metrics,
    gatk_collect_illumina_lane_metrics,
    gatk_collect_independent_replicate_metrics,
    gatk_collect_jumping_library_metrics,
    gatk_collect_multiple_metrics,
    gatk_collect_oxog_metrics,
    gatk_collect_quality_yield_metrics,
    gatk_collect_rnaseq_metrics,
    gatk_collect_rrbs_metrics,
    gatk_collect_sequencing_artifact_metrics,
    gatk_collect_targeted_pcr_metrics,
    gatk_collect_umi_prevalence_metrics,
    gatk_collect_variant_calling_metrics,
    gatk_collect_wgs_metrics,
    gatk_compare_base_qualities,
    gatk_compare_duplicates_spark,
    gatk_compare_sams,
    gatk_convert_haplotype_database_to_vcf,
    gatk_crosscheck_fingerprints,
    gatk_dump_tabix_index,
    gatk_extract_fingerprint,
    gatk_extract_illumina_barcodes,
    gatk_get_sample_name,
    gatk_htsget_download,
    gatk_identify_contaminant,
    gatk_illumina_basecalls_to_reads,
    gatk_lift_over_haplotype_map,
    gatk_mark_illumina_adapters,
    gatk_mean_quality_by_cycle,
    gatk_quality_score_distribution,
    gatk_split_cram,
    gatk_validate_sam_file,
    gatk_view_sam,
)
from .nanoplot_p12_cli_functions import (
    nanoplot_plots,
)
from .p16_readqc_cli_functions import (
    demultiplex_nanopore_barcodes,
    filter_long_reads_with_reference,
    filter_nanopore_reads_by_gc_summary,
    filter_nanopore_reads_by_start_time,
    trim_nanopore_adapters_discard_middle,
)
from .picard_metrics_p08_cli_functions import (
    picard_accumulate_metrics,
    picard_arrays_variant_calling_metrics,
    picard_base_distribution_by_cycle,
    picard_check_fingerprint,
    picard_check_illumina_directory,
    picard_check_terminator_block,
    picard_cluster_crosscheck_metrics,
    picard_convert_artifact_to_oxog,
    picard_crosscheck_fingerprints,
    picard_duplicate_metrics,
    picard_extract_fingerprint,
    picard_fingerprint_metrics,
    picard_format_contamination_metrics,
    picard_oxog_metrics,
    picard_read_group_checksum,
    picard_rnaseq_metrics,
    picard_sam_error_metrics,
    picard_sequencing_artifact_metrics,
    picard_snvq_yield_metrics,
    picard_variant_calling_metrics,
    picard_wgs_metrics_nonzero_coverage,
)
from .pyslow5_p41_functions import (
    slow5_read_records,
)
from .trimmomatic_p15_cli_functions import (
    trimmomatic_pair_reads,
    trimmomatic_process_reads,
)

__all__ = [
    *__all__,
    "afterqc_debubble",
    "afterqc_move_barcodes_to_read_names",
    "afterqc_overlap_analysis",
    "demultiplex_nanopore_barcodes",
    "fastp_process_reads",
    "filter_long_reads_with_reference",
    "filter_nanopore_reads_by_gc_summary",
    "filter_nanopore_reads_by_start_time",
    "gatk_calculate_fingerprint_metrics",
    "gatk_calculate_read_group_checksum",
    "gatk_check_duplicate_marking",
    "gatk_check_fingerprint",
    "gatk_check_pileup",
    "gatk_check_terminator_block",
    "gatk_collect_alignment_summary_metrics",
    "gatk_collect_base_distribution_by_cycle",
    "gatk_collect_hiseqx_pf_fail_metrics",
    "gatk_collect_hs_metrics",
    "gatk_collect_illumina_basecalling_metrics",
    "gatk_collect_illumina_lane_metrics",
    "gatk_collect_independent_replicate_metrics",
    "gatk_collect_jumping_library_metrics",
    "gatk_collect_multiple_metrics",
    "gatk_collect_oxog_metrics",
    "gatk_collect_quality_yield_metrics",
    "gatk_collect_rnaseq_metrics",
    "gatk_collect_rrbs_metrics",
    "gatk_collect_sequencing_artifact_metrics",
    "gatk_collect_targeted_pcr_metrics",
    "gatk_collect_umi_prevalence_metrics",
    "gatk_collect_variant_calling_metrics",
    "gatk_collect_wgs_metrics",
    "gatk_compare_base_qualities",
    "gatk_compare_duplicates_spark",
    "gatk_compare_sams",
    "gatk_convert_haplotype_database_to_vcf",
    "gatk_crosscheck_fingerprints",
    "gatk_dump_tabix_index",
    "gatk_extract_fingerprint",
    "gatk_extract_illumina_barcodes",
    "gatk_get_sample_name",
    "gatk_htsget_download",
    "gatk_identify_contaminant",
    "gatk_illumina_basecalls_to_reads",
    "gatk_lift_over_haplotype_map",
    "gatk_mark_illumina_adapters",
    "gatk_mean_quality_by_cycle",
    "gatk_quality_score_distribution",
    "gatk_split_cram",
    "gatk_validate_sam_file",
    "gatk_view_sam",
    "nanoplot_plots",
    "picard_accumulate_metrics",
    "picard_arrays_variant_calling_metrics",
    "picard_base_distribution_by_cycle",
    "picard_check_fingerprint",
    "picard_check_illumina_directory",
    "picard_check_terminator_block",
    "picard_cluster_crosscheck_metrics",
    "picard_convert_artifact_to_oxog",
    "picard_crosscheck_fingerprints",
    "picard_duplicate_metrics",
    "picard_extract_fingerprint",
    "picard_fingerprint_metrics",
    "picard_format_contamination_metrics",
    "picard_oxog_metrics",
    "picard_read_group_checksum",
    "picard_rnaseq_metrics",
    "picard_sam_error_metrics",
    "picard_sequencing_artifact_metrics",
    "picard_snvq_yield_metrics",
    "picard_variant_calling_metrics",
    "picard_wgs_metrics_nonzero_coverage",
    "slow5_read_records",
    "trim_nanopore_adapters_discard_middle",
    "trimmomatic_pair_reads",
    "trimmomatic_process_reads",
]

_register_public_verbose(globals())
