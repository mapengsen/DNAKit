"""GATK reference, FASTA, interval and feature-file command-line facades (pack P01).

Each public function runs one GATK 4.7.0.0 tool (or one family of tools selected by a mode argument) through
``dnakit.backends.external.execute_bounded_command`` with a fixed option whitelist and no shell. The user installs GATK
(conda ``gatk4`` or the release zip) with Java 17 and either puts the ``gatk`` launcher on PATH, sets ``DNAKIT_GATK`` or
passes ``gatk_path``. Outputs are written to explicit new paths that must not exist; nothing is overwritten. Results are
JSON-serialisable dictionaries with status, backend, tool, version, command, outputs, a bounded log tail and provenance.
15 public functions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gatk_p01_core as _core
from dnakit.bioinformatics._gatk_p01_core import Param as _Param
from dnakit.bioinformatics._gatk_p01_core import Spec as _Spec

_SPECS: dict[str, _Spec] = {
    'gatk_index_feature_file': _Spec(
        'gatk_index_feature_file',
        {'default': 'IndexFeatureFile'},
        (
            _Param('feature_path', 'in', 'input', req=True),
            _Param('index_path', 'out', 'output'),
        ),
        derived=(('index', 'feature_path', ('.idx', '.tbi', '.csi'), False),),
    ),
    'gatk_sort_gff': _Spec(
        'gatk_sort_gff',
        {'default': 'SortGff'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('sequence_dictionary', 'in', 'SEQUENCE_DICTIONARY'),
        ),
    ),
    'gatk_gtf_to_bed': _Spec(
        'gatk_gtf_to_bed',
        {'default': 'GtfToBed'},
        (
            _Param('gtf_path', 'in', 'gtf-path', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('sequence_dictionary', 'in', 'sequence-dictionary', req=True),
            _Param('sort_by_transcript', 'bool', 'sort-by-transcript'),
            _Param('use_basic_transcript', 'bool', 'use-basic-transcript'),
        ),
    ),
    'gatk_bait_designer': _Spec(
        'gatk_bait_designer',
        {'default': 'BaitDesigner'},
        (
            _Param('targets', 'in', 'TARGETS', req=True),
            _Param('reference_path', 'in', 'REFERENCE_SEQUENCE', req=True),
            _Param('design_name', 'str', 'DESIGN_NAME', req=True),
            _Param('output_dir', 'outd', 'OUTPUT_DIRECTORY', req=True),
            _Param('bait_size', 'int', 'BAIT_SIZE', minimum=1),
            _Param('bait_offset', 'int', 'BAIT_OFFSET', minimum=1),
            _Param('design_strategy', 'choice', 'DESIGN_STRATEGY', choices=('CenteredConstrained', 'Fixed', 'Flanking')),
        ),
    ),
    'gatk_bwa_mem_index_image': _Spec(
        'gatk_bwa_mem_index_image',
        {'default': 'BwaMemIndexImageCreator'},
        (
            _Param('reference_path', 'in', 'input', req=True),
        ),
        derived=(('image', 'reference_path', ('.img',), True),),
    ),
    'gatk_check_reference_compatibility': _Spec(
        'gatk_check_reference_compatibility',
        {'default': 'CheckReferenceCompatibility'},
        (
            _Param('references_to_compare', 'ins', 'references-to-compare', req=True),
            _Param('alignment_path', 'in', 'input'),
            _Param('variant_path', 'in', 'variant'),
            _Param('output_path', 'out', 'output'),
        ),
        any_of=((None, ('alignment_path', 'variant_path')),),
    ),
    'gatk_compare_references': _Spec(
        'gatk_compare_references',
        {'default': 'CompareReferences'},
        (
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('references_to_compare', 'ins', 'references-to-compare', req=True),
            _Param('display_only_differing_sequences', 'bool', 'display-only-differing-sequences'),
        ),
        hook='tab_table',
    ),
    'gatk_compose_str_table': _Spec(
        'gatk_compose_str_table',
        {'default': 'ComposeSTRTableFile'},
        (
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('max_period', 'int', 'max-period', minimum=1),
            _Param('max_repeats', 'int', 'max-repeats', minimum=1),
            _Param('decimation', 'str', 'decimation'),
        ),
    ),
    'gatk_count_reference_bases': _Spec(
        'gatk_count_reference_bases',
        {'count_bases': 'CountBasesInReference', 'non_n_size': 'NonNFastaSize'},
        (
            _Param('fasta_path', 'in', {'count_bases': 'reference', 'non_n_size': 'INPUT'}, req=True),
            _Param('output_path', 'out', 'OUTPUT', v=('non_n_size',), req=('non_n_size',)),
            _Param('intervals', 'in', 'INTERVALS', v=('non_n_size',)),
        ),
        selector='mode',
        hook='base_counts',
    ),
    'gatk_extract_reference_sequences': _Spec(
        'gatk_extract_reference_sequences',
        {'interval_list': 'ExtractSequences', 'intervals': 'FastaReferenceMaker'},
        (
            _Param('reference_path', 'in', {'interval_list': 'REFERENCE_SEQUENCE', 'intervals': 'reference'}, req=True),
            _Param('output_path', 'out', {'interval_list': 'OUTPUT', 'intervals': 'output'}, req=True),
            _Param('interval_list', 'in', 'INTERVAL_LIST', v=('interval_list',), req=('interval_list',)),
            _Param('intervals', 'strs', 'intervals', v=('intervals',), req=('intervals',)),
            _Param('line_length', 'int', {'interval_list': 'LINE_LENGTH', 'intervals': 'line-width'}, minimum=1),
        ),
        selector='source',
    ),
    'gatk_fasta_alternate_reference_maker': _Spec(
        'gatk_fasta_alternate_reference_maker',
        {'default': 'FastaAlternateReferenceMaker'},
        (
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('snp_mask', 'in', 'snp-mask'),
            _Param('use_iupac_sample', 'str', 'use-iupac-sample'),
            _Param('line_width', 'int', 'line-width', minimum=1),
        ),
    ),
    'gatk_find_bad_genomic_kmers_spark': _Spec(
        'gatk_find_bad_genomic_kmers_spark',
        {'default': 'FindBadGenomicKmersSpark'},
        (
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('k_size', 'int', 'k-size', minimum=1),
            _Param('kmer_max_dust_score', 'int', 'kmer-max-dust-score', minimum=0),
            _Param('high_copy_fasta', 'in', 'high-copy-fasta'),
        ),
    ),
    'gatk_normalize_fasta': _Spec(
        'gatk_normalize_fasta',
        {'default': 'NormalizeFasta'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('line_length', 'int', 'LINE_LENGTH', minimum=1),
            _Param('truncate_sequence_names_at_whitespace', 'bool', 'TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE'),
        ),
    ),
    'gatk_scatter_intervals_by_ns': _Spec(
        'gatk_scatter_intervals_by_ns',
        {'default': 'ScatterIntervalsByNs'},
        (
            _Param('reference_path', 'in', 'REFERENCE', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('output_type', 'choice', 'OUTPUT_TYPE', choices=('N', 'ACGT', 'BOTH')),
            _Param('max_to_merge', 'int', 'MAX_TO_MERGE', minimum=0),
        ),
    ),
    'gatk_shift_fasta': _Spec(
        'gatk_shift_fasta',
        {'default': 'ShiftFasta'},
        (
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('shift_back_output', 'out', 'shift-back-output', req=True),
            _Param('shift_offset_list', 'strs', 'shift-offset-list'),
            _Param('interval_file_name', 'out', 'interval-file-name'),
            _Param('line_width', 'int', 'line-width', minimum=1),
        ),
    ),
}


def gatk_index_feature_file(
    feature_path: str,
    *,
    index_path: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a tribble index (.idx, or .tbi for bgzipped VCF) for a feature file.

    Run `gatk IndexFeatureFile`; without `index_path` the index file next to the input is detected after the run and reported as `outputs.index`.

    Parameters:
    - ``feature_path`` (required): Feature file (VCF, BED, GFF ...) in a tribble-supported format.
    - ``index_path`` (optional): Explicit index file; by default the index is written next to the feature file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.index` is the index.
    """
    return _core.execute(
        _SPECS['gatk_index_feature_file'],
        {'feature_path': feature_path, 'index_path': index_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sort_gff(
    input_path: str,
    output_path: str,
    *,
    sequence_dictionary: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort a GFF3 file and add the flush (`###`) directives.

    Run Picard `gatk SortGff`; contigs are ordered by `sequence_dictionary` when given.

    Parameters:
    - ``input_path`` (required): GFF3 file.
    - ``output_path`` (required): New sorted GFF3 file (with flush directives).
    - ``sequence_dictionary`` (optional): Sequence dictionary that defines the contig order.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the sorted GFF3.
    """
    return _core.execute(
        _SPECS['gatk_sort_gff'],
        {'input_path': input_path, 'output_path': output_path, 'sequence_dictionary': sequence_dictionary},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_gtf_to_bed(
    gtf_path: str,
    output_path: str,
    sequence_dictionary: str,
    *,
    sort_by_transcript: bool | None = None,
    use_basic_transcript: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a Gencode GTF file to a BED file.

    Run `gatk GtfToBed`.

    Parameters:
    - ``gtf_path`` (required): Gencode GTF file.
    - ``output_path`` (required): New BED file.
    - ``sequence_dictionary`` (required): Sequence dictionary (.dict) of the reference, used to order contigs.
    - ``sort_by_transcript`` (optional): Sort the output by transcript.
    - ``use_basic_transcript`` (optional): Keep only transcripts tagged `basic`.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the BED file.
    """
    return _core.execute(
        _SPECS['gatk_gtf_to_bed'],
        {'gtf_path': gtf_path, 'output_path': output_path, 'sequence_dictionary': sequence_dictionary, 'sort_by_transcript': sort_by_transcript, 'use_basic_transcript': use_basic_transcript},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_bait_designer(
    targets: str,
    reference_path: str,
    design_name: str,
    output_dir: str,
    *,
    bait_size: int | None = None,
    bait_offset: int | None = None,
    design_strategy: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design oligonucleotide baits for hybrid selection (exome capture) from target intervals.

    Run Picard `gatk BaitDesigner`; the files written to `output_dir` are listed in `files`.

    Parameters:
    - ``targets`` (required): Interval list of the target regions.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``design_name`` (required): Name of the bait design.
    - ``output_dir`` (required): New directory for the design files.
    - ``bait_size`` (optional): Length of each bait.
    - ``bait_offset`` (optional): Offset between consecutive baits.
    - ``design_strategy`` (optional): How baits are placed over a target.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `files` (file names inside `output_dir`).
    """
    return _core.execute(
        _SPECS['gatk_bait_designer'],
        {'targets': targets, 'reference_path': reference_path, 'design_name': design_name, 'output_dir': output_dir, 'bait_size': bait_size, 'bait_offset': bait_offset, 'design_strategy': design_strategy},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_bwa_mem_index_image(
    reference_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a BWA-MEM index image (`.img`) of a reference for the GATK BWA tools.

    Run `gatk BwaMemIndexImageCreator`; the image `<fasta>.img` is written next to the FASTA and must not exist yet; it is reported as `outputs.image`.

    Parameters:
    - ``reference_path`` (required): Reference FASTA.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.image` is the index image.
    """
    return _core.execute(
        _SPECS['gatk_bwa_mem_index_image'],
        {'reference_path': reference_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_check_reference_compatibility(
    references_to_compare: Sequence[str],
    *,
    alignment_path: str | None = None,
    variant_path: str | None = None,
    output_path: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a BAM/CRAM or VCF for compatibility with one or more reference FASTA files.

    Run `gatk CheckReferenceCompatibility` (experimental); at least one of `alignment_path` and `variant_path` is needed.

    Parameters:
    - ``references_to_compare`` (required): Reference FASTA files to test against.
    - ``alignment_path`` (optional): BAM/CRAM whose sequence dictionary is checked.
    - ``variant_path`` (optional): VCF whose contig header is checked.
    - ``output_path`` (optional): New report file (default: standard output).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields (the table is in `log_tail` unless `output_path` is given).
    """
    return _core.execute(
        _SPECS['gatk_check_reference_compatibility'],
        {'references_to_compare': references_to_compare, 'alignment_path': alignment_path, 'variant_path': variant_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_compare_references(
    reference_path: str,
    references_to_compare: Sequence[str],
    *,
    display_only_differing_sequences: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare reference FASTA files contig by contig (MD5, length, base differences) and summarise the differences.

    Run `gatk CompareReferences` (experimental) and parse the table it prints.

    Parameters:
    - ``reference_path`` (required): Main reference FASTA.
    - ``references_to_compare`` (required): Reference FASTA files to compare with.
    - ``display_only_differing_sequences`` (optional): List only sequences that differ.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields plus `comparison_table` (rows of the printed table).
    """
    return _core.execute(
        _SPECS['gatk_compare_references'],
        {'reference_path': reference_path, 'references_to_compare': references_to_compare, 'display_only_differing_sequences': display_only_differing_sequences},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_compose_str_table(
    reference_path: str,
    output_path: str,
    *,
    max_period: int | None = None,
    max_repeats: int | None = None,
    decimation: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compose a genome-wide short-tandem-repeat (STR) location table for DRAGstr model calibration.

    Run `gatk ComposeSTRTableFile`.

    Parameters:
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New STR table (binary).
    - ``max_period`` (optional): Maximum repeat unit length.
    - ``max_repeats`` (optional): Maximum repeat count.
    - ``decimation`` (optional): Decimation table (`period:factor,...`).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the STR table.
    """
    return _core.execute(
        _SPECS['gatk_compose_str_table'],
        {'reference_path': reference_path, 'output_path': output_path, 'max_period': max_period, 'max_repeats': max_repeats, 'decimation': decimation},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_count_reference_bases(
    fasta_path: str,
    *,
    mode: str = 'count_bases',
    output_path: str | None = None,
    intervals: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the bases of a FASTA file: the number of A, C, G, T and N, or the number of non-N bases.

    Run `gatk CountBasesInReference` (`mode=count_bases`, parsed into `base_counts`) or Picard `NonNFastaSize` (`mode=non_n_size`, result in the output file, optionally limited to an interval list).

    ``mode``: `count_bases` = CountBasesInReference, `non_n_size` = NonNFastaSize.

    Parameters:
    - ``fasta_path`` (required): FASTA file.
    - ``output_path`` (required for some modes): New file for the non-N size.
    - ``intervals`` (optional): Interval list limiting the count.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `base_counts` and `total_bases` for `count_bases`, `outputs.output_path` for `non_n_size`.
    """
    return _core.execute(
        _SPECS['gatk_count_reference_bases'],
        {'fasta_path': fasta_path, 'output_path': output_path, 'intervals': intervals},
        mode=mode,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_extract_reference_sequences(
    reference_path: str,
    output_path: str,
    *,
    source: str = 'interval_list',
    interval_list: str | None = None,
    intervals: Sequence[str] | None = None,
    line_length: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the sequence of intervals from a reference into a new FASTA file.

    Run Picard `ExtractSequences` (`source=interval_list`, intervals from an interval-list file) or `gatk FastaReferenceMaker` (`source=intervals`, intervals written as `chr:start-end`).

    ``source``: `interval_list` = ExtractSequences, `intervals` = FastaReferenceMaker.

    Parameters:
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New FASTA file.
    - ``interval_list`` (required for some modes): Interval list to extract (Picard format).
    - ``intervals`` (required for some modes): Intervals to extract, for example `chr1:100-200`.
    - ``line_length`` (optional): Line width of the output FASTA.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the FASTA.
    """
    return _core.execute(
        _SPECS['gatk_extract_reference_sequences'],
        {'reference_path': reference_path, 'output_path': output_path, 'interval_list': interval_list, 'intervals': intervals, 'line_length': line_length},
        mode=source,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_fasta_alternate_reference_maker(
    reference_path: str,
    variant_path: str,
    output_path: str,
    *,
    snp_mask: str | None = None,
    use_iupac_sample: str | None = None,
    line_width: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create an alternative reference by applying the variants of a VCF to a reference FASTA.

    Run `gatk FastaAlternateReferenceMaker`; restrict the regions with the `intervals` option.

    Parameters:
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``variant_path`` (required): Indexed VCF whose alleles are applied.
    - ``output_path`` (required): New alternate-reference FASTA.
    - ``snp_mask`` (optional): VCF of positions to mask with N.
    - ``use_iupac_sample`` (optional): Sample whose heterozygous genotypes are written as IUPAC codes.
    - ``line_width`` (optional): Line width of the output FASTA.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the FASTA.
    """
    return _core.execute(
        _SPECS['gatk_fasta_alternate_reference_maker'],
        {'reference_path': reference_path, 'variant_path': variant_path, 'output_path': output_path, 'snp_mask': snp_mask, 'use_iupac_sample': use_iupac_sample, 'line_width': line_width},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_find_bad_genomic_kmers_spark(
    reference_path: str,
    output_path: str,
    *,
    k_size: int | None = None,
    kmer_max_dust_score: int | None = None,
    high_copy_fasta: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the k-mers that occur at high frequency in a reference (for example to mask them in SV discovery).

    Run `gatk FindBadGenomicKmersSpark` (beta, Spark local mode).

    Parameters:
    - ``reference_path`` (required): Reference FASTA (a .2bit file also works).
    - ``output_path`` (required): New text file of over-represented k-mers.
    - ``k_size`` (optional): Length of the k-mers.
    - ``kmer_max_dust_score`` (optional): Maximum DUST score of a reported k-mer.
    - ``high_copy_fasta`` (optional): FASTA of known high-copy sequences to include.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the k-mer list.
    """
    return _core.execute(
        _SPECS['gatk_find_bad_genomic_kmers_spark'],
        {'reference_path': reference_path, 'output_path': output_path, 'k_size': k_size, 'kmer_max_dust_score': kmer_max_dust_score, 'high_copy_fasta': high_copy_fasta},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_normalize_fasta(
    input_path: str,
    output_path: str,
    *,
    line_length: int | None = None,
    truncate_sequence_names_at_whitespace: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rewrite a FASTA file so that all sequence lines have the same length.

    Run Picard `gatk NormalizeFasta`.

    Parameters:
    - ``input_path`` (required): FASTA file.
    - ``output_path`` (required): New normalised FASTA.
    - ``line_length`` (optional): Sequence line length (default 100).
    - ``truncate_sequence_names_at_whitespace`` (optional): Cut sequence names at the first whitespace.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new FASTA.
    """
    return _core.execute(
        _SPECS['gatk_normalize_fasta'],
        {'input_path': input_path, 'output_path': output_path, 'line_length': line_length, 'truncate_sequence_names_at_whitespace': truncate_sequence_names_at_whitespace},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_scatter_intervals_by_ns(
    reference_path: str,
    output_path: str,
    *,
    output_type: str | None = None,
    max_to_merge: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write an interval list made by splitting a reference at its runs of N.

    Run Picard `gatk ScatterIntervalsByNs`.

    Parameters:
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New interval list.
    - ``output_type`` (optional): Which intervals to write.
    - ``max_to_merge`` (optional): Merge N-runs of at most this length into the neighbouring ACGT interval.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the interval list.
    """
    return _core.execute(
        _SPECS['gatk_scatter_intervals_by_ns'],
        {'reference_path': reference_path, 'output_path': output_path, 'output_type': output_type, 'max_to_merge': max_to_merge},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_shift_fasta(
    reference_path: str,
    output_path: str,
    shift_back_output: str,
    *,
    shift_offset_list: Sequence[str] | None = None,
    interval_file_name: str | None = None,
    line_width: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a shifted FASTA (circular contigs such as chrM rotated by an offset) and the file that shifts coordinates back.

    Run `gatk ShiftFasta` (beta).

    Parameters:
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New shifted FASTA.
    - ``shift_back_output`` (required): New chain file that maps shifted coordinates back.
    - ``shift_offset_list`` (optional): Offset per contig, in the order of the dictionary.
    - ``interval_file_name`` (optional): New interval file for the shifted contigs.
    - ``line_width`` (optional): Line width of the output FASTA.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs` lists the FASTA and the shift-back file.
    """
    return _core.execute(
        _SPECS['gatk_shift_fasta'],
        {'reference_path': reference_path, 'output_path': output_path, 'shift_back_output': shift_back_output, 'shift_offset_list': shift_offset_list, 'interval_file_name': interval_file_name, 'line_width': line_width},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


__all__ = [
    "gatk_index_feature_file",
    "gatk_sort_gff",
    "gatk_gtf_to_bed",
    "gatk_bait_designer",
    "gatk_bwa_mem_index_image",
    "gatk_check_reference_compatibility",
    "gatk_compare_references",
    "gatk_compose_str_table",
    "gatk_count_reference_bases",
    "gatk_extract_reference_sequences",
    "gatk_fasta_alternate_reference_maker",
    "gatk_find_bad_genomic_kmers_spark",
    "gatk_normalize_fasta",
    "gatk_scatter_intervals_by_ns",
    "gatk_shift_fasta",
]

_register_public_verbose(globals(), tuple(__all__))
